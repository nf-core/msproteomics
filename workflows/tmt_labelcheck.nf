/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TMT LABEL CHECK WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { checkWorkflowDatabase; resolveDatabase; generateSdrfBookkeeping; resolveProteome; resolveOrganism; resolveFragPipeDatabase; resolveFragPipeDatabaseChannels; getTmtLabelcheckWorkflow } from '../subworkflows/local/utils_nfcore_msproteomics_pipeline'

// SDRF generation (bookkeeping only)
include { GENERATESDRFFROMSAMPLESHEET } from '../modules/local/generatesdrffromsamplesheet/main'

// TMT Label Check subworkflow
include { TMT_LABELCHECK } from '../subworkflows/local/tmt_labelcheck/main'
include { FRAGPIPESDRF } from '../modules/local/fragpipesdrf/main'
include { FRAGPIPESAVEWORKFLOW } from '../modules/local/fragpipesaveworkflow/main'
include { FRAGPIPESECONDPASS } from '../modules/local/fragpipesecondpass/main'

// Samplesheet parsing (shared with the FragPipe workflow)
include { msfraggerConsoleLogs; subMzmlWeights; parseSamplesheet; uniqueSpectraFiles; requireFragpipeLcmsFiles; speclibOrDiannRuns; fragpipeWorkflowProperties; generateFragpipeManifest; fragpipeManifestEntries } from '../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MSPROTEOMICS_TMT_LABELCHECK {

    main:
    // FragPipe runs every tool from one licensed distribution; there is no public default image (a stub run runs no tool)
    if (workflow.containerEngine && !workflow.stubRun && !params.fragpipe_container) {
        error("--mode fragpipe with ${workflow.containerEngine} needs --fragpipe_container: the licensed FragPipe image that runs every FragPipe tool")
    }
    // One search per spectra file: TMT per-channel samplesheets have one row per channel
    // (FragpipeRun.java:1258-1262 takes each LC-MS file once).
    // Auto-select TMT labelcheck workflow file based on tmt_type if not explicitly provided
    def tmt_workflow = file(params.fragpipe_workflow ?: getTmtLabelcheckWorkflow(params.tmt_type, projectDir), checkIfExists: true)
    // SpecLibGen/DIA-NN on: every file is in FragPipe's single group "" (RC-04).
    def workflow_props = fragpipeWorkflowProperties(tmt_workflow)
    def samplesheet_rows = parseSamplesheet(params.input, "${projectDir}/assets/schema_input.json", speclibOrDiannRuns(workflow_props))
    // Second pass of Write sub mzML: _sub.mzML inputs with their _percolator.weights (CmdMsfragger.java:612-626, CmdPercolator.java:185-201)
    def spectra_files = subMzmlWeights(requireFragpipeLcmsFiles(uniqueSpectraFiles(samplesheet_rows)))
    ch_input = channel.fromList(spectra_files)

    // iProteome (--proteome): fasta and hosted FragPipe database; --database wins
    def proteome = resolveProteome(params)
    checkWorkflowDatabase(fragpipeWorkflowProperties(tmt_workflow), [resolveDatabase(proteome, params), resolveFragPipeDatabase(proteome, params)])

    GENERATESDRFFROMSAMPLESHEET(
        channel.fromPath(params.input, checkIfExists: true),
        generateSdrfBookkeeping(params, resolveOrganism(proteome, params)),
        file(resolveFragPipeDatabase(proteome, params), glob: false, checkIfExists: true)
    )

    def db_channels = resolveFragPipeDatabaseChannels(proteome, params)

    ch_workflow_file = channel.of([[id: 'labelcheck_workflow'], tmt_workflow])

    TMT_LABELCHECK(
        ch_input,
        db_channels.ch_fasta,
        ch_workflow_file,
        params.tmt_type,
        db_channels.ch_prebuilt_db
    )

    // MSFragger console text in the order FragPipe runs MSFragger (FRAGPIPESDRF, FRAGPIPESECONDPASS)
    def ch_msfragger_logs = TMT_LABELCHECK.out.msfragger_logs
        .toList()
        .map { items -> msfraggerConsoleLogs(items, spectra_files) }

    // FragPipe writes fragpipe-files.fp-manifest and fragpipe.workflow (its saved UI state) into the work dir
    // of every run (FragpipeRun.java:371-375), as in MSPROTEOMICS_FRAGPIPE; the SDRF and the write-sub-mzML
    // finalizer are written from them
    def fragpipe_workdir = file("${params.outdir}/fragpipe").toString()
    FRAGPIPESAVEWORKFLOW(
        channel.of(generateFragpipeManifest(fragpipeManifestEntries(spectra_files, 'spectra')))
            .collectFile(name: 'staged.fp-manifest')
            .map { manifest ->
                [[id: 'fragpipe'], tmt_workflow, manifest, file("${projectDir}/assets/fragpipe_ui_state", checkIfExists: true),
                 file(resolveFragPipeDatabase(proteome, params), glob: false, checkIfExists: true), fragpipe_workdir, spectra_files.collect { _meta, spectra_file -> spectra_file }]
            }
    )

    // Write sub mzML: FragPipe's finalizer writes the second-pass workflow and manifest from fragpipe.workflow,
    // the manifest and the console text (FragpipeRun.java:589-700)
    if (workflow_props.getProperty('tab-run.write_sub_mzml', 'false').toLowerCase() == 'true') {
        FRAGPIPESECONDPASS(
            FRAGPIPESAVEWORKFLOW.out.workflow
                .combine(channel.of(generateFragpipeManifest(fragpipeManifestEntries(spectra_files))).collectFile(name: 'fragpipe-files.fp-manifest'))
                .combine(ch_msfragger_logs.map { logs -> [logs] })
                .map { _meta, saved_workflow, manifest, logs -> [[id: 'fragpipe'], saved_workflow, manifest, logs, fragpipe_workdir] }
        )
    }

    // FragPipe saves the SDRF template after the run when workflow.misc.save-sdrf=true (FragpipeRun.java:582-587;
    // the checkbox defaults to on, TabRun.java:429-430), as in workflows/fragpipe.nf.
    if (workflow_props.getProperty('workflow.misc.save-sdrf', 'true') == 'true') {
        def tmt_label = workflow_props.getProperty('tmtintegrator.run-tmtintegrator') == 'true'
            ? workflow_props.getProperty('tmtintegrator.channel_num')
            : ''
        def ch_manifest = channel.of(generateFragpipeManifest(fragpipeManifestEntries(spectra_files)))
            .collectFile(name: 'fragpipe-files.fp-manifest')
        // FragPipe writes the SDRF from the MSFragger console text in MSFragger's run order (msfraggerConsoleLogs)
        // and keys instruments by the file MSFragger searched (TabMsfragger.java:2388), so the manifest names are
        // mapped to the searched names
        def sdrf_manifest_names = spectra_files.collectEntries { meta, spectra_file -> [meta.id, spectra_file.name] }
        def ch_sdrf_searched_names = TMT_LABELCHECK.out.searched_lcms
            .map { meta, searched -> "${sdrf_manifest_names[meta.id]}\t${searched.name}" }
            .collectFile(name: 'searched_files.tsv', newLine: true, sort: true)
        FRAGPIPESDRF(
            ch_manifest
                .combine(TMT_LABELCHECK.out.database)
                .combine(ch_sdrf_searched_names)
                .combine(ch_msfragger_logs.map { logs -> [logs] })
                .combine(channel.value(file("${projectDir}/assets/unimod.obo", checkIfExists: true)))
                .combine(FRAGPIPESAVEWORKFLOW.out.workflow.map { _meta, saved_workflow -> saved_workflow })
                .map { manifest, fasta, searched_names, logs, unimod_obo, workflow_file -> [[id: 'fragpipe'], workflow_file, manifest, searched_names, fasta, unimod_obo, logs, tmt_label] }
        )
    }

    emit:
    versions      = channel.empty() // every module reports its version on the 'versions' topic
    multiqc_files = channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
