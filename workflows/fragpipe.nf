/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE WORKFLOW (LFQ, TMT quant, TMT label check, any FragPipe workflow)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Unified FragPipe workflow; MSPROTEOMICS_TMT_LABELCHECK runs it and analyses its PSMs.
    When samples have labels (label column in samplesheet), TMT annotation is
    built automatically. When no labels are present, annotation is empty and
    FragPipe runs in LFQ mode.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { checkWorkflowDatabase; resolveDatabase; generateSdrfBookkeeping; resolveProteome; resolveOrganism; resolveFragPipeDatabase; resolveFragPipeDatabaseChannels; getTmtLabelcheckWorkflow } from '../subworkflows/local/utils_nfcore_msproteomics_pipeline'

// SDRF generation (bookkeeping only)
include { GENERATESDRFFROMSAMPLESHEET } from '../modules/local/generatesdrffromsamplesheet/main'

// FragPipe subworkflows
include { FRAGPIPE as FRAGPIPE_PIPELINE } from '../subworkflows/local/fragpipe/main'
include { FRAGPIPE_HEADLESS_WF             } from '../subworkflows/local/fragpipe_headless_wf/main'
include { FRAGPIPESDRF } from '../modules/local/fragpipesdrf/main'
include { FRAGPIPESAVEWORKFLOW } from '../modules/local/fragpipesaveworkflow/main'
include { WRITETEXTFILE as WRITE_FRAGPIPE_MANIFEST } from '../modules/local/writetextfile/main'
include { FRAGPIPESECONDPASS } from '../modules/local/fragpipesecondpass/main'
include { PMULTIQC as PMULTIQC_FRAGPIPE } from '../modules/bigbio/pmultiqc/main'

// FragPipe utilities
include { msfraggerConsoleLogs; subMzmlWeights; generateFragpipeManifest; fragpipeManifestEntries; parseSamplesheet; speclibOrDiannRuns; checkDiaPreprocessing; uniqueSpectraFiles; requireFragpipeLcmsFiles; tmtAnnotationLines; fragpipeWorkflowProperties; pmultiqcFragpipeInputs } from '../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MSPROTEOMICS_FRAGPIPE {

    main:
    // FragPipe runs every tool from one licensed distribution; there is no public default image (a stub run runs no tool)
    if (workflow.containerEngine && !workflow.stubRun && !params.fragpipe_container) {
        error("--mode fragpipe with ${workflow.containerEngine} needs --fragpipe_container: the licensed FragPipe image that runs every FragPipe tool")
    }
    // The label check is its own defined workflow: the shipped label check workflow of --tmt_type, no user workflow.
    // Every other run needs --fragpipe_workflow, there is no default (nf-core/msproteomics#2, user decisions 2026-10-16)
    if (params.tmt_mode == 'labelcheck' && params.fragpipe_workflow) {
        error("--tmt_mode labelcheck runs the shipped label check workflow of --tmt_type and takes no --fragpipe_workflow (see docs/usage.md)")
    }
    if (params.tmt_mode != 'labelcheck' && !params.fragpipe_workflow) {
        error("--mode fragpipe requires --fragpipe_workflow (a FragPipe .workflow file, e.g. the shipped assets/LFQ-MBR.workflow or assets/TMT6.workflow, see docs/usage.md)")
    }
    def fragpipe_workflow_file = file(params.fragpipe_workflow ?: getTmtLabelcheckWorkflow(params.tmt_type, projectDir), checkIfExists: true)
    def workflow_props = fragpipeWorkflowProperties(fragpipe_workflow_file)

    // One row per samplesheet line (TMT: one per channel); one item per spectra file.
    // SpecLibGen/DIA-NN on: every file is in FragPipe's single group "".
    def samplesheet_rows = parseSamplesheet(params.input, "${projectDir}/assets/schema_input.json", speclibOrDiannRuns(workflow_props))
    // Second pass of Write sub mzML: _sub.mzML inputs with their _percolator.weights (CmdMsfragger.java:612-626, CmdPercolator.java:185-201)
    def spectra_files    = subMzmlWeights(checkDiaPreprocessing(
        requireFragpipeLcmsFiles(uniqueSpectraFiles(samplesheet_rows)),
        workflow_props.getProperty('diaumpire.run-diaumpire') == 'true',
        workflow_props.getProperty('diatracer.run-diatracer') == 'true'
    ))

    // iProteome (--proteome): fasta and hosted FragPipe database; --database wins
    def proteome = resolveProteome(params)
    // the database FragPipe searches, resolved once (--database > iProteome)
    def fragpipe_database = resolveFragPipeDatabase(proteome, params)
    checkWorkflowDatabase(workflow_props, [resolveDatabase(proteome, params), fragpipe_database])
    // checks that can stop the run come before the first channel operator (see workflows/diann.nf)
    def db_channels = resolveFragPipeDatabaseChannels(proteome, params)

    GENERATESDRFFROMSAMPLESHEET(
        channel.fromPath(params.input, checkIfExists: true),
        generateSdrfBookkeeping(params, resolveOrganism(proteome, params)),
        file(fragpipe_database, glob: false, checkIfExists: true)
    )

    ch_workflow_file = channel.of([[id: 'fragpipe_config'], fragpipe_workflow_file])

    // TMT annotation (one plex per FragPipe group, see tmtAnnotationLines) in the multi-plex
    // TMT-Integrator format (plex\tchannel\tsample\tsample_name\tcondition\treplicate).
    // Empty ([]) when TMT-Integrator does not run, as FragPipe then uses no annotation.
    def tmt_annotation_lines = tmtAnnotationLines(samplesheet_rows, spectra_files, workflow_props)
    ch_annotation = channel.fromList(tmt_annotation_lines)
        .collectFile(name: 'tmt_annotation.tsv', newLine: true, sort: false)
        .ifEmpty([])

    ch_unique_files = channel.fromList(spectra_files)

    if (params.fragpipe_mode == 'headless') {
        //
        // ALL-IN-ONE MODE: Run FragPipe headless in a single process
        //
        // FragPipe headless reads every LC-MS file as given (.raw natively).
        // Generates a manifest from the samplesheet for FragPipe's --manifest flag.
        // For TMT: generates per-experiment annotation content for auto-discovery.
        //
        def headless_files = spectra_files.collect { meta, spectra_file -> [meta, spectra_file.name] }

        // Manifest rows: each spectra file with its experiment, bioreplicate, data type and FragPipe group;
        // FRAGPIPEHEADLESS puts each TMT group's files in one directory with its annotation (headlessLayout)
        ch_manifest_rows = channel.value(
            headless_files.collect { meta, fname -> [fname, meta.experiment, meta.bioreplicate, meta.data_type, meta.group] }
        )

        // TMT annotation for headless mode: group, channel, sample_name per channel, the same rows as the
        // pipeline-mode annotation; the blank group "" stays "". [] for LFQ.
        ch_annotation_rows = channel.value(
            tmt_annotation_lines.drop(1).collect { line ->
                def cols = line.split('\t')
                [cols[0], cols[1], cols[3]]
            }
        )

        // Workflow file (plain path, not tuple)
        ch_wf_file = channel.of(fragpipe_workflow_file)

        // FragPipe headless never generates decoys; FRAGPIPE_HEADLESS_WF checks the database
        // (and runs the database update first with --add_philosopher_decoy / --add_philosopher_contaminants).
        FRAGPIPE_HEADLESS_WF(
            ch_unique_files,
            db_channels.ch_fasta,
            db_channels.ch_prebuilt_db,
            ch_wf_file,
            ch_manifest_rows,
            ch_annotation_rows
        )
        // FragPipe's results: each group's Report psm.tsv (<group>/psm.tsv, psm.tsv for the group "") and
        // IonQuant's combined_modified_peptide.tsv in the work dir root
        ch_headless_results = FRAGPIPE_HEADLESS_WF.out.all_results.flatten()
        ch_group_psms = ch_headless_results
            .filter { result -> result.name == 'psm.tsv' }
            .toSortedList { a_psm, b_psm -> a_psm.parent.name <=> b_psm.parent.name }
        ch_combined_modified_peptide = ch_headless_results
            .filter { result -> result.name == 'combined_modified_peptide.tsv' && result.parent.name == 'results' }
            .map { result -> [[id: 'fragpipe'], result] }
        // pmultiqc: FragPipe's work dir root files and each group's directory
        ch_pmultiqc_files = FRAGPIPE_HEADLESS_WF.out.all_results.map { results ->
            def files = [results].flatten()
            pmultiqcFragpipeInputs(
                files.findAll { result -> result.parent.name == 'results' },
                files.findAll { result -> result.name == 'psm.tsv' && result.parent.name != 'results' }.collect { psm -> psm.parent }
            )
        }
    } else {
        //
        // PIPELINE MODE: Modular FragPipe subworkflows with Nextflow orchestration
        //
        // FragPipe run-folder files (FragpipeRun.java:371-376): the manifest of the LC-MS files
        // (TabWorkflow.manifestSave: absolute path, experiment, bioreplicate, data type; no trailing
        // newline) and FragPipe's saved UI state. Published to <outdir>/fragpipe/ like headless mode.
        // Read by DIA-NN's post-steps and MBG too.
        WRITE_FRAGPIPE_MANIFEST(channel.of([[id: 'fragpipe'], 'fragpipe-files.fp-manifest', generateFragpipeManifest(fragpipeManifestEntries(spectra_files))]))
        ch_fragpipe_manifest = WRITE_FRAGPIPE_MANIFEST.out.file.map { _meta, manifest -> manifest }.first()
        // fragpipe.workflow: FragPipe's saved UI state for this run, written by FragPipe's own code.
        // FragPipe loads the manifest and drops the LC-MS files that do not exist (TabWorkflow.manifestLoad),
        // so the task loads a manifest naming the files it stages (spectra/<name>), not the input paths
        // (an S3 input is not a path inside the task). The database is staged too.
        def fragpipe_workdir = file("${params.outdir}/fragpipe").toString()
        def fragpipe_db = file(fragpipe_database, glob: false, checkIfExists: true)
        def ui_state_dir = file("${projectDir}/assets/fragpipe_ui_state", checkIfExists: true)
        def ch_staged_manifest = channel.of(generateFragpipeManifest(fragpipeManifestEntries(spectra_files, 'spectra')))
            .collectFile(name: 'staged.fp-manifest')
        FRAGPIPESAVEWORKFLOW(
            ch_staged_manifest.map { manifest ->
                [[id: 'fragpipe'], fragpipe_workflow_file, manifest, ui_state_dir, fragpipe_db, fragpipe_workdir, spectra_files.collect { _meta, spectra_file -> spectra_file }]
            }
        )

        FRAGPIPE_PIPELINE(
            ch_unique_files,
            db_channels.ch_fasta,
            ch_workflow_file,
            db_channels.ch_prebuilt_db,
            ch_annotation,
            ch_fragpipe_manifest,
            // Metaproteomics taxonomy dumps; [] = FragPipe's tools/metaproteomics/{names,nodes}.dmp (CmdMetaproteomics.java:39-42)
            channel.value([
                params.metaproteomics_names_dmp ? file(params.metaproteomics_names_dmp, checkIfExists: true) : [],
                params.metaproteomics_nodes_dmp ? file(params.metaproteomics_nodes_dmp, checkIfExists: true) : []
            ])
        )

        // each group's final Report psm.tsv, groups by name; IonQuant's combined_modified_peptide.tsv
        ch_group_psms = FRAGPIPE_PIPELINE.out.results_dir
            .map { _meta, dir -> dir }
            .collect(sort: true)
            .map { dirs -> dirs.sort(false) { dir -> dir.name }.collect { dir -> dir.resolve('psm.tsv') } }
        ch_combined_modified_peptide = FRAGPIPE_PIPELINE.out.combined_modified_peptide

        // pmultiqc: the whole work dir of the last IonQuant pass (its combined_*.tsv, the blank group's tables and each
        // group's rewritten psm.tsv and ion.tsv, in FragPipe's layout), else each group's final Report directory;
        // with fragpipe.workflow and the manifest
        ch_pmultiqc_files = FRAGPIPE_PIPELINE.out.quant_dir.map { dir -> [1, [dir]] }
            .mix(FRAGPIPE_PIPELINE.out.results_dir.map { _meta, dir -> dir }.collect(sort: true).map { dirs -> [0, dirs] })
            .toList()
            .map { sources -> sources.max { source -> source[0] } }
            .combine(FRAGPIPESAVEWORKFLOW.out.workflow.map { _meta, saved_workflow -> saved_workflow })
            .combine(ch_fragpipe_manifest)
            .map { _source, dirs, saved_workflow, manifest -> pmultiqcFragpipeInputs([saved_workflow, manifest], dirs) }

        // MSFragger console text in the order FragPipe runs MSFragger (FRAGPIPESDRF, FRAGPIPESECONDPASS)
        def ch_msfragger_logs = FRAGPIPE_PIPELINE.out.msfragger_logs
            .toList()
            .map { items -> msfraggerConsoleLogs(items, spectra_files) }

        // Write sub mzML: FragPipe's finalizer writes the second-pass workflow and manifest from
        // fragpipe.workflow, the manifest and the console text (FragpipeRun.java:589-700); the
        // _sub.mzML files are published to the same fragpipe/ folder (WRITESUBMZML)
        if (workflow_props.getProperty('tab-run.write_sub_mzml', 'false').toLowerCase() == 'true') {
            FRAGPIPESECONDPASS(
                FRAGPIPESAVEWORKFLOW.out.workflow
                    .combine(ch_fragpipe_manifest)
                    .combine(ch_msfragger_logs.map { logs -> [logs] })
                    .map { _meta, saved_workflow, manifest, logs -> [[id: 'fragpipe'], saved_workflow, manifest, logs, fragpipe_workdir] }
            )
        }

        // FragPipe saves the SDRF template after the run when workflow.misc.save-sdrf=true
        // (FragpipeRun.java:582-587; the checkbox defaults to on, TabRun.java:429) from the run
        // database (FragpipeRun.java:439-442), the MSFragger console text in MSFragger's run order
        // (msfraggerConsoleLogs) and the TMT-Integrator channel type when
        // TMT-Integrator runs. Instruments are keyed by the file MSFragger searched
        // (TabMsfragger.java:2388), so the manifest names are mapped to the searched names.
        if (workflow_props.getProperty('workflow.misc.save-sdrf', 'true') == 'true') {
            def tmt_label = workflow_props.getProperty('tmtintegrator.run-tmtintegrator') == 'true'
                ? workflow_props.getProperty('tmtintegrator.channel_num')
                : ''
            def sdrf_manifest_names = spectra_files.collectEntries { meta, spectra_file -> [meta.id, spectra_file.name] }
            def ch_sdrf_searched_names = FRAGPIPE_PIPELINE.out.searched_lcms
                // DIA-Umpire / diaTracer pseudo spectra have no manifest name of their own: they map to themselves
                .map { meta, searched -> "${sdrf_manifest_names[meta.id] ?: searched.name}\t${searched.name}" }
                .collectFile(name: 'searched_files.tsv', newLine: true, sort: true)
            FRAGPIPESDRF(
                ch_fragpipe_manifest
                    .combine(FRAGPIPE_PIPELINE.out.database)
                    .combine(ch_sdrf_searched_names)
                    .combine(ch_msfragger_logs.map { logs -> [logs] })
                    .combine(channel.value(file("${projectDir}/assets/unimod.obo", checkIfExists: true)))
                    .combine(FRAGPIPESAVEWORKFLOW.out.workflow.map { _meta, saved_workflow -> saved_workflow })
                    .map { manifest, fasta, searched_names, logs, unimod_obo, workflow_file -> [[id: 'fragpipe'], workflow_file, manifest, searched_names, fasta, unimod_obo, logs, tmt_label] }
            )
        }
    }

    // Proteomics QC report (pmultiqc FragPipe plugin); multiqc_config.yml is staged with the results (ext.args --config)
    if (params.enable_pmultiqc) {
        def pmultiqc_config = file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        PMULTIQC_FRAGPIPE(ch_pmultiqc_files.map { inputs -> inputs + [pmultiqc_config] })
    }

    emit:
    group_psms                = ch_group_psms                 // channel: [ path(psm.tsv), ... ] - each group's Report psm.tsv, groups by name
    combined_modified_peptide = ch_combined_modified_peptide  // channel: [ val(meta), path(combined_modified_peptide.tsv) ] - IonQuant's, when it runs
    multiqc_files             = channel.empty()               // every module reports its version on the 'versions' topic (main.nf collects them)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
