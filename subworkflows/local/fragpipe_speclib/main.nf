/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_SPECLIB SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Spectral library generation and DIA-NN quantification (FragPipe's SpecLibGen and
    DIA-Quant steps, FragpipeRun.java:2240-2339; DIA-NN runs after SpecLibGen, :2486).

    Modules:
    - SPECLIBGEN     (FragPipe's tools/speclib/gen_con_spec_lib.py, CmdSpecLibGen.java)
    - DIANNFRAGPIPE (FragPipe's bundled DIA-NN and its post-steps, CmdDiann.java)

    Execution: Aggregate (one FragPipe group "" when either runs, RC-04)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPECLIBGEN     } from '../../../modules/local/speclibgen/main'
include { DIANNFRAGPIPE } from '../../../modules/local/diannfragpipe/main'
include { WRITETEXTFILE as WRITE_DIANN_EXPERIMENT_ANNOTATION } from '../../../modules/local/writetextfile/main'
include { shouldRunTool; getToolEntry; speclibSpectraFiles; diannSpectraFiles; fragpipeManifestEntries; lfqExperimentAnnotation; fpopCoadaptrMeta } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'
include { diannSizeOrder } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

workflow FRAGPIPE_SPECLIB {
    take:
    ch_psms          // channel: [ val(meta), path(psm.tsv) ] - the group's psm.tsv (one group)
    ch_peptides      // channel: [ val(meta), path(peptide.tsv) ] - the group's peptide.tsv (one group)
    ch_proteins      // channel: [ val(meta), path(protein.tsv) ] - the group's protein.tsv (FragReporter)
    ch_pepxml        // channel: [ val(meta), path(pepxml) ] - validated pepXML per searched file
    ch_searched      // channel: [ val(meta), path(spectra) ] - files MSFragger searched (after DIA-Umpire / diaTracer)
    ch_uncalibrated  // channel: [ val(meta), path(*_uncalibrated.mzML) ] - MSFragger's uncalibrated mzML
    ch_files         // channel: [ val(meta), path(spectra) ] - LC-MS table (samplesheet files, meta.data_type)
    ch_manifest      // channel: path(fragpipe-files.fp-manifest) - value channel
    ch_fasta         // channel: path(fasta) - value channel for database
    ch_tool_configs  // channel: val(tool_configs_json) - parsed JSON config (value channel)
    aggregate_meta   // val(meta) for aggregate outputs

    main:

    def wrap = { items -> [items] }  // keep collected lists intact through combine

    // Only checked when SpecLibGen / DIA-NN run: they force FragPipe's single group (RC-04)
    def single = { files, name ->
        if (files.size() != 1) {
            error("SpecLibGen / DIA-NN use FragPipe's single group (RC-04) but ${files.size()} ${name} files were produced")
        }
        files[0]
    }
    ch_group_psm = ch_psms.map { _meta, psm -> psm }.toList().map(wrap)
    ch_group_peptide = ch_peptides.map { _meta, peptide -> peptide }.toList().map(wrap)
    ch_group_protein = ch_proteins.map { _meta, protein -> protein }.toList().map(wrap)
    // LC-MS table in samplesheet order (channel.fromList emits it in order), so task inputs and hashes are stable for -resume
    ch_table = ch_files.toList().map { files -> wrap.call(files) }
    ch_uncalibrated_by_id = ch_uncalibrated.toList().map { items -> wrap.call(items.collectEntries { meta, mzml -> [(meta.id): mzml] }) }

    //
    // SPECLIBGEN: group work dir with psm.tsv, peptide.tsv, pepXMLs and spectra (CmdSpecLibGen.java:114-216)
    //
    ch_speclib_spectra = ch_searched.toSortedList { a, b -> a[0].id <=> b[0].id }.map(wrap)
        .combine(ch_uncalibrated_by_id)
        .combine(ch_table)
        .map { searched, uncalibrated, files ->
            def spectra = speclibSpectraFiles(searched, uncalibrated, files)
            wrap.call(spectra)
        }

    ch_for_speclibgen = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'speclibgen') }
        .combine(ch_group_psm)
        .combine(ch_group_peptide)
        .combine(ch_pepxml.flatMap { _meta, pepxml -> [pepxml].flatten() }.collect(sort: true).map(wrap))
        .combine(ch_speclib_spectra)
        .map { configs, psms, peptides, pepxmls, spectra ->
            def psm = single.call(psms, 'psm.tsv')
            def peptide = single.call(peptides, 'peptide.tsv')
            def entry = getToolEntry(configs, 'speclibgen')
            def rt_cal_file = entry.rt_cal_file ? file(entry.rt_cal_file, checkIfExists: true) : []
            def im_cal_file = entry.im_cal_file ? file(entry.im_cal_file, checkIfExists: true) : []
            [aggregate_meta, psm, peptide, pepxmls, spectra, rt_cal_file, im_cal_file, entry.argv_parts]
        }

    SPECLIBGEN(ch_for_speclibgen, ch_fasta)

    //
    // DIANNFRAGPIPE: runs when DIA-NN is on and the LC-MS table has DIA or DIA-Quant files
    // (FragpipeRun.java:2301-2306), on those files (CmdDiann.java:209-238).
    //
    ch_diann_spectra = ch_table
        .combine(ch_uncalibrated_by_id)
        .combine(ch_tool_configs)
        .map { files, uncalibrated, configs ->
            def spectra = diannSpectraFiles(files, uncalibrated, shouldRunTool(configs, 'diaumpire'))
            wrap.call(diannSizeOrder(spectra))
        }

    ch_experiment_annotation = ch_table.map { wrapped ->
        def entries = fragpipeManifestEntries(wrapped[0])
        lfqExperimentAnnotation(entries)
    }

    // --lib: the workflow's diann.library, else SpecLibGen's library.tsv (CmdDiann.java:315-326)
    ch_speclib_library = SPECLIBGEN.out.library_tsv.map { _meta, library -> library }.ifEmpty([]).map(wrap)

    ch_for_diann = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'diann') }
        .combine(ch_diann_spectra)
        .filter { _configs, spectra -> !spectra.isEmpty() }
        .combine(ch_speclib_library)
        .combine(ch_group_psm)
        .combine(ch_group_peptide)
        .combine(ch_group_protein)
        .combine(ch_manifest)
        .combine(ch_experiment_annotation)
        .map { configs, spectra, speclib_library, psms, peptides, proteins, manifest, experiment_annotation ->
            def psm = single.call(psms, 'psm.tsv')
            def peptide = single.call(peptides, 'peptide.tsv')
            def protein = single.call(proteins, 'protein.tsv')
            def entry = getToolEntry(configs, 'diann')
            // External library (diann.library, DiannPanel "Spectral library (optional)"): readability
            // check and message of CmdDiann.java:305-311; MSstats is then off (:137-140, parser)
            def library = speclib_library
            if (entry.library) {
                library = file(entry.library)
                if (!library.exists() || !java.nio.file.Files.isRegularFile(library) || !java.nio.file.Files.isReadable(library)) {
                    error("The custom library file does not exist or is not readable: ${entry.library}\nPlease double check the spectral library (optinal) box in the Quant (DIA) tab.")
                }
                log.info("There are external spectral library. Will not generate the MSstats file.")
            }
            if (!library) {
                error("DIA-NN reads SpecLibGen's library.tsv (CmdDiann.java:313-326) but SpecLibGen does not run and diann.library is empty")
            }
            // FPOP coADAPTr's DIA mode edits DIA-NN's report.tsv in place (meta.fpop_coadaptr, conf/modules.config)
            def coadaptr_meta = fpopCoadaptrMeta(aggregate_meta, configs, ['diann'])
            [coadaptr_meta, library, spectra, psm, peptide, protein, manifest, experiment_annotation, groovy.json.JsonOutput.toJson(entry)]
        }

    DIANNFRAGPIPE(ch_for_diann.map { item -> item.take(7) + [item[8]] }, ch_fasta)

    // experiment_annotation.tsv of the DIA-Quant work dir (generateLFQExperimentAnnotation(wd, 1), CmdDiann.java:610-617)
    WRITE_DIANN_EXPERIMENT_ANNOTATION(ch_for_diann.map { item -> [item[0], 'experiment_annotation.tsv', item[7]] })

    emit:
    library_tsv       = SPECLIBGEN.out.library_tsv         // channel: [ val(meta), path(library.tsv) ]
    easypqp_files     = SPECLIBGEN.out.easypqp_files       // channel: [ val(meta), path(easypqp_files) ]
    rt_reference      = SPECLIBGEN.out.rt_reference        // channel: [ val(meta), path(easypqp_rt_reference_run.tsv) ]
    diann_results     = DIANNFRAGPIPE.out.results_files   // channel: [ val(meta), [ path(dia-quant-output/*) ] ]
    diann_report      = DIANNFRAGPIPE.out.report          // channel: [ val(meta), path(report.tsv) ]
    diann_pg_matrix   = DIANNFRAGPIPE.out.pg_matrix       // channel: [ val(meta), path(report.pg_matrix.tsv) ]
    diann_pr_matrix   = DIANNFRAGPIPE.out.pr_matrix       // channel: [ val(meta), path(report.pr_matrix.tsv) ]
    diann_speclib     = DIANNFRAGPIPE.out.speclib         // channel: [ val(meta), path(report-tsv.speclib) ]
    diann_msstats     = DIANNFRAGPIPE.out.msstats         // channel: [ val(meta), path(msstats.csv) ]
    diann_experiment_annotation = WRITE_DIANN_EXPERIMENT_ANNOTATION.out.file  // channel: [ val(meta), path(experiment_annotation.tsv) ]
    versions          = SPECLIBGEN.out.versions_speclibgen
        .mix(SPECLIBGEN.out.versions_fragpipe)
        .mix(DIANNFRAGPIPE.out.versions_diann)
        .mix(DIANNFRAGPIPE.out.versions_fragpipe)
}
