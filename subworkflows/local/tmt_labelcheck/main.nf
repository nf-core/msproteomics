/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TMT LABEL CHECK SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A nf-core compliant subworkflow for TMT labeling efficiency QC.

    This subworkflow uses a FragPipe workflow file parsed by PARSEFRAGPIPEWORKFLOW
    and reuses generic FragPipe subworkflows for maximum code reuse:

    Pipeline Steps:
    1. PARSEFRAGPIPEWORKFLOW: Parse workflow file to get tool configs
    2. FRAGPIPE_DATABASE: FragPipe's decoy check (decoys/contaminants only with --add_philosopher_decoy / --add_philosopher_contaminants)
    3. CHECKCENTROID + FRAGPIPE_CONVERT: centroid check, optional DIA preprocessing (no conversion)
    4. FRAGPIPE_SEARCH: MSFragger (+ Crystal-C / MSBooster when the workflow enables them)
       - Single-pass search (no three-stage calibration) for faster QC
       - TMT is set as VARIABLE modification to detect unlabeled peptides
       - Protein N-terminal acetylation included to account for biological blocking
    5. FRAGPIPE_VALIDATE: Percolator PSM-level FDR control (parallel per sample)
    6. FRAGPIPE_INFERENCE: ProteinProphet + Philosopher Filter (aggregate → per-sample)
    7. FRAGPIPE_QUANT: IonQuant aggregate per-sample PSMs (aggregate)
       - Produces combined_modified_peptide.tsv for combined analysis
    8a. TMTLABELCHECKANALYZE (mode=psm): Per-sample labeling efficiency (aggregate)
    8b. TMTLABELCHECKANALYZE (mode=ionquant): Combined labeling efficiency (aggregate)

    Labeling Efficiency Formula:
        Total_Sites = peptide.count('K') + 1  (lysines + N-terminus)
        Labeled_Sites = count(TMT modifications)
        Efficiency = sum(Labeled_Sites) / sum(Total_Sites) * 100
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import config parser
include { PARSEFRAGPIPEWORKFLOW } from '../../../modules/local/parsefragpipeworkflow/main'

// Import modules used directly (not through subworkflows)
include { CHECKCENTROID } from '../../../modules/local/checkcentroid/main'
include { TMTLABELCHECKANALYZE } from '../../../modules/local/tmtlabelcheckanalyze/main'
include { TMTLABELCHECKANALYZE as TMTLABELCHECKANALYZE_IONQUANT } from '../../../modules/local/tmtlabelcheckanalyze/main'

// Import generic FragPipe subworkflows
include { FRAGPIPE_DATABASE  } from '../fragpipe_database/main'
include { FRAGPIPE_CONVERT   } from '../fragpipe_convert/main'
include { FRAGPIPE_SEARCH    } from '../fragpipe_search/main'
include { FRAGPIPE_VALIDATE  } from '../fragpipe_validate/main'
include { FRAGPIPE_INFERENCE } from '../fragpipe_inference/main'
include { FRAGPIPE_QUANT     } from '../fragpipe_quant/main'
include { FRAGPIPE_MATCHED_FRAGMENTS_INPUTS } from '../fragpipe_matched_fragments/main'
include { FRAGPIPE_PHILOSOPHER_QUANT } from '../fragpipe_philosopher_quant/main'
include { WRITESUBMZML     } from '../../../modules/local/writesubmzml/main'
include { checkCentroidFile; shouldRunTool; requireToolField; msboosterFraggerParamsType; getMsfraggerNumSlices } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TMT_LABELCHECK {
    take:
    ch_file            // channel: [ val(sample_meta), path(mzml/raw) ] - sample files with metadata
    ch_database        // channel: [ val(meta), path(fasta) ] - FASTA to update with decoys (empty unless --add_philosopher_decoy / --add_philosopher_contaminants)
    ch_workflow_file   // channel: [ val(meta), path(workflow_file) ] - FragPipe workflow file for label check
    tmt_type           // val: TMT type string (TMT6, TMT10, TMT16, TMTPRO, etc.)
    ch_decoy_db        // channel: [ val(meta), path(fasta) ] - FASTA searched as given (empty when ch_database is used)

    main:

    //
    // STEP 1: Parse workflow config to get tool run flags and args
    //
    // FragPipe adjusts the tools to the data types of the manifest files (TabWorkflow.adjustToolsBasedOnDataTypes)
    ch_data_types = ch_file
        .map { meta, _spectra -> meta.data_type }
        .unique()
        .toSortedList()
    PARSEFRAGPIPEWORKFLOW(ch_workflow_file.combine(ch_data_types.map { data_types -> [data_types] }))

    // Read tool configs JSON as a string value channel.
    // Strings avoid the Groovy Map-as-named-parameter bug in channel operators.
    // Use splitText to avoid blocking I/O inside map{} on cloud executors.
    ch_tool_configs = PARSEFRAGPIPEWORKFLOW.out.tool_configs_json
        .map { _meta, json_file -> json_file }
        .splitText()
        .collect()
        .map { lines -> lines.join('') }
        .first()

    // MSFragger params per data type (CmdMsfragger.java:444-457)
    ch_fragger_params = PARSEFRAGPIPEWORKFLOW.out.msfragger_config.map { _meta, params -> ['DDA', params] }
        .mix(PARSEFRAGPIPEWORKFLOW.out.msfragger_dia_config.map { _meta, params -> ['DIA', params] })
        .mix(PARSEFRAGPIPEWORKFLOW.out.msfragger_gpfdia_config.map { _meta, params -> ['GPF-DIA', params] })
        .mix(PARSEFRAGPIPEWORKFLOW.out.msfragger_dda_plus_config.map { _meta, params -> ['DDA+', params] })
        .toList()


    //
    // STEP 2: FragPipe's database decoy check (CHECKDBCONFIG); ch_database only holds a
    // FASTA when --add_philosopher_decoy / --add_philosopher_contaminants asks for FragPipe's database update first.
    //
    FRAGPIPE_DATABASE(ch_database, ch_decoy_db, ch_tool_configs)
    ch_fasta = FRAGPIPE_DATABASE.out.fasta_path

    //
    // STEP 3: FragPipe's CheckCentroid on the first LC-MS file (FragpipeRun.java:1333-1340), then
    // optional DIA preprocessing; LC-MS files are searched as given (MSFragger reads .raw natively).
    // DIA-Quant files are not searched (FragpipeRun.java:1273-1280)
    //
    // Skipped when no file is searched (only DIA-Quant files), as FragPipe (FragpipeRun.java:1335)
    CHECKCENTROID(
        ch_file.toList()
            .filter { files -> checkCentroidFile(files) != null }
            .map { files -> checkCentroidFile(files) }
    )

    FRAGPIPE_CONVERT(
        ch_file
            .filter { meta, _spectra -> meta.data_type != 'DIA-Quant' }
            .combine(CHECKCENTROID.out.log.map { _meta, log -> log })
            .map { meta, spectra, _log -> [meta, spectra] },
        ch_tool_configs
    )
    ch_mzml = FRAGPIPE_CONVERT.out.mzml

    //
    // STEP 4: MSFragger search, Crystal-C and MSBooster as the workflow file sets them: FragPipe
    // routes MSFragger output through both when enabled (FragpipeRun.java:2454-2457), the label
    // check workflows only switch them off.
    //
    ch_msbooster_params = PARSEFRAGPIPEWORKFLOW.out.msbooster_config
        .map { _meta, params -> params }
        .first()

    // MSBooster reads one MSFragger params file chosen from the data types of all files and
    // whether DIA-Umpire / diaTracer run (CmdMSBooster.java:134-149)
    ch_msbooster_fragger_params = ch_file
        .map { meta, _spectra -> meta.data_type }
        .toList()
        .map { data_types -> [data_types] }
        .combine(ch_tool_configs)
        .combine(ch_fragger_params.map { sets -> [sets] })
        .map { data_types, configs, sets ->
            if (!shouldRunTool(configs, 'msbooster')) {
                return []
            }
            def params_type = msboosterFraggerParamsType(data_types, shouldRunTool(configs, 'diaumpire'), shouldRunTool(configs, 'diatracer'))
            sets.find { set -> set[0] == params_type }[1]
        }
        .first()

    FRAGPIPE_SEARCH(
        ch_mzml,
        ch_fasta,
        ch_tool_configs,
        ch_fragger_params,
        ch_msbooster_fragger_params,
        ch_msbooster_params,
        ch_tool_configs.map { configs -> getMsfraggerNumSlices(configs) }
    )

    // FragPipe's LC-MS directory after the search; Percolator names one spectra file per run
    ch_lcms = FRAGPIPE_SEARCH.out.lcms

    //
    // STEP 5: PSM validation via Percolator (parallel per sample)
    // FRAGPIPE_VALIDATE conditionally runs Percolator/PeptideProphet/PTMProphet
    // based on tool_configs parsed from the workflow file.
    //
    FRAGPIPE_VALIDATE(
        FRAGPIPE_SEARCH.out.pin,
        FRAGPIPE_SEARCH.out.pin_rescored,
        FRAGPIPE_SEARCH.out.pepxml,
        FRAGPIPE_SEARCH.out.lcms.map { meta, files -> [meta, files[0]] },   // searched LC-MS file
        FRAGPIPE_SEARCH.out.calibrated_mzml,
        FRAGPIPE_SEARCH.out.uncalibrated_mzml,
        ch_fasta,
        ch_tool_configs
    )

    //
    // STEP 6: Protein inference + FDR filtering (aggregate → per-sample)
    // FRAGPIPE_INFERENCE runs ProteinProphet (aggregate) then Philosopher Filter (per-sample).
    // Filter flags (--sequential --picked --prot 0.01 --razor) come from tool_configs.
    //
    def aggregate_meta = [id: 'tmt_labelcheck']

    FRAGPIPE_INFERENCE(
        FRAGPIPE_VALIDATE.out.pepxml,
        FRAGPIPE_VALIDATE.out.pepxml_group,
        ch_fasta,
        ch_tool_configs,
        aggregate_meta
    )

    //
    // STEP 6b: FreeQuant / LabelQuant + final Report when configured, then Export matched fragments
    // (optional, tab-run.export_matched_fragments): FragPipe runs it right after PhilosopherReport,
    // before IonQuant (FragpipeRun.java:2411, 2463-2471, 2490).
    //
    ch_search_lcms = FRAGPIPE_SEARCH.out.lcms.flatMap { meta, files -> files.collect { f -> [meta, f] } }

    FRAGPIPE_MATCHED_FRAGMENTS_INPUTS(
        ch_file,
        ch_search_lcms,
        FRAGPIPE_CONVERT.out.umpire_mzml
    )

    FRAGPIPE_PHILOSOPHER_QUANT(
        FRAGPIPE_INFERENCE.out.results_dir,
        ch_search_lcms,
        ch_file,
        channel.value([]),  // No TMT annotation: label check runs no TMT-Integrator
        ch_tool_configs,
        FRAGPIPE_MATCHED_FRAGMENTS_INPUTS.out.inputs,
        aggregate_meta
    )
    ch_results_dir = FRAGPIPE_PHILOSOPHER_QUANT.out.final_results_dir

    //
    // STEP 6c: Write sub mzML (tab-run.write_sub_mzml) per LC-MS file after the final report, as in
    // FRAGPIPE_WF (FragpipeRun.java:2388-2404, :2489; CmdWriteSubMzml.java:97-130)
    //
    ch_for_write_sub_mzml = ch_file
        .filter { meta, _spectra -> meta.data_type != 'DIA-Quant' }
        .map { meta, spectra -> [meta.id, meta, spectra] }
        .join(FRAGPIPE_SEARCH.out.mzbin_all.map { meta, mzbin -> [meta.id, mzbin] }, failOnDuplicate: true, remainder: true)
        .combine(ch_tool_configs)
        .filter { _key, _meta, _spectra, _mzbin, configs -> shouldRunTool(configs, 'writesubmzml') }
        .map { key, meta, spectra, mzbin, configs ->
            if (mzbin == null) {
                error("WRITESUBMZML: MSFragger wrote no .mzBIN_all for ${key}, but Write sub mzML is on (write_mzbin_all = 1, FragpipeRun.java:1481)")
            }
            [meta.group, meta, spectra, mzbin, requireToolField(configs, 'writesubmzml', 'prob_threshold')]
        }
        .combine(ch_results_dir.map { meta, dir -> [meta.group, dir.resolve('psm.tsv')] }, by: 0)
        .map { _group, meta, spectra, mzbin, threshold, psm -> [meta, spectra, mzbin, psm, threshold] }

    WRITESUBMZML(ch_for_write_sub_mzml)

    //
    // STEP 7: IonQuant - Aggregate per-sample PSMs into combined peptide table
    // FRAGPIPE_QUANT conditionally runs IonQuant/TMTIntegrator/FreeQuant
    // based on tool_configs. For label check, only IonQuant runs.
    //
    FRAGPIPE_QUANT(
        ch_results_dir,
        ch_search_lcms,
        ch_file,
        channel.value([]),  // No TMT annotation: label check runs no TMT-Integrator
        ch_tool_configs,
        [],  // no AppendFile: label-check workflows run no PTM-Shepherd glycan assignment
        [id: 'tmt_labelcheck_combined']
    )

    //
    // STEP 8a: Combined TMT QC analysis (from per-sample psm.tsv files)
    // Collects all per-sample PSMs and generates a single combined report.
    // Use results_dir (which contains psm.tsv in sample-named subdirs) to avoid
    // Nextflow file staging collision: all per-sample files are named psm.tsv.
    //
    ch_psm_for_analyze = ch_results_dir
        .map { _meta, dir -> dir }
        .collect(sort: true)
        .map { dirs -> [[id: 'tmt_labelcheck'], dirs.sort(false) { dir -> dir.name }.collect { dir -> dir.resolve('psm.tsv') }] }  // each group's Report psm.tsv, groups by name

    TMTLABELCHECKANALYZE(
        ch_psm_for_analyze,
        tmt_type,
        'psm'
    )

    //
    // STEP 8b: Combined TMT QC analysis (from IonQuant combined_modified_peptide.tsv)
    // Provides aggregate view with proper peptide grouping across samples
    //
    TMTLABELCHECKANALYZE_IONQUANT(
        FRAGPIPE_QUANT.out.combined_modified_peptide,
        tmt_type,
        'ionquant'
    )

    emit:
    mzml                 = ch_lcms                                            // channel: [ val(meta), [ path(lcms), path(msfragger side files)... ] ]
    pepxml               = FRAGPIPE_VALIDATE.out.pepxml                       // channel: [ val(meta), path(pepxml) ]
    pin                  = FRAGPIPE_SEARCH.out.pin                                  // channel: [ val(meta), path(pin) ]
    protxml              = FRAGPIPE_INFERENCE.out.protxml                     // channel: [ val(meta), path(protxml) ]
    psms                 = ch_results_dir.map { meta, dir -> [meta, dir.resolve('psm.tsv')] }  // channel: [ val(meta), path(psm.tsv) ] - per sample
    combined_peptide     = FRAGPIPE_QUANT.out.combined_modified_peptide       // channel: [ val(meta), path(combined_modified_peptide.tsv) ]
    html_report          = TMTLABELCHECKANALYZE.out.html_report             // channel: [ val(meta), path(report.html) ] - combined
    md_report            = TMTLABELCHECKANALYZE.out.md_report               // channel: [ val(meta), path(report.md) ] - combined
    summary              = TMTLABELCHECKANALYZE.out.summary                 // channel: [ val(meta), path(labeling_summary.tsv) ] - combined
    per_sample           = TMTLABELCHECKANALYZE.out.per_sample              // channel: [ val(meta), path(per_sample_efficiency.csv) ] - combined
    html_report_combined = TMTLABELCHECKANALYZE_IONQUANT.out.html_report    // channel: combined report
    summary_combined     = TMTLABELCHECKANALYZE_IONQUANT.out.summary        // channel: combined summary
    searched_lcms        = FRAGPIPE_SEARCH.out.lcms.map { meta, files -> [meta, files[0]] }  // channel: [ val(meta), path(file) ] - the LC-MS file MSFragger searched
    database             = ch_fasta                                           // channel: path(fasta) - the run database (decoys, contaminants)
    msfragger_logs       = FRAGPIPE_SEARCH.out.msfragger_logs                 // channel: [ val(meta), path(log) ] - MSFragger console output
    versions             = channel.empty()                                    // versions collected via topic channels
}
