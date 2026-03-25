/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TMT LABEL CHECK SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A nf-core compliant subworkflow for TMT labeling efficiency QC.

    This subworkflow uses a FragPipe workflow file parsed by PARSE_FRAGPIPE_WORKFLOW
    and reuses generic FragPipe subworkflows for maximum code reuse:

    Pipeline Steps:
    1. PARSE_FRAGPIPE_WORKFLOW: Parse workflow file to get tool configs
    2. FRAGPIPE_DATABASE: Add decoys/contaminants to FASTA (once)
    3. FRAGPIPE_CONVERT: RAW → mzML conversion (parallel per sample)
    4. MSFRAGGER: Database search with calibrate_mass=0 (parallel per sample)
       - Single-pass search (no three-stage calibration) for faster QC
       - TMT is set as VARIABLE modification to detect unlabeled peptides
       - Protein N-terminal acetylation included to account for biological blocking
    5. FRAGPIPE_VALIDATE: Percolator PSM-level FDR control (parallel per sample)
    6. FRAGPIPE_INFERENCE: ProteinProphet + Philosopher Filter (aggregate → per-sample)
    7. FRAGPIPE_QUANT: IonQuant aggregate per-sample PSMs (aggregate)
       - Produces combined_modified_peptide.tsv for combined analysis
    8a. TMT_LABELCHECK_ANALYZE (mode=psm): Per-sample labeling efficiency (aggregate)
    8b. TMT_LABELCHECK_ANALYZE (mode=ionquant): Combined labeling efficiency (aggregate)

    Labeling Efficiency Formula:
        Total_Sites = peptide.count('K') + 1  (lysines + N-terminus)
        Labeled_Sites = count(TMT modifications)
        Efficiency = sum(Labeled_Sites) / sum(Total_Sites) * 100
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import config parser
include { PARSE_FRAGPIPE_WORKFLOW } from '../../../modules/local/parse_fragpipe_workflow/main'

// Import modules used directly (not through subworkflows)
include { MSFRAGGER                             } from '../../../modules/local/msfragger/main'
include { TMT_LABELCHECK_ANALYZE                                    } from '../../../modules/local/tmtlabelcheck/main'
include { TMT_LABELCHECK_ANALYZE as TMT_LABELCHECK_ANALYZE_IONQUANT } from '../../../modules/local/tmtlabelcheck/main'

// Import generic FragPipe subworkflows
include { FRAGPIPE_DATABASE  } from '../fragpipe_database/main'
include { FRAGPIPE_CONVERT   } from '../fragpipe_convert/main'
include { FRAGPIPE_VALIDATE  } from '../fragpipe_validate/main'
include { FRAGPIPE_INFERENCE } from '../fragpipe_inference/main'
include { FRAGPIPE_QUANT     } from '../fragpipe_quant/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TMT_LABELCHECK {
    take:
    ch_file            // channel: [ val(sample_meta), path(mzml/raw) ] - sample files with metadata
    ch_database        // channel: [ val(meta), path(fasta) ] - protein database FASTA file
    ch_workflow_file   // channel: [ val(meta), path(workflow_file) ] - FragPipe workflow file for label check
    tmt_type           // val: TMT type string (TMT6, TMT10, TMT16, TMTPRO, etc.)
    ch_decoy_db        // channel: [ val(meta), path(fasta) ] - prebuilt philosopher database (empty if not provided)

    main:

    //
    // STEP 1: Parse workflow config to get tool run flags and args
    //
    PARSE_FRAGPIPE_WORKFLOW(ch_workflow_file)

    // Read tool configs JSON as a string value channel.
    // Strings avoid the Groovy Map-as-named-parameter bug in channel operators.
    // Use splitText to avoid blocking I/O inside map{} on cloud executors.
    ch_tool_configs = PARSE_FRAGPIPE_WORKFLOW.out.tool_configs_json
        .map { _meta, json_file -> json_file }
        .splitText()
        .collect()
        .map { lines -> lines.join('') }
        .first()

    // Extract config files as value channels
    ch_fragger_params = PARSE_FRAGPIPE_WORKFLOW.out.msfragger_config
        .map { _meta, params -> params }
        .first()

    ch_peptideprophet_config = PARSE_FRAGPIPE_WORKFLOW.out.peptideprophet_config
        .map { _meta, config -> config }
        .first()

    ch_tmtintegrator_config = PARSE_FRAGPIPE_WORKFLOW.out.tmtintegrator_config
        .map { _meta, config -> config }
        .first()

    //
    // STEP 2: Prepare database with decoys and contaminants
    // If prebuilt database (ch_decoy_db) is provided, use it directly.
    // Otherwise, run PHILOSOPHER_DATABASE to add decoys and contaminants.
    // concat().first(): prebuilt items come first; if empty, falls through to generated.
    //
    FRAGPIPE_DATABASE(ch_database)
    ch_fasta = ch_decoy_db
        .map { _meta, fasta -> fasta }
        .concat(FRAGPIPE_DATABASE.out.fasta_path)
        .first()

    //
    // STEP 3: Convert RAW to mzML if needed (parallel per sample)
    // FRAGPIPE_CONVERT handles branching by file type (.raw, .mzml, .d)
    //
    FRAGPIPE_CONVERT(ch_file, ch_tool_configs)
    ch_mzml = FRAGPIPE_CONVERT.out.mzml

    //
    // STEP 4: MSFragger search (single pass, calibrate_mass=0)
    // Called directly (not through FRAGPIPE_SEARCH) because label check
    // uses single-pass search for speed, not three-stage calibration.
    // TMT is set as VARIABLE modification via the workflow file.
    //
    MSFRAGGER(ch_mzml, ch_fasta, ch_fragger_params, [], [])

    //
    // STEP 5: PSM validation via Percolator (parallel per sample)
    // FRAGPIPE_VALIDATE conditionally runs Percolator/PeptideProphet/PTMProphet
    // based on tool_configs parsed from the workflow file.
    //
    FRAGPIPE_VALIDATE(
        MSFRAGGER.out.pin,
        MSFRAGGER.out.pepxml,
        ch_mzml,
        ch_fasta,
        ch_tool_configs,
        ch_peptideprophet_config,
        'DDA'
    )

    //
    // STEP 6: Protein inference + FDR filtering (aggregate → per-sample)
    // FRAGPIPE_INFERENCE runs ProteinProphet (aggregate) then Philosopher Filter (per-sample).
    // Filter flags (--sequential --picked --prot 0.01 --razor) come from tool_configs.
    //
    def aggregate_meta = [id: 'tmt_labelcheck']

    FRAGPIPE_INFERENCE(
        FRAGPIPE_VALIDATE.out.pepxml,
        ch_fasta,
        ch_tool_configs,
        aggregate_meta
    )

    //
    // STEP 7: IonQuant - Aggregate per-sample PSMs into combined peptide table
    // FRAGPIPE_QUANT conditionally runs IonQuant/TMTIntegrator/FreeQuant
    // based on tool_configs. For label check, only IonQuant runs.
    //
    FRAGPIPE_QUANT(
        FRAGPIPE_INFERENCE.out.results_dir,
        ch_mzml,
        [],  // No annotation file needed for label check
        ch_fasta,
        ch_tool_configs,
        ch_tmtintegrator_config,
        [id: 'tmt_labelcheck_combined']
    )

    //
    // STEP 8a: Combined TMT QC analysis (from per-sample psm.tsv files)
    // Collects all per-sample PSMs and generates a single combined report.
    // Use results_dir (which contains psm.tsv in sample-named subdirs) to avoid
    // Nextflow file staging collision: all per-sample files are named psm.tsv.
    //
    ch_psm_for_analyze = FRAGPIPE_INFERENCE.out.results_dir
        .map { meta, dir -> dir }
        .collect(sort: true)
        .map { dirs -> [[id: 'tmt_labelcheck'], dirs] }

    TMT_LABELCHECK_ANALYZE(
        ch_psm_for_analyze,
        tmt_type,
        'psm'
    )

    //
    // STEP 8b: Combined TMT QC analysis (from IonQuant combined_modified_peptide.tsv)
    // Provides aggregate view with proper peptide grouping across samples
    //
    TMT_LABELCHECK_ANALYZE_IONQUANT(
        FRAGPIPE_QUANT.out.combined_modified_peptide,
        tmt_type,
        'ionquant'
    )

    emit:
    mzml                 = ch_mzml                                            // channel: [ val(meta), path(mzml) ]
    pepxml               = FRAGPIPE_VALIDATE.out.pepxml                       // channel: [ val(meta), path(pepxml) ]
    pin                  = MSFRAGGER.out.pin                                  // channel: [ val(meta), path(pin) ]
    protxml              = FRAGPIPE_INFERENCE.out.protxml                     // channel: [ val(meta), path(protxml) ]
    psms                 = FRAGPIPE_INFERENCE.out.psms                        // channel: [ val(meta), path(psm.tsv) ] - per sample
    combined_peptide     = FRAGPIPE_QUANT.out.combined_modified_peptide       // channel: [ val(meta), path(combined_modified_peptide.tsv) ]
    html_report          = TMT_LABELCHECK_ANALYZE.out.html_report             // channel: [ val(meta), path(report.html) ] - combined
    md_report            = TMT_LABELCHECK_ANALYZE.out.md_report               // channel: [ val(meta), path(report.md) ] - combined
    summary              = TMT_LABELCHECK_ANALYZE.out.summary                 // channel: [ val(meta), path(labeling_summary.tsv) ] - combined
    per_sample           = TMT_LABELCHECK_ANALYZE.out.per_sample              // channel: [ val(meta), path(per_sample_efficiency.csv) ] - combined
    html_report_combined = TMT_LABELCHECK_ANALYZE_IONQUANT.out.html_report    // channel: combined report
    summary_combined     = TMT_LABELCHECK_ANALYZE_IONQUANT.out.summary        // channel: combined summary
    versions             = Channel.empty()                                    // versions collected via topic channels
}
