/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_HEADLESS_WF SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs FragPipe in headless mode (all-in-one process) on the LC-MS files as given
    (FragPipe reads .raw natively), after FragPipe's database decoy check.

    Modules:
    - PARSEFRAGPIPEWORKFLOW (workflow decoy tag and tool run flags)
    - FRAGPIPE_DATABASE (CHECKDBCONFIG, plus PHILOSOPHERDATABASE with --add_philosopher_decoy / --add_philosopher_contaminants;
      FragPipe headless hard-fails if the FASTA decoys do not match the workflow decoy tag)
    - FRAGPIPEHEADLESS (all-in-one FragPipe headless execution)

    Execution: a single aggregate process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARSEFRAGPIPEWORKFLOW } from '../../../modules/local/parsefragpipeworkflow/main'
include { parseFragpipeWorkflowTemplates } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'
include { FRAGPIPE_DATABASE   } from '../fragpipe_database/main'
include { FRAGPIPEHEADLESS } from '../../../modules/local/fragpipeheadless/main'

workflow FRAGPIPE_HEADLESS_WF {
    take:
    ch_files                // channel: [ val(meta), path(spectra) ] - sample files
    ch_fasta                // channel: [ val(meta), path(fasta) ] - FASTA to update with decoys and contaminants (empty unless --add_philosopher_decoy / --add_philosopher_contaminants)
    ch_database             // channel: [ val(meta), path(fasta) ] - FASTA searched as given (empty when ch_fasta is used)
    ch_workflow_file        // channel: path(workflow) - FragPipe .workflow file
    ch_manifest_rows        // channel: val([[filename, experiment, bioreplicate, data_type, group], ...]) - value channel, one row per LC-MS file
    ch_annotation_rows      // channel: val([[group, channel, sample_name], ...]) - value channel, TMT annotation, [] for LFQ
    main:

    //
    // STEP 1: FragPipe's database decoy check (FragpipeRun.checkDbConfig) so the run stops here, before the
    // headless process, wherever FragPipe headless would stop.
    // Must run as a separate process (not inside headless script) to avoid Fusion mount corruption.
    // The decoy tag and run flags come from the same .workflow FragPipe headless reads.
    //
    // Data types of the manifest files: FragPipe adjusts the tools to them (TabWorkflow.adjustToolsBasedOnDataTypes)
    ch_data_types = ch_files
        .map { meta, _spectra -> meta.data_type }
        .unique()
        .toSortedList()
    PARSEFRAGPIPEWORKFLOW(
        ch_workflow_file
            .map { wf -> [[id: wf.baseName], wf] }
            .combine(ch_data_types.map { data_types -> [data_types] }),
        parseFragpipeWorkflowTemplates(projectDir)
    )
    ch_tool_configs = PARSEFRAGPIPEWORKFLOW.out.tool_configs_json
        .map { _meta, json_file -> json_file }
        .splitText()
        .collect()
        .map { lines -> lines.join('') }

    FRAGPIPE_DATABASE(ch_fasta, ch_database, ch_tool_configs)
    ch_prepared_db = FRAGPIPE_DATABASE.out.fasta_path

    //
    // STEP 2: Collect all files (strip meta) into a single list for headless mode
    //
    ch_collected_files = ch_files
        .map { _meta, file -> file }
        .collect()

    //
    // STEP 3: Run FragPipe headless
    //
    FRAGPIPEHEADLESS(
        ch_collected_files,
        ch_prepared_db,
        ch_workflow_file.map { wf -> wf.text },  // value channel: the .workflow is read once, FRAGPIPEHEADLESS writes it into its work dir
        ch_manifest_rows,
        ch_annotation_rows
    )

    emit:
    all_results      = FRAGPIPEHEADLESS.out.all_results      // path: results/**
    combined_protein = FRAGPIPEHEADLESS.out.combined_protein  // path: combined_protein.tsv
    combined_peptide = FRAGPIPEHEADLESS.out.combined_peptide  // path: combined_peptide.tsv
    combined_ion     = FRAGPIPEHEADLESS.out.combined_ion      // path: combined_ion.tsv
}
