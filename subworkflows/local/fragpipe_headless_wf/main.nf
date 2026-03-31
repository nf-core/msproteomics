/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_HEADLESS_WF SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs FragPipe in headless mode (all-in-one process).
    Pre-converts .raw files to .mzML via ThermoRawFileParser, prepares the
    database with decoys via Philosopher, then passes all spectra files to
    FRAGPIPE_HEADLESS as a single batch.

    Modules:
    - THERMORAWFILEPARSER (RAW to mzML)
    - PHILOSOPHER_DATABASE (decoy generation — required because FragPipe headless
      hard-fails if FASTA has no decoys, unlike GUI mode which prompts the user)
    - FRAGPIPE_HEADLESS (all-in-one FragPipe headless execution)

    Execution: Conversion is per-sample parallel; headless is a single aggregate process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { THERMORAWFILEPARSER } from '../../../modules/nf-core/thermorawfileparser/main'
include { PHILOSOPHER_DATABASE } from '../../../modules/local/philosopher/database/main'
include { FRAGPIPE_HEADLESS   } from '../../../modules/local/fragpipe_headless/main'

workflow FRAGPIPE_HEADLESS_WF {
    take:
    ch_files                // channel: [ val(meta), path(spectra) ] - sample files
    ch_database             // channel: path(fasta) - FASTA database
    ch_workflow_file        // channel: path(workflow) - FragPipe .workflow file
    ch_manifest_content     // channel: val(manifest_string) - TSV manifest content
    ch_annotation_content   // channel: val(annotation_string) - TMT annotation (experiment\tchannel\tsample_name), empty for LFQ
    ch_file_experiment_map  // channel: val(map_string) - file-to-experiment mapping (filename\texperiment), empty for LFQ
    main:

    //
    // Branch input files by extension: .raw files need conversion, others pass through
    //
    ch_branched = ch_files
        .branch {
            _meta, file ->
            raw:   file.name.toLowerCase().endsWith('.raw')
            ready: true  // .mzML, .d, .mzXML pass through
        }

    //
    // STEP 1: Convert RAW files to mzML using ThermoRawFileParser
    //
    THERMORAWFILEPARSER(ch_branched.raw)

    //
    // STEP 2: Combine converted mzML with already-ready files
    //
    ch_all_spectra = THERMORAWFILEPARSER.out.spectra
        .mix(ch_branched.ready)

    //
    // STEP 3: Prepare database with decoys via Philosopher
    // FragPipe headless hard-fails if no decoys are present (FragpipeRun.java:2724-2727).
    // PHILOSOPHER_DATABASE handles both cases: adds decoys if missing, annotates if present.
    // Must run as a separate process (not inside headless script) to avoid Fusion mount corruption.
    //
    ch_db_with_meta = ch_database.map { db -> [[id: db.baseName], db] }
    PHILOSOPHER_DATABASE(ch_db_with_meta)
    ch_prepared_db = PHILOSOPHER_DATABASE.out.fasta.map { _meta, db -> db }

    //
    // STEP 4: Collect all files (strip meta) into a single list for headless mode
    //
    ch_collected_files = ch_all_spectra
        .map { _meta, file -> file }
        .collect()

    //
    // STEP 5: Run FragPipe headless
    //
    FRAGPIPE_HEADLESS(
        ch_collected_files,
        ch_prepared_db,
        ch_workflow_file,
        ch_manifest_content,
        ch_annotation_content,
        ch_file_experiment_map
    )

    emit:
    all_results      = FRAGPIPE_HEADLESS.out.all_results      // path: results/**
    combined_protein = FRAGPIPE_HEADLESS.out.combined_protein  // path: combined_protein.tsv
    combined_peptide = FRAGPIPE_HEADLESS.out.combined_peptide  // path: combined_peptide.tsv
    combined_ion     = FRAGPIPE_HEADLESS.out.combined_ion      // path: combined_ion.tsv
    versions         = THERMORAWFILEPARSER.out.versions_thermorawfileparser
        .mix(PHILOSOPHER_DATABASE.out.versions_philosopher)
        .mix(FRAGPIPE_HEADLESS.out.versions)
}
