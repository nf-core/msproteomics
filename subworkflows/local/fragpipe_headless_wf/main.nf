/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_HEADLESS_WF SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs FragPipe in headless mode (all-in-one process).
    Pre-converts .raw files to .mzML via ThermoRawFileParser, then passes
    all spectra files to FRAGPIPE_HEADLESS as a single batch.

    Modules:
    - THERMORAWFILEPARSER (RAW to mzML)
    - FRAGPIPE_HEADLESS (all-in-one FragPipe headless execution)

    Execution: Conversion is per-sample parallel; headless is a single aggregate process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { THERMORAWFILEPARSER } from '../../../modules/nf-core/thermorawfileparser/main'
include { FRAGPIPE_HEADLESS   } from '../../../modules/local/fragpipe_headless/main'

workflow FRAGPIPE_HEADLESS_WF {
    take:
    ch_files             // channel: [ val(meta), path(spectra) ] - sample files
    ch_database          // channel: path(fasta) - FASTA database
    ch_workflow_file     // channel: path(workflow) - FragPipe .workflow file
    ch_manifest_content  // channel: val(manifest_string) - TSV manifest content

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
    // STEP 3: Collect all files (strip meta) into a single list for headless mode
    //
    ch_collected_files = ch_all_spectra
        .map { _meta, file -> file }
        .collect()

    //
    // STEP 4: Run FragPipe headless
    //
    FRAGPIPE_HEADLESS(
        ch_collected_files,
        ch_database,
        ch_workflow_file,
        ch_manifest_content
    )

    emit:
    all_results      = FRAGPIPE_HEADLESS.out.all_results      // path: results/**
    combined_protein = FRAGPIPE_HEADLESS.out.combined_protein  // path: combined_protein.tsv
    combined_peptide = FRAGPIPE_HEADLESS.out.combined_peptide  // path: combined_peptide.tsv
    combined_ion     = FRAGPIPE_HEADLESS.out.combined_ion      // path: combined_ion.tsv
    versions         = THERMORAWFILEPARSER.out.versions_thermorawfileparser
        .mix(FRAGPIPE_HEADLESS.out.versions)
}
