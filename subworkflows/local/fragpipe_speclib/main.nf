/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_SPECLIB SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Spectral library generation for DIA analysis.

    Modules:
    - SPECLIBGEN (generates spectral libraries using EasyPQP/FragPipe-SpecLib)

    Execution: Aggregate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPECLIBGEN } from '../../../modules/local/speclibgen/main'
include { shouldRunTool; getToolArgs } from '../fragpipe_utils'

workflow FRAGPIPE_SPECLIB {
    take:
    ch_psms          // channel: [ val(meta), path(psm.tsv) ] - per-sample PSM files
    ch_peptides      // channel: [ val(meta), path(peptide.tsv) ] - per-sample peptide files
    ch_mzml          // channel: [ val(meta), path(mzml) ] - mzML files
    ch_fasta         // channel: path(fasta) - value channel for database
    ch_tool_configs  // channel: val(tool_configs_map) - parsed JSON config
    aggregate_meta   // val(meta) for aggregate outputs

    main:

    //
    // SPECLIBGEN (spectral library generation, aggregate)
    // Generates spectral libraries for DIA analysis.
    // Collects per-sample psm.tsv and peptide.tsv files (avoiding Fusion
    // directory staging issues where directories may be staged as file symlinks).
    //
    ch_all_psms = ch_psms
        .map { meta, psm -> psm }
        .collect(sort: true)
        .map { items -> [items] }  // Wrap to prevent combine from flattening

    ch_all_peptides = ch_peptides
        .map { meta, pep -> pep }
        .collect(sort: true)
        .map { items -> [items] }  // Wrap to prevent combine from flattening

    ch_all_mzml = ch_mzml
        .map { meta, mzml -> mzml }
        .collect(sort: true)
        .map { items -> [items] }  // Wrap to prevent combine from flattening

    // Combine collected files with tool_configs
    ch_for_speclibgen = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'speclibgen') }
        .combine(ch_all_psms)
        .combine(ch_all_peptides)
        .combine(ch_all_mzml)
        .map { configs, psms_wrapped, peps_wrapped, mzml_wrapped ->
            def args = getToolArgs(configs, 'speclibgen')
            [aggregate_meta, psms_wrapped, peps_wrapped, mzml_wrapped, args]
        }

    SPECLIBGEN(ch_for_speclibgen, ch_fasta)

    emit:
    library_tsv      = SPECLIBGEN.out.library_tsv              // channel: [ val(meta), path(library.tsv) ]
    library_speclib  = SPECLIBGEN.out.library_speclib          // channel: [ val(meta), path(library.speclib) ]
    library_parquet  = SPECLIBGEN.out.library_parquet          // channel: [ val(meta), path(library.parquet) ]
    results_dir      = SPECLIBGEN.out.results_dir              // channel: [ val(meta), path(results_dir) ]
    versions         = SPECLIBGEN.out.versions_speclibgen
        .mix(SPECLIBGEN.out.versions_fragpipe)
}
