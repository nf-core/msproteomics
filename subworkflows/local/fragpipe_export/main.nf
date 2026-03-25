/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_EXPORT SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Export and specialized analysis modules.

    Modules:
    - SKYLINE (create Skyline documents)
    - SAINTEXPRESS (AP-MS scoring for protein-protein interactions)
    - FPOP (Fast Photochemical Oxidation of Proteins analysis)
    - METAPROTEOMICS (taxonomic assignment and analysis)

    Execution: Aggregate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SKYLINE        } from '../../../modules/local/skyline/main'
include { SAINTEXPRESS   } from '../../../modules/local/saintexpress/main'
include { FPOP           } from '../../../modules/local/fpop/main'
include { METAPROTEOMICS } from '../../../modules/local/metaproteomics/main'
include { shouldRunTool; getToolArgs } from '../fragpipe_utils'

workflow FRAGPIPE_EXPORT {
    take:
    ch_results_dir   // channel: [ val(meta), path(results_dir) ] - PHILOSOPHER_FILTER output
    ch_mzml          // channel: [ val(meta), path(mzml) ] - mzML files
    ch_speclib       // channel: [ val(meta), path(speclib) ] - spectral library (optional)
    ch_fasta         // channel: path(fasta) - value channel for database
    ch_taxon_files   // channel: tuple path(names.dmp), path(nodes.dmp) - for metaproteomics
    ch_saint_files   // channel: tuple path(inter), path(bait), path(prey) - for SAINT
    ch_tool_configs  // channel: val(tool_configs_map) - parsed JSON config
    aggregate_meta   // val(meta) for aggregate outputs

    main:

    //
    // SKYLINE (create Skyline documents, aggregate)
    //
    // Wrap collected lists in an outer list to prevent .combine() from flattening
    // empty lists ([] + configs would produce [configs] instead of [[], configs]).
    ch_all_results = ch_results_dir
        .map { meta, dir -> dir }
        .collect(sort: true)
        .map { items -> [items] }

    ch_all_mzml = ch_mzml
        .map { meta, mzml -> mzml }
        .collect(sort: true)
        .map { items -> [items] }

    // Speclib is optional — use a sentinel file so combine() always emits a 4-element
    // tuple. Plain [] gets flattened by combine(), collapsing the tuple and causing
    // a MissingMethodException in downstream filter/map closures.
    ch_speclib_path = ch_speclib
        .map { meta, lib -> lib }
        .ifEmpty(file('NO_SPECLIB'))

    ch_for_skyline = ch_all_results
        .combine(ch_all_mzml)
        .combine(ch_speclib_path)
        .combine(ch_tool_configs)
        .filter { results_wrapped, mzml_wrapped, speclib, configs -> shouldRunTool(configs, 'skyline') }
        .map { results_wrapped, mzml_wrapped, speclib, configs ->
            def args = getToolArgs(configs, 'skyline')
            [aggregate_meta, results_wrapped, mzml_wrapped, speclib, args]
        }

    // Skyline path - would be provided by external config in real usage
    SKYLINE(ch_for_skyline, 'skyline')

    //
    // SAINTEXPRESS (AP-MS scoring, aggregate)
    //
    ch_for_saintexpress = ch_saint_files
        .combine(ch_tool_configs)
        .filter { inter, bait, prey, configs -> shouldRunTool(configs, 'saintexpress') }
        .map { inter, bait, prey, configs ->
            def args = getToolArgs(configs, 'saintexpress')
            def mode = args.contains('--mode int') ? 'int' : 'spc'
            [aggregate_meta, inter, bait, prey, mode, args]
        }

    SAINTEXPRESS(ch_for_saintexpress)

    //
    // FPOP (Fast Photochemical Oxidation analysis, aggregate)
    //
    ch_peptide_files = ch_results_dir
        .map { meta, dir -> dir.resolve('combined_modified_peptide.tsv') }
        .collect(sort: true)
        .filter { files -> files.size() > 0 }
        .map { files -> files.first() }

    ch_for_fpop = ch_peptide_files
        .combine(ch_tool_configs)
        .filter { peptide_file, configs -> shouldRunTool(configs, 'fpop') }
        .map { peptide_file, configs ->
            def args = getToolArgs(configs, 'fpop')
            [aggregate_meta, peptide_file, file('NO_FILE'), args]
        }

    FPOP(ch_for_fpop)

    //
    // METAPROTEOMICS (taxonomic analysis, aggregate)
    //
    // Use first() without ifEmpty() - if ch_results_dir is empty, the channel
    // stays empty and combine() correctly produces nothing (skipping metaproteomics).
    ch_project_dir = ch_results_dir
        .map { meta, dir -> dir }
        .first()

    ch_taxon_names = ch_taxon_files.map { names, nodes -> names }
    ch_taxon_nodes = ch_taxon_files.map { names, nodes -> nodes }

    ch_for_metaproteomics = ch_project_dir
        .combine(ch_tool_configs)
        .filter { project, configs -> shouldRunTool(configs, 'metaproteomics') }
        .map { project, configs ->
            def args = getToolArgs(configs, 'metaproteomics')
            [aggregate_meta, project, args]
        }

    METAPROTEOMICS(ch_for_metaproteomics, ch_fasta, ch_taxon_names, ch_taxon_nodes)

    emit:
    skyline_doc      = SKYLINE.out.skyline_document            // channel: [ val(meta), path(fragpipe.sky) ]
    skyline_reports  = SKYLINE.out.reports                     // channel: [ val(meta), path(*.csv) ]
    saint_results    = SAINTEXPRESS.out.results_list           // channel: [ val(meta), path(list.txt) ]
    fpop_results     = FPOP.out.results_dir                    // channel: [ val(meta), path(results_dir) ]
    taxonomy_results = METAPROTEOMICS.out.results_dir          // channel: [ val(meta), path(results_dir) ]
    versions         = SKYLINE.out.versions_skyline
        .mix(SKYLINE.out.versions_fragpipe)
        .mix(SAINTEXPRESS.out.versions_saintexpress)
        .mix(SAINTEXPRESS.out.versions_fragpipe)
        .mix(FPOP.out.versions_fpop)
        .mix(FPOP.out.versions_fragpipe)
        .mix(METAPROTEOMICS.out.versions_metaproteomics)
        .mix(METAPROTEOMICS.out.versions_fragpipe)
}
