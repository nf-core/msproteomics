/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_PTM SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PTM analysis and characterization.

    Modules:
    - PTMSHEPHERD (PTM mass shift analysis, localization, and profiling)

    Execution: Aggregate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PTMSHEPHERD } from '../../../modules/local/ptmshepherd/main'
include { shouldRunTool; getToolArgs } from '../fragpipe_utils'

workflow FRAGPIPE_PTM {
    take:
    ch_results_dir   // channel: [ val(meta), path(results_dir) ] - per-sample from PHILOSOPHER_FILTER
    ch_protxml       // channel: [ val(meta), path(protxml) ] - from PROTEINPROPHET
    ch_mzml          // channel: [ val(meta), path(mzml) ] - mzML files
    ch_tool_configs  // channel: val(tool_configs_map) - parsed JSON config
    aggregate_meta   // val(meta) for aggregate outputs
    ch_fasta         // channel: path(fasta) - protein database FASTA

    main:

    //
    // PTMSHEPHERD (aggregate PTM analysis)
    // Analyzes mass shifts, PTM localization, and modification profiles
    //
    ch_all_results = ch_results_dir
        .map { meta, dir -> dir }
        .collect(sort: true)
        .map { items -> [items] }  // Wrap to prevent combine from flattening

    ch_protxml_path = ch_protxml
        .map { meta, protxml -> protxml }

    ch_all_mzml = ch_mzml
        .map { meta, mzml -> mzml }
        .collect(sort: true)
        .map { items -> [items] }  // Wrap to prevent combine from flattening

    ch_ptmshepherd_pre = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'ptmshepherd') }
        .combine(ch_all_results)
        .combine(ch_protxml_path)
        .combine(ch_all_mzml)
        .map { configs, results_wrapped, protxml, mzml_wrapped ->
            def args = getToolArgs(configs, 'ptmshepherd')
            [aggregate_meta, results_wrapped, protxml, mzml_wrapped, args ?: '# PTMShepherd default config']
        }

    // Create config file via collectFile (avoids writing files in workflow scope)
    ch_ptmshepherd_config = ch_ptmshepherd_pre
        .map { meta, _results, _protxml, _mzml, config_content ->
            ["shepherd_${meta.id}.config", config_content]
        }
        .collectFile()
        .map { config_file -> config_file }

    ch_for_ptmshepherd = ch_ptmshepherd_pre
        .map { meta, results, protxml, mzml, _config_content -> [meta, results, protxml, mzml] }
        .combine(ch_ptmshepherd_config)
        .map { meta, results, protxml, mzml, config_file -> [meta, results, protxml, mzml, config_file] }

    PTMSHEPHERD(ch_for_ptmshepherd, ch_fasta)

    emit:
    results_dir      = PTMSHEPHERD.out.results_dir             // channel: [ val(meta), path(results_dir) ]
    global_profile   = PTMSHEPHERD.out.global_profile          // channel: [ val(meta), path(global.profile.tsv) ]
    global_modsummary = PTMSHEPHERD.out.global_modsummary      // channel: [ val(meta), path(global.modsummary.tsv) ]
    diagmine         = PTMSHEPHERD.out.diagmine                // channel: [ val(meta), path(*diagmine.tsv) ]
    localization     = PTMSHEPHERD.out.localization            // channel: [ val(meta), path(*localization.tsv) ]
    glycoprofile     = PTMSHEPHERD.out.glycoprofile            // channel: [ val(meta), path(*glycoprofile.tsv) ]
    versions         = PTMSHEPHERD.out.versions_ptmshepherd
        .mix(PTMSHEPHERD.out.versions_fragpipe)
}
