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
include { shouldRunTool; getToolArgs; getToolEntry } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

workflow FRAGPIPE_PTM {
    take:
    ch_results_dir   // channel: [ val(meta), path(group_dir) ] - one per FragPipe group (meta.group), psm.tsv inside
    ch_mzml          // channel: [ val(meta), path(mzml) ] - LC-MS files
    ch_tool_configs  // channel: val(tool_configs_json) - parsed workflow config
    aggregate_meta   // val(meta) for aggregate outputs
    ch_fasta         // channel: path(fasta) - the searched database

    main:

    //
    // PTMSHEPHERD (one run over all groups, CmdPtmshepherd.java). It rewrites each group's psm.tsv
    // (glycan assignment), and FragPipe runs MBG, IonQuant and SpecLibGen after it on that file
    // (FragpipeRun.java:2474-2483), so the rewritten group directories replace the input ones.
    //
    ch_groups = ch_results_dir
        .collect(flat: false)
        .map { items -> [items.sort { item -> item[0].group }] }  // wrapped: one list of [meta, dir]

    ch_all_mzml = ch_mzml
        .map { _meta, mzml -> mzml }
        .collect(sort: true)
        .map { items -> [items] }

    // custom annotation file of ptmshepherd.annotation_file, '' for the keywords (PtmshepherdPanel.java:225-232);
    // a declared file that does not exist stops the run; wrapped in a list so combine adds [] as one element
    ch_annotation_file = ch_tool_configs.map { configs ->
        def annotation = getToolEntry(configs, 'ptmshepherd').annotation_file
        [annotation ? file(annotation, checkIfExists: true) : []]
    }

    ch_ptmshepherd_pre = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'ptmshepherd') }
        .combine(ch_annotation_file)
        .combine(ch_groups)
        .combine(ch_all_mzml)
        .map { configs, annotation_file, groups, mzml ->
            // dataset name: the group, "dataset01" for the blank group (PtmshepherdParams.java:86)
            def datasets = groups.collect { item -> [item[1].name, item[0].group == '' ? 'dataset01' : item[0].group] }
            [aggregate_meta, groups.collect { item -> item[1] }, datasets, mzml, getToolArgs(configs, 'ptmshepherd'), annotation_file]
        }

    // property block of shepherd.config (parser text); the module adds the database and dataset lines
    PTMSHEPHERD(
        ch_ptmshepherd_pre
            .map { meta, dirs, datasets, mzml, config_content, annotation_file -> [meta, dirs, datasets, mzml, config_content + '\n', annotation_file] },
        ch_fasta
    )

    // Rewritten group directories back to their group metas (directory name = meta.id)
    ch_ptms_psm_dirs = PTMSHEPHERD.out.psm_dirs
        .combine(ch_groups)
        .flatMap { _meta, dirs, groups ->
            def by_id = groups.collectEntries { item -> [item[0].id, item[0]] }
            [dirs].flatten().collect { dir ->
                if (!by_id.containsKey(dir.name)) {
                    error("PTMSHEPHERD returned group directory ${dir.name}, which is not one of ${by_id.keySet()}")
                }
                [by_id[dir.name], dir]
            }
        }

    // PTM-Shepherd off: the groups go on unchanged
    ch_unchanged = ch_results_dir
        .combine(ch_tool_configs)
        .filter { _meta, _dir, configs -> !shouldRunTool(configs, 'ptmshepherd') }
        .map { meta, dir, _configs -> [meta, dir] }

    emit:
    psm_dirs     = ch_ptms_psm_dirs.mix(ch_unchanged)          // channel: [ val(meta), path(group_dir) ] - per group
    results_dir  = PTMSHEPHERD.out.results_dir                // channel: [ val(meta), path(ptm-shepherd-output) ]
    glyco_masses = PTMSHEPHERD.out.glyco_masses               // channel: [ val(meta), path(glyco_masses_list.txt) ]
    config       = PTMSHEPHERD.out.config                     // channel: [ val(meta), path(shepherd.config) ]
    versions     = PTMSHEPHERD.out.versions_ptmshepherd
        .mix(PTMSHEPHERD.out.versions_fragpipe)
}
