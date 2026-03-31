/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_GLYCO SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Glycoproteomics analysis modules.

    Modules:
    - MBG (Mass-Based Glycoproteomics matching)
    - OPAIR (O-glycoproteomics analysis)

    Execution: Per-sample or aggregate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MBG   } from '../../../modules/local/mbg/main'
include { OPAIR } from '../../../modules/local/opair/main'
include { shouldRunTool; getToolArgs } from '../fragpipe_utils'

workflow FRAGPIPE_GLYCO {
    take:
    ch_psm_files     // channel: [ val(meta), path(psm.tsv) ] - from PHILOSOPHER_FILTER
    ch_mzml          // channel: [ val(meta), path(mzml) ] - mzML files
    ch_manifest      // channel: path(manifest.fp-manifest) - FragPipe manifest
    ch_glycan_dbs    // channel: tuple path(residue_db), path(glycan_mod_db), path(oglycan_db)
    ch_tool_configs  // channel: val(tool_configs_map) - parsed JSON config

    main:

    //
    // MBG (Mass-Based Glycoproteomics, per-sample)
    // Performs glycan matching on PSM data
    //
    ch_psm_keyed = ch_psm_files
        .map { meta, psm -> [meta.id, meta, psm] }

    ch_for_mbg = ch_psm_keyed
        .combine(ch_manifest)
        .combine(ch_tool_configs)
        .filter { key, meta, psm, manifest, configs -> shouldRunTool(configs, 'mbg') }
        .map { key, meta, psm, manifest, configs ->
            def args = getToolArgs(configs, 'mbg')
            [meta, psm, manifest, args]
        }

    // Extract glycan database files
    ch_residue_db = ch_glycan_dbs.map { residue, glycan_mod, oglycan -> residue }
    ch_glycan_mod_db = ch_glycan_dbs.map { residue, glycan_mod, oglycan -> glycan_mod }
    ch_oglycan_db = ch_glycan_dbs.map { residue, glycan_mod, oglycan -> oglycan }

    MBG(ch_for_mbg, ch_residue_db, ch_glycan_mod_db)

    //
    // O-Pair (O-glycoproteomics, per-sample)
    // Performs O-glycan site localization
    //
    ch_mzml_keyed = ch_mzml
        .map { meta, mzml -> [meta.id, mzml] }

    ch_for_opair = ch_psm_keyed
        .join(ch_mzml_keyed, by: 0, failOnMismatch: false, failOnDuplicate: false)
        .combine(ch_tool_configs)
        .filter { key, meta, psm, mzml, configs -> shouldRunTool(configs, 'opair') }
        .map { key, meta, psm, mzml, configs ->
            def args = getToolArgs(configs, 'opair')
            [meta, psm, mzml, args]
        }

    OPAIR(ch_for_opair, ch_oglycan_db)

    emit:
    mbg_results      = MBG.out.results                         // channel: [ val(meta), path(*_mbg*.tsv) ]
    mbg_glycan       = MBG.out.glycan_results                  // channel: [ val(meta), path(*_glycan*.tsv) ]
    opair_results    = OPAIR.out.results                       // channel: [ val(meta), path(*_opair_results.tsv) ]
    opair_glycoforms = OPAIR.out.glycoforms                    // channel: [ val(meta), path(*_opair_glycoforms.tsv) ]
    mbg_dir          = MBG.out.results_dir                     // channel: [ val(meta), path(results_dir) ]
    opair_dir        = OPAIR.out.results_dir                   // channel: [ val(meta), path(results_dir) ]
    versions         = MBG.out.versions_mbg
        .mix(MBG.out.versions_fragpipe)
        .mix(OPAIR.out.versions_opair)
        .mix(OPAIR.out.versions_fragpipe)
}
