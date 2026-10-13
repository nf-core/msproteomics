/*
 * Test harness: FRAGPIPE_PHILOSOPHER_QUANT followed by FRAGPIPE_QUANT, wired as FRAGPIPE_WF wires
 * them without the glyco chain (subworkflows/local/fragpipe/main.nf, STEP 6b and 7b).
 */
include { FRAGPIPE_PHILOSOPHER_QUANT } from '../main'
include { FRAGPIPE_QUANT             } from '../../fragpipe_quant/main'

workflow QUANT_CHAIN {
    take:
    ch_results_dir                // channel: [ val(meta), path(results_dir) ] - per-group PHILOSOPHERFILTER output
    ch_mzml                       // channel: [ val(meta), path(mzml) ] - searched spectra
    ch_lcms_files                 // channel: [ val(meta), path(spectra) ] - FragPipe's LC-MS file list
    ch_annotation                 // channel: path(annotation_file) - [] without TMT-Integrator
    ch_tool_configs               // channel: val(tool_configs_json)
    ch_emf_inputs                 // channel: [ path(manifest), [spectra], [umpire_spectra] ]
    ch_glyco_masses               // channel: path(glyco_masses_list.txt) - [] without AppendFile
    aggregate_meta                // val(meta) for aggregate outputs

    main:
    FRAGPIPE_PHILOSOPHER_QUANT(ch_results_dir, ch_mzml, ch_lcms_files, ch_annotation, ch_tool_configs, ch_emf_inputs, aggregate_meta)
    FRAGPIPE_QUANT(FRAGPIPE_PHILOSOPHER_QUANT.out.final_results_dir, ch_mzml, ch_lcms_files, ch_annotation, ch_tool_configs,
        ch_glyco_masses, aggregate_meta)
    ch_report_results_dir = FRAGPIPE_PHILOSOPHER_QUANT.out.final_results_dir.filter { meta, _dir -> meta.report_after_quant }

    emit:
    final_results_dir     = FRAGPIPE_PHILOSOPHER_QUANT.out.final_results_dir
    report_results_dir    = ch_report_results_dir
    report_psms           = ch_report_results_dir.map { meta, dir -> [meta, dir.resolve('psm.tsv')] }
    experiment_annotation = FRAGPIPE_QUANT.out.experiment_annotation.mix(FRAGPIPE_PHILOSOPHER_QUANT.out.experiment_annotation)
    combined_protein      = FRAGPIPE_QUANT.out.combined_protein
    ionquant_dir          = FRAGPIPE_QUANT.out.ionquant_dir
    tmt_abundance         = FRAGPIPE_QUANT.out.tmt_abundance
    tmt_dir               = FRAGPIPE_QUANT.out.tmt_dir
}
