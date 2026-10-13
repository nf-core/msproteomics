/*
 * Test harness: FRAGPIPE_GLYCO followed by FRAGPIPE_QUANT, wired as FRAGPIPE_WF wires them
 * (subworkflows/local/fragpipe/main.nf, STEP 7a-7b), for the FragPipe golden-run comparison.
 */
include { FRAGPIPE_GLYCO } from '../main'
include { FRAGPIPE_QUANT } from '../../fragpipe_quant/main'

workflow GLYCO_CHAIN {
    take:
    ch_results_dir   // channel: [ val(meta), path(group_dir) ] - FragPipe Report group dirs
    ch_lcms          // channel: [ val(meta), path(spectra) ] - LC-MS files, file-table order
    ch_fasta         // channel: path(fasta)
    ch_tool_configs  // channel: val(tool_configs_json)

    main:
    FRAGPIPE_GLYCO(ch_results_dir, ch_lcms, ch_lcms, channel.empty(), channel.empty(), ch_fasta, ch_tool_configs, [id: 'fragpipe'])
    // FRAGPIPE_QUANT take: results, searched spectra, LC-MS file list, TMT annotation, tool configs,
    // AppendFile masses, aggregate meta
    FRAGPIPE_QUANT(FRAGPIPE_GLYCO.out.psm_dirs, ch_lcms, ch_lcms, channel.value([]), ch_tool_configs, FRAGPIPE_GLYCO.out.modlist_append, [id: 'fragpipe'])

    emit:
    psm_dirs     = FRAGPIPE_GLYCO.out.psm_dirs
    ptm_results  = FRAGPIPE_GLYCO.out.ptm_results
    config       = FRAGPIPE_GLYCO.out.shepherd_config
    mbg_manifest = FRAGPIPE_GLYCO.out.mbg_manifest
    mbg_log      = FRAGPIPE_GLYCO.out.mbg_log
    ionquant_dir = FRAGPIPE_QUANT.out.ionquant_dir
}
