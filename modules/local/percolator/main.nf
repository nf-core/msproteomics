include { percolatorCommands; shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'
include { shellArgs } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { fragpipeLcmsName } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils'

process PERCOLATOR {
    tag "${meta.id}"
    label 'process_medium'

    input:
    // pin: MSFragger PIN (<stem>.pin); PercolatorOutputToPepXML always reads it (CmdPercolator.java:292)
    // edited_pin: MSBooster <stem>_edited.pin, or [] when MSBooster did not run; Percolator reads it
    //   instead of the MSFragger PIN when MSBooster ran (CmdPercolator.java:209-213)
    // pepxml_files: MSFragger pepXML(s) of the file (<stem>.pepXML, or <stem>_rank<N>.pepXML for
    //   DDA+/DIA; Crystal-C output carries the _c suffix)
    // spectra: the searched LC-MS file; msfragger_mzml: MSFragger's <stem>_calibrated.mzML /
    //   <stem>_uncalibrated.mzML ([] when none). FragPipe names <stem>_calibrated.mzML when
    //   has_calibrated_mzml, <stem>_uncalibrated.mzML for .raw/.d, else the searched file
    //   (CmdPercolator.java:219-224); PercolatorOutputToPepXML falls back to the files that
    //   exist next to it (PercolatorOutputToPepXML.java:266-295), so all of them are staged.
    // percolator_cli: workflow percolator.cmd-opts
    // min_prob: PercolatorOutputToPepXML min probability as FragPipe prints it (Double.toString)
    // keep_target_tsv / keep_decoy_tsv: false deletes the TSV like CmdPercolator.java:233-247
    // has_crystalc: Crystal-C ran; the PIN/TSV basename drops the _c suffix (CmdPercolator.java:145-150)
    // write_weights: Write sub mzML is on (tabRun.isWriteSubMzml(), FragpipeRun.java:1642); Percolator saves
    //   <stripped>_percolator.weights for the second pass (CmdPercolator.java:180-183)
    // init_weights: second pass, the first pass's <stem>_percolator.weights of the _sub.mzML input, [] when none
    //   (CmdPercolator.java:185-201); staged in weights/
    tuple val(meta), path(pin), path(edited_pin), path(pepxml_files), path(spectra), path(msfragger_mzml), val(has_calibrated_mzml), val(percolator_cli), val(min_prob), val(keep_target_tsv), val(keep_decoy_tsv), val(has_crystalc), val(write_weights), path(init_weights, stageAs: 'weights/*')
    val decoy_tag    // decoy protein prefix from workflow database.decoy-tag (FragpipeRun.java:1442)
    // meta.data_type: FragPipe data type of the searched file. PercolatorOutputToPepXML gets
    // "DDA" only for DDA files, "DIA" for every other type (CmdPercolator.java:227, :297).

    output:
    tuple val(meta), path("${prefix}")                     , emit: results_dir
    tuple val(meta), path("${prefix}/interact-*.pep.xml")  , emit: pepxml
    tuple val(meta), path("${prefix}/*_target_psms.tsv")   , emit: target_psms, optional: true
    tuple val(meta), path("${prefix}/*_decoy_psms.tsv")    , emit: decoy_psms, optional: true
    // in the task dir, not in ${prefix}: a file inside a published directory output is not published again on its own,
    // and the weights are published to fragpipe/<group>/ (conf/modules.config)
    tuple val(meta), path("*_percolator.weights")          , emit: weights, optional: true
    tuple val("${task.process}"), val('percolator'), val(fragpipeToolVersions().percolator), emit: versions_percolator, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = percolatorCommands(meta, pin, edited_pin, pepxml_files, percolator_cli, args, min_prob, has_crystalc, decoy_tag, write_weights, init_weights, task)
    def gdir = shellQuote(cmds.group_dir)
    def target_tsv = "${cmds.stripped}_percolator_target_psms.tsv"
    def decoy_tsv = "${cmds.stripped}_percolator_decoy_psms.tsv"
    // the group dir holds the PIN(s) and pepXML(s) like FragPipe's (the task dir itself for the blank group)
    def inputs = ([pin] + [edited_pin].flatten().findAll { f -> f } + [pepxml_files].flatten()).collect { f -> shellQuote(f.name) }
    def link_inputs = cmds.group_dir == '.' ? '' : inputs.collect { name -> "ln -s \"\$PWD\"/${name} ${gdir}/${name}" }.join('\n    ')
    def unlink_inputs = cmds.group_dir == '.' ? '' : "rm ${inputs.collect { name -> "${gdir}/${name}" }.join(' ')}"
    def weights = "${cmds.stripped}_percolator.weights"
    def outputs = cmds.interact.collect { name -> shellQuote(name) } + (keep_target_tsv ? [shellQuote(target_tsv)] : []) + (keep_decoy_tsv ? [shellQuote(decoy_tsv)] : [])
    def move_outputs = cmds.group_dir == prefix ? '' : "mv ${outputs.collect { name -> "${gdir}/${name}" }.join(' ')} ${shellQuote(prefix)}/"
    def drop_tsv = (keep_target_tsv ? [] : [target_tsv]) + (keep_decoy_tsv ? [] : [decoy_tsv])
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix} ${gdir}
    ${link_inputs}

    for part in ${shellQuote(cmds.percolator[0])} ${shellQuote(cmds.fragpipe_lib)}; do
        [ -e "\$part" ] || { echo "ERROR: FragPipe component \$part not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    done
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # Percolator, argv in CmdPercolator.java:152-213 order, run in the group dir (CmdPercolator.java:214)
    perc_cmd=(${shellArgs(cmds.percolator)})
    printf '%s\\n' "\${perc_cmd[*]}" > ${prefix}/percolator.log
    (cd ${gdir} && "\${perc_cmd[@]}") 2>&1 | tee -a ${prefix}/percolator.log

    # PercolatorOutputToPepXML, argv in CmdPercolator.java:280-300 order, run in the group dir (CmdPercolator.java:228)
    # -Duser.home: the container user has no passwd entry, so the JVM home would be "?" in the task dir
    conv_cmd=(${shellQuote(cmds.converter[0])} -Duser.home="\$PWD" ${cmds.converter.drop(1).collect { token -> shellQuote(token) }.join(' ')} "\$PWD"/${shellQuote(fragpipeLcmsName(spectra.name, has_calibrated_mzml))})
    printf '%s\\n' "\${conv_cmd[*]}" >> ${prefix}/percolator.log
    (cd ${gdir} && "\${conv_cmd[@]}") 2>&1 | tee -a ${prefix}/percolator.log

    ${unlink_inputs}
    ${drop_tsv ? "rm ${drop_tsv.collect { name -> "${gdir}/${shellQuote(name)}" }.join(' ')}" : ''}
    ${move_outputs}
    ${write_weights && cmds.group_dir != '.' ? "mv ${gdir}/${shellQuote(weights)} ." : ''}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = percolatorCommands(meta, pin, edited_pin, pepxml_files, percolator_cli, task.ext.args ?: '', min_prob, has_crystalc, decoy_tag, write_weights, init_weights, task)
    """
    mkdir -p ${prefix}
    touch ${cmds.interact.collect { name -> "${prefix}/${name}" }.join(' ')}
    ${keep_target_tsv ? "touch ${prefix}/${cmds.stripped}_percolator_target_psms.tsv" : ''}
    ${keep_decoy_tsv ? "touch ${prefix}/${cmds.stripped}_percolator_decoy_psms.tsv" : ''}
    ${write_weights ? "touch ${cmds.stripped}_percolator.weights" : ''}
    printf '%s\\n' ${shellQuote(cmds.percolator.join(' '))} ${shellQuote(cmds.converter.join(' ') + " {CWD}/${fragpipeLcmsName(spectra.name, has_calibrated_mzml)}")} ${shellQuote("cwd ${cmds.group_dir}")} > ${prefix}/percolator.log
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
