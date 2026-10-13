include { peptideProphetCommands; peptideProphetDirLines; shellQuote; fragpipeToolVersions } from '../../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process PEPTIDEPROPHET {
    tag "${meta.id}"
    label 'process_medium'

    input:
    // pepxml_files: one MSFragger pepXML (per-file mode), or every pepXML of one FragPipe group in
    //   FragPipe order (combine mode, CmdPeptideProphet.java:302-336)
    // spectra: the searched LC-MS files of those pepXMLs plus MSFragger's _calibrated/_uncalibrated
    //   mzML; RewritePepxml falls back to the files that exist (RewritePepxml.java:31-62)
    // lcms_names: per LC-MS file the path name FragPipe passes to RewritePepxml
    //   (CmdPeptideProphet.java:439-450), distinct, in FragPipe order
    // peptideprophet_cli: cmd-opts with FragPipe's enzyme flags (CmdPeptideProphet.java:375-404)
    // combine: one --combine run for the group (CmdPeptideProphet.java:190-202)
    tuple val(meta), path(pepxml_files), path(spectra), val(lcms_names), val(peptideprophet_cli), val(combine)
    path(fasta)
    val decoy_tag    // decoy protein prefix from workflow database.decoy-tag (FragpipeRun.java:1442)

    output:
    tuple val(meta), path("${prefix}")                   , emit: results_dir
    tuple val(meta), path("${prefix}/interact*.pep.xml") , emit: pepxml
    tuple val("${task.process}"), val('philosopher'), val(fragpipeToolVersions().philosopher), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = peptideProphetCommands(meta, pepxml_files, lcms_names, peptideprophet_cli, args, combine, fasta, decoy_tag, task)
    def gdir = shellQuote(cmds.group_dir)
    def lcms_args = cmds.lcms.collect { name -> "\"\$PWD\"/${shellQuote(name)}" }.join(' ')
    // the group dir holds the pepXMLs like FragPipe's (the task dir itself for the blank group)
    def inputs = [pepxml_files].flatten().collect { f -> shellQuote(f.name) }
    def link_inputs = cmds.group_dir == '.' ? '' : inputs.collect { name -> "ln -s \"\$PWD\"/${name} ${gdir}/${name}" }.join('\n    ')
    def unlink_inputs = cmds.group_dir == '.' ? '' : "rm ${inputs.collect { name -> "${gdir}/${name}" }.join(' ')}"
    def move_interact = cmds.group_dir == prefix ? '' : "mv ${shellQuote(cmds.interact_path)} ${shellQuote(prefix)}/"
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix} ${gdir}
    ${link_inputs}
    [ -x ${shellQuote(cmds.philosopher)} ] || { echo "ERROR: ${cmds.philosopher} not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    [ -d ${shellQuote(cmds.fragpipe_lib)} ] || { echo "ERROR: ${cmds.fragpipe_lib} not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # Philosopher workspace with its own temp dir (CmdPeptideProphet.java:248-262,
    # CmdPhilosopherWorkspaceCleanInit.java:66-74); per-file mode runs in <group dir>/fragpipe-<pepxml>-temp
    run_dir=${shellQuote(cmds.run_dir)}
    phi_temp="\$PWD/phi_temp/\$(cat /proc/sys/kernel/random/uuid)"
    mkdir -p "\$run_dir" "\$phi_temp"
    # InteractParser's 256-byte <pepXML dir>/<base name> buffer: a short /tmp link for a long work dir
${peptideProphetDirLines(cmds.group_dir, cmds.run_names)}
    (cd "\$PP_GROUP_DIR"/${shellQuote(cmds.run_sub)} && ${shellQuote(cmds.philosopher)} workspace --init --nocheck --temp "\$phi_temp")

    # PeptideProphet (CmdPeptideProphet.java:264-275 per file, :319-336 combined); WEBSERVER_ROOT as in setupEnv (:365-368)
    pp_cmd=(${cmds.peptideprophet.collect { token -> token.startsWith('{CWD}/') ? "\"\$PWD\"/" + shellQuote(token.substring('{CWD}/'.length())) : shellQuote(token) }.join(' ')})
    printf '%s\\n' "\${pp_cmd[*]}" > ${prefix}/peptideprophet.log
    (cd "\$PP_GROUP_DIR"/${shellQuote(cmds.run_sub)} && WEBSERVER_ROOT=fake-WEBSERVER_ROOT-value "\${pp_cmd[@]}") 2>&1 | tee -a ${prefix}/peptideprophet.log

    # RewritePepxml always follows PeptideProphet (CmdPeptideProphet.java:340-355, :428-451)
    rewrite_cmd=(${shellQuote(cmds.rewrite[0])} -Duser.home="\$PWD" ${cmds.rewrite.drop(1).collect { token -> shellQuote(token) }.join(' ')} "\$PWD"/${shellQuote(cmds.interact_path)} ${lcms_args})
    printf '%s\\n' "\${rewrite_cmd[*]}" >> ${prefix}/peptideprophet.log
    "\${rewrite_cmd[@]}" 2>&1 | tee -a ${prefix}/peptideprophet.log

    ${unlink_inputs}
    ${move_interact}
    rm -rf ${combine ? '"\$run_dir"/.meta' : '"\$run_dir"'} phi_temp
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = peptideProphetCommands(meta, pepxml_files, lcms_names, peptideprophet_cli, task.ext.args ?: '', combine, fasta, decoy_tag, task)
    def pp_line = "cd ${cmds.run_dir}; " + cmds.peptideprophet.join(' ')
    def rewrite_line = (cmds.rewrite + ["{CWD}/${cmds.interact_path}"] + cmds.lcms.collect { name -> "{CWD}/${name}" }).join(' ')
    """
    mkdir -p ${prefix}
    touch ${prefix}/${cmds.interact}
    printf '%s\\n' ${shellQuote(pp_line)} ${shellQuote(rewrite_line)} > ${prefix}/peptideprophet.log
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
