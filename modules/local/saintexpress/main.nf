include { saintexpressCommand; shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process SAINTEXPRESS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    // FragPipe's Downstream SAINTexpress run for one mode (TabDownstream.java:114-160, CmdSaintExpress.java:56-96):
    // saint_input (staged files): inter.dat, bait.dat and prey.dat of REPRINTTOSAINTINPUT (CmdSaintExpress.reprintTsvToSaintInput);
    //   a mode without a reprint.<mode>.tsv never reaches this process (REPRINTTOSAINTINPUT reports the skip)
    // mode: 'spc' (SAINTexpress-spc) or 'int' (SAINTexpress-int) (TabDownstream.java:70-83)
    // saint_cli: "-R<max-replicates> -L<virtual-controls>"
    // cmd_opts: saintexpress.cmd-opts, passed as ONE argument when not empty (CmdSaintExpress.java:88-90)
    tuple val(meta), path(saint_input, stageAs: 'saint_input/*'), val(mode), val(saint_cli), val(cmd_opts)

    output:
    tuple val(meta), path("saintexpress-${mode}")            , emit: results_dir
    tuple val(meta), path("saintexpress-${mode}/list.txt")   , emit: results_list
    tuple val("${task.process}"), val('saintexpress'), val(fragpipeToolVersions().saintexpress), emit: versions_saintexpress, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cmds = saintexpressCommand(saint_input, mode, saint_cli, cmd_opts, args, task)
    """
    export HOME=\$(pwd)
    [ -x ${shellQuote(cmds.argv[0])} ] || { echo "ERROR: ${cmds.argv[0]} not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # SAINTexpress reads inter/bait/prey.dat from, and writes list.txt to, saintexpress-<mode> (CmdSaintExpress.java:84-93)
    mkdir -p saintexpress-${mode}
    cp saint_input/* saintexpress-${mode}/
    saint_cmd=(${cmds.argv.collect { t -> shellQuote(t) }.join(' ')})
    printf '%s\\n' "\${saint_cmd[*]}" > saintexpress-${mode}.log
    saint_status=0
    (cd saintexpress-${mode} && "\${saint_cmd[@]}") 2>&1 | tee -a saintexpress-${mode}.log || saint_status=\$?
    if [ "\$saint_status" -ne 0 ]; then
        echo "ERROR: SAINTexpress-${mode}, FragPipe's optional Downstream step (TabDownstream.java:114-145), exited \$saint_status; its output is above (saintexpress-${mode}.log). The run's other results do not depend on it." >&2
        exit "\$saint_status"
    fi
    mv saintexpress-${mode}.log saintexpress-${mode}/saintexpress.log
    """

    stub:
    def cmds = saintexpressCommand(saint_input, mode, saint_cli, cmd_opts, task.ext.args ?: '', task)
    """
    mkdir -p saintexpress-${mode}
    cp saint_input/* saintexpress-${mode}/
    touch saintexpress-${mode}/list.txt
    printf '%s\\n' ${shellQuote(cmds.argv.join(' '))} > saintexpress-${mode}/saintexpress.log
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
