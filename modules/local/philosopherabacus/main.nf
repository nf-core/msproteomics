include { abacusCommand; shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process PHILOSOPHERABACUS {
    tag "${meta.id}"
    label 'process_high'

    input:
    // Multi-experiment report (CmdPhilosopherAbacus.java:73-146), run in the protXML's directory
    // (FragPipe's work dir) on the group directories.
    // group_workspaces: each group's .meta after filter/report, in group order (TreeMap of LcmsFileGroup by name)
    // group_names: the group directory names in the same order; Abacus reads the database of the first
    //   (philosopher lib/aba/pro.go database.RestoreWithPath(args[0]))
    // protxml: combined.prot.xml, or [] when ProteinProphet did not run (philosopher lib/aba/pro.go reads it from the cwd)
    // combined_pepxml: iProphet combined.pep.xml for --peptide (philosopher lib/aba/pep.go), else []
    // abacus_cli: "<--picked/--uniqueonly from the filter text> --razor --reprint --tag <tag> [--protein] [--peptide]"
    // plex: TMT plex for "--labels --plex <n>" with Philosopher TMT extraction, 0 otherwise (CmdPhilosopherAbacus.java:125-129)
    tuple val(meta), path(group_workspaces, stageAs: 'group_meta?/*'), val(group_names), path(protxml), path(combined_pepxml), val(abacus_cli), val(plex)

    output:
    tuple val(meta), path("${prefix}")                         , emit: results_dir
    tuple val(meta), path("${prefix}/reprint.spc.tsv")         , emit: reprint_spc, optional: true
    tuple val(meta), path("${prefix}/reprint.int.tsv")         , emit: reprint_int, optional: true
    tuple val(meta), path("${prefix}/combined_protein.tsv")    , emit: combined_protein, optional: true
    tuple val(meta), path("${prefix}/combined_peptide.tsv")    , emit: combined_peptide, optional: true
    tuple val("${task.process}"), val('philosopher'), val(fragpipeToolVersions().philosopher), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = abacusCommand(group_workspaces, group_names, protxml, abacus_cli, args, plex, task)
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    [ -x ${shellQuote(cmds.argv[0])} ] || { echo "ERROR: ${cmds.argv[0]} not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # Group directories with their workspaces, as in FragPipe's work dir
    ${cmds.stage}

    # Work-dir workspace (CmdPhilosopherWorkspaceCleanInit.java:66-74)
    phi_temp="\$PWD/phi_temp/\$(cat /proc/sys/kernel/random/uuid)"
    mkdir -p "\$phi_temp"
    ${shellQuote(cmds.argv[0])} workspace --init --nocheck --temp "\$phi_temp"

    abacus_cmd=(${cmds.argv.collect { t -> shellQuote(t) }.join(' ')})
    printf '%s\\n' "\${abacus_cmd[*]}" > ${prefix}/abacus.log
    "\${abacus_cmd[@]}" 2>&1 | tee -a ${prefix}/abacus.log

    for f in *.tsv *.csv; do
        if [ -e "\$f" ]; then mv "\$f" ${prefix}/; fi
    done
    rm -rf phi_temp
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = abacusCommand(group_workspaces, group_names, protxml, abacus_cli, task.ext.args ?: '', plex, task)
    def opts = cmds.argv.drop(2)
    """
    mkdir -p ${prefix}
    touch ${prefix}/reprint.spc.tsv ${prefix}/reprint.int.tsv
    ${opts.contains('--protein') ? "touch ${prefix}/combined_protein.tsv" : ''}
    ${opts.contains('--peptide') ? "touch ${prefix}/combined_peptide.tsv" : ''}
    printf '%s\\n' ${shellQuote(cmds.argv.join(' '))} > ${prefix}/abacus.log
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
