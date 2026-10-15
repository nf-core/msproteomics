include { philosopherSkipReport } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { philosopherFilterCommands; workdirWord; goMemLimitGb; shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process PHILOSOPHERFILTER {
    tag "${meta.id}"
    label 'process_high'

    input:
    // One FragPipe group. The task dir plays FragPipe's work dir: the group runs in
    // <task>/<group> (the blank group "" in <task> itself), the protXML is <task>/<name>.prot.xml.
    // pepxml_files: the group's validated pepXMLs (interact and PTMProphet .mod files); filter reads
    //   the group dir (CmdPhilosopherFilter.java:99-100), Philosopher prefers .mod. files (pep.go:315-322)
    // protxml: combined protXML, or [] when it is not used (dont_use_protxml)
    // first_group_meta: the .meta workspace of FragPipe's first group (sharedPepxmlFiles.firstKey()),
    //   [] for the first group itself; later groups read its db.bin, protxml.bin and razor.bin via
    //   --dbbin/--probin (CmdPhilosopherFilter.java:103-113)
    // first_group: name of the first group ("" = work-dir root)
    // filter_cli: filter options after FragPipe's rewrite (CmdPhilosopherFilter.java:82-97)
    // dont_use_protxml: no --protxml/--probin/--razor (FragpipeRun.java:1716-1755)
    // report_cli: report flags (CmdPhilosopherReport.java:54-62)
    tuple val(meta), path(pepxml_files, stageAs: 'pepxml_inputs/*'), path(protxml), path(first_group_meta, stageAs: 'first_group_meta'), val(first_group), val(filter_cli), val(dont_use_protxml), val(report_cli)
    path fasta
    val decoy_tag    // decoy protein prefix from workflow database.decoy-tag (FragpipeRun.java:1442)

    output:
    tuple val(meta), path("${prefix}/*")                     , emit: results_files
    tuple val(meta), path("${prefix}/.meta")                 , emit: workspace
    // psm.tsv / peptide.tsv: absent when Report runs later in PHILOSOPHERREPORT (ext.skip_report)
    tuple val(meta), path("${prefix}/psm.tsv")               , emit: psms, optional: true
    tuple val(meta), path("${prefix}/peptide.tsv")           , emit: peptides, optional: true
    tuple val(meta), path("${prefix}/protein.tsv")           , emit: proteins, optional: true
    tuple val(meta), path("${prefix}/ion.tsv")               , emit: ions, optional: true
    tuple val(meta), path("${prefix}/philosopher_filter.log"), emit: log
    tuple val("${task.process}"), val('philosopher'), val(fragpipeToolVersions().philosopher), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = goMemLimitGb(task)
    // skip_report: FreeQuant / LabelQuant update the workspace before Report (FragpipeRun.java:2463-2471),
    // so Report runs later in PHILOSOPHERREPORT and this task writes no tables
    def skip_report = philosopherSkipReport(task)
    def cmds = philosopherFilterCommands(meta, pepxml_files, protxml, first_group_meta, first_group, filter_cli, args, dont_use_protxml, report_cli, fasta, decoy_tag, prefix, task)
    """
    export HOME=\$(pwd)
    [ -x ${shellQuote(cmds.philosopher)} ] || { echo "ERROR: ${cmds.philosopher} not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # FragPipe work-dir layout: group dir with its pepXMLs, first group's workspace for --dbbin/--probin
    mkdir -p ${shellQuote(cmds.group_dir)}
    ${cmds.link_pepxml}
    ${cmds.link_first}
    : > ${prefix}.log

    # Group workspace (CmdPhilosopherWorkspaceCleanInit.java:66-74)
    phi_temp="\$PWD/phi_temp/\$(cat /proc/sys/kernel/random/uuid)"
    mkdir -p "\$phi_temp"
    (cd ${shellQuote(cmds.group_dir)} && ${shellQuote(cmds.philosopher)} workspace --init --nocheck --temp "\$phi_temp") 2>&1 | tee -a ${prefix}.log

    # Database annotation in the first group only (CmdPhilosopherDbAnnotate.java:58-74)
    ${cmds.annotate_lines}

    # Filter (CmdPhilosopherFilter.java:79-121); its stdout+stderr is the group's filter.log
    filter_cmd=(${cmds.filter.collect { t -> workdirWord(t) }.join(' ')})
    printf '%s\\n' "\${filter_cmd[*]}" >> ${prefix}.log
    (cd ${shellQuote(cmds.group_dir)} && GOMEMLIMIT=${ram}GiB GOGC=200 GOMAXPROCS=${task.cpus} "\${filter_cmd[@]}") > ${shellQuote(cmds.group_dir)}/filter.log 2>&1 || { cat ${shellQuote(cmds.group_dir)}/filter.log >&2; exit 1; }
    cat ${shellQuote(cmds.group_dir)}/filter.log >> ${prefix}.log

    # Report (CmdPhilosopherReport.java:47-69), unless it runs after FreeQuant / LabelQuant
    if [ "${skip_report}" = "false" ]; then
        report_cmd=(${cmds.report.collect { t -> shellQuote(t) }.join(' ')})
        printf '%s\\n' "\${report_cmd[*]}" >> ${prefix}.log
        (cd ${shellQuote(cmds.group_dir)} && GOMEMLIMIT=${ram}GiB GOGC=200 GOMAXPROCS=${task.cpus} "\${report_cmd[@]}") 2>&1 | tee -a ${prefix}.log
    fi

    rm -rf phi_temp
    ${cmds.move_results}
    # the pepXML links staged into the group dir are inputs, not results
    rm -f ${cmds.pepxml_names.collect { name -> shellQuote("${prefix}/${name}") }.join(' ')}
    mv ${prefix}.log ${prefix}/philosopher_filter.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = philosopherFilterCommands(meta, pepxml_files, protxml, first_group_meta, first_group, filter_cli, task.ext.args ?: '', dont_use_protxml, report_cli, fasta, decoy_tag, prefix, task)
    def skip_report = philosopherSkipReport(task)
    def lines = ["cwd ${cmds.group_dir}"] + (cmds.annotate ? [cmds.annotate.join(' ')] : []) + [cmds.filter.join(' ')] + (skip_report ? [] : [cmds.report.join(' ')])
    """
    # .meta holds Philosopher's workspace (its binary state file, ...); a placeholder file keeps the directory on
    # object storage, which has no object for an empty directory
    mkdir -p ${prefix}/.meta
    touch ${prefix}/.meta/meta.bin
    touch ${prefix}/filter.log
    ${skip_report ? '' : "touch ${prefix}/psm.tsv ${prefix}/peptide.tsv ${prefix}/ion.tsv"}
    ${skip_report || dont_use_protxml ? '' : "touch ${prefix}/protein.tsv"}
    printf '%s\\n' ${lines.collect { line -> shellQuote(line) }.join(' ')} > ${prefix}/philosopher_filter.log
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
