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
    tuple val(meta), path("${prefix}/psm.tsv")               , emit: psms
    tuple val(meta), path("${prefix}/peptide.tsv")           , emit: peptides
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
    // so Report runs later in PHILOSOPHERREPORT; the tables here are placeholders
    def skip_report = philosopherSkipReport(task)
    def cmds = philosopherFilterCommands(meta, pepxml_files, protxml, first_group_meta, first_group, filter_cli, args, dont_use_protxml, report_cli, fasta, decoy_tag, task)
    def link_pepxml = cmds.pepxml_names.collect { name -> "ln -s \"\$PWD\"/pepxml_inputs/${shellQuote(name)} ${shellQuote(cmds.group_dir)}/${shellQuote(name)}" }.join('\n    ')
    def link_first = cmds.is_first ? '' : "mkdir -p ${shellQuote(cmds.first_dir)} && ln -s \"\$PWD\"/first_group_meta ${shellQuote(cmds.first_dir)}/.meta"
    def annotate = cmds.annotate ? """annotate_cmd=(${cmds.annotate.collect { t -> workdirWord(t) }.join(' ')})
    printf '%s\\n' "\${annotate_cmd[*]}" >> ${prefix}.log
    (cd ${shellQuote(cmds.group_dir)} && GOMEMLIMIT=${ram}GiB GOGC=100 GOMAXPROCS=${task.cpus} "\${annotate_cmd[@]}") 2>&1 | tee -a ${prefix}.log""" : ''
    """
    export HOME=\$(pwd)
    [ -x ${shellQuote(cmds.philosopher)} ] || { echo "ERROR: ${cmds.philosopher} not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # FragPipe work-dir layout: group dir with its pepXMLs, first group's workspace for --dbbin/--probin
    mkdir -p ${shellQuote(cmds.group_dir)}
    ${link_pepxml}
    ${link_first}
    : > ${prefix}.log

    # Group workspace (CmdPhilosopherWorkspaceCleanInit.java:66-74)
    phi_temp="\$PWD/phi_temp/\$(cat /proc/sys/kernel/random/uuid)"
    mkdir -p "\$phi_temp"
    (cd ${shellQuote(cmds.group_dir)} && ${shellQuote(cmds.philosopher)} workspace --init --nocheck --temp "\$phi_temp") 2>&1 | tee -a ${prefix}.log

    # Database annotation in the first group only (CmdPhilosopherDbAnnotate.java:58-74)
    ${annotate}

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
    else
        (cd ${shellQuote(cmds.group_dir)} && touch psm.tsv peptide.tsv)
    fi

    rm -rf phi_temp
    ${cmds.group_dir == '.' ? "mkdir -p ${shellQuote(prefix)}; for f in psm.tsv peptide.tsv ion.tsv protein.tsv protein.fas filter.log; do if [ -e \"\$f\" ]; then mv \"\$f\" ${shellQuote(prefix)}/; fi; done; mv .meta ${shellQuote(prefix)}/" : (cmds.group_dir == prefix ? '' : "mv ${shellQuote(cmds.group_dir)} ${shellQuote(prefix)}")}
    # the pepXML links staged into the group dir are inputs, not results
    rm -f ${cmds.pepxml_names.collect { name -> shellQuote("${prefix}/${name}") }.join(' ')}
    mv ${prefix}.log ${prefix}/philosopher_filter.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = philosopherFilterCommands(meta, pepxml_files, protxml, first_group_meta, first_group, filter_cli, task.ext.args ?: '', dont_use_protxml, report_cli, fasta, decoy_tag, task)
    def skip_report = philosopherSkipReport(task)
    def lines = ["cwd ${cmds.group_dir}"] + (cmds.annotate ? [cmds.annotate.join(' ')] : []) + [cmds.filter.join(' ')] + (skip_report ? [] : [cmds.report.join(' ')])
    """
    # .meta holds Philosopher's workspace (its binary state file, ...); a placeholder file keeps the directory on
    # object storage, which has no object for an empty directory
    mkdir -p ${prefix}/.meta
    touch ${prefix}/.meta/${'meta'}.bin
    touch ${prefix}/psm.tsv ${prefix}/peptide.tsv ${prefix}/ion.tsv ${prefix}/filter.log
    ${dont_use_protxml ? '' : "touch ${prefix}/protein.tsv"}
    printf '%s\\n' ${lines.collect { line -> shellQuote(line) }.join(' ')} > ${prefix}/philosopher_filter.log
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
