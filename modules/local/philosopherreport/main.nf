include { fragpipeTools } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process PHILOSOPHERREPORT {
    tag "${meta.id}"
    label 'process_medium'

    input:
    // results_dir: the group's FREEQUANT / PHILOSOPHERLABELQUANT directory with the .meta/ workspace
    //              (Filter ran without report/clean); copied to ${prefix}/, the group work dir
    // report_cli: philosopher report options (CmdPhilosopherReport)
    tuple val(meta), path(results_dir, stageAs: 'input_group/*'), val(report_cli)

    output:
    tuple val(meta), path("${prefix}/*")              , emit: results_files
    tuple val(meta), path("${prefix}/psm.tsv")      , emit: psms
    tuple val(meta), path("${prefix}/peptide.tsv")   , emit: peptides
    tuple val(meta), path("${prefix}/protein.tsv")   , emit: proteins
    tuple val(meta), path("${prefix}/ion.tsv")       , emit: ions, optional: true
    tuple val("${task.process}"), val('philosopher'), val(fragpipeToolVersions().philosopher), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!task.memory) {
        error("${task.process} needs task.memory (JVM / Go memory limit)")
    }
    def ram = Math.max(1, task.memory.toGiga() - 2)
    def tools = fragpipeTools(task)

    """
    export HOME=\$(pwd)
    WORK_DIR=\$(pwd)
    mkdir -p ${prefix}

    PHILOSOPHER=${shellQuote(tools.philosopher)}
    [ -x "\$PHILOSOPHER" ] || { echo "ERROR: \$PHILOSOPHER not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }

    # Version capture
    echo ${shellQuote(tools.fragpipe_version)} > .fragpipe_version

    export GOMEMLIMIT="${ram}GiB"
    export GOGC=200
    export GOMAXPROCS=${task.cpus}

    # Report and the final workspace clean run in the group work dir on the workspace Filter,
    # FreeQuant and LabelQuant built (FragpipeRun.java:2463-2471, 2573-2581): no re-init or
    # database annotation (the workspace already holds the annotated database).
    cp -rL ${results_dir}/. ${prefix}/
    # Philosopher keeps the workspace's absolute paths in the .meta workspace binary; point them at this copy, and its
    # Temp dir (in the task that made the workspace) at one in this task
    mkdir -p phi_temp
    philosopher_relocate_workspace.py "\$PWD"/${prefix} "\$PWD/phi_temp"
    cd ${prefix}
    REPORT_CMD=(\$PHILOSOPHER report ${report_cli} ${args})
    printf '%s\\n' "\${REPORT_CMD[*]}" > \${WORK_DIR}/philosopher_report.log
    "\${REPORT_CMD[@]}" 2>&1 | tee -a \${WORK_DIR}/philosopher_report.log
    \$PHILOSOPHER workspace --clean --nocheck 2>&1
    cd \${WORK_DIR}
    rm -rf phi_temp
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    """
    mkdir -p ${prefix}
    cp -rL ${results_dir}/. ${prefix}/
    rm -rf ${prefix}/.meta
    touch ${prefix}/protein.tsv
    touch ${prefix}/peptide.tsv
    touch ${prefix}/psm.tsv
    printf '%s\\n' "philosopher report ${report_cli}" > philosopher_report.log
    echo ${shellQuote(tools.fragpipe_version)} > .fragpipe_version
    """
}
