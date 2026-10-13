include { fragpipeTools } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process FREEQUANT {
    tag "${meta.id}"
    label 'process_medium'

    input:
    // FreeQuant performs label-free quantification using Philosopher
    // results_dir: the group's PHILOSOPHERFILTER directory with its .meta/ workspace (filter ran
    //              without report/clean); copied to ${prefix}/, FragPipe's group work dir
    // lcms_files: LC-MS files of the FragPipe group, staged into lcms_files/ (the --dir FragPipe passes)
    // config_cli: FreeQuant options (--ptw <rt-tol> --tol <mz-tol>, QuantPanelLabelfree.getFreequantOptsAsText)
    tuple val(meta), path(results_dir, stageAs: 'input_group/*'), path(lcms_files, stageAs: 'lcms_files/*'), val(config_cli)

    output:
    tuple val(meta), path("${prefix}")                          , emit: results_dir
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

    # The image's binary runs in place, as FragPipe runs it (golden FreeQuant argv)
    PHILOSOPHER="${tools.philosopher}"
    echo "${tools.fragpipe_version}" > .fragpipe_version

    export GOMAXPROCS=${task.cpus}
    export GOMEMLIMIT="${ram}GiB"
    export GOGC=200

    # The group work dir with the workspace Filter left (FragPipe runs FreeQuant there, no re-init)
    mkdir -p ${prefix}
    cp -rL ${results_dir}/. ${prefix}/
    # Philosopher keeps the workspace's absolute paths in the .meta workspace binary; point them at this copy
    philosopher_relocate_workspace.py ${prefix}

    # philosopher freequant <opts> --dir <LC-MS dir> in the group dir (CmdFreequant.java:108-128).
    # FragPipe's --raw is unreachable: FreeQuant takes only mzML (CmdFreequant.java:44, checked in
    # FRAGPIPE_QUANT).
    CMD=(\$PHILOSOPHER freequant ${config_cli} --dir "\${WORK_DIR}/lcms_files" ${args})
    printf '%s\\n' "\${CMD[*]}" > \${WORK_DIR}/freequant.log
    (cd ${prefix} && "\${CMD[@]}") 2>&1 | tee -a \${WORK_DIR}/freequant.log

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    cp -rL ${results_dir}/. ${prefix}/
    printf '%s\\n' "philosopher freequant ${config_cli} --dir lcms_files" > freequant.log
    echo "${fragpipeTools(task).fragpipe_version}" > .fragpipe_version
    """
}
