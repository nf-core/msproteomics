include { fragpipeTools; labelquantAnnotationName } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process PHILOSOPHERLABELQUANT {
    tag "${meta.id}"
    label 'process_high'

    input:
    // Philosopher labelquant extracts TMT/isobaric reporter ion intensities (CmdLabelquant.java:62-165).
    // results_dir: the group's FREEQUANT directory (.meta/ workspace); copied to ${prefix}/, the group work dir
    // lcms_files: mzML files of the FragPipe group, staged into lcms_files/ (the --dir FragPipe passes)
    // annotation: <group>_annotation.txt ("channel sample", labelquantAnnotationText), staged into the LC-MS
    //             directory, where FragPipe requires it (CmdLabelquant.java:86-91), and copied into the group dir for Report
    // config_cli: --tol --level --minprob --purity --removelow --plex (parser, CmdLabelquant.java:98-116)
    // brand: label type in lower case (--brand, CmdLabelquant.java:144-145)
    tuple val(meta), path(results_dir, stageAs: 'input_group/*'), path(lcms_files, stageAs: 'lcms_files/*'), path(annotation, stageAs: 'lcms_files/*'), val(config_cli), val(brand)

    output:
    tuple val(meta), path("${prefix}")              , emit: results_dir
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
    def annotation_name = labelquantAnnotationName(meta)
    if (annotation.toString().tokenize("/").last() != annotation_name) {
        error("${task.process} ${meta.id}: the TMT annotation must be ${annotation_name} (CmdLabelquant.java:86-91), got ${annotation}")
    }

    """
    export HOME=\$(pwd)
    WORK_DIR=\$(pwd)
    PHILOSOPHER=${shellQuote(tools.philosopher)}
    [ -x "\$PHILOSOPHER" ] || { echo "ERROR: \$PHILOSOPHER not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }

    # Version capture
    echo ${shellQuote(tools.fragpipe_version)} > .fragpipe_version

    export GOMAXPROCS=${task.cpus}
    export GOMEMLIMIT="${ram}GiB"
    export GOGC=200

    mkdir -p ${prefix}
    cp -rL ${results_dir}/. ${prefix}/
    # Philosopher keeps the workspace's absolute paths in the .meta workspace binary; point them at this copy
    philosopher_relocate_workspace.py ${prefix}
    cp lcms_files/${annotation_name} ${prefix}/

    # philosopher labelquant <opts> --annot <abs> --brand <type> --dir <LC-MS dir> in the group dir
    # (CmdLabelquant.java:94-160; --raw is unreachable: labelquant takes only mzML, CmdLabelquant.java:51)
    CMD=("\$PHILOSOPHER" labelquant ${config_cli} --annot "\${WORK_DIR}"/lcms_files/${shellQuote(annotation_name)} --brand ${shellQuote(brand)} --dir "\${WORK_DIR}/lcms_files" ${args})
    printf '%s\\n' "\${CMD[*]}" > labelquant.log
    (cd ${prefix} && "\${CMD[@]}") 2>&1 | tee -a labelquant.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    def annotation_name = labelquantAnnotationName(meta)
    """
    mkdir -p ${prefix}
    cp -rL ${results_dir}/. ${prefix}/
    cp lcms_files/${annotation_name} ${prefix}/
    printf '%s\\n' "philosopher labelquant ${config_cli} --annot lcms_files/${annotation_name} --brand ${brand} --dir lcms_files" > labelquant.log
    echo ${shellQuote(tools.fragpipe_version)} > .fragpipe_version
    """
}
