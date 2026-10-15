/*
 * PHILOSOPHERWORKSPACECLEAN: FragPipe's "Phi wrk clean" of one group work dir (CmdPhilosopherWorkspaceClean)
 *
 * FragPipe runs `philosopher workspace --clean --nocheck` in every Philosopher work dir after the
 * last Philosopher-dependent command (FragpipeRun.java:2557-2583, CmdPhilosopherWorkspaceClean.java:43-48);
 * Philosopher then removes the .meta workspace (lib/wrk/wrk.go:133-142).
 *
 * Input:
 *   - meta: Group metadata map
 *   - results_dir: the group directory with its .meta workspace (Filter + Report done)
 *
 * Output:
 *   - results_dir: the group directory without the .meta workspace
 *   - versions_philosopher / versions_fragpipe: versions (topic: versions)
 */
include { fragpipeTools } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process PHILOSOPHERWORKSPACECLEAN {
    tag "${meta.id}"
    label 'process_single'

    input:
    // results_dir: the group's final Philosopher directory with the .meta/ workspace; copied to ${prefix}/,
    //              the group work dir (the staged input is a link to the upstream task's directory)
    tuple val(meta), path(results_dir, stageAs: 'input_group/*')

    output:
    tuple val(meta), path("${prefix}")                          , emit: results_dir
    tuple val("${task.process}"), val('philosopher'), val(fragpipeToolVersions().philosopher), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    """
    export HOME=\$(pwd)
    PHILOSOPHER=${shellQuote(tools.philosopher)}
    [ -x "\$PHILOSOPHER" ] || { echo "ERROR: \$PHILOSOPHER not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    echo ${shellQuote(tools.fragpipe_version)} > .fragpipe_version

    mkdir -p ${prefix}
    cp -rL ${results_dir}/. ${prefix}/
    # philosopher workspace --clean --nocheck in the group work dir (CmdPhilosopherWorkspaceClean.java:43-48)
    (cd ${prefix} && "\$PHILOSOPHER" workspace --clean --nocheck ${args}) 2>&1
    # Clean removes .meta only when it holds meta.bin (wrk.go:133-142); a workspace left behind is an error
    [ ! -e ${prefix}/.meta ] || { echo "ERROR: philosopher workspace --clean left ${prefix}/.meta" >&2; exit 1; }
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    cp -rL ${results_dir}/. ${prefix}/
    rm -rf ${prefix}/.meta
    echo ${shellQuote(fragpipeTools(task).fragpipe_version)} > .fragpipe_version
    """
}
