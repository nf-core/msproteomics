/*
 * FRAGPIPESECONDPASS: FragPipe's second-pass workflow and manifest (Write sub mzML)
 *
 * When tab-run.write_sub_mzml=true, FragPipe's finalizer writes, after every tool has run
 * (FragpipeRun.java:589-700):
 *   - <wd>/fragpipe-second-pass.workflow: the run's fragpipe.workflow with Write sub mzML off,
 *     Percolator on, PeptideProphet off, Print decoys off, no mass calibration, MSFragger's optimized
 *     values parsed from the console text, the FragPipe databases folder and an empty decoy tag and workdir
 *   - <wd>/fragpipe-files-second-pass.fp-manifest: the run's manifest naming <wd>/<stem>_sub.mzML
 * (TabWorkflow.workflowExt = ".workflow", manifestExt = ".fp-manifest", TabWorkflow.java:138-140).
 * bin/fragpipe_second_pass.py writes both from the staged files; the databases folder is
 * FragpipeLocations.getDirDatabases() of the FragPipe install (fragpipeTools).
 *
 * Input:
 *   - meta: run metadata map
 *   - workflow: the run's fragpipe.workflow (FRAGPIPESAVEWORKFLOW)
 *   - manifest: the run's fragpipe-files.fp-manifest
 *   - logs: MSFragger console output of the run, in FragPipe's order
 *   - workdir: FragPipe's work dir (the published fragpipe/ folder), where the _sub.mzML files are
 *
 * Output:
 *   - workflow: fragpipe-second-pass.workflow
 *   - manifest: fragpipe-files-second-pass.fp-manifest
 *   - versions: python version
 */
include { fragpipeTools } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process FRAGPIPESECONDPASS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    tuple val(meta), path(workflow_file), path(manifest), path(logs, stageAs: 'logs/*/*'), val(workdir)

    output:
    tuple val(meta), path("fragpipe-second-pass.workflow")         , emit: workflow
    tuple val(meta), path("fragpipe-files-second-pass.fp-manifest"), emit: manifest
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def tools = fragpipeTools(task)
    def log_args = [logs].flatten().collect { log_file -> "--log ${log_file}" }.join(' ')
    """
    fragpipe_second_pass.py \\
        --workflow ${workflow_file} \\
        --manifest ${manifest} \\
        ${log_args} \\
        --workdir '${workdir}' \\
        --databases-dir '${tools.databases_dir}' \\
        ${args}
    """

    stub:
    // Stand-ins: the run's workflow unchanged, and its manifest naming <workdir>/<stem>_sub.mzML per LC-MS file
    """
    cp ${workflow_file} fragpipe-second-pass.workflow
    awk 'BEGIN { FS = OFS = "\\t" } { n = \$1; sub(/.*[/]/, "", n); if (n ~ /[.]/) sub(/[.][^.]*\$/, "", n); \$1 = "${workdir}/" n "_sub.mzML"; print }' ${manifest} > fragpipe-files-second-pass.fp-manifest
    """
}
