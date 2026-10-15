/*
 * FRAGPIPESDRF: Write the FragPipe SDRF template (sdrf.tsv) for a run.
 *
 * Port of the FragPipe 24.0 writer triggered by workflow.misc.save-sdrf=true
 * (FragpipeRun.java, TabWorkflow.sdrfSave, SDRFtable.java, TabMsfragger SDRF helpers).
 *
 * Input:
 *   - meta: Run metadata map
 *   - workflow_file: FragPipe's saved UI state fragpipe.workflow (FRAGPIPESAVEWORKFLOW: every msfragger.* key, workflow.misc.sdrf-type)
 *   - manifest: fragpipe-files.fp-manifest (path, experiment, bioreplicate, data type)
 *   - searched_names: manifest file name -> name of the file MSFragger searched (TSV)
 *   - fasta: Protein database searched (organism guess from UniProt OS= headers)
 *   - unimod: unimod.obo shipped in the FragPipe tools/ folder
 *   - logs: MSFragger log(s), concatenated in the given order (calibrated fragment tolerance, instrument)
 *   - tmt_label: TMT-Integrator QuantLabel name (e.g. 'TMT-16') when TMT-Integrator runs, else ''
 *
 * Output:
 *   - sdrf: sdrf.tsv
 *   - versions: Software versions (topic channel)
 */
process FRAGPIPESDRF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    tuple val(meta), path(workflow_file), path(manifest), path(searched_names), path(fasta), path(unimod), path(logs, stageAs: 'logs/*/*'), val(tmt_label)

    output:
    tuple val(meta), path("sdrf.tsv"), emit: sdrf
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def log_args = [logs].flatten().collect { log_file -> "--log ${log_file}" }.join(' ')
    """
    fragpipe_sdrf.py \\
        --workflow ${workflow_file} \\
        --manifest ${manifest} \\
        --searched-names ${searched_names} \\
        --fasta ${fasta} \\
        --unimod ${unimod} \\
        --tmt-label '${tmt_label ?: ''}' \\
        ${log_args} \\
        --output sdrf.tsv \\
        ${args}
    """

    stub:
    """
    touch sdrf.tsv
    """
}
