/*
 * CHECKDBCONFIG: FragPipe's database decoy check (FragpipeRun.checkDbConfig)
 *
 * Input:
 *   - meta: Database metadata map
 *   - fasta: Protein database that downstream tools will search
 *   - decoy_tag: Workflow database.decoy-tag
 *
 * Output:
 *   - fasta: The unchanged input database, emitted only when it passes the check
 *   - versions_python: Python version (topic: versions)
 */
process CHECKDBCONFIG {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    tuple val(meta), path(fasta)
    val decoy_tag    // decoy protein prefix from workflow database.decoy-tag (FragpipeRun.java:1442)

    output:
    tuple val(meta), path(fasta, includeInputs: true), emit: fasta
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    check_db_config.py \\
        --fasta ${fasta} \\
        --decoy-tag '${decoy_tag}' \\
        ${args}
    """

    stub:
    """
    echo "stub: decoy check skipped for ${fasta} (decoy tag '${decoy_tag}')"
    """
}
