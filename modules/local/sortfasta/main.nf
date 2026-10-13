/*
 * SORTFASTA: the protein-sorted FASTA MSFragger calibrates against in a split database search
 * (FragPipe 24.0 msfragger_pep_split.py sample_fasta(), lines 614-620, 632-633).
 *
 * Input:
 *   - fasta: Protein database
 *
 * Output:
 *   - fasta: sorted_<fasta name>, the proteins sorted bytewise
 *   - versions: Python version (topic: versions)
 */
process SORTFASTA {
    tag "${fasta.name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    path fasta

    output:
    path "sorted_${fasta.name}", emit: fasta
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sort_fasta.py ${fasta} sorted_${fasta.name} ${args}
    """

    stub:
    """
    touch sorted_${fasta.name}
    """
}
