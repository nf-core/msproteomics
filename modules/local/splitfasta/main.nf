/*
 * SPLITFASTA: Split a FASTA file into M chunks at protein boundaries.
 *
 * Replicates get_fasta_offsets() and set_up_directories() from
 * FragPipe's msfragger_pep_split.py for split-database search.
 *
 * Input:
 *   - meta: Sample/experiment metadata map
 *   - fasta: Protein sequence database in FASTA format
 *   - num_chunks: Number of chunks to split into (integer >= 1)
 *
 * Output:
 *   - fasta_chunks: Tuple of (meta, chunk FASTA files) in split_db/{0,1,...}/
 *   - versions: Software versions
 */
process SPLITFASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    tuple val(meta), path(fasta)
    val num_chunks

    output:
    tuple val(meta), path("split_db/*/*.fasta"), emit: fasta_chunks
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    split_fasta.py \\
        --fasta ${fasta} \\
        --num_chunks ${num_chunks} \\
        --outdir split_db \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    for i in \$(seq 0 \$((${num_chunks} - 1))); do
        mkdir -p split_db/\$i
        echo ">dummy_protein_chunk_\${i}" > split_db/\$i/${fasta.name}
        echo "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERM" >> split_db/\$i/${fasta.name}
    done
    """
}
