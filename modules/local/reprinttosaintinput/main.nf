/*
 * REPRINTTOSAINTINPUT: SAINTexpress inputs from a REPRINT table, as CmdSaintExpress.configure
 * prepares them before running SAINTexpress (CmdSaintExpress.java:56-81, reprintTsvToSaintInput :98-167).
 *
 * Input:
 *   - meta: Aggregate metadata map
 *   - reprint_tsv: reprint.<mode>.tsv of the work dir (Abacus --reprint or IonQuant), or [] when it does not exist
 *   - mode: 'spc' or 'int'
 *
 * Output:
 *   - saint_input (the emit): saintexpress-<mode>/{inter,bait,prey}.dat
 *   - skipped: saintexpress-<mode>/skipped_mode.txt, FragPipe's error for a missing table (the mode is skipped)
 *   - versions: Python version
 */
process REPRINTTOSAINTINPUT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    tuple val(meta), path(reprint_tsv), val(mode)

    output:
    tuple val(meta), path("saintexpress-${mode}/{inter,bait,prey}.dat"), val(mode), emit: saint_input, optional: true
    tuple val(meta), path("saintexpress-${mode}/skipped_mode.txt")     , emit: skipped, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // the task dir stands for FragPipe's work dir: the table is read as reprint.<mode>.tsv there (CmdSaintExpress.java:59)
    """
    reprint_to_saint_input.py reprint.${mode}.tsv saintexpress-${mode}
    """

    stub:
    def outputs = reprint_tsv
        ? "touch saintexpress-${mode}/inter.dat saintexpress-${mode}/bait.dat saintexpress-${mode}/prey.dat"
        : "echo \"Input file \$PWD/reprint.${mode}.tsv is not accessible.\" > saintexpress-${mode}/skipped_mode.txt"
    """
    mkdir -p saintexpress-${mode}
    ${outputs}
    """
}
