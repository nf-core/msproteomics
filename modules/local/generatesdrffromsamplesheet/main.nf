/*
 * GENERATESDRFFROMSAMPLESHEET: Bookkeeping SDRF from the CSV samplesheet (bin/generate_sdrf.py)
 *
 * Input:
 *   - samplesheet: the pipeline's CSV samplesheet
 *   - params_json: experiment-level values (generateSdrfBookkeeping), passed as one shell word
 *   - database: searched FASTA (organism guess when the organism is empty)
 *
 * Output:
 *   - sdrf: samplesheet.sdrf.tsv
 *   - versions_python: Python version (topic: versions)
 */
include { shellQuote } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process GENERATESDRFFROMSAMPLESHEET {
    tag "generate_sdrf"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    path(samplesheet)
    val(params_json)
    path(database)

    output:
    path("*.sdrf.tsv"), emit: sdrf
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    generate_sdrf.py \\
        --input ${samplesheet} \\
        --params ${shellQuote(params_json)} \\
        --output samplesheet.sdrf.tsv \\
        --database ${database} \\
        ${args}
    """

    stub:
    """
    touch samplesheet.sdrf.tsv
    """
}
