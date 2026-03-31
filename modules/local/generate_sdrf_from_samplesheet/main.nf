process GENERATE_SDRF_FROM_SAMPLESHEET {
    tag "generate_sdrf"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    path(samplesheet)
    val(params_json)

    output:
    path("*.sdrf.tsv"), emit: sdrf
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    generate_sdrf.py \\
        --input ${samplesheet} \\
        --params '${params_json}' \\
        --output samplesheet.sdrf.tsv \\
        ${args}
    """

    stub:
    """
    touch samplesheet.sdrf.tsv
    """
}
