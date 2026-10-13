process PMULTIQC {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pmultiqc:0.0.39--pyhdfd78af_0' :
        'quay.io/biocontainers/pmultiqc:0.0.39--pyhdfd78af_0' }"

    input:
    tuple val(meta), path('results/*'), path(multiqc_config)

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*.db")  , emit: db, optional: true
    tuple val(meta), path("*_data"), emit: data
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    multiqc \\
        -f \\
        --config ${multiqc_config} \\
        ${args} \\
        ./results \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pmultiqc: \$(multiqc --pmultiqc_version 2>&1 | sed -e "s/pmultiqc, version //g")
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_report_data
    touch multiqc_report_data/multiqc.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pmultiqc: stub
    END_VERSIONS
    """
}
