process THERMORAWFILEPARSER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/thermorawfileparser:2.0.0.dev--h9ee0642_1' :
        'quay.io/biocontainers/thermorawfileparser:2.0.0.dev--h9ee0642_1' }"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("*.{mzML,mzML.gz,mgf,mgf.gz,parquet,parquet.gz}"), emit: spectra
    tuple val("${task.process}"), val('thermorawfileparser'), eval(workflow.stubRun ? "echo 2.0.0.0" : "thermorawfileparser --version"), emit: versions_thermorawfileparser, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--format 0") || args.contains("-f 0") ? "mgf" :
                args.contains("--format 1") || args.contains("-f 1") ? "mzML" :
                args.contains("--format 2") || args.contains("-f 2") ? "mzML" :
                args.contains("--format 3") || args.contains("-f 3") ? "parquet" :
                "mzML"
    suffix = args.contains("--gzip")? "${suffix}.gz" : "${suffix}"

    """
    thermorawfileparser \\
        -i $raw \\
        -b ${prefix}.${suffix} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--format 0") || args.contains("-f 0") ? "mgf" :
                args.contains("--format 1") || args.contains("-f 1") ? "mzML" :
                args.contains("--format 2") || args.contains("-f 2") ? "mzML" :
                args.contains("--format 3") || args.contains("-f 3") ? "parquet" :
                "mzML"
    suffix = args.contains("--gzip")? "${suffix}.gz" : "${suffix}"

    """
    touch ${prefix}.${suffix}
    """
}
