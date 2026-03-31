process MZML_STATISTICS {
    tag "$meta.mzml_id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quantms-utils:0.0.24--pyh7e72e81_0' :
        'biocontainers/quantms-utils:0.0.24--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(ms_file)

    output:
    path "*_ms_info.parquet", emit: ms_statistics
    tuple val(meta), path("*_ms2_info.parquet"), emit: ms2_statistics, optional: true
    path "*_feature_info.parquet", emit: feature_statistics, optional: true
    path "versions.yml", emit: versions
    path "*.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"
    def string_ms2_file = task.ext.id_only == true || task.ext.mzml_features == true ? "--ms2_file" : ""
    def string_features_file = task.ext.mzml_features == true ? "--feature_detection" : ""

    """
    quantmsutilsc mzmlstats --ms_path "${ms_file}" \\
        ${string_ms2_file} \\
        ${string_features_file} \\
        2>&1 | tee ${ms_file.baseName}_mzml_statistics.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quantms-utils: \$(pip show quantms-utils 2>/dev/null | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"

    """
    touch ${prefix}_ms_info.parquet
    touch ${prefix}_ms2_info.parquet
    touch ${prefix}_feature_info.parquet
    echo "Stub execution" > ${prefix}_mzml_statistics.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quantms-utils: stub
    END_VERSIONS
    """
}
