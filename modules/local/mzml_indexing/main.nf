process MZML_INDEXING {
    tag "$meta.mzml_id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/bigbio/openms-tools-thirdparty-sif:2025.04.14' :
        'ghcr.io/bigbio/openms-tools-thirdparty:2025.04.14' }"

    input:
    tuple val(meta), path(mzmlfile)

    output:
    tuple val(meta), path("out/*.mzML"), emit: mzmls_indexed
    path "versions.yml", emit: versions
    path "*.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p out
    FileConverter -in ${mzmlfile} -out out/${mzmlfile.baseName}.mzML $args 2>&1 | tee ${mzmlfile.baseName}_mzmlindexing.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FileConverter: \$(FileConverter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p out
    touch out/${mzmlfile.baseName}.mzML
    echo "Stub execution" > ${mzmlfile.baseName}_mzmlindexing.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FileConverter: stub
    END_VERSIONS
    """
}
