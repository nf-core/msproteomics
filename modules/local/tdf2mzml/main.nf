process TDF2MZML {
    tag "$meta.mzml_id"
    label 'process_single'
    label 'error_retry'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/bigbio/tdf2mzml@sha256:98a38d2e9d7803e75dfdfdd7282567fec6a4f57df6859c67e0307483c9a710e9' :
        'quay.io/bigbio/tdf2mzml@sha256:98a38d2e9d7803e75dfdfdd7282567fec6a4f57df6859c67e0307483c9a710e9' }"

    input:
    tuple val(meta), path(rawfile)

    output:
    tuple val(meta), path("*.mzML"), emit: mzmls_converted
    tuple val(meta), path("*.d"),   emit: dotd_files
    path "versions.yml",   emit: versions
    path "*.log",   emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def target_name = file(rawfile.baseName).baseName

    """
    echo "Converting..." | tee --append ${rawfile.baseName}_conversion.log
    tdf2mzml.py -i *.d $args 2>&1 | tee --append ${rawfile.baseName}_conversion.log
    mv *.mzml ${target_name}.mzML
    mv *.d ${target_name}.d

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tdf2mzml.py: \$(tdf2mzml.py --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"

    """
    touch ${prefix}.mzML
    mkdir -p ${prefix}.d
    echo "Stub execution" > ${rawfile.baseName}_conversion.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tdf2mzml.py: stub
    END_VERSIONS
    """
}
