process TDF2MZML {
    tag "$meta.id"
    label 'process_single'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/bigbio/tdf2mzml@sha256:98a38d2e9d7803e75dfdfdd7282567fec6a4f57df6859c67e0307483c9a710e9' :
        'quay.io/bigbio/tdf2mzml@sha256:98a38d2e9d7803e75dfdfdd7282567fec6a4f57df6859c67e0307483c9a710e9' }"

    input:
    tuple val(meta), path(rawfile)

    output:
    tuple val(meta), path("${rawfile.baseName}.mzML"), emit: mzmls_converted
    path "versions.yml",   emit: versions, topic: versions
    path "*.log",   emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    // No Conda package: the image bundles Bruker's proprietary libtimsdata.so (TDF SDK), and tdf2mzml is
    // neither on Bioconda nor on PyPI
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("TDF2MZML does not support Conda. Use Docker, Singularity or Podman instead.")
    }
    def args = task.ext.args ?: ''

    """
    echo "Converting..." | tee --append ${rawfile.baseName}_conversion.log
    tdf2mzml.py -i ${rawfile} $args 2>&1 | tee --append ${rawfile.baseName}_conversion.log
    mv *.mzml ${rawfile.baseName}.mzML

    # tdf2mzml.py has no --version: the image's tool writes MAJOR_VERSION + "-" + MINOR_VERSION
    # ('0.3', '20200211.1', tdf2mzml.py:29-30,504) as its mzML software version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tdf2mzml.py: 0.3-20200211.1
    END_VERSIONS
    """

    stub:
    """
    touch ${rawfile.baseName}.mzML
    echo "Stub execution" > ${rawfile.baseName}_conversion.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tdf2mzml.py: stub
    END_VERSIONS
    """
}
