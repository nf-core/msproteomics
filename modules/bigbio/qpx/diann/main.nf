process QPX_DIANN {
    tag "qpx_diann"
    label 'process_medium'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qpx:1.1.5--pyhdfd78af_0' :
        'quay.io/biocontainers/qpx:1.1.5--pyhdfd78af_0' }"

    input:
    path(diann_report)
    path(pg_matrix)
    path(sdrf)
    path(diann_log)
    val(project_accession)
    // Optional: the FASTA used for the search, or [] when it is not available.
    // qpx fills null pg.sequence_coverage / pg.molecular_weight and
    // feature.pg_positions from it, for target rows only, never overwriting a
    // producer value; proteins absent from it stay null.
    path(fasta)

    output:
    path "qpx_output/*", emit: qpx_dataset
    // Optional: the parquet views are the dataset's source of truth, so a MuData
    // view that cannot be built must not destroy an otherwise complete run.
    path "*.h5mu"      , emit: mudata, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = project_accession ?: 'diann'
    def pg_arg  = pg_matrix ? "--pg-matrix-path ${pg_matrix}" : ''
    def log_arg = diann_log ? "--diann-log ${diann_log}" : ''
    def acc_arg = project_accession ? "--project-accession ${project_accession}" : ''
    def fasta_arg = fasta ? "--fasta ${fasta}" : ''
    // No --qvalue-threshold here: a cutoff, if any, comes from ext.args
    """
    set -o pipefail
    qpxc convert diann \\
        --report-path ${diann_report} \\
        --sdrf-file ${sdrf} \\
        ${pg_arg} \\
        ${log_arg} \\
        ${acc_arg} \\
        ${fasta_arg} \\
        --output-folder qpx_output \\
        --output-prefix ${prefix} \\
        --standardized-intensities \\
        --max-cpus ${task.cpus} \\
        --max-memory ${task.memory ? task.memory.toGiga() : 4}GB \\
        --compression zstd \\
        ${args}

    # qpx >= 1.1.5 writes the MuData view itself at the end of `convert`
    # (--mudata, the default), via write_dataset_mudata: it refuses a MuData
    # missing a required quantification modality (bigbio/qpx#316) and is
    # best-effort, so a view that cannot be built leaves the parquet views
    # intact. Building it again here would double the time/memory on large
    # datasets, so only move it out of qpx_output (keeping it out of the
    # qpx_dataset emit). `--no-mudata` in ext.args skips the view entirely.
    if [ -f "qpx_output/${prefix}.h5mu" ]; then
        mv "qpx_output/${prefix}.h5mu" "${prefix}.h5mu"
        echo "MuData -> ${prefix}.h5mu"
    else
        echo "WARNING: no MuData view was written; see the log above. The qpx_output parquet views are complete and authoritative."
    fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    qpx: \$(qpxc --version 2>&1 | sed 's/^qpx //')
    mudata: \$(python -c 'import mudata; print(mudata.__version__)')
END_VERSIONS
    """

    stub:
    def prefix = project_accession ?: 'diann'
    """
    mkdir -p qpx_output
    touch qpx_output/stub.parquet
    touch ${prefix}.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qpx: stub
        mudata: stub
    END_VERSIONS
    """
}
