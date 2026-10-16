/*
 * QPXDIANN: Convert DIA-NN results to the QPX Parquet dataset, .h5mu and .h5ad
 *
 * Adapted from bigbio quantmsdiann modules/bigbio/qpx/diann (v2.3.0), with
 * pg_matrix required (without it qpxc writes no pg view and the MuData has no
 * proteins modality) and the inline MuData step replaced by the guarded
 * bin/qpx_export.py (silent-failure guards, NaN for missing values, dense float32).
 * The q-value filter (--qvalue-threshold) comes from ext.args (conf/modules.config).
 *
 * Input:
 *   - meta: experiment x search-database metadata map (meta.id is the output prefix)
 *   - report: DIA-NN main report (report.tsv or report.parquet)
 *   - pg_matrix: DIA-NN protein-group matrix (pg_matrix.tsv)
 *   - diann_log: DIA-NN log of the final quantification (version detection)
 *   - sdrf: SDRF of the run (run set and obs metadata)
 *
 * Output:
 *   - qpx_dataset: QPX Parquet files
 *   - h5mu: MuData with the precursors and proteins modalities
 *   - h5ad: AnnData of the proteins modality
 *   - summary: key<TAB>value summary (count.* rows are data-dependent)
 *   - versions: qpx, anndata and mudata versions (topic)
 */
process QPXDIANN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qpx:1.1.4--pyhdfd78af_0' :
        'quay.io/biocontainers/qpx:1.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(report), path(pg_matrix), path(diann_log)
    path(sdrf)

    output:
    tuple val(meta), path("qpx_output/*")            , emit: qpx_dataset
    tuple val(meta), path("${prefix}.h5mu")          , emit: h5mu
    tuple val(meta), path("${prefix}.proteins.h5ad") , emit: h5ad
    tuple val(meta), path("${prefix}.summary.tsv")   , emit: summary
    tuple val("${task.process}"), val('qpx'), eval("qpxc --version | sed 's/^qpx //'"), emit: versions_qpx, topic: versions
    tuple val("${task.process}"), val('anndata'), eval("python3 -c \"import importlib.metadata as m; print(m.version('anndata'))\""), emit: versions_anndata, topic: versions
    tuple val("${task.process}"), val('mudata'), eval("python3 -c \"import importlib.metadata as m; print(m.version('mudata'))\""), emit: versions_mudata, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!task.memory) {
        error("QPXDIANN: task.memory is not set, but qpxc convert diann needs --max-memory. Set memory for process_medium (or QPXDIANN) in your config.")
    }
    """
    qpxc convert diann \\
        --report-path ${report} \\
        --sdrf-file ${sdrf} \\
        --pg-matrix-path ${pg_matrix} \\
        --diann-log ${diann_log} \\
        --output-folder qpx_output \\
        --output-prefix ${prefix} \\
        --standardized-intensities \\
        --max-cpus ${task.cpus} \\
        --max-memory ${task.memory.toGiga()}GB \\
        --compression zstd \\
        ${args}

    qpx_export.py \\
        --engine diann \\
        --sdrf-file ${sdrf} \\
        --pg-matrix-file ${pg_matrix} \\
        --qpx-folder qpx_output \\
        --prefix ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p qpx_output
    touch qpx_output/${prefix}.dataset.parquet
    touch qpx_output/${prefix}.feature.parquet
    touch qpx_output/${prefix}.ontology.parquet
    touch qpx_output/${prefix}.pg.parquet
    touch qpx_output/${prefix}.provenance.parquet
    touch qpx_output/${prefix}.run.parquet
    touch qpx_output/${prefix}.sample.parquet
    touch ${prefix}.h5mu
    touch ${prefix}.proteins.h5ad
    printf 'engine\\tdiann\\nmodalities\\tprecursors,proteins\\nobs_key\\trun_stem\\nprefix\\t${prefix}\\n' > ${prefix}.summary.tsv
    """
}
