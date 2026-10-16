/*
 * QPXEXPORT: .h5mu and proteins .h5ad from a DIA-NN QPX Parquet dataset (QPX_DIANN, bigbio/nf-modules)
 *
 * Runs bin/qpx_export.py on the dataset qpxc convert diann wrote: silent-failure guards
 * (G1, G2, G4-G10, G16, G17: modalities, obs = the SDRF runs, values against pg_matrix), then the
 * MuData with missing values as NaN in a dense float32 matrix (qpx's own MuData view stores them as
 * sparse zeros) and the proteins AnnData. The pr_matrix filter (--matrix-qvalue) comes from ext.args.
 *
 * Input:
 *   - meta: Metadata map (meta.id is the dataset prefix, the qpxc --output-prefix)
 *   - qpx_dataset: QPX Parquet files of the dataset
 *   - pg_matrix: DIA-NN protein-group matrix (pg_matrix.tsv)
 *   - sdrf: SDRF of the run (run set and obs metadata)
 *
 * Output:
 *   - h5mu: MuData with the precursors and proteins modalities
 *   - h5ad: AnnData of the proteins modality
 *   - summary: key<TAB>value summary (count.* rows are data-dependent)
 *   - versions: qpx, anndata and mudata versions (topic)
 */
process QPXEXPORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qpx:1.1.5--pyhdfd78af_0' :
        'quay.io/biocontainers/qpx:1.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(qpx_dataset, stageAs: 'qpx_output/*'), path(pg_matrix)
    path(sdrf)

    output:
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
    """
    qpx_export.py \\
        --engine diann \\
        --sdrf-file ${sdrf} \\
        --pg-matrix-file ${pg_matrix} \\
        --qpx-folder qpx_output \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5mu
    touch ${prefix}.proteins.h5ad
    printf 'engine\\tdiann\\nmodalities\\tprecursors,proteins\\nobs_key\\trun_stem\\nprefix\\t${prefix}\\n' > ${prefix}.summary.tsv
    """
}
