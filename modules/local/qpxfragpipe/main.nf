/*
 * QPXFRAGPIPE: FragPipe LFQ results -> QPX Parquet, .h5mu and proteins .h5ad
 *
 * Runs three steps and fails the task on any error:
 *   1. qpx_export.py --precheck: guard G3 (MaxLFQ / Total columns for every experiment),
 *      because qpxc exits 0 while silently dropping the protein view without them.
 *   2. qpxc convert fragpipe: QPX Parquet dataset in qpx_output/.
 *   3. qpx_export.py: guards G1-G11, NaN for missing values, then writes the outputs.
 * Runs map to IonQuant experiments through FragPipe's experiment_annotation.tsv, which both
 * pipeline mode (IONQUANT) and headless mode (FragPipe's results/) write.
 *
 * Input:
 *   - meta: Aggregate metadata map (meta.id is the output prefix)
 *   - combined_protein: IonQuant combined_protein.tsv
 *   - combined_ion: IonQuant combined_ion.tsv
 *   - experiment_annotation: FragPipe experiment_annotation.tsv (file, sample, sample_name, condition, replicate)
 *   - sdrf: SDRF of the run (expected run set, obs metadata)
 *
 * Output:
 *   - qpx_dataset: QPX Parquet dataset (qpx_output/*)
 *   - h5mu: MuData with the precursors and proteins modalities
 *   - h5ad: AnnData of the proteins modality
 *   - summary: sorted key<TAB>value summary TSV
 *   - versions: qpx, anndata and mudata versions (topic)
 */
process QPXFRAGPIPE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qpx:1.1.5--pyhdfd78af_0' :
        'quay.io/biocontainers/qpx:1.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(combined_protein), path(combined_ion), path(experiment_annotation)
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
    """
    qpx_export.py \\
        --engine fragpipe \\
        --precheck \\
        --pg-file ${combined_protein} \\
        --sdrf-file ${sdrf} \\
        --experiment-annotation-file ${experiment_annotation}

    qpxc convert fragpipe \\
        --ion-file ${combined_ion} \\
        --pg-file ${combined_protein} \\
        --sdrf-file ${sdrf} \\
        --experiment-annotation-file ${experiment_annotation} \\
        --max-cpus ${task.cpus} \\
        --output-folder qpx_output \\
        --output-prefix ${prefix} \\
        ${args}

    qpx_export.py \\
        --engine fragpipe \\
        --sdrf-file ${sdrf} \\
        --qpx-folder qpx_output \\
        --prefix ${prefix} \\
        --pg-file ${combined_protein} \\
        --experiment-annotation-file ${experiment_annotation}
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
    printf 'engine\\tfragpipe\\nprefix\\t${prefix}\\n' > ${prefix}.summary.tsv
    """
}
