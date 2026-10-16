/*
 * TMTANNDATA: TMT-Integrator abundance reports -> .h5mu and .h5ad
 *
 * Runs tmt_to_anndata.py on every abundance_<level>_<norm>.tsv TMT-Integrator wrote.
 * Guards G4, G8, G9, G12-G15 fail the task before anything is written.
 *
 * Input:
 *   - meta: Aggregate metadata map (meta.id is the output prefix)
 *   - abundance: TMT-Integrator abundance_<level>_<norm>.tsv files
 *   - annotation: multi-plex TMT annotation (plex channel sample sample_name condition replicate)
 *   - log2transformed: TMT-Integrator log2transformed setting ('true' / 'false')
 *
 * Output:
 *   - h5mu: one MuData per normalization, one modality per level
 *   - h5ad: one AnnData per level x normalization
 *   - summary: sorted key<TAB>value summary TSV
 *   - versions: anndata and mudata versions (topic)
 */
process TMTANNDATA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qpx:1.1.4--pyhdfd78af_0' :
        'quay.io/biocontainers/qpx:1.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(abundance), path(annotation), val(log2transformed)

    output:
    tuple val(meta), path("${prefix}.tmt_*.h5mu")     , emit: h5mu
    tuple val(meta), path("${prefix}.tmt_*.h5ad")     , emit: h5ad
    tuple val(meta), path("${prefix}.tmt.summary.tsv"), emit: summary
    tuple val("${task.process}"), val('anndata'), eval("python3 -c \"import importlib.metadata as m; print(m.version('anndata'))\""), emit: versions_anndata, topic: versions
    tuple val("${task.process}"), val('mudata'), eval("python3 -c \"import importlib.metadata as m; print(m.version('mudata'))\""), emit: versions_mudata, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    tmt_to_anndata.py \\
        --abundance ${abundance} \\
        --annotation ${annotation} \\
        --log2transformed ${log2transformed} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tmt_None.h5mu
    touch ${prefix}.tmt_protein_None.h5ad
    printf 'log2transformed\\t${log2transformed}\\nprefix\\t${prefix}\\n' > ${prefix}.tmt.summary.tsv
    """
}
