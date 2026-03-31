process DIANNR {
    tag "diannr"
    label 'process_high'
    // NOTE: diann-rpackage is not available on conda; environment.yml covers R base only
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://dongzehe/diann-r:0.0.1' :
        'docker.io/dongzehe/diann-r:0.0.1' }"

    input:
        path(report_path)
        val(q)
        val(pgq)
        val(contaminant_pattern)
    output:
        path "contaminants_removed_*.tsv"       , emit: contaminants_removed_report
        path "precursors.tsv"                   , emit: precursors
        path "peptides.tsv"                     , emit: peptides
        path "peptides_maxlfq.tsv"              , emit: peptides_maxlfq
        path "unique_genes.tsv"                 , emit: unique_genes
        path "protein_groups_maxlfq.tsv"        , emit: protein_groups_maxlfq
        path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'diannr.r'

    stub:
    """
    touch contaminants_removed_report.tsv
    touch precursors.tsv
    touch peptides.tsv
    touch peptides_maxlfq.tsv
    touch unique_genes.tsv
    touch protein_groups_maxlfq.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        diann-rpackage: stub
    END_VERSIONS
    """
}
