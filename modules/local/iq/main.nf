process IQ {
    tag "iq"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://dongzehe/iq:0.0.1' :
        'docker.io/dongzehe/iq:0.0.1' }"

    input:
        path(report_path)
        val(q)
        val(pgq)
        val(contaminant_pattern)
    output:
        path "contaminants_removed_*.tsv"        , emit: contaminants_removed_report
        path "maxlfq.tsv"                       , emit: maxlfq
        path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'iq.r'

    stub:
    """
    touch contaminants_removed_report.tsv
    touch maxlfq.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        iq: stub
    END_VERSIONS
    """
}
