/*
 * IQ: MaxLFQ protein quantification of a DIA-NN report with the iq R package (templates/iq.r)
 *
 * Input:
 *   - report_path: DIA-NN report.tsv
 *   - q / pgq: precursor and protein-group q-value cut-offs
 *   - contaminant_pattern: protein prefix of contaminants ('' for none)
 *
 * Output:
 *   - contaminants_removed_report: the report without contaminant rows (when a pattern is given)
 *   - maxlfq: maxlfq.tsv
 *   - versions: versions.yml (topic: versions)
 */
process IQ {
    tag "iq"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed074b4541efbe9f6ccf072705df59ae1e3e962e46e66d4b8275c4f5bfe3645f/data' :
        'community.wave.seqera.io/library/r-iq_r-readr:2434db545350cb42' }"

    input:
        path(report_path)
        val(q)
        val(pgq)
        val(contaminant_pattern)
    output:
        path "contaminants_removed_*.tsv"        , emit: contaminants_removed_report, optional: true
        path "maxlfq.tsv"                       , emit: maxlfq
        path "versions.yml"                     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'iq.r'

    stub:
    def contaminants_cmd = contaminant_pattern ? 'touch contaminants_removed_report.tsv' : ''
    """
    ${contaminants_cmd}
    touch maxlfq.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "\$(Rscript -e 'cat(R.version\$major, R.version\$minor, sep = ".")')"
        iq: "\$(Rscript -e 'cat(as.character(packageVersion("iq")))')"
    END_VERSIONS
    """
}
