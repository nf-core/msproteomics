process PMULTIQC {
    tag "pmultiqc"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pmultiqc:0.0.39--pyhdfd78af_0' :
        'biocontainers/pmultiqc:0.0.39--pyhdfd78af_0' }"

    input:
    path 'results/*'
    path quantms_log

    output:
    path "*.html", emit: ch_pmultiqc_report
    path "*.db", optional: true, emit: ch_pmultiqc_db
    path "versions.yml", emit: versions
    path "*_data", emit: data

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def enable_pmultiqc = task.ext.enable_pmultiqc != null ? task.ext.enable_pmultiqc : true
    def export_mztab = task.ext.export_mztab != null ? task.ext.export_mztab : true
    def skip_table = task.ext.skip_table_plots != null ? task.ext.skip_table_plots : false
    def idxml_skip = task.ext.pmultiqc_idxml_skip != null ? task.ext.pmultiqc_idxml_skip : true
    def contaminant_str = task.ext.contaminant_string ?: "CONT"
    def quant_method = task.ext.quantification_method ?: "feature_intensity"
    def disable_pmultiqc = (enable_pmultiqc) && (export_mztab) ? "--quantms_plugin" : ""
    def disable_table_plots = (enable_pmultiqc) && (skip_table) ? "--disable_table" : ""
    def disable_idxml_index = (enable_pmultiqc) && (idxml_skip) ? "--ignored_idxml" : ""
    def contaminant_affix = contaminant_str ? "--contaminant_affix ${contaminant_str}" : ""

    """
    multiqc \\
        -f \\
        ${disable_pmultiqc} \\
        --config ./results/multiqc_config.yml \\
        ${args} \\
        ${disable_table_plots} \\
        ${disable_idxml_index} \\
        ${contaminant_affix} \\
        --quantification_method ${quant_method} \\
        ./results \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pmultiqc: \$(multiqc --pmultiqc_version 2>&1 | sed -e "s/pmultiqc, version //g")
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_report_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pmultiqc: stub
    END_VERSIONS
    """
}
