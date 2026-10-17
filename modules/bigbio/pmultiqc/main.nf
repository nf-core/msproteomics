process PMULTIQC {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // pmultiqc is published to bioconda/biocontainers on release; a single image
    // serves Docker (native) and Singularity (via the galaxy depot mirror).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pmultiqc:0.0.48--pyhdfd78af_0' :
        'biocontainers/pmultiqc:0.0.48--pyhdfd78af_0' }"

    input:
    // Everything the report needs, staged flat under `results/`. Consumers collect
    // their inputs (config, expdesign, quant tables, qpx dataset, logs, ...) into a
    // single channel; all pmultiqc/plugin flags are supplied via `task.ext.args`.
    path 'results/*'

    output:
    path "*.html",       emit: ch_pmultiqc_report
    path "*.db",         optional: true, emit: ch_pmultiqc_db
    path "versions.yml", emit: versions, topic: versions
    path "*_data",       emit: data

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    multiqc \\
        -f \\
        ${args} \\
        ./results \\
        -o . 2>&1 | tee pmultiqc.log
    # pmultiqc 0.0.48 logs these and exits 0 with an empty report: a plugin error (pmultiqc/modules/core/core.py:107),
    # no plugin (core.py:113), and a plugin that found nothing to report (fragpipe/fragpipe.py:135,155; diann/diann.py:102
    # for --diann-plugin). The quantms plugin has no such exit: its get_data() always returns True (quantms/quantms.py:191)
    if grep -F \\
        -e "encountered an error and was skipped" \\
        -e "No pmultiqc plugin selected; skipping." \\
        -e "No FragPipe files found." \\
        -e "Required input not found: psm.tsv" \\
        -e "DIANN report not found. Please check your data!" \\
        pmultiqc.log; then
        echo "pmultiqc: the plugin failed or found nothing to report (see the line above)" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pmultiqc: \$(multiqc --pmultiqc-version | sed -e "s/pmultiqc, version //g")
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_report_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pmultiqc: \$(multiqc --pmultiqc-version | sed -e "s/pmultiqc, version //g" 2>/dev/null || echo "0.0.48")
    END_VERSIONS
    """
}
