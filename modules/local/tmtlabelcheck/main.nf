process TMT_LABELCHECK_ANALYZE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_pandas:76f3c63ec67531f0' :
        'community.wave.seqera.io/library/matplotlib_pandas:76f3c63ec67531f0' }"

    input:
    tuple val(meta), path(input_files)
    val tmt_type
    val mode  // 'psm' for per-sample PSM files, 'ionquant' for combined_modified_peptide.tsv

    output:
    tuple val(meta), path("report.html")              , emit: html_report
    tuple val(meta), path("report.md")                , emit: md_report
    tuple val(meta), path("labeling_summary.tsv")     , emit: summary
    tuple val(meta), path("per_sample_efficiency.csv"), emit: per_sample
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def pass_threshold = task.ext.pass_threshold ?: 0.95
    def warn_threshold = task.ext.warn_threshold ?: 0.85
    def prob_threshold = task.ext.prob_threshold ?: 0.95
    if (mode == 'psm') {
        // input_files are sample result directories (e.g., sample1/, sample2/)
        // each containing psm.tsv. Find and concatenate them.
        """
        # Find all psm.tsv files in input directories
        PSM_FILES=\$(find -L ${input_files instanceof List ? input_files.join(' ') : input_files} -name 'psm.tsv' -type f | sort)
        FIRST=\$(echo "\$PSM_FILES" | head -1 || true)
        head -1 "\$FIRST" > combined_psm.tsv
        for f in \$PSM_FILES; do
            tail -n +2 "\$f" >> combined_psm.tsv
        done
        tmt_qc.py analyze \\
            --psm-file combined_psm.tsv \\
            --tmt-type ${tmt_type} \\
            --prob-threshold ${prob_threshold} \\
            --pass-threshold ${pass_threshold} \\
            --warn-threshold ${warn_threshold}
        """
    } else {
        """
        tmt_qc_ionquant.py \\
            --combined-modified-peptide ${input_files} \\
            --tmt-type ${tmt_type} \\
            --pass-threshold ${pass_threshold} \\
            --warn-threshold ${warn_threshold} \\
            --task-process ${task.process}
        """
    }

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "<html><body>STUB</body></html>" > report.html
    echo "STUB" > report.md
    echo -e "metric\\tvalue\\nqc_status\\tSTUB" > labeling_summary.tsv
    echo "sample,psm_count" > per_sample_efficiency.csv
    """
}
