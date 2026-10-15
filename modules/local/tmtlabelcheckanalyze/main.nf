include { tmtLabelcheckThresholds } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process TMTLABELCHECKANALYZE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seaborn:0.13.2' :
        'quay.io/biocontainers/seaborn:0.13.2' }"

    input:
    // input_files: 'psm' mode, every group's psm.tsv in the subworkflow's order; 'ionquant' mode, combined_modified_peptide.tsv
    tuple val(meta), path(input_files, stageAs: 'input?/*')
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
    def thresholds = tmtLabelcheckThresholds(task)
    def pass_threshold = thresholds.pass
    def warn_threshold = thresholds.warn
    def prob_threshold = thresholds.prob
    if (mode == 'psm') {
        // tmt_qc.py concatenates the groups' psm.tsv into combined_psm.tsv (header once)
        """
        tmt_qc.py analyze \\
            --psm-file ${[input_files].flatten().join(' ')} \\
            --tmt-type ${tmt_type} \\
            --prob-threshold ${prob_threshold} \\
            --pass-threshold ${pass_threshold} \\
            --warn-threshold ${warn_threshold} \\
            ${args}
        """
    } else {
        """
        tmt_qc_ionquant.py \\
            --combined-modified-peptide ${input_files} \\
            --tmt-type ${tmt_type} \\
            --pass-threshold ${pass_threshold} \\
            --warn-threshold ${warn_threshold} \\
            --task-process ${task.process} \\
            ${args}
        """
    }

    stub:
    """
    echo "<html><body>STUB</body></html>" > report.html
    echo "STUB" > report.md
    echo -e "metric\\tvalue\\nqc_status\\tSTUB" > labeling_summary.tsv
    echo "sample,psm_count" > per_sample_efficiency.csv
    """
}
