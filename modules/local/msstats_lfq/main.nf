process MSSTATS_LFQ {
    tag "$msstats_csv_input.Name"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-msstats:4.14.0--r44he5774e6_0' :
        'biocontainers/bioconductor-msstats:4.14.0--r44he5774e6_0' }"

    input:
    path msstats_csv_input

    output:
    // The generation of the PDFs from MSstats are very unstable, especially with auto-contrasts.
    // And users can easily fix anything based on the csv and the included script -> make optional
    path "*.pdf", optional: true
    path "*.csv", emit: msstats_csv
    path "*.log", emit: log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ref_con = task.ext.ref_condition ?: ""
    def contrasts = task.ext.contrasts ?: "pairwise"
    def remove_one_feat = task.ext.remove_one_feat_prot != null ? task.ext.remove_one_feat_prot : true
    def remove_few = task.ext.removeFewMeasurements != null ? task.ext.removeFewMeasurements : true
    def feature_subset = task.ext.feature_subset_protein ?: "top3"
    def quant_method = task.ext.quant_summary_method ?: "TMP"
    def threshold = task.ext.msstats_threshold ?: 0.05

    """
    Rscript ${moduleDir}/resources/msstats_plfq.R \\
        ${msstats_csv_input} \\
        "${contrasts}" \\
        "${ref_con}" \\
        ${remove_one_feat} \\
        ${remove_few} \\
        ${feature_subset} \\
        ${quant_method} \\
        ${msstats_csv_input.baseName} \\
        ${threshold} \\
        $args \\
        2>&1 | tee msstats.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-msstats: \$(Rscript -e "library(MSstats); cat(as.character(packageVersion('MSstats')))" 2>/dev/null)
    END_VERSIONS
    """

    stub:
    """
    touch msstats_results.csv
    echo "Stub execution" > msstats.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        bioconductor-msstats: stub
    END_VERSIONS
    """
}
