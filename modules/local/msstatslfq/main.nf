/*
 * MSSTATSLFQ: MSstats label-free statistics on an MSstats input table (bin/msstats_plfq.R)
 *
 * Input:
 *   - msstats_csv_input: MSstats input CSV
 *   - contrasts / ref_condition: contrasts ('pairwise' or a list) and reference condition
 *   - remove_one_feat_prot / remove_few_measurements / feature_subset_protein / quant_summary_method / threshold: MSstats options
 *
 * Output:
 *   - pdf: MSstats plots (optional)
 *   - msstats_csv: result tables
 *   - log: msstats.log
 *   - versions: versions.yml (topic: versions)
 */
process MSSTATSLFQ {
    tag "${msstats_csv_input.name}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-msstats:4.14.0--r44he5774e6_0' :
        'quay.io/biocontainers/bioconductor-msstats:4.14.0--r44he5774e6_0' }"

    input:
    path msstats_csv_input
    val contrasts                // "pairwise" or an explicit contrast list
    val ref_condition            // reference (control) condition, '' for none
    val remove_one_feat_prot     // boolean: removeOneFeatProts
    val remove_few_measurements  // boolean: removeFewMeasurements
    val feature_subset_protein   // featureSubset ("top3", "all", ...)
    val quant_summary_method     // summaryMethod ("TMP", "linear")
    val threshold                // adjusted p-value threshold

    output:
    // The generation of the PDFs from MSstats are very unstable, especially with auto-contrasts.
    // And users can easily fix anything based on the csv and the included script -> make optional
    path "*.pdf", emit: pdf, optional: true
    path "*.csv", emit: msstats_csv
    path "*.log", emit: log
    path "versions.yml" , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    msstats_plfq.R \\
        ${msstats_csv_input} \\
        "${contrasts}" \\
        "${ref_condition}" \\
        ${remove_one_feat_prot} \\
        ${remove_few_measurements} \\
        ${feature_subset_protein} \\
        ${quant_summary_method} \\
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
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-msstats: \$(Rscript -e "library(MSstats); cat(as.character(packageVersion('MSstats')))" 2>/dev/null)
    END_VERSIONS
    """
}
