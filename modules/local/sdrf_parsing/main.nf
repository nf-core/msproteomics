process SDRF_PARSING {
    tag "$sdrf.Name"
    label 'process_tiny'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sdrf-pipelines:0.0.33--pyhdfd78af_0' :
        'biocontainers/sdrf-pipelines:0.0.33--pyhdfd78af_0' }"

    input:
    path sdrf

    output:
    path "${sdrf.baseName}_openms_design.tsv", emit: ch_expdesign
    path "${sdrf.baseName}_config.tsv"       , emit: ch_sdrf_config_file
    path "*.log"                             , emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def extensionconversions
    if (task.ext.convert_dotd) {
        extensionconversions = ",.d.gz:.mzML,.d.tar.gz:.mzML,d.tar:.mzML,.d.zip:.mzML,.d:.mzML"
    } else {
        extensionconversions = ",.gz:,.tar.gz:,.tar:,.zip:"
    }

    """
    ## -t2 since the one-table format parser is broken in OpenMS2.5
    ## -l for legacy behavior to always add sample columns

    parse_sdrf convert-openms \\
        -t2 -l \\
        --extension_convert raw:mzML$extensionconversions \\
        -s ${sdrf} \\
        $args \\
        2>&1 | tee ${sdrf.baseName}_parsing.log

    mv openms.tsv ${sdrf.baseName}_config.tsv
    mv experimental_design.tsv ${sdrf.baseName}_openms_design.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sdrf-pipelines: \$(parse_sdrf --version 2>/dev/null | awk -F ' ' '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    echo -e "Fraction_Group\\tFraction\\tSpectra_Filepath\\tLabel\\tSample" > ${sdrf.baseName}_openms_design.tsv

    # Create config TSV with columns expected by create_input_channel
    HEADER="Filename\\tURI\\tProteomics Data Acquisition Method\\tDissociationMethod\\tLabel\\tFixedModifications\\tVariableModifications\\tPrecursorMassTolerance\\tPrecursorMassToleranceUnit\\tFragmentMassTolerance\\tFragmentMassToleranceUnit\\tEnzyme"
    echo -e "\$HEADER" > ${sdrf.baseName}_config.tsv
    # Extract data file names from SDRF (skip header), create config rows
    tail -n +2 ${sdrf} | while IFS=\$'\\t' read -r line; do
        FN=\$(echo "\$line" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i!="") fn=\$i; print fn}' | head -1 || true)
        # Use column 11 (comment[data file]) from standard SDRF
        FN=\$(echo "\$line" | cut -f11)
        FN=\${FN/.raw/.mzML}
        echo -e "\${FN}\\t\${FN}\\tData-independent acquisition\\tHCD\\tlabel free sample\\t\\t\\t20\\tppm\\t20\\tppm\\tTrypsin" >> ${sdrf.baseName}_config.tsv
    done

    echo "Stub execution" > ${sdrf.baseName}_parsing.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sdrf-pipelines: stub
    END_VERSIONS
    """
}
