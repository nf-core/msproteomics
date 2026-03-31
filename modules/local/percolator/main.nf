process PERCOLATOR {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    tuple val(meta), path(pin_file), path(pepxml_files), path(mzml_file), val(percolator_cli)
    val data_type  // 'DDA' or 'DIA' - affects PercolatorOutputToPepXML output format
                   // Auto-overridden to 'DIA' when ranked pepXMLs detected (DDA+ ion mobility)

    // percolator_cli: CLI string for Percolator (e.g., "--only-psms --no-terminate --post-processing-tdc --trainFDR 0.01 --testFDR 0.01")
    // See percolator --help for all options. Common flags:
    //   --only-psms: Only report PSM-level results (no peptide/protein inference)
    //   --no-terminate: Don't stop on convergence (run all iterations)
    //   --post-processing-tdc: Use target-decoy competition for final FDR estimation
    //   --trainFDR <value>: FDR threshold for training (default 0.01)
    //   --testFDR <value>: FDR threshold for testing/output (default 0.01)

    output:
    tuple val(meta), path("${prefix}")                     , emit: results_dir
    tuple val(meta), path("${prefix}/interact-*.pep.xml")  , emit: pepxml
    tuple val(meta), path("${prefix}/*_target_psms.tsv")   , emit: target_psms
    tuple val(meta), path("${prefix}/*_decoy_psms.tsv")    , emit: decoy_psms
    tuple val("${task.process}"), val('percolator'), eval("cat .percolator_version"), emit: versions_percolator, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def basename = pin_file.baseName  // Use actual PIN filename (MSBOOSTER outputs clean names)
    def min_prob = task.ext.min_prob ?: '0.5'
    def decoy_tag = task.ext.decoy_tag ? "--protein-decoy-pattern ${task.ext.decoy_tag}" : ''
    // percolator_cli passed directly from input - no defaults, caller must specify all flags
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    PERCOLATOR_SRC=\$(find "\$TOOLS_DIR" -path '*/percolator_*/linux/percolator' -type f 2>/dev/null | head -1 || true)
    PERCOLATOR="\$(pwd)/percolator"
    cp "\$PERCOLATOR_SRC" "\$PERCOLATOR"
    chmod +x "\$PERCOLATOR"
    # FragPipe JAR is in lib/ (sibling of tools/), not inside tools/ — search directly
    FP_LIB="${tools_dir}/../lib"
    FP_JAR=\$(find "\$FP_LIB" -name 'fragpipe*.jar' -type f 2>/dev/null | head -1 || true)

    _ver=\$(\$PERCOLATOR --version 2>&1 | grep -oP 'version \\K[\\d.]+' | head -1 || true)

    echo "\$_ver" > .percolator_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    # Determine output flags based on --only-psms presence (matching FragPipe CmdPercolator.java:157-175).
    # When --only-psms is set, Percolator outputs PSM-level results (--results-psms/--decoy-results-psms).
    # Otherwise, it outputs peptide-level results (--results-peptides/--decoy-results-peptides).
    # Note: output filenames use _psms.tsv suffix regardless for downstream compatibility.
    if echo "${percolator_cli} ${args}" | grep -q -- '--only-psms'; then
        RESULTS_FLAG="--results-psms"
        DECOY_FLAG="--decoy-results-psms"
    else
        RESULTS_FLAG="--results-peptides"
        DECOY_FLAG="--decoy-results-peptides"
    fi

    # Run Percolator with user-specified flags
    # --num-threads and --protein-decoy-pattern are added by the module
    PERC_CMD="\$PERCOLATOR --num-threads ${task.cpus} ${decoy_tag} \$RESULTS_FLAG ${basename}_percolator_target_psms.tsv \$DECOY_FLAG ${basename}_percolator_decoy_psms.tsv ${percolator_cli} ${args} ${pin_file}"
    printf '%s\\n' "\$PERC_CMD" > ${prefix}/percolator.log
    \$PERC_CMD 2>&1 | tee -a ${prefix}/percolator.log

    # PercolatorOutputToPepXML expects pepXML at ${basename}.pepXML (or ${basename}_rank1.pepXML
    # for DDA+/DIA). Nextflow stages files with their original names, which may differ from
    # basename (derived from PIN filename). Create symlinks so the Java code finds them.
    for f in ${pepxml_files}; do
        if [ ! -e "${basename}.pepXML" ] && [[ "\$f" == *.pepXML ]] && [[ "\$f" != *_rank*.pepXML ]]; then
            ln -sf "\$f" "${basename}.pepXML"
        fi
        # Handle ranked pepXMLs (DDA+ ion mobility)
        if [[ "\$f" == *_rank*.pepXML ]]; then
            rank=\$(echo "\$f" | sed -n 's/.*_rank\\([0-9]*\\)\\.pepXML/\\1/p')
            if [ -n "\$rank" ] && [ ! -e "${basename}_rank\${rank}.pepXML" ]; then
                ln -sf "\$f" "${basename}_rank\${rank}.pepXML"
            fi
        fi
    done

    # Auto-detect effective data_type for perc2pepxml.
    # DDA+ (ion mobility, e.g. Bruker timsTOF .d) produces ranked pepXMLs (_rank1.pepXML etc.)
    # and requires DIA mode in perc2pepxml for per-rank processing
    # (CmdPercolator.java: getDataType().contentEquals("DDA") returns false for "DDA+").
    EFFECTIVE_DATA_TYPE="${data_type}"
    if [ -f "${basename}_rank1.pepXML" ]; then
        EFFECTIVE_DATA_TYPE="DIA"
    fi

    # Convert Percolator output to pepXML format using FragPipe's converter
    # Arguments: pin_file basename target_psms decoy_psms output_prefix data_type min_prob mzml_path
    PEPXML_CMD="java \${JAVA_OPTS:-} -cp \"\$FP_JAR:\$FP_LIB/*\" org.nesvilab.fragpipe.tools.percolator.PercolatorOutputToPepXML ${pin_file} ${basename} ${basename}_percolator_target_psms.tsv ${basename}_percolator_decoy_psms.tsv interact-${basename} \$EFFECTIVE_DATA_TYPE ${min_prob} \$(pwd)/${mzml_file}"
    printf '\\n%s\\n' "\$PEPXML_CMD" >> ${prefix}/percolator.log
    \$PEPXML_CMD 2>&1 | tee -a ${prefix}/percolator.log

    shopt -s nullglob
    for f in interact-*.pep.xml *_psms.tsv; do mv "\$f" "${prefix}/"; done
    shopt -u nullglob

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def basename = pin_file.baseName  // Use actual PIN filename
    """
    mkdir -p ${prefix}
    # Create ranked interact files if ranked pepXMLs are present (DDA+ ion mobility)
    if compgen -G "${basename}_rank*.pepXML" > /dev/null; then
        for f in ${basename}_rank*.pepXML; do
            rank=\$(echo "\$f" | sed -n 's/.*_rank\\([0-9]*\\).*/\\1/p')
            touch "${prefix}/interact-${basename}_rank\${rank}.pep.xml"
        done
    else
        touch ${prefix}/interact-${basename}.pep.xml
    fi
    touch ${prefix}/${basename}_percolator_target_psms.tsv
    touch ${prefix}/${basename}_percolator_decoy_psms.tsv
    touch ${prefix}/percolator.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .percolator_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
