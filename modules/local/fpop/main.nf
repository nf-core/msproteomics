process FPOP {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // FPOP (Fast Photochemical Oxidation of Proteins) analysis
    // input_file: combined_modified_peptide.tsv (LFQ) or TMT abundance files
    // secondary_file: optional secondary TMT file for TMT mode
    // config_cli: bash-sourceable KEY='value' lines (REGION_SIZE, CONTROL_LABEL,
    //             FPOP_LABEL, SUBTRACT_CONTROL, IS_TMT)
    tuple val(meta), path(input_file), path(secondary_file), val(config_cli)

    output:
    tuple val(meta), path("${prefix}")                      , emit: results_dir
    tuple val(meta), path("${prefix}/*_fpop*.tsv")          , emit: results, optional: true
    tuple val(meta), path("${prefix}/*_fpop*.csv")          , emit: results_csv, optional: true
    tuple val("${task.process}"), val('fpop'), eval("cat .fpop_version"), emit: versions_fpop, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    def secondary_arg = secondary_file.name != 'NO_FILE' ? secondary_file : ''

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    # Source algorithm config (sets REGION_SIZE, CONTROL_LABEL, FPOP_LABEL,
    # SUBTRACT_CONTROL, IS_TMT)
    cat > .fpop.env << 'FPOP_CONFIG'
${config_cli}
FPOP_CONFIG
    . .fpop.env

    TOOLS_DIR="${tools_dir}"

    # Discover FPOP script
    FPOP_SCRIPT="${tools_dir}/fpop/FragPipe_FPOP_Analysis.py"
    [ ! -f "\$FPOP_SCRIPT" ] && echo "WARNING: FPOP script not found at \$FPOP_SCRIPT" >&2
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fpop_version

    echo "\$_fp_ver" > .fragpipe_version

    # Run FPOP analysis Python script (positional args)
    CMD="python3 \$FPOP_SCRIPT ${input_file} \$REGION_SIZE \$CONTROL_LABEL \$FPOP_LABEL \$SUBTRACT_CONTROL \$IS_TMT"
    if [[ -n "${secondary_arg}" && "\$IS_TMT" == "true" ]]; then
        CMD="\$CMD ${secondary_file}"
    fi
    printf '%s\\n' "\$CMD" > ${prefix}/fpop.log

    python3 \$FPOP_SCRIPT \\
        ${input_file} \\
        \$REGION_SIZE \\
        \$CONTROL_LABEL \\
        \$FPOP_LABEL \\
        \$SUBTRACT_CONTROL \\
        \$IS_TMT \\
        ${secondary_arg ? secondary_file : ''} \\
        ${args} \\
        2>&1 | tee -a ${prefix}/fpop.log

    # Move output files to results directory
    shopt -s nullglob
    for f in *_fpop*.tsv *_fpop*.csv; do mv "\$f" ${prefix}/; done
    shopt -u nullglob

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/sample_fpop_results.tsv
    touch ${prefix}/sample_fpop_results.csv
    touch ${prefix}/fpop.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fpop_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
