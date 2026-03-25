process SAINTEXPRESS {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // SAINTexpress performs AP-MS scoring for protein-protein interactions
    // inter_file: interaction file (inter.dat)
    // bait_file: bait definitions (bait.dat)
    // prey_file: prey definitions (prey.dat)
    // mode: "spc" (spectral count) or "int" (intensity)
    // config_cli: SAINTexpress parameters as CLI string (coupled with sample inputs)
    tuple val(meta), path(inter_file), path(bait_file), path(prey_file), val(mode), val(config_cli)

    output:
    tuple val(meta), path("${prefix}")                          , emit: results_dir
    tuple val(meta), path("${prefix}/list.txt")                 , emit: results_list
    tuple val("${task.process}"), val('saintexpress'), eval("cat .saintexpress_version"), emit: versions_saintexpress, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    def saint_suffix = mode == 'int' ? 'SAINTexpress-int' : 'SAINTexpress-spc'

    // SAINTexpress CLI parameters:
    // -R: max replicates per bait
    // -L: virtual controls
    // Input files must be in current directory: inter.dat, bait.dat, prey.dat

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    SAINTEXPRESS_SPC_SRC=\$(find "\$TOOLS_DIR" -path '*/SAINTexpress/SAINTexpress-spc' -type f 2>/dev/null | head -1 || true)
    SAINTEXPRESS_SPC="\$(pwd)/SAINTexpress-spc"
    cp "\$SAINTEXPRESS_SPC_SRC" "\$SAINTEXPRESS_SPC"
    chmod +x "\$SAINTEXPRESS_SPC"
    SAINTEXPRESS_INT_SRC=\$(find "\$TOOLS_DIR" -path '*/SAINTexpress/SAINTexpress-int' -type f 2>/dev/null | head -1 || true)
    SAINTEXPRESS_INT="\$(pwd)/SAINTexpress-int"
    cp "\$SAINTEXPRESS_INT_SRC" "\$SAINTEXPRESS_INT"
    chmod +x "\$SAINTEXPRESS_INT"
    # Version capture
    _se_ver=\$("\$SAINTEXPRESS_SPC" --version 2>&1 | grep -oP 'version \\K[\\d.]+' | head -1 || true)

    echo "\$_se_ver" > .saintexpress_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    # Select binary based on mode
    if [ "${saint_suffix}" = "SAINTexpress-int" ]; then
        SAINT_BIN="\$SAINTEXPRESS_INT"
    else
        SAINT_BIN="\$SAINTEXPRESS_SPC"
    fi

    # SAINTexpress requires input files in working directory
    cp ${inter_file} ${prefix}/inter.dat
    cp ${bait_file} ${prefix}/bait.dat
    cp ${prey_file} ${prefix}/prey.dat

    cd ${prefix}

    # Run SAINTexpress
    CMD="\$SAINT_BIN ${config_cli} ${args}"
    printf '%s\\n' "\$CMD" > saintexpress.log
    \$SAINT_BIN ${config_cli} ${args} 2>&1 | tee -a saintexpress.log

    cd ..

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/list.txt
    touch ${prefix}/saintexpress.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _se_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_se_ver" > .saintexpress_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    """
}
