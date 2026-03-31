process OPAIR {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // O-Pair performs O-glycoproteomics analysis
    // psm_file: psm.tsv from PHILOSOPHER_FILTER
    // mzml_files: mass spec files
    // config_cli: O-Pair parameters as CLI string (coupled with sample inputs)
    tuple val(meta), path(psm_file), path(mzml_files, stageAs: 'spectra/*'), val(config_cli)
    path(glycan_db)  // glycan database file (shared resource)

    output:
    tuple val(meta), path("${prefix}")                              , emit: results_dir
    tuple val(meta), path("${prefix}/*_opair_results.tsv")          , emit: results, optional: true
    tuple val(meta), path("${prefix}/*_opair_glycoforms.tsv")       , emit: glycoforms, optional: true
    tuple val("${task.process}"), val('opair'), eval("cat .opair_version"), emit: versions_opair, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    // O-Pair CLI parameters:
    // -b: product PPM tolerance
    // -c: precursor PPM tolerance
    // -f: oxonium filter rules file
    // -m: min oxonium intensity
    // -g: O-glycan database file
    // -x: glycan residues file
    // -y: glycan mods file
    // -n: max number of glycans
    // -z: allowed sites (e.g., STY)
    // -t: number of threads
    // -i: min isotope
    // -j: max isotope
    // -r: file list
    // -s: psm file
    // -o: output directory

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    # Discover O-Pair DLL
    OPAIR_DLL="${tools_dir}/opair/CMD.dll"
    [ ! -f "\$OPAIR_DLL" ] && echo "WARNING: O-Pair DLL not found at \$OPAIR_DLL" >&2
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .opair_version

    echo "\$_fp_ver" > .fragpipe_version

    # O-Pair expects Glycan_Mods/glycan_residues.txt and glycan_mods.txt relative to its
    # DLL directory (GlobalVariables.cs:61). The container stores these files in
    # tools/Glycan_Databases/ instead of tools/opair/Glycan_Mods/. Rather than relying on
    # symlinks (which can fail with read-only container filesystems), pass the paths
    # explicitly via -x (glycan residues) and -y (glycan mods) CLI flags.
    GLYCAN_SRC="\$TOOLS_DIR/Glycan_Databases"
    GLYCAN_RESIDUES_FLAG=""
    GLYCAN_MODS_FLAG=""
    if [ -d "\$GLYCAN_SRC" ]; then
        [ -f "\$GLYCAN_SRC/glycan_residues.txt" ] && GLYCAN_RESIDUES_FLAG="-x \$GLYCAN_SRC/glycan_residues.txt"
        [ -f "\$GLYCAN_SRC/glycan_mods.txt" ] && GLYCAN_MODS_FLAG="-y \$GLYCAN_SRC/glycan_mods.txt"
    fi

    # Create file list for O-Pair
    find spectra -type f -name "*.mzML" > ${prefix}/filelist_opair.txt

    # Run O-Pair
    CMD="dotnet \$OPAIR_DLL -t ${task.cpus} -s ${psm_file} -r ${prefix}/filelist_opair.txt -o ${prefix} -g ${glycan_db} \$GLYCAN_RESIDUES_FLAG \$GLYCAN_MODS_FLAG ${config_cli} ${args}"
    printf '%s\\n' "\$CMD" > ${prefix}/opair.log
    dotnet "\$OPAIR_DLL" \\
        -t ${task.cpus} \\
        -s ${psm_file} \\
        -r ${prefix}/filelist_opair.txt \\
        -o ${prefix} \\
        -g ${glycan_db} \\
        \$GLYCAN_RESIDUES_FLAG \\
        \$GLYCAN_MODS_FLAG \\
        ${config_cli} \\
        ${args} \\
        2>&1 | tee -a ${prefix}/opair.log

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/sample_opair_results.tsv
    touch ${prefix}/sample_opair_glycoforms.tsv
    touch ${prefix}/opair.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .opair_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
