process CRYSTALC {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // Crystal-C removes chimeric artifact PSMs from pepXML files
    // pepxml: search results from MSFragger
    // mzml: mass spec file for spectrum access
    // config_file: Crystal-C params file (coupled with sample inputs)
    tuple val(meta), path(pepxml_file), path(mzml_file), path(config_file)
    path(fasta)  // protein database FASTA (shared resource, required by CrystalC)

    output:
    tuple val(meta), path("${prefix}")                  , emit: results_dir
    tuple val(meta), path("${prefix}/*_c.pepXML")       , emit: pepxml_filtered
    tuple val("${task.process}"), val('crystalc'), eval("cat .crystalc_version"), emit: versions_crystalc, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    CRYSTALC_JAR=\$(find "\$TOOLS_DIR" -name '*crystal*.jar' -type f 2>/dev/null | head -1 || true)
    BATMASS_JAR=\$(find "\$TOOLS_DIR" -name 'batmass-io*.jar' -type f 2>/dev/null | head -1 || true)
    GRPPR_JAR=\$(find "\$TOOLS_DIR" -name 'grppr*.jar' -type f 2>/dev/null | head -1 || true)
    THERMO_DIR="${tools_dir}/../ext/thermo"
    _ver=\$(ls "\$TOOLS_DIR"/original-crystalc-*.jar 2>/dev/null | grep -oP 'crystalc-\\K[\\d.]+' | head -1 || true)

    echo "\$_ver" > .crystalc_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version
    NATIVE_FLAGS=""
    [ -d "\$THERMO_DIR" ] && NATIVE_FLAGS="-Dlibs.thermo.dir=\$THERMO_DIR"

    # Crystal-C requires: fasta, raw_file_location, raw_file_extension, output_location
    # Inject runtime paths into a copy of the config (config_file may be a symlink)
    # Determine file extension from mzML filename
    MZML_NAME="${mzml_file}"
    if echo "\$MZML_NAME" | grep -qi '\\.mzml\$'; then
        RAW_EXT="mzML"
    elif echo "\$MZML_NAME" | grep -qi '\\.raw\$'; then
        RAW_EXT="raw"
    else
        RAW_EXT="mzXML"
    fi
    if [ -s "${config_file}" ]; then
        cp ${config_file} crystalc_run.params

        # Validate required analytical parameters are present in config
        # (CrystalcPanel.java:67-70, CrystalcParams.java defaults)
        for required_param in precursor_charge isotope_number precursor_mass precursor_isolation_window correct_isotope_error; do
            if ! grep -q "^\${required_param}" crystalc_run.params; then
                echo "ERROR: Required parameter '\${required_param}' not found in Crystal-C config file" >&2
                exit 1
            fi
        done

        # Inject runtime paths (ensure newline before appending in case config lacks trailing newline)
        echo "" >> crystalc_run.params
        echo "fasta = ${fasta}" >> crystalc_run.params
        echo "raw_file_location = ." >> crystalc_run.params
        echo "raw_file_extension = \$RAW_EXT" >> crystalc_run.params
        echo "output_location = ." >> crystalc_run.params
        sed -i '/^thread/d' crystalc_run.params
        echo "thread = ${task.cpus}" >> crystalc_run.params
    else
        # No config provided - create minimal params with runtime params only
        # CrystalC will use compiled-in defaults for algorithm parameters
        cat > crystalc_run.params <<PARAMS
# CrystalC defaults - no algorithm config provided
fasta = ${fasta}
raw_file_location = .
raw_file_extension = \$RAW_EXT
output_location = .
thread = ${task.cpus}
PARAMS
    fi

    # Run Crystal-C
    export JAVA_OPTS="-Xmx${ram}G"
    CMD="java \$JAVA_OPTS \$NATIVE_FLAGS -cp \\"\$CRYSTALC_JAR:\$BATMASS_JAR:\$GRPPR_JAR\\" crystalc.Run crystalc_run.params ${pepxml_file} ${args}"
    printf '%s\\n' "\$CMD" > ${prefix}/crystalc.log
    java \$JAVA_OPTS \$NATIVE_FLAGS -cp "\$CRYSTALC_JAR:\$BATMASS_JAR:\$GRPPR_JAR" crystalc.Run \\
        crystalc_run.params \\
        ${pepxml_file} \\
        ${args} \\
        2>&1 | tee -a ${prefix}/crystalc.log

    # Move output files to prefix directory
    shopt -s nullglob
    for f in *_c.pepXML; do mv "\$f" "${prefix}/"; done
    shopt -u nullglob

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def basename = pepxml_file.baseName
    """
    mkdir -p ${prefix}
    touch ${prefix}/${basename}_c.pepXML
    touch ${prefix}/crystalc.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .crystalc_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
