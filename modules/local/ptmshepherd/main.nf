process PTMSHEPHERD {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // psm_dirs: directories containing psm.tsv files from PHILOSOPHER_FILTER
    // protxml: combined protein inference file from PROTEINPROPHET
    // mzml_files: mass spec files for intensity extraction
    // config_file: shepherd.config - PTMShepherd requires a config file (coupled with sample inputs)
    tuple val(meta), path(psm_dirs), path(protxml), path(mzml_files, stageAs: 'spectra/*'), path(config_file)
    path fasta  // protein database FASTA (shared resource, required by PTMShepherd config)

    output:
    tuple val(meta), path("${prefix}")                                  , emit: results_dir
    tuple val(meta), path("${prefix}/global.profile.tsv")               , emit: global_profile, optional: true
    tuple val(meta), path("${prefix}/global.modsummary.tsv")            , emit: global_modsummary, optional: true
    tuple val(meta), path("${prefix}/*diagmine.tsv")                    , emit: diagmine, optional: true
    tuple val(meta), path("${prefix}/*localization.tsv")                , emit: localization, optional: true
    tuple val(meta), path("${prefix}/*glycoprofile.tsv")                , emit: glycoprofile, optional: true
    tuple val("${task.process}"), val('ptmshepherd'), eval("cat .ptmshepherd_version"), emit: versions_ptmshepherd, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    """
    
    # Clean up files from previous attempts (Nextflow retries in same work dir,
    # and bash -C noclobber prevents overwriting existing files)
    rm -f .ptmshepherd_version .fragpipe_version shepherd.config 2>/dev/null || true
    rm -rf ${prefix} ptm-shepherd-output 2>/dev/null || true

    export HOME=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    PTMS_JAR=\$(find "\$TOOLS_DIR" -name 'ptmshepherd*.jar' -type f 2>/dev/null | head -1 || true)
    BATMASS_JAR=\$(find "\$TOOLS_DIR" -name 'batmass-io*.jar' -type f 2>/dev/null | head -1 || true)
    COMMONS_MATH_JAR=\$(find "\$TOOLS_DIR" -name 'commons-math3*.jar' -type f 2>/dev/null | head -1 || true)
    HIPP_CORE_JAR=\$(find "\$TOOLS_DIR" -path '*/hipparchus-*/hipparchus-core*.jar' -type f 2>/dev/null | head -1 || true)
    HIPP_STAT_JAR=\$(find "\$TOOLS_DIR" -path '*/hipparchus-*/hipparchus-stat*.jar' -type f 2>/dev/null | head -1 || true)
    # IonQuant is optional for PTMShepherd (only needed for some features)
    IONQUANT_JAR=\$(find "\$TOOLS_DIR" -name 'IonQuant*.jar' -type f 2>/dev/null | head -1 || true)
    THERMO_DIR="${tools_dir}/../ext/thermo"
    _ver=\$(ls "\$TOOLS_DIR"/ptmshepherd-*.jar 2>/dev/null | grep -oP 'ptmshepherd-\\K[\\d.]+' | head -1 || true)

    echo "\$_ver" >| .ptmshepherd_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" >| .fragpipe_version
    CP="\$PTMS_JAR:\$BATMASS_JAR:\$COMMONS_MATH_JAR:\$HIPP_CORE_JAR:\$HIPP_STAT_JAR"
    [ -n "\$IONQUANT_JAR" ] && CP="\$CP:\$IONQUANT_JAR"
    NATIVE_FLAGS=""
    [ -d "\$THERMO_DIR" ] && NATIVE_FLAGS="-Dlibs.thermo.dir=\$THERMO_DIR"

    if [ -s "${config_file}" ]; then
        # Config provided - copy and inject runtime params
        cp "${config_file}" shepherd.config
        sed -i '/^threads/d' shepherd.config
        sed -i '/^database /d' shepherd.config
        sed -i '/^dataset /d' shepherd.config
        # Ensure file ends with a newline before appending runtime params
        sed -i -e '\$a\\' shepherd.config
    else
        # No config provided - create minimal config with runtime params only
        # PTMShepherd will use compiled-in defaults for algorithm parameters
        cat > shepherd.config <<CONFIG
# PTMShepherd defaults - no algorithm config provided
CONFIG
    fi
    echo "threads = ${task.cpus}" >> shepherd.config

    # PTMShepherd 3.0.11+ requires msfragger_massdiff_to_varmod (no compiled-in default).
    # Inject 0 (disabled) when the config does not already set it.
    if ! grep -q '^msfragger_massdiff_to_varmod' shepherd.config; then
        echo "msfragger_massdiff_to_varmod = 0" >> shepherd.config
    fi

    # Add database path (PtmshepherdParams.java:72)
    echo "database = \$(pwd)/${fasta}" >> shepherd.config

    # Build dataset lines from per-sample PSM directories (PtmshepherdParams.java:86-90)
    # Format: dataset = <sample_name> <psm.tsv_path> <spectra_dir>
    for PSM_DIR in ${psm_dirs}; do
        SAMPLE_NAME=\$(basename "\$PSM_DIR")
        if [ -f "\$PSM_DIR/psm.tsv" ]; then
            echo "dataset = \$SAMPLE_NAME \$(pwd)/\$PSM_DIR/psm.tsv \$(pwd)/spectra" >> shepherd.config
        fi
    done

    # Run PTMShepherd
    export JAVA_OPTS="-Xmx${ram}G"
    CMD="java \$JAVA_OPTS \$NATIVE_FLAGS -cp \\"\$CP\\" edu.umich.andykong.ptmshepherd.PTMShepherd shepherd.config ${args}"
    printf '%s\\n' "\$CMD" >| ${prefix}/ptmshepherd.log
    eval "\$CMD" 2>&1 | tee -a ${prefix}/ptmshepherd.log

    # Move output files to prefix directory (all are optional depending on analysis type)
    # PTMShepherd writes outputs to ptm-shepherd-output/ subdirectory
    shopt -s nullglob
    for f in ptm-shepherd-output/global.profile.tsv ptm-shepherd-output/global.modsummary.tsv \
             ptm-shepherd-output/*diagmine.tsv ptm-shepherd-output/*localization.tsv ptm-shepherd-output/*glycoprofile.tsv \
             global.profile.tsv global.modsummary.tsv *diagmine.tsv *localization.tsv *glycoprofile.tsv; do
        [ -f "\$f" ] && mv "\$f" ${prefix}/
    done
    shopt -u nullglob

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/global.profile.tsv
    touch ${prefix}/global.modsummary.tsv
    touch ${prefix}/ptmshepherd.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" >| .ptmshepherd_version
    echo "\$_fp_ver" >| .fragpipe_version

    """
}
