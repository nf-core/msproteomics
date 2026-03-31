process MSFRAGGER {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(mzml_files)
    path fasta
    path params_file  // fragger.params file (native MSFragger format). Pass [] when not using.
    path pepindex      // Prebuilt pepindex files (staged alongside FASTA). Pass [] when not using.
    path msfragger_dir // Unzipped MSFragger tool directory (optional). Pass [] when not using.

    output:
    // Search outputs
    tuple val(meta), path("${prefix}/*.pepXML")              , emit: pepxml, optional: true
    tuple val(meta), path("${prefix}/*.pin")                 , emit: pin, optional: true
    tuple val(meta), path("${prefix}/*.tsv")                 , emit: tsv, optional: true

    // Calibrated mzML (when write_calibrated_mzml=true)
    tuple val(meta), path("*_calibrated.mzML")               , emit: calibrated_mzml, optional: true

    // Uncalibrated mzML (produced when processing .d or .raw files)
    tuple val(meta), path("*_uncalibrated.mzML")             , emit: uncalibrated_mzml, optional: true

    // Results directory and logs
    tuple val(meta), path("${prefix}")                       , emit: results_dir
    tuple val(meta), path("${prefix}/*.log")                 , emit: log, optional: true

    // License agreement
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT", emit: license_agreement

    // Version tracking
    tuple val("${task.process}"), val('msfragger'), eval("cat .msfragger_version"), emit: versions_msfragger, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def mem = task.ext.java_xmx ? task.ext.java_xmx as int : (task.memory ? task.memory.toGiga() : 8)
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    def agree_license = task.ext.agree_fragpipe_license_agreement ?: false
    // Only override calibrate_mass/write_calibrated_mzml when ext is explicitly set.
    // Otherwise, the params file value (from workflow file) is preserved.
    def calibrate_mass = task.ext.calibrate_mass != null ? task.ext.calibrate_mass : ''
    def write_calibrated = task.ext.write_calibrated_mzml != null ? (task.ext.write_calibrated_mzml ? '1' : '0') : ''
    // CLI mode defaults (when no params file and no ext override)
    def calibrate_mass_cli = calibrate_mass ?: '2'
    def write_calibrated_cli = write_calibrated ?: '0'
    def mzml_input = mzml_files instanceof List ? mzml_files.join(' ') : mzml_files

    def partial_flag = meta.chunk_id != null ? "--partial ${meta.chunk_id}" : ''

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    export JAVA_OPTS="-Xmx${mem}G"

    # License gate
    if [ "${agree_license}" != "true" ]; then
        echo "ERROR: You must agree to the FragPipe license agreement before using this tool." >&2
        echo "Set ext.agree_fragpipe_license_agreement = true in your Nextflow config." >&2
        exit 1
    fi
    echo "INFO: User has set agree_fragpipe_license_agreement = true. Proceeding with licensed tool." >&2
    echo "User agreed to FragPipe license agreement (agree_fragpipe_license_agreement=true)" > I_AGREE_FRAGPIPE_LICENSE_AGREEMENT

    TOOLS_DIR="${tools_dir}"

    # JAR + native lib discovery: prefer provided dir, fall back to tools_dir
    if [ -d "${msfragger_dir}" ] && [ "\$(ls -A ${msfragger_dir} 2>/dev/null)" ]; then
        MSFRAGGER_JAR=\$(find ${msfragger_dir} -maxdepth 1 -name "MSFragger*.jar" -not -name "*-sources*" -type f 2>/dev/null | head -1 || true)
        BRUKER_DIR="${msfragger_dir}/ext/bruker"
        THERMO_DIR="${msfragger_dir}/ext/thermo"
    else
        MSFRAGGER_JAR=\$(find "\$TOOLS_DIR" -name 'MSFragger*.jar' -not -name '*-sources*' -type f 2>/dev/null | head -1 || true)
        BRUKER_DIR="${tools_dir}/../ext/bruker"
        THERMO_DIR="${tools_dir}/../ext/thermo"
    fi
    NATIVE_FLAGS=""
    [ -d "\$BRUKER_DIR" ] && NATIVE_FLAGS="\$NATIVE_FLAGS -Dlibs.bruker.dir=\$BRUKER_DIR"
    [ -d "\$THERMO_DIR" ] && NATIVE_FLAGS="\$NATIVE_FLAGS -Dlibs.thermo.dir=\$THERMO_DIR"

    # Version capture
    _ms_ver=\$(java -jar "\$MSFRAGGER_JAR" --help 2>&1 | grep -oP 'MSFragger-\\K[\\d.]+' | head -1 || true)

    echo "\$_ms_ver" > .msfragger_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    if [ -s "${params_file}" ]; then
        # Params file mode: matches FragPipe calling convention
        # Copy params file and override runtime parameters
        cp ${params_file} ${prefix}/fragger.params
        sed -i '/^database_name/d' ${prefix}/fragger.params
        sed -i '/^num_threads/d' ${prefix}/fragger.params
        sed -i '/^write_mzbin_all/d' ${prefix}/fragger.params
        echo "database_name = \$(pwd)/${fasta}" >> ${prefix}/fragger.params
        echo "num_threads = ${task.cpus}" >> ${prefix}/fragger.params

        # Only override calibrate_mass if explicitly set via ext
        if [ -n "${calibrate_mass}" ]; then
            sed -i '/^calibrate_mass/d' ${prefix}/fragger.params
            echo "calibrate_mass = ${calibrate_mass}" >> ${prefix}/fragger.params
        fi

        # Only override write_calibrated_mzml if explicitly set via ext
        if [ -n "${write_calibrated}" ]; then
            sed -i '/^write_calibrated_mzml/d' ${prefix}/fragger.params
            echo "write_calibrated_mzml = ${write_calibrated}" >> ${prefix}/fragger.params
        fi

        CMD="java \$JAVA_OPTS -Dfile.encoding=UTF-8 \$NATIVE_FLAGS -jar \\"\$MSFRAGGER_JAR\\" ${prefix}/fragger.params ${mzml_input} ${partial_flag} ${args}"
    else
        # CLI mode (backward compatible): ext.args provides all parameters
        CMD="java \$JAVA_OPTS -Dfile.encoding=UTF-8 \$NATIVE_FLAGS -jar \\"\$MSFRAGGER_JAR\\" --database_name \$(pwd)/${fasta} --num_threads ${task.cpus} --calibrate_mass ${calibrate_mass_cli} --write_calibrated_mzml ${write_calibrated_cli} ${mzml_input} ${partial_flag} ${args}"
    fi

    printf '%s\\n' "\$CMD" > ${prefix}/msfragger.log
    eval "\$CMD" 2>&1 | tee -a ${prefix}/msfragger.log

    shopt -s nullglob
    for f in *.pepXML *.pin *.tsv; do mv "\$f" "${prefix}/"; done
    shopt -u nullglob

    # Rename outputs with chunk prefix for split mode to avoid file name collisions
    # when Nextflow stages files from multiple chunks into MERGE_SPLIT_SEARCH
    if [ -n "${partial_flag}" ]; then
        cd "${prefix}"
        shopt -s nullglob
        for f in *.pepXML; do
            [ -f "\$f" ] && mv "\$f" "chunk_${meta.chunk_id}_\${f}"
        done
        for f in *.pin; do
            [ -f "\$f" ] && mv "\$f" "chunk_${meta.chunk_id}_\${f}"
        done
        for f in *_scores_histogram.tsv; do
            [ -f "\$f" ] && mv "\$f" "chunk_${meta.chunk_id}_\${f}"
        done
        shopt -u nullglob
        cd ..
    fi

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def write_calibrated = task.ext.write_calibrated_mzml ?: false
    def mzml_input = mzml_files instanceof List ? mzml_files.join(' ') : mzml_files
    def partial_flag = meta.chunk_id != null ? "--partial ${meta.chunk_id}" : ''

    """
    mkdir -p ${prefix}

    # Create per-file outputs for all input files (supports aggregate mode)
    for f in ${mzml_input}; do
        [[ ! -e "\$f" ]] && touch "\$f"
        fname=\$(basename "\$f")
        stem="\${fname%.*}"
        # In split mode, prefix outputs with chunk_id to avoid name collisions
        if [ -n "${partial_flag}" ]; then
            chunk_prefix="chunk_${meta.chunk_id}_"
        else
            chunk_prefix=""
        fi
        touch "${prefix}/\${chunk_prefix}\${stem}.pepXML"
        touch "${prefix}/\${chunk_prefix}\${stem}.pin"
        touch "${prefix}/\${chunk_prefix}\${stem}_scores_histogram.tsv"
        # MSFragger produces _uncalibrated.mzML for .d and .raw inputs
        lower=\$(echo "\$fname" | tr '[:upper:]' '[:lower:]')
        if [[ "\$lower" == *.d ]] || [[ "\$lower" == *.raw ]]; then
            touch "\${stem}_uncalibrated.mzML"
        fi
        if [[ "${write_calibrated}" == "true" ]]; then
            touch "\${stem}_calibrated.mzML"
        fi
    done
    touch ${prefix}/msfragger.log
    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .msfragger_version
    echo "\$_fp_ver" > .fragpipe_version
    """
}
