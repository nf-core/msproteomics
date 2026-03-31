process DIATRACER {
    tag "${meta.id}"
    label 'process_high'

    input:
    // diaTracer processes Bruker .d files for DIA analysis
    // d_file: Bruker .d directory
    // config_cli: diaTracer parameters as CLI string (coupled with sample inputs)
    tuple val(meta), path(d_file), val(config_cli)
    path diatracer_dir // Unzipped diaTracer tool directory (optional). Pass [] when not using.

    output:
    tuple val(meta), path("${prefix}")                          , emit: results_dir
    tuple val(meta), path("${prefix}/*_diatracer.mzML")         , emit: mzml
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT", emit: license_agreement
    tuple val("${task.process}"), val('diatracer'), eval("cat .diatracer_version"), emit: versions_diatracer, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def basename = d_file.baseName.replaceAll(/\.d$/, '')
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    def agree_license = task.ext.agree_fragpipe_license_agreement ?: false

    // diaTracer CLI parameters:
    // --dFilePath: input .d file
    // --workDir: output directory
    // --threadNum: number of threads
    // --writeInter: write intermediate files (0/1)
    // --deltaApexIM: ion mobility tolerance
    // --deltaApexRT: retention time tolerance (scans)
    // --massDefectFilter: enable mass defect filter (0/1)
    // --massDefectOffset: mass defect offset
    // --ms1MS2Corr: MS1-MS2 correlation threshold
    // --RFMax: max RF value

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    # License gate
    if [ "${agree_license}" != "true" ]; then
        echo "ERROR: You must agree to the FragPipe license agreement before using this tool." >&2
        echo "Set ext.agree_fragpipe_license_agreement = true in your Nextflow config." >&2
        exit 1
    fi
    echo "INFO: User has set agree_fragpipe_license_agreement = true. Proceeding with licensed tool." >&2
    echo "User agreed to FragPipe license agreement (agree_fragpipe_license_agreement=true)" > I_AGREE_FRAGPIPE_LICENSE_AGREEMENT

    TOOLS_DIR="${tools_dir}"

    # JAR discovery: prefer provided dir, fall back to tools_dir
    if [ -d "${diatracer_dir}" ] && [ "\$(ls -A ${diatracer_dir} 2>/dev/null)" ]; then
        DIATRACER_JAR=\$(find ${diatracer_dir} -maxdepth 1 -name "diaTracer*.jar" -type f 2>/dev/null | head -1 || true)
    else
        DIATRACER_JAR=\$(find "\$TOOLS_DIR" -name 'diaTracer*.jar' -type f 2>/dev/null | head -1 || true)
    fi

    # Version capture
    _dt_ver=\$(java -jar "\$DIATRACER_JAR" --version 2>&1 | grep -oP 'diaTracer-\\K[\\d.]+' | head -1 || true)

    echo "\$_dt_ver" > .diatracer_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    BRUKER_DIR="${tools_dir}/../ext/bruker"
    NATIVE_FLAGS=""
    [ -d "\$BRUKER_DIR" ] && NATIVE_FLAGS="-Dlibs.bruker.dir=\$BRUKER_DIR"

    export JAVA_OPTS="-Xmx${ram}G"
    CMD="java \$JAVA_OPTS \$NATIVE_FLAGS -jar \\"\$DIATRACER_JAR\\" --dFilePath ${d_file} --workDir ${prefix} --threadNum ${task.cpus} ${config_cli} ${args}"
    printf '%s\\n' "\$CMD" > ${prefix}/diatracer.log
    eval "\$CMD" 2>&1 | tee -a ${prefix}/diatracer.log

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def basename = d_file.baseName.replaceAll(/\.d$/, '')
    """
    mkdir -p ${prefix}
    touch ${prefix}/${basename}_diatracer.mzML
    touch ${prefix}/diatracer.log
    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .diatracer_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
