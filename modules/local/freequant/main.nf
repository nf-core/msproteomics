process FREEQUANT {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // FreeQuant performs label-free quantification using Philosopher
    // results_dir: directory containing psm.tsv from PHILOSOPHER_FILTER
    // mzml_dir: directory containing mass spec files
    // config_cli: FreeQuant parameters as CLI string (coupled with sample inputs)
    tuple val(meta), path(results_dir), path(mzml_dir), val(config_cli)

    output:
    tuple val(meta), path("${prefix}")                          , emit: results_dir
    tuple val(meta), path("${prefix}/ion.tsv")                  , emit: ions, optional: true
    tuple val("${task.process}"), val('philosopher'), eval("cat .philosopher_version"), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    // FreeQuant (Philosopher freequant) CLI parameters:
    // --ptw: peak width in minutes
    // --tol: mass tolerance in ppm
    // --isolated: use isolated peaks only
    // --dir: directory containing LCMS files
    // --raw: use raw files (for Thermo .raw)

    """
    
    export HOME=\$(pwd)
    WORK_DIR=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    PHILOSOPHER_SRC=\$(find "\$TOOLS_DIR" -path '*/Philosopher/philosopher*' -type f 2>/dev/null | head -1 || true)
    PHILOSOPHER="\$WORK_DIR/philosopher"
    cp "\$PHILOSOPHER_SRC" "\$PHILOSOPHER"
    chmod +x "\$PHILOSOPHER"
    _ver=\$(\$PHILOSOPHER version 2>&1 | grep -oP 'version=v\\K[\\d.]+' | head -1 || true)

    echo "\$_ver" > .philosopher_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    export GOMAXPROCS=${task.cpus}
    export GOMEMLIMIT="${ram}GiB"
    export GOGC=200

    # Detect .raw files and add --raw flag (CmdFreequant.java:116-120)
    RAW_FLAG=""
    if compgen -G "${mzml_dir}/*.raw" > /dev/null; then
        RAW_FLAG="--raw"
    fi

    # Initialize Philosopher workspace
    cd ${results_dir}
    \$PHILOSOPHER workspace --init --nocheck 2>&1

    # Run FreeQuant (philosopher freequant, see PhilosopherProps.java:31)
    CMD="\$PHILOSOPHER freequant ${config_cli} \$RAW_FLAG --dir \${WORK_DIR}/${mzml_dir} ${args}"
    printf '%s\\n' "\$CMD" > \${WORK_DIR}/${prefix}/freequant.log
    eval "\$CMD" 2>&1 | tee -a \${WORK_DIR}/${prefix}/freequant.log

    # Do NOT clean workspace - preserve .meta/ binaries for downstream PHILOSOPHER_REPORT.
    # FreeQuant updates .meta/ with quantification data; report needs these to generate
    # TSVs with correct intensities.

    cd \${WORK_DIR}

    # Copy results and .meta/ workspace to output directory
    [ -f "${results_dir}/ion.tsv" ] && cp "${results_dir}/ion.tsv" ${prefix}/
    cp -a ${results_dir}/.meta ${prefix}/

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/ion.tsv
    touch ${prefix}/freequant.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .philosopher_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
