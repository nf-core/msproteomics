process PHILOSOPHER_LABELQUANT {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // Philosopher labelquant extracts TMT/isobaric reporter ion intensities from mzML spectra.
    // results_dir: directory containing psm.tsv and .meta/ workspace from PHILOSOPHER_FILTER or FREEQUANT
    // mzml_dir: directory containing mzML files for the experiment/group
    // annot: TMT annotation file mapping channels to samples
    // config_cli: labelquant CLI flags (--tol, --level, --plex, --brand, --minprob, --purity, --removelow)
    // NOTE: Philosopher labelquant runs INSIDE the results_dir workspace (CmdLabelquant.java:159).
    //       It reads .meta/ binaries and mzML files to extract reporter ion intensities.
    //       Downstream PHILOSOPHER_REPORT re-generates TSVs with TMT intensity columns.
    tuple val(meta), path(results_dir), path(mzml_dir), path(annot), val(config_cli)

    output:
    tuple val(meta), path("${prefix}")              , emit: results_dir
    tuple val("${task.process}"), val('philosopher'), eval("cat .philosopher_version"), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    """
    export HOME=\$(pwd)
    WORK_DIR=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    discover_tool() {
        local name="\$1"; shift
        local matches count
        matches=\$(find "\$TOOLS_DIR" "\$@" -type f 2>/dev/null)
        count=\$(echo "\$matches" | grep -c . 2>/dev/null)
        if [ "\$count" -eq 0 ]; then
            echo "ERROR: \$name not found in \$TOOLS_DIR" >&2; exit 1
        elif [ "\$count" -gt 1 ]; then
            echo "ERROR: Multiple \$name matches in \$TOOLS_DIR:" >&2
            echo "\$matches" >&2; exit 1
        fi
        echo "\$matches"
    }

    # Tool discovery
    PHILOSOPHER=\$(discover_tool 'philosopher' -path '*/Philosopher/philosopher*')

    # Version capture
    _ver=\$(\$PHILOSOPHER version 2>&1 | grep -oP 'version=v\\K[\\d.]+' | head -1 || true)

    echo "\$_ver" > .philosopher_version
    _ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_ver" > .fragpipe_version

    export GOMAXPROCS=${task.cpus}
    export GOMEMLIMIT="${ram}GiB"
    export GOGC=200

    # Enter results directory which contains .meta/ workspace from filter/freequant
    cd ${results_dir}
    \$PHILOSOPHER workspace --init --nocheck 2>&1

    # Run labelquant - extracts TMT reporter ion intensities from mzML
    # (CmdLabelquant.java:94-161: runs in group working directory with --annot and --dir flags)
    CMD="\$PHILOSOPHER labelquant ${config_cli} --annot \${WORK_DIR}/${annot} --dir \${WORK_DIR}/${mzml_dir} ${args}"
    printf '%s\\n' "\$CMD" > \${WORK_DIR}/${prefix}/labelquant.log
    \$CMD 2>&1 | tee -a \${WORK_DIR}/${prefix}/labelquant.log

    # Don't clean workspace - downstream PHILOSOPHER_REPORT needs .meta/
    cd \${WORK_DIR}

    # Copy results to output directory (preserve .meta/ for downstream report)
    cp -r ${results_dir}/. ${prefix}/

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/psm.tsv
    touch ${prefix}/peptide.tsv
    touch ${prefix}/protein.tsv
    touch ${prefix}/labelquant.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    PHILOSOPHER=\$(find "\$TOOLS_DIR" -name 'philosopher' -type f 2>/dev/null | head -1) || true
    if [ -z "\$PHILOSOPHER" ]; then
        PHILOSOPHER=\$(which philosopher 2>/dev/null || true)
    fi
    _ver=""
    if [ -n "\$PHILOSOPHER" ]; then
        _ver=\$("\$PHILOSOPHER" version 2>&1 | grep -oP 'version=v\\K[\\d.]+' | head -1) || true
    fi

    echo "\$_ver" > .philosopher_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1) || true

    echo "\$_fp_ver" > .fragpipe_version


    """
}
