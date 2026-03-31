process PHILOSOPHER_PIPELINE {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    tuple val(meta), path(pepxml_files), path(protxml), path(mzml_files), val(filter_flags), val(report_flags)
    path fasta
    // fasta: Philosopher database with decoys/contaminants from PHILOSOPHER_DATABASE
    // filter_flags: e.g., "--sequential --razor --prot 0.01 --picked"
    // report_flags: e.g., "--msstats"

    output:
    tuple val(meta), path("${prefix}")               , emit: results_dir
    tuple val(meta), path("${prefix}/protein.tsv")   , emit: proteins
    tuple val(meta), path("${prefix}/peptide.tsv")   , emit: peptides
    tuple val(meta), path("${prefix}/psm.tsv")       , emit: psms
    tuple val("${task.process}"), val('philosopher'), eval("cat .philosopher_version"), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def decoy_tag = task.ext.decoy_tag ?: 'rev_'
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    // Config passed as CLI strings directly (required inputs, no defaults)
    def filter_cli = filter_flags
    def report_cli = report_flags
    tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    """
    export HOME=\$(pwd)
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

    export GOMEMLIMIT="${ram}GiB"
    export GOGC=200
    export GOMAXPROCS=${task.cpus}

    # Initialize workspace and annotate database
    \$PHILOSOPHER workspace --init --nocheck 2>&1
    \$PHILOSOPHER database --annotate ${fasta} --prefix ${decoy_tag} 2>&1 | tee ${prefix}/philosopher.log

    FILTER_CMD="\$PHILOSOPHER filter --tag ${decoy_tag} --pepxml . --protxml ${protxml} ${filter_cli} ${args}"
    printf '\\n%s\\n' "\$FILTER_CMD" >> ${prefix}/philosopher.log
    \$FILTER_CMD 2>&1 | tee -a ${prefix}/philosopher.log

    REPORT_CMD="\$PHILOSOPHER report ${report_cli}"
    printf '\\n%s\\n' "\$REPORT_CMD" >> ${prefix}/philosopher.log
    \$REPORT_CMD 2>&1 | tee -a ${prefix}/philosopher.log

    shopt -s nullglob
    for f in *.tsv; do mv "\$f" ${prefix}/; done
    shopt -u nullglob
    \$PHILOSOPHER workspace --clean --nocheck 2>&1

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/protein.tsv
    touch ${prefix}/peptide.tsv
    touch ${prefix}/psm.tsv
    touch ${prefix}/philosopher.log
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
