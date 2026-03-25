process PHILOSOPHER_DATABASE {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}_philosopher.fasta"), emit: fasta
    tuple val(meta), path("philosopher_database.log")   , emit: log
    tuple val("${task.process}"), val('philosopher'), eval("cat .philosopher_version"), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def decoy_tag = task.ext.decoy_tag ?: 'rev_'
    def mem_bytes = task.memory ? (long)(task.memory.toBytes() * 0.9) : 4294967296L
    tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    """
    export HOME=\$(pwd)

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

    # Go runtime tuning (CmdPhilosopherDbAnnotate.java:68-70)
    export GOMEMLIMIT="${mem_bytes}"
    export GOGC=100
    export GOMAXPROCS=${task.cpus}

    # Check if FASTA already has decoy entries (e.g., from FragPipe database preparation).
    # philosopher database --custom adds decoys + contaminants to a raw FASTA.
    # If decoys already exist, --custom would create double-decoys (rev_rev_),
    # which breaks all downstream tools. Use the database as-is in that case.
    # FragPipe's CmdPhilosopherDbAnnotate.java always uses --annotate (not --custom)
    # because FragPipe handles decoy generation separately.
    \$PHILOSOPHER workspace --init --nocheck 2>&1

    if grep -qm1 "^>${decoy_tag}" ${fasta}; then
        # Database already has decoys — just annotate (CmdPhilosopherDbAnnotate.java).
        # --custom would add double-decoys (rev_rev_), breaking all downstream tools.
        CMD="\$PHILOSOPHER database --annotate ${fasta} --prefix ${decoy_tag}"
        printf '%s\\n' "\$CMD" > philosopher_database.log
        \$CMD 2>&1 | tee -a philosopher_database.log
        cp ${fasta} ${prefix}_philosopher.fasta
    else
        # Raw target-only database — add decoys and contaminants
        CMD="\$PHILOSOPHER database --custom ${fasta} ${args}"
        printf '%s\\n' "\$CMD" > philosopher_database.log
        \$CMD 2>&1 | tee -a philosopher_database.log
        mv *-${fasta}.fas ${prefix}_philosopher.fasta
    fi

    \$PHILOSOPHER workspace --clean --nocheck 2>&1

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_philosopher.fasta
    touch philosopher_database.log
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
