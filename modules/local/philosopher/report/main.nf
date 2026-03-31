process PHILOSOPHER_REPORT {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // results_dir: directory containing .meta/ workspace binaries from PHILOSOPHER_FILTER (skip_report=true)
    //              or from FREEQUANT (which preserves .meta/ after updating psm.tsv with intensities)
    // report_cli: e.g., "--msstats" - CLI flags for philosopher report command
    tuple val(meta), path(results_dir), val(report_cli)
    path fasta
    // fasta: Philosopher database with decoys/contaminants from PHILOSOPHER_DATABASE
    // NOTE: philosopher report requires an annotated database in the workspace.
    // Since the workspace may come from a different task's work directory,
    // we re-init and re-annotate to ensure the database is accessible.

    output:
    tuple val(meta), path("${prefix}")              , emit: results_dir
    tuple val(meta), path("${prefix}/psm.tsv")      , emit: psms
    tuple val(meta), path("${prefix}/peptide.tsv")   , emit: peptides
    tuple val(meta), path("${prefix}/protein.tsv")   , emit: proteins
    tuple val(meta), path("${prefix}/ion.tsv")       , emit: ions, optional: true
    tuple val("${task.process}"), val('philosopher'), eval("cat .philosopher_version"), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def decoy_prefix = task.ext.decoy_tag ?: 'rev_'
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

    export GOMEMLIMIT="${ram}GiB"
    export GOGC=200
    export GOMAXPROCS=${task.cpus}

    # Work inside results_dir which contains .meta/ workspace binaries
    cd ${results_dir}

    # Re-initialize workspace and annotate database
    # The .meta/ binaries from filter/freequant are present, but we need
    # the database annotation accessible in this task's context
    \$PHILOSOPHER workspace --init --nocheck 2>&1
    \$PHILOSOPHER database --annotate \${WORK_DIR}/${fasta} --prefix ${decoy_prefix} 2>&1 | tee \${WORK_DIR}/${prefix}/philosopher_report.log

    # Run report to generate TSV files with correct intensities
    # After FreeQuant has updated .meta/ binaries with quantification data,
    # report reads the updated binaries and produces TSVs with intensities filled in
    REPORT_CMD="\$PHILOSOPHER report ${report_cli} ${args}"
    printf '\\n%s\\n' "\$REPORT_CMD" >> \${WORK_DIR}/${prefix}/philosopher_report.log
    \$REPORT_CMD 2>&1 | tee -a \${WORK_DIR}/${prefix}/philosopher_report.log

    # Move generated TSV files to output directory
    shopt -s nullglob
    for f in *.tsv; do mv "\$f" \${WORK_DIR}/${prefix}/; done
    shopt -u nullglob

    # Clean workspace
    \$PHILOSOPHER workspace --clean --nocheck 2>&1

    cd \${WORK_DIR}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/protein.tsv
    touch ${prefix}/peptide.tsv
    touch ${prefix}/psm.tsv
    touch ${prefix}/philosopher_report.log
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
