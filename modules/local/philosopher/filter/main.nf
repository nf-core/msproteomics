process PHILOSOPHER_FILTER {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    tuple val(meta), path(pepxml_files), path(protxml), val(filter_flags), val(report_flags)
    path fasta
    // fasta: Philosopher database with decoys/contaminants from PHILOSOPHER_DATABASE
    // filter_flags: e.g., "--sequential --razor --prot 0.01 --picked"
    // report_flags: e.g., "--msstats"
    // NOTE: Philosopher filter+report do NOT read mzML files. Per-sample psm.tsv Intensity
    // starts at 0.0 after report. IonQuant updates psm.tsv in-place with intensities from mzML
    // (IonQuant bypasses .meta/ binaries, reads/writes TSV directly).
    // NOTE: FragPipe passes --dbbin for multi-group runs (CmdPhilosopherFilter.java:103-106),
    // pointing to the first group's directory. This module initializes its own workspace per sample,
    // so --dbbin is not needed for single-experiment runs. Multi-group support would require
    // passing a shared database workspace path.

    output:
    tuple val(meta), path("${prefix}")              , emit: results_dir
    tuple val(meta), path("${prefix}/psm.tsv")      , emit: psms
    tuple val(meta), path("${prefix}/peptide.tsv")  , emit: peptides
    tuple val(meta), path("${prefix}/protein.tsv")  , emit: proteins
    tuple val(meta), path("${prefix}/ion.tsv")      , emit: ions, optional: true
    tuple val("${task.process}"), val('philosopher'), eval("cat .philosopher_version"), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def decoy_prefix = task.ext.decoy_tag ?: 'rev_'
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def skip_report = task.ext.skip_report ?: false
    tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    // Config passed as CLI strings directly
    def filter_cli = filter_flags
    def report_cli = report_flags

    // When skip_report is true, the .meta/ workspace binaries are preserved in the output
    // so that downstream tools (FreeQuant, LabelQuant) can update them before a separate
    // PHILOSOPHER_REPORT step generates TSVs with correct intensities.
    // When skip_report is false (default), behavior is unchanged: filter + report + clean.

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
    \$PHILOSOPHER database --annotate ${fasta} --prefix ${decoy_prefix} 2>&1 | tee ${prefix}/philosopher_filter.log

    # Run filter and report - fail immediately on error (matching FragPipe behavior)
    # --pepxml .: FragPipe passes the group working directory (CmdPhilosopherFilter.java:101-102).
    # In Nextflow, staged files are in the current directory, so "." is equivalent.
    FILTER_CMD="\$PHILOSOPHER filter --tag ${decoy_prefix} --pepxml . --protxml ${protxml} ${filter_cli} ${args}"
    printf '\\n%s\\n' "\$FILTER_CMD" >> ${prefix}/philosopher_filter.log
    \$FILTER_CMD 2>&1 | tee -a ${prefix}/philosopher_filter.log

    if [ "${skip_report}" = "false" ]; then
        REPORT_CMD="\$PHILOSOPHER report ${report_cli}"
        printf '\\n%s\\n' "\$REPORT_CMD" >> ${prefix}/philosopher_filter.log
        \$REPORT_CMD 2>&1 | tee -a ${prefix}/philosopher_filter.log

        for f in psm.tsv peptide.tsv protein.tsv ion.tsv; do
            [ -f "\$f" ] && mv "\$f" "${prefix}/"
        done

        \$PHILOSOPHER workspace --clean --nocheck 2>&1
    else
        # skip_report=true: preserve .meta/ workspace for downstream FreeQuant/LabelQuant.
        # philosopher filter only updates .meta/ binaries; it does NOT produce TSV files.
        # Create placeholder TSVs so Nextflow output declarations are satisfied.
        # PHILOSOPHER_REPORT will regenerate these with correct intensities after FreeQuant.
        touch ${prefix}/psm.tsv
        touch ${prefix}/peptide.tsv
        touch ${prefix}/protein.tsv

        # Copy .meta/ workspace to output directory so downstream tools can use it
        cp -a .meta ${prefix}/
    fi

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/protein.tsv
    touch ${prefix}/peptide.tsv
    touch ${prefix}/psm.tsv
    touch ${prefix}/philosopher_filter.log
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
