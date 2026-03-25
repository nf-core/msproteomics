process PEPTIDEPROPHET {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    tuple val(meta), path(pepxml), path(mzml), path(config_file)
    path(fasta)

    output:
    tuple val(meta), path("${prefix}")                    , emit: results_dir
    tuple val(meta), path("${prefix}/interact-*.pep.xml") , emit: pepxml
    tuple val("${task.process}"), val('philosopher'), eval("cat .philosopher_version"), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def decoy_tag = task.ext.decoy_tag ?: 'rev_'
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    // Config file parsing is done in bash for compatibility with nf-test
    // Note: Enzyme-conditional flags (--nontt --nonmc, --enzyme nonspecific) for
    // nonspecific/dual-enzyme/nocleavage searches (CmdPeptideProphet.java:391-404)
    // are already included in config_file by parse_fragpipe_workflow.py via
    // _get_peptideprophet_enzyme_flags().

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

    export GOMAXPROCS=${task.cpus}
    export GOMEMLIMIT="${ram}GiB"
    export GOGC=200

    # Parse config file for peptideprophet flags (format: peptideprophet=flags)
    CMD_FLAGS=\$(grep -E '^peptideprophet=' ${config_file} | sed 's/^peptideprophet=//')
    if [ -z "\$CMD_FLAGS" ]; then
        echo "ERROR: No 'peptideprophet=' key found in config file ${config_file}" >&2
        exit 1
    fi

    # Initialize philosopher workspace
    \$PHILOSOPHER workspace --init --nocheck 2>&1

    CMD="\$PHILOSOPHER peptideprophet --decoy ${decoy_tag} --database ${fasta} \$CMD_FLAGS ${args} ${pepxml}"
    printf '%s\\n' "\$CMD" > ${prefix}/peptideprophet.log
    \$CMD 2>&1 | tee -a ${prefix}/peptideprophet.log

    shopt -s nullglob
    for f in interact-*.pep.xml; do mv "\$f" "${prefix}/"; done
    shopt -u nullglob

    # Clean up philosopher workspace
    \$PHILOSOPHER workspace --clean --nocheck 2>&1

    # Normalize base_name attribute in pepXML: PeptideProphet embeds absolute paths from
    # the Nextflow work directory (e.g., "/work/ab/123.../sample1") which become invalid
    # after files are moved. Downstream tools (ProteinProphet, Philosopher filter/report)
    # use this attribute to associate peptides with source files, so we simplify it to
    # just the base name (e.g., "sample1").
    for f in ${prefix}/interact-*.pep.xml; do
        if [[ -f "\$f" ]]; then
            sed -i "s|<msms_run_summary base_name=\\"[^\\"]*\\"|<msms_run_summary base_name=\\"${mzml.baseName}\\"|g" "\$f"
        fi
    done

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def basename = pepxml.baseName.replaceAll('\\.pepXML$', '')
    """
    mkdir -p ${prefix}
    touch ${prefix}/interact-${basename}.pep.xml
    touch ${prefix}/peptideprophet.log
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
