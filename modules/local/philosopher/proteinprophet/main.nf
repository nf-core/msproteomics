process PROTEINPROPHET {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    tuple val(meta), path(pepxml_files), val(config_cli)  // pepXML files + config as CLI string (e.g., "--maxppmdiff 2000000")

    output:
    tuple val(meta), path("${prefix}")                    , emit: results_dir
    tuple val(meta), path("${prefix}/combined.prot.xml")  , emit: protxml
    tuple val(meta), path(pepxml_files)                   , emit: pepxml  // Pass through
    tuple val("${task.process}"), val('philosopher'), eval("cat .philosopher_version"), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def pepxml_list = pepxml_files instanceof List ? pepxml_files : [pepxml_files]
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
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

    export GOMAXPROCS=${task.cpus}
    export GOMEMLIMIT="${ram}GiB"
    export GOGC=200

    > filelist_proteinprophet.txt
    for f in ${pepxml_list.join(' ')}; do
        echo "\$(pwd)/\$f" >> filelist_proteinprophet.txt
    done

    export XML_ONLY=1
    \$PHILOSOPHER workspace --init --nocheck 2>&1

    CMD="\$PHILOSOPHER proteinprophet ${config_cli} --output combined filelist_proteinprophet.txt ${args}"
    printf '%s\\n' "\$CMD" > ${prefix}/proteinprophet.log

    # Run ProteinProphet
    \$CMD 2>&1 | tee -a ${prefix}/proteinprophet.log

    mv combined.prot.xml ${prefix}/

    \$PHILOSOPHER workspace --clean --nocheck 2>&1

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/combined.prot.xml
    touch ${prefix}/proteinprophet.log
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
