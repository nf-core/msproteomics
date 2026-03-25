process PTMPROPHET {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // PTMProphet performs PTM site localization
    // pepxml: search results with PTMs
    // config_cli: PTMProphet parameters as CLI string (coupled with sample inputs)
    tuple val(meta), path(pepxml_file), val(config_cli)

    output:
    tuple val(meta), path("${prefix}")                          , emit: results_dir
    tuple val(meta), path("${prefix}/*.mod.pep.xml")            , emit: mod_pepxml
    tuple val("${task.process}"), val('ptmprophet'), eval("cat .ptmprophet_version"), emit: versions_ptmprophet, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def decoy_tag = task.ext.decoy_tag ?: 'rev_'
    def basename = pepxml_file.baseName
    // Replace .pep with .mod.pep in output filename
    def output_file = basename.replaceAll(/\.pep$/, '') + ".mod.pep.xml"
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    // PTMProphet uses UPPERCASE=value format for parameters
    // Common params: MINPROB=0.5, MAXTHREADS=1, STATIC, KEEPOLD, LABILEMODS, NOSTACK
    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    PTMPROPHET_SRC=\$(find "\$TOOLS_DIR" -path '*/PTMProphet/PTMProphetParser*' -type f 2>/dev/null | head -1 || true)
    PTMPROPHET="\$(pwd)/PTMProphetParser"
    cp "\$PTMPROPHET_SRC" "\$PTMPROPHET"
    chmod +x "\$PTMPROPHET"
    # Version capture
    _ver=\$(ls "\$TOOLS_DIR"/PTMProphet/PTMProphetParser-* 2>/dev/null | grep -oP 'Parser-\\K[\\d.]+' | head -1 || true)

    echo "\$_ver" > .ptmprophet_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    # Copy input pepxml to output directory (PTMProphet writes output alongside input)
    cp ${pepxml_file} "${prefix}/"

    # FragPipe always sets MAXTHREADS=1 (CmdPtmProphet.java:118-119) because PTMProphet
    # is run as many parallel instances (one per pepXML). Using MAXTHREADS>1 within a single
    # PTMProphet instance provides minimal benefit and can cause thread contention.
    # Nextflow handles parallelism at the process level instead.
    # Strip any existing MAXTHREADS= from config_cli to avoid duplication
    CLEAN_CLI=\$(echo "${config_cli}" | sed 's/MAXTHREADS=[0-9]*//g' | tr -s ' ')
    PTMPROPHET_CMD="\$PTMPROPHET MAXTHREADS=1 \$CLEAN_CLI ${args} ${prefix}/${pepxml_file.name} ${prefix}/${output_file}"
    printf '%s\\n' "\$PTMPROPHET_CMD" > "${prefix}/ptmprophet.log"
    \$PTMPROPHET MAXTHREADS=1 \$CLEAN_CLI ${args} "${prefix}/${pepxml_file.name}" "${prefix}/${output_file}" 2>&1 | tee -a "${prefix}/ptmprophet.log"

    # Clean up copied input file
    rm -f "${prefix}/${pepxml_file.name}"

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def basename = pepxml_file.baseName
    def output_file = basename.replaceAll(/\.pep$/, '') + ".mod.pep.xml"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${output_file}
    touch ${prefix}/ptmprophet.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .ptmprophet_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
