include { iprophetCommand; shellQuote; fragpipeToolVersions } from '../../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process IPROPHET {
    tag "${meta.id}"
    label 'process_high'

    input:
    // pepxml_files: validated pepXMLs of every group (sharedPepxmlFiles, CmdIprophet.java:53)
    // pepxml_paths: per pepXML its FragPipe location relative to the work dir ("<group>/<file>",
    //   or "<file>" for the blank group), in FragPipe order (InputLcmsFile.compareTo)
    // iprophet_cli: "--decoy <tag> --nonsp --output combined" (CmdIprophet.java:44-50)
    tuple val(meta), path(pepxml_files, stageAs: 'pepxml_inputs/?/*'), val(pepxml_paths), val(iprophet_cli)

    output:
    tuple val(meta), path("${prefix}")                   , emit: results_dir
    tuple val(meta), path("${prefix}/combined.pep.xml")  , emit: pepxml
    tuple val("${task.process}"), val('philosopher'), val(fragpipeToolVersions().philosopher), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = iprophetCommand(pepxml_files, pepxml_paths, iprophet_cli, args, task)
    def stage = cmds.links.collect { staged, rel ->
        def dir = rel.contains('/') ? rel[0..<rel.lastIndexOf('/')] : ''
        (dir ? "mkdir -p ${shellQuote(dir)}\n    " : '') + "ln -s \"\$PWD\"/${shellQuote(staged)} ${shellQuote(rel)}"
    }.join('\n    ')
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    [ -x ${shellQuote(cmds.argv[0])} ] || { echo "ERROR: ${cmds.argv[0]} not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # pepXMLs at their FragPipe locations
    ${stage}

    # Work-dir workspace (CmdPhilosopherWorkspaceCleanInit.java:66-74)
    phi_temp="\$PWD/phi_temp/\$(cat /proc/sys/kernel/random/uuid)"
    mkdir -p "\$phi_temp"
    ${shellQuote(cmds.argv[0])} workspace --init --nocheck --temp "\$phi_temp"

    # iProphet in the work dir (CmdIprophet.java:41-58)
    iprophet_cmd=(${cmds.argv.collect { t -> shellQuote(t) }.join(' ')} ${cmds.paths.collect { rel -> "\"\$PWD\"/${shellQuote(rel)}" }.join(' ')})
    printf '%s\\n' "\${iprophet_cmd[*]}" > ${prefix}/iprophet.log
    "\${iprophet_cmd[@]}" 2>&1 | tee -a ${prefix}/iprophet.log

    mv combined.pep.xml ${prefix}/
    rm -rf phi_temp
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = iprophetCommand(pepxml_files, pepxml_paths, iprophet_cli, task.ext.args ?: '', task)
    """
    mkdir -p ${prefix}
    touch ${prefix}/combined.pep.xml
    printf '%s\\n' ${shellQuote((cmds.argv + cmds.paths.collect { rel -> "{CWD}/${rel}" }).join(' '))} > ${prefix}/iprophet.log
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
