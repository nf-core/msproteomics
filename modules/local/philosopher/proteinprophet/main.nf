include { proteinProphetCommand; shellQuote; fragpipeToolVersions } from '../../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process PROTEINPROPHET {
    tag "${meta.id}"
    label 'process_high'

    input:
    // pepxml_files: plain interact pepXMLs of every group (FragpipeRun.java:1677-1680 feeds
    //   sharedPepxmlFiles, never PTMProphet's .mod.pep.xml)
    // filelist_paths: per pepXML its FragPipe location relative to the work dir ("<group>/<file>",
    //   or "<file>" for the blank group), in FragPipe order (InputLcmsFile.compareTo); the
    //   filelist holds their absolute paths (CmdProteinProphet.java:150-166), and protXML
    //   source_files record them
    // proteinprophet_cli: workflow protein-prophet.cmd-opts
    tuple val(meta), path(pepxml_files, stageAs: 'pepxml_inputs/?/*'), val(filelist_paths), val(proteinprophet_cli)
    // database: the run's FASTA. ProteinProphet opens the file named in the pepXML header (search_database
    //   local_path, a symlink to this FASTA in the work dir of the search or of the split-database merge) for
    //   protein length, coverage and mu; staging it makes its real location visible in the container. Without it ProteinProphet warns
    //   "Cannot open DB", writes prot_length 0 and no percent_coverage
    path(database, stageAs: 'database/*')

    output:
    tuple val(meta), path("${prefix}")              , emit: results_dir
    tuple val(meta), path("${prefix}/*.prot.xml")   , emit: protxml
    tuple val("${task.process}"), val('philosopher'), val(fragpipeToolVersions().philosopher), emit: versions_philosopher, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = proteinProphetCommand(pepxml_files, filelist_paths, proteinprophet_cli, args, task)
    def stage = cmds.links.collect { staged, rel ->
        def dir = rel.contains('/') ? rel[0..<rel.lastIndexOf('/')] : ''
        (dir ? "mkdir -p ${shellQuote(dir)}\n    " : '') +
            "ln -s \"\$PWD\"/${shellQuote(staged)} ${shellQuote(rel)}\n    printf '%s\\n' \"\$PWD\"/${shellQuote(rel)} >> filelist_proteinprophet.txt"
    }.join('\n    ')
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    [ -x ${shellQuote(cmds.philosopher)} ] || { echo "ERROR: ${cmds.philosopher} not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # pepXMLs at their FragPipe locations, listed in FragPipe order (CmdProteinProphet.java:150-166)
    : > filelist_proteinprophet.txt
    ${stage}

    # Work-dir workspace (CmdPhilosopherWorkspaceCleanInit.java:66-74)
    phi_temp="\$PWD/phi_temp/\$(cat /proc/sys/kernel/random/uuid)"
    mkdir -p "\$phi_temp"
    ${shellQuote(cmds.philosopher)} workspace --init --nocheck --temp "\$phi_temp"

    # ProteinProphet in the work dir with XML_ONLY=1 and Philosopher's dir appended to PATH as
    # FragPipe does (CmdProteinProphet.java:172-216)
    pp_cmd=(${cmds.argv.collect { token -> shellQuote(token) }.join(' ')} "\$PWD"/filelist_proteinprophet.txt)
    printf '%s\\n' "\${pp_cmd[*]}" > ${prefix}/proteinprophet.log
    XML_ONLY=1 PATH="\$PATH;;${cmds.philosopher_dir}" "\${pp_cmd[@]}" 2>&1 | tee -a ${prefix}/proteinprophet.log

    mv ${shellQuote(cmds.output)} ${prefix}/
    rm -rf phi_temp
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = proteinProphetCommand(pepxml_files, filelist_paths, proteinprophet_cli, task.ext.args ?: '', task)
    """
    mkdir -p ${prefix}
    touch ${prefix}/${cmds.output}
    printf '%s\\n' ${shellQuote(cmds.argv.join(' ') + ' {CWD}/filelist_proteinprophet.txt')} > ${prefix}/proteinprophet.log
    printf '{CWD}/%s\\n' ${cmds.filelist.collect { rel -> shellQuote(rel) }.join(' ')} > ${prefix}/filelist_proteinprophet.txt
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
