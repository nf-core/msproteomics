include { ptmProphetCommand; shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process PTMPROPHET {
    tag "${meta.id}"
    label 'process_medium'

    input:
    // pepxml_file: ONE validated interact pepXML; ranked lists are split by the caller
    //   (CmdPtmProphet.java:107-124 runs one PTMProphet per pepXML)
    // ptmprophet_cli: PtmProphetPanel.getCmdLineOpts(), "--" style already translated
    //   (FragPipe's plan, bin/plan_to_tool_configs.py)
    // spectra: the LC-MS files the pepXML refers to (searched file + MSFragger mzML), lcms_names: the
    //   names FragPipe writes as base_name; RewritePepxml points base_name at the staged copies
    tuple val(meta), path(pepxml_file, stageAs: 'input/*'), path(spectra), val(lcms_names), val(ptmprophet_cli)

    output:
    tuple val(meta), path("${prefix}")                 , emit: results_dir
    tuple val(meta), path("${prefix}/*.mod.pep.xml")   , emit: mod_pepxml
    tuple val("${task.process}"), val('ptmprophet'), val(fragpipeToolVersions().ptmprophet), emit: versions_ptmprophet, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = ptmProphetCommand(pepxml_file, ptmprophet_cli, args, task)
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    [ -x ${shellQuote(cmds.argv[0])} ] || { echo "ERROR: ${cmds.argv[0]} not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # Not a FragPipe step: FragPipe runs PTMProphet alone in the group dir that holds the pepXML (CmdPtmProphet.java:107-127).
    # PTMProphet reads scans from <base_name>.mzML and fails without it ("cannot read scan in data file"); the
    # input pepXML's base_name names another task's work dir, so RewritePepxml (argv as CmdPeptideProphet.java:426-451)
    # points it at this task's staged spectra. Only base_name changes: the golden mod.pep.xml matches either way.
    cp -L input/${shellQuote(cmds.input)} .
    rewrite_cmd=(java -Duser.home="\$PWD" -cp ${shellQuote(cmds.fragpipe_lib + '/*')} org.nesvilab.fragpipe.util.RewritePepxml "\$PWD"/${shellQuote(cmds.input)} ${[lcms_names].flatten().collect { n -> '"$PWD"/' + shellQuote(n) }.join(' ')})
    printf '%s\\n' "\${rewrite_cmd[*]}" > ${prefix}/ptmprophet.log
    "\${rewrite_cmd[@]}" 2>&1 | tee -a ${prefix}/ptmprophet.log

    # PTMProphet in the pepXML directory, argv as CmdPtmProphet.java:111-122
    ptm_cmd=(${cmds.argv.collect { token -> shellQuote(token) }.join(' ')})
    printf '%s\\n' "\${ptm_cmd[*]}" >> ${prefix}/ptmprophet.log
    "\${ptm_cmd[@]}" 2>&1 | tee -a ${prefix}/ptmprophet.log

    mv ${shellQuote(cmds.output)} ${prefix}/
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = ptmProphetCommand(pepxml_file, ptmprophet_cli, task.ext.args ?: '', task)
    """
    mkdir -p ${prefix}
    touch ${prefix}/${cmds.output}
    printf '%s\\n' ${shellQuote((['java', '-cp', cmds.fragpipe_lib + '/*', 'org.nesvilab.fragpipe.util.RewritePepxml', '{CWD}/' + cmds.input] + [lcms_names].flatten().collect { n -> '{CWD}/' + n }).join(' '))} ${shellQuote(cmds.argv.join(' '))} > ${prefix}/ptmprophet.log
    printf '%s\\n' ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
