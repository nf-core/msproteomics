include { fragpipeJarLines; fragpipeTools; shellArgs; rawLinkCopies } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { ptmShepherdCommands; taskFileCommand; shellQuote; fragpipeToolVersions; javaXmxGb } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

/*
 * PTMSHEPHERD: PTM-Shepherd mass-shift profiling and glycan assignment (CmdPtmshepherd.java)
 *
 * Input:
 *   - meta: aggregate metadata map
 *   - psm_dirs: one directory per FragPipe group (named meta.id of the group), each with psm.tsv
 *   - datasets: [dir_name, dataset_name] per group, in FragPipe's group order; the blank group's
 *     dataset is "dataset01" (PtmshepherdParams.java:86)
 *   - mzml_files: LC-MS files of all groups (one directory, as FragPipe requires)
 *   - config_text: sorted property block from PARSEFRAGPIPEWORKFLOW (@THREADS@/@TOOLS_DIR@/@ANNOTATION_FILE@ tokens)
 *   - annotation_file: the custom annotation file of ptmshepherd.annotation_file (annotation-custom/annotation-glyco
 *     radio, PtmshepherdPanel.java:229-232), staged as annotation/<name>; [] when annotation_file is a keyword
 *   - fasta: the database FragPipe searched
 *
 * Output:
 *   - psm_dirs: copies of the group directories, psm.tsv as PTM-Shepherd rewrote it
 *   - results_dir: ptm-shepherd-output/ (all of it)
 *   - glyco_masses: ptm-shepherd-output/glyco_masses_list.txt (glycan assignment only), for AppendFile
 *   - config: shepherd.config as run
 */
process PTMSHEPHERD {
    tag "${meta.id}"
    label 'process_high'

    input:
    // config_text: shepherd.config property block from the workflow parser (@THREADS@, @TOOLS_DIR@, @ANNOTATION_FILE@ tokens)
    tuple val(meta), path(psm_dirs, stageAs: 'input/*'), val(datasets), path(mzml_files, stageAs: 'spectra/*'), val(config_text), path(annotation_file, stageAs: 'annotation/*')
    path fasta

    output:
    tuple val(meta), path("groups/*", type: 'dir')                         , emit: psm_dirs
    tuple val(meta), path("ptm-shepherd-output")                           , emit: results_dir
    tuple val(meta), path("ptm-shepherd-output/glyco_masses_list.txt")     , emit: glyco_masses, optional: true
    tuple val(meta), path("shepherd.config")                               , emit: config
    tuple val(meta), path("ptmshepherd.log")                               , emit: log
    tuple val("${task.process}"), val('ptmshepherd'), val(fragpipeToolVersions().ptmshepherd), emit: versions_ptmshepherd, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = javaXmxGb(task, 2)
    def cmds = ptmShepherdCommands(datasets, config_text, fasta.name, annotation_file, task)
    // IonQuant joins the class path only for kl glyco features; the Thermo libraries are optional
    def jar_lines = fragpipeJarLines(fragpipeTools(task), ['batmass_io'] + (cmds.class_path.contains('{IONQUANT_JAR}') ? ['ionquant'] : []), ['msfragger'])
    """
    ${rawLinkCopies(mzml_files)}
    export HOME=\$(pwd)
    TASK_DIR=\$(pwd)

    # PTM-Shepherd rewrites psm.tsv in place; work on copies so upstream task dirs stay untouched
    mkdir -p groups
    ${cmds.groups.collect { g -> "cp -rL ${shellQuote("input/${g}")} ${shellQuote("groups/${g}")}\n    test -f ${shellQuote("groups/${g}/psm.tsv")}" }.join('\n    ')}

    # shepherd.config as PtmshepherdParams.createConfig() (CmdPtmshepherd runs it from the work dir)
    ${taskFileCommand('shepherd.config', cmds.config)}

    # Glycan definition files FragPipe takes from its tools dir (PTMSGlycanAssignPanel.java:131-134)
    for GLYCO_FILE in ${cmds.glyco_files.collect { f -> shellQuote(f) }.join(' ')}; do
        [ -f "\$GLYCO_FILE" ] || { echo "ERROR: glycan definition file \$GLYCO_FILE not found" >&2; exit 1; }
    done
    ${jar_lines}
    for JAR in ${cmds.class_path.collect { j -> shellQuote(j) }.join(' ')}; do
        [ -f "\$JAR" ] || { echo "ERROR: PTM-Shepherd class path jar \$JAR missing (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    done
    JAVA_ARGS=( "-Duser.home=\$PWD" "-Xmx${ram}G" )
    # Thermo libraries ship next to the MSFragger jar (FragpipeRun.java:1348-1356)
    if [ -d ${shellQuote(cmds.thermo_dir)} ]; then
        JAVA_ARGS+=( ${shellQuote("-Dlibs.thermo.dir=${cmds.thermo_dir}")} )
    fi

    echo ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # CmdPtmshepherd.java:163 passes the config path wrapped in literal quotes (a ProcessBuilder, no shell: the
    # element holds the quote characters; PTMShepherd.parseParamFile: replaceAll("['\"]", "") strips them)
    CMD=( java "\${JAVA_ARGS[@]}" -cp ${shellQuote(cmds.class_path.join(':'))} edu.umich.andykong.ptmshepherd.PTMShepherd ${shellArgs([cmds.config_arg])} ${args} )
    printf '%q ' "\${CMD[@]}" > ptmshepherd.log
    printf '\\n' >> ptmshepherd.log
    "\${CMD[@]}" 2>&1 | tee -a ptmshepherd.log
    test -d ptm-shepherd-output
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def cmds = ptmShepherdCommands(datasets, config_text, fasta.name, annotation_file, task)
    // the stub writes the config with paths relative to the task dir (deterministic snapshots)
    """
    mkdir -p groups ptm-shepherd-output
    ${cmds.groups.collect { g -> "cp -rL ${shellQuote("input/${g}")} ${shellQuote("groups/${g}")}" }.join('\n    ')}
    ${taskFileCommand('shepherd.config', cmds.config.replace('{TASK_DIR}/', ''))}
    touch ptm-shepherd-output/global.profile.tsv ptm-shepherd-output/global.modsummary.tsv ptmshepherd.log
    ${cmds.glyco_mode ? 'touch ptm-shepherd-output/glyco_masses_list.txt' : ''}
    echo ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
