include { opairCommands; taskDirArg; taskFileCommand; shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'
include { textFileCommand } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/write_text_file'

/*
 * OPAIR: O-Pair O-glycan localization for one FragPipe group (CmdOPair.java)
 *
 * Input:
 *   - meta: group metadata map (meta.id = group directory name)
 *   - psm_dir: the group directory with psm.tsv (copied: O-Pair writes its results into psm.tsv)
 *   - spectra: the group's LC-MS files as FragPipe lists them (<base>_calibrated.mzML,
 *     <base>_uncalibrated.mzML or the input file, CmdOPair.java:107-114)
 *   - pairs: PairScans <base>.pairs files ([] for single scan type); O-Pair finds them next to
 *     the listed spectra after stripping _calibrated/_uncalibrated
 *   - oxonium_file: the user's oxonium ion definitions file (opair.oxonium_filtering_file, "-f", CmdOPair.java:178-180)
 *     or [] for FragPipe's default in the tools dir (CmdOPair.java:171-177)
 *   - filelist: spectra file names, in FragPipe's order, for filelist_opair.txt
 *   - argv: O-Pair argv from PARSEFRAGPIPEWORKFLOW (one token per flag, @TOKENS@ filled here (@OXONIUM_FILE@ = the staged oxonium file))
 *
 * Output:
 *   - psm_dir: the group directory: psm.tsv with O-Pair results, filelist_opair.txt, O-Pair output
 *   - log: opair.log
 */
process OPAIR {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(psm_dir, stageAs: 'input/*'), path(spectra, stageAs: 'spectra/*'), path(pairs, stageAs: 'spectra/*'), path(oxonium_file, stageAs: 'oxonium/*'), val(filelist), val(argv)

    output:
    tuple val(meta), path("${meta.id}", type: 'dir'), emit: psm_dir
    tuple val(meta), path("opair.log")               , emit: log
    tuple val("${task.process}"), val('opair'), val(fragpipeToolVersions().opair), emit: versions_opair, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // staged as input/<name>: compare the last path element
    if (psm_dir.toString().tokenize("/").last() != meta.id) {
        error("OPAIR: group directory ${psm_dir} does not match meta.id ${meta.id}")
    }
    def cmds = opairCommands(meta.id, filelist, argv, oxonium_file, task)
    """
    export HOME=\$(pwd)
    TASK_DIR=\$(pwd)

    # O-Pair rewrites psm.tsv: run it on a copy of the group directory
    cp -rL ${shellQuote("input/${meta.id}")} ${shellQuote(meta.id)}
    test -f ${shellQuote("${meta.id}/psm.tsv")}

    # filelist_opair.txt in the group directory, absolute paths (CmdOPair.java:100-117)
    for NAME in ${cmds.spectra.collect { n -> shellQuote(n) }.join(' ')}; do
        test -e "spectra/\$NAME"
    done
    ${taskFileCommand("${meta.id}/filelist_opair.txt", cmds.filelist)}

    [ -f ${shellQuote(cmds.dll)} ] || { echo ${shellQuote("ERROR: Could not find O-Pair executable file ${cmds.dll}")} >&2; exit 1; }
    # FragPipe stops when the glycan definition files are missing (CmdOPair.java:150-190)
    for TOOL_FILE in ${cmds.tool_files.collect { f -> shellQuote(f) }.join(' ')}; do
        [ -f "\$TOOL_FILE" ] || { echo "ERROR: O-Pair file \$TOOL_FILE not found" >&2; exit 1; }
    done

    echo ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # FragPipe runs O-Pair from its work dir with one argv token per flag ("-b 20.0")
    CMD=( dotnet ${shellQuote(cmds.dll)} ${cmds.argv.collect { t -> taskDirArg(t) }.join(' ')} ${args} )
    printf '%q ' "\${CMD[@]}" > opair.log
    printf '\\n' >> opair.log
    "\${CMD[@]}" 2>&1 | tee -a opair.log
    """

    stub:
    def cmds = opairCommands(meta.id, filelist, argv, oxonium_file, task)
    // relative paths in the stub (deterministic snapshots); the log holds the filled argv, {TASK_DIR} unexpanded
    """
    cp -rL ${shellQuote("input/${meta.id}")} ${shellQuote(meta.id)}
    ${taskFileCommand("${meta.id}/filelist_opair.txt", cmds.filelist.replace('{TASK_DIR}/', ''))}
    ${textFileCommand('opair.log', cmds.argv.collect { t -> "${t}\n" }.join(''))}
    echo ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
