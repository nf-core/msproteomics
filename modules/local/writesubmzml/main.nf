/*
 * WRITESUBMZML: FragPipe "Write sub mzML" (CmdWriteSubMzml / util.WriteSubMzml)
 *
 * FragPipe 24.0 runs it per LC-MS file after the final Philosopher report when
 * tab-run.write_sub_mzml=true (FragpipeRun.java:2388-2404, graph :2489). In <wd>/<group> it runs
 * (CmdWriteSubMzml.java:97-130):
 *   java -Xmx<ram>G -Dlibs.bruker.dir=<dir> -Dlibs.thermo.dir=<dir> -cp <batmass-io jar>:<install>/lib/*
 *     org.nesvilab.fragpipe.util.WriteSubMzml <LC-MS file> <wd>/<group>/psm.tsv <wd>/<stem>_sub.mzML <threshold> 1
 * WriteSubMzml reads <LC-MS stem>.mzBIN_all next to the LC-MS file (MSFragger's write_mzbin_all output),
 * drops the MS2 scans whose PSM probability in psm.tsv is above the threshold (and the MS3 scans of
 * their parents), writes the other scans with MSFragger's calibrated peaks into <stem>_sub.mzML, and
 * deletes the .mzBIN_all on exit (WriteSubMzml.java:72-100, :102-161, :176-215).
 * The argv comes from writeSubMzmlArgv (fragpipe_utils_a.nf); the task dir stands for FragPipe's work dir.
 *
 * Input:
 *   - meta: metadata of the LC-MS file (meta.group: its FragPipe group)
 *   - lcms: the LC-MS file as given (FragPipe passes the input file), staged in lcms/
 *   - mzbin_all: MSFragger's <stem>.mzBIN_all of the file, staged next to it
 *   - psm: psm.tsv of the file's group (final Philosopher report)
 *   - prob_threshold: tab-run.sub_mzml_prob_threshold (tool_configs writesubmzml.prob_threshold)
 *
 * Output:
 *   - sub_mzml: <stem>_sub.mzML
 *   - log: write_sub_mzml.log, the command line followed by the tool output
 *   - versions: FragPipe version (topic: versions)
 */
include { fragpipeTools; writeSubMzmlArgv; shellArgs; rawSymlinksToFiles; fragpipeJarLines; hashSetClassPathLines } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { javaXmxGb; shellQuote } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process WRITESUBMZML {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(lcms, stageAs: 'lcms/*'), path(mzbin_all, stageAs: 'lcms/*'), path(psm, stageAs: 'group_psm/psm.tsv'), val(prob_threshold)

    output:
    tuple val(meta), path("*_sub.mzML")        , emit: sub_mzml
    tuple val(meta), path("write_sub_mzml.log"), emit: log
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // -Xmx<ram>G from the task memory; FragPipe always passes its RAM setting (CmdWriteSubMzml.java:103)
    def ram = javaXmxGb(task, 0)
    def tools = fragpipeTools(task)
    def cmd = writeSubMzmlArgv(lcms.name.tokenize('/').last(), 'lcms', meta, prob_threshold, ram, tools)
    """
    export HOME=\$(pwd)
    echo ${shellQuote(tools.fragpipe_version)} > .fragpipe_version

    # FragPipe's work-dir layout: the group's psm.tsv in <wd>/<group>, the LC-MS file with its .mzBIN_all
    mkdir -p ${shellQuote(cmd.group_dir)}
    ln -s "\$PWD/group_psm/psm.tsv" ${shellQuote(cmd.psm)}
    ${rawSymlinksToFiles("lcms/${lcms.name.tokenize('/').last()}")}
    ${fragpipeJarLines(tools, ['msfragger', 'batmass_io'])}
    ${hashSetClassPathLines(cmd.classpath_words)}

    CMD=(${shellArgs(cmd.argv)} ${args})
    printf '%s\\n' "\${CMD[*]}" > write_sub_mzml.log
    (cd ${shellQuote(cmd.group_dir)} && "\${CMD[@]}") 2>&1 | tee -a write_sub_mzml.log
    """

    stub:
    def tools = fragpipeTools(task)
    def cmd = writeSubMzmlArgv(lcms.name.tokenize('/').last(), 'lcms', meta, prob_threshold, javaXmxGb(task, 0), tools)
    // the log keeps the jar tokens literal (fragpipeJarLines does not run in the stub)
    """
    printf '%s\\n' ${cmd.argv.collect { w -> "'${w}'" }.join(' ')} ${shellQuote("cwd ${cmd.group_dir}")} > write_sub_mzml.log
    touch ${shellQuote(cmd.sub_mzml)}
    echo ${shellQuote(tools.fragpipe_version)} > .fragpipe_version
    """
}
