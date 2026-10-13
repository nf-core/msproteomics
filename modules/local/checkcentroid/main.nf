/*
 * CHECKCENTROID: FragPipe's centroid check (CmdCheckCentroid / util.CheckCentroid)
 *
 * FragPipe runs it once, on the first LC-MS file, before DIA-Umpire, diaTracer and MSFragger
 * (FragpipeRun.java:1333-1340, CmdCheckCentroid.java:52-65). It exits 1 when an mzML/mzXML file
 * has profile-mode scans; other formats pass (CheckCentroid.java:41-58).
 * argv: java -Xmx<ram>G -cp <fragpipe jar>:<batmass-io jar> org.nesvilab.fragpipe.util.CheckCentroid <file> <threads>
 *
 * Input:
 *   - meta: Metadata of the checked LC-MS file
 *   - lcms: The first LC-MS file of the run
 *
 * Output:
 *   - log: check_centroid.log, the command line followed by the tool output
 *   - versions: FragPipe version (topic: versions)
 */
include { fragpipeTools; fragpipeJarLines; shellArgs } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process CHECKCENTROID {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(lcms)

    output:
    tuple val(meta), path("check_centroid.log"), emit: log
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ram = task.memory ? task.memory.toGiga() : 8
    def tools = fragpipeTools(task)
    // Class path: lib/fragpipe-<version>.jar and the batmass-io jar its build names (ToolingUtils.BATMASS_IO_JAR:
    // batmass-io-1.36.5 in the FragPipe 24.0 release, 1.36.6 in FragPipePlus; fragpipeJarLines)
    """
    export HOME=\$(pwd)
    echo "${tools.fragpipe_version}" > .fragpipe_version
    ${fragpipeJarLines(tools, ['batmass_io'])}

    CMD=(java -Duser.home=\$PWD -Xmx${ram}G -cp ${shellArgs(["${tools.fragpipe_jar}:${tools.batmass_io_jar}"])} org.nesvilab.fragpipe.util.CheckCentroid "\$(pwd)/${lcms}" ${task.cpus} ${args})
    printf '%s\\n' "\${CMD[*]}" > check_centroid.log
    "\${CMD[@]}" 2>&1 | tee -a check_centroid.log
    """

    stub:
    def tools = fragpipeTools(task)
    """
    echo "org.nesvilab.fragpipe.util.CheckCentroid \$(pwd)/${lcms} ${task.cpus}" > check_centroid.log
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}
