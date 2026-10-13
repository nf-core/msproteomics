include { fragpipeJarLines; fragpipeTools; shellArgs; rawLinkCopies } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { pairScansCommands; shellQuote; javaXmxGb } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'
include { textFileCommand } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/write_text_file'

/*
 * PAIRSCANS: pair parent/child scans of one LC-MS file for O-Pair (CmdPairScans.java, util/PairScans.java)
 *
 * Input:
 *   - meta: spectra metadata map
 *   - spectra: the LC-MS file FragPipe lists (mzML, mzXML or Thermo .raw)
 *   - pairscans_args: [activation1, activation2, reverse_scan_order, single_scan_type]
 *
 * Output:
 *   - pairs: <basename>.pairs, written next to the spectra file as PairScans does (none when
 *     PairScans skips the file)
 *   - log: <basename>.pairscans.log
 */
process PAIRSCANS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(spectra, stageAs: 'spectra/*'), val(pairscans_args)

    output:
    // optional: PairScans writes nothing for an unknown activation type or format (PairScans.java:80-88);
    // O-Pair then treats the file as single scan type, as in FragPipe
    tuple val(meta), path("spectra/${spectra.baseName}.pairs"), emit: pairs, optional: true
    tuple val(meta), path("${spectra.baseName}.pairscans.log"), emit: log
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ram = javaXmxGb(task, 2)
    def cmds = pairScansCommands(spectra.toString(), pairscans_args, task)
    def is_raw = spectra.name.toLowerCase().endsWith('.raw')
    """
    ${rawLinkCopies(spectra)}
    export HOME=\$(pwd)
    ${fragpipeJarLines(fragpipeTools(task), ['batmass_io'], ['msfragger'])}
    for JAR in ${cmds.jars.collect { j -> shellQuote(j) }.join(' ')}; do
        [ -f "\$JAR" ] || { echo "ERROR: Could not find the FragPipe jar and batmass-io for scan pairing (\$JAR)" >&2; exit 1; }
    done
    ${is_raw ? "[ -d ${shellQuote(cmds.thermo_dir)} ] || { echo 'ERROR: When processing .RAW files PairScans requires native Thermo libraries' >&2; exit 1; }" : ''}

    # java -Xmx<ram>G [bruker] [thermo] -cp <fragpipe jar>:<batmass-io> PairScans <file> <threads> ... (CmdPairScans.java:108-127);
    # the Thermo / Bruker libraries ship next to the MSFragger jar (FragpipeRun.java:1348-1356)
    JAVA_ARGS=( "-Duser.home=\$PWD" "-Xmx${ram}G" )
    if [ -d ${shellQuote(cmds.bruker_dir)} ]; then
        JAVA_ARGS+=( ${shellQuote("-Dlibs.bruker.dir=${cmds.bruker_dir}")} )
    fi
    if [ -d ${shellQuote(cmds.thermo_dir)} ]; then
        JAVA_ARGS+=( ${shellQuote("-Dlibs.thermo.dir=${cmds.thermo_dir}")} )
    fi
    echo ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # PairScans writes <basename>.pairs next to its input, i.e. into this task's spectra/ dir
    CMD=( java "\${JAVA_ARGS[@]}" ${shellArgs(cmds.argv)} ${args} )
    printf '%q ' "\${CMD[@]}" > ${spectra.baseName}.pairscans.log
    printf '\\n' >> ${spectra.baseName}.pairscans.log
    "\${CMD[@]}" 2>&1 | tee -a ${spectra.baseName}.pairscans.log
    """

    stub:
    def cmds = pairScansCommands(spectra.toString(), pairscans_args, task)
    // the log holds the argv after the JVM options, {TASK_DIR} unexpanded (deterministic snapshots)
    """
    ${textFileCommand("${spectra.baseName}.pairscans.log", cmds.argv.collect { t -> "${t}\n" }.join(''))}
    touch spectra/${spectra.baseName}.pairs
    echo ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
