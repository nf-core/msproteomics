/*
 * MBG: match-between-glycans glycoform inference for one FragPipe group (CmdMBGMatch.java)
 *
 * Input:
 *   - meta: group metadata map (meta.id = group directory name)
 *   - psm_dir: the group directory with the PTM-Shepherd glyco psm.tsv (copied: MBG replaces psm.tsv)
 *   - spectra: every LC-MS file of the run (MBG reads them through the manifest)
 *   - manifest_rows: [file name, experiment, bioreplicate, data type] per LC-MS file, in
 *     fragpipe-files.fp-manifest order (TabWorkflow.manifestSave)
 *   - argv: MBG argv from PARSEFRAGPIPEWORKFLOW (@TOKENS@ filled here, @NOPASEF@ by the workflow)
 *
 * Output:
 *   - psm_dir: the group directory with MBG's updated psm.tsv
 *   - manifest: fragpipe-files.fp-manifest given to MBG
 *   - log: mbg.log
 */
include { fragpipeTools; fragpipeManifestText; writeTextFile; fillArgvTokens; shellArgs; fragpipeJarLines } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process MBG {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(psm_dir, stageAs: 'input/*'), path(spectra, stageAs: 'spectra/*'), val(manifest_rows), val(argv)

    output:
    tuple val(meta), path("${meta.id}", type: 'dir')     , emit: psm_dir
    tuple val(meta), path("fragpipe-files.fp-manifest")  , emit: manifest
    tuple val(meta), path("mbg.log")                     , emit: log
    tuple val("${task.process}"), val('mbg'), eval("cat .mbg_version"), emit: versions_mbg, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    // staged as input/<name>: compare the last path element
    if (psm_dir.toString().tokenize("/").last() != meta.id) {
        error("MBG: group directory ${psm_dir} does not match meta.id ${meta.id}")
    }
    // fragpipe-files.fp-manifest with absolute LC-MS paths (TabWorkflow.java:1241-1250)
    def manifest = writeTextFile('fragpipe-files.fp-manifest', fragpipeManifestText(manifest_rows, '{TASK_DIR}/spectra'))
    def mbg_args = shellArgs(fillArgvTokens(argv, [
        '@TOOLS_DIR@': tools.argv_tools_dir,
        '@THREADS@'  : task.cpus,
        '@PSM@'      : "{TASK_DIR}/${meta.id}/psm.tsv",
        '@MANIFEST@' : '{TASK_DIR}/fragpipe-files.fp-manifest',
    ]))
    def residues = argv[argv.indexOf('--toaddresiduals') + 1]
    // java [bruker] [thermo] -cp <MBG>:<batmass-io>:<IonQuant> com.mbg.MBG ... (CmdMBGMatch.java:127-136), with
    // the Thermo / Bruker libraries next to the MSFragger jar (FragpipeRun.java:1348-1356).
    // FragPipe passes no -Xmx, so MBG gets the JVM default (1/4 of the host RAM); inside a task's
    // memory limit that default is too small for IonQuant's index (OutOfMemoryError on 2 Q Exactive
    // runs at 64 GB), so the heap is the task memory - a resource setting, not an algorithm one.
    """
    export HOME=\$(pwd)
    echo "${tools.mbg_version}" > .mbg_version
    echo "${tools.fragpipe_version}" > .fragpipe_version

    # MBG replaces psm.tsv: run it on a copy of the group directory
    cp -rL "input/${meta.id}" "${meta.id}"

    ${manifest}

    # Every residue/mod to add must parse as a glycan (CmdMBGMatch.verifyResidues)
    java -Duser.home="\$PWD" -cp "${tools.fragpipe_lib}/*" "\$(command -v MbgVerifyResidues.java)" \\
        "${tools.glycan_residues}" "${tools.glycan_mods}" '${residues}'

    ${fragpipeJarLines(tools, ['ionquant', 'msfragger', 'batmass_io'])}
    CMD=( java -Duser.home=\$PWD -Xmx${ram}G ${shellArgs(["-Dlibs.bruker.dir=${tools.bruker_dir}", "-Dlibs.thermo.dir=${tools.thermo_dir}", '-cp', "${tools.mbg_jar}:${tools.batmass_io_jar}:${tools.ionquant_jar}"])} com.mbg.MBG ${mbg_args} ${args} )
    printf '%q ' "\${CMD[@]}" > mbg.log
    printf '\\n' >> mbg.log
    # FragPipe runs MBG from its work dir
    "\${CMD[@]}" 2>&1 | tee -a mbg.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    """
    cp -rL "input/${meta.id}" "${meta.id}"
    ${writeTextFile('fragpipe-files.fp-manifest', fragpipeManifestText(manifest_rows, 'spectra'))}
    # the argv the real run would pass, before run-time tokens are filled
    ${writeTextFile('mbg.log', argv.join('\n') + '\n')}
    echo "${tools.mbg_version}" > .mbg_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}
