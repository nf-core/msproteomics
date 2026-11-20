/*
 * CRYSTALC: Crystal-C chimeric-spectrum correction of one MSFragger pepXML (CmdCrystalc)
 *
 * Params file as CmdCrystalc.java:183-196: the jar's crystalc.params algorithm settings
 * (config_file, from CrystalcPanel.toParams) plus the runtime keys fasta, raw_file_location,
 * raw_file_extension (the LC-MS file's own extension) and output_location.
 * argv as CmdCrystalc.java:205-215:
 *   java -Dlibs.thermo.dir=<MSFragger dir>/ext/thermo -Xmx<ram>G -cp <crystalc>:<batmass-io>:<grppr> crystalc.Run <params> <pepXML>
 *
 * Input:
 *   - meta: Sample metadata map (crystalc_index: index of the pepXML among all Crystal-C runs of the run, names the params file)
 *   - pepxml_file: MSFragger pepXML
 *   - lcms_file: The LC-MS file the pepXML was searched from (.mzML, .mzXML or .raw)
 *   - config_file: Crystal-C algorithm settings (FragPipe's plan, bin/plan_to_tool_configs.py)
 *   - fasta: Protein database
 *
 * Output:
 *   - results_dir / pepxml_filtered: <pepXML base>_c.pepXML
 *   - versions: Crystal-C and FragPipe versions (topic: versions)
 */
include { fragpipeTools; crystalcParamsName; crystalcRuntimeParams; heredocText; fragpipeJarLines; shellArgs; rawLinkCopies } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { javaXmxGb } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process CRYSTALC {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(pepxml_file), path(lcms_file), path(config_file)
    path(fasta)  // protein database FASTA (shared resource, required by CrystalC)

    output:
    tuple val(meta), path("${prefix}")                  , emit: results_dir
    tuple val(meta), path("${prefix}/*_c.pepXML")       , emit: pepxml_filtered
    tuple val("${task.process}"), val('crystalc'), eval("cat .crystalc_version"), emit: versions_crystalc, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba (copied from modules/nf-core/diann/main.nf)
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CRYSTALC module does not support Conda: it runs in the FragPipe image given with --fragpipe_container. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = javaXmxGb(task, 0)
    def tools = fragpipeTools(task)
    // <ccParamsFilePrefix>-<index>-<pepXML name>.params; fragpipe_search counts the index over the whole run as CmdCrystalc.java:176-185
    def params_file = crystalcParamsName(meta, pepxml_file)
    def runtime_params = heredocText(crystalcRuntimeParams(task.cpus, fasta.name, lcms_file.name, prefix))
    // FragPipe 24.0 class path (CmdCrystalc.JAR_CRYSTALC_NAME, JAR_DEPS) and the Thermo libraries next to the
    // MSFragger jar (FragpipeRun.java:1348-1356, CmdCrystalc.java:205-208), all in the FragPipe image
    """
    ${rawLinkCopies(lcms_file)}
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    echo "${tools.fragpipe_version}" > .fragpipe_version
    echo "${tools.crystalc_version}" > .crystalc_version

    # Algorithm settings, then the runtime keys (CmdCrystalc.java:190-196)
    cp ${config_file} ${params_file}
    cat >> ${params_file} <<FP_EOF
${runtime_params}
FP_EOF

    ${fragpipeJarLines(tools, ['msfragger', 'batmass_io'])}

    CMD=(java -Duser.home=\$PWD ${shellArgs(["-Dlibs.thermo.dir=${tools.thermo_dir}"])} -Xmx${ram}G -cp ${shellArgs(["${tools.crystalc_jar}:${tools.batmass_io_jar}:${tools.grppr_jar}"])} crystalc.Run ${params_file} ${pepxml_file} ${args})
    printf '%s\\n' "\${CMD[*]}" > ${prefix}/crystalc.log
    "\${CMD[@]}" 2>&1 | tee -a ${prefix}/crystalc.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def basename = pepxml_file.baseName
    def tools = fragpipeTools(task)
    """
    mkdir -p ${prefix}
    touch ${prefix}/${basename}_c.pepXML
    touch ${prefix}/crystalc.log
    echo "${tools.crystalc_version}" > .crystalc_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}
