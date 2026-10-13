/*
 * DIAUMPIRE: DIA-Umpire SE pseudo-spectra from one DIA LC-MS file (CmdUmpireSe)
 *
 * Input:
 *   - meta: Sample metadata map
 *   - dia_file: DIA LC-MS file (mzML, mzXML or raw)
 *   - umpire_params: umpire-se.params text with Thread = THREADS
 *
 * Output:
 *   - results_dir: <run>_Q1/Q2/Q3 mzML (and mgf), umpire-se.params, diaumpire.log
 *   - q1_spectra / q2_spectra / q3_spectra / mgf_spectra: the pseudo-spectra
 *   - versions_diaumpire / versions_fragpipe: versions (topic: versions)
 */
include { fragpipeTools; setFragpipeParamValues; writeTextFile; fragpipeJarLines; shellArgs; rawLinkCopies } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { javaXmxGb } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process DIAUMPIRE {
    tag "${meta.id}"
    label 'process_high'

    input:
    // DIA-Umpire generates pseudo-spectra from DIA data
    // dia_file: DIA mass spec file (mzML, mzXML, or raw)
    // umpire_params: umpire-se.params text with Thread = THREADS (parse_fragpipe_workflow.py generate_diaumpire_config)
    tuple val(meta), path(dia_file), val(umpire_params)

    output:
    tuple val(meta), path("${prefix}")                          , emit: results_dir
    tuple val(meta), path("${prefix}/*_Q1.mzML")                , emit: q1_spectra, optional: true
    tuple val(meta), path("${prefix}/*_Q2.mzML")                , emit: q2_spectra, optional: true
    tuple val(meta), path("${prefix}/*_Q3.mzML")                , emit: q3_spectra, optional: true
    tuple val(meta), path("${prefix}/*_Q*.mgf")                 , emit: mgf_spectra, optional: true
    // DIA-Umpire's own log, written in its cwd; FragPipe moves it with the Q files into the work dir
    // (CmdUmpireSe.java:151-157, UmpireSeGarbageFiles.java:30-32) and keeps it
    tuple val(meta), path("diaumpire_se.log")                   , emit: se_log
    tuple val("${task.process}"), val('diaumpire'), eval("cat .diaumpire_version"), emit: versions_diaumpire, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba (copied from modules/nf-core/diann/main.nf)
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DIAUMPIRE module does not support Conda: it runs in the FragPipe image given with --fragpipe_container. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = javaXmxGb(task, 0)
    def tools = fragpipeTools(task)
    // FragPipe's template with the workflow values; Thread is the run's thread count (UmpirePanel.java:427)
    def params_text = setFragpipeParamValues(umpire_params, [Thread: task.cpus])
    // DIA-Umpire SE generates pseudo-spectra from DIA data; it reads .raw itself, with the Thermo libraries next to
    // the MSFragger jar (FragpipeRun.java:1348-1356). Argv as CmdUmpireSe.java:130-141:
    // java -Xmx<ram>G -Dlibs.thermo.dir=<dir> -cp <jars> dia_umpire_se.DIA_Umpire_SE <lcms file> <params file>
    // (batmass-io-1.36.6 of the FragPipe 24.0 release jar; the tag-24.0 source says 1.36.5)
    """
    ${rawLinkCopies(dia_file)}
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    echo "${tools.diaumpire_version}" > .diaumpire_version
    echo "${tools.fragpipe_version}" > .fragpipe_version

    ${writeTextFile('umpire-se.params', params_text)}

    ${fragpipeJarLines(tools, ['msfragger', 'batmass_io'])}

    CMD=(java -Duser.home=\$PWD -Xmx${ram}G ${shellArgs(["-Dlibs.thermo.dir=${tools.thermo_dir}"])} ${args} -cp ${shellArgs(["${tools.diaumpire_jar}:${tools.batmass_io_jar}"])} dia_umpire_se.DIA_Umpire_SE ${dia_file} umpire-se.params)
    printf '%s\\n' "\${CMD[*]}" > ${prefix}/diaumpire.log
    "\${CMD[@]}" 2>&1 | tee -a ${prefix}/diaumpire.log

    # Move output files to prefix directory
    shopt -s nullglob
    for f in *_Q1.mzML *_Q2.mzML *_Q3.mzML *_Q1.mgf *_Q2.mgf *_Q3.mgf; do cp -L "\$f" ${prefix}/; done
    shopt -u nullglob
    # FragPipe writes umpire-se.params into its work dir (CmdUmpireSe.java:100), so it is published with the results
    cp umpire-se.params ${prefix}/
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def basename = dia_file.baseName
    def tools = fragpipeTools(task)
    """
    mkdir -p ${prefix}
    touch ${prefix}/${basename}_Q1.mzML
    touch ${prefix}/${basename}_Q2.mzML
    touch ${prefix}/${basename}_Q3.mzML
    ${writeTextFile("${prefix}/umpire-se.params", setFragpipeParamValues(umpire_params, [Thread: task.cpus]))}
    touch ${prefix}/diaumpire.log diaumpire_se.log
    echo "${tools.diaumpire_version}" > .diaumpire_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}
