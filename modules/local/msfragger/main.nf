/*
 * MSFRAGGER: MSFragger database search of one or more LC-MS files (CmdMsfragger)
 *
 * Input:
 *   - meta: Sample metadata map
 *   - mzml_files: LC-MS files
 *   - params_name / params_text: FragPipe's params file name and text for the data type
 *   - fasta: protein database
 *
 * Output:
 *   - pepxml / pin / tsv: search results
 *   - lcms_files / mzbin_all / pepindex: files MSFragger writes next to its inputs
 *   - results_dir / log: results directory and msfragger.log
 *   - license_agreement: I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
 *   - versions_msfragger / versions_fragpipe: versions (topic: versions)
 */
include { fragpipeJavaMemGb; msfraggerPartialArgs; fragpipeTools; requireFragpipeLicenseAgreement; writeTextFile; setFragpipeParamValues; shellArgs; msfraggerStubFiles; msfraggerTrailingArgv; fragpipeJarLines; rawLinkCopies } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process MSFRAGGER {
    tag "${meta.id}"
    label 'process_high'

    input:
    // params_name: FragPipe's name of the file's MSFragger params set (fragger.params / fragger_dda.params /
    //   fragger_dia.params / fragger_gpfdia.params / fragger_dda_plus.params, CmdMsfragger.java:444-447, 541-552)
    // params_text: text of that params set
    tuple val(meta), path(mzml_files), val(params_name), val(params_text), path(fasta)

    output:
    // Search outputs
    tuple val(meta), path("${prefix}/*.pepXML")              , emit: pepxml, optional: true
    tuple val(meta), path("${prefix}/*.pin")                 , emit: pin, optional: true
    tuple val(meta), path("${prefix}/*.tsv")                 , emit: tsv, optional: true

    // Files MSFragger writes next to the LC-MS file (FragPipe's LC-MS directory after the search):
    // <base>_uncalibrated.mzML (.raw/.d input), <base>_calibrated.mzML (write_calibrated_mzml = 1),
    // <base>_uncalibrated.mzBIN and <base>.mzBIN_calibrated. Downstream tools read them from there.
    tuple val(meta), path("{*_uncalibrated.mzML,*_calibrated.mzML,*_uncalibrated.mzBIN,*.mzBIN_calibrated}"), emit: lcms_files, optional: true

    // <base>.mzBIN_all next to the LC-MS file with write_mzbin_all = 1 (Write sub mzML, FragpipeRun.java:1481);
    // only WriteSubMzml reads it (WriteSubMzml.java:73), so it is not part of the LC-MS directory above
    tuple val(meta), path("*.mzBIN_all")                     , emit: mzbin_all, optional: true

    // Peptide index MSFragger builds next to the FASTA during the search, one <fasta>.<N>.pepindex
    // per distinct digestion parameter set of the call (with calibrate_mass 1/2: the first-search
    // set, then the main-search set). As in FragPipe, there is no separate index step.
    tuple val(meta), path("*.pepindex")                      , emit: pepindex, optional: true

    // Results directory and logs
    tuple val(meta), path("${prefix}")                       , emit: results_dir
    tuple val(meta), path("${prefix}/*.log")                 , emit: log, optional: true

    // License agreement
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT", emit: license_agreement

    // Version tracking
    tuple val("${task.process}"), val('msfragger'), eval("cat .msfragger_version"), emit: versions_msfragger, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba (copied from modules/nf-core/diann/main.nf)
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MSFRAGGER module does not support Conda: it runs in the FragPipe image given with --fragpipe_container. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def mem = fragpipeJavaMemGb(task)
    def tools = fragpipeTools(task)
    def agreement = requireFragpipeLicenseAgreement(task)
    def spectra = mzml_files instanceof List ? mzml_files : [mzml_files]
    // FragPipe fills database_name and num_threads from its run inputs (CmdMsfragger.java:441,
    // TabMsfragger.java:1657), in place, so the file keeps FragPipe's layout (Props.java:240-273)
    def fragger_params = setFragpipeParamValues(params_text, [database_name: "{TASK_DIR}/${fasta}", num_threads: task.cpus])
    // argv as CmdMsfragger.java:516-560: java -jar -Dfile.encoding=UTF-8 -Xmx<ram>G <jar> <params> [--license <file>] <lcms files>;
    // MSFragger finds ext/thermo and ext/bruker next to its jar, FragPipe passes no -Dlibs.* (CmdMsfragger.java:390-403).
    // The license is the one FragPipe finds in <tools>/.. (FragpipeLocations.locateLicense); the log shows <license>.
    def partial = msfraggerPartialArgs(meta)
    // --no-ptm-score-in-pin after the LC-MS files with Write sub mzML or in its second pass (CmdMsfragger.java:586, :612-630)
    def argv = ['java', '-Duser.home={TASK_DIR}', '-jar', '-Dfile.encoding=UTF-8', "-Xmx${mem}G", tools.msfragger_jar,
        "{TASK_DIR}/${prefix}/${params_name}", '--license', tools.license] + spectra.collect { f -> f.name } + msfraggerTrailingArgv(params_text, meta) + partial
    def log_argv = argv.collect { word -> word == tools.license ? '<license>' : word }
    """
    ${rawLinkCopies(spectra)}
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    ${writeTextFile('I_AGREE_FRAGPIPE_LICENSE_AGREEMENT', agreement)}
    ${fragpipeJarLines(tools, ['msfragger'])}
    echo "${tools.fragpipe_version}" > .fragpipe_version


    ${writeTextFile("${prefix}/${params_name}", fragger_params)}

    CMD=(${shellArgs(argv)} ${args})
    LOG=(${shellArgs(log_argv)} ${args})
    printf '%s\\n' "\${LOG[*]}" > ${prefix}/msfragger.log
    "\${CMD[@]}" 2>&1 | tee -a ${prefix}/msfragger.log

    shopt -s nullglob
    for f in *.pepXML *.pin *.tsv; do mv "\$f" "${prefix}/"; done
    shopt -u nullglob
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    def spectra = mzml_files instanceof List ? mzml_files : [mzml_files]
    """
    ${rawLinkCopies(spectra)}
    mkdir -p ${prefix}
    ${writeTextFile("${prefix}/${params_name}", params_text)}
    touch ${msfraggerStubFiles(spectra, fasta.name, params_text, prefix).collect { f -> "'${f}'" }.join(' ')}
    touch ${prefix}/msfragger.log
    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
    echo "${tools.msfragger_version}" > .msfragger_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}
