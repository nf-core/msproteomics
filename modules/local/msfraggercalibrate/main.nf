/*
 * MSFRAGGERCALIBRATE: MSFragger calibration search (--split1) of one file for a split-database search
 *
 * The MSFragger run of calibrate() in FragPipe 24.0 tools/msfragger_pep_split.py (lines 498-503) for one
 * spectra file and its data type's params set. FragPipe runs the split script per params set; the
 * per-file patched FragPipe used for golden comparison runs it per file, so calibration is per file.
 * Run only when the params have calibrate_mass = 1 or 2 (msfragger_pep_split.py:631-632).
 * MSFRAGGERAPPLYCALIBRATION turns its output into the search params and spectra (lines 524-611).
 *
 * Input:
 *   - meta: Sample metadata map
 *   - spectra: one mass spectrometry data file
 *   - params_text: text of the file's MSFragger params (fragger.params / fragger_dia.params / ...)
 *   - fasta: the protein-sorted FASTA (SORTFASTA, msfragger_pep_split.py sample_fasta)
 *
 * Output:
 *   - calibration: MSFragger console output and the .mzBIN_calibrated it wrote (none when it wrote none)
 *   - uncalibrated_mzml: <stem>_uncalibrated.mzML for .raw / .d input
 *   - log: MSFragger console output (msfragger_calibrate.log)
 *   - versions: Software versions
 */
include { fragpipeJavaMemGb; fragpipeTools; requireFragpipeLicenseAgreement; writeTextFile; setFragpipeParamValues; shellArgs; fragpipeJarLines; rawLinkCopies } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process MSFRAGGERCALIBRATE {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(spectra), val(params_text)
    path fasta

    output:
    tuple val(meta), path("msfragger_calibrate.log"), path("*.mzBIN_calibrated", arity: '0..*'), emit: calibration
    // MSFragger writes <stem>_uncalibrated.mzML next to .raw / .d input (write_uncalibrated_mzml = 1)
    tuple val(meta), path("*_uncalibrated.mzML")                                         , emit: uncalibrated_mzml, optional: true
    // and <stem>_uncalibrated.mzBIN, which a slice-db FragPipe run keeps in raw_files/ next to the input
    tuple val(meta), path("*_uncalibrated.mzBIN")                                        , emit: uncalibrated_mzbin, optional: true
    // Console log: read by FRAGPIPESDRF for the instruments (FragpipeRun.java:586)
    tuple val(meta), path("msfragger_calibrate.log")                                     , emit: log
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT"  , emit: license_agreement
    tuple val("${task.process}"), val('msfragger'), eval("cat .msfragger_version"), emit: versions_msfragger, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba (copied from modules/nf-core/diann/main.nf)
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MSFRAGGERCALIBRATE module does not support Conda: it runs in the FragPipe image given with --fragpipe_container. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def mem = fragpipeJavaMemGb(task)
    def tools = fragpipeTools(task)
    def agreement = requireFragpipeLicenseAgreement(task)
    // The search params pointing at the sorted FASTA (msfragger_pep_split.py:499-501), values replaced in place
    def calibrate_params = setFragpipeParamValues(params_text, [database_name: "{TASK_DIR}/${fasta}", num_threads: task.cpus])
    // argv as msfragger_pep_split.py:502: msfragger_cmd + [params, --license <file>, --split1, files]
    // (msfragger_cmd = "java -jar -Dfile.encoding=UTF-8 -Xmx<ram>G <jar>", CmdMsfragger.java:516-528).
    // MSFragger finds ext/thermo and ext/bruker next to its jar; FragPipe passes no -Dlibs.* (msfragger_pep_split.py:502).
    def argv = ['java', '-Duser.home={TASK_DIR}', '-jar', '-Dfile.encoding=UTF-8', "-Xmx${mem}G", tools.msfragger_jar,
        '{TASK_DIR}/calibrate_fragger.params', '--license', tools.license, '--split1', spectra.name]
    def log_argv = argv.collect { word -> word == tools.license ? '<license>' : word }
    """
    ${rawLinkCopies(spectra)}
    export HOME=\$(pwd)
    ${writeTextFile('I_AGREE_FRAGPIPE_LICENSE_AGREEMENT', agreement)}
    ${fragpipeJarLines(tools, ['msfragger'])}
    echo "${tools.fragpipe_version}" > .fragpipe_version

    ${writeTextFile('calibrate_fragger.params', calibrate_params)}


    CMD=(${shellArgs(argv)} ${args})
    LOG=(${shellArgs(log_argv)} ${args})
    printf '%s\\n' "\${LOG[*]}" > msfragger_calibrate.log
    "\${CMD[@]}" 2>&1 | tee -a msfragger_calibrate.log
    """

    stub:
    def tools = fragpipeTools(task)
    """
    touch "${spectra.baseName}.mzBIN_calibrated" "${spectra.baseName}_uncalibrated.mzBIN"
    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT msfragger_calibrate.log
    echo "${tools.msfragger_version}" > .msfragger_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}
