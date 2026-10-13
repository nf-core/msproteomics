/*
 * MSBOOSTER: MSBooster rescoring features (CmdMSBooster)
 *
 * msbooster_params.txt holds the workflow-derived keys from PARSEFRAGPIPEWORKFLOW plus the
 * runtime keys CmdMSBooster.java:160-212 writes: DiaNN (FragPipe's bundled
 * DIA-NN, Diann.fallbackDiannPath), fragger, unimodObo, mzmlDirectory, pinPepXMLDirectory.
 * useIM and findBestImModel are on only when the run has timsTOF data (.d or _diatracer.mzML;
 * CmdMSBooster.java:151-157, 170, 178). msboosterParamsText renders the file. Diann.getLDPRELOAD sets no
 * LD_PRELOAD when DIA-NN runs as it is, which it does in the FragPipe image.
 * argv: java -Djava.awt.headless=true -Xmx<ram>G -cp MSBooster-1.4.17.jar:batmass-io-1.36.6.jar mainsteps.MainClass --paramsList <file>
 *
 * Input:
 *   - meta: Sample (or, when MSBooster pools runs, run) metadata
 *   - pin_files: MSFragger PIN file(s)
 *   - mzml_files: FragPipe's LC-MS directory content for those PINs (LC-MS files plus the
 *     MSFragger mzML outputs); MSBooster matches <base>.mzML / <base>_uncalibrated.mzML / .mgf itself
 *   - params_text: workflow-derived msbooster params (PARSEFRAGPIPEWORKFLOW msbooster.config text)
 *   - fragger_params: the MSFragger params file MSBooster reads (CmdMSBooster.java:134-149)
 *   - library: the file of msbooster.spectral-library-path (spectraPredFile/RTPredFile/IMPredFile; CmdMSBooster.java:180-184), staged as spectral_library/<name>, [] without one
 *   - has_ion_mobility: true when any LC-MS file of the run is timsTOF data
 *   - threads: numThreads of params_text (workflow.threads / --fragpipe_threads); the task asks for that many cores
 *   - pin_dirs: FragPipe group of each PIN (same order as pin_files; '' for the blank group); MSBooster is given <group>/<pin>
 *
 * Output:
 *   - results_dir / pin_edited: <base>_edited.pin, the rescored PIN MSBooster writes
 *   - predictions: MSBooster/ (spectraRT.tsv, spectraRT_full.tsv, spectraRT.predicted.bin for FragPipe-PDV,
 *     MSBooster_plots/), under the PINs' deepest common directory (msboosterOutputDir)
 *   - versions: MSBooster and FragPipe versions (topic: versions)
 */
include { fragpipeTools; msboosterParamsText; msboosterPinPaths; msboosterOutputDir; msboosterPinStaging; writeTextFile; fragpipeJarLines; shellArgs } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process MSBOOSTER {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(pin_files), path(mzml_files)
    val params_text
    path fragger_params
    path library, stageAs: 'spectral_library/*'
    val has_ion_mobility
    val threads
    val pin_dirs

    output:
    tuple val(meta), path("${prefix}")              , emit: results_dir
    tuple val(meta), path("${prefix}/*.pin")        , emit: pin_edited
    tuple val(meta), path("${booster_dir}")         , emit: predictions
    tuple val("${task.process}"), val('msbooster'), eval("cat .msbooster_version"), emit: versions_msbooster, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? task.memory.toGiga() : 8
    def pins = pin_files instanceof List ? pin_files : [pin_files]
    def pin_paths = msboosterPinPaths(pins.collect { p -> p.name }, pin_dirs)
    booster_dir = msboosterOutputDir(pin_paths)
    def tools = fragpipeTools(task)
    def msbooster_params = msboosterParamsText(params_text, has_ion_mobility, threads, tools, fragger_params.name, pin_paths, library)
    // MSBooster and batmass-io: the jars the image's FragPipe build names (FragPipe 24.0 release: MSBooster-1.4.14,
    // batmass-io-1.36.5; FragPipePlus: 1.4.17, 1.36.6; fragpipeJarLines)
    // renamePin = 1 writes <base>_edited.pin next to each PIN. It keeps that name: Percolator reads
    // <base>_edited.pin and PercolatorOutputToPepXML reads <base>.pin (CmdPercolator.java:209-213, :292).
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    ${fragpipeJarLines(tools, ['msbooster', 'batmass_io'])}
    echo "${tools.fragpipe_version}" > .fragpipe_version

    ${msboosterPinStaging(pin_paths)}
    ${writeTextFile("${prefix}/msbooster_params.txt", msbooster_params)}

    CMD=(java -Duser.home=\$PWD -Djava.awt.headless=true -Xmx${ram}G -cp ${shellArgs(["${tools.msbooster_jar}:${tools.batmass_io_jar}"])} mainsteps.MainClass --paramsList \$(pwd)/${prefix}/msbooster_params.txt ${args})
    printf '%s\\n' "\${CMD[*]}" > ${prefix}/${prefix}_msbooster.log
    "\${CMD[@]}" 2>&1 | tee -a ${prefix}/${prefix}_msbooster.log

    mv ${pin_paths.collect { p -> "'${p.replaceAll(/\.pin$/, '')}_edited.pin'" }.join(' ')} ${prefix}/
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def pins = pin_files instanceof List ? pin_files : [pin_files]
    def pin_paths = msboosterPinPaths(pins.collect { p -> p.name }, pin_dirs)
    booster_dir = msboosterOutputDir(pin_paths)
    def tools = fragpipeTools(task)
    """
    mkdir -p ${prefix} '${booster_dir}/MSBooster_plots'
    touch '${booster_dir}/spectraRT.tsv' '${booster_dir}/spectraRT_full.tsv' '${booster_dir}/spectraRT.predicted.bin'
    # the real run's params with the task directory as <task dir>, so the stub output does not depend on the work dir
    ${writeTextFile("${prefix}/msbooster_params.txt", msboosterParamsText(params_text, has_ion_mobility, threads, tools, fragger_params.name, pin_paths, library).replace('{TASK_DIR}', '<task dir>'))}
    touch ${pins.collect { p -> "'${prefix}/${p.baseName}_edited.pin'" }.join(' ')}
    touch ${prefix}/${prefix}_msbooster.log
    echo "${tools.msbooster_version}" > .msbooster_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}
