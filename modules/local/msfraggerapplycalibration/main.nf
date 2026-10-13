/*
 * MSFRAGGERAPPLYCALIBRATION: search params and spectra of a split-database search after MSFragger's
 * calibration search
 *
 * Runs parse_calibration_output.py, the port of calibrate() in FragPipe 24.0 tools/msfragger_pep_split.py
 * (lines 524-611): the optimized values MSFragger printed are applied to the search params (calibrate_mass
 * = 2 only), check_spectral_files is set to 0, and the spectra file is replaced by its .mzBIN_calibrated
 * (or copied when MSFragger wrote none).
 *
 * Input:
 *   - meta: Sample metadata map
 *   - spectra: the spectra file MSFRAGGERCALIBRATE calibrated
 *   - calibrate_log: MSFRAGGERCALIBRATE console output
 *   - calibrated_mzbin: the <stem>.mzBIN_calibrated MSFRAGGERCALIBRATE wrote ([] when none)
 *   - params_file: the file's MSFragger search params
 *
 * Output:
 *   - calibrated: the calibrated spectra and the search params after calibration
 *   - versions: Python version (topic: versions)
 */
process MSFRAGGERAPPLYCALIBRATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    tuple val(meta), path(spectra), path(calibrate_log), path(calibrated_mzbin), path(params_file)

    output:
    tuple val(meta), path("calibrated_spectra/*"), path("calibrated/${params_file.name}"), emit: calibrated
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p calibrated
    parse_calibration_output.py \\
        --log_file ${calibrate_log} \\
        --params_file ${params_file} \\
        --output calibrated/${params_file.name} \\
        --spectra_files ${spectra} \\
        --output_dir calibrated_spectra ${args}
    """

    stub:
    """
    mkdir -p calibrated calibrated_spectra
    touch "calibrated_spectra/${spectra.baseName}.mzBIN_calibrated"
    cp ${params_file} calibrated/${params_file.name}
    """
}
