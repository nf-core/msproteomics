/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_CONVERT SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Converts raw files to mzML and optionally preprocesses DIA data.

    Modules:
    - THERMORAWFILEPARSER (RAW to mzML)
    - DIAUMPIRE (DIA preprocessing, generates pseudo-spectra)
    - DIATRACER (Bruker DIA preprocessing)

    Execution: Per-sample parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { THERMORAWFILEPARSER } from '../../../modules/nf-core/thermorawfileparser/main'
include { DIAUMPIRE           } from '../../../modules/local/diaumpire/main'
include { DIATRACER           } from '../../../modules/local/diatracer/main'
include { shouldRunTool; getToolArgs } from '../fragpipe_utils'

workflow FRAGPIPE_CONVERT {
    take:
    ch_file          // channel: [ val(meta), path(raw/mzml/d_file) ] - sample files
    ch_tool_configs  // channel: val(tool_configs_map) - parsed JSON config

    main:

    //
    // Branch input files by type
    //
    ch_files = ch_file
        .branch {
            _meta, file ->
            raw:   file.name.toLowerCase().endsWith('.raw')
            mzml:  file.name.toLowerCase().endsWith('.mzml')
            d:     file.name.toLowerCase().endsWith('.d')
            other: true
        }

    //
    // STEP 1: Convert RAW files to mzML using ThermoRawFileParser
    //
    THERMORAWFILEPARSER(ch_files.raw)

    // Combine converted mzML with already-mzML files
    ch_mzml_from_raw = THERMORAWFILEPARSER.out.spectra
    ch_mzml_direct = ch_files.mzml

    //
    // STEP 2: DIA-Umpire for DIA preprocessing (optional)
    // Generates pseudo-spectra (Q1, Q2, Q3) from DIA data
    //
    ch_diaumpire_pre = ch_mzml_from_raw
        .mix(ch_mzml_direct)
        .combine(ch_tool_configs)
        .filter { meta, mzml, configs -> shouldRunTool(configs, 'diaumpire') }
        .map { meta, mzml, configs ->
            def args = getToolArgs(configs, 'diaumpire')
            [meta, mzml, args ?: '# DIA-Umpire default params']
        }

    // Create config files via collectFile (avoids writing files in workflow scope)
    ch_diaumpire_configs = ch_diaumpire_pre
        .map { meta, _mzml, config_content ->
            ["diaumpire_${meta.id}.params", config_content]
        }
        .collectFile()
        .map { config_file ->
            def id = config_file.name.replaceAll(/^diaumpire_/, '').replaceAll(/\.params$/, '')
            [id, config_file]
        }

    ch_for_diaumpire = ch_diaumpire_pre
        .map { meta, mzml, _config_content -> [meta.id, meta, mzml] }
        .join(ch_diaumpire_configs)
        .map { key, meta, mzml, config_file -> [meta, mzml, config_file] }

    DIAUMPIRE(ch_for_diaumpire)

    //
    // STEP 3: diaTracer for Bruker DIA preprocessing (optional)
    //
    ch_for_diatracer = ch_files.d
        .combine(ch_tool_configs)
        .filter { meta, d_file, configs -> shouldRunTool(configs, 'diatracer') }
        .map { meta, d_file, configs ->
            def args = getToolArgs(configs, 'diatracer')
            [meta, d_file, args]
        }

    DIATRACER(ch_for_diatracer, [])

    //
    // .d files that don't go through diaTracer pass through directly
    // (e.g., DDA workflows with Bruker .d input - MSFragger handles .d natively)
    //
    ch_d_passthrough = ch_files.d
        .combine(ch_tool_configs)
        .filter { _meta, _d_file, configs -> !shouldRunTool(configs, 'diatracer') }
        .map { meta, d_file, _configs -> [meta, d_file] }

    //
    // Collect all spectra outputs
    // Includes: converted mzML, direct mzML, diaTracer output, .d passthrough
    //
    ch_mzml_all = ch_mzml_from_raw
        .mix(ch_mzml_direct)
        .mix(DIATRACER.out.mzml)
        .mix(ch_d_passthrough)

    // Pseudo-spectra from DIA-Umpire (optional)
    ch_pseudo_spectra = DIAUMPIRE.out.q1_spectra
        .mix(DIAUMPIRE.out.q2_spectra)
        .mix(DIAUMPIRE.out.q3_spectra)

    emit:
    mzml           = ch_mzml_all           // channel: [ val(meta), path(spectra) ] - mzML or .d files
    pseudo_spectra = ch_pseudo_spectra     // channel: [ val(meta), path(q*_mzml) ] - DIA-Umpire output
    versions       = THERMORAWFILEPARSER.out.versions_thermorawfileparser
        .mix(DIAUMPIRE.out.versions_diaumpire)
        .mix(DIATRACER.out.versions_diatracer)
}
