/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_CONVERT SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Optional DIA preprocessing of the LC-MS files FragPipe searches.

    LC-MS files are never converted: MSFragger reads .raw, .mzML, .mzXML, .mgf, .mzBIN and
    Bruker .d natively (CmdMsfragger.java:390-403, 764-766), and so do DIA-Umpire and diaTracer.

    Modules:
    - DIAUMPIRE (DIA preprocessing, generates pseudo-spectra)
    - DIATRACER (Bruker DIA preprocessing)

    Execution: Per-sample parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DIAUMPIRE           } from '../../../modules/local/diaumpire/main'
include { DIATRACER           } from '../../../modules/local/diatracer/main'
include { ddaPseudoSpectraMeta; shouldRunTool; getToolArgs; umpireQTiers } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

workflow FRAGPIPE_CONVERT {
    take:
    ch_file          // channel: [ val(meta), path(raw/mzml/d_file) ] - sample files
    ch_tool_configs  // channel: val(tool_configs_map) - parsed JSON config

    main:

    //
    // Branch input files by type: MSFragger's LC-MS formats (CmdMsfragger.java:764-766)
    //
    ch_files = ch_file
        .branch {
            _meta, file ->
            d:     file.name.toLowerCase().endsWith('.d')
            lcms:  file.name.toLowerCase() ==~ /.+\.(raw|mzml|mzxml|mgf|mzbin)/
            other: true
        }

    // No silent drop: anything else is not an input FragPipe can search
    ch_files.other.subscribe { _meta, file ->
        error("FRAGPIPE_CONVERT: unsupported LC-MS file ${file.name} (FragPipe reads .raw, .mzML, .mzXML, .mgf, .mzBIN and Bruker .d)")
    }

    //
    // STEP 1: DIA-Umpire (optional) on every file that is not DDA, DIA-Quant or DDA+
    // (CmdUmpireSe.java:59-78, :120-123); it reads .raw itself (CmdUmpireSe.java:129-141).
    // FragpipeRun.java:1360-1366 refuses .d files, which workflows/fragpipe.nf checks before
    // the run (checkDiaPreprocessing).
    //
    ch_lcms = ch_files.lcms
        .combine(ch_tool_configs)
        .branch { meta, _mzml, configs ->
            umpire: shouldRunTool(configs, 'diaumpire') && !(meta.data_type in ['DDA', 'DIA-Quant', 'DDA+'])
            other:  true
        }

    // umpire-se.params text (UmpirePanel.collect(), parse_fragpipe_workflow.py), the same for every file
    ch_diaumpire_params = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'diaumpire') }
        .map { configs -> getToolArgs(configs, 'diaumpire') }

    DIAUMPIRE(
        ch_lcms.umpire
            .map { meta, mzml, _configs -> [meta, mzml] }
            .combine(ch_diaumpire_params)
    )

    // The enabled Q tiers replace the DIA file, typed DDA (CmdUmpireSe.java:66-75, FragpipeRun.java:1380-1382).
    ch_umpire_spectra = DIAUMPIRE.out.q1_spectra.map { meta, f -> [meta, 'Q1', f] }
        .mix(DIAUMPIRE.out.q2_spectra.map { meta, f -> [meta, 'Q2', f] })
        .mix(DIAUMPIRE.out.q3_spectra.map { meta, f -> [meta, 'Q3', f] })
        .combine(ch_tool_configs)
        .filter { _meta, tier, _f, configs -> umpireQTiers(configs).contains(tier) }
        .map { meta, _tier, f, configs ->
            [ddaPseudoSpectraMeta(meta, f.baseName, shouldRunTool(configs, 'speclibgen') || shouldRunTool(configs, 'diann')), f]
        }

    //
    // STEP 2: diaTracer (optional) on DIA / DIA-Lib / GPF-DIA .d files (CmdDiaTracer.java:131-145);
    // other .d files are searched as they are.
    //
    ch_d = ch_files.d
        .combine(ch_tool_configs)
        .branch { meta, _d_file, configs ->
            diatracer: shouldRunTool(configs, 'diatracer') && meta.data_type in ['DIA', 'GPF-DIA', 'DIA-Lib']
            other:     true
        }

    DIATRACER(
        ch_d.diatracer.map { meta, d_file, configs -> [meta, d_file, getToolArgs(configs, 'diatracer')] }
    )

    // The _diatracer.mzML replaces the .d, typed DDA (CmdDiaTracer.java:139)
    ch_diatracer_spectra = DIATRACER.out.mzml
        .combine(ch_tool_configs)
        .map { meta, f, configs ->
            [ddaPseudoSpectraMeta(meta, f.baseName, shouldRunTool(configs, 'speclibgen') || shouldRunTool(configs, 'diann')), f]
        }

    //
    // LC-MS files FragPipe searches after DIA preprocessing
    //
    ch_mzml_all = ch_lcms.other.map { meta, mzml, _configs -> [meta, mzml] }
        .mix(ch_umpire_spectra)
        .mix(ch_d.other.map { meta, d_file, _configs -> [meta, d_file] })
        .mix(ch_diatracer_spectra)

    emit:
    mzml        = ch_mzml_all        // channel: [ val(meta), path(spectra) ] - LC-MS files to search (meta.data_type per file)
    umpire_mzml = ch_umpire_spectra  // channel: [ val(meta), path(mzml) ] - DIA-Umpire pseudo-spectra (FragPipe writes them to its work dir)
    versions    = DIAUMPIRE.out.versions_diaumpire
        .mix(DIATRACER.out.versions_diatracer)
}
