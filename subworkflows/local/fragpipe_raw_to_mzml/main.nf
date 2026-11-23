/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_RAW_TO_MZML SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Opt-in conversion of the Thermo .raw files of a FragPipe search to mzML, before the
    manifest is written and before FRAGPIPE / FRAGPIPE_HEADLESS_WF.

    FragPipe itself reads .raw (FRAGPIPE_CONVERT converts nothing, as FragPipe 24.0). A pipeline
    that wants FragPipe to read the ThermoRawFileParser mzML instead calls this subworkflow on
    its LC-MS files and gives `spectra` to FRAGPIPE / FRAGPIPE_HEADLESS_WF, so the manifest names
    the files FragPipe gets. The msproteomics workflows do not call it.

    - Thermo .raw (file name ends with .raw, any case): THERMORAWFILEPARSER -> mzML, same meta
    - Every other file (Bruker .d, mzML, ...): as is (ThermoRawFileParser reads Thermo files only)

    Modules:
    - THERMORAWFILEPARSER (nf-core)

    Execution: Per-file parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { THERMORAWFILEPARSER } from '../../../modules/nf-core/thermorawfileparser/main'

workflow FRAGPIPE_RAW_TO_MZML {
    take:
    ch_files    // channel: [ val(meta), path(spectra) ] - LC-MS files of a FragPipe search

    main:
    ch_branched = ch_files.branch { _meta, spectra ->
        raw:   spectra.name.toLowerCase().endsWith('.raw')
        other: true
    }

    THERMORAWFILEPARSER(ch_branched.raw)

    emit:
    spectra  = ch_branched.other.mix(THERMORAWFILEPARSER.out.spectra)  // channel: [ val(meta), path(spectra) ] - the mzML in place of each .raw
    raw      = ch_branched.raw                                         // channel: [ val(meta), path(raw) ] - the Thermo .raw inputs
    versions = THERMORAWFILEPARSER.out.versions_thermorawfileparser    // channel: [ val(process), val(tool), val(version) ]
}
