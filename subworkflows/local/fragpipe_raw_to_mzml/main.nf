/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_RAW_TO_MZML SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The one conversion of Thermo .raw files to mzML for every FragPipe route: FragPipe never gets a
    .raw. Thermo RawFileReader under Mono misreads a staged symlink ("No open raw file", 0 scans),
    and one conversion gives every route the same mzML.

    A caller runs it on its LC-MS files before it writes the manifest and gives `spectra` to
    FRAGPIPE / FRAGPIPE_HEADLESS_WF, so the manifest names the files FragPipe searches (FRAGPIPEHEADLESS
    checks that each manifest row names a staged file). workflows/fragpipe.nf (both modes and the TMT
    label check) calls it; tests/python/test_fragpipe_raw_conversion.py fails if a route does not.

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
