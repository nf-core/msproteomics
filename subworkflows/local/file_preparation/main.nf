//
// Decompression, raw file conversion, mzML indexing and optional peak picking for DIA-NN
//

include { THERMORAWFILEPARSER } from '../../../modules/nf-core/thermorawfileparser/main'
include { TDF2MZML            } from '../../../modules/local/tdf2mzml/main'
include { DECOMPRESS          } from '../../../modules/local/decompress/main'
include { MZMLINDEXING } from '../../../modules/local/mzmlindexing/main'
include { OPENMSPEAKPICKER } from '../../../modules/local/openmspeakpicker/main'

workflow FILE_PREPARATION {
    take:
    ch_rawfiles            // channel: [ val(meta), .raw | .mzML | .d | .dia, optionally .gz/.tar/.tar.gz/.zip ]
    val_reindex_mzml       // boolean: index the mzML files (MZMLINDEXING)
    val_openms_peakpicking // boolean: pick peaks in the mzML files (OPENMSPEAKPICKER)
    val_convert_dotd       // boolean: convert Bruker .d to mzML (TDF2MZML)

    main:
    // Decompress archives: X.mzML.gz -> X.mzML, X.raw.gz -> X.raw, X.d.tar.gz -> X.d
    ch_rawfiles
        .branch { _meta, ms_file ->
            compressed: ['.tar', '.zip', '.gz'].any { ext -> hasExtension(ms_file, ext) }
            uncompressed: true
        }
        .set { ch_compression }

    DECOMPRESS(ch_compression.compressed)
    ch_ms_files = ch_compression.uncompressed.mix(DECOMPRESS.out.decompressed_files)

    ch_ms_files
        .branch { _meta, ms_file ->
            raw: hasExtension(ms_file, '.raw')
            mzml: hasExtension(ms_file, '.mzML')
            dotd: hasExtension(ms_file, '.d')
            dia: hasExtension(ms_file, '.dia')
            unsupported: true
        }
        .set { ch_by_format }

    ch_by_format.unsupported.map { _meta, ms_file ->
        error("Unsupported mass spectrometry file for DIA-NN: ${ms_file}. Supported formats: .raw, .mzML, .d (Bruker), .dia, optionally compressed (.gz, .tar, .tar.gz, .zip)")
    }

    // Indexing is done here when requested; search engines index on demand otherwise
    // (peeking into remote files costs a lot of RAM, see https://github.com/bigbio/quantms/issues/61).
    if (val_reindex_mzml) {
        MZMLINDEXING(ch_by_format.mzml)
        ch_mzml = MZMLINDEXING.out.mzmls_indexed
    } else {
        ch_mzml = ch_by_format.mzml
    }

    THERMORAWFILEPARSER(ch_by_format.raw)
    ch_mzml = ch_mzml.mix(THERMORAWFILEPARSER.out.spectra)

    // Peak picking applies to the mzML files only; .d and .dia files pass through.
    if (val_openms_peakpicking) {
        OPENMSPEAKPICKER(ch_mzml)
        ch_mzml = OPENMSPEAKPICKER.out.mzmls_picked
    }

    if (val_convert_dotd) {
        TDF2MZML(ch_by_format.dotd)
        ch_dotd = TDF2MZML.out.mzmls_converted
    } else {
        ch_dotd = ch_by_format.dotd
    }

    emit:
    results = ch_mzml.mix(ch_dotd, ch_by_format.dia) // channel: [ val(meta), .mzML | .d | .dia ]
}

//
// check file extension
//
def hasExtension(file, extension) {
    return file.toString().toLowerCase().endsWith(extension.toLowerCase())
}
