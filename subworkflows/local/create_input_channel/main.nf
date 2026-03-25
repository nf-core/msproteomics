//
// Create channel for input file
//
include { SDRF_PARSING } from '../../../modules/local/sdrf_parsing/main'


workflow CREATE_INPUT_CHANNEL {
    take:
    ch_sdrf_or_design
    is_sdrf

    main:
    ch_versions = Channel.empty()

    SDRF_PARSING(ch_sdrf_or_design)
    ch_versions = ch_versions.mix(SDRF_PARSING.out.versions)
    ch_config = SDRF_PARSING.out.ch_sdrf_config_file
    ch_expdesign = SDRF_PARSING.out.ch_expdesign

    def Set enzymes = []
    def Set files = []

    // Wrapper holds shared state across map invocations for labelling/acquisition consistency checks
    def wrapper = [
        labelling_type: "",
        acquisition_method: "",
        experiment_id: ch_sdrf_or_design,
    ]

    ch_config
        .splitCsv(header: true, sep: '\t')
        .map { create_meta_channel(it, enzymes, files, wrapper) }
        .branch {
            ch_meta_config_dia: it[0].acquisition_method.contains("dia")
            ch_meta_config_iso: it[0].labelling_type.contains("tmt") || it[0].labelling_type.contains("itraq")
            ch_meta_config_lfq: it[0].labelling_type.contains("label free")
        }
        .set { result }
    ch_meta_config_iso = result.ch_meta_config_iso
    ch_meta_config_lfq = result.ch_meta_config_lfq
    ch_meta_config_dia = result.ch_meta_config_dia

    emit:
    ch_meta_config_iso // [meta, [spectra_files ]]
    ch_meta_config_lfq // [meta, [spectra_files ]]
    ch_meta_config_dia // [meta, [spectra files ]]
    ch_expdesign
    versions = ch_versions
}

// Function to get list of [meta, [ spectra_files ]]
def create_meta_channel(LinkedHashMap row, enzymes, files, wrapper) {
    def meta = [:]
    def filestr

    if (!params.sample_dir) {
        filestr = row.URI.toString()
    }
    else {
        filestr = row.Filename.toString()
    }

    meta.mzml_id = file(filestr).name.take(file(filestr).name.lastIndexOf('.'))
    meta.id = meta.mzml_id  // alias for nf-core module compatibility
    meta.experiment_id = file(wrapper.experiment_id.toString()).baseName

    // apply transformations given by specified sample_dir and type
    if (params.sample_dir) {
        filestr = params.sample_dir + File.separator + filestr
        filestr = (params.local_input_type
            ? filestr.take(filestr.lastIndexOf('.')) + '.' + params.local_input_type
            : filestr)
    }

    // existence check
    if (!file(filestr).exists()) {
        error("ERROR: Please check input file -> File Uri does not exist!\n${filestr}")
    }

    // Parse metadata from SDRF config
    if (row["Proteomics Data Acquisition Method"].toString().toLowerCase().contains("data-dependent acquisition")) {
        meta.acquisition_method = "dda"
    }
    else if (row["Proteomics Data Acquisition Method"].toString().toLowerCase().contains("data-independent acquisition")) {
        meta.acquisition_method = "dia"
    }
    else {
        error("Currently DIA and DDA are supported for the pipeline. Check and Fix your SDRF.")
    }

    // dissociation method conversion
    if (row.DissociationMethod == "COLLISION-INDUCED DISSOCIATION") {
        meta.dissociationmethod = "CID"
    }
    else if (row.DissociationMethod == "HIGHER ENERGY BEAM-TYPE COLLISION-INDUCED DISSOCIATION") {
        meta.dissociationmethod = "HCD"
    }
    else if (row.DissociationMethod == "ELECTRON TRANSFER DISSOCIATION") {
        meta.dissociationmethod = "ETD"
    }
    else if (row.DissociationMethod == "ELECTRON CAPTURE DISSOCIATION") {
        meta.dissociationmethod = "ECD"
    }
    else {
        meta.dissociationmethod = row.DissociationMethod
    }

    wrapper.acquisition_method = meta.acquisition_method
    meta.labelling_type = row.Label
    meta.fixedmodifications = row.FixedModifications
    meta.variablemodifications = row.VariableModifications
    meta.precursormasstolerance = Double.parseDouble(row.PrecursorMassTolerance)
    meta.precursormasstoleranceunit = row.PrecursorMassToleranceUnit
    meta.fragmentmasstolerance = Double.parseDouble(row.FragmentMassTolerance)
    meta.fragmentmasstoleranceunit = row.FragmentMassToleranceUnit
    meta.enzyme = row.Enzyme

    enzymes += row.Enzyme
    if (enzymes.size() > 1) {
        error("Currently only one enzyme is supported for the whole experiment. Specified was '${enzymes}'. Check or split your SDRF.\n${filestr}")
    }
    // Nothing to determine for dia. Only LFQ allowed there.
    if (!meta.acquisition_method.equals("dia")) {
        if (wrapper.labelling_type.equals("")) {
            if (meta.labelling_type.contains("tmt") || meta.labelling_type.contains("itraq") || meta.labelling_type.contains("label free")) {
                wrapper.labelling_type = meta.labelling_type
            }
            else {
                error("Unsupported quantification type '${meta.labelling_type}'.")
            }
        }
        else {
            if (meta.labelling_type != wrapper.labelling_type) {
                error("Currently, only one label type per design is supported: was '${wrapper.labelling_type}', now is '${meta.labelling_type}'.")
            }
        }
    }

    if (wrapper.labelling_type.contains("label free") || meta.acquisition_method == "dia") {
        if (filestr in files) {
            error("Currently only one search engine setting/DIA-NN setting per file is supported for the whole experiment. ${filestr} has multiple entries in your SDRF. Maybe you have a (isobaric) labelled experiment? Otherwise, consider splitting your design into multiple experiments.")
        }
        files += filestr
    }

    return [meta, filestr]
}
