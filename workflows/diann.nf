/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DIA WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { generateSdrfBookkeeping; resolveProteome; resolveOrganism; resolveDatabase; resolveDiannSpeclib } from '../subworkflows/local/utils_nfcore_msproteomics_pipeline'

// SDRF generation (bookkeeping only)
include { GENERATESDRFFROMSAMPLESHEET } from '../modules/local/generatesdrffromsamplesheet/main'

// DIA workflow modules
include { FILE_PREPARATION       } from '../subworkflows/local/file_preparation/main'
include { DIA_PROTEOMICS_ANALYSIS } from '../subworkflows/nf-core/dia_proteomics_analysis/main'
include { MSSTATSLFQ as MSSTATS } from '../modules/local/msstatslfq/main'
include { IQ                     } from '../modules/local/iq/main'
include { PMULTIQC               } from '../modules/local/pmultiqc/main'

// Samplesheet parsing (shared with the FragPipe workflows)
include { parseSamplesheet       } from '../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MSPROTEOMICS_DIANN {

    main:
    // DIA-NN runs, so FragPipe's getGroup() rule puts every file in the single group "".
    ch_input = channel.fromList(parseSamplesheet(params.input, "${projectDir}/assets/schema_input.json", true))

    // iProteome (--proteome): fasta and hosted DIA-NN library; --database / --diann_speclib win
    def proteome = resolveProteome(params)
    def database = resolveDatabase(proteome, params)

    GENERATESDRFFROMSAMPLESHEET(
        channel.fromPath(params.input, checkIfExists: true),
        generateSdrfBookkeeping(params, resolveOrganism(proteome, params)),
        file(database, glob: false, checkIfExists: true)
    )

    // Convert raw files if needed
    FILE_PREPARATION(ch_input, params.reindex_mzml, params.openms_peakpicking, params.convert_dotd)

    // DIA-NN input: one item per MS file with the search settings (precursor = MS1, fragment = MS2 tolerance)
    ch_dia_input = FILE_PREPARATION.out.results.map { meta, ms_file ->
        [
            meta,
            ms_file,
            params.enzyme,
            params.fixed_mods ?: '',
            params.variable_mods ?: '',
            params.precursor_mass_tolerance,
            params.fragment_mass_tolerance,
            params.precursor_mass_tolerance_unit,
            params.fragment_mass_tolerance_unit
        ]
    }

    // glob: false keeps the query string of a UniProt REST URL (fromPath treats '?' as a glob; nf-core/msproteomics e5ce83c)
    ch_searchdb = channel.fromPath(database, glob: false, checkIfExists: true).map { fasta -> [[id: fasta.baseName], fasta] }

    // Experimental design in the quantms two-section format read by QUANTMSUTILS_DIANN2MZTAB:
    // file table (Fraction_Group, Fraction, Spectra_Filepath, Label, Sample), an empty line,
    // then the sample table (Sample, MSstats_Condition, MSstats_BioReplicate).
    // Spectra_Filepath is the name of the file DIA-NN reads (after FILE_PREPARATION), Sample the sample's number;
    // condition defaults to the sample name (parseSamplesheet), bioreplicate to 1.
    ch_expdesign = FILE_PREPARATION.out.results
        .map { meta, ms_file -> [meta, ms_file.name] }
        .toSortedList { a, b -> a[1] <=> b[1] }
        .map { rows ->
            // Sample is the quantms study-variable number (1..n in file order): diann2mztab reads it as an
            // integer (quantmsutils diann2mztab.py convert_to_mztab: index_ref["study_variable"].astype("int"))
            def sample_index = rows.collect { meta, _name -> meta.sample }.unique().withIndex().collectEntries { name, i -> [name, i + 1] }
            def file_lines   = rows.collect { meta, name -> "1\t${meta.fraction}\t${name}\t1\t${sample_index[meta.sample]}" }
            def sample_lines = rows.collect { meta, _name -> "${sample_index[meta.sample]}\t${meta.condition}\t${meta.bioreplicate ?: 1}" }.unique()
            ((["Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample"] + file_lines + [""] +
                ["Sample\tMSstats_Condition\tMSstats_BioReplicate"] + sample_lines).join('\n') + '\n')
        }
        .collectFile(name: 'experimental_design.tsv')
        .map { expdesign -> [[id: 'experiment'], expdesign] }

    // Pre-built in-silico spectral library (skips QUANTMSUTILS_DIANNCFG and in-silico generation):
    // --diann_speclib > the iProteome's hosted library > built by DIANN_INSILICOLIBRARYGENERATION
    def diann_speclib = resolveDiannSpeclib(proteome, params)
    ch_speclib = diann_speclib
        ? channel.of(file(diann_speclib, checkIfExists: true))
        : channel.empty()

    // Pre-built empirical library and its assembly log (skips preliminary analysis and assembly)
    if (params.diann_skip_preliminary_analysis) {
        if (!params.diann_empirical_assembly_lib || !params.diann_empirical_assembly_log) {
            error("--diann_skip_preliminary_analysis needs --diann_empirical_assembly_lib and --diann_empirical_assembly_log")
        }
        ch_empirical_lib = channel.of(file(params.diann_empirical_assembly_lib, checkIfExists: true))
        ch_empirical_log = channel.of(file(params.diann_empirical_assembly_log, checkIfExists: true))
    } else {
        ch_empirical_lib = channel.empty()
        ch_empirical_log = channel.empty()
    }

    DIA_PROTEOMICS_ANALYSIS(
        ch_dia_input,
        ch_searchdb,
        ch_expdesign,
        [params.diann_random_preanalysis, params.diann_empirical_assembly_ms_n, params.diann_random_preanalysis_seed],
        params.diann_scan_window,
        params.diann_mass_acc_automatic,
        params.diann_scan_window_automatic,
        params.diann_debug,
        params.diann_pg_level,
        ch_speclib,
        ch_empirical_lib,
        ch_empirical_log
    )

    // Post-processing
    ch_diann_report = DIA_PROTEOMICS_ANALYSIS.out.diann_report.map { _meta, report -> report }
    IQ(ch_diann_report, params.diann_maxlfq_q, params.diann_maxlfq_pgq, params.contaminant_pattern ?: '')

    if (!params.skip_post_msstats) {
        MSSTATS(
            DIA_PROTEOMICS_ANALYSIS.out.msstats_in.map { _meta, msstats -> msstats },
            params.contrasts,
            params.ref_condition ?: '',
            params.msstats_remove_one_feat_prot,
            params.msstatslfq_removeFewMeasurements,
            params.msstatslfq_feature_subset_protein,
            params.msstatslfq_quant_summary_method,
            params.msstats_threshold
        )
    }

    // Proteomics QC report (pmultiqc) over the mzTab, experimental design, mzML statistics and DIA-NN report
    if (params.enable_pmultiqc) {
        def pmultiqc_config = file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        PMULTIQC(
            DIA_PROTEOMICS_ANALYSIS.out.final_result.map { _meta, mztab -> mztab }
                .mix(
                    ch_expdesign.map { _meta, expdesign -> expdesign },
                    DIA_PROTEOMICS_ANALYSIS.out.mzml_statistics.map { _meta, stats -> stats },
                    ch_diann_report
                )
                .collect()
                .map { files -> [[id: 'experiment'], files, pmultiqc_config] }
        )
    }

    emit:
    multiqc_files = channel.empty() // every module reports its version on the 'versions' topic (main.nf collects them)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
