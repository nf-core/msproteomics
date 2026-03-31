/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DIA WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_msproteomics_pipeline'

// SDRF generation (bookkeeping only)
include { GENERATE_SDRF_FROM_SAMPLESHEET } from '../modules/local/generate_sdrf_from_samplesheet/main'

// DIA workflow modules
include { FILE_PREPARATION       } from '../subworkflows/local/file_preparation/main'
include { DIA_PROTEOMICS_ANALYSIS } from '../subworkflows/nf-core/dia_proteomics_analysis/main'
include { MSSTATS_LFQ as MSSTATS } from '../modules/local/msstats_lfq/main'
include { DIANNR                 } from '../modules/local/diannr/main'
include { IQ                     } from '../modules/local/iq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def parse_samplesheet() {
    return Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample    = row.sample
            def filepath  = row.spectra
            def condition = row.condition ?: sample
            def fraction  = row.fraction ?: null
            def label     = row.label ?: null

            def meta = [
                id:        fraction ? "${sample}_f${fraction}" : sample,
                group:     condition,
                sample:    sample,
                condition: condition,
                label:     label,
                fraction:  fraction ? fraction as int : 1
            ]
            [ meta, file(filepath) ]
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MSPROTEOMICS_DIANN {

    main:
    ch_input = parse_samplesheet()

    GENERATE_SDRF_FROM_SAMPLESHEET(
        Channel.fromPath(params.input, checkIfExists: true),
        WorkflowUtils.generateSdrfBookkeeping(params)
    )

    def organism = WorkflowUtils.convertOrganismToStandardName(params.organism ?: 'Homo sapiens')

    // Convert raw files if needed
    FILE_PREPARATION(ch_input)

    // Resolve database
    def database
    if (params.database) {
        database = params.database
    } else if (params.databases && params.databases[organism]?.database) {
        database = params.databases[organism].database
    } else {
        error "Neither --database is set nor a default database is found for ${organism}"
    }

    // Build DIA input channel with params from nextflow.config
    ch_dia_input = FILE_PREPARATION.out.results.map { meta, ms_file ->
        [
            meta,
            ms_file,
            params.enzyme ?: 'Trypsin/P',
            params.fixed_mods ?: '',
            params.variable_mods ?: 'Oxidation (M)',
            params.diann_library_ms1_acc,
            params.diann_library_mass_acc,
            'ppm',
            'ppm'
        ]
    }

    Channel.fromPath(database).map { fasta ->
        [[id: fasta.baseName], fasta]
    }.set { ch_searchdb }

    // Build experimental design in quantms two-section format (required by QUANTMSUTILS_DIANN2MZTAB).
    // Section 1: file table (Fraction_Group, Fraction, Spectra_Filepath, Label, Sample)
    // Section 2: sample table (Sample, MSstats_Condition, MSstats_BioReplicate)
    // Sections separated by an empty line.
    // Only generated when condition column provides distinct values.
    ch_expdesign_lines = ch_input
        .filter { meta, _file -> meta.condition != meta.sample }

    ch_expdesign_file_section = ch_expdesign_lines
        .map { meta, ms_file ->
            def mzml_name = ms_file.name.replaceAll(/\.(raw|RAW|d)$/, '.mzML')
            "1\t${meta.fraction}\t${mzml_name}\t1\t${meta.sample}"
        }
        .collectFile(name: 'expdesign_files.tsv', sort: true, newLine: true)

    ch_expdesign_sample_section = ch_expdesign_lines
        .map { meta, _ms_file ->
            "${meta.sample}\t${meta.condition}\t1"
        }
        .unique()
        .collectFile(name: 'expdesign_samples.tsv', sort: true, newLine: true)

    // Merge the two sections into a single experimental design file.
    // Use collectFile to avoid blocking file I/O inside a .map{} closure.
    ch_expdesign_file = ch_expdesign_file_section
        .combine(ch_expdesign_sample_section)
        .flatMap { file_section, sample_section ->
            // Read file contents outside of a closure that returns a channel value.
            // collectFile handles the I/O natively via Nextflow's file system.
            def file_lines = file_section.text.trim()
            def sample_lines = sample_section.text.trim()
            def content = "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n"
            content += file_lines + "\n\n"
            content += "Sample\tMSstats_Condition\tMSstats_BioReplicate\n"
            content += sample_lines + "\n"
            [content]
        }
        .collectFile(name: 'experimental_design.tsv', storeDir: "${workDir}/tmp")

    ch_expdesign_meta = ch_expdesign_file
        .map { expdesign -> [ [id: 'experiment'], expdesign ] }
        .ifEmpty( [ [id: 'experiment'], [] ] )

    // Spectral library selection
    if (params.diann_speclib != null && params.diann_speclib.toString() != "") {
        ch_speclib = Channel.of(file(params.diann_speclib, checkIfExists: true))
    } else {
        ch_speclib = Channel.empty()
    }

    // Handle diann_skip_preliminary_analysis mode
    if (params.diann_skip_preliminary_analysis) {
        ch_empirical_lib = Channel.of(file(params.diann_speclib))
        ch_empirical_log = Channel.fromPath(params.diann_empirical_assembly_log)
    } else {
        ch_empirical_lib = Channel.empty()
        ch_empirical_log = Channel.empty()
    }

    DIA_PROTEOMICS_ANALYSIS(
        ch_dia_input,
        ch_searchdb,
        ch_expdesign_meta,
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
    DIANNR(ch_diann_report, params.diann_maxlfq_q, params.diann_maxlfq_pgq, params.contaminant_pattern)
    IQ(ch_diann_report, params.diann_maxlfq_q, params.diann_maxlfq_pgq, params.contaminant_pattern)

    if (!params.skip_post_msstats) {
        MSSTATS(DIA_PROTEOMICS_ANALYSIS.out.msstats_in.map { _meta, msstats -> msstats })
    }

    emit:
    versions = Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
