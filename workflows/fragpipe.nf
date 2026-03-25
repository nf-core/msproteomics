/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE WORKFLOW (LFQ + TMT Quant)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Unified FragPipe workflow that handles both DDA LFQ and TMT quantification.
    When samples have labels (label column in samplesheet), TMT annotation is
    built automatically. When no labels are present, annotation is empty and
    FragPipe runs in LFQ mode.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_msproteomics_pipeline'

// SDRF generation (bookkeeping only)
include { GENERATE_SDRF_FROM_SAMPLESHEET } from '../modules/local/generate_sdrf_from_samplesheet/main'

// FragPipe subworkflows
include { FRAGPIPE_WF as FRAGPIPE_PIPELINE } from '../subworkflows/local/fragpipe/main'
include { FRAGPIPE_HEADLESS_WF             } from '../subworkflows/local/fragpipe_headless_wf/main'

// FragPipe utilities
include { generateFragpipeManifest } from '../subworkflows/local/fragpipe_utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MSPROTEOMICS_FRAGPIPE {

    take:
    ch_input    // channel: [ meta, spectra_file ]

    main:

    GENERATE_SDRF_FROM_SAMPLESHEET(
        Channel.fromPath(params.input, checkIfExists: true),
        WorkflowUtils.generateSdrfBookkeeping(params)
    )

    def organism = WorkflowUtils.convertOrganismToStandardName(params.organism ?: 'Homo sapiens')
    def db_channels = WorkflowUtils.resolveFragPipeDatabaseChannels(organism, params)

    ch_workflow_file = Channel.of(
        [[id: 'fragpipe_config'], file(params.fragpipe_workflow, checkIfExists: true)]
    )

    // TMT annotation: built from samplesheet label column.
    // Uses multi-plex TMTIntegrator format (plex\tchannel\tsample\tsample_name\tcondition\treplicate)
    // so the IonQuant module can extract per-plex 2-column annotations.
    // When no labels are present, ch_annotation is empty ([]) and FragPipe runs in LFQ mode.
    ch_annotation = ch_input
        .filter { meta, _file -> meta.label != null }
        .map { meta, _ms_file ->
            def plex = meta.group ?: meta.sample
            "plex\tchannel\tsample\tsample_name\tcondition\treplicate\n${plex}\t${meta.label.replaceAll(/^TMT/, '')}\t${meta.sample}\t${meta.sample}\t${meta.condition ?: meta.sample}\t1\n"
        }
        .collectFile(name: 'tmt_annotation.tsv', sort: true, keepHeader: true)
        .ifEmpty([])

    if (params.fragpipe_mode == 'allinone') {
        //
        // ALL-IN-ONE MODE: Run FragPipe headless in a single process
        //
        // Pre-converts .raw to .mzML before running FragPipe headless.
        // Generates a manifest from the samplesheet for FragPipe's --manifest flag.
        //

        // Collect samplesheet entries for manifest generation.
        // Each entry is [group, filename] where group = condition for experiment grouping.
        // After conversion, filenames will be .mzML — update extensions accordingly.
        ch_manifest = ch_input
            .map { meta, spectra_file ->
                def fname = spectra_file.name.toLowerCase().endsWith('.raw')
                    ? spectra_file.baseName + '.mzML'
                    : spectra_file.name
                [ meta.condition, fname ]
            }
            .collect(flat: false)
            .map { entries -> generateFragpipeManifest(entries, 'DDA') }

        // Resolve database path for headless mode
        ch_database_file = db_channels.ch_fasta
            .map { _meta, fasta -> fasta }

        // Workflow file (plain path, not tuple)
        ch_wf_file = Channel.of(file(params.fragpipe_workflow, checkIfExists: true))

        FRAGPIPE_HEADLESS_WF(
            ch_input,
            ch_database_file,
            ch_wf_file,
            ch_manifest
        )
    } else {
        //
        // PIPELINE MODE: Modular FragPipe subworkflows with Nextflow orchestration
        //
        FRAGPIPE_PIPELINE(
            ch_input,
            db_channels.ch_fasta,
            ch_workflow_file,
            db_channels.ch_prebuilt_db,
            ch_annotation,
            []
        )
    }

    emit:
    multiqc_report = Channel.empty()
    versions = Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
