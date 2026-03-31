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

workflow MSPROTEOMICS_FRAGPIPE {

    main:
    ch_input = parse_samplesheet()

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

    // Deduplicate input files: TMT samplesheets have multiple rows per file (one per channel).
    // Only convert each .raw file once. Applies to both headless and pipeline modes.
    // Override meta.id with spectra file basename so ThermoRawFileParser outputs unique filenames.
    // Without this, all plexes get meta.id="sample-01" (first channel) causing staging collisions.
    ch_unique_files = ch_input
        .unique { _meta, spectra_file -> spectra_file.name }
        .map { meta, spectra_file ->
            [meta + [id: spectra_file.baseName], spectra_file]
        }

    if (params.fragpipe_mode == 'headless') {
        //
        // ALL-IN-ONE MODE: Run FragPipe headless in a single process
        //
        // Pre-converts .raw to .mzML before running FragPipe headless.
        // Generates a manifest from the samplesheet for FragPipe's --manifest flag.
        // For TMT: generates per-experiment annotation content for auto-discovery.
        //

        // Collect samplesheet entries for manifest generation.
        // Each entry is [group, filename] where group = condition for experiment grouping.
        // After conversion, filenames will be .mzML — update extensions accordingly.
        // TMT samplesheets have multiple rows per file (one per channel).
        // Deduplicate by [condition, filename] so each file appears once in the manifest.
        ch_manifest = ch_input
            .map { meta, spectra_file ->
                def fname = spectra_file.name.toLowerCase().endsWith('.raw')
                    ? spectra_file.baseName + '.mzML'
                    : spectra_file.name
                [ meta.condition, fname ]
            }
            .unique()
            .collect(flat: false)
            .map { entries -> generateFragpipeManifest(entries, 'DDA') }

        // TMT annotation for headless mode: 3-column format (experiment\tchannel\tsample_name).
        // FragPipe headless auto-discovers *annotation.txt files in each experiment directory.
        // When no labels are present, annotation_content is empty and headless runs in LFQ mode.
        // FragPipe 24.0 requires globally unique sample names across all plexes
        // (CmdTmtIntegrator.java:205). Prefix with short plex ID from experiment name.
        ch_annotation_content = ch_input
            .filter { meta, _file -> meta.label != null }
            .map { meta, _ms_file ->
                def experiment = meta.condition ?: meta.sample
                def channel = meta.label.replaceAll(/^TMT/, '')
                def plexId = experiment.split('_')[-1]
                "${experiment}\t${channel}\t${plexId}_${meta.sample}"
            }
            .collect()
            .map { lines -> lines.join('\n') }
            .ifEmpty('')

        // File-to-experiment mapping using spectra file basename (matches ThermoRawFileParser output).
        // Dedup by spectra filename keeps 1 entry per unique .raw file.
        // Used by FRAGPIPE_HEADLESS to organize multi-plex files into per-experiment directories.
        ch_file_experiment_map = ch_input
            .unique { _meta, spectra_file -> spectra_file.name }
            .map { meta, spectra_file ->
                def experiment = meta.condition ?: meta.sample
                "${spectra_file.baseName}.mzML\t${experiment}"
            }
            .collect()
            .map { lines -> lines.join('\n') }
            .ifEmpty('')

        // Resolve database path for headless mode.
        // Handles both raw FASTA (ch_fasta) and prebuilt philosopher database (ch_prebuilt_db).
        // FragPipe headless handles decoy generation internally via Philosopher.
        ch_database_file = db_channels.ch_fasta
            .mix(db_channels.ch_prebuilt_db)
            .map { _meta, db -> db }

        // Workflow file (plain path, not tuple)
        ch_wf_file = Channel.of(file(params.fragpipe_workflow, checkIfExists: true))

        FRAGPIPE_HEADLESS_WF(
            ch_unique_files,
            ch_database_file,
            ch_wf_file,
            ch_manifest,
            ch_annotation_content,
            ch_file_experiment_map
        )
    } else {
        //
        // PIPELINE MODE: Modular FragPipe subworkflows with Nextflow orchestration
        //
        FRAGPIPE_PIPELINE(
            ch_unique_files,
            db_channels.ch_fasta,
            ch_workflow_file,
            db_channels.ch_prebuilt_db,
            ch_annotation,
            []
        )
    }

    emit:
    versions = Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
