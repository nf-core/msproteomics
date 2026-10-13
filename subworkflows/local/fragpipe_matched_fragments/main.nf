/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_MATCHED_FRAGMENTS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FragPipe "Export matched fragments" (tab-run.export_matched_fragments).

    FragPipe 24.0 orders CmdExportMatchedFragments right after PhilosopherReport
    (FragpipeRun.java:2411-2417 and :2490), before IonQuant, PTM-Shepherd, TMT-Integrator,
    SpecLibGen and every other tool that reads psm.tsv, and the export rewrites each group's
    psm.tsv in place. Groups whose meta.export_matched_fragments is set go through one
    aggregate EXPORTMATCHEDFRAGMENTS task; the others pass through unchanged.

    FRAGPIPE_MATCHED_FRAGMENTS_INPUTS builds the export's manifest and spectra once per run.

    Modules:
    - EXPORTMATCHEDFRAGMENTS (aggregate)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPORTMATCHEDFRAGMENTS } from '../../../modules/local/exportmatchedfragments/main'
include { fpopCoadaptrMeta; fragpipeRootGroupDir; matchedFragmentsManifest; shouldRunTool } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'
include { matchedFragmentsExportPsms } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

workflow FRAGPIPE_MATCHED_FRAGMENTS {
    take:
    ch_results_dir  // channel: [ val(meta), path(results_dir) ] - per-group Philosopher report output (meta.export_matched_fragments)
    ch_emf_inputs   // channel: [ path(manifest), [spectra], [umpire_spectra] ] - value channel (see EXPORTMATCHEDFRAGMENTS)
    ch_tool_configs // channel: val(tool_configs) - DIA-Umpire (umpire-se.params) and SpecLibGen (fragpipe.workflow) run flags
    aggregate_meta  // val(meta) for the aggregate task

    main:

    ch_groups = ch_results_dir
        .branch { meta, _dir ->
            export:      meta.export_matched_fragments
            passthrough: true
        }

    ch_for_export = ch_groups.export
        .toList()
        .filter { groups -> groups.size() > 0 }
        .map { groups -> [groups] }  // wrap so combine() does not flatten the group list
        .combine(ch_emf_inputs)
        .combine(ch_tool_configs)
        .multiMap { groups, manifest, spectra, umpire_spectra, configs ->
            def sorted = groups.toSorted { group -> group[0].id }
            // meta.fpop_coadaptr: FPOP coADAPTr edits the annotated psm.tsv afterwards (publish config)
            files:     [fpopCoadaptrMeta(aggregate_meta, configs, ['filter', 'report']), sorted.collect { group -> group[1] }, spectra, umpire_spectra]
            manifest:  manifest
            root:      fragpipeRootGroupDir(sorted.collect { group -> group[0] })
            flags:     [shouldRunTool(configs, 'diaumpire'), shouldRunTool(configs, 'speclibgen'), matchedFragmentsExportPsms(manifest.text, shouldRunTool(configs, 'speclibgen'))]
        }

    EXPORTMATCHEDFRAGMENTS(
        ch_for_export.files,
        ch_for_export.manifest,
        ch_for_export.root,
        ch_for_export.flags
    )

    // results/<group dir>/ back to the group's meta (the dir name is meta.id)
    ch_exported = EXPORTMATCHEDFRAGMENTS.out.results_dirs
        .flatMap { _meta, dirs -> [dirs].flatten().collect { dir -> [dir.name, dir] } }
        .join(ch_groups.export.map { meta, _dir -> [meta.id, meta] }, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, dir, meta -> [meta, dir] }

    emit:
    results_dir = ch_exported.mix(ch_groups.passthrough)  // channel: [ val(meta), path(results_dir) ] - psm.tsv with matched fragments when exported
    export_log  = EXPORTMATCHEDFRAGMENTS.out.log         // channel: [ val(meta), path(export_matched_fragments.log) ]
    versions    = channel.empty()  // channel: empty - every module reports its version on the 'versions' topic
}

workflow FRAGPIPE_MATCHED_FRAGMENTS_INPUTS {
    take:
    ch_file            // channel: [ val(meta), path(spectra) ] - every LC-MS file of the run (meta.experiment, meta.bioreplicate, meta.data_type)
    ch_spectra         // channel: [ val(meta), path(spectra) ] - spectra MSFragger searched and its _uncalibrated.mzML / _calibrated.mzML files
    ch_umpire_spectra  // channel: [ val(meta), path(mzml) ] - DIA-Umpire pseudo-spectra

    main:

    // Bruker .d directories and .raw files are never read: the export opens their mzML files
    ch_mzml = ch_spectra
        .map { _meta, spectra -> spectra }
        .flatten()
        .filter { spectra -> !(spectra.name.toLowerCase().endsWith('.d') || spectra.name.toLowerCase().endsWith('.raw')) }
        .collect(sort: true)
        .ifEmpty([])

    ch_manifest = ch_file
        .toList()
        .map { files -> [files] }  // wrap so combine() does not flatten the file list
        .combine(ch_mzml.map { files -> [files] })
        .map { files, spectra -> matchedFragmentsManifest(files, spectra.collect { f -> f.name }) }
        .collectFile(name: 'fragpipe-files.fp-manifest')

    ch_umpire = ch_umpire_spectra
        .map { _meta, spectra -> spectra }
        .collect(sort: true)
        .ifEmpty([])
        .map { files -> [files] }

    emit:
    // ch_manifest holds exactly one file (collectFile of one toList() item) and the others are
    // collect() results, so the combination has one item; first() only turns it into a value channel.
    inputs = ch_manifest.combine(ch_mzml.map { files -> [files] }).combine(ch_umpire).first()  // channel: [ path(manifest), [spectra], [umpire_spectra] ] - value channel
}
