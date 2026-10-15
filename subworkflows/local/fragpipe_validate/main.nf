/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_VALIDATE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PSM validation and site localization, as FragpipeRun.java:1603-1668 configures them.

    Modules:
    - PERCOLATOR (per searched pepXML basename, CmdPercolator)
    - PHILOSOPHERPEPTIDEPROPHET (per pepXML, or one --combine run per group, CmdPeptideProphet)
    - PTMPROPHET (per validated pepXML, CmdPtmProphet)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PERCOLATOR     } from '../../../modules/local/percolator/main'
include { PHILOSOPHERPEPTIDEPROPHET } from '../../../modules/local/philosopherpeptideprophet/main'
include { PTMPROPHET     } from '../../../modules/local/ptmprophet/main'
include { percolatorInitWeights; shouldRunTool; getToolArgs; getDecoyTag; requireToolField; fragpipeGroupDir; fragpipeLcmsOrderKey; compareFragpipeLcmsKeys; sortByRank; fragpipeLcmsName } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

workflow FRAGPIPE_VALIDATE {
    take:
    ch_pin                // channel: [ val(meta), path(pin) ] - MSFragger PIN (<stem>.pin)
    ch_pin_rescored       // channel: [ val(meta), path(pin) ] - MSBooster <stem>_edited.pin, or the MSFragger PIN
    ch_pepxml             // channel: [ val(meta), path(pepxml(s)) ] - MSFragger/Crystal-C pepXML, a list for ranked output
    ch_spectra            // channel: [ val(meta), path(spectra) ] - the searched LC-MS file
    ch_calibrated_mzml    // channel: [ val(meta), path(<stem>_calibrated.mzML) ] - when MSFragger writes it
    ch_uncalibrated_mzml  // channel: [ val(meta), path(<stem>_uncalibrated.mzML) ] - .raw/.d searches
    ch_fasta              // channel: path(fasta) - value channel for database
    ch_tool_configs       // channel: val(tool_configs_map) - parsed JSON config

    main:

    // Decoy tag (workflow database.decoy-tag) - single source for every decoy-aware tool
    ch_decoy_tag = ch_tool_configs.map { configs -> getDecoyTag(configs) }

    ch_pepxml_keyed = ch_pepxml.map { meta, pepxml -> [meta.id, meta, sortByRank(pepxml)] }

    //
    // Spectra each pepXML refers to: the searched file plus MSFragger's _calibrated/_uncalibrated mzML.
    // FragPipe names <stem>_calibrated.mzML when calibrated mzML is written, <stem>_uncalibrated.mzML for
    // .raw/.d, else the file (CmdPercolator.java:219-224, CmdPeptideProphet.java:439-450); the converters
    // fall back to the files that exist next to it, so all of them are staged.
    //
    ch_spectra_all = ch_spectra
        .map { meta, spectra -> [meta.id, meta, spectra] }
        .join(ch_calibrated_mzml.map { meta, mzml -> [meta.id, mzml] }, remainder: true)
        .join(ch_uncalibrated_mzml.map { meta, mzml -> [meta.id, mzml] }, remainder: true)
        .combine(ch_tool_configs)
        .map { key, _meta, spectra, calibrated, uncalibrated, configs ->
            if (spectra == null) {
                error("FRAGPIPE_VALIDATE: MSFragger mzML for ${key} has no searched spectra file")
            }
            def has_calibrated = requireToolField(configs, 'msfragger', 'has_calibrated_mzml')
            if (has_calibrated && calibrated == null) {
                error("FRAGPIPE_VALIDATE: ${key}: msfragger.write_calibrated_mzml is on with calibration, but MSFragger wrote no _calibrated.mzML")
            }
            [key, spectra, [calibrated, uncalibrated].findAll { f -> f != null }, has_calibrated]
        }

    // Every validated pepXML carries the name of its LC-MS file, FragPipe's order key (InputLcmsFile.compareTo)
    ch_pepxml_named = ch_pepxml_keyed
        .join(ch_spectra.map { meta, spectra -> [meta.id, spectra.name] }, failOnMismatch: true, failOnDuplicate: true)
        .map { key, meta, pepxml, spectra_name -> [key, meta + [spectra_name: spectra_name], pepxml] }

    //
    // STEP 1: Percolator + PercolatorOutputToPepXML (CmdPercolator.java:117-261)
    //
    // the run decision, then the joins: when Percolator runs, every searched file must have all its
    // inputs; when it does not, MSFragger may write no PIN at all (output_format tsv_pepXML)
    ch_percolator_on = ch_tool_configs.filter { configs -> shouldRunTool(configs, 'percolator') }
    ch_for_percolator = ch_pin
        .combine(ch_percolator_on)
        .map { meta, pin, _configs -> [meta.id, pin] }
        .join(ch_pin_rescored.combine(ch_percolator_on).map { meta, pin, _configs -> [meta.id, pin] }, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_pepxml_named.combine(ch_percolator_on).map { key, meta, pepxml, _configs -> [key, meta, pepxml] }, failOnMismatch: true, failOnDuplicate: true)
        .map { key, pin, rescored, meta, pepxml -> [key, meta, pin, rescored, pepxml] }
        .join(ch_spectra_all.combine(ch_percolator_on).map { key, spectra, msfragger_mzml, has_calibrated, _configs -> [key, spectra, msfragger_mzml, has_calibrated] }, failOnMismatch: true, failOnDuplicate: true)
        .combine(ch_tool_configs)
        .map { _key, meta, pin, rescored, pepxml, spectra, msfragger_mzml, has_calibrated, configs ->
            // MSBooster ran for this file when its rescored PIN is MSBooster's <stem>_edited.pin; the
            // TMT label check runs no MSBooster and passes the MSFragger PIN as both inputs
            def msbooster = rescored.name != pin.name
            def keep_tsv = requireToolField(configs, 'percolator', 'keep_tsv_files')
            // Metaproteomics reads the target TSV, so FragPipe keeps it (CmdPercolator.java:233-237)
            def keep_target = keep_tsv || shouldRunTool(configs, 'metaproteomics')
            [meta, pin, msbooster ? rescored : [], pepxml, spectra, msfragger_mzml, has_calibrated,
             getToolArgs(configs, 'percolator'), requireToolField(configs, 'percolator', 'min_prob'),
             keep_target, keep_tsv, shouldRunTool(configs, 'crystalc'),
             // Write sub mzML: Percolator saves its weights (CmdPercolator.java:180-183; tabRun.isWriteSubMzml(), FragpipeRun.java:1642)
             shouldRunTool(configs, 'writesubmzml'),
             // second pass: the weights found for a _sub.mzML input (subMzmlWeights; CmdPercolator.java:185-201)
             percolatorInitWeights(meta)]
        }

    PERCOLATOR(ch_for_percolator, ch_decoy_tag)

    //
    // STEP 2: PeptideProphet when Percolator does not run (FragpipeRun.java:1604-1623) + RewritePepxml
    //
    ch_pp_files = ch_pepxml_named
        .join(ch_spectra_all, failOnMismatch: true, failOnDuplicate: true)
        .combine(ch_tool_configs)
        // With both on, FragPipe runs PeptideProphet and Percolator then overwrites the same interact
        // files (FragpipeRun.java:1603-1652), so only Percolator's result matters
        .filter { row -> shouldRunTool(row[-1], 'peptideprophet') && !shouldRunTool(row[-1], 'percolator') }
        .map { _key, meta, pepxml, spectra, msfragger_mzml, has_calibrated, configs ->
            [meta, pepxml, [spectra] + msfragger_mzml, fragpipeLcmsName(spectra.name, has_calibrated),
             getToolArgs(configs, 'peptideprophet'), requireToolField(configs, 'peptideprophet', 'combine')]
        }
        .branch { row ->
            combined: row[5]
            per_file: true
        }

    // Per-file mode: one PeptideProphet per pepXML, ranked files included (CmdPeptideProphet.java:226-288)
    ch_pp_per_file = ch_pp_files.per_file
        .flatMap { meta, pepxmls, spectra, lcms_name, cli, combine ->
            def files = [pepxmls].flatten()
            files.collect { f -> [meta + [pp_files: files.size()], f, spectra, [lcms_name], cli, combine] }
        }

    // Combine mode: one --combine run per group over the group's pepXMLs in FragPipe order (CmdPeptideProphet.java:302-336)
    ch_pp_combined = ch_pp_files.combined
        .toList()
        .flatMap { rows ->
            rows.groupBy { row -> row[0].group }.collect { group, group_rows ->
                def ordered = group_rows.sort(false) { a, b -> compareFragpipeLcmsKeys(fragpipeLcmsOrderKey(a[0]), fragpipeLcmsOrderKey(b[0])) }
                def group_meta = ordered[0][0] + [id: fragpipeGroupDir(group), group: group, pp_combined_ids: ordered.collect { row -> row[0].id }]
                [group_meta, ordered.collectMany { row -> [row[1]].flatten() }, ordered.collectMany { row -> row[2] },
                 ordered.collect { row -> row[3] }.unique(), ordered[0][4], true]
            }
        }

    PHILOSOPHERPEPTIDEPROPHET(ch_pp_per_file.mix(ch_pp_combined), ch_fasta, ch_decoy_tag)

    // Per-file outputs back to one item per LC-MS file
    ch_pp_out = PHILOSOPHERPEPTIDEPROPHET.out.pepxml
        .branch { meta, _pepxml ->
            combined: meta.pp_combined_ids != null
            per_file: true
        }
    ch_pp_validated = ch_pp_out.per_file
        .map { meta, pepxml -> [groupKey(meta.id, meta.pp_files), meta.findAll { k, _v -> k != 'pp_files' }, pepxml] }
        .groupTuple()
        .map { _key, metas, pepxmls -> [metas[0], sortByRank(pepxmls)] }
        .mix(ch_pp_out.combined.map { meta, pepxml -> [meta.findAll { k, _v -> k != 'pp_combined_ids' } + [lcms_ids: meta.pp_combined_ids], pepxml] })

    ch_validated = PERCOLATOR.out.pepxml.map { meta, pepxml -> [meta, sortByRank(pepxml)] }
        .mix(ch_pp_validated)

    //
    // STEP 3: PTMProphet, one run per validated pepXML (CmdPtmProphet.java:90-124). It does not
    // replace the validated pepXMLs FragPipe passes on (FragpipeRun.java:1656-1667).
    //
    // PTMProphet reads the spectra its pepXML names (msms_run_summary base_name); they are staged with
    // the pepXML, whose base_name the module points at the staged copy (works without a shared filesystem).
    ch_spectra_by_id = ch_spectra_all
        .map { key, spectra, msfragger_mzml, has_calibrated -> [key, [spectra] + msfragger_mzml, fragpipeLcmsName(spectra.name, has_calibrated)] }
        .toList()
        .map { rows -> [rows.collectEntries { row -> [(row[0].toString()): [row[1], row[2]]] }] }

    ch_for_ptmprophet = ch_validated
        .combine(ch_tool_configs)
        .filter { _meta, _pepxml, configs -> shouldRunTool(configs, 'ptmprophet') }
        .combine(ch_spectra_by_id)
        .flatMap { meta, pepxml, configs, by_id ->
            def ids = (meta.lcms_ids ?: [meta.id]).collect { id -> id.toString() }
            def missing = ids.findAll { id -> !by_id.containsKey(id) }
            if (missing) {
                error("FRAGPIPE_VALIDATE: no spectra for ${missing} of ${meta.id}")
            }
            def files = [pepxml].flatten()
            files.collect { f ->
                [meta + [ptm_files: files.size()], f, ids.collectMany { id -> by_id[id][0] }, ids.collect { id -> by_id[id][1] }.unique(),
                 getToolArgs(configs, 'ptmprophet')]
            }
        }

    PTMPROPHET(ch_for_ptmprophet)

    ch_mod_pepxml = PTMPROPHET.out.mod_pepxml
        .map { meta, mod -> [groupKey(meta.id, meta.ptm_files), meta.findAll { k, _v -> k != 'ptm_files' }, mod] }
        .groupTuple()
        .map { _key, metas, mods -> [metas[0], sortByRank(mods)] }

    // Without Percolator and PeptideProphet the search pepXMLs go on unvalidated (FragpipeRun.java:1604-1652)
    ch_unvalidated = ch_pepxml_named
        .combine(ch_tool_configs)
        .filter { _key, _meta, _pepxml, configs -> !shouldRunTool(configs, 'percolator') && !shouldRunTool(configs, 'peptideprophet') }
        .map { _key, meta, pepxml, _configs -> [meta, pepxml] }

    ch_pepxml_validated = ch_validated.mix(ch_unvalidated)

    emit:
    pepxml       = ch_pepxml_validated                          // channel: [ val(meta), path(pepxml(s)) ] - sharedPepxmlFiles: ProteinProphet, iProphet
    pepxml_group = ch_pepxml_validated.mix(ch_mod_pepxml)       // channel: [ val(meta), path(pepxml(s)) ] - everything in a group dir for the filter
    target_psms  = PERCOLATOR.out.target_psms                   // channel: [ val(meta), path(target_psms.tsv) ]
    decoy_psms   = PERCOLATOR.out.decoy_psms                    // channel: [ val(meta), path(decoy_psms.tsv) ]
    mod_pepxml   = ch_mod_pepxml                                // channel: [ val(meta), path(mod.pep.xml(s)) ] - PTM localized
    results_dir  = PERCOLATOR.out.results_dir                   // channel: [ val(meta), path(results_dir) ]
    versions     = PERCOLATOR.out.versions_percolator
        .mix(PERCOLATOR.out.versions_fragpipe)
        .mix(PHILOSOPHERPEPTIDEPROPHET.out.versions_philosopher)
        .mix(PHILOSOPHERPEPTIDEPROPHET.out.versions_fragpipe)
        .mix(PTMPROPHET.out.versions_ptmprophet)
        .mix(PTMPROPHET.out.versions_fragpipe)
}
