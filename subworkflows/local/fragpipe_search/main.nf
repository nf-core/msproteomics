/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_SEARCH SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Database search with MSFragger, optional AI rescoring with MSBooster,
    and optional chimeric artifact removal with Crystal-C.

    Supports two modes:
    - Standard (num_slices <= 1): search per-sample on the full FASTA
    - Split database (num_slices > 1): calibrate each file (calibrate_mass 1/2),
      split FASTA into M chunks, search M chunks x N files in parallel with each
      file's params, merge results per file

    MSFragger builds its pepindex inside every search call, as in FragPipe (FragPipe never
    prebuilds one for a search; CmdMsfraggerDigest only runs for transfer learning,
    FragpipeRun.java:2255). The pepindex header stores the absolute database_name, so an
    index built in another task directory is never reused by MSFragger.

    Post-Search (PARALLEL branches):
    - CrystalC (optional) - Modifies pepXML (chimeric artifact removal)
    - MSBooster (optional) - Modifies PIN (AI-based feature rescoring), one run over all files

    Note: Crystal-C reads mzML, mzXML and RAW only; other formats stop the run, as in FragPipe.

    Modules:
    - MSFRAGGERCALIBRATE (per-file mass calibration + param optimization for split search)
    - SPLITFASTA (FASTA splitting for split search)
    - MSFRAGGER (database search)
    - MERGESPLITSEARCH (merge split search results)
    - MSBOOSTER (AI rescoring, optional)
    - CRYSTALC (chimeric artifact removal, optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MSFRAGGER            } from '../../../modules/local/msfragger/main'
include { MSFRAGGERCALIBRATE } from '../../../modules/local/msfraggercalibrate/main'
include { SPLITFASTA } from '../../../modules/local/splitfasta/main'
include { MERGESPLITSEARCH } from '../../../modules/local/mergesplitsearch/main'
include { MSBOOSTER           } from '../../../modules/local/msbooster/main'
include { CRYSTALC            } from '../../../modules/local/crystalc/main'
include { SORTFASTA } from '../../../modules/local/sortfasta/main'
include { MSFRAGGERAPPLYCALIBRATION } from '../../../modules/local/msfraggerapplycalibration/main'
include { WRITETEXTFILE as WRITE_FRAGGER_PARAMS } from '../../../modules/local/writetextfile/main'
include { msboosterThreads; setFragpipeParamValues } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { shouldRunTool; getToolEntry; getToolArgs; crystalcRunIndex; fragpipeLcmsOrderKey; compareFragpipeLcmsKeys; fraggerParamsType; fraggerCalibrateMass; fragpipeValidationLcms; fragpipeLcmsDir; fraggerWritesCalibratedMzml } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

workflow FRAGPIPE_SEARCH {
    take:
    ch_mzml              // channel: [ val(meta), path(mzml/d) ] - spectra files
    ch_fasta             // channel: path(fasta) - value channel for database
    ch_tool_configs      // channel: val(tool_configs_map) - parsed JSON config
    ch_fragger_params    // channel: val([[params_type, path(fragger*.params)], ...]) - value channel, one MSFragger
                         //   params file per FragPipe data type ('DDA', 'DIA', 'GPF-DIA', 'DDA+'; CmdMsfragger.java:444-457)
    ch_msbooster_fragger_params // channel: path(fragger*.params) - value channel, the one params file MSBooster reads (CmdMSBooster.java:134-149)
    ch_msbooster_params  // channel: path(msbooster.config) - value channel, workflow-derived MSBooster params
    ch_num_slices        // channel: val(int) - value channel, msfragger.misc.slice-db of the workflow
                         //   (1 = one search on the full FASTA, >1 = split database search)

    main:

    // Each MSFragger params set read once: [[params_type, path, text], ...] (value channel)
    ch_fragger_sets = ch_fragger_params.map { sets -> sets.collect { set -> [set[0], set[1], set[1].text] } }

    // Filter to samples that should run MSFragger; each file gets the params set of its data type
    // (CmdMsfragger.java:541-552: DDA, DIA + DIA-Lib, GPF-DIA, DDA+).
    ch_for_msfragger = ch_mzml
        .combine(ch_tool_configs)
        .filter { _meta, _mzml, configs -> shouldRunTool(configs, 'msfragger') }
        .map { meta, mzml, _configs -> [meta, mzml] }
        .combine(ch_fragger_sets.map { sets -> [sets] })
        .map { meta, mzml, sets ->
            def params_type = fraggerParamsType(meta.data_type)
            def set = sets.find { s -> s[0] == params_type }
            [meta, mzml, set[1], set[2]]
        }

    // Create a meta-wrapped fasta for modules that need tuple input
    ch_fasta_tuple = ch_fasta.map { fasta -> [[id: 'database'], fasta] }

    // Number of database slices: msfragger.misc.slice-db of the workflow (TabMsfragger.getNumDbSlices(),
    // passed to CmdMsfragger.configure at FragpipeRun.java:1470). 1 = one search on the full FASTA,
    // >1 = msfragger_pep_split.py (CmdMsfragger.java:337, :519-528).
    ch_slices = ch_num_slices

    ch_search_mode = ch_for_msfragger
        .combine(ch_slices)
        .branch { _meta, _spectra, _params_file, _params_text, slices ->
            split: slices > 1
            standard: true
        }

    // ---------------------------------------------------------------------
    // SPLIT DATABASE (slices > 1)
    // Replicates FragPipe 24.0 tools/msfragger_pep_split.py, which FragPipe runs once per
    // data-type params set (CmdMsfragger.java:555-611); here per file, as the per-file patched
    // FragPipe used for golden comparison does:
    //   1. calibrate_mass 1/2: calibrate the file against the full sorted FASTA (--split1) and
    //      apply the optimized values to the file's params (msfragger_pep_split.py:498-611);
    //      calibrate_mass 0: search the file and its params as they are (:630-631)
    //   2. split the FASTA into M chunks
    //   3. search every chunk with --partial (:146-153)
    //   4. merge the chunk results per file (:401-496)
    // ---------------------------------------------------------------------

    // FragPipe's split database checks (CmdMsfragger.java:337-357) run in PARSEFRAGPIPEWORKFLOW,
    // which knows the manifest data types and fails the run as FragPipe does. calibrate_mass of
    // each file's params set comes from tool_configs (no params file read on the head node).
    ch_split_input = ch_search_mode.split
        .combine(ch_tool_configs)
        .map { meta, spectra, params_file, params_text, _slices, configs ->
            def calibrate_mass = fraggerCalibrateMass(configs, fraggerParamsType(meta.data_type))
            if (!(calibrate_mass in [0, 1, 2])) {
                error("msfragger_pep_split.py needs calibrate_mass 0, 1 or 2 in ${params_file.name}, found ${calibrate_mass}")
            }
            [meta, spectra, params_file, params_text, calibrate_mass]
        }
        .branch { _meta, _spectra, _params_file, _params_text, calibrate_mass ->
            calibrate: calibrate_mass in [1, 2]
            direct: true
        }

    // The calibration search runs against the protein-sorted FASTA (msfragger_pep_split.py:614-620, 632-633)
    SORTFASTA(ch_fasta.combine(ch_split_input.calibrate.first()).map { items -> items[0] })

    MSFRAGGERCALIBRATE(
        ch_split_input.calibrate.map { meta, spectra, _params_file, params_text, _cal -> [meta, spectra, params_text] },
        SORTFASTA.out.fasta.first()
    )

    // Optimized values into the search params, calibrated spectra collected (msfragger_pep_split.py:524-611)
    MSFRAGGERAPPLYCALIBRATION(
        ch_split_input.calibrate
            .map { meta, spectra, params_file, _params_text, _cal -> [meta.id, meta, spectra, params_file] }
            .join(MSFRAGGERCALIBRATE.out.calibration.map { meta, log_file, mzbin -> [meta.id, log_file, mzbin] }, failOnMismatch: true, failOnDuplicate: true)
            .map { _id, meta, spectra, params_file, log_file, mzbin -> [meta, spectra, log_file, mzbin, params_file] }
    )

    // [meta, spectra to search, search params, its text]; the calibrated params file is a small task
    // output with no text emit, read once here
    ch_split_search_input = MSFRAGGERAPPLYCALIBRATION.out.calibrated
        .map { meta, spectra, params_file -> [meta, spectra, params_file, params_file.text] }
        .mix(ch_split_input.direct.map { meta, spectra, params_file, params_text, _cal -> [meta, spectra, params_file, params_text] })

    // The FASTA is split only when the workflow asks for slices (msfragger_pep_split.py:96-137).
    SPLITFASTA(
        ch_fasta_tuple.combine(ch_slices).filter { _meta, _fasta, slices -> slices > 1 }.map { meta, fasta, _slices -> [meta, fasta] },
        ch_slices
    )

    // SPLITFASTA outputs split_db/{0,1,...}/<fasta name>
    ch_chunk_fastas = SPLITFASTA.out.fasta_chunks
        .flatMap { _meta, fastas ->
            def fasta_list = fastas instanceof List ? fastas : [fastas]
            fasta_list.collect { fasta ->
                def chunk_idx = fasta.parent.name
                [[id: "chunk_${chunk_idx}", chunk_id: chunk_idx as int], fasta]
            }
        }

    // Every (file, chunk) pair with the file's params; MSFRAGGER adds --partial <chunk_id>
    // (msfragger_pep_split.py:146-153).
    ch_for_split_search = ch_split_search_input
        .combine(ch_chunk_fastas.map { chunk_meta, fasta -> [chunk_meta.chunk_id, fasta] })
        .map { meta, spectra, params_file, params_text, chunk_id, chunk_fasta ->
            def new_meta = meta.clone()
            new_meta.original_id = meta.id
            new_meta.chunk_id = chunk_id
            new_meta.id = "${meta.id}_chunk_${chunk_id}"
            [new_meta, spectra, params_file, params_text, chunk_fasta]
        }

    ch_for_standard_search = ch_search_mode.standard
        .combine(ch_fasta)
        .map { meta, spectra, params_file, params_text, _slices, fasta -> [meta, spectra, params_file, params_text, fasta] }

    MSFRAGGER(
        ch_for_standard_search.mix(ch_for_split_search)
            .map { meta, spectra, params_file, params_text, fasta -> [meta, spectra, params_file.name, params_text, fasta] }
    )

    // FragPipe writes each params set its files use into its work dir, with the run's FASTA as database_name
    // (CmdMsfragger.java:441-552), also for a split database: msfragger_pep_split.py calibrates and searches in a
    // temporary dir it deletes (:597-603, :645) and leaves that file as it is
    WRITE_FRAGGER_PARAMS(
        ch_for_msfragger
            .map { _meta, _mzml, params_file, params_text -> [params_file, params_text] }
            .unique { item -> item[0].name }
            .combine(ch_fasta)
            .map { params_file, params_text, fasta ->
                [[id: params_file.baseName], params_file.name, setFragpipeParamValues(params_text, [database_name: fasta.toUriString()])]
            }
    )

    ch_msfragger_out = [
        pepxml: MSFRAGGER.out.pepxml.filter { meta, _f -> meta.chunk_id == null },
        pin: MSFRAGGER.out.pin.filter { meta, _f -> meta.chunk_id == null },
        tsv: MSFRAGGER.out.tsv.filter { meta, _f -> meta.chunk_id == null },
        results_dir: MSFRAGGER.out.results_dir.filter { meta, _d -> meta.chunk_id == null },
        lcms_files: MSFRAGGER.out.lcms_files.branch { meta, _f -> split: meta.chunk_id != null; standard: true },
        log: MSFRAGGER.out.log.branch { meta, _f -> split: meta.chunk_id != null; standard: true },
    ]

    // Merge the chunks of each file (msfragger_pep_split.py:401-496); groupKey carries the slice count.
    ch_split_keyed = ch_split_search_input
        .map { meta, _spectra, _params_file, params_text -> [meta.id, meta, params_text] }

    // Each file's chunk result directories in chunk order (msfragger_pep_split.py merges chunks 0..N-1)
    ch_split_dirs = MSFRAGGER.out.results_dir
        .filter { meta, _dir -> meta.chunk_id != null }
        .combine(ch_slices)
        .map { meta, dir, slices -> [groupKey(meta.original_id, slices), [meta.chunk_id, dir]] }
        .groupTuple()
        .map { key, chunks -> [key.toString(), chunks.sort { c -> c[0] }.collect { c -> c[1] }] }

    ch_for_merge = ch_split_keyed
        .join(ch_split_dirs, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, params_text, dirs -> [meta, dirs, params_text] }

    MERGESPLITSEARCH(ch_for_merge, ch_fasta)

    ch_msfragger_pepxml = ch_msfragger_out.pepxml.mix(MERGESPLITSEARCH.out.pepxml)
    ch_msfragger_pin = ch_msfragger_out.pin.mix(MERGESPLITSEARCH.out.pin)
    ch_msfragger_tsv = ch_msfragger_out.tsv  // no per-file TSV from the merge
    ch_msfragger_results_dir = ch_msfragger_out.results_dir
    // Files MSFragger writes next to the LC-MS file (_uncalibrated.mzML for .raw/.d, _calibrated.mzML,
    // mzBIN caches) wherever it reads the original file: in the search, in split calibration, or
    // (calibrate_mass 0) in every chunk search, of which chunk 0 is taken.
    ch_split_meta = ch_split_keyed.map { id, meta, _params_text -> [id, meta] }
    ch_msfragger_lcms_files = ch_msfragger_out.lcms_files.standard
        .mix(MSFRAGGERCALIBRATE.out.uncalibrated_mzml)
        .mix(
            ch_msfragger_out.lcms_files.split
                .filter { meta, _files -> meta.chunk_id == 0 }
                .map { meta, files -> [meta.original_id, files] }
                // no failOnMismatch: lcms_files is optional (MSFragger writes none for an mzML input without calibration output)
                .join(ch_split_meta)
                .map { _id, files, meta -> [meta, files] }
        )
    // MSFragger console output per file for FRAGPIPESDRF (FragpipeRun.java:586): the search log;
    // in split mode the calibration log (calibrate_mass 1/2, it read the original file) or chunk 0's
    // search log (calibrate_mass 0).
    ch_msfragger_logs = ch_msfragger_out.log.standard
        .mix(MSFRAGGERCALIBRATE.out.log)
        .mix(
            ch_msfragger_out.log.split
                .filter { meta, _log_file -> meta.chunk_id == 0 }
                .map { meta, log_file -> [meta.original_id, log_file] }
                // no failOnMismatch: chunk 0 logs of calibrated files are dropped here (their log is MSFRAGGERCALIBRATE's)
                .join(ch_split_input.direct.map { meta, _spectra, _params_file, _params_text, _cal -> [meta.id, meta] })
                .map { _id, log_file, meta -> [meta, log_file] }
        )

    //
    // =====================================================================
    // POST-SEARCH: CrystalC + MSBooster (both paths converge here)
    // =====================================================================
    //

    //
    // Key per-sample search outputs for downstream joining
    //
    ch_pin_keyed = ch_msfragger_pin
        .map { meta, pin -> [meta.id, meta, pin] }

    // pepXML may be a list for ranked output (output_report_topN > 1)
    ch_pepxml_keyed = ch_msfragger_pepxml
        .map { meta, pepxml -> [meta.id, meta, pepxml] }

    ch_mzml_keyed = ch_mzml
        .map { meta, mzml -> [meta.id, meta, mzml] }

    //
    // CrystalC (optional, per-sample, PARALLEL with MSBooster) on the LC-MS file as given.
    // FragPipe refuses formats Crystal-C cannot read: mzML, mzXML, and RAW with the Thermo
    // libraries (CmdCrystalc.java:58-60, 138-145; CmdBase.getNotSupportedExts1).
    //
    ch_crystalc_params = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'crystalc') }
        .map { configs -> getToolArgs(configs, 'crystalc') }
        .collectFile(name: 'crystalc.params', newLine: true)

    ch_for_crystalc = ch_pepxml_keyed
        .join(ch_mzml_keyed.map { key, _meta, mzml -> [key, mzml] }, failOnMismatch: true, failOnDuplicate: true)
        .combine(ch_tool_configs)
        .filter { _key, _meta, _pepxml, _lcms, configs -> shouldRunTool(configs, 'crystalc') }
        .flatMap { _key, meta, pepxml, lcms, _configs ->
            def ext = lcms.name.tokenize('.').last()
            if (!(ext.toLowerCase() in ['mzml', 'mzxml', 'raw'])) {
                error("Crystal-C can't work with '.${ext.toLowerCase()}' files. You can convert files using msconvert from ProteoWizard.")
            }
            // One Crystal-C run per pepXML, ranked pepXMLs included (CmdCrystalc.java:179-219)
            def pepxmls = pepxml instanceof List ? pepxml : [pepxml]
            pepxmls.withIndex().collect { f, rank -> [meta + [crystalc_count: pepxmls.size(), spectra_name: lcms.name], f, lcms, rank] }
        }
        .toList()
        .flatMap { items -> crystalcRunIndex(items) }
        .combine(ch_crystalc_params)

    CRYSTALC(ch_for_crystalc, ch_fasta)

    //
    // MSBooster (optional, PARALLEL with CrystalC). The run has timsTOF data when any searched
    // file is .d or _diatracer.mzML (CmdMSBooster.java:151-157).
    //
    ch_has_im = ch_mzml
        .map { _meta, lcms -> lcms.name.toLowerCase() ==~ /.+(\.d|_diatracer\.mzml)/ }
        .collect()
        .map { flags -> flags.any { has_im -> has_im } }

    //
    // FragPipe's LC-MS directory after the search: each LC-MS file plus what MSFragger wrote next
    // to it (_uncalibrated.mzML for .raw/.d, _calibrated.mzML, mzBIN caches). Every later tool
    // finds its spectra there.
    //
    ch_lcms = fragpipeLcmsDir(ch_mzml, ch_msfragger_lcms_files)

    // The one spectra file Percolator / PeptideProphet name in their outputs (CmdPercolator.java:219-224,
    // CmdPeptideProphet.java:442-447): <base>_calibrated.mzML when MSFragger wrote it, else
    // <base>_uncalibrated.mzML for .raw/.d input, else the LC-MS file itself.
    // hasCalibratedMzml from the MSFragger settings, the same in every data-type params set (FragpipeRun.java:1642)
    ch_has_calibrated_mzml = ch_fragger_sets.map { sets -> fraggerWritesCalibratedMzml(sets[0][2]) }
    ch_validation_lcms = ch_lcms
        .combine(ch_has_calibrated_mzml)
        .map { meta, files, has_cal -> [meta, fragpipeValidationLcms(files, has_cal)] }

    ch_spectra_for_boost = ch_lcms
        .map { meta, files -> [meta.id, files] }

    // One MSBooster run over every file, as FragPipe (CmdMSBooster: one --paramsList call for all PINs).
    // Per-file runs are not equivalent: DIA-NN's RT predictions depend on the peptide batch, so
    // delta_RT_loess differs in the 4th decimal (measured on a FragPipe .raw run).
    // The run decision comes before the strict join: MSFragger writes no PIN when its output_format
    // has none (tsv_pepXML, PeptideProphet runs), and every searched file would be unmatched.
    ch_msbooster_on = ch_tool_configs.filter { configs -> shouldRunTool(configs, 'msbooster') }
    ch_msbooster_in = ch_pin_keyed
        .combine(ch_msbooster_on)
        .map { key, meta, pin, _configs -> [key, meta, pin] }
        .join(ch_spectra_for_boost.combine(ch_msbooster_on).map { key, files, _configs -> [key, files] }, failOnMismatch: true, failOnDuplicate: true)
        .map { _key, meta, pin, spectra -> [pin, spectra, meta.group, meta + [spectra_name: spectra[0].name]] }
        // FragPipe adds the PINs in the TreeMap order of the searched LC-MS files (FragpipeRun.java:1458,
        // InputLcmsFile.compareTo); MSBooster's HashSet keeps that order within a bucket (CmdMSBooster.java:200-212)
        .toSortedList { a, b -> compareFragpipeLcmsKeys(fragpipeLcmsOrderKey(a[3]), fragpipeLcmsOrderKey(b[3])) }
        .filter { items -> items }
        .multiMap { items ->
            msbooster: [[id: 'msbooster'], items.collect { item -> item[0] }, items.collectMany { item -> item[1] }]
            pin_dirs: items.collect { item -> item[2] }
        }

    // numThreads the parser wrote (workflow.threads / --fragpipe_threads); it sizes the task (conf/modules.config)
    // value channel: the MSBooster params are read once and the text reused below
    ch_msbooster_text = ch_msbooster_params.map { f -> f.text }
    ch_msbooster_threads = ch_msbooster_in.msbooster.combine(ch_msbooster_text).map { _meta, _pins, _spectra, text -> msboosterThreads(text) }

    // msbooster.spectral-library-path is staged; a declared file that does not exist stops the run
    ch_msbooster_library = ch_tool_configs.map { configs ->
        def library = getToolEntry(configs, 'msbooster').spectral_library
        library ? file(library, checkIfExists: true) : []
    }

    MSBOOSTER(ch_msbooster_in.msbooster, ch_msbooster_text, ch_msbooster_fragger_params, ch_msbooster_library, ch_has_im, ch_msbooster_threads, ch_msbooster_in.pin_dirs)

    //
    // Merge outputs
    //
    // pepXML: Use CrystalC output for samples that went through it,
    // MSFragger output for all others.
    //
    ch_pepxml_from_crystalc_keyed = CRYSTALC.out.pepxml_filtered
        .map { meta, pepxml -> [groupKey(meta.id, meta.crystalc_count), pepxml] }
        .groupTuple(sort: true)
        .map { key, pepxmls -> [key.getGroupTarget(), pepxmls.size() == 1 ? pepxmls[0] : pepxmls] }

    ch_pepxml_from_msfragger_keyed = ch_pepxml_keyed
        .map { key, meta, pepxml -> [key, meta, pepxml] }

    // Left join: all search pepXML samples, with CrystalC output where available.
    // remainder: true is required because CrystalC only runs on a subset of samples
    // (it is off, or the file is not searched). Unmatched left entries get null crystalc_pepxml,
    // and the ternary below falls back to msfragger_pepxml.
    ch_pepxml_final = ch_pepxml_from_msfragger_keyed
        .join(ch_pepxml_from_crystalc_keyed, remainder: true)
        .map { _key, meta, msfragger_pepxml, crystalc_pepxml ->
            def pepxml = crystalc_pepxml ?: msfragger_pepxml
            [meta, pepxml]
        }

    // PIN: Use MSBooster output if available, else search output.
    // MSBooster output is already per-sample.
    // Back to one PIN per file, keyed by the PIN name (the run returns every PIN)
    ch_pin_from_msbooster_keyed = MSBOOSTER.out.pin_edited
        .flatMap { _meta, pins -> [pins].flatten() }
        .map { pin -> [pin.name.replaceAll(/_edited\.pin$/, ".pin"), pin] }  // MSBooster keeps <stem>_edited.pin
        // no failOnMismatch: without MSBooster the left side is empty and every PIN is unmatched
        .join(ch_pin_keyed.map { key, _meta, pin -> [pin.name, key] }, failOnDuplicate: true)
        .map { _name, pin, key -> [key, pin] }

    // remainder: true is required because MSBooster only runs on a subset of samples.
    // Unmatched left entries get null msbooster_pin, ternary falls back to msfragger_pin.
    ch_pin_final = ch_pin_keyed
        .join(ch_pin_from_msbooster_keyed, remainder: true)
        .map { _key, meta, msfragger_pin, msbooster_pin ->
            def pin = msbooster_pin ?: msfragger_pin
            [meta, pin]
        }

    emit:
    pepxml            = ch_pepxml_final                    // channel: [ val(meta), path(pepxml(s)) ] - CrystalC or MSFragger; may be list for ranked
    pin               = ch_pin_keyed.map { _key, meta, pin -> [meta, pin] }  // channel: [ val(meta), path(pin) ] - raw search per-sample
    pin_rescored      = ch_pin_final                       // channel: [ val(meta), path(pin) ] - MSBooster edited or raw
    tsv               = ch_msfragger_tsv                   // channel: [ val(meta), path(tsv) ]
    results_dir       = ch_msfragger_results_dir           // channel: [ val(meta), path(results_dir) ] - per-sample
    msfragger_logs    = ch_msfragger_logs                  // channel: [ val(meta), path(log) ] - per-sample logs, or the split-mode calibration log
    msfragger_split_logs = ch_msfragger_out.log.split     // channel: [ val(meta), path(log) ] - split-mode chunk search logs (meta.id <file>_chunk_<i>)
    lcms              = ch_lcms                            // channel: [ val(meta), [ path(lcms), path(msfragger side files)... ] ] - FragPipe's LC-MS directory after the search
    // Write sub mzML (write_mzbin_all = 1): split database searches are refused with it (parse_fragpipe_workflow.py check_write_sub_mzml)
    mzbin_all         = MSFRAGGER.out.mzbin_all.filter { meta, _f -> meta.chunk_id == null }  // channel: [ val(meta), path(<stem>.mzBIN_all) ] - read by WRITESUBMZML
    validation_lcms   = ch_validation_lcms                 // channel: [ val(meta), path(spectra) ] - spectra file named by Percolator / PeptideProphet
    uncalibrated_mzml = ch_lcms.flatMap { meta, files -> files.findAll { f -> f.name.endsWith('_uncalibrated.mzML') }.collect { f -> [meta, f] } }  // channel: [ val(meta), path(*_uncalibrated.mzML) ] - for .d/.raw
    calibrated_mzml   = ch_lcms.flatMap { meta, files -> files.findAll { f -> f.name.endsWith('_calibrated.mzML') && !f.name.endsWith('_uncalibrated.mzML') }.collect { f -> [meta, f] } }  // channel: [ val(meta), path(*_calibrated.mzML) ] - write_calibrated_mzml with calibration
}
