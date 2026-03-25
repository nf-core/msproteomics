/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_SEARCH SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Database search with MSFragger, optional AI rescoring with MSBooster,
    and optional chimeric artifact removal with Crystal-C.

    Supports two modes:
    - Standard (num_slices <= 1): Prebuild pepindex once, search per-sample
    - Split database (num_slices > 1): Calibrate, split FASTA into M chunks,
      prebuild pepindex per chunk, search M chunks x N samples in parallel,
      merge results per sample

    Post-Search (PARALLEL branches):
    - CrystalC (optional) - Modifies pepXML (chimeric artifact removal)
    - MSBooster (optional) - Modifies PIN (AI-based feature rescoring)

    Note: CrystalC does NOT support Bruker .d files (only mzML, mzXML, RAW).

    Modules:
    - MSFRAGGER_CALIBRATE (mass calibration + param optimization for split search)
    - MSFRAGGER_INDEX (pepindex prebuilding)
    - SPLIT_FASTA (FASTA splitting for split search)
    - MSFRAGGER (database search)
    - MERGE_SPLIT_SEARCH (merge split search results)
    - MSBOOSTER (AI rescoring, optional)
    - CRYSTALC (chimeric artifact removal, optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MSFRAGGER            } from '../../../modules/local/msfragger/main'
include { MSFRAGGER_CALIBRATE } from '../../../modules/local/msfragger_calibrate/main'
include { MSFRAGGER_INDEX     } from '../../../modules/local/msfragger_index/main'
include { SPLIT_FASTA         } from '../../../modules/local/split_fasta/main'
include { MERGE_SPLIT_SEARCH  } from '../../../modules/local/merge_split_search/main'
include { MSBOOSTER           } from '../../../modules/local/msbooster/main'
include { CRYSTALC            } from '../../../modules/local/crystalc/main'
include { shouldRunTool; getToolArgs } from '../fragpipe_utils'

workflow FRAGPIPE_SEARCH {
    take:
    ch_mzml              // channel: [ val(meta), path(mzml/d) ] - spectra files
    ch_fasta             // channel: path(fasta) - value channel for database
    ch_tool_configs      // channel: val(tool_configs_map) - parsed JSON config
    ch_fragger_params    // channel: path(fragger.params) - value channel, [] for CLI mode
    ch_msbooster_params  // channel: path(msbooster.config) - value channel, [] for CLI mode
    ch_num_slices        // val: integer number of database slices (1 = standard, >1 = split)
                         // IMPORTANT: must be a plain integer, not a channel (used in if-condition and groupTuple size)

    main:

    // Filter to samples that should run MSFragger
    ch_for_msfragger = ch_mzml
        .combine(ch_tool_configs)
        .filter { _meta, _mzml, configs -> shouldRunTool(configs, 'msfragger') }
        .map { meta, mzml, _configs -> [meta, mzml] }

    // Create a meta-wrapped fasta for modules that need tuple input
    ch_fasta_tuple = ch_fasta.map { fasta -> [[id: 'database'], fasta] }

    if (ch_num_slices > 1) {
        //
        // =====================================================================
        // SPLIT DATABASE PATH (num_slices > 1)
        // Replicates FragPipe's msfragger_pep_split.py calibration flow:
        //   1. Calibrate ALL samples against full FASTA (--split1)
        //   2. Split FASTA into M chunks
        //   3. Prebuild pepindex per chunk (using calibrated params)
        //   4. Search M chunks x N samples in parallel (calibrated spectra)
        //   5. Merge results per sample across chunks
        // Per-file scoring is independent (confirmed from FragPipe source),
        // so M×N parallelism produces identical results to FragPipe's M tasks.
        // =====================================================================
        //

        //
        // STEP 0: Calibrate all samples against full FASTA (1 AGGREGATE task)
        // Runs MSFragger with --split1 to perform mass calibration and parameter
        // optimization. Produces .mzBIN_calibrated spectra and optimized params.
        // Matches FragPipe's calibrate() from msfragger_pep_split.py.
        //
        ch_all_spectra = ch_for_msfragger
            .map { _meta, spectra -> spectra }
            .collect()

        MSFRAGGER_CALIBRATE(ch_fasta, ch_all_spectra, ch_fragger_params, [])

        // Calibrated params (value channel): calibrate_mass=0, check_spectral_files=0,
        // optimized tolerances from --split1 output
        ch_calibrated_params = MSFRAGGER_CALIBRATE.out.params.first()

        // Calibrated spectra: .mzBIN_calibrated files (or original copies)
        // Re-associate with sample metadata by matching on stem name.
        ch_calibrated_spectra = MSFRAGGER_CALIBRATE.out.calibrated_spectra
            .flatten()
            .map { spectra ->
                // Strip .mzBIN_calibrated suffix to get original sample stem
                def stem = spectra.baseName.replaceAll(/\.mzBIN_calibrated$/, '')
                    .replaceAll(/\.mzBIN$/, '')
                [stem, spectra]
            }

        ch_for_msfragger_calibrated = ch_for_msfragger
            .map { meta, spectra ->
                [spectra.baseName, meta]
            }
            .join(ch_calibrated_spectra)
            .map { _stem, meta, calibrated_spectra ->
                [meta, calibrated_spectra]
            }


        //
        // STEP 1: Split FASTA into M chunks (1 task)
        //
        SPLIT_FASTA(ch_fasta_tuple, ch_num_slices)


        // Flatten chunk FASTAs into individual (chunk_meta, fasta) tuples.
        // SPLIT_FASTA outputs split_db/{0,1,...}/database.fasta
        ch_chunk_fastas = SPLIT_FASTA.out.fasta_chunks
            .flatMap { _meta, fastas ->
                def fasta_list = fastas instanceof List ? fastas : [fastas]
                fasta_list.collect { fasta ->
                    // Extract chunk index from parent directory name
                    def chunk_idx = fasta.parent.name
                    [[id: "chunk_${chunk_idx}", chunk_id: chunk_idx as int], fasta]
                }
            }

        //
        // STEP 2: Prebuild pepindex per chunk (M parallel tasks)
        // Uses calibrated params (calibrate_mass=0, optimized tolerances)
        //
        MSFRAGGER_INDEX(ch_chunk_fastas, ch_calibrated_params, [])


        //
        // STEP 3: Cartesian product: M chunks x N samples (M*N parallel search tasks)
        // Per-file scoring is independent in MSFragger (confirmed from FragPipe source:
        // CmdMsfragger.java splits files across separate invocations when cmdline > 30K).
        // M×N parallelism produces identical results to FragPipe's sequential per-chunk runs.
        //
        ch_indexed_chunks = MSFRAGGER_INDEX.out.indexed_fasta
            .map { chunk_meta, fasta, pepindex ->
                [chunk_meta.chunk_id, fasta, pepindex]
            }

        ch_for_split_search = ch_for_msfragger_calibrated
            .combine(ch_indexed_chunks)
            .map { meta, spectra, chunk_id, chunk_fasta, chunk_pepindex ->
                def new_meta = meta.clone()
                new_meta.original_id = meta.id
                new_meta.chunk_id = chunk_id
                new_meta.id = "${meta.id}_chunk_${chunk_id}"
                [new_meta, spectra, chunk_fasta, chunk_pepindex]
            }

        // Run MSFragger for each (sample, chunk) pair.
        // MSFRAGGER reads meta.chunk_id to add --partial flag automatically.
        // Calibrated params have calibrate_mass=0 (already calibrated by --split1).
        MSFRAGGER(
            ch_for_split_search.map { meta, spectra, _fasta, _pepindex -> [meta, spectra] },
            ch_for_split_search.map { _meta, _spectra, fasta, _pepindex -> fasta },
            ch_calibrated_params,
            ch_for_split_search.map { _meta, _spectra, _fasta, pepindex -> pepindex },
            []
        )


        //
        // STEP 4: Group by sample and merge across chunks (N tasks)
        // Collect pepXML, PIN, and histogram files from all chunks for each sample,
        // then merge via MERGE_SPLIT_SEARCH.
        //
        ch_split_pepxml = MSFRAGGER.out.pepxml
            .map { meta, pepxml -> [meta.original_id, pepxml] }
            .groupTuple(size: ch_num_slices, sort: 'hash')

        ch_split_pin = MSFRAGGER.out.pin
            .map { meta, pin -> [meta.original_id, pin] }
            .groupTuple(size: ch_num_slices, sort: 'hash')

        // Histograms are in the tsv output (_scores_histogram.tsv files)
        ch_split_histograms = MSFRAGGER.out.tsv
            .filter { meta, tsv ->
                def files = tsv instanceof List ? tsv : [tsv]
                files.any { it.name.contains('_scores_histogram') }
            }
            .map { meta, tsv ->
                def files = tsv instanceof List ? tsv : [tsv]
                def histograms = files.findAll { it.name.contains('_scores_histogram') }
                [meta.original_id, histograms]
            }
            .groupTuple(size: ch_num_slices, sort: 'hash')
            .map { id, histogram_lists -> [id, histogram_lists.flatten()] }

        // Join all grouped outputs by original sample ID
        ch_for_merge = ch_split_pepxml
            .join(ch_split_pin)
            .join(ch_split_histograms)
            .map { id, pepxmls, pins, histograms ->
                def meta = [id: id]
                // Flatten lists of lists from groupTuple
                def flat_pepxmls = pepxmls instanceof List ? pepxmls.flatten() : [pepxmls]
                def flat_pins = pins instanceof List ? pins.flatten() : [pins]
                def flat_histograms = histograms instanceof List ? histograms.flatten() : [histograms]
                [meta, flat_pepxmls, flat_pins, flat_histograms]
            }

        MERGE_SPLIT_SEARCH(ch_for_merge, ch_num_slices, ch_calibrated_params, [])


        // Set merged outputs as the search results for downstream steps
        ch_msfragger_pepxml = MERGE_SPLIT_SEARCH.out.pepxml
        ch_msfragger_pin = MERGE_SPLIT_SEARCH.out.pin
        ch_msfragger_tsv = Channel.empty()  // No per-sample TSV from merge
        ch_msfragger_results_dir = Channel.empty()  // No per-sample results_dir from merge
        ch_msfragger_uncalibrated_mzml = Channel.empty()  // Not available in split mode

        // Versions from split-path processes
        ch_search_versions = MSFRAGGER.out.versions_msfragger
            .mix(MSFRAGGER.out.versions_fragpipe)
            .mix(MSFRAGGER_CALIBRATE.out.versions_msfragger)
            .mix(MSFRAGGER_CALIBRATE.out.versions_fragpipe)
            .mix(MSFRAGGER_INDEX.out.versions_msfragger)
            .mix(MSFRAGGER_INDEX.out.versions_fragpipe)
            .mix(MERGE_SPLIT_SEARCH.out.versions_msfragger)
            .mix(MERGE_SPLIT_SEARCH.out.versions_fragpipe)

    } else {
        //
        // =====================================================================
        // STANDARD PATH (num_slices <= 1)
        // Prebuild pepindex once on full FASTA, share with all sample searches.
        // Identical to previous behavior except MSFragger skips digest step.
        // =====================================================================
        //

        //
        // STEP 1: Prebuild pepindex on full FASTA (1 task)
        // MSFragger auto-detects pepindex co-located with FASTA and reuses it.
        //
        MSFRAGGER_INDEX(ch_fasta_tuple, ch_fragger_params, [])


        // Extract pepindex as a value channel (broadcast to all samples)
        ch_pepindex = MSFRAGGER_INDEX.out.indexed_fasta
            .map { _meta, _fasta, pepindex -> pepindex }
            .first()

        //
        // STEP 2: MSFragger search (per-sample, parallel with prebuilt pepindex)
        //
        MSFRAGGER(ch_for_msfragger, ch_fasta, ch_fragger_params, ch_pepindex, [])


        ch_msfragger_pepxml = MSFRAGGER.out.pepxml
        ch_msfragger_pin = MSFRAGGER.out.pin
        ch_msfragger_tsv = MSFRAGGER.out.tsv
        ch_msfragger_results_dir = MSFRAGGER.out.results_dir
        ch_msfragger_uncalibrated_mzml = MSFRAGGER.out.uncalibrated_mzml

        // Versions from standard-path processes
        ch_search_versions = MSFRAGGER.out.versions_msfragger
            .mix(MSFRAGGER.out.versions_fragpipe)
            .mix(MSFRAGGER_INDEX.out.versions_msfragger)
            .mix(MSFRAGGER_INDEX.out.versions_fragpipe)
    }

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
    // CrystalC (optional, per-sample, PARALLEL with MSBooster)
    // CrystalC does NOT support .d files - only mzML, mzXML, RAW.
    // Creates per-file config with runtime paths.
    //
    ch_crystalc_pre = ch_pepxml_keyed
        .join(ch_mzml_keyed.map { key, meta, mzml -> [key, mzml] }, failOnMismatch: false, failOnDuplicate: false)
        .combine(ch_tool_configs)
        .filter { key, meta, pepxml, mzml, configs ->
            shouldRunTool(configs, 'crystalc') &&
            !mzml.name.toLowerCase().endsWith('.d')  // CrystalC does not support .d files
        }
        .map { key, meta, pepxml, mzml, configs ->
            def args = getToolArgs(configs, 'crystalc')
            // Determine file extension for Crystal-C config
            def ext = mzml.name.toLowerCase().endsWith('.mzml') ? 'mzML' :
                      mzml.name.toLowerCase().endsWith('.raw') ? 'raw' : 'mzXML'
            // Crystal-C params file needs runtime paths
            def config_content = args ? args + '\n' : ''
            config_content += "raw_file_location = .\n"
            config_content += "raw_file_extension = ${ext}\n"
            config_content += "output_location = .\n"
            [meta, pepxml, mzml, config_content]
        }

    // Create config files via collectFile (avoids writing files in workflow scope).
    ch_crystalc_configs = ch_crystalc_pre
        .map { meta, _pepxml, _mzml, config_content ->
            ["crystalc_${meta.id}.params", config_content]
        }
        .collectFile()
        .map { config_file ->
            def id = config_file.name.replaceAll(/^crystalc_/, '').replaceAll(/\.params$/, '')
            [id, config_file]
        }

    ch_for_crystalc = ch_crystalc_pre
        .map { meta, pepxml, mzml, _config_content -> [meta.id, meta, pepxml, mzml] }
        .join(ch_crystalc_configs)
        .map { key, meta, pepxml, mzml, config_file -> [meta, pepxml, mzml, config_file] }

    CRYSTALC(ch_for_crystalc, ch_fasta)

    //
    // MSBooster (optional, per-sample, PARALLEL with CrystalC)
    // MSBooster requires mzML/MGF files and cannot read .d directories directly.
    // For Bruker .d input, MSFragger produces _uncalibrated.mzML during processing,
    // which MSBooster can use (handles _uncalibrated suffix for matching).
    //

    // Detect ion mobility data from original input files (Bruker timsTOF = .d directories)
    ch_has_im = ch_mzml
        .map { _meta, mzml -> mzml.name.toLowerCase().endsWith('.d') }
        .collect()
        .map { flags -> flags.any { it } }

    // For MSBooster: prefer uncalibrated mzML if available (Bruker .d), else use original spectra.
    // In split mode, uncalibrated_mzml is empty — use original spectra.
    ch_uncalibrated_mzml_keyed = ch_msfragger_uncalibrated_mzml
        .map { meta, mzml -> [meta.id, mzml] }

    ch_original_spectra_keyed = ch_mzml
        .map { meta, spectra -> [meta.id, spectra] }

    ch_spectra_for_boost = ch_uncalibrated_mzml_keyed
        .join(ch_original_spectra_keyed, remainder: true)
        .map { key, mzml, spectra ->
            [key, mzml ?: spectra]
        }

    // Per-sample MSBooster: join PIN with spectra, one invocation per sample
    ch_for_msbooster = ch_pin_keyed
        .join(ch_spectra_for_boost)
        .combine(ch_tool_configs)
        .filter { _key, _meta, _pin, _spectra, configs -> shouldRunTool(configs, 'msbooster') }
        .map { _key, meta, pin, spectra, _configs -> [meta, pin, spectra] }

    MSBOOSTER(ch_for_msbooster, ch_msbooster_params, ch_fragger_params, ch_has_im)

    //
    // Merge outputs
    //
    // pepXML: Use CrystalC output for samples that went through it,
    // MSFragger output for all others.
    //
    ch_pepxml_from_crystalc_keyed = CRYSTALC.out.pepxml_filtered
        .map { meta, pepxml -> [meta.id, pepxml] }

    ch_pepxml_from_msfragger_keyed = ch_pepxml_keyed
        .map { key, meta, pepxml -> [key, meta, pepxml] }

    // Left join: all search pepXML samples, with CrystalC output where available.
    // remainder: true is required because CrystalC only runs on a subset of samples
    // (e.g., .d files are excluded). Unmatched left entries get null crystalc_pepxml,
    // and the ternary below falls back to msfragger_pepxml.
    ch_pepxml_final = ch_pepxml_from_msfragger_keyed
        .join(ch_pepxml_from_crystalc_keyed, remainder: true)
        .map { key, meta, msfragger_pepxml, crystalc_pepxml ->
            def pepxml = crystalc_pepxml ?: msfragger_pepxml
            [meta, pepxml]
        }

    // PIN: Use MSBooster output if available, else search output.
    // MSBooster output is already per-sample.
    ch_pin_from_msbooster_keyed = MSBOOSTER.out.pin_edited
        .map { meta, pin -> [meta.id, pin] }

    // remainder: true is required because MSBooster only runs on a subset of samples.
    // Unmatched left entries get null msbooster_pin, ternary falls back to msfragger_pin.
    ch_pin_final = ch_pin_keyed
        .join(ch_pin_from_msbooster_keyed, remainder: true)
        .map { key, meta, msfragger_pin, msbooster_pin ->
            def pin = msbooster_pin ?: msfragger_pin
            [meta, pin]
        }

    emit:
    pepxml            = ch_pepxml_final                    // channel: [ val(meta), path(pepxml(s)) ] - CrystalC or MSFragger; may be list for ranked
    pin               = ch_pin_keyed.map { _key, meta, pin -> [meta, pin] }  // channel: [ val(meta), path(pin) ] - raw search per-sample
    pin_rescored      = ch_pin_final                       // channel: [ val(meta), path(pin) ] - MSBooster edited or raw
    tsv               = ch_msfragger_tsv                   // channel: [ val(meta), path(tsv) ]
    results_dir       = ch_msfragger_results_dir           // channel: [ val(meta), path(results_dir) ] - per-sample
    uncalibrated_mzml = ch_msfragger_uncalibrated_mzml     // channel: [ val(meta), path(*_uncalibrated.mzML) ] - for .d/.raw
    versions          = ch_search_versions
        .mix(MSBOOSTER.out.versions_msbooster)
        .mix(MSBOOSTER.out.versions_fragpipe)
        .mix(CRYSTALC.out.versions_crystalc)
        .mix(CRYSTALC.out.versions_fragpipe)
}
