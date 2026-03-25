/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_QUANT SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Quantification module orchestration.

    Modules:
    - IONQUANT (label-free quantification with MaxLFQ, aggregate)
    - TMTINTEGRATOR (TMT quantification, aggregate)
    - FREEQUANT (simple label-free using Philosopher, per-sample)
    - PHILOSOPHER_LABELQUANT (TMT reporter ion extraction via Philosopher, per-sample)
    - PHILOSOPHER_REPORT (generate TSVs after FreeQuant/Labelquant, per-sample)

    Execution: Aggregate (mostly)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IONQUANT                      } from '../../../modules/local/ionquant/main'
include { IONQUANT as IONQUANT_MS1      } from '../../../modules/local/ionquant/main'
include { IONQUANT as IONQUANT_ISOBARIC } from '../../../modules/local/ionquant/main'
include { TMTINTEGRATOR            } from '../../../modules/local/tmtintegrator/main'
include { FREEQUANT                } from '../../../modules/local/freequant/main'
include { PHILOSOPHER_LABELQUANT   } from '../../../modules/local/philosopher/labelquant/main'
include { PHILOSOPHER_REPORT       } from '../../../modules/local/philosopher/report/main'
include { shouldRunTool; getToolArgs; getToolModmasses; getToolReportArgs } from '../fragpipe_utils'

workflow FRAGPIPE_QUANT {
    take:
    ch_results_dir          // channel: [ val(meta), path(results_dir) ] - per-sample PHILOSOPHER_FILTER output
    ch_mzml                 // channel: [ val(meta), path(mzml) ] - mzML files
    ch_annotation           // channel: path(annotation_file) - optional annotation for TMT/IsoQuant
    ch_fasta                // channel: path(fasta) - protein database for PHILOSOPHER_REPORT
    ch_tool_configs         // channel: val(tool_configs_map) - parsed JSON config
    ch_tmtintegrator_config // channel: path(tmt-integrator.yml) - value channel, [] when not used
    aggregate_meta          // val(meta) for aggregate outputs

    main:

    //
    // STEP 1: IonQuant (label-free quantification, aggregate)
    // Combines all samples for MaxLFQ and MBR
    //
    ch_all_results = ch_results_dir
        .map { meta, dir -> dir }
        .collect(sort: true)
        .map { items -> [items] }  // Wrap in list to prevent combine from flattening

    ch_all_mzml = ch_mzml
        .map { meta, mzml -> mzml }
        .collect(sort: true)
        .map { items -> [items] }  // Wrap in list to prevent combine from flattening

    ch_for_ionquant = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'ionquant') }
        .combine(ch_all_results)
        .combine(ch_all_mzml)
        .map { configs, results_wrapped, mzml_wrapped ->
            def args = getToolArgs(configs, 'ionquant')
            def modmasses = getToolModmasses(configs, 'ionquant')
            [aggregate_meta, results_wrapped, mzml_wrapped, args, modmasses]
        }

    // Pass annotation through - [] signals "not provided" to modules
    ch_annotation_file = ch_annotation

    IONQUANT(
        ch_for_ionquant.map { meta, results, mzml, _args, _modmasses -> [meta, results, mzml] },
        ch_annotation_file,
        ch_for_ionquant.map { _meta, _results, _mzml, args, _modmasses -> args },
        ch_for_ionquant.map { _meta, _results, _mzml, _args, modmasses -> modmasses },
        []
    )

    //
    // STEP 2: TMT Two-Pass IonQuant (when TMTIntegrator is enabled)
    // Pass 1 (MS1): Adds precursor intensity columns to psm.tsv
    // Pass 2 (Isobaric): Adds TMT reporter ion intensities to psm.tsv (in-place)
    // Replicates FragPipe's FragpipeRun.java:1978-2031
    //

    // --- TMT IonQuant MS1 Pass ---
    ch_for_ionquant_ms1 = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'ionquant_ms1') }
        .combine(ch_all_results)
        .combine(ch_all_mzml)
        .map { configs, results_wrapped, mzml_wrapped ->
            def args = getToolArgs(configs, 'ionquant_ms1')
            def modmasses = getToolModmasses(configs, 'ionquant_ms1')
            [aggregate_meta, results_wrapped, mzml_wrapped, args, modmasses]
        }

    IONQUANT_MS1(
        ch_for_ionquant_ms1.map { meta, results, mzml, _args, _modmasses -> [meta, results, mzml] },
        ch_annotation_file,
        ch_for_ionquant_ms1.map { _meta, _results, _mzml, args, _modmasses -> args },
        ch_for_ionquant_ms1.map { _meta, _results, _mzml, _args, modmasses -> modmasses },
        []
    )

    // --- TMT IonQuant Isobaric Pass ---
    // IonQuant MS1 pass copies per-sample directories and adds MS1 intensity columns
    // to psm.tsv within the copies. The isobaric pass MUST use the MS1-updated dirs
    // (from IONQUANT_MS1 output), not the original philosopher dirs.
    //
    // IONQUANT_MS1.out.results_dir emits [meta, prefix_dir] where prefix_dir contains
    // per-sample subdirectories with MS1-updated psm.tsv files. We extract those subdirs
    // to pass as psm_dirs input to the isobaric pass.
    //
    // When no annotation file is provided, the ionquant module script handles this gracefully:
    // it disables isoquant (--perform-isoquant 0) and exits cleanly when both quantification
    // modes are disabled (both --perform-isoquant 0 and --perform-ms1quant 0).
    ch_ms1_results = IONQUANT_MS1.out.results_dir
        .map { _meta, dir ->
            // Extract per-sample subdirectories from MS1 output (dirs containing psm.tsv)
            // Use NIO Files.list() instead of toFile().listFiles() — toFile() throws
            // UnsupportedOperationException on S3-backed paths (AWS Batch / Tower).
            def subdirs = []
            java.nio.file.Files.list(dir).withCloseable { stream ->
                subdirs = stream.collect(java.util.stream.Collectors.toList())
                    .findAll { java.nio.file.Files.isDirectory(it) && java.nio.file.Files.exists(it.resolve('psm.tsv')) }
                    .sort { it.toString() }
            }
            [subdirs]  // Wrap in list to prevent combine from flattening
        }

    ch_for_ionquant_isobaric = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'ionquant_isobaric') }
        .combine(ch_ms1_results)  // Data dependency + MS1-updated dirs
        .combine(ch_all_mzml)
        .map { configs, results_wrapped, mzml_wrapped ->
            def args = getToolArgs(configs, 'ionquant_isobaric')
            def modmasses = getToolModmasses(configs, 'ionquant_isobaric')
            [aggregate_meta, results_wrapped, mzml_wrapped, args, modmasses]
        }

    IONQUANT_ISOBARIC(
        ch_for_ionquant_isobaric.map { meta, results, mzml, _args, _modmasses -> [meta, results, mzml] },
        ch_annotation_file,
        ch_for_ionquant_isobaric.map { _meta, _results, _mzml, args, _modmasses -> args },
        ch_for_ionquant_isobaric.map { _meta, _results, _mzml, _args, modmasses -> modmasses },
        []
    )

    //
    // STEP 3: FreeQuant (simple label-free or TMT Philosopher mode, per-sample)
    // In LFQ workflows: simple label-free quantification
    // In TMT Philosopher mode (intensity_extraction_tool=1): adds MS1 intensities
    //   before Labelquant extracts TMT reporter ions (FragpipeRun.java:2049-2055)
    //
    ch_results_keyed = ch_results_dir
        .map { meta, dir -> [meta.id, meta, dir] }

    ch_mzml_keyed = ch_mzml
        .map { meta, mzml -> [meta.id, mzml] }

    // FreeQuant runs when either freequant is enabled OR labelquant is enabled
    // (Philosopher TMT mode requires FreeQuant before Labelquant)
    ch_for_freequant = ch_results_keyed
        .join(ch_mzml_keyed, by: 0, failOnMismatch: false, failOnDuplicate: false)
        .combine(ch_tool_configs)
        .filter { key, meta, dir, mzml, configs ->
            shouldRunTool(configs, 'freequant') || shouldRunTool(configs, 'labelquant')
        }
        .map { key, meta, dir, mzml, configs ->
            def args = getToolArgs(configs, 'freequant')
            // FreeQuant needs directory containing mzML, not the file directly
            def mzml_dir = mzml.parent
            [meta, dir, mzml_dir, args]
        }

    FREEQUANT(ch_for_freequant)

    //
    // STEP 4: Philosopher Labelquant (TMT Philosopher mode, per-sample)
    // Extracts TMT reporter ion intensities from mzML into .meta/ workspace.
    // Runs AFTER FreeQuant (which adds MS1 intensities) and BEFORE Report.
    // Only active when intensity_extraction_tool=1 (CmdLabelquant.java).
    //
    ch_for_labelquant = FREEQUANT.out.results_dir
        .map { meta, dir -> [meta.id, meta, dir] }
        .join(
            ch_mzml.map { meta, mzml -> [meta.id, mzml] },
            by: 0, failOnMismatch: false, failOnDuplicate: false
        )
        .combine(ch_tool_configs)
        .filter { key, meta, dir, mzml, configs -> shouldRunTool(configs, 'labelquant') }
        .combine(ch_annotation)
        .map { key, meta, dir, mzml, configs, annot ->
            def args = getToolArgs(configs, 'labelquant')
            // Labelquant needs directory containing mzML files
            def mzml_dir = mzml.parent
            [meta, dir, mzml_dir, annot, args]
        }

    // Always invoke PHILOSOPHER_LABELQUANT (empty input = no tasks executed, but .out is defined)
    PHILOSOPHER_LABELQUANT(ch_for_labelquant)

    //
    // STEP 5: Philosopher Report after FreeQuant/Labelquant (per-sample)
    // Regenerates TSVs from .meta/ binaries with correct intensity columns.
    // After FreeQuant-only: produces TSVs with MS1 intensities.
    // After FreeQuant+Labelquant: produces TSVs with TMT reporter ion columns.
    //
    // Use Labelquant output when available, otherwise FreeQuant output
    ch_report_from_labelquant = PHILOSOPHER_LABELQUANT.out.results_dir

    ch_report_from_freequant_only = FREEQUANT.out.results_dir
        .combine(ch_tool_configs)
        .filter { meta, dir, configs ->
            shouldRunTool(configs, 'freequant') && !shouldRunTool(configs, 'labelquant')
        }
        .map { meta, dir, _configs -> [meta, dir] }

    ch_for_report = ch_report_from_labelquant
        .mix(ch_report_from_freequant_only)
        .combine(ch_tool_configs)
        .map { meta, dir, configs ->
            def report_args = getToolReportArgs(configs, 'filter')
            [meta, dir, report_args]
        }

    // Always invoke PHILOSOPHER_REPORT (empty input = no tasks executed, but .out is defined)
    PHILOSOPHER_REPORT(ch_for_report, ch_fasta)

    //
    // STEP 6: TMTIntegrator (TMT quantification, aggregate)
    // Uses per-sample psm.tsv/protein.tsv from philosopher directories, which now
    // contain TMT reporter ion columns added by:
    //   - IonQuant isobaric pass (intensity_extraction_tool=0), OR
    //   - Philosopher Labelquant + Report (intensity_extraction_tool=1)
    // Gated on the appropriate upstream completion to ensure TMT columns are present.
    //

    // Gate on IonQuant isobaric pass completion (IonQuant extraction mode)
    ch_isobaric_done = IONQUANT_ISOBARIC.out.results_dir.map { _meta, _dir -> true }.collect()

    // Gate on Labelquant -> Report completion (Philosopher extraction mode)
    ch_labelquant_report_done = PHILOSOPHER_REPORT.out.results_dir.map { _meta, _dir -> true }.collect()

    // Pass directories WITH TMT intensity columns to TMTIntegrator:
    // For IonQuant extraction mode, use IONQUANT_ISOBARIC output (contains per-sample dirs
    //   with "Intensity " columns added by IonQuant isobaric pass).
    // For Philosopher extraction mode, use PHILOSOPHER_REPORT output dirs (contain TMT columns).
    //
    // IonQuant path: use IONQUANT_ISOBARIC output directory (single dir containing per-sample subdirs).
    // TMTIntegrator's script uses find to locate psm.tsv at any depth.
    ch_tmt_results_from_ionquant = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'ionquant_isobaric') }
        .combine(IONQUANT_ISOBARIC.out.results_dir.map { _meta, dir -> dir })
        .map { _configs, ionquant_dir ->
            [[ionquant_dir]]  // Wrap: outer list for concat, inner list for path collection
        }

    // Labelquant path: collect PHILOSOPHER_REPORT output dirs, wrapped in outer list.
    // filter removes empty emission from .collect() when PHILOSOPHER_REPORT didn't run.
    ch_tmt_results_from_labelquant = PHILOSOPHER_REPORT.out.results_dir
        .map { _meta, dir -> dir }
        .collect(sort: true)
        .filter { items -> items.size() > 0 }
        .map { items -> [items] }  // Wrap to prevent combine from flattening

    // Use whichever source provides results (only one path will be active per run).
    // concat+first ensures exactly one emission: ionquant results take priority,
    // and labelquant results are used only if ionquant path was not active.
    ch_tmt_results = ch_tmt_results_from_ionquant
        .concat(ch_tmt_results_from_labelquant)
        .first()

    // TMTIntegrator gate: wait for whichever extraction mode was used.
    // Collapse collected versions to a boolean sentinel (true) so .combine() adds
    // exactly one scalar element instead of flattening a multi-element list.
    // filter removes empty emission from .collect() when the inactive path didn't run.
    ch_extraction_done = ch_isobaric_done
        .mix(ch_labelquant_report_done)
        .filter { items -> items.size() > 0 }
        .map { _items -> true }
        .first()

    ch_for_tmtintegrator = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'tmtintegrator') }
        .combine(ch_extraction_done)
        .combine(ch_tmt_results)
        .map { configs, _gate, results_dirs ->
            [aggregate_meta, results_dirs]
        }

    TMTINTEGRATOR(
        ch_for_tmtintegrator,
        ch_annotation_file,
        ch_tmtintegrator_config
    )

    emit:
    combined_protein = IONQUANT.out.combined_protein           // channel: [ val(meta), path(combined_protein.tsv) ]
    combined_peptide = IONQUANT.out.combined_peptide           // channel: [ val(meta), path(combined_peptide.tsv) ]
    combined_ion     = IONQUANT.out.combined_ions              // channel: [ val(meta), path(combined_ion.tsv) ]
    combined_modified_peptide = IONQUANT.out.combined_modified_peptide  // channel: [ val(meta), path(combined_modified_peptide.tsv) ]
    tmt_abundance    = TMTINTEGRATOR.out.abundance             // channel: [ val(meta), path(abundance_*.tsv) ]
    tmt_ratio        = TMTINTEGRATOR.out.ratio                 // channel: [ val(meta), path(ratio_*.tsv) ]
    ionquant_dir     = IONQUANT.out.results_dir                // channel: [ val(meta), path(results_dir) ]
    tmt_dir          = TMTINTEGRATOR.out.results_dir           // channel: [ val(meta), path(tmt-report) ]
    freequant_ions   = FREEQUANT.out.ions                      // channel: [ val(meta), path(ion.tsv) ]
    labelquant_results_dir = PHILOSOPHER_LABELQUANT.out.results_dir  // channel: [ val(meta), path(results_dir) ] - with TMT intensities in .meta/
    report_results_dir = PHILOSOPHER_REPORT.out.results_dir    // channel: [ val(meta), path(results_dir) ] - TSVs with correct intensities
    report_psms      = PHILOSOPHER_REPORT.out.psms             // channel: [ val(meta), path(psm.tsv) ]
    report_peptides  = PHILOSOPHER_REPORT.out.peptides         // channel: [ val(meta), path(peptide.tsv) ]
    report_proteins  = PHILOSOPHER_REPORT.out.proteins         // channel: [ val(meta), path(protein.tsv) ]
    versions         = IONQUANT.out.versions_ionquant
        .mix(IONQUANT.out.versions_fragpipe)
        .mix(IONQUANT_MS1.out.versions_ionquant)
        .mix(IONQUANT_MS1.out.versions_fragpipe)
        .mix(IONQUANT_ISOBARIC.out.versions_ionquant)
        .mix(IONQUANT_ISOBARIC.out.versions_fragpipe)
        .mix(TMTINTEGRATOR.out.versions_tmtintegrator)
        .mix(TMTINTEGRATOR.out.versions_fragpipe)
        .mix(FREEQUANT.out.versions_philosopher)
        .mix(FREEQUANT.out.versions_fragpipe)
        .mix(PHILOSOPHER_LABELQUANT.out.versions_philosopher)
        .mix(PHILOSOPHER_LABELQUANT.out.versions_fragpipe)
        .mix(PHILOSOPHER_REPORT.out.versions_philosopher)
        .mix(PHILOSOPHER_REPORT.out.versions_fragpipe)
}
