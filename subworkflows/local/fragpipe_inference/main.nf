/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_INFERENCE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Protein inference and FDR filtering.

    Modules:
    - PROTEINPROPHET (aggregate protein inference)
    - PHILOSOPHER_FILTER (per-sample FDR filtering and reporting)

    Execution: Mixed (aggregate → per-sample)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROTEINPROPHET    } from '../../../modules/local/philosopher/proteinprophet/main'
include { PHILOSOPHER_FILTER } from '../../../modules/local/philosopher/filter/main'
include { shouldRunTool; getToolArgs; getToolReportArgs } from '../fragpipe_utils'

workflow FRAGPIPE_INFERENCE {
    take:
    ch_pepxml        // channel: [ val(meta), path(pepxml) ] - all samples
    ch_fasta         // channel: path(fasta) - value channel for database
    ch_tool_configs  // channel: val(tool_configs_map) - parsed JSON config
    aggregate_meta   // val(meta) for aggregate outputs

    main:

    //
    // STEP 1: ProteinProphet (aggregate protein inference)
    // Combines all pepXML files into a single protein-level analysis
    //
    ch_all_pepxml = ch_pepxml
        .map { meta, pepxml -> pepxml }
        .collect(sort: true)
        .map { items -> [items] }  // Wrap to prevent combine from flattening

    ch_for_proteinprophet = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'proteinprophet') }
        .combine(ch_all_pepxml)
        .map { configs, pepxml_wrapped ->
            def args = getToolArgs(configs, 'proteinprophet')
            [aggregate_meta, pepxml_wrapped, args]
        }

    PROTEINPROPHET(ch_for_proteinprophet)

    //
    // STEP 2: Philosopher Filter (per-group FDR filtering)
    // Groups pepXMLs by meta.group for fractionation support:
    //   Non-fractionated: meta.group == meta.id → each group has one pepXML (same as per-sample)
    //   Fractionated: meta.group groups fractions → filter runs on all fractions together
    // This matches FragPipe's CmdPhilosopherFilter.java:60-122 which iterates
    // mapGroupsToProtxml.entrySet() (one filter invocation per experiment group).
    //
    ch_pepxml_grouped = ch_pepxml
        .map { meta, pepxml -> [meta.group ?: meta.id, pepxml] }
        .groupTuple(sort: true)
        .map { group, pepxmls ->
            // Flatten nested lists: each sample may contribute a list of ranked pepXMLs
            // (e.g., DDA+ ion mobility: interact-sample_rank1.pep.xml through _rank5.pep.xml).
            // groupTuple wraps each emission in a list, so [[rank1, rank2, ...]] needs flattening.
            def flat_pepxmls = pepxmls.flatten()
            [[id: group, group: group], flat_pepxmls]
        }

    ch_for_filter = ch_pepxml_grouped
        .combine(PROTEINPROPHET.out.protxml.map { _meta, protxml -> protxml })
        .combine(ch_tool_configs)
        .filter { meta, pepxmls, protxml, configs -> shouldRunTool(configs, 'filter') }
        .map { meta, pepxmls, protxml, configs ->
            def filter_args = getToolArgs(configs, 'filter')
            // --razor is enforced by parse_fragpipe_workflow.py (matching CmdPhilosopherFilter.java:84-114)
            def report_args = getToolReportArgs(configs, 'filter')
            [meta, pepxmls, protxml, filter_args, report_args]
        }

    PHILOSOPHER_FILTER(ch_for_filter, ch_fasta)

    emit:
    protxml      = PROTEINPROPHET.out.protxml                  // channel: [ val(meta), path(combined.prot.xml) ]
    results_dir  = PHILOSOPHER_FILTER.out.results_dir          // channel: [ val(meta), path(results_dir) ] per-sample
    psms         = PHILOSOPHER_FILTER.out.psms                 // channel: [ val(meta), path(psm.tsv) ]
    peptides     = PHILOSOPHER_FILTER.out.peptides             // channel: [ val(meta), path(peptide.tsv) ]
    proteins     = PHILOSOPHER_FILTER.out.proteins             // channel: [ val(meta), path(protein.tsv) ]
    ions         = PHILOSOPHER_FILTER.out.ions                 // channel: [ val(meta), path(ion.tsv) ] optional
    pepxml       = ch_pepxml                                   // channel: [ val(meta), path(pepxml) ] - ProteinProphet doesn't modify pepXML
    versions     = PROTEINPROPHET.out.versions_philosopher
        .mix(PROTEINPROPHET.out.versions_fragpipe)
        .mix(PHILOSOPHER_FILTER.out.versions_philosopher)
        .mix(PHILOSOPHER_FILTER.out.versions_fragpipe)
}
