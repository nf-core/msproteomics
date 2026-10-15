process PARSEFRAGPIPEWORKFLOW {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    // data_types: FragPipe data types of the manifest files (DDA, DDA+, DIA, DIA-Quant, DIA-Lib, GPF-DIA).
    //   FragPipe adjusts tools to them after loading the manifest (TabWorkflow.adjustToolsBasedOnDataTypes).
    tuple val(meta), path(workflow_file), val(data_types)
    // templates: FragPipe 24.0's fragger_closed.params and diaumpire_se.params (this module's assets/);
    //   staged because a container sees only the module's resources/usr/bin
    path(templates, stageAs: 'parser_templates/*')

    output:
    // JSON file with all tool configs - primary output for dynamic workflows
    tuple val(meta), path("tool_configs.json"), emit: tool_configs_json
    // MSFragger params per FragPipe data type (CmdMsfragger.java:444-457): DDA, DIA (+DIA-Lib), GPF-DIA, DDA+;
    // the DDA file is fragger_dda.params when the manifest also has DIA, GPF-DIA or DIA-Lib files (:444);
    // always written: without msfragger.* keys they hold FragPipe's closed-search defaults (TabMsfragger.java:515-517)
    tuple val(meta), path("configs/{fragger,fragger_dda}.params"), emit: msfragger_config
    tuple val(meta), path("configs/fragger_dia.params")      , emit: msfragger_dia_config
    tuple val(meta), path("configs/fragger_gpfdia.params")   , emit: msfragger_gpfdia_config
    tuple val(meta), path("configs/fragger_dda_plus.params") , emit: msfragger_dda_plus_config
    tuple val(meta), path("configs/msbooster.config")      , emit: msbooster_config, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def data_types_arg = data_types.toUnique().toSorted().join(',')
    if ([templates].flatten().collect { f -> f.toString().tokenize('/').last() }.sort() != ['diaumpire_se.params', 'fragger_closed.params']) {
        error("PARSEFRAGPIPEWORKFLOW: templates must be fragger_closed.params and diaumpire_se.params, got ${templates}")
    }
    """
    # tool_configs.json (run flags and args of every tool) and the tool config files in configs/, one run
    parse_fragpipe_workflow.py \\
        --workflow ${workflow_file} \\
        --data-types '${data_types_arg}' \\
        --output-json \\
        --outdir . \\
        --configs-outdir configs \\
        --fragger-template parser_templates/fragger_closed.params \\
        --diaumpire-template parser_templates/diaumpire_se.params \\
        ${args}
    """

    stub:
    // The parser is a pure, fast Python script: the stub runs it for real so stub runs
    // follow the workflow file (tools enabled, TMT label type) instead of a canned config.
    def args = task.ext.args ?: ''
    def data_types_arg = data_types.toUnique().toSorted().join(',')
    if ([templates].flatten().collect { f -> f.toString().tokenize('/').last() }.sort() != ['diaumpire_se.params', 'fragger_closed.params']) {
        error("PARSEFRAGPIPEWORKFLOW: templates must be fragger_closed.params and diaumpire_se.params, got ${templates}")
    }
    """
    # tool_configs.json (run flags and args of every tool) and the tool config files in configs/, one run
    parse_fragpipe_workflow.py \\
        --workflow ${workflow_file} \\
        --data-types '${data_types_arg}' \\
        --output-json \\
        --outdir . \\
        --configs-outdir configs \\
        --fragger-template parser_templates/fragger_closed.params \\
        --diaumpire-template parser_templates/diaumpire_se.params \\
        ${args}
    """

}
