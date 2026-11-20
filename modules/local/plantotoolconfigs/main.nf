/*
 * PLANTOTOOLCONFIGS: tool configs from FragPipe's configured run (FRAGPIPEPLAN fp_work/)
 *
 * Runs bin/plan_to_tool_configs.py, which writes the tool configs the FragPipe
 * subworkflows read: tool_configs.json and the per-data-type MSFragger params and MSBooster params in configs/.
 *
 * Input:
 *   - meta: plan id map
 *   - workdir: FragPipe's work dir after --write-plan (FRAGPIPEPLAN): run log, configs, run files
 *
 * Output:
 *   - tool_configs_json: tool run flags, args and fields (one entry per tool: run, args, config_type and tool-specific fields)
 *   - msfragger_*_config: MSFragger params of each data type the run has (FragPipe's files)
 *   - msbooster_config:   MSBooster params without the run-time keys (empty when MSBooster does not run)
 *   - versions_python:    Python version
 */
process PLANTOTOOLCONFIGS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/370b76e81fbff1ab2427ad3a246bca8a4aa578a68216a5e800dceaedf7bcfa40/data' :
        'community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658' }"

    input:
    tuple val(meta), path(workdir, stageAs: 'fp_work')

    output:
    tuple val(meta), path("tool_configs.json")                    , emit: tool_configs_json
    tuple val(meta), path("configs/{fragger,fragger_dda}.params") , emit: msfragger_config, optional: true
    tuple val(meta), path("configs/fragger_dia.params")           , emit: msfragger_dia_config, optional: true
    tuple val(meta), path("configs/fragger_gpfdia.params")        , emit: msfragger_gpfdia_config, optional: true
    tuple val(meta), path("configs/fragger_dda_plus.params")      , emit: msfragger_dda_plus_config, optional: true
    tuple val(meta), path("configs/msbooster.config")             , emit: msbooster_config
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    plan_to_tool_configs.py \\
        --workdir ${workdir} \\
        --outdir . \\
        ${args}
    """

    stub:
    // plan_to_tool_configs.py is a pure, fast Python script: the stub runs it for real
    def args = task.ext.args ?: ''
    """
    plan_to_tool_configs.py \\
        --workdir ${workdir} \\
        --outdir . \\
        ${args}
    """
}
