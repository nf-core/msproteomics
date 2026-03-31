process PARSE_FRAGPIPE_WORKFLOW {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(workflow_file)

    output:
    // JSON file with all tool configs - primary output for dynamic workflows
    tuple val(meta), path("tool_configs.json"), emit: tool_configs_json
    // Legacy config file outputs for backwards compatibility
    tuple val(meta), path("configs/fragger.params")        , emit: msfragger_config, optional: true
    tuple val(meta), path("configs/ionquant.config")       , emit: ionquant_config, optional: true
    tuple val(meta), path("configs/msbooster.config")      , emit: msbooster_config, optional: true
    tuple val(meta), path("configs/percolator.config")     , emit: percolator_config, optional: true
    tuple val(meta), path("configs/phi-report.config")     , emit: phi_report_config, optional: true
    tuple val(meta), path("configs/peptide-prophet.config"), emit: peptideprophet_config, optional: true
    tuple val(meta), path("configs/protein-prophet.config"), emit: proteinprophet_config, optional: true
    tuple val(meta), path("configs/database.config")       , emit: database_config, optional: true
    tuple val(meta), path("configs/tmtintegrator.yml")     , emit: tmtintegrator_config, optional: true
    tuple val(meta), path("configs")                       , emit: all_configs
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash

    mkdir -p configs

    # Generate JSON output with run flags and args for all tools
    parse_fragpipe_workflow.py \\
        --workflow ${workflow_file} \\
        --output-json \\
        --outdir .

    # Also generate legacy config files for backwards compatibility
    parse_fragpipe_workflow.py \\
        --workflow ${workflow_file} \\
        --outdir configs \\
        ${args}

    # Python now outputs fragger.params directly (no rename needed)

    # Ensure all expected config files exist (create empty ones if missing)
    touch configs/fragger.params
    touch configs/ionquant.config
    touch configs/msbooster.config
    touch configs/percolator.config
    touch configs/phi-report.config
    touch configs/peptide-prophet.config
    touch configs/protein-prophet.config
    touch configs/database.config
    touch configs/tmtintegrator.yml

    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p configs

    # Create stub JSON file with example tool configs
    cat <<-'EOF' > tool_configs.json
{
  "msfragger": {"run": true, "args": "precursor_mass_lower = -20\nprecursor_mass_upper = 20", "config_type": "params_file"},
  "msbooster": {"run": true, "args": "useRT = true\nrtModel = DIA-NN", "config_type": "params_file"},
  "percolator": {"run": true, "args": "--only-psms --no-terminate", "config_type": "cli"},
  "peptideprophet": {"run": false, "args": "", "config_type": "cli"},
  "proteinprophet": {"run": true, "args": "proteinprophet=--maxppmdiff 2000000", "config_type": "cli"},
  "filter": {"run": true, "args": "filter=--sequential --prot 0.01 --picked", "config_type": "cli"},
  "ionquant": {"run": true, "args": "--mbr 1 --maxlfq 1 --perform-ms1quant 1 --perform-isoquant 0 --site-reports 1 --msstats 1 --ionmobility 0", "config_type": "cli", "modmasses": "15.9949,42.0106"},
  "tmtintegrator": {"run": false, "args": "", "config_type": "params_file"},
  "freequant": {"run": false, "args": "", "config_type": "cli"},
  "diaumpire": {"run": false, "args": "", "config_type": "cli"},
  "diatracer": {"run": false, "args": "", "config_type": "cli"},
  "diann": {"run": false, "args": "", "config_type": "cli"},
  "speclibgen": {"run": false, "args": "", "config_type": "cli"},
  "ptmprophet": {"run": false, "args": "", "config_type": "cli"},
  "ptmshepherd": {"run": false, "args": "", "config_type": "cli"},
  "crystalc": {"run": false, "args": "", "config_type": "cli"},
  "opair": {"run": false, "args": "", "config_type": "cli"},
  "mbg": {"run": false, "args": "", "config_type": "cli"},
  "fpop": {"run": false, "args": "", "config_type": "cli"},
  "saintexpress": {"run": false, "args": "", "config_type": "cli"},
  "skyline": {"run": false, "args": "", "config_type": "cli"},
  "metaproteomics": {"run": false, "args": "", "config_type": "cli"},
  "transferlearning": {"run": false, "args": "", "config_type": "cli"},
  "database": {"run": true, "args": "database=--prefix rev_", "config_type": "cli"}
}
EOF

    # Create stub config files (params file format for MSFragger/MSBooster, CLI for others)
    printf "precursor_mass_lower = -20\nprecursor_mass_upper = 20\n" > configs/fragger.params
    echo "--mbr 1 --maxlfq 1 --requantify 1 --mztol 10" > configs/ionquant.config
    printf "useDetect = false\nuseRT = true\nrtModel = DIA-NN\nspectraModel = DIA-NN\n" > configs/msbooster.config
    echo "--only-psms --no-terminate" > configs/percolator.config
    echo "filter=--sequential --prot 0.01 --picked" > configs/phi-report.config
    echo "peptideprophet=--decoyprobs --ppm" > configs/peptide-prophet.config
    echo "proteinprophet=--maxppmdiff 2000000" > configs/protein-prophet.config
    echo "database=--prefix rev_" > configs/database.config

    """
}
