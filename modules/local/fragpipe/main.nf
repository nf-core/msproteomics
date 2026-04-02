process FRAGPIPE {
    tag "${meta.id}"
    label 'process_high'

    // This module needs the full fragpipe container.

    input:
    tuple val(meta), path(mzml_files)
    path(fasta)
    path(workflow_file)
    val(data_type)  // 'DDA', 'DIA', 'GPF-DIA', 'DIA-Quant', 'DIA-Lib', 'DDA+'

    output:
    tuple val(meta), path("${prefix}")                 , emit: results_dir
    tuple val(meta), path("${prefix}/*/protein.tsv")   , emit: proteins
    tuple val(meta), path("${prefix}/*/peptide.tsv")   , emit: peptides
    tuple val(meta), path("${prefix}/*/psm.tsv")       , emit: psms
    tuple val(meta), path("${prefix}/*/ion.tsv")       , emit: ions        , optional: true
    tuple val(meta), path("${prefix}/*/*.pepXML")      , emit: pepxml      , optional: true
    tuple val(meta), path("${prefix}/combined.prot.xml"), emit: protxml    , optional: true
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def threads = task.cpus ?: 16

    // Get mzML file list - handle both single file and list
    def mzml_list = mzml_files instanceof List ? mzml_files : [mzml_files]
    def num_samples = mzml_list.size()

    """
    
    #!/bin/bash

    echo "============================================================"
    echo "FragPipe All-in-One Workflow"
    echo "============================================================"
    echo "Experiment: ${prefix}"
    echo "Samples: ${num_samples}"
    echo "Data type: ${data_type}"
    echo "Workflow: ${workflow_file}"
    echo "RAM: ${ram}G"
    echo "Threads: ${threads}"
    echo "============================================================"
    echo ""

    # Set HOME to current directory for FragPipe config
    export HOME=\$(pwd)
    mkdir -p \$HOME/.config

    # Create results directory
    mkdir -p ${prefix}

    # ========================================================================
    # STEP 1: Copy input files to local directory
    # Required for mzBIN cache writes on network mounts
    # ========================================================================
    echo "Copying input files to local directory..."
    mkdir -p mzml_local fasta_local

    # Copy mzML files
    for f in ${mzml_list.join(' ')}; do
        cp -L "\$f" mzml_local/
    done
    echo "  Copied \$(find mzml_local -maxdepth 1 -name '*.mzML' | wc -l) mzML files"

    # Copy FASTA to local (MSFragger writes index files next to FASTA)
    cp -L ${fasta} fasta_local/database.fasta
    echo "  Copied FASTA database"
    echo ""

    # ========================================================================
    # STEP 2: Create manifest file
    # ========================================================================
    echo "Creating manifest file..."
    python3 /usr/local/bin/create_manifest.py \\
        --mzml_files mzml_local/*.mzML \\
        --experiment "${prefix}" \\
        --data_type "${data_type}" \\
        --output manifest.tsv

    echo ""
    cat manifest.tsv
    echo ""

    # ========================================================================
    # STEP 3: Prepare workflow file
    # ========================================================================
    echo "Preparing workflow file..."
    python3 /usr/local/bin/prepare_workflow.py \\
        --workflow ${workflow_file} \\
        --fasta fasta_local/database.fasta \\
        --output fragpipe.workflow \\
        --threads ${threads}
    echo ""

    # ========================================================================
    # STEP 4: Run FragPipe
    # ========================================================================
    echo "Running FragPipe..."
    echo ""

    # Initialize log file
    LOG_FILE="${prefix}/${prefix}_fragpipe.log"
    echo "=== FragPipe ===" > \$LOG_FILE
    echo "Timestamp: \$(date)" >> \$LOG_FILE
    echo "Experiment: ${prefix}" >> \$LOG_FILE
    echo "Workflow: ${workflow_file}" >> \$LOG_FILE
    echo "Samples: ${num_samples}" >> \$LOG_FILE
    echo "" >> \$LOG_FILE

    # Version capture
    _fp_ver=\$(fragpipe --version 2>&1 | grep -oP 'FragPipe version \\K[0-9.]+' || true)

    echo "\$_fp_ver" > .fragpipe_version

    # Define the main command
    CMD="/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/bin/fragpipe --headless --workflow fragpipe.workflow --manifest manifest.tsv --workdir ${prefix} --config-tools-folder /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools --config-diann /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 --config-python /usr/bin/python3 --ram ${ram} --threads ${threads} ${args}"

    echo "=== Command ===" >> \$LOG_FILE
    echo "\$CMD" >> \$LOG_FILE
    echo "" >> \$LOG_FILE
    echo "=== Output ===" >> \$LOG_FILE

    # Run FragPipe
    eval "\$CMD" 2>&1 | tee -a \$LOG_FILE

    echo "" >> \$LOG_FILE
    echo "=== Results ===" >> \$LOG_FILE
    find ${prefix} -name "*.tsv" -o -name "*.xml" | head -20 >> \$LOG_FILE

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/${prefix}
    touch ${prefix}/${prefix}/protein.tsv
    touch ${prefix}/${prefix}/peptide.tsv
    touch ${prefix}/${prefix}/psm.tsv
    touch ${prefix}/${prefix}/ion.tsv
    touch ${prefix}/fragpipe.log
    echo "24.0" > .fragpipe_version
    """
}
