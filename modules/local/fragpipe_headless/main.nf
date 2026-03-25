/*
 * FRAGPIPE_HEADLESS: Run FragPipe in headless mode (all-in-one process)
 *
 * Runs the entire FragPipe pipeline (MSFragger, MSBooster, Percolator,
 * ProteinProphet, Philosopher Filter, IonQuant, etc.) in a single process
 * using FragPipe's built-in headless mode.
 *
 * This produces results identical to running FragPipe from the GUI.
 *
 * Input:
 *   - raw_files:         All raw data files (.d, .raw, .mzML)
 *   - database:          FASTA database file
 *   - workflow_file:     FragPipe .workflow configuration file
 *   - manifest_content:  Tab-separated manifest content (filename\texperiment\tbioreplicate\tdata_type)
 *
 * Output:
 *   - all_results:       All FragPipe output files
 *   - combined_protein:  Combined protein report
 *   - combined_peptide:  Combined peptide report
 *   - combined_ion:      Combined ion report
 *   - versions:          Software versions
 */
process FRAGPIPE_HEADLESS {
    tag "fragpipe_headless"
    label 'process_high'

    // FragPipe is proprietary software — not available via conda
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    path raw_files, stageAs: "raw_files/*"
    path database
    path workflow_file
    val  manifest_content

    output:
    path "results/**",                   emit: all_results
    path "results/combined_protein.tsv", emit: combined_protein, optional: true
    path "results/combined_peptide.tsv", emit: combined_peptide, optional: true
    path "results/combined_ion.tsv",     emit: combined_ion,     optional: true
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mem_gb = (task.memory.toGiga() * 0.9).intValue()
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    def lib_dir = tools_dir.replaceAll('/tools$', '/lib')
    """
    export JAVA_OPTS="-Xmx${mem_gb}G"

    mkdir -p results

    # Generate manifest with absolute paths to staged raw files
    while IFS=\$'\\t' read -r filename experiment biorep dtype; do
        echo -e "\$(pwd)/raw_files/\${filename}\\t\${experiment}\\t\${biorep}\\t\${dtype}"
    done <<< "${manifest_content}" > manifest.fp-manifest

    # Update workflow file for container environment
    cp ${workflow_file} run.workflow
    # Set database path to staged file (replace if exists, append if not)
    if grep -q '^database.db-path=' run.workflow; then
        sed -i "s|^database.db-path=.*|database.db-path=\$(pwd)/${database}|" run.workflow
    else
        echo "database.db-path=\$(pwd)/${database}" >> run.workflow
    fi
    # Add decoys to FASTA if not already present
    DECOY_TAG=\$(grep -oP '^database.decoy-tag=\\K.*' run.workflow || echo "rev_")
    PHILOSOPHER=\$(find ${tools_dir} -name "philosopher*" -type f | head -1)
    if ! grep -q ">\${DECOY_TAG}" ${database} && [ -n "\${PHILOSOPHER}" ]; then
        echo "Adding decoy sequences with tag \${DECOY_TAG}..."
        "\${PHILOSOPHER}" workspace --init --nocheck
        "\${PHILOSOPHER}" database --custom ${database} --contam "" --prefix "\${DECOY_TAG}"
        # philosopher generates a new FASTA with decoys; find and use it
        DECOY_DB=\$(find . -maxdepth 1 -name "*.fas" | head -1)
        if [ -n "\${DECOY_DB}" ]; then
            sed -i "s|^database.db-path=.*|database.db-path=\$(pwd)/\${DECOY_DB}|" run.workflow
        fi
    fi
    # Remove developer-machine tool paths (fragpipe-config.bin-* and legacy keys);
    # --config-tools-folder handles tool discovery inside the container.
    sed -i '/^fragpipe-config\\.bin-/d' run.workflow
    sed -i '/^philosopher\\.exe=/d' run.workflow
    sed -i '/^philospher\\.exe=/d' run.workflow
    sed -i '/^msfragger\\.ext-thermo=/d' run.workflow
    sed -i '/^diann\\.exec-path=/d' run.workflow
    sed -i '/^diann\\.exe=/d' run.workflow

    # Run FragPipe headless
    # FragPipe uses Gradle application plugin — launch via main class + classpath, not -jar.
    java -Djava.awt.headless=true \${JAVA_OPTS} \\
        -cp "${lib_dir}/*" \\
        org.nesvilab.fragpipe.FragPipeMain \\
        --headless \\
        --workflow run.workflow \\
        --manifest manifest.fp-manifest \\
        --workdir \$(pwd)/results \\
        --threads ${task.cpus} \\
        --ram ${mem_gb} \\
        --config-tools-folder ${tools_dir}

    # --help prints "FragPipePlus v24.0" on the first line
    FRAGPIPE_VERSION=\$(java -cp "${lib_dir}/*" org.nesvilab.fragpipe.FragPipeMain --help 2>&1 | grep -oP 'v\\K[0-9.]+' | head -1 || true)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragpipe: "\${FRAGPIPE_VERSION}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p results/sample1
    touch results/combined_protein.tsv
    touch results/combined_peptide.tsv
    touch results/combined_ion.tsv
    touch results/sample1/psm.tsv
    touch results/sample1/protein.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragpipe: "24.0"
    END_VERSIONS
    """
}
