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
 *   - raw_files:            All raw data files (.d, .raw, .mzML)
 *   - database:             FASTA database file
 *   - workflow_file:        FragPipe .workflow configuration file
 *   - manifest_content:     Tab-separated manifest content (filename\texperiment\tbioreplicate\tdata_type)
 *   - annotation_content:   TMT annotation content (experiment\tchannel\tsample_name), empty for LFQ
 *   - file_experiment_map:  File-to-experiment mapping (filename\texperiment), empty for LFQ
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

    // No default container — users must provide a licensed FragPipe image
    // via process.withName or process.container in their config

    input:
    path raw_files, stageAs: "raw_files/*"
    path database
    path workflow_file
    val  manifest_content
    val  annotation_content
    val  file_experiment_map

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

    ANNOTATION_CONTENT="${annotation_content}"
    FILE_EXP_MAP="${file_experiment_map}"

    # Determine TMT mode
    IS_MULTIPLEX=false
    if [ -n "\${ANNOTATION_CONTENT}" ]; then
        NUM_TMT_EXP=\$(echo "\${ANNOTATION_CONTENT}" | cut -f1 | sort -u | wc -l)
        [ "\${NUM_TMT_EXP}" -gt 1 ] && IS_MULTIPLEX=true
    fi

    if [ "\${IS_MULTIPLEX}" = true ]; then
        #
        # MULTI-PLEX TMT: Per-experiment subdirectories with annotation files.
        # FragPipe auto-discovers *annotation.txt in each mzML's parent directory.
        # Requires exactly 1 annotation file per directory (TmtiPanel.java:868-880).
        #
        EXPERIMENTS=\$(echo "\${ANNOTATION_CONTENT}" | cut -f1 | sort -u)

        # Create per-experiment subdirectories and move files
        for EXP in \${EXPERIMENTS}; do
            mkdir -p "raw_files/\${EXP}"
        done

        for f in raw_files/*.mzML raw_files/*.d; do
            [ -e "\$f" ] || continue
            fname=\$(basename "\$f")
            # Look up experiment from file-experiment map
            TARGET_EXP=\$(echo "\${FILE_EXP_MAP}" | grep -F "\${fname}" | head -1 | cut -f2)
            if [ -n "\${TARGET_EXP}" ] && [ -d "raw_files/\${TARGET_EXP}" ]; then
                cp "\$f" "raw_files/\${TARGET_EXP}/\${fname}"
            fi
        done

        # Create per-experiment annotation files (space-separated: channel sample_name)
        for EXP in \${EXPERIMENTS}; do
            echo "\${ANNOTATION_CONTENT}" | awk -F'\\t' -v expname="\${EXP}" '\$1==expname {print \$2" "\$3}' > "raw_files/\${EXP}/\${EXP}_annotation.txt"
        done

        # Generate manifest pointing into subdirectories
        >| manifest.fp-manifest
        for EXP in \${EXPERIMENTS}; do
            for f in raw_files/\${EXP}/*.mzML raw_files/\${EXP}/*.d; do
                [ -e "\$f" ] || continue
                fname=\$(basename "\$f")
                [[ "\$fname" == *annotation.txt ]] && continue
                echo -e "\$(pwd)/raw_files/\${EXP}/\${fname}\\t\${EXP}\\t1\\tDDA"
            done
        done >> manifest.fp-manifest

    elif [ -n "\${ANNOTATION_CONTENT}" ]; then
        #
        # SINGLE-PLEX TMT: One annotation file in raw_files/ (no subdirs needed).
        # All files from one experiment share the directory.
        #
        EXP=\$(echo "\${ANNOTATION_CONTENT}" | cut -f1 | sort -u | head -1)
        echo "\${ANNOTATION_CONTENT}" | awk -F'\\t' '{print \$2" "\$3}' > "raw_files/\${EXP}_annotation.txt"

        >| manifest.fp-manifest
        for f in raw_files/*.mzML raw_files/*.d; do
            [ -e "\$f" ] || continue
            fname=\$(basename "\$f")
            [[ "\$fname" == *annotation.txt ]] && continue
            echo -e "\$(pwd)/raw_files/\${fname}\\t\${EXP}\\t1\\tDDA"
        done >> manifest.fp-manifest

    else
        #
        # LFQ MODE: No annotation, flat directory.
        #
        FILE_EXP_MAP_FOR_LFQ="${file_experiment_map}"
        >| manifest.fp-manifest
        for f in raw_files/*.mzML raw_files/*.d; do
            [ -e "\$f" ] || continue
            fname=\$(basename "\$f")
            MATCHED_EXP=""
            if echo "\${FILE_EXP_MAP_FOR_LFQ}" | grep -qF "\${fname}"; then
                MATCHED_EXP=\$(echo "\${FILE_EXP_MAP_FOR_LFQ}" | grep -F "\${fname}" | head -1 | cut -f2)
            fi
            echo -e "\$(pwd)/raw_files/\${fname}\\t\${MATCHED_EXP:-experiment1}\\t1\\tDDA"
        done >> manifest.fp-manifest
    fi

    # Update workflow file for container environment
    cp ${workflow_file} run.workflow
    # Set database path to staged file (replace if exists, append if not)
    if grep -q '^database.db-path=' run.workflow; then
        sed -i "s|^database.db-path=.*|database.db-path=\$(pwd)/${database}|" run.workflow
    else
        echo "database.db-path=\$(pwd)/${database}" >> run.workflow
    fi
    # FragPipe headless handles decoy generation internally via Philosopher.
    # Do NOT run philosopher here — it corrupts FASTA on Fusion-mounted filesystems.
    # Set explicit tool paths in workflow file for headless mode.
    # Remove any existing tool path configs first, then discover and inject actual paths.
    sed -i '/^fragpipe-config\\.bin-/d' run.workflow
    sed -i '/^philosopher\\.exe=/d' run.workflow
    sed -i '/^philospher\\.exe=/d' run.workflow
    sed -i '/^msfragger\\.ext-thermo=/d' run.workflow
    sed -i '/^diann\\.exec-path=/d' run.workflow
    sed -i '/^diann\\.exe=/d' run.workflow

    # Discover and inject tool paths explicitly (avoids async discovery issues in headless mode)
    MSFRAGGER_JAR=\$(find ${tools_dir} -name "MSFragger*.jar" -not -name "*original*" | head -1)
    IONQUANT_JAR=\$(find ${tools_dir} -name "IonQuant*.jar" | head -1)
    DIATRACER_JAR=\$(find ${tools_dir} -name "diaTracer*.jar" -o -name "DiaTracer*.jar" | head -1)
    DIANN_BIN=\$(find ${tools_dir} -path "*/diann/*" -name "diann-*" -type f | head -1)
    PHILOSOPHER=\$(find ${tools_dir} -path "*/Philosopher/*" -name "philosopher*" -type f | head -1)
    [ -n "\${MSFRAGGER_JAR}" ] && echo "fragpipe-config.bin-msfragger=\${MSFRAGGER_JAR}" >> run.workflow
    [ -n "\${IONQUANT_JAR}" ] && echo "fragpipe-config.bin-ionquant=\${IONQUANT_JAR}" >> run.workflow
    [ -n "\${DIATRACER_JAR}" ] && echo "fragpipe-config.bin-diatracer=\${DIATRACER_JAR}" >> run.workflow
    [ -n "\${PHILOSOPHER}" ] && echo "philosopher.exe=\${PHILOSOPHER}" >> run.workflow
    [ -n "\${DIANN_BIN}" ] && echo "diann.exec-path=\${DIANN_BIN}" >> run.workflow

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
    cat <<-END_VERSIONS >| versions.yml
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

    cat <<-END_VERSIONS >| versions.yml
    "${task.process}":
        fragpipe: "24.0"
    END_VERSIONS
    """
}
