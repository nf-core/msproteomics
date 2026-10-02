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
 *   - raw_files:            All raw data files (.d, .mzML) staged into raw_files/
 *   - database:             FASTA database file
 *   - workflow_file:        FragPipe .workflow configuration file
 *   - manifest_content:     Manifest rows (file name\texperiment\tbioreplicate\tdata_type).
 *                           Authoritative when non-empty; bioreplicate may be empty (fractions)
 *   - annotation_content:   TMT annotation content (experiment\tchannel\tsample_name), empty for LFQ
 *   - file_experiment_map:  File-to-experiment mapping (file name\texperiment), used only when
 *                           manifest_content is empty (bioreplicate 1, data type DDA)
 *
 * File names are matched exactly against the first column. The task fails if a staged
 * .mzML/.d file has no row, a row names a file that is not staged, a file name is listed
 * twice, a TMT file's experiment is not in annotation_content, or manifest_content and
 * file_experiment_map are both empty.
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
    def annotation_b64 = annotation_content.toString().getBytes('UTF-8').encodeBase64().toString()
    """
    export JAVA_OPTS="-Xmx${mem_gb}G"

    mkdir -p results

    ${fragpipeLookupScript(manifest_content, file_experiment_map)}

    printf '%s' '${annotation_b64}' | base64 -d > input_annotation.tsv

    if [ -s input_annotation.tsv ]; then
        #
        # TMT: one subdirectory per plex (experiment) holding its mzML files and annotation.
        # FragPipe auto-discovers *annotation.txt in each mzML's parent directory and
        # requires exactly 1 annotation file per directory (TmtiPanel.java:868-880).
        #
        cut -f1 input_annotation.tsv | sort -u > tmt_experiments.txt
        UNKNOWN=\$(awk -F'\\t' 'NR == FNR { plex[\$0] = 1; next } !(\$2 in plex) { print \$1 " -> " \$2 }' tmt_experiments.txt lookup.tsv)
        if [ -n "\${UNKNOWN}" ]; then
            echo "ERROR: FRAGPIPE_HEADLESS: \${LOOKUP_SOURCE} experiment(s) not in annotation_content: \${UNKNOWN}" >&2
            exit 1
        fi

        while IFS= read -r EXP; do
            mkdir -p "raw_files/\${EXP}"
            # space-separated: channel sample_name
            EXP="\${EXP}" awk -F'\\t' '\$1 == ENVIRON["EXP"] { print \$2" "\$3 }' input_annotation.tsv > "raw_files/\${EXP}/\${EXP}_annotation.txt"
        done < tmt_experiments.txt

        cut -f1,2 lookup.tsv | while IFS=\$'\\t' read -r fname EXP; do
            cp -R "raw_files/\${fname}" "raw_files/\${EXP}/\${fname}"
        done

        RAW_DIR="\$(pwd)/raw_files" awk -F'\\t' -v OFS='\\t' '{ print ENVIRON["RAW_DIR"] "/" \$2 "/" \$1, \$2, \$3, \$4 }' lookup.tsv > manifest.fp-manifest
    else
        #
        # LFQ: flat directory.
        #
        RAW_DIR="\$(pwd)/raw_files" awk -F'\\t' -v OFS='\\t' '{ print ENVIRON["RAW_DIR"] "/" \$1, \$2, \$3, \$4 }' lookup.tsv > manifest.fp-manifest
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
    ${fragpipeLookupScript(manifest_content, file_experiment_map)}

    mkdir -p results/sample1
    touch results/combined_protein.tsv
    touch results/combined_peptide.tsv
    touch results/combined_ion.tsv
    touch results/sample1/psm.tsv
    touch results/sample1/protein.tsv

    # Manifest and LFQ experiment annotation as FragPipe writes them into the work dir
    # (ToolingUtils.generateLFQExperimentAnnotation); relative paths keep stub output deterministic
    awk -F'\\t' -v OFS='\\t' '{ print "raw_files/" \$1, \$2, \$3, \$4 }' lookup.tsv > results/fragpipe-files.fp-manifest
    awk -F'\\t' -v OFS='\\t' '
        BEGIN { print "file", "sample", "sample_name", "condition", "replicate" }
        {
            sample = (\$3 == "") ? \$2 : \$2 "_" \$3
            split(\$2, parts, "_")
            print "raw_files/" \$1, sample, sample, parts[1], ((\$3 == "") ? 1 : \$3)
        }
    ' lookup.tsv > results/experiment_annotation.tsv

    cat <<-END_VERSIONS >| versions.yml
    "${task.process}":
        fragpipe: "24.0"
    END_VERSIONS
    """
}

/*
 * Bash shared by script: and stub: that writes lookup.tsv (file name, experiment,
 * bioreplicate, data type; one row per staged .mzML/.d file) and sets LOOKUP_SOURCE.
 * manifest_content is authoritative; file_experiment_map is the fallback.
 * Values are base64-encoded so '$', backticks and quotes in them never expand in bash.
 */
def fragpipeLookupScript(manifest_content, file_experiment_map) {
    def manifest_b64 = manifest_content.toString().getBytes('UTF-8').encodeBase64().toString()
    def map_b64      = file_experiment_map.toString().getBytes('UTF-8').encodeBase64().toString()
    return """
    printf '%s' '${manifest_b64}' | base64 -d > input_manifest.tsv
    printf '%s' '${map_b64}' | base64 -d > input_file_experiment_map.tsv

    if [ -s input_manifest.tsv ]; then
        LOOKUP_SOURCE=manifest_content
        awk -F'\\t' -v OFS='\\t' '
            NF == 0 { next }
            NF != 4 || \$1 == "" || \$2 == "" || \$4 == "" {
                print "ERROR: FRAGPIPE_HEADLESS: manifest_content line " NR " is not file<TAB>experiment<TAB>bioreplicate<TAB>data_type: " \$0 > "/dev/stderr"
                bad = 1
                next
            }
            { print }
            END { exit bad }
        ' input_manifest.tsv > lookup.tsv
    elif [ -s input_file_experiment_map.tsv ]; then
        LOOKUP_SOURCE=file_experiment_map
        awk -F'\\t' -v OFS='\\t' '
            NF == 0 { next }
            NF != 2 || \$1 == "" || \$2 == "" {
                print "ERROR: FRAGPIPE_HEADLESS: file_experiment_map line " NR " is not file<TAB>experiment: " \$0 > "/dev/stderr"
                bad = 1
                next
            }
            { print \$1, \$2, "1", "DDA" }
            END { exit bad }
        ' input_file_experiment_map.tsv > lookup.tsv
    else
        echo "ERROR: FRAGPIPE_HEADLESS: manifest_content and file_experiment_map are both empty" >&2
        exit 1
    fi

    for f in raw_files/*.mzML raw_files/*.d; do
        [ -e "\$f" ] || continue
        basename "\$f"
    done > staged_files.txt

    if [ ! -s staged_files.txt ] || [ ! -s lookup.tsv ]; then
        echo "ERROR: FRAGPIPE_HEADLESS: no staged .mzML/.d files or no \${LOOKUP_SOURCE} rows" >&2
        exit 1
    fi

    DUPLICATES=\$(cut -f1 lookup.tsv | sort | uniq -d)
    if [ -n "\${DUPLICATES}" ]; then
        echo "ERROR: FRAGPIPE_HEADLESS: file name(s) listed more than once in \${LOOKUP_SOURCE}: \${DUPLICATES}" >&2
        exit 1
    fi

    MISSING=\$(awk -F'\\t' 'NR == FNR { listed[\$1] = 1; next } !(\$0 in listed)' lookup.tsv staged_files.txt)
    if [ -n "\${MISSING}" ]; then
        echo "ERROR: FRAGPIPE_HEADLESS: no \${LOOKUP_SOURCE} row for staged file(s): \${MISSING}" >&2
        exit 1
    fi

    UNSTAGED=\$(awk -F'\\t' 'NR == FNR { staged[\$0] = 1; next } !(\$1 in staged) { print \$1 }' staged_files.txt lookup.tsv)
    if [ -n "\${UNSTAGED}" ]; then
        echo "ERROR: FRAGPIPE_HEADLESS: \${LOOKUP_SOURCE} row(s) for file(s) not staged in raw_files/: \${UNSTAGED}" >&2
        exit 1
    fi
    """
}
