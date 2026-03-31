/*
 * MSFRAGGER_CALIBRATE: Mass calibration + parameter optimization for split-database search
 *
 * Replicates calibrate() from FragPipe's msfragger_pep_split.py.
 * Runs MSFragger with --split1 flag on ALL samples to perform mass calibration,
 * then extracts optimized parameters (fragment_mass_tolerance, etc.).
 *
 * This is an AGGREGATE process: it runs once with all samples collected together,
 * producing calibrated spectra (.mzBIN_calibrated) and optimized params for downstream search.
 *
 * Input:
 *   - fasta: Protein sequence database (FASTA)
 *   - mzml_files: ALL mass spectrometry data files collected (mzML format)
 *   - params_file: MSFragger params file (native fragger.params format)
 *   - msfragger_dir: Unzipped MSFragger tool directory (optional)
 *
 * Output:
 *   - calibrated_spectra: .mzBIN_calibrated files (or original copies if not calibrated)
 *   - params: Updated fragger.params with optimized tolerances + check_spectral_files=0 + calibrate_mass=0
 *   - versions: Software versions
 */
process MSFRAGGER_CALIBRATE {
    tag 'calibrate'
    label 'process_high'

    input:
    path fasta
    path mzml_files
    path params_file  // fragger.params file (native MSFragger format)
    path msfragger_dir // Unzipped MSFragger tool directory (optional). Pass [] when not using.

    output:
    path "calibrated_spectra/*"                , emit: calibrated_spectra
    path "calibrated_fragger.params"           , emit: params
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT"  , emit: license_agreement
    tuple val("${task.process}"), val('msfragger'), eval("cat .msfragger_version"), emit: versions_msfragger, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mem = task.ext.java_xmx ? task.ext.java_xmx as int : (task.memory ? task.memory.toGiga() : 8)
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    def agree_license = task.ext.agree_fragpipe_license_agreement ?: false
    def mzml_input = mzml_files instanceof List ? mzml_files.join(' ') : mzml_files

    """
    
    export HOME=\$(pwd)
    export JAVA_OPTS="-Xmx${mem}G"

    # License gate
    if [ "${agree_license}" != "true" ]; then
        echo "ERROR: You must agree to the FragPipe license agreement before using this tool." >&2
        echo "Set ext.agree_fragpipe_license_agreement = true in your Nextflow config." >&2
        exit 1
    fi
    echo "INFO: User has set agree_fragpipe_license_agreement = true. Proceeding with licensed tool." >&2
    echo "User agreed to FragPipe license agreement (agree_fragpipe_license_agreement=true)" > I_AGREE_FRAGPIPE_LICENSE_AGREEMENT

    TOOLS_DIR="${tools_dir}"

    # JAR + native lib discovery: prefer provided dir, fall back to tools_dir
    if [ -d "${msfragger_dir}" ] && [ "\$(ls -A ${msfragger_dir} 2>/dev/null)" ]; then
        MSFRAGGER_JAR=\$(find ${msfragger_dir} -maxdepth 1 -name "MSFragger*.jar" -not -name "*-sources*" -type f 2>/dev/null | head -1 || true)
        BRUKER_DIR="${msfragger_dir}/ext/bruker"
        THERMO_DIR="${msfragger_dir}/ext/thermo"
    else
        MSFRAGGER_JAR=\$(find "\$TOOLS_DIR" -name 'MSFragger*.jar' -not -name '*-sources*' -type f 2>/dev/null | head -1 || true)
        BRUKER_DIR="${tools_dir}/../ext/bruker"
        THERMO_DIR="${tools_dir}/../ext/thermo"
    fi
    NATIVE_FLAGS=""
    [ -d "\$BRUKER_DIR" ] && NATIVE_FLAGS="\$NATIVE_FLAGS -Dlibs.bruker.dir=\$BRUKER_DIR"
    [ -d "\$THERMO_DIR" ] && NATIVE_FLAGS="\$NATIVE_FLAGS -Dlibs.thermo.dir=\$THERMO_DIR"

    # Version capture
    _ms_ver=\$(java -jar "\$MSFRAGGER_JAR" --help 2>&1 | grep -oP 'MSFragger-\\K[\\d.]+' | head -1 || true)

    echo "\$_ms_ver" > .msfragger_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    # =========================================================================
    # Step 1: Sort FASTA proteins lexicographically
    # Replicates sample_fasta() from msfragger_pep_split.py with stride=1
    # (min(num_parts, 1) = 1, so all proteins are included, just sorted)
    # =========================================================================
    echo "Sorting FASTA proteins lexicographically..." >&2
    python3 -c "
import pathlib, sys
fasta = pathlib.Path('${fasta}')
prots = [b'>' + e + b'\\n' for e in fasta.read_bytes()[1:].split(b'\\n>')]
sorted_prots = sorted(prots)
with pathlib.Path('sorted_${fasta}').open('wb') as f:
    f.writelines(sorted_prots)
print(f'Sorted {len(sorted_prots)} proteins', file=sys.stderr)
"

    # =========================================================================
    # Step 2: Create calibration params file
    # Points to sorted FASTA, preserves original calibrate_mass value
    # =========================================================================
    if [ -s "${params_file}" ]; then
        cp ${params_file} calibrate_fragger.params
        sed -i '/^database_name/d' calibrate_fragger.params
        sed -i '/^num_threads/d' calibrate_fragger.params
        # Ensure file ends with newline before appending
        sed -i -e '\$a\\' calibrate_fragger.params
    else
        # No params file provided — create minimal with runtime params only
        touch calibrate_fragger.params
    fi
    echo "database_name = \$(pwd)/sorted_${fasta}" >> calibrate_fragger.params
    echo "num_threads = ${task.cpus}" >> calibrate_fragger.params

    # =========================================================================
    # Step 3: Run MSFragger calibration with --split1
    # Command format: java -jar MSFragger.jar <params_path> --split1 <file1> <file2> ...
    # params BEFORE --split1 flag (matches FragPipe calling convention)
    # =========================================================================
    echo "Running MSFragger calibration (--split1)..." >&2
    CMD="java \$JAVA_OPTS -Dfile.encoding=UTF-8 \$NATIVE_FLAGS -jar \\"\$MSFRAGGER_JAR\\" \$(pwd)/calibrate_fragger.params --split1 ${mzml_input} ${args}"
    printf '%s\\n' "\$CMD" > msfragger_calibrate.log
    eval "\$CMD" 2>&1 | tee -a msfragger_calibrate.log

    # =========================================================================
    # Step 4: Parse calibration output and collect calibrated spectra
    # Uses parse_calibration_output.py to extract optimized params from stdout
    # and collect mzBIN_calibrated files (or original copies)
    # =========================================================================
    echo "Parsing calibration output and collecting calibrated spectra..." >&2
    parse_calibration_output.py \\
        --log_file msfragger_calibrate.log \\
        --params_file calibrate_fragger.params \\
        --output calibrated_fragger.params \\
        --spectra_files ${mzml_input} \\
        --output_dir calibrated_spectra

    echo "Calibration complete." >&2
    """

    stub:
    def mzml_input = mzml_files instanceof List ? mzml_files.join(' ') : mzml_files

    """
    mkdir -p calibrated_spectra

    # Create mzBIN_calibrated files for each input spectra file
    for f in ${mzml_input}; do
        fname=\$(basename "\$f")
        stem="\${fname%.*}"
        touch "calibrated_spectra/\${stem}.mzBIN_calibrated"
    done

    # Copy params file as calibrated output
    cp ${params_file} calibrated_fragger.params

    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .msfragger_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
