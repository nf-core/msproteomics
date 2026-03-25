process SPECLIBGEN {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // SpecLibGen generates spectral libraries using EasyPQP/FragPipe-SpecLib.
    // Calls fragpipe-speclib convertpsm + library directly (bypasses gen_con_spec_lib.py
    // to avoid path resolution issues on Fusion filesystem).
    // psm_files: per-sample psm.tsv files from PHILOSOPHER_FILTER
    // peptide_files: per-sample peptide.tsv files from PHILOSOPHER_FILTER
    // mzml_files: mass spec files for spectrum extraction
    // config_cli: bash-sourceable KEY='value' lines (CONVERT_ARGS, FRAGMENT_TYPES,
    //             LIBRARY_ARGS, RT_CAL, IM_CAL, KEEP_INTERMEDIATE)
    tuple val(meta), path(psm_files, stageAs: 'psms/sample_*'), path(peptide_files, stageAs: 'peps/sample_*'), path(mzml_files, stageAs: 'spectra/*'), val(config_cli)
    path(fasta)  // protein database (shared resource)

    output:
    tuple val(meta), path("${prefix}")                          , emit: results_dir
    tuple val(meta), path("${prefix}/library.tsv")              , emit: library_tsv, optional: true
    tuple val(meta), path("${prefix}/library.speclib")          , emit: library_speclib, optional: true
    tuple val(meta), path("${prefix}/library.parquet")          , emit: library_parquet, optional: true
    tuple val("${task.process}"), val('speclibgen'), eval("cat .speclibgen_version"), emit: versions_speclibgen, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'
    def decoy_tag = task.ext.decoy_tag ?: 'rev_'

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    # Source algorithm config (sets CONVERT_ARGS, FRAGMENT_TYPES, LIBRARY_ARGS,
    # RT_CAL, IM_CAL, KEEP_INTERMEDIATE)
    cat > .speclibgen.env << 'SPECLIBGEN_CONFIG'
${config_cli}
SPECLIBGEN_CONFIG
    . .speclibgen.env

    # Build convert command: algorithm args + runtime args (decoy_prefix, enable_unannotated)
    CONVERT_EXTRA_ARGS="\$CONVERT_ARGS --fragment_types '\$FRAGMENT_TYPES' --decoy_prefix ${decoy_tag}"

    # Discover SpecLibGen directory (for iRT reference files)
    SPECLIBGEN_DIR="${tools_dir}/speclib"
    _fp_ver=\$(echo "${tools_dir}" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .speclibgen_version

    echo "\$_fp_ver" > .fragpipe_version

    # Discover fragpipe-speclib (Python entry point for EasyPQP)
    SPECLIB=\$(which fragpipe-speclib 2>/dev/null || find /usr/local/bin -name "fragpipe-speclib" -type f 2>/dev/null | head -1 || true)
    echo "Using fragpipe-speclib: \$SPECLIB"

    # Patch easypqp bug: np.array_split converts DataFrame chunks to numpy arrays,
    # breaking .itertuples() calls in the thread function. Fixed in FragPipe-SpecLib
    # commit e035d1df (2026-02-27) but not yet in the container image.
    # Replace np.array_split(psms, nthreads) with a DataFrame-preserving split.
    CONVERT_PY=\$(python3 -c "import easypqp.convert; print(easypqp.convert.__file__)")
    if grep -q 'np.array_split(psms, nthreads)' "\$CONVERT_PY" 2>/dev/null; then
        echo "Patching easypqp convert.py to fix np.array_split DataFrame bug..."
        EASYPQP_DIR=\$(pwd)/easypqp_patched/easypqp
        mkdir -p "\$EASYPQP_DIR"
        cp -r "\$(dirname \$CONVERT_PY)"/* "\$EASYPQP_DIR/"
        sed -i 's|exe.map(f, np.array_split(psms, nthreads))|exe.map(f, [psms.iloc[idx] for idx in np.array_split(range(len(psms)), nthreads)])|' "\$EASYPQP_DIR/convert.py"
        export PYTHONPATH="\$(pwd)/easypqp_patched:\${PYTHONPATH:-}"
    fi

    # Merge per-sample psm.tsv files into a single psm.tsv
    HEADER_DONE=false
    for f in psms/sample_*; do
        [ -f "\$f" ] || continue
        if [ "\$HEADER_DONE" = false ]; then
            cp "\$f" psm.tsv
            HEADER_DONE=true
        else
            tail -n +2 "\$f" >> psm.tsv
        fi
    done
    [ "\$HEADER_DONE" = false ] && { echo "ERROR: No psm.tsv files found in psms/" >&2; exit 1; }
    echo "Merged psm.tsv: \$(wc -l < psm.tsv) lines"

    # Merge per-sample peptide.tsv files into a single peptide.tsv
    HEADER_DONE=false
    for f in peps/sample_*; do
        [ -f "\$f" ] || continue
        if [ "\$HEADER_DONE" = false ]; then
            cp "\$f" peptide.tsv
            HEADER_DONE=true
        else
            tail -n +2 "\$f" >> peptide.tsv
        fi
    done
    if [ "\$HEADER_DONE" = false ]; then
        echo "WARNING: No peptide.tsv files found — library step may use PSM-only mode" >&2
    else
        echo "Merged peptide.tsv: \$(wc -l < peptide.tsv) lines"
    fi

    # Create filelist of spectra files
    find -L spectra -type f \\( -name "*.mzML" -o -name "*.mzXML" \\) | sort > filelist_speclibgen.txt
    echo "Spectra files: \$(wc -l < filelist_speclibgen.txt)"

    # --- Per-sample EasyPQP convertpsm ---
    cd ${prefix}

    CONVERT_COUNT=0
    for mzml in ../spectra/*.mzML ../spectra/*.mzXML; do
        [ -f "\$mzml" ] || continue
        # Get basename without calibration suffix
        BASENAME=\$(basename "\$mzml" | sed 's/\\.[^.]*\$//' | sed 's/_\\(un\\)\\?calibrated\$//')

        # Extract PSMs for this spectra file (Spectrum column starts with basename.)
        head -1 ../psm.tsv > "\${BASENAME}_temp-psm.tsv"
        grep "^\${BASENAME}\\." ../psm.tsv >> "\${BASENAME}_temp-psm.tsv" || true

        # Count PSMs (excluding header)
        PSM_COUNT=\$(tail -n +2 "\${BASENAME}_temp-psm.tsv" | wc -l)
        echo "Sample \${BASENAME}: \${PSM_COUNT} PSMs"

        if [ "\$PSM_COUNT" -gt 0 ]; then
            echo "Running convertpsm for \${BASENAME}..."
            eval python3 "\$SPECLIB" convertpsm \$CONVERT_EXTRA_ARGS --enable_unannotated \\
                --psm "\${BASENAME}_temp-psm.tsv" --spectra "\$mzml" \\
                --exclude-range -1.5,3.5 \\
                --psms "\${BASENAME}.psmpkl" --peaks "\${BASENAME}.peakpkl" \\
                || { echo "WARNING: convertpsm failed for \${BASENAME}" >&2; }
            CONVERT_COUNT=\$((CONVERT_COUNT + 1))
        fi
    done

    echo "Converted \${CONVERT_COUNT} samples"

    # Create filelist of all psmpkl and peakpkl files for library step
    ls -1 *.psmpkl *.peakpkl 2>/dev/null > filelist_easypqp_library.txt || true
    if [ ! -s filelist_easypqp_library.txt ]; then
        echo "ERROR: No psmpkl files were generated by convertpsm" >&2
        exit 1
    fi

    # --- EasyPQP library step ---
    RT_REF_ARG=""
    IM_REF_ARG=""

    # Handle RT calibration
    case "\$RT_CAL" in
        noiRT) ;;  # No reference, auto-select a run
        ciRT)
            if [ -f "\$SPECLIBGEN_DIR/hela_irtkit.tsv" ]; then
                cp "\$SPECLIBGEN_DIR/hela_irtkit.tsv" irt.tsv
                RT_REF_ARG="--rt_reference irt.tsv"
            fi
            ;;
        Pierce_iRT)
            if [ -f "\$SPECLIBGEN_DIR/Pierce_iRT.tsv" ]; then
                cp "\$SPECLIBGEN_DIR/Pierce_iRT.tsv" irt.tsv
                RT_REF_ARG="--rt_reference irt.tsv"
            fi
            ;;
        Biognosys_iRT)
            echo "Using Biognosys iRT for alignment"
            ;;
        *)
            # File path — user-provided RT reference
            if [ -f "\$RT_CAL" ]; then
                RT_REF_ARG="--rt_reference \$RT_CAL"
            fi
            ;;
    esac

    # Handle IM calibration
    case "\$IM_CAL" in
        noIM) ;;
        *)
            if [ -f "\$IM_CAL" ]; then
                IM_REF_ARG="--im_reference \$IM_CAL"
            fi
            ;;
    esac

    # Build peptide args (only if peptide.tsv was merged)
    PEPTIDE_ARG=""
    if [ -f ../peptide.tsv ]; then
        PEPTIDE_ARG="--peptidetsv ../peptide.tsv"
    fi

    echo "Running EasyPQP library..."
    LIBRARY_RC=0
    python3 "\$SPECLIB" library \\
        --psmtsv ../psm.tsv \$PEPTIDE_ARG \\
        --out easypqp_lib_openswath.tsv \\
        \$LIBRARY_ARGS \\
        \$RT_REF_ARG \$IM_REF_ARG \\
        filelist_easypqp_library.txt \\
        || LIBRARY_RC=\$?

    # Fallback: if noiRT failed, retry with ciRT (matching gen_con_spec_lib.py behavior)
    if [ "\$LIBRARY_RC" -ne 0 ] && [ -z "\$RT_REF_ARG" ]; then
        echo "Not enough peptides for alignment with automatic selection. Retrying with ciRT..."
        if [ -f "\$SPECLIBGEN_DIR/hela_irtkit.tsv" ]; then
            cp "\$SPECLIBGEN_DIR/hela_irtkit.tsv" irt.tsv
        fi
        python3 "\$SPECLIB" library \\
            --psmtsv ../psm.tsv \$PEPTIDE_ARG \\
            --out easypqp_lib_openswath.tsv \\
            \$LIBRARY_ARGS \\
            --rt_reference irt.tsv \$IM_REF_ARG \\
            filelist_easypqp_library.txt \\
            || { echo "ERROR: Library generation failed even with ciRT" >&2; exit 1; }
    elif [ "\$LIBRARY_RC" -ne 0 ]; then
        echo "ERROR: Library generation failed" >&2
        exit 1
    fi

    # --- Post-processing: convert library formats ---
    if [ -f easypqp_lib_openswath.tsv ]; then
        cp easypqp_lib_openswath.tsv library.tsv
        python3 "\$SPECLIB" export --in easypqp_lib_openswath.tsv --type Spectronaut --out library.speclib 2>/dev/null || true
        python3 "\$SPECLIB" export --in easypqp_lib_openswath.tsv --type Parquet --out library.parquet 2>/dev/null || true
    fi

    # Clean up intermediate files if configured
    if [ "\$KEEP_INTERMEDIATE" != "true" ]; then
        rm -f *.psmpkl *.peakpkl *_temp-psm.tsv *_run_peaks.tsv filelist_easypqp_library.txt 2>/dev/null || true
    fi

    cd ..

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/library.tsv
    touch ${prefix}/library.speclib
    touch ${prefix}/speclibgen.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .speclibgen_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
