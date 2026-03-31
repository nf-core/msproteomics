#!/usr/bin/env python3
"""Parse MSFragger calibration output and generate calibrated params file.

Replicates calibrate() lines 498-611 from FragPipe's msfragger_pep_split.py.
Extracts optimized parameters from MSFragger --split1 calibration log, applies
them to the original fragger.params file, and collects calibrated spectra.

Usage:
    parse_calibration_output.py \
        --log_file msfragger_calibrate.log \
        --params_file fragger.params \
        --output calibrated_fragger.params \
        --spectra_files sample1.mzML sample2.mzML \
        --output_dir calibrated_spectra
"""

import argparse
import os
import pathlib
import re
import shutil
import sys


def parse_optimized_params(log_text: str) -> dict:
    """Extract optimized parameters from MSFragger calibration log.

    Matches the regex patterns from msfragger_pep_split.py lines 526-544.
    Note: the log prints 'remove_precursor_peaks' (plural) but the param
    name is 'remove_precursor_peak' (singular).

    Args:
        log_text: Full text of the MSFragger calibration stdout/stderr.

    Returns:
        Dictionary of parameter name to optimized value string.
    """
    patterns = {
        "fragment_mass_tolerance": r"New fragment_mass_tolerance = (.+) PPM",
        "precursor_true_tolerance": r"New precursor_true_tolerance = (.+) PPM",
        "use_topN_peaks": r"New use_topN_peaks = (.+)",
        "minimum_ratio": r"New minimum_ratio = (.+)",
        "intensity_transform": r"New intensity_transform = (.+)",
        "remove_precursor_peak": r"New remove_precursor_peaks = (.+)",
    }
    result = {}
    for param, pattern in patterns.items():
        matches = re.compile(pattern).findall(log_text)
        if matches:
            result[param] = matches[0].strip()
    return result


def apply_optimized_params(params_txt: str, optimized: dict) -> str:
    """Apply optimized parameters to the fragger.params text.

    Replicates msfragger_pep_split.py lines 566-587: for calibrate_mass=2,
    update tolerance values (setting units to PPM=1) and spectral processing
    parameters.

    Args:
        params_txt: Original fragger.params file content.
        optimized: Dictionary from parse_optimized_params().

    Returns:
        Updated params text with optimized values applied.
    """
    if "precursor_true_tolerance" in optimized:
        params_txt = re.sub(
            r"^precursor_true_tolerance\s*=\s*[0-9.]+",
            f"precursor_true_tolerance = {optimized['precursor_true_tolerance']}",
            params_txt,
            flags=re.MULTILINE,
        )
        params_txt = re.sub(
            r"^precursor_true_units\s*=\s*[01]",
            "precursor_true_units = 1",
            params_txt,
            flags=re.MULTILINE,
        )

    if "fragment_mass_tolerance" in optimized:
        params_txt = re.sub(
            r"^fragment_mass_tolerance\s*=\s*[0-9.]+",
            f"fragment_mass_tolerance = {optimized['fragment_mass_tolerance']}",
            params_txt,
            flags=re.MULTILINE,
        )
        params_txt = re.sub(
            r"^fragment_mass_units\s*=\s*[01]",
            "fragment_mass_units = 1",
            params_txt,
            flags=re.MULTILINE,
        )

    if "use_topN_peaks" in optimized:
        params_txt = re.sub(
            r"^use_topN_peaks\s*=\s*[0-9]+",
            f"use_topN_peaks = {optimized['use_topN_peaks']}",
            params_txt,
            flags=re.MULTILINE,
        )

    if "minimum_ratio" in optimized:
        params_txt = re.sub(
            r"^minimum_ratio\s*=\s*[0-9.]+",
            f"minimum_ratio = {optimized['minimum_ratio']}",
            params_txt,
            flags=re.MULTILINE,
        )

    if "intensity_transform" in optimized:
        params_txt = re.sub(
            r"^intensity_transform\s*=\s*[0-9]",
            f"intensity_transform = {optimized['intensity_transform']}",
            params_txt,
            flags=re.MULTILINE,
        )

    if "remove_precursor_peak" in optimized:
        params_txt = re.sub(
            r"^remove_precursor_peak\s*=\s*[0-9]",
            f"remove_precursor_peak = {optimized['remove_precursor_peak']}",
            params_txt,
            flags=re.MULTILINE,
        )

    return params_txt


def finalize_params(params_txt: str) -> str:
    """Apply final calibration settings to the params text.

    Replicates msfragger_pep_split.py lines 589-593 plus additional
    downstream requirements:
    - Set check_spectral_files = 0 (split search requirement)
    - Set calibrate_mass = 0 (calibration already done)
    - Remove database_name line (downstream modules set their own)

    Args:
        params_txt: Params text after optimized parameters are applied.

    Returns:
        Final params text ready for downstream search.
    """
    # Disable check_spectral_files for split search (msfragger_pep_split.py lines 589-593)
    params_txt, n = re.subn(
        r"^check_spectral_files\s*=\s*[0-9]",
        "check_spectral_files = 0",
        params_txt,
        flags=re.MULTILINE,
    )
    if n == 0:
        params_txt += "\ncheck_spectral_files = 0\n"

    # Set calibrate_mass = 0 (calibration is done, downstream search should not re-calibrate)
    params_txt = re.sub(
        r"^calibrate_mass\s*=\s*[0-9]",
        "calibrate_mass = 0",
        params_txt,
        flags=re.MULTILINE,
    )

    # Remove database_name line (downstream modules set their own database path)
    params_txt = re.sub(
        r"^database_name\s*=\s*.*\n?",
        "",
        params_txt,
        flags=re.MULTILINE,
    )

    return params_txt


def collect_calibrated_spectra(spectra_files: list, output_dir: str) -> None:
    """Collect calibrated spectra into the output directory.

    Replicates msfragger_pep_split.py lines 595-610: for each input spectra
    file, moves the .mzBIN_calibrated version if it exists, otherwise copies
    the original file/directory (e.g., Bruker .d directories).

    Args:
        spectra_files: List of original spectra file paths.
        output_dir: Directory to collect calibrated spectra into.
    """
    os.makedirs(output_dir, exist_ok=True)
    for spec_file in spectra_files:
        spec_path = pathlib.Path(spec_file)
        calibrated = spec_path.with_suffix(".mzBIN_calibrated")
        if calibrated.exists():
            dest = pathlib.Path(output_dir) / calibrated.name
            print(f"  Moving calibrated: {calibrated.name}")
            shutil.move(str(calibrated), str(dest))
        else:
            dest = pathlib.Path(output_dir) / spec_path.name
            print(f"  Copying original (not calibrated): {spec_path.name}")
            if spec_path.is_dir():
                shutil.copytree(str(spec_path), str(dest))
            else:
                shutil.copy(str(spec_path), str(dest))


def main():
    parser = argparse.ArgumentParser(
        description="Parse MSFragger calibration output and generate calibrated params."
    )
    parser.add_argument(
        "--log_file",
        required=True,
        help="Path to MSFragger calibration log file",
    )
    parser.add_argument(
        "--params_file",
        required=True,
        help="Path to original fragger.params file",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to write calibrated fragger.params",
    )
    parser.add_argument(
        "--spectra_files",
        nargs="+",
        required=True,
        help="List of original spectra files",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to collect calibrated spectra into",
    )
    args = parser.parse_args()

    # Read inputs
    log_path = pathlib.Path(args.log_file)
    params_path = pathlib.Path(args.params_file)

    if not log_path.exists():
        print(f"ERROR: Log file not found: {args.log_file}", file=sys.stderr)
        sys.exit(1)
    if not params_path.exists():
        print(f"ERROR: Params file not found: {args.params_file}", file=sys.stderr)
        sys.exit(1)

    log_text = log_path.read_text()
    params_txt = params_path.read_text()

    # Parse optimized parameters from calibration log
    optimized = parse_optimized_params(log_text)
    if optimized:
        print(f"Found optimized parameters: {', '.join(f'{k}={v}' for k, v in optimized.items())}")
    else:
        print(
            "WARNING: No optimized parameters found in calibration log. "
            "Only check_spectral_files and calibrate_mass will be set.",
            file=sys.stderr,
        )

    # Apply optimized params and finalize
    params_txt = apply_optimized_params(params_txt, optimized)
    params_txt = finalize_params(params_txt)
    pathlib.Path(args.output).write_text(params_txt)
    print(f"Calibrated params written to {args.output}")

    # Collect calibrated spectra
    print("Collecting calibrated spectra:")
    collect_calibrated_spectra(args.spectra_files, args.output_dir)
    print(f"Calibrated spectra collected in {args.output_dir}")


if __name__ == "__main__":
    main()
