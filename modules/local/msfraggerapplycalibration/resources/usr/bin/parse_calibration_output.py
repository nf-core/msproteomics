#!/usr/bin/env python3
"""Turn an MSFragger --split1 calibration run into the params and spectra of a split search.

Replicates calibrate() of FragPipe 24.0 tools/msfragger_pep_split.py (lines 498-611):
the optimized values MSFragger prints are applied to the search params (calibrate_mass = 2 only),
check_spectral_files is set to 0, calibrate_mass and everything else are kept, and each input
file is replaced by its .mzBIN_calibrated (or copied when MSFragger wrote none).

Usage:
    parse_calibration_output.py \
        --log_file msfragger_calibrate.log \
        --params_file fragger.params \
        --output calibrated_fragger.params \
        --spectra_files sample1.mzML \
        --output_dir calibrated_spectra
"""

import argparse
import pathlib
import re
import shutil

# msfragger_pep_split.py:79-82
RE_CALIBRATE_MASS = re.compile(r"^calibrate_mass *= *([012])\b", re.MULTILINE)


def calibrate_mass_of(params_txt: str) -> int:
    """calibrate_mass of a params file; msfragger_pep_split.py fails without one."""
    match = RE_CALIBRATE_MASS.search(params_txt)
    if match is None:
        raise ValueError("The MSFragger params file has no 'calibrate_mass = 0|1|2' line")
    return int(match.group(1))


def _single(pattern: str, out: str):
    """`value, = re.findall(...)` of msfragger_pep_split.py: the value when there is exactly
    one match, else None (the script's `except ValueError`)."""
    found = re.compile(pattern).findall(out)
    return found[0] if len(found) == 1 else None


def calibrated_params(params_txt: str, out: str) -> str:
    """Search params after calibration (msfragger_pep_split.py:524-593); `out` is MSFragger's output."""
    calibrate_mass = calibrate_mass_of(params_txt)
    if calibrate_mass not in (1, 2):
        raise ValueError(f"calibrate_mass = {calibrate_mass}: msfragger_pep_split.py does not calibrate")
    params_txt_new = params_txt
    if calibrate_mass == 2:
        new_precursor_true_tol = _single("New precursor_true_tolerance = (.+) PPM", out)
        new_ms2_tol = re.compile("New fragment_mass_tolerance = (.+) PPM").findall(out)
        if len(new_ms2_tol) != 1:
            raise ValueError(
                f"Expected one 'New fragment_mass_tolerance = ... PPM' line in the MSFragger "
                f"calibration output, found {len(new_ms2_tol)} (msfragger_pep_split.py:527 fails too)"
            )
        new_ms2_tol = new_ms2_tol[0]
        new_use_topN_peak = _single("New use_topN_peaks = (.+)", out)
        new_minimum_ratio = _single("New minimum_ratio = (.+)", out)
        new_intensity_transform = _single("New intensity_transform = (.+)", out)
        new_remove_precursor_peak = _single("New remove_precursor_peaks = (.+)", out)

        def sub(pattern, repl, text):
            return re.compile(pattern, re.MULTILINE).sub(repl, text)

        if new_precursor_true_tol is not None:
            params_txt_new = sub(r"^precursor_true_tolerance\s*=\s*[0-9.]+",
                                 f"precursor_true_tolerance = {new_precursor_true_tol}", params_txt_new)
            params_txt_new = sub(r"^precursor_true_units\s*=\s*[01]", "precursor_true_units = 1", params_txt_new)
        params_txt_new = sub(r"^fragment_mass_tolerance\s*=\s*[0-9.]+",
                             f"fragment_mass_tolerance = {new_ms2_tol}", params_txt_new)
        params_txt_new = sub(r"^fragment_mass_units\s*=\s*[01]", "fragment_mass_units = 1", params_txt_new)
        if new_use_topN_peak is not None:
            params_txt_new = sub(r"^use_topN_peaks\s*=\s*[0-9]+", f"use_topN_peaks = {new_use_topN_peak}",
                                 params_txt_new)
        if new_minimum_ratio is not None:
            params_txt_new = sub(r"^minimum_ratio\s*=\s*[0-9.]+", f"minimum_ratio = {new_minimum_ratio}",
                                 params_txt_new)
        if new_intensity_transform is not None:
            params_txt_new = sub(r"^intensity_transform\s*=\s*[0-9]",
                                 f"intensity_transform = {new_intensity_transform}", params_txt_new)
        if new_remove_precursor_peak is not None:
            params_txt_new = sub(r"^remove_precursor_peak\s*=\s*[0-9]",
                                 f"remove_precursor_peak = {new_remove_precursor_peak}", params_txt_new)

    # disable check_spectral_files in split search (msfragger_pep_split.py:589-593)
    params_txt_new, n = re.compile(r"^check_spectral_files\s*=\s*[0-9]", re.MULTILINE).subn(
        "check_spectral_files = 0", params_txt_new)
    if n == 0:
        params_txt_new += "\ncheck_spectral_files = 0"
    return params_txt_new


def collect_calibrated_spectra(spectra_files: list, output_dir: pathlib.Path) -> None:
    """msfragger_pep_split.py:595-610: move <stem>.mzBIN_calibrated, or copy the input file."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for spec_file in spectra_files:
        spec_path = pathlib.Path(spec_file)
        calibrated = spec_path.with_suffix(".mzBIN_calibrated")
        if calibrated.exists():
            shutil.move(str(calibrated), str(output_dir / calibrated.name))
        elif spec_path.is_dir():
            shutil.copytree(str(spec_path), str(output_dir / spec_path.name))
        else:
            shutil.copy(str(spec_path), str(output_dir / spec_path.name))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--log_file", required=True, help="MSFragger --split1 output")
    parser.add_argument("--params_file", required=True, help="The file's MSFragger search params")
    parser.add_argument("--output", required=True, help="Calibrated params to write")
    parser.add_argument("--spectra_files", nargs="+", required=True, help="Spectra MSFragger calibrated")
    parser.add_argument("--output_dir", required=True, help="Directory for the calibrated spectra")
    args = parser.parse_args()

    params_txt = pathlib.Path(args.params_file).read_bytes().decode()
    out = pathlib.Path(args.log_file).read_bytes().decode()
    pathlib.Path(args.output).write_bytes(calibrated_params(params_txt, out).encode())
    collect_calibrated_spectra(args.spectra_files, pathlib.Path(args.output_dir))


if __name__ == "__main__":
    main()
