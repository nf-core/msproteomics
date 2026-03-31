#!/usr/bin/env python3
"""
Generate SDRF-proteomics TSV from a simple CSV samplesheet and experiment-level parameters.

Takes a minimal CSV samplesheet (sample, spectra, condition) plus a JSON string of
experiment-level parameters and produces a standards-compliant SDRF file.

SDRF-Proteomics specification:
  https://github.com/bigbio/proteomics-sample-metadata/blob/master/sdrf-proteomics/README.adoc
Quick start guide:
  https://github.com/bigbio/proteomics-sample-metadata/blob/master/sdrf-proteomics/quickstart.adoc
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Ontology lookup tables
# ---------------------------------------------------------------------------

ENZYME_ONTOLOGY = {
    "Trypsin": "NT=Trypsin;AC=MS:1001251",
    "Trypsin/P": "NT=Trypsin/P;AC=MS:1001313",
    "Lys-C": "NT=Lys-C;AC=MS:1001309",
    "Lys-C/P": "NT=Lys-C/P;AC=MS:1001310",
    "Arg-C": "NT=Arg-C;AC=MS:1001303",
    "Asp-N": "NT=Asp-N;AC=MS:1001304",
    "Glu-C": "NT=Glu-C;AC=MS:1001917",
    "Chymotrypsin": "NT=Chymotrypsin;AC=MS:1001306",
    "CNBr": "NT=CNBr;AC=MS:1001307",
    "nonspecific": "NT=unspecific cleavage;AC=MS:1001956",
}

DISSOCIATION_ONTOLOGY = {
    "HCD": "NT=HCD;AC=MS:1000422",
    "CID": "NT=CID;AC=MS:1000133",
    "ETD": "NT=ETD;AC=MS:1000598",
    "EThcD": "NT=EThcD;AC=MS:1002631",
}

ACQUISITION_METHOD = {
    "DIA": "NT=Data-independent acquisition;AC=MS:1003215",
    "DDA": "NT=Data-dependent acquisition;AC=MS:1003214",
}

# SDRF column order
SDRF_COLUMNS = [
    "source name",
    "characteristics[organism]",
    "characteristics[organism part]",
    "characteristics[disease]",
    "characteristics[cell type]",
    "characteristics[biological replicate]",
    "assay name",
    "technology type",
    "comment[label]",
    "comment[data file]",
    "comment[file uri]",
    "comment[technical replicate]",
    "comment[fraction identifier]",
    "comment[cleavage agent details]",
    "comment[instrument]",
    "comment[modification parameters]",
    "comment[modification parameters]",
    "comment[proteomics data acquisition method]",
    "comment[dissociation method]",
    "comment[precursor mass tolerance]",
    "comment[fragment mass tolerance]",
]


def map_enzyme(enzyme_name):
    """Map enzyme name to MS ontology string."""
    if enzyme_name in ENZYME_ONTOLOGY:
        return ENZYME_ONTOLOGY[enzyme_name]
    return f"NT={enzyme_name};AC=MS:1001045"


def map_dissociation(method):
    """Map dissociation method to MS ontology string."""
    if method in DISSOCIATION_ONTOLOGY:
        return DISSOCIATION_ONTOLOGY[method]
    return f"NT={method};AC=MS:1000044"


def map_acquisition(method):
    """Map acquisition method to MS ontology string."""
    if method in ACQUISITION_METHOD:
        return ACQUISITION_METHOD[method]
    return f"NT={method}"


def format_modification(mod_string, mod_type):
    """
    Format modification string for SDRF.

    Args:
        mod_string: Comma-separated modification names or empty string.
        mod_type: 'fixed' or 'variable'.

    Returns:
        Formatted SDRF modification string, or empty string if no mods.
    """
    if not mod_string or mod_string.strip() == "":
        return ""
    mods = [m.strip() for m in mod_string.split(";") if m.strip()]
    if not mods:
        return ""
    return ";".join(f"NT={m};MT={mod_type}" for m in mods)


def assign_replicates(rows):
    """
    Auto-assign biological replicate numbers within each condition group.
    Numbers sequentially within each condition group.
    """
    counters = defaultdict(int)
    for row in rows:
        condition = row.get("condition", "")
        counters[condition] += 1
        row["_replicate"] = str(counters[condition])


def generate_sdrf(input_csv, params, output_path):
    """
    Generate SDRF file from CSV samplesheet and experiment parameters.

    Args:
        input_csv: Path to input CSV samplesheet.
        params: Dict of experiment-level parameters.
        output_path: Path to write output SDRF TSV.
    """
    # Read input CSV
    with open(input_csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        print("ERROR: Empty samplesheet", file=sys.stderr)
        sys.exit(1)

    # Validate required columns
    required = {"sample", "spectra"}
    header_set = set(rows[0].keys())
    missing = required - header_set
    if missing:
        print(f"ERROR: Missing required columns: {missing}", file=sys.stderr)
        sys.exit(1)

    # Auto-assign replicates if not provided
    assign_replicates(rows)

    # Extract experiment params
    organism = params.get("organism", "Homo sapiens")
    enzyme = params.get("enzyme", "Trypsin/P")
    fixed_mods = params.get("fixed_mods", "")
    variable_mods = params.get("variable_mods", "")
    instrument = params.get("instrument", "")
    acquisition = params.get("acquisition_method", "DDA")
    precursor_tol = params.get("precursor_mass_tolerance", "20 ppm")
    fragment_tol = params.get("fragment_mass_tolerance", "20 ppm")
    dissociation = params.get("dissociation_method", "HCD")

    # Format modifications
    fixed_mod_str = format_modification(fixed_mods, "fixed")
    variable_mod_str = format_modification(variable_mods, "variable")

    # Write SDRF
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(SDRF_COLUMNS)

        for row in rows:
            sample = row["sample"]
            filepath = row["spectra"]
            label = row.get("label", "")
            fraction = row.get("fraction", "1")
            replicate = row.get("_replicate", "1")

            # Determine label
            if label:
                label_str = f"NT={label}"
            else:
                label_str = "NT=label free sample;AC=MS:1002038"

            sdrf_row = [
                sample,                                         # source name
                organism,                                       # characteristics[organism]
                "not available",                                # characteristics[organism part]
                "not available",                                # characteristics[disease]
                "not applicable",                               # characteristics[cell type]
                replicate,                                      # characteristics[biological replicate]
                sample,                                         # assay name
                "proteomic profiling by mass spectrometry",     # technology type
                label_str,                                      # comment[label]
                os.path.basename(filepath),                     # comment[data file]
                filepath,                                       # comment[file uri]
                "1",                                            # comment[technical replicate]
                fraction,                                       # comment[fraction identifier]
                map_enzyme(enzyme),                             # comment[cleavage agent details]
                instrument or "NT=mass spectrometer;AC=MS:1000031",  # comment[instrument]
                fixed_mod_str,                                  # comment[modification parameters]
                variable_mod_str,                               # comment[modification parameters]
                map_acquisition(acquisition),                   # comment[proteomics data acquisition method]
                map_dissociation(dissociation),                 # comment[dissociation method]
                precursor_tol,                                  # comment[precursor mass tolerance]
                fragment_tol,                                   # comment[fragment mass tolerance]
            ]
            writer.writerow(sdrf_row)

    print(f"Generated SDRF with {len(rows)} samples: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate SDRF-proteomics TSV from CSV samplesheet and experiment parameters."
    )
    parser.add_argument(
        "--input", required=True, help="Path to input CSV samplesheet"
    )
    parser.add_argument(
        "--params", required=True, help="JSON string of experiment-level parameters"
    )
    parser.add_argument(
        "--output", required=True, help="Path to output SDRF TSV file"
    )
    args = parser.parse_args()

    params = json.loads(args.params)
    generate_sdrf(args.input, params, args.output)


if __name__ == "__main__":
    main()
