#!/usr/bin/env python3
"""
Create FragPipe manifest file from mzML files.

FragPipe manifest format (tab-separated):
<path_to_file> <experiment> <bioreplicate> <data_type>
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create FragPipe manifest file from mzML files"
    )
    parser.add_argument(
        "--mzml_files", nargs="+", required=True, help="Input mzML files"
    )
    parser.add_argument(
        "--experiment",
        default="experiment",
        help="Experiment name (default: experiment)",
    )
    parser.add_argument(
        "--output",
        default="manifest.tsv",
        help="Output manifest file (default: manifest.tsv)",
    )
    parser.add_argument(
        "--data_type",
        default="DDA",
        choices=["DDA", "DIA", "DDA+", "GPF-DIA", "DIA-Quant", "DIA-Lib"],
        help="Data type (default: DDA)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.output, "w") as f:
        for i, mzml_file in enumerate(args.mzml_files, 1):
            # Get absolute path
            mzml_path = Path(mzml_file).resolve()

            if not mzml_path.exists():
                print(f"Warning: File not found: {mzml_path}", file=sys.stderr)
                continue

            # Format: path, experiment, bioreplicate, data_type
            experiment_name = (
                f"{args.experiment}_{i}"
                if len(args.mzml_files) > 1
                else args.experiment
            )
            bioreplicate = 1

            line = f"{mzml_path}\t{experiment_name}\t{bioreplicate}\t{args.data_type}\n"
            f.write(line)

    print(f"Created manifest file: {args.output}")
    print(f"  - {len(args.mzml_files)} files")
    print(f"  - Experiment: {args.experiment}")
    print(f"  - Data type: {args.data_type}")


if __name__ == "__main__":
    main()
