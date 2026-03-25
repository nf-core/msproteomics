#!/usr/bin/env python3
"""
Prepare FragPipe workflow file with updated paths.

Updates the workflow file with:
- Database path
- Thread count
- Other runtime configurations
"""

import argparse
import os
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description='Prepare FragPipe workflow file'
    )
    parser.add_argument(
        '--workflow',
        required=True,
        help='Input workflow file'
    )
    parser.add_argument(
        '--fasta',
        required=True,
        help='FASTA database file'
    )
    parser.add_argument(
        '--output',
        default='fragpipe.workflow',
        help='Output workflow file (default: fragpipe.workflow)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=0,
        help='Number of threads (default: auto)'
    )
    parser.add_argument(
        '--enable_msbooster',
        action='store_true',
        default=True,
        help='Enable MSBooster (default: True)'
    )
    parser.add_argument(
        '--enable_ionquant',
        action='store_true',
        default=True,
        help='Enable IonQuant (default: True)'
    )
    return parser.parse_args()


def update_workflow(input_file, output_file, fasta_path, threads,
                    enable_msbooster=True, enable_ionquant=True):
    """Update workflow file with new paths and settings."""

    fasta_abs = Path(fasta_path).resolve()

    # Read the input workflow
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Update lines
    updated_lines = []
    for line in lines:
        # Update database path
        if line.startswith('database.db-path='):
            line = f'database.db-path={fasta_abs}\n'

        # Update thread count
        elif line.startswith('msfragger.num_threads=') and threads > 0:
            line = f'msfragger.num_threads={threads}\n'
        elif line.startswith('diann.threads=') and threads > 0:
            line = f'diann.threads={threads}\n'

        # Update MSBooster setting
        elif line.startswith('msbooster.run-msbooster='):
            value = 'true' if enable_msbooster else 'false'
            line = f'msbooster.run-msbooster={value}\n'

        # Update IonQuant setting
        elif line.startswith('ionquant.run-ionquant='):
            value = 'true' if enable_ionquant else 'false'
            line = f'ionquant.run-ionquant={value}\n'

        # Update paths for Docker container
        elif line.startswith('diann.exec-path='):
            line = 'diann.exec-path=/opt/tools/diann/diann-1.8.1\n'
        elif line.startswith('diann.exe='):
            line = 'diann.exe=/opt/tools/diann/diann-1.8.1\n'

        # Disable decoy generation (assume database already has decoys)
        elif line.startswith('database.generate-decoys='):
            line = 'database.generate-decoys=false\n'

        updated_lines.append(line)

    # Write the output workflow
    with open(output_file, 'w') as f:
        f.writelines(updated_lines)

    print(f"Created workflow file: {output_file}")
    print(f"  - Database: {fasta_abs}")
    print(f"  - Threads: {threads if threads > 0 else 'auto'}")
    print(f"  - MSBooster: {'enabled' if enable_msbooster else 'disabled'}")
    print(f"  - IonQuant: {'enabled' if enable_ionquant else 'disabled'}")


def main():
    args = parse_args()

    if not os.path.exists(args.workflow):
        print(f"Error: Workflow file not found: {args.workflow}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.fasta):
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    update_workflow(
        args.workflow,
        args.output,
        args.fasta,
        args.threads,
        args.enable_msbooster,
        args.enable_ionquant
    )


if __name__ == '__main__':
    main()
