#!/usr/bin/env python3
"""Write a FASTA file with its proteins sorted, as FragPipe 24.0 msfragger_pep_split.py sample_fasta() does
for the calibration search (lines 614-620, called with n = min(num_parts, 1) = 1, so every protein is kept).

Each protein is '>' + its bytes up to the next '\\n>' + '\\n', and the list is sorted bytewise.

Usage:
    sort_fasta.py <input.fasta> <output.fasta>
"""

import pathlib
import sys


def sorted_fasta(data: bytes) -> bytes:
    """msfragger_pep_split.py:615-620: proteins split on b'\\n>', each re-prefixed with '>' and suffixed '\\n'."""
    prots = [b">" + e + b"\n" for e in data[1:].split(b"\n>")]
    return b"".join(sorted(prots))


def main(argv: list) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    src, dst = (pathlib.Path(a) for a in argv)
    out = sorted_fasta(src.read_bytes())
    dst.write_bytes(out)
    print(f"Sorted {out.count(b'>')} proteins", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
