#!/usr/bin/env python3
"""Port of FragPipe 24.0 CmdSaintExpress.configure's input steps (CmdSaintExpress.java:56-81, 98-167).

Reads a REPRINT table (reprint.spc.tsv or reprint.int.tsv) and writes SAINTexpress's
inter.dat, bait.dat and prey.dat into the output directory, byte for byte as FragPipe does
(reprintTsvToSaintInput). When the table is not a readable regular file, FragPipe logs an error
and skips the mode (CmdSaintExpress.java:72-79):
    log.error("Input file " + inputPath.toAbsolutePath() + " is not accessible.");
The same message goes to stderr and to <output dir>/skipped_mode.txt, and the exit status is 0.

Usage: reprint_to_saint_input.py <reprint.<mode>.tsv> <output dir>
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path


def java_split_tab(line: str) -> list[str]:
    """String.split("\\t"): trailing empty strings are removed."""
    parts = line.split("\t")
    while parts and parts[-1] == "":
        parts.pop()
    return parts


def java_parse_double_is_zero(text: str) -> bool:
    """Double.parseDouble(e) != 0 is false (Java trims whitespace; NaN != 0 is true)."""
    return float(text.strip()) == 0


def convert(input_path: Path, out_path: Path) -> None:
    inter_dat: list[str] = []
    bait_dat: list[str] = []
    prey_dat: list[str] = []
    # BufferedReader.readLine(): a line ends at \n, \r or \r\n; no empty line after the final terminator
    text = input_path.read_text(encoding="utf-8")
    rows = re.split(r"\r\n|\r|\n", text)
    if rows and rows[-1] == "":
        rows.pop()
    lines = iter(rows)
    first_row = java_split_tab(next(lines))
    has_protlen = first_row[2].upper() == "PROTLEN"
    has_geneid = first_row[1].upper() == "GENEID"
    column_offset = 3 if has_protlen else 2 if has_geneid else 1
    ip = []
    for s in first_row[column_offset:]:
        if not (s.endswith("_INT") or s.endswith("_SPC")):
            raise AssertionError(s)
        ip.append(s[:-4])
    second_row = java_split_tab(next(lines))
    bait_names = [ip[i] if second_row[i + column_offset].upper() == "CONTROL" else second_row[i + column_offset] for i in range(len(ip))]
    for i in range(len(ip)):
        bait_dat.append(f"{ip[i]}\t{bait_names[i]}\t{'C' if bait_names[i].upper().startswith('CONTROL') else 'T'}")
    for line in lines:
        ls = java_split_tab(line)
        prey = ls[0]
        prey_dat.append(f"{ls[0]}\t{ls[2] if has_protlen else -1}\t{ls[1] if has_geneid else ls[0]}")
        for i in range(len(ip)):
            e = ls[i + column_offset]
            if not java_parse_double_is_zero(e):
                inter_dat.append(f"{ip[i]}\t{bait_names[i]}\t{prey}\t{e}")
    out_path.mkdir(parents=True, exist_ok=True)
    (out_path / "inter.dat").write_text("".join(e + "\n" for e in inter_dat) + "\n", encoding="utf-8")
    (out_path / "bait.dat").write_text("".join(e + "\n" for e in bait_dat), encoding="utf-8")
    (out_path / "prey.dat").write_text("".join(e + "\n" for e in prey_dat), encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    input_path, out_path = Path(argv[0]), Path(argv[1])
    if not (input_path.is_file() and os.access(input_path, os.R_OK)):
        message = f"Input file {input_path.absolute()} is not accessible."
        print(message, file=sys.stderr)
        out_path.mkdir(parents=True, exist_ok=True)
        (out_path / "skipped_mode.txt").write_text(message + "\n", encoding="utf-8")
        return 0
    convert(input_path, out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
