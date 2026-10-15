#!/usr/bin/env python3
"""Replicate FragPipe 24.0's database decoy check (FragpipeRun.checkDbConfig).

FragPipe counts FASTA entries whose first header field starts with the workflow
decoy tag and refuses to run (headless: log.error + abort) when the decoy share
is 0, 100 %, or outside 40-60 %. Databases of 1 GiB or more are not inspected.

Sources (FragPipe-GUI/src/main/java/org/nesvilab/):
  fragpipe/FragpipeRun.java:2703-2765   checkDbConfig (messages, thresholds)
  fragpipe/tabs/TabDatabase.java:79     databaseSizeLimit = 1 << 30
  fragpipe/tabs/TabDatabase.java:278-285 validateFasta (size gate, counting)
  utils/FastaUtils.java:44-85           getDecoysCnt / getProtsTotal / readFasta
"""

import argparse
import sys
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

DATABASE_SIZE_LIMIT = 1 << 30  # TabDatabase.java:79


def _java_trim(s: str) -> str:
    """String.trim(): strip leading/trailing chars <= U+0020."""
    start, end = 0, len(s)
    while start < end and s[start] <= " ":
        start += 1
    while end > start and s[end - 1] <= " ":
        end -= 1
    return s[start:end]


def first_header_fields(fasta: Path) -> list:
    """FastaUtils.readFasta(p).ordered.get(0): trimmed first '|' field per header."""
    fields = []
    # InputStreamReader(..., UTF_8) replaces malformed bytes with U+FFFD and never throws
    with open(fasta, "r", encoding="utf-8", errors="replace", newline=None) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.startswith(">"):
                continue
            # readFasta emits a first descriptor iff a '|' follows or text remains after '>'
            if line.find("|", 1) < 0 and len(line) <= 1:
                continue
            nxt = line.find("|", 1)
            fields.append(_java_trim(line[1 : nxt if nxt >= 0 else len(line)]))
    return fields


def java_percent(fraction: float) -> str:
    """DecimalFormat("##.#'%'").format(fraction * 100) with HALF_EVEN rounding."""
    q = Decimal(fraction * 100).quantize(Decimal("0.1"), rounding=ROUND_HALF_EVEN)
    text = format(q, "f")
    if text.endswith(".0"):
        text = text[:-2]
    if text.startswith("0."):  # '##' allows zero integer digits
        text = text[1:]
    return text + "%"


def check_db_config(fasta: Path, decoy_tag: str) -> str:
    """Return FragPipe's headless error message, or '' when the database passes."""
    if fasta.stat().st_size >= DATABASE_SIZE_LIMIT:
        return ""  # isBigDatabase: decoy share is not checked
    fields = first_header_fields(fasta)
    if not fields:
        return "Database looks to be empty"
    decoys = sum(1 for f in fields if f.startswith(decoy_tag))
    fraction = decoys / float(len(fields))
    if fraction <= 0:
        return "No decoys found in the FASTA file."
    if fraction >= 1:
        return "All FASTA entries seem to be decoys."
    if fraction < 0.4 or fraction > 0.6:
        return f"FASTA file contains {java_percent(fraction)} decoys."
    return ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--decoy-tag", required=True)
    args = parser.parse_args()
    message = check_db_config(args.fasta, args.decoy_tag)
    if message:
        print(f"ERROR: {message}", file=sys.stderr)
        return 1
    print(f"Decoy check passed for {args.fasta} (decoy tag '{args.decoy_tag}')")
    return 0


if __name__ == "__main__":
    sys.exit(main())
