#!/usr/bin/env python3
"""Compare MSFragger pepindex files after normalizing the path-dependent header.

Layout (measured on MSFragger 4.4.1): 4-byte header hash, 4-byte big-endian header length,
header text (params, including the absolute database_name), then the index body.
The hash and length change with the database path, so they are skipped; the header text is
compared with database_name replaced, and the body byte for byte.

Usage: compare_pepindex.py <expected.pepindex> <actual.pepindex> [...pairs]
"""
import hashlib
import re
import struct
import sys


def read(path):
    with open(path, "rb") as f:
        head = f.read(8)
        if len(head) != 8:
            raise ValueError(f"{path}: shorter than the 8-byte pepindex prefix")
        (length,) = struct.unpack(">i", head[4:8])
        header = f.read(length).decode("latin-1")
        if len(header) != length:
            raise ValueError(f"{path}: truncated header ({len(header)} of {length} bytes)")
        body = hashlib.sha256()
        for chunk in iter(lambda: f.read(1 << 24), b""):
            body.update(chunk)
    header = re.sub(r"(?m)^database_name = .*$", "database_name = <DB>", header)
    return header, body.hexdigest()


def main(argv):
    if len(argv) < 2 or len(argv) % 2:
        sys.exit(__doc__)
    ok = True
    for expected, actual in zip(argv[0::2], argv[1::2]):
        e, a = read(expected), read(actual)
        same = e == a
        ok &= same
        print(f"{'IDENTICAL' if same else 'DIFFERENT'} header={e[0] == a[0]} body={e[1] == a[1]} {expected} {actual}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
