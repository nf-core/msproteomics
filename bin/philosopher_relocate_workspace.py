#!/usr/bin/env python3
"""Move a Philosopher workspace (.meta/) to the directory it now lives in.

Philosopher stores absolute paths in .meta/meta.bin (met.Data: Home, MetaFile, MetaDir, DB;
philosopher lib/met/met.go New()) and writes reports and the workspace through them
(lib/rep/rep.go: PSMReport(m.Home, ...); Serialize() to d.MetaFile). FragPipe runs Filter,
FreeQuant, LabelQuant and Report in one group directory; a Nextflow task gets a copy of that
directory elsewhere, so every string that starts with the old Home is rewritten to the new one.

The msgpack stream is copied element by element; only str elements change.

Usage: philosopher_relocate_workspace.py <group_dir>
"""

import struct
import sys
from pathlib import Path


def _str_header(n: int) -> bytes:
    if n < 32:
        return bytes([0xA0 | n])
    if n < 0x100:
        return bytes([0xD9, n])
    if n < 0x10000:
        return b"\xda" + struct.pack(">H", n)
    return b"\xdb" + struct.pack(">I", n)


# fixed-size types: first byte -> number of payload bytes
FIXED = {0xC0: 0, 0xC2: 0, 0xC3: 0, 0xCA: 4, 0xCB: 8, 0xCC: 1, 0xCD: 2, 0xCE: 4, 0xCF: 8,
         0xD0: 1, 0xD1: 2, 0xD2: 4, 0xD3: 8, 0xD4: 2, 0xD5: 3, 0xD6: 5, 0xD7: 9, 0xD8: 17}


class Relocator:
    def __init__(self, data: bytes, old: str, new: str):
        self.data, self.old, self.new = data, old.encode(), new.encode()
        self.out = bytearray()
        self.changed = 0

    def _take(self, pos: int, n: int) -> int:
        self.out += self.data[pos:pos + n]
        return pos + n

    def _len(self, pos: int, size: int) -> int:
        return int.from_bytes(self.data[pos:pos + size], "big")

    def element(self, pos: int) -> int:
        b = self.data[pos]
        if b <= 0x7F or b >= 0xE0:  # fixint
            return self._take(pos, 1)
        if 0x80 <= b <= 0x8F:  # fixmap
            return self._items(self._take(pos, 1), 2 * (b & 0x0F))
        if 0x90 <= b <= 0x9F:  # fixarray
            return self._items(self._take(pos, 1), b & 0x0F)
        if 0xA0 <= b <= 0xBF:  # fixstr
            return self._str(pos + 1, b & 0x1F)
        if b in FIXED:
            return self._take(pos, 1 + FIXED[b])
        if b in (0xD9, 0xDA, 0xDB):  # str8/16/32
            size = {0xD9: 1, 0xDA: 2, 0xDB: 4}[b]
            return self._str(pos + 1 + size, self._len(pos + 1, size))
        if b in (0xC4, 0xC5, 0xC6):  # bin8/16/32
            size = {0xC4: 1, 0xC5: 2, 0xC6: 4}[b]
            return self._take(pos, 1 + size + self._len(pos + 1, size))
        if b in (0xC7, 0xC8, 0xC9):  # ext8/16/32
            size = {0xC7: 1, 0xC8: 2, 0xC9: 4}[b]
            return self._take(pos, 2 + size + self._len(pos + 1, size))
        if b in (0xDC, 0xDD):  # array16/32
            size = 2 if b == 0xDC else 4
            n = self._len(pos + 1, size)
            return self._items(self._take(pos, 1 + size), n)
        if b in (0xDE, 0xDF):  # map16/32
            size = 2 if b == 0xDE else 4
            n = self._len(pos + 1, size)
            return self._items(self._take(pos, 1 + size), 2 * n)
        raise ValueError(f"unknown msgpack type byte 0x{b:02x} at {pos}")

    def _items(self, pos: int, n: int) -> int:
        for _ in range(n):
            pos = self.element(pos)
        return pos

    def _str(self, start: int, n: int) -> int:
        value = self.data[start:start + n]
        if value == self.old or value.startswith(self.old + b"/"):
            value = self.new + value[len(self.old):]
            self.changed += 1
        self.out += _str_header(len(value)) + value
        return start + n


def home_of(data: bytes) -> str:
    """The Home value of the met.Data map (first str value after the key "Home")."""
    key = _str_header(4) + b"Home"
    i = data.find(key)
    if i < 0:
        raise ValueError("meta.bin has no Home field")
    r = Relocator(data, "\0", "\0")
    pos = i + len(key)
    end = r.element(pos)
    value = bytes(r.out)
    hdr = value[0]
    skip = 1 if 0xA0 <= hdr <= 0xBF else {0xD9: 2, 0xDA: 3, 0xDB: 5}[hdr]
    assert end - pos == len(value)
    return value[skip:].decode()


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    group = Path(sys.argv[1]).resolve()
    meta = group / ".meta" / "meta.bin"
    data = meta.read_bytes()
    old = home_of(data)
    r = Relocator(data, old, str(group))
    if r.element(0) != len(data):
        raise ValueError(f"{meta}: trailing bytes after the workspace record")
    meta.write_bytes(bytes(r.out))
    print(f"Philosopher workspace {old} -> {group} ({r.changed} paths)", file=sys.stderr)


if __name__ == "__main__":
    main()
