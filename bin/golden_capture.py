#!/usr/bin/env python3
"""Capture a finished FragPipe run as a golden data set.

Layout written (nothing is overwritten; an existing target is an error)::

    <root>/<dataset>/<workflow>/
        manifest.json               step manifest (bin/fragpipe_log_to_manifest.py)
                                    plus a "golden" block: layout + sha256 per file
        fragpipe_log.redacted.txt   the FragPipe log with license material removed
        files/work/...              FragPipe work dir   ({WORKDIR})
        files/input/...             spectra dir         ({LCMS_DIR}, symlinks resolved)
        files/db/...                FASTA + pepindex    ({DB_DIR})

The FragPipe ``log_*.txt`` is not copied verbatim (it contains license
material); any file named like a license is refused, and every copied file is
scanned for license file names / ``--license`` arguments before the capture is
accepted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fragpipe_log_to_manifest as fpm  # noqa: E402

LAYOUT = {"{WORKDIR}": "files/work", "{LCMS_DIR}": "files/input", "{DB_DIR}": "files/db"}
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
# A license file or license argument. Plain "licens" is not used: protein
# descriptions such as "DNA replication licensing factor" contain it.
LICENSE_REF_RE = re.compile(rb"licen[cs]e\.dat|(?<![\w-])-+licen[cs]e|licen[cs]e[ _-]?(key|file|path)", re.IGNORECASE)
# Redaction placeholders (bin/fragpipe_log_to_manifest.py): a license flag whose value is a
# placeholder, or a bare placeholder token.
LICENSE_PLACEHOLDER_RE = re.compile(rb"(?<![\w-])-+licen[cs]e[\w.-]*(?:=| +)\{LICENSE(?:_PATH)?\}|\{LICENSE(?:_PATH)?\}", re.IGNORECASE)


class CaptureError(Exception):
    """The run cannot be captured as golden data."""


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


# Binary spectra/index/image formats: byte patterns in them are not text (a .raw can
# contain "license"-like byte runs). Everything else is scanned unless it has NUL bytes.
BINARY_SUFFIXES = {".raw", ".mzbin", ".mzbin_calibrated", ".pepindex", ".tdf", ".tdf_bin", ".bin", ".png", ".pdf", ".gz", ".zip", ".wiff", ".wiff2", ".scan", ".speclib", ".parquet"}
NUL_PROBE_BYTES = 8192


def _is_binary(p: Path) -> bool:
    """Binary by extension (also inside a Bruker .d dir) or a NUL byte in the first 8 KiB."""
    if p.suffix.lower() in BINARY_SUFFIXES or any(part.lower().endswith(".d") for part in p.parts[:-1]):
        return True
    with p.open("rb") as fh:
        return b"\0" in fh.read(NUL_PROBE_BYTES)


def _mentions_license(p: Path) -> bool:
    """True if a text file references a license file, flag value or notice.

    The agreed redaction placeholders (``--license {LICENSE}``, ``--license={LICENSE}``,
    ``{LICENSE_PATH}``) are removed first; anything left that matches is a hit.
    Scanned line by line, so a placeholder is never split across reads.
    """
    with p.open("rb") as fh:
        for line in fh:
            if LICENSE_REF_RE.search(LICENSE_PLACEHOLDER_RE.sub(b"", line)):
                return True
    return False


def check_root(root: Path) -> None:
    """The storage root must exist (or be creatable) and be writable. No fallback."""
    try:
        root.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        raise CaptureError(f"golden root {root} cannot be created: {e}") from e
    if not os.access(root, os.W_OK):
        raise CaptureError(f"golden root {root} is not writable")


def _copy_run(log: Path, roots: dict, manifest: dict, dest: Path, lcms_source: Path | None = None) -> None:
    entries = []
    lcms_names = {re.fullmatch(r"\{INPUT:(.+)\}", f["path"]).group(1) for f in manifest["lcms_files"]}
    for key, rel in LAYOUT.items():
        src_root = Path(roots[key])
        for src in sorted(src_root.rglob("*")):
            if src.is_dir():
                continue
            r = src.relative_to(src_root)
            if re.search(r"licen[cs]", src.name, re.IGNORECASE):
                raise CaptureError(f"refusing to capture a license-named file: {key}/{r}")
            if key == "{WORKDIR}" and len(r.parts) == 1 and re.fullmatch(r"log_.*\.txt", src.name):
                continue  # replaced by fragpipe_log.redacted.txt
            if fpm.EXCLUDED_DIRS & set(r.parts[:-1]):
                continue  # run environment (Philosopher state, JVM user home, "?")
            if not src.exists():
                raise CaptureError(f"dangling symlink: {src}")
            dst = dest / rel / r
            dst.parent.mkdir(parents=True, exist_ok=True)
            origin = src
            if key == "{LCMS_DIR}" and src.name in lcms_names and src.stat().st_size == 0:
                # placeholder file (e.g. a Docker bind-mount target): the real spectra live elsewhere
                if lcms_source is None:
                    raise CaptureError(f"LCMS input {src.name} is empty in the run dir; pass --lcms-source DIR with the real file")
                origin = lcms_source / src.name
                if not origin.is_file() or origin.stat().st_size == 0:
                    raise CaptureError(f"--lcms-source has no non-empty {src.name}")
            shutil.copyfile(origin, dst)  # follows symlinks: inputs are stored as real files
            scanned = not _is_binary(dst)
            if scanned and _mentions_license(dst):
                raise CaptureError(f"captured file references a license file or argument, aborting: {rel}/{r}")
            entry = {"path": f"{key}/{r}", "golden": f"{rel}/{r}", "size": dst.stat().st_size, "sha256": _sha256(dst)}
            if not scanned:
                entry["content_scan"] = "skipped: binary"
            if origin != src:
                entry["copied_from"] = str(origin)
            entries.append(entry)

    redacted = fpm.redact_text(log.read_text(encoding="utf-8"))
    if fpm.contains_license_material(redacted):
        raise CaptureError("redacted log still mentions a license")
    (dest / "fragpipe_log.redacted.txt").write_text(redacted, encoding="utf-8")
    manifest["golden"] = {"layout": LAYOUT, "files": entries, "log": "fragpipe_log.redacted.txt"}
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def capture(log: Path, root: Path, dataset: str, workflow: str, lcms_source: Path | None = None) -> Path:
    for n, v in (("dataset", dataset), ("workflow", workflow)):
        if not NAME_RE.match(v):
            raise CaptureError(f"{n} name {v!r} must match {NAME_RE.pattern}")
    check_root(root)
    target = root / dataset / workflow
    if target.exists():
        raise CaptureError(f"{target} already exists; remove it first to recapture")
    manifest = fpm.build_manifest(log)
    roots = manifest["roots"]
    work = Path(roots["{WORKDIR}"]).resolve()
    if log.resolve().parent != work:
        raise CaptureError(f"log {log} is not in the run's work dir {work}")

    partial = target.parent / (target.name + ".partial")
    if partial.exists():
        raise CaptureError(f"{partial} exists (interrupted capture?); remove it first")
    partial.mkdir(parents=True)
    done = False
    try:
        _copy_run(log, roots, manifest, partial, lcms_source)
        partial.rename(target)
        done = True
    finally:
        if not done:
            shutil.rmtree(partial)
    return target


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--log", type=Path, required=True, help="FragPipe log_<timestamp>.txt inside the run's work dir")
    ap.add_argument("--dataset", required=True, help="data set name, e.g. tmt_C02_D02")
    ap.add_argument("--workflow", required=True, help="workflow name, e.g. TMT-labelcheck-304_perfile")
    ap.add_argument("--root", type=Path, default=os.environ.get("GOLDEN_ROOT"), required="GOLDEN_ROOT" not in os.environ, help="golden storage root (default: $GOLDEN_ROOT)")
    ap.add_argument("--lcms-source", type=Path, help="dir with the real spectra files when the run dir holds empty placeholders")
    a = ap.parse_args(argv)
    try:
        t = capture(a.log, a.root, a.dataset, a.workflow, a.lcms_source)
    except CaptureError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(t)
    return 0


if __name__ == "__main__":
    sys.exit(main())
