#!/usr/bin/env python3
"""Portable copies of the files FRAGPIPEPLAN's FragPipe run wrote: input paths back to the originals.

FRAGPIPEPLAN plans on staged files, so FragPipe's run log, its saved fragpipe.workflow and its
fragpipe-files.fp-manifest name paths inside the task. This writes copies where:

- every staged input (LC-MS file, FASTA, TMT annotation) is its original path, from the exact map the
  pipeline gives (``--map``: ``<kind>\\t<path relative to the task dir>\\t<original path>`` per line,
  kind ``lcms`` or ``file``, s3:// kept);
- a file FragPipe would write in a staged LC-MS directory (``<stem>.pepXML``, ``diaumpire_se.log``, ...) is
  in the original directory of that directory's files; when they come from several, the entry whose stem
  starts the name decides, and a name that no input's stem starts is an error. A staged LC-MS directory named on its own
  is the original directory of its files (an error when they come from several);
- FragPipe's work dir (``fp_work``) is ``@WORKDIR@`` and the FragPipe tools dir is ``@TOOLS_DIR@``.

The run log (log_*.txt) has its licence material redacted first (fragpipe_log_to_manifest.redact_text,
the golden logs' redaction). The run fails if any task path or licence material is left in an output.

Usage: fragpipe_portable.py --task-dir DIR --tools-dir DIR --map map.tsv --outdir OUT FILE...
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path, PurePosixPath

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fragpipe_log_to_manifest as fpm  # noqa: E402

WORKDIR_TOKEN = "@WORKDIR@"
TOOLS_TOKEN = "@TOOLS_DIR@"


class PortableError(Exception):
    """A path cannot be mapped to its original without guessing."""


def stem(name: str) -> str:
    for ext in (".mzML", ".mzml", ".raw", ".d", ".mgf", ".mzXML"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def boundary(path: str) -> re.Pattern:
    """A path, matched only as a whole path or as a directory prefix ("~" ends the FASTA in the run log's
    "~~~~~~Sample of <fasta>~~~~~~~" header, FragpipeRun.java:438)."""
    return re.compile(re.escape(path) + r"(?=$|[/\s\"'<:;,=)\]~])")


class PathMap:

    def __init__(self, task_dir: str, tools_dir: str, rows: list[tuple[str, str, str]]):
        self.task_dir = task_dir.rstrip("/")
        self.tools_dir = tools_dir.rstrip("/")
        self.exact = {f"{self.task_dir}/{rel}": orig for _kind, rel, orig in rows}
        # staged LC-MS dirs: staged dir -> [(stem, original dir)]
        self.dirs: dict[str, list[tuple[str, str]]] = {}
        for kind, rel, orig in rows:
            if kind != "lcms":
                continue
            staged = f"{self.task_dir}/{rel}"
            d = str(PurePosixPath(staged).parent)
            o = orig.rsplit("/", 1)[0]
            self.dirs.setdefault(d, []).append((stem(PurePosixPath(staged).name), o))

    def derived(self, staged_dir: str, name: str) -> str:
        """Original location of a file FragPipe names in a staged dir: the dir its inputs come from, or, when
        they come from several, the dir of the input whose stem starts the name."""
        origins = {o for _st, o in self.dirs[staged_dir]}
        if len(origins) == 1:
            return f"{origins.pop()}/{name}"
        matches = sorted((s for s in self.dirs[staged_dir] if name.startswith(s[0])), key=lambda s: len(s[0]), reverse=True)
        if not matches:
            raise PortableError(f"{staged_dir}/{name}: no staged input's stem starts the name")
        dirs = {o for st, o in matches if st == matches[0][0]}
        if len(dirs) != 1:
            raise PortableError(f"{staged_dir}/{name}: inputs with the same stem come from {sorted(dirs)}")
        return f"{dirs.pop()}/{name}"

    def apply(self, text: str) -> str:
        for staged in sorted(self.exact, key=len, reverse=True):
            text = boundary(staged).sub(lambda _m, o=self.exact[staged]: o, text)
        for staged_dir in sorted(self.dirs, key=len, reverse=True):
            text = re.sub(re.escape(staged_dir) + r"/([^/\s\"'<:;,=)\]]+)", lambda m, d=staged_dir: self.derived(d, m.group(1)), text)
            origins = {o for _st, o in self.dirs[staged_dir]}
            if boundary(staged_dir).search(text):
                if len(origins) != 1:
                    raise PortableError(f"{staged_dir} is named on its own, but its files come from {sorted(origins)}")
                text = boundary(staged_dir).sub(lambda _m, o=origins.pop(): o, text)
        text = boundary(f"{self.task_dir}/fp_work").sub(WORKDIR_TOKEN, text)
        text = boundary(self.tools_dir).sub(TOOLS_TOKEN, text)
        return text

    def check(self, name: str, text: str) -> None:
        at = text.find(self.task_dir)
        if at >= 0:
            raise PortableError(f"{name}: a task path is left: ...{text[max(0, at - 60):at + 120]}...")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--task-dir", required=True)
    ap.add_argument("--tools-dir", required=True)
    ap.add_argument("--map", required=True, type=Path, help="<lcms|file>\\t<path relative to the task dir>\\t<original path> per line")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("files", nargs="+", type=Path)
    a = ap.parse_args(argv)
    rows = [tuple(line.split("\t", 2)) for line in a.map.read_text().splitlines() if line]
    # the task dir as given: FragPipe prints its paths with toAbsolutePath().normalize(), which follows no link
    paths = PathMap(a.task_dir, a.tools_dir, rows)
    a.outdir.mkdir(parents=True, exist_ok=True)
    try:
        for f in a.files:
            text = f.read_text()
            if f.name.startswith("log_"):
                text = fpm.redact_text(text)
            text = paths.apply(text)
            paths.check(f.name, text)
            if fpm.contains_license_material(text):
                raise PortableError(f"{f.name}: licence material is left")
            (a.outdir / f.name).write_text(text)
    except PortableError as ex:
        print(f"fragpipe_portable.py: {ex}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
