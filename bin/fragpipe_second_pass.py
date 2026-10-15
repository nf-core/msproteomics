#!/usr/bin/env python3
"""FragPipe 24.0's second-pass workflow and manifest (Write sub mzML finalizer, FragpipeRun.java:589-700).

When tab-run.write_sub_mzml=true, FragPipe's finalizer writes, after every tool has run:
  <wd>/fragpipe-second-pass.workflow       from <wd>/fragpipe.workflow and the console text (:636-682)
  <wd>/fragpipe-files-second-pass.fp-manifest   from <wd>/fragpipe-files.fp-manifest (:684-697)

Java semantics replicated:
  - console text: lines ending with U+200B are dropped first (:599), each value is the FIRST match of its
    pattern (matcher.find(), :600-633), parsed with Float.parseFloat / Integer.parseInt and printed with
    "" + value (Float.toString / Integer.toString); \\d is ASCII as in java.util.regex
  - files are read with BufferedReader.readLine (\\n, \\r or \\r\\n end a line), every line is String.trim()med
    (chars <= U+0020) and written with "\\n"
  - manifest: String.split("\\t") drops trailing empty fields; the LC-MS path becomes
    wd.resolve(upToLastDot(fileName) + "_sub.mzML").toAbsolutePath().normalize()

Float.toString is the shortest decimal that reads back as the same float, in Java's layout
(decimal for 1e-3 <= |v| < 1e7, else d.dddE<n>), which is what JDK 17's FloatingDecimal prints for the
values MSFragger reports (fixed 6-decimal numbers); JDK-4511638-type corner cases are not reproduced.
"""

import argparse
import posixpath
import re
import struct
import sys
from decimal import Decimal

JAVA_TRIM = "".join(chr(c) for c in range(0x21))

# FragpipeRun.java:600-633 (pattern, kind)
NEW_VALUES = {
    "fragment_mass_tolerance": (r"New fragment_mass_tolerance = ([\d.]+) PPM", "float"),
    "use_topN_peaks": (r"New use_topN_peaks = (\d+)", "int"),
    "minimum_ratio": (r"New minimum_ratio = ([\d.]+)", "float"),
    "intensity_transform": (r"New intensity_transform = (\d+)", "int"),
    "remove_precursor_peak": (r"New remove_precursor_peak = (\d+)", "int"),
    "max_fragment_charge": (r"New max_fragment_charge = (\d+)", "int"),
    "precursor_true_tolerance": (r"New precursor_true_tolerance = ([\d.]+) PPM", "float"),
}


def to_float32(value: float) -> float:
    return struct.unpack("f", struct.pack("f", value))[0]


def java_float_to_string(text: str) -> str:
    """Float.toString(Float.parseFloat(text)) for a non-negative decimal without exponent."""
    v = to_float32(float(text))
    if v == 0.0:
        return "0.0"
    if v != v or v in (float("inf"), float("-inf")):
        raise ValueError(f"not a finite float: {text}")
    for digits in range(1, 10):
        candidate = f"{v:.{digits}g}"
        if to_float32(float(candidate)) == v:
            break
    d = Decimal(candidate).normalize()
    sign, digit_tuple, exponent = d.as_tuple()
    digit_str = "".join(map(str, digit_tuple))
    sci_exp = len(digit_str) - 1 + exponent  # value = d.ddd * 10^sci_exp
    if 1e-3 <= abs(v) < 1e7:
        if exponent >= 0:
            plain = digit_str + "0" * exponent + ".0"
        elif sci_exp >= 0:
            plain = digit_str[: sci_exp + 1] + "." + digit_str[sci_exp + 1 :]
        else:
            plain = "0." + "0" * (-sci_exp - 1) + digit_str
        return ("-" if sign else "") + plain
    mantissa = digit_str[0] + "." + (digit_str[1:] or "0")
    return ("-" if sign else "") + f"{mantissa}E{sci_exp}"


def java_lines(text: str) -> list:
    """BufferedReader.readLine over the whole text."""
    if text == "":
        return []
    lines = re.split(r"\r\n|\r|\n", text)
    if text.endswith(("\n", "\r")):
        lines.pop()
    return lines


def console_values(console_text: str) -> dict:
    text = re.sub("[^\n]+​\n", "", console_text)
    values = {}
    for key, (pattern, kind) in NEW_VALUES.items():
        m = re.search(pattern, text, re.ASCII)
        if m is None:
            values[key] = None
        elif kind == "float":
            values[key] = java_float_to_string(m.group(1))
        else:
            values[key] = str(int(m.group(1)))
    return values


def second_pass_workflow(workflow_text: str, console_text: str, databases_dir: str) -> str:
    """FragpipeRun.java:636-682; the first rule that applies rewrites the line."""
    v = console_values(console_text)
    rules = [
        ("tab-run.write_sub_mzml", True, "tab-run.write_sub_mzml=false"),
        ("percolator.run-percolator", True, "percolator.run-percolator=true"),
        ("peptide-prophet.run-peptide-prophet", True, "peptide-prophet.run-peptide-prophet=false"),
        ("phi-report.print-decoys", True, "phi-report.print-decoys=false"),
        ("msfragger.calibrate_mass", True, "msfragger.calibrate_mass=0"),
        ("msfragger.fragment_mass_tolerance", v["fragment_mass_tolerance"] is not None, f"msfragger.fragment_mass_tolerance={v['fragment_mass_tolerance']}"),
        ("msfragger.fragment_mass_units", v["fragment_mass_tolerance"] is not None, "msfragger.fragment_mass_units=1"),
        ("msfragger.use_topN_peaks", v["use_topN_peaks"] is not None, f"msfragger.use_topN_peaks={v['use_topN_peaks']}"),
        ("msfragger.minimum_ratio", v["minimum_ratio"] is not None, f"msfragger.minimum_ratio={v['minimum_ratio']}"),
        ("msfragger.intensity_transform", v["intensity_transform"] is not None, f"msfragger.intensity_transform={v['intensity_transform']}"),
        ("msfragger.remove_precursor_peak", v["remove_precursor_peak"] is not None, f"msfragger.remove_precursor_peak={v['remove_precursor_peak']}"),
        ("msfragger.max_fragment_charge", v["max_fragment_charge"] is not None, f"msfragger.max_fragment_charge={v['max_fragment_charge']}"),
        ("msfragger.precursor_true_tolerance", v["precursor_true_tolerance"] is not None, f"msfragger.precursor_true_tolerance={v['precursor_true_tolerance']}"),
        ("msfragger.precursor_true_units", v["precursor_true_tolerance"] is not None, "msfragger.precursor_true_units=1"),
        ("database.db-path", True, "database.db-path=" + databases_dir.replace("\\", "\\\\")),
        ("database.decoy-tag", True, "database.decoy-tag="),
        ("workdir=", True, "workdir="),
    ]
    out = []
    for raw in java_lines(workflow_text):
        line = raw.strip(JAVA_TRIM)
        for prefix, applies, replacement in rules:
            if applies and line.startswith(prefix):
                line = replacement
                break
        out.append(line + "\n")
    return "".join(out)


def resolve_normalize(workdir: str, name: str) -> str:
    """wd.resolve(name).toAbsolutePath().normalize().toString() (FragpipeRun.java:694) for an absolute local
    path, and for a URI work dir (s3://bucket/path, an S3 --outdir) the same on its path part: scheme and
    authority are kept, the path is resolved and normalised."""
    m = re.match(r"^([A-Za-z][A-Za-z0-9+.-]*://[^/]*)(/.*)?$", workdir)
    if m:
        return m.group(1) + posixpath.normpath(posixpath.join(m.group(2) or "/", name))
    if not workdir.startswith("/"):
        raise ValueError(f"--workdir must be an absolute path or a URI, got {workdir}")
    return posixpath.normpath(posixpath.join(workdir, name))


def second_pass_manifest(manifest_text: str, workdir: str) -> str:
    """FragpipeRun.java:684-697."""
    out = []
    for raw in java_lines(manifest_text):
        parts = raw.strip(JAVA_TRIM).split("\t")
        while len(parts) > 1 and parts[-1] == "":
            parts.pop()
        name = posixpath.basename(parts[0].rstrip("/"))
        stem = name[: name.rfind(".")] if "." in name else name
        parts[0] = resolve_normalize(workdir, stem + "_sub.mzML")
        out.append("\t".join(parts) + "\n")
    return "".join(out)


def read_text(path: str) -> str:
    with open(path, "rb") as fh:
        return fh.read().decode("utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workflow", required=True, help="the run's fragpipe.workflow")
    ap.add_argument("--manifest", required=True, help="the run's fragpipe-files.fp-manifest")
    ap.add_argument("--log", action="append", default=[], help="MSFragger console output, in FragPipe's order (repeatable)")
    ap.add_argument("--workdir", required=True, help="FragPipe's work dir (absolute path or URI, e.g. s3://bucket/out/fragpipe), where the _sub.mzML files are")
    ap.add_argument("--databases-dir", required=True, help="FragpipeLocations.getDirDatabases() of the FragPipe install")
    ap.add_argument("--outdir", default=".", help="where to write the two files")
    a = ap.parse_args(argv)
    console = "".join(read_text(p) for p in a.log)
    with open(posixpath.join(a.outdir, "fragpipe-second-pass.workflow"), "w", encoding="utf-8", newline="") as fh:
        fh.write(second_pass_workflow(read_text(a.workflow), console, a.databases_dir))
    with open(posixpath.join(a.outdir, "fragpipe-files-second-pass.fp-manifest"), "w", encoding="utf-8", newline="") as fh:
        fh.write(second_pass_manifest(read_text(a.manifest), a.workdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
