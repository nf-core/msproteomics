#!/usr/bin/env python3
"""Parse a FragPipe 24.0 headless log into a golden-test step manifest.

The log is the ``log_<timestamp>.txt`` FragPipe writes into its work directory
(or the console output of ``fragpipe --headless``). Its structure, from the
FragPipe 24.0 source (``FragpipeRun.java``, ``ProcessBuilderInfo.java``):

* ``LCMS files:`` block, one ``Experiment/Group: <name>`` per group followed by
  ``  - <path>\\t<data type>`` lines (``FragpipeRun.java:1160-1166``);
* ``~~~~~~~~~ fragpipe.config ~~~~~~~~~`` block with ``key=value`` properties;
* ``<N> commands to execute:`` followed by N pairs of lines: a header
  ``<name>[ [Work dir: <dir>]]`` and the command line
  (``FragpipeRun.java:412``, ``printProcessDescription`` at ``1170-1182``);
* the same header + command pairs again when each process runs, each closed by
  ``Process '<name>' finished, exit code: <n>`` (``ProcessBuilderInfo.java:202``);
* ``ALL JOBS DONE IN ...`` at the end of a complete run.

Command lines are ``String.join(" ", argv)`` with arguments containing
whitespace wrapped in double quotes (``FragpipeRun.java:1179-1181``), so the
argv is recovered by splitting on single spaces and re-joining quoted runs.

Anything that does not fit this structure raises ``LogFormatError``: the
manifest is the basis for coverage, so a silently skipped step is a bug.

License handling: the MSFragger/IonQuant license is confidential. Every
``--license <arg>`` pair is rewritten to ``--license {LICENSE}``, any
``license.dat`` path to ``{LICENSE_PATH}``, any other command token mentioning a
license to ``{LICENSE}``; command lines are never dropped, only non-command
informational lines mentioning a license. Logs already scrubbed this way are
accepted as is. The same ``redact_lines`` function is used for the manifest,
the captured log and log fixtures.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

MANIFEST_SCHEMA_VERSION = 1
LICENSE_ARG = "{LICENSE}"
LICENSE_PATH = "{LICENSE_PATH}"
LICENSE_RE = re.compile(r"licen[cs]", re.IGNORECASE)
# A license-named flag: any number of dashes, licence/license (any case), optional
# suffix (-key, file, ...), optional "=value".
LICENSE_FLAG_RE = re.compile(r"(?P<flag>-+licen[cs]e[\w.-]*)(?P<eq>=(?P<val>.*))?", re.IGNORECASE)
# A bare license flag name (no value), e.g. "--license" or "-licence-file": not confidential.
LICENSE_FLAG_NAME_RE = re.compile(r'(?<![\w-])-+licen[cs]e[\w.-]*(?=$|[\s"=,\]])', re.IGNORECASE)
LICENSE_PATH_RE = re.compile(r'[^ "\'=:]*licen[cs]e\.dat', re.IGNORECASE)
# First token of a command: a path (absolute, relative, home, quoted, drive letter).
EXEC_PATH_RE = re.compile(r'^("?(/|\./|\.\./|~/)|"?[A-Za-z]:\\)')
FLAG_TOKEN_RE = re.compile(r"-{1,3}[A-Za-z][\w.-]*(=.*)?")
INTERPRETERS = {"java", "python", "python3", "bash", "sh", "perl", "rscript", "mono", "dotnet", "env", "docker", "philosopher", "percolator", "diann"}

COMMANDS_HEADER_RE = re.compile(r"^(\d+) commands to execute:$")
STEP_HEADER_RE = re.compile(r"^(?P<name>.+?)(?: \[Work dir: (?P<wd>.+)\])?$")
FINISHED_RE = re.compile(r"^Process '(?P<name>.+)' finished, exit code: (?P<code>-?\d+)$")
ALL_DONE_RE = re.compile(r"ALL JOBS DONE IN ")
SECTION_RULE_RE = re.compile(r"^~{3,}$")
CONFIG_START = "~~~~~~~~~ fragpipe.config ~~~~~~~~~"
UUID_RE = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")
FRAGPIPE_JAR_RE = re.compile(r"(?P<root>/[^:\s\"]+)/lib/fragpipe-[0-9][0-9.]*\.jar")

# FragPipe process name (ProcessBuilderInfo.name) -> (Cmd class, tool).
# Names and classes from FragPipe 24.0 `cmd/Cmd*.java` NAME constants and the
# PbiBuilder.setName(...) calls (e.g. CmdPercolator.java:229).
KNOWN_STEPS: dict[str, tuple[str, str]] = {
    "CheckCentroid": ("CmdCheckCentroid", "CheckCentroid"),
    "WorkspaceCleanInit": ("CmdPhilosopherWorkspaceCleanInit", "PhilosopherWorkspace"),
    "WorkspaceClean": ("CmdPhilosopherWorkspaceClean", "PhilosopherWorkspace"),
    "MSFragger": ("CmdMsfragger", "MSFragger"),
    "MSFragger move pepxml": ("CmdMsfragger", "FileMove"),
    "MSFragger move pin": ("CmdMsfragger", "FileMove"),
    "MSFragger move tsv": ("CmdMsfragger", "FileMove"),
    "Percolator": ("CmdPercolator", "Percolator"),
    "Percolator: Convert to pepxml": ("CmdPercolator", "PercolatorOutputToPepXML"),
    "Percolator delete temp": ("CmdPercolator", "FileDelete"),
    "Rewrite pepxml": ("CmdPeptideProphet", "RewritePepxml"),
    "AppendFile": ("CmdAppendFile", "AppendFile"),
    "ProteinProphet": ("CmdProteinProphet", "ProteinProphet"),
    "PhilosopherDbAnnotate": ("CmdPhilosopherDbAnnotate", "PhilosopherDbAnnotate"),
    "PhilosopherFilter": ("CmdPhilosopherFilter", "PhilosopherFilter"),
    "PhilosopherReport": ("CmdPhilosopherReport", "PhilosopherReport"),
    "IonQuant": ("CmdIonquant", "IonQuant"),
    # Known FragPipe steps whose input/output conventions are not implemented yet:
    # their outputs are recorded as empty with an explicit flag.
    "MSBooster": ("CmdMSBooster", "MSBooster"),
    "Crystal-C": ("CmdCrystalc", "CrystalC"),
    "PeptideProphet": ("CmdPeptideProphet", "PeptideProphet"),
    "PTMProphet": ("CmdPtmProphet", "PTMProphet"),
    "PTMShepherd": ("CmdPtmshepherd", "PTMShepherd"),
    "TmtIntegrator": ("CmdTmtIntegrator", "TMTIntegrator"),
    "FreeQuant": ("CmdFreequant", "FreeQuant"),
    "SpecLibGen": ("CmdSpecLibGen", "SpecLibGen"),
    "diaTracer": ("CmdDiaTracer", "diaTracer"),
    "UmpireSe": ("CmdUmpireSe", "DIA-Umpire"),
    # CmdDiann.java: getCmdName() + " <step>" (lines 409, 487, 520, 546), FileDelete/FileMove for the .quant and speclib clean-up
    "DIA-Quant run DIA-NN": ("CmdDiann", "DIA-NN"),
    "DIA-Quant propagate information": ("CmdDiann", "DIA-NN propagate information"),
    "DIA-Quant convert DIA-NN output to MSstats.csv": ("CmdDiann", "DIA-NN MSstats"),
    "DIA-Quant generate site reports": ("CmdDiann", "DIA-NN site reports"),
    "DIA-Quant delete .quant file": ("CmdDiann", "FileDelete"),
    "DIA-Quant move and rename speclib for skyline": ("CmdDiann", "FileMove"),
    "SaintExpress": ("CmdSaintExpress", "SAINTexpress"),
    "PhilosopherAbacus": ("CmdPhilosopherAbacus", "PhilosopherAbacus"),
    "MBG": ("CmdMBGMatch", "MBG"),
    "OPair": ("CmdOPair", "OPair"),
    "PairScans": ("CmdPairScans", "PairScans"),
    # Remaining FragPipe 24.0 Cmd NAME constants (git grep 'String NAME = ' at tag 24.0, cmd/*.java)
    "Batched Run": ("CmdBatch", "BatchedRun"),
    "DbDownload": ("CmdDatabaseDownload", "DbDownload"),
    "DbUpdate": ("CmdDatabaseUpdate", "DbUpdate"),
    "ExportMatchedFragments": ("CmdExportMatchedFragments", "ExportMatchedFragments"),
    "FPOP": ("CmdFpopQuant", "FPOP"),
    "FPOP-adapter": ("CmdFPOPcoadaptr", "FPOP-adapter"),
    "MSFragger Digest": ("CmdMsfraggerDigest", "MSFraggerDigest"),
    "Metaproteomics": ("CmdMetaproteomics", "Metaproteomics"),
    "Quant (Isobaric)": ("CmdLabelquant", "Labelquant"),
    "Skyline": ("CmdSkyline", "Skyline"),
    "Transfer Learning": ("CmdTransferLearning", "TransferLearning"),
    "WriteSubMzml": ("CmdWriteSubMzml", "WriteSubMzml"),
    "iProphet": ("CmdIprophet", "iProphet"),
}
# Sub-steps named "<Cmd NAME> <action>" (CmdDiann.java:195-546: getCmdName() + " run DIA-NN", ...).
KNOWN_STEP_PREFIXES: dict[str, str] = {"DIA-Quant ": "CmdDiann"}


def step_kind(name: str) -> tuple[str, str] | None:
    """(Cmd class, tool) for a FragPipe process name, or None if unknown."""
    if name in KNOWN_STEPS:
        return KNOWN_STEPS[name]
    for prefix, cls in KNOWN_STEP_PREFIXES.items():
        if name.startswith(prefix):
            return cls, "DIA-Quant: " + name[len(prefix) :]
    # Sub-steps named "<Cmd NAME><sep><action>" (e.g. "MSFragger move pepxml"): longest NAME prefix.
    for known in sorted(KNOWN_STEPS, key=len, reverse=True):
        if name.startswith(known) and name[len(known) : len(known) + 1] in (" ", ":"):
            cls, tool = KNOWN_STEPS[known]
            return cls, f"{tool}: {name[len(known):].lstrip(' :')}"
    return None


LCMS_EXTS = (".mzml", ".mzxml", ".raw", ".d", ".mgf")  # compared lower-case: FragPipe accepts any case (.RAW, .mzML)

CONVENTIONS_IMPLEMENTED = {
    "CheckCentroid",
    "PhilosopherWorkspace",
    "MSFragger",
    "FileMove",
    "FileDelete",
    "Percolator",
    "PercolatorOutputToPepXML",
    "ProteinProphet",
    "PhilosopherDbAnnotate",
    "PhilosopherFilter",
    "PhilosopherReport",
    "IonQuant",
    "PeptideProphet",
    "RewritePepxml",
    "AppendFile",
    "PTMShepherd",
    "MBG",
    "PairScans",
    "OPair",
}

# Files FragPipe itself writes into the work dir (no logged command), with the
# FragPipe 24.0 source location that writes them.
FRAGPIPE_GENERATED = {
    "fragger.params": "CmdMsfragger.configure (params written per data type, CmdMsfragger.java:444-447)",
    "fragger_dda.params": "CmdMsfragger.java:444-447",
    "fragger_dia.params": "CmdMsfragger.java:444-447",
    "fragger_gpfdia.params": "CmdMsfragger.java:444-447",
    "fragger_dda_plus.params": "CmdMsfragger.java:444-447",
    "fragpipe-files.fp-manifest": "FragpipeRun.java:372",
    "fragpipe.workflow": "FragpipeRun.java:380",
    "fragpipe.job": "FragpipeRun.java:388",
    "filelist_proteinprophet.txt": "CmdProteinProphet.java:157",
    "filelist_ionquant.txt": "CmdIonquant.java:310",
    "modmasses_ionquant.txt": "CmdIonquant.java:347",
    "experiment_annotation.tsv": "ToolingUtils.java:223/301",
    "sdrf.tsv": "FragpipeRun.java:585 (written after all processes finish)",
    "shepherd.config": "CmdPtmshepherd.java:50 (CONFIG_FN, written at configure)",
}
# Directories that are run environment, not FragPipe results: Philosopher
# workspace state, the JVM user home FragPipe gives child processes, and the
# literal "?" dir created when user.home resolves to "?" in a container.
EXCLUDED_DIRS = {".meta", "fragpipe_child_jvm_user_home", "?"}
PERCOLATOR_OUTPUT_FLAGS = {
    "--results-psms",
    "--decoy-results-psms",
    "--results-peptides",
    "--decoy-results-peptides",
    "--results-proteins",
    "--decoy-results-proteins",
    "--weights",
}
PHILOSOPHER_REPORT_FILES = ("psm.tsv", "ion.tsv", "peptide.tsv", "protein.tsv", "protein.fas")
IONQUANT_INPLACE_FILES = ("psm.tsv", "ion.tsv", "peptide.tsv", "protein.tsv")
IONQUANT_OUTPUT_GLOBS = ("combined_*.tsv", "msstats.csv", "msstats_ptm.csv")


class LogFormatError(Exception):
    """The log does not have the FragPipe 24.0 headless structure."""


# ---------------------------------------------------------------- redaction


def _redact_token(tok: str) -> str:
    """Redact a non-flag token that mentions a license, keeping its non-license prefix."""
    if not contains_license_material(tok):
        return tok
    m = LICENSE_PATH_RE.search(tok)
    if m:  # e.g. -Dkey=/opt/fp/license.dat -> -Dkey={LICENSE_PATH}
        return LICENSE_PATH_RE.sub(LICENSE_PATH, tok)
    return LICENSE_ARG


def redact_command(line: str) -> str:
    """Redact license material in a command line token by token; argv shape is kept.

    ``-license X`` / ``--Licence X`` / ``--license-key X`` (any dashes, case,
    spelling): the flag is kept and the value token becomes ``{LICENSE}``.
    ``--license=X`` becomes ``--license={LICENSE}``. An empty value (two spaces)
    stays an empty argument and the following token is redacted too, so nothing
    is merged, dropped or leaked. The value is redacted even when it starts with
    ``-``: a license value cannot be told apart from an option, and leaking one
    is worse than masking an option of a malformed command. Any other token mentioning a license
    is replaced in place (``license.dat`` paths by ``{LICENSE_PATH}``).
    """
    toks = line.split(" ")
    out = list(toks)
    i = 0
    while i < len(toks):
        t = toks[i]
        m = LICENSE_FLAG_RE.fullmatch(t)
        if m:
            if m.group("eq") is not None:
                out[i] = f"{m.group('flag')}={LICENSE_ARG if m.group('val') not in (LICENSE_ARG, LICENSE_PATH) else m.group('val')}"
                i += 1
                continue
            j = i + 1
            while j < len(toks) and toks[j] == "":
                j += 1  # empty value argument(s): keep them, redact the next value too
            if j < len(toks):
                out[j] = toks[j] if toks[j] in (LICENSE_ARG, LICENSE_PATH) else LICENSE_ARG
            i = j + 1
            continue
        if t not in (LICENSE_ARG, LICENSE_PATH):
            out[i] = _redact_token(t)
        i += 1
    return " ".join(out)


def is_command_line(line: str, after_step_header: bool = False) -> bool:
    """A command line: the line after a step header, a first token that is a path or a
    known interpreter/script, or a plain word followed by ``-flag`` arguments."""
    if after_step_header:
        return True
    first = line.split(" ", 1)[0].strip('"')
    if not first:
        return False
    base = first.rsplit("/", 1)[-1].rsplit("\\", 1)[-1].lower()
    has_flag = any(FLAG_TOKEN_RE.fullmatch(t) for t in line.split(" ")[1:])
    return bool(
        EXEC_PATH_RE.match(line)
        or "/" in first
        or (has_flag and re.fullmatch(r"[\w.+-]+", first) is not None)
        or base in INTERPRETERS
        or re.search(r"\.(py|sh|jar|exe|pl|r)$", base)
        or re.match(r"(python|perl|java)[\d.]*$", base)
    )


def redact_line(line: str, after_step_header: bool = False) -> str | None:
    """Redact one log line. Command lines are never dropped (tokens are replaced);
    a non-command informational line that mentions a license is dropped (``None``).
    Already-scrubbed placeholders are kept as they are.
    """
    if not LICENSE_RE.search(line.replace(LICENSE_ARG, "").replace(LICENSE_PATH, "")):
        return line
    if is_command_line(line, after_step_header):
        red = redact_command(line)
        if contains_license_material(red):
            raise RuntimeError("license redaction left material in a command line")
        return red
    return None


def redact_lines(lines: list[str]) -> list[str]:
    """Redact a whole log; lines right after a known step header are commands."""
    out = []
    prev = None
    for ln in lines:
        hdr = prev is not None and _is_step_header(prev)
        r = redact_line(ln, after_step_header=hdr)
        if r is not None:
            out.append(r)
        prev = ln
    return out


def _is_step_header(line: str) -> bool:
    m = STEP_HEADER_RE.match(line)
    return bool(m) and step_kind(m.group("name")) is not None


def redact_text(text: str) -> str:
    return "\n".join(redact_lines(text.split("\n")))


def contains_license_material(text: str) -> bool:
    """True if text mentions a license other than placeholders and bare license flag names."""
    text = LICENSE_FLAG_NAME_RE.sub("", text)
    for allowed in (LICENSE_ARG, LICENSE_PATH):
        text = text.replace(allowed, "")
    return bool(LICENSE_RE.search(text))


# ------------------------------------------------------------------- argv


def split_command(line: str) -> list[str]:
    """Invert FragPipe's ``join(" ")`` with whitespace args wrapped in double quotes.

    An empty argument is not quoted (FragpipeRun.java:1180 quotes only args that
    match ``.*\\s.*``), so two consecutive spaces mean an empty argument and a
    leading/trailing space an empty first/last one; they are kept as ``""``.
    The wrap is only added around whitespace-containing arguments, so a quoted token without a space
    is an argument with literal quote characters (CmdPtmshepherd.java:163 passes ``"<config>"``
    through a ProcessBuilder, no shell) and keeps them.
    """
    out: list[str] = []
    buf: list[str] | None = None
    for tok in line.split(" "):
        if buf is not None:
            buf.append(tok)
            if tok.endswith('"'):
                out.append(" ".join(buf)[1:-1])
                buf = None
        elif tok.startswith('"') and not (len(tok) > 1 and tok.endswith('"')):
            buf = [tok]
        else:
            out.append(tok)
    if buf is not None:
        raise LogFormatError(f"unterminated quote in command line: {line!r}")
    return out


class Normalizer:
    """Replace absolute paths with placeholders, longest root first.

    ``inputs`` maps an exact external input file path to ``{INPUT:<basename>}``;
    ``roots`` maps a directory to a placeholder such as ``{WORKDIR}``.
    A root only matches at a path boundary (followed by ``/``, ``:``, ``"`` or end).
    """

    def __init__(self, roots: dict[str, str], inputs: dict[str, str] | None = None):
        self.inputs = dict(inputs or {})
        self.roots = sorted(roots.items(), key=lambda kv: len(kv[0]), reverse=True)
        self._input_re = (
            re.compile("|".join(re.escape(p) for p in sorted(self.inputs, key=len, reverse=True)) + r"(?=$|[:\"\s])")
            if self.inputs
            else None
        )

    def path(self, s: str) -> str:
        if self._input_re is not None:
            s = self._input_re.sub(lambda m: self.inputs[m.group(0)], s)
        for root, token in self.roots:
            s = re.sub(re.escape(root.rstrip("/")) + r"(?=$|[/:\"\s])", lambda _m, t=token: t, s)
        # Philosopher workspaces get a fresh java.io.tmpdir/<uuid> per run.
        return re.sub(r"\{TMP\}/" + UUID_RE.pattern, "{TMP}/{UUID}", s)

    def argv(self, argv: list[str]) -> list[str]:
        return [self.path(a) for a in argv]


# ------------------------------------------------------------------ parsing


def _parse_header(line: str) -> tuple[str, str | None]:
    m = STEP_HEADER_RE.match(line)
    if not m:
        raise LogFormatError(f"bad step header: {line!r}")
    return m.group("name"), m.group("wd")


def parse_log(text: str) -> dict:
    """Parse the log text into a raw structure (no filesystem access)."""
    lines = redact_lines(text.replace("\r\n", "\n").split("\n"))

    # LCMS files
    try:
        i = lines.index("LCMS files:")
    except ValueError as e:
        raise LogFormatError("no 'LCMS files:' section") from e
    lcms: list[dict] = []
    group = None
    i += 1
    while i < len(lines) and lines[i].strip():
        s = lines[i]
        if s.startswith("  Experiment/Group: "):
            group = s[len("  Experiment/Group: ") :]
        elif s.startswith("  - "):
            if group is None:
                raise LogFormatError(f"LCMS file before any Experiment/Group: {s!r}")
            path, sep, dtype = s[4:].rpartition("\t")
            if not sep:
                raise LogFormatError(f"LCMS line without data type: {s!r}")
            lcms.append({"path": path, "experiment": group, "data_type": dtype})
        elif not s.startswith("  (if "):
            raise LogFormatError(f"unexpected line in LCMS files section: {s!r}")
        i += 1
    if not lcms:
        raise LogFormatError("LCMS files section is empty")

    # fragpipe.config properties
    try:
        c = lines.index(CONFIG_START)
    except ValueError as e:
        raise LogFormatError("no fragpipe.config section") from e
    config: dict[str, str] = {}
    j = c + 1
    while j < len(lines) and not SECTION_RULE_RE.match(lines[j]):
        s = lines[j]
        if s and not s.startswith("#") and "=" in s:
            k, _, v = s.partition("=")
            config[k] = v
        j += 1
    if j >= len(lines):
        raise LogFormatError("fragpipe.config section is not terminated")
    exec_start = j + 1

    # planned commands
    plan_idx = next((k for k, s in enumerate(lines) if COMMANDS_HEADER_RE.match(s)), None)
    if plan_idx is None:
        raise LogFormatError("no '<N> commands to execute:' line")
    n_planned = int(COMMANDS_HEADER_RE.match(lines[plan_idx]).group(1))
    planned = []
    k = plan_idx + 1
    while not SECTION_RULE_RE.match(lines[k]):
        if k + 1 >= len(lines) or SECTION_RULE_RE.match(lines[k + 1]):
            raise LogFormatError(f"planned step header without command at log line {k + 1}")
        name, wd = _parse_header(lines[k])
        planned.append((name, wd, lines[k + 1]))
        k += 2
    if len(planned) != n_planned:
        raise LogFormatError(f"log announces {n_planned} commands but lists {len(planned)}")
    if plan_idx > exec_start:
        raise LogFormatError("command plan found after the execution section start")

    # Executed steps. FragPipe runs steps of a parallel group concurrently, so
    # their start blocks (header + command, printed atomically by the
    # synchronized printProcessDescription) and their "Process '<name>'
    # finished" markers interleave with each other and with tool output. A step
    # starts at the first not-yet-started planned step whose header and command
    # appear on two consecutive lines; a finish marker closes the earliest
    # started running step of that name. Tool output lines are not attributed.
    pending: dict[tuple[str, str], list[int]] = {}
    for idx, (name, wd, cmd) in enumerate(planned, start=1):
        header = name + (f" [Work dir: {wd}]" if wd is not None else "")
        pending.setdefault((header, cmd), []).append(idx)
    started: dict[int, int] = {}
    finished: dict[int, tuple[int, int]] = {}
    running: list[int] = []
    ambiguous: dict[int, list[int]] = {}
    pos = exec_start
    while pos < len(lines):
        line = lines[pos]
        if pos + 1 < len(lines) and pending.get((line, lines[pos + 1])):
            idx = pending[(line, lines[pos + 1])].pop(0)
            started[idx] = pos
            running.append(idx)
            pos += 2
            continue
        m = FINISHED_RE.match(line)
        if m:
            cands = [i for i in running if planned[i - 1][0] == m.group("name")]
            if not cands:
                raise LogFormatError(f"log line {pos + 1}: '{m.group('name')}' finished but no such step is running")
            idx = cands[0]
            if len(cands) > 1:
                for i in cands:
                    ambiguous.setdefault(i, sorted(set(cands) - {i}))
            running.remove(idx)
            finished[idx] = (pos, int(m.group("code")))
        elif ALL_DONE_RE.search(line):
            break
        pos += 1
    never = [i for i in range(1, len(planned) + 1) if i not in started]
    if never:
        raise LogFormatError(f"steps planned but never started: {never} ({planned[never[0] - 1][0]} first)")
    if running:
        raise LogFormatError(f"steps with no 'finished, exit code' line: {running} ({planned[running[0] - 1][0]} first)")

    steps = []
    for idx, (name, wd, cmd) in enumerate(planned, start=1):
        st, (fin, code) = started[idx], finished[idx]
        steps.append(
            {
                "step_index": idx,
                "name": name,
                "cwd_raw": wd,
                "command_line": cmd,
                "argv_raw": split_command(cmd),
                "exit_code": code,
                "log_line": st + 1,
                "finished_log_line": fin + 1,
            }
        )
    # Concurrency: steps whose [start, finish] log intervals overlap ran at the same time.
    for s1 in steps:
        s1["concurrent_with"] = [
            s2["step_index"]
            for s2 in steps
            if s2 is not s1 and s2["log_line"] < s1["finished_log_line"] and s1["log_line"] < s2["finished_log_line"]
        ]
        s1["finish_attribution_ambiguous_with"] = ambiguous.get(s1["step_index"], [])
    batch = 0
    seen: set[int] = set()
    by_idx = {x["step_index"]: x for x in steps}
    for s1 in steps:
        if s1["step_index"] in seen:
            continue
        batch += 1
        stack = [s1["step_index"]]
        while stack:
            i = stack.pop()
            if i in seen:
                continue
            seen.add(i)
            by_idx[i]["parallel_batch"] = batch
            stack += by_idx[i]["concurrent_with"]
    if not any(ALL_DONE_RE.search(s) for s in lines[pos:]):
        raise LogFormatError("no 'ALL JOBS DONE' line after the last step: run incomplete")
    for s in steps:
        if step_kind(s["name"]) is None:
            raise LogFormatError(
                f"unrecognized FragPipe step name {s['name']!r} (step {s['step_index']}); add it to KNOWN_STEPS"
            )
    return {"lcms_files": lcms, "config": config, "steps": steps}


# ---------------------------------------------------------- roots/listing


def discover_roots(parsed: dict) -> dict[str, str]:
    """Absolute directories the run used, keyed by placeholder."""
    cfg = parsed["config"]
    if "workdir" not in cfg:
        raise LogFormatError("fragpipe.config has no 'workdir' property")
    if "database.db-path" not in cfg:
        raise LogFormatError("fragpipe.config has no 'database.db-path' property")
    roots = {"{WORKDIR}": cfg["workdir"].rstrip("/")}
    lcms_dirs = {str(Path(f["path"]).parent) for f in parsed["lcms_files"]}
    if len(lcms_dirs) != 1:
        raise LogFormatError(f"LCMS files are in {len(lcms_dirs)} directories; only one is supported: {sorted(lcms_dirs)}")
    roots["{LCMS_DIR}"] = lcms_dirs.pop()
    roots["{DB_DIR}"] = str(Path(cfg["database.db-path"]).parent)
    fp_roots = {
        m.group("root") for s in parsed["steps"] for a in s["argv_raw"] for m in FRAGPIPE_JAR_RE.finditer(a)
    }
    if len(fp_roots) != 1:
        raise LogFormatError(f"expected one FragPipe install root (dir with lib/fragpipe-*.jar), found {sorted(fp_roots)}")
    roots["{FRAGPIPE}"] = fp_roots.pop()
    tmp = set()
    for s in parsed["steps"]:
        a = s["argv_raw"]
        for t in range(len(a) - 1):
            if a[t] == "--temp":
                p = Path(a[t + 1])
                if not UUID_RE.fullmatch(p.name):
                    raise LogFormatError(f"--temp argument does not end in a UUID: {a[t + 1]}")
                tmp.add(str(p.parent))
    if len(tmp) > 1:
        raise LogFormatError(f"several java.io.tmpdir roots: {sorted(tmp)}")
    if tmp:
        roots["{TMP}"] = tmp.pop()
    vals = list(roots.values())
    if len(set(vals)) != len(vals):
        raise LogFormatError(f"two placeholders share a directory: {roots}")
    return roots


def list_files(roots: dict[str, str]) -> dict[str, dict]:
    """Final state of every file under WORKDIR, LCMS_DIR and DB_DIR (symlinks reported, not followed)."""
    out = {}
    for key in ("{WORKDIR}", "{LCMS_DIR}", "{DB_DIR}"):
        base = Path(roots[key])
        if not base.is_dir():
            raise FileNotFoundError(f"{key} directory does not exist: {base}")
        for dp, dns, fns in os.walk(base):
            dns[:] = sorted(d for d in dns if d not in EXCLUDED_DIRS)
            for fn in sorted(fns):
                p = Path(dp) / fn
                st = p.lstat()
                out[str(p)] = {"size": st.st_size, "mtime_ns": st.st_mtime_ns, "symlink": p.is_symlink()}
    return out


# ----------------------------------------------------- per-tool conventions


def _abs(path: str, cwd: str | None) -> str:
    if os.path.isabs(path):
        return os.path.normpath(path)
    if cwd is None:
        raise LogFormatError(f"relative path {path!r} in a step without a logged work dir")
    return os.path.normpath(os.path.join(cwd, path))


def _read_params(path: str) -> dict[str, str]:
    params = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if "=" in line:
                k, _, v = line.partition("=")
                params[k.strip()] = v.strip()
    return params


def derive_io(step: dict, files: dict[str, dict], parsed: dict, roots: dict[str, str], state: dict) -> None:
    """Fill step['inputs'], step['outputs'], step['flags'] (absolute paths at this point)."""
    tool = step["tool"]
    a = step["argv_raw"]
    cwd = step["cwd_raw"]
    inputs: list[tuple[str, str]] = []
    outputs: list[tuple[str, str]] = []
    flags: list[str] = []

    def opt(flag: str) -> str | None:
        return a[a.index(flag) + 1] if flag in a and a.index(flag) + 1 < len(a) else None

    if tool not in CONVENTIONS_IMPLEMENTED:
        flags.append(f"io_convention_not_implemented:{tool}")
    elif tool == "CheckCentroid":
        inputs += [(x, "argv") for x in a if x.lower().endswith(LCMS_EXTS)]
    elif tool == "PhilosopherWorkspace":
        flags.append("philosopher_workspace_state:.meta (directory state, not a file output)")
    elif tool == "MSFragger":
        jar = next(t for t, x in enumerate(a) if x.endswith(".jar"))
        params_file = a[jar + 1]
        inputs.append((params_file, "argv"))
        params = _read_params(params_file)
        ranked = params.get("data_type", "0") != "0"
        top_n = int(params["output_report_topN"]) if ranked else 0
        # Every spectra file in argv: the manifest's LC-MS files and files an earlier step wrote (DIA-Umpire's
        # <name>_Q1.mzML pseudo-spectra, searched with the DDA params, and diaTracer's <name>_diatracer.mzML are not in the manifest)
        lcms = [x for x in a[jar + 2 :] if x.lower().endswith(LCMS_EXTS)]
        if not lcms:
            raise LogFormatError(f"MSFragger step {step['step_index']} has no LCMS file in argv")
        inputs += [(x, "argv") for x in lcms]
        fasta = params.get("database_name")
        if not fasta:
            raise LogFormatError(f"{params_file} has no database_name")
        inputs.append((fasta, "params:database_name"))
        fmt = params.get("output_format", "").lower()
        exts = {"pepxml": ".pepXML", "pin": ".pin", "tsv": ".tsv"}
        for x in lcms:
            stem = os.path.splitext(x)[0]
            for key, ext in exts.items():
                if key not in fmt.split("_"):
                    continue
                if key == "pepxml" and ranked:
                    # non-DDA data types write <stem>_rank1..N.pepXML (CmdMsfragger.outputs, N = output_report_topN)
                    outputs += [(f"{stem}_rank{k}{ext}", "convention:ranked pepXML") for k in range(1, top_n + 1)]
                else:
                    outputs.append((stem + ext, "convention:output_format"))
            for p in files:
                if p.startswith(stem + "_") and os.path.dirname(p) == os.path.dirname(x):
                    outputs.append((p, "convention:written_next_to_input"))
        if not state.get("pepindex_claimed"):
            idx = sorted(p for p in files if p.startswith(fasta + ".") and p.endswith(".pepindex"))
            outputs += [(p, "convention:pepindex") for p in idx]
            if idx:
                flags.append("pepindex_attributed_to_first_msfragger_step_inferred")
            state["pepindex_claimed"] = True
    elif tool == "FileMove":
        src, dst = a[-2], a[-1]
        inputs.append((src, "argv"))
        outputs.append((dst, "argv"))
        state.setdefault("moved", {})[os.path.normpath(src)] = step["step_index"]
    elif tool == "FileDelete":
        inputs.append((a[-1], "argv"))
        state.setdefault("deleted", {})[os.path.normpath(a[-1])] = step["step_index"]
    elif tool == "Percolator":
        inputs.append((_abs(a[-1], cwd), "argv"))
        for f in PERCOLATOR_OUTPUT_FLAGS:
            v = opt(f)
            if v is not None:
                outputs.append((_abs(v, cwd), "argv"))
    elif tool == "PercolatorOutputToPepXML":
        cls = next(t for t, x in enumerate(a) if x.endswith("PercolatorOutputToPepXML"))
        rest = a[cls + 1 :]
        if len(rest) != 8:
            raise LogFormatError(f"PercolatorOutputToPepXML expects 8 args, got {len(rest)}: {rest}")
        pin, basename, target, decoy, outprefix, dtype, _minprob, mzml = rest
        inputs += [(_abs(pin, cwd), "argv"), (_abs(target, cwd), "argv"), (_abs(decoy, cwd), "argv"), (mzml, "argv")]
        if dtype == "DDA":
            inputs.append((_abs(basename + ".pepXML", cwd), "convention:<basename>.pepXML"))
            outputs.append((_abs(outprefix + ".pep.xml", cwd), "convention:<outprefix>.pep.xml"))
        else:
            rank_re = re.compile(re.escape(basename) + r"_rank(\d+)\.pepXML")
            ranks = sorted(int(m.group(1)) for p in files if os.path.dirname(p) == cwd and (m := rank_re.fullmatch(os.path.basename(p))))
            if not ranks:
                raise LogFormatError(f"{dtype} PercolatorOutputToPepXML step {step['step_index']}: no {basename}_rank*.pepXML in {cwd}")
            for k in ranks:
                inputs.append((_abs(f"{basename}_rank{k}.pepXML", cwd), "convention:<basename>_rank<k>.pepXML"))
                outputs.append((_abs(f"{outprefix}_rank{k}.pep.xml", cwd), "convention:<outprefix>_rank<k>.pep.xml"))
    elif tool == "PeptideProphet":
        if opt("--database"):
            inputs.append((opt("--database"), "argv"))
        pepxmls = [x for x in a if x.endswith(".pepXML")]
        inputs += [(_abs(x, cwd), "argv") for x in pepxmls]
        if "--combine" in a:
            outputs.append((_abs("interact.pep.xml", cwd), "convention:--combine interact.pep.xml"))
        else:
            outputs += [(_abs("interact-" + os.path.basename(x)[: -len(".pepXML")] + ".pep.xml", cwd), "convention:interact-<stem>.pep.xml") for x in pepxmls]
        flags.append("peptideprophet_model_files_not_attributed")
    elif tool == "RewritePepxml":
        cls = next(t for t, x in enumerate(a) if x.endswith("RewritePepxml"))
        target, *spectra = a[cls + 1 :]
        inputs.append((target, "argv"))
        inputs += [(x, "argv") for x in spectra]
        outputs.append((target, "convention:in-place rewrite"))
    elif tool == "AppendFile":
        cls = next(t for t, x in enumerate(a) if x.endswith("AppendToFile"))
        target, source = a[cls + 1 : cls + 3]
        inputs += [(target, "argv"), (source, "argv")]
        outputs.append((target, "convention:in-place append"))
    elif tool == "PTMShepherd":
        inputs.append((a[-1], "argv"))
        cfg = _read_params(a[-1])
        with open(a[-1], encoding="utf-8") as fh:
            datasets = [ln.split("=", 1)[1].split() for ln in fh if ln.split("=", 1)[0].strip() == "dataset"]
        for ds in datasets:
            if len(ds) != 3:
                raise LogFormatError(f"unexpected PTM-Shepherd dataset line: {ds}")
            inputs.append((ds[1], "shepherd.config:dataset"))
            if cfg.get("run_glyco_mode", "false").lower() == "true":
                outputs.append((ds[1], "convention:glyco mode rewrites psm.tsv in place"))
                flags.append("ptmshepherd_glyco_psm_rewrite_inferred")
        outdir = os.path.join(cwd, "ptm-shepherd-output")
        outputs += [(p, "convention:ptm-shepherd-output/") for p in sorted(files) if p.startswith(outdir + "/")]
    elif tool == "MBG":
        for f in ("--manifest", "--residuedb", "--glycanmoddb"):
            if opt(f):
                inputs.append((os.path.normpath(_abs(opt(f), cwd)), "argv"))
        psm = _abs(opt("--psm"), cwd)
        inputs.append((psm, "argv"))
        outputs.append((psm, "convention:MBG --match rewrites --psm in place"))
        flags.append("mbg_psm_rewrite_inferred")
    elif tool == "PairScans":
        # CmdPairScans.java:120-126 (24.0): PairScans <spectra> <threads> <act1> <act2> <reverse> <single>;
        # writes <stem>.pairs next to the spectra file (observed in the O-Pair run).
        cls = next(t for t, x in enumerate(a) if x.endswith(".PairScans"))
        spectra = a[cls + 1]
        inputs.append((spectra, "argv"))
        outputs.append((os.path.splitext(spectra)[0] + ".pairs", "convention:<stem>.pairs next to spectra"))
    elif tool == "OPair":
        # CmdOPair.java (24.0): single-token "-x value" arguments; -r filelist (CmdOPair.java:102), -s psm.tsv, -o output dir.
        val = {t[1]: t[3:] for t in a if re.fullmatch(r"-[a-z] .*|-[a-z]", t) and len(t) >= 2}
        for k in ("r", "s", "f", "x", "y"):
            if val.get(k):
                inputs.append((os.path.normpath(val[k]), f"argv:-{k}"))
        if val.get("r"):
            with open(val["r"], encoding="utf-8") as fh:
                base = os.path.dirname(val["r"])
                inputs += [(_abs(x.strip(), base) if os.path.isabs(x.strip()) else os.path.join(roots["{LCMS_DIR}"], x.strip()), "filelist_opair") for x in fh if x.strip()]
        flags.append("opair_outputs_not_determined:-o " + (val.get("o") or "?"))
    elif tool == "ProteinProphet":
        filelist = _abs(a[-1], cwd)
        inputs.append((filelist, "argv"))
        with open(filelist, encoding="utf-8") as fh:
            inputs += [(_abs(x.strip(), cwd), "filelist") for x in fh if x.strip()]
        name = opt("--output")
        if name is None:
            raise LogFormatError("ProteinProphet without --output")
        outputs.append((_abs(name + ".prot.xml", cwd), "convention:--output <name>.prot.xml"))
    elif tool == "PhilosopherDbAnnotate":
        inputs.append((opt("--annotate"), "argv"))
        flags.append("philosopher_workspace_state:.meta (database annotation, not a file output)")
    elif tool == "PhilosopherFilter":
        pdir = _abs(opt("--pepxml"), cwd)
        inputs += [(p, "convention:--pepxml dir") for p in sorted(files) if os.path.dirname(p) == pdir and p.endswith(".pep.xml")]
        if opt("--protxml"):
            inputs.append((_abs(opt("--protxml"), cwd), "argv"))
        for f in ("--dbbin", "--probin"):
            if opt(f):
                flags.append(f"reads_philosopher_workspace_state:{f}")
        outputs.append((_abs("filter.log", cwd), "convention:stdout capture (CmdPhilosopherFilter.java:41)"))
        flags.append("philosopher_workspace_state:.meta (filter results, not a file output)")
    elif tool == "PhilosopherReport":
        outputs += [(_abs(f, cwd), "convention:philosopher report") for f in PHILOSOPHER_REPORT_FILES if _abs(f, cwd) in files]
        flags.append("reads_philosopher_workspace_state:.meta")
    elif tool == "IonQuant":
        fl = _abs(opt("--filelist"), cwd)
        inputs.append((fl, "argv"))
        if opt("--modlist"):
            inputs.append((_abs(opt("--modlist"), cwd), "argv"))
        with open(fl, encoding="utf-8") as fh:
            rows = [x.rstrip("\n").split("\t") for x in fh if x.strip()]
        if rows[0] != ["flag", "value"]:
            raise LogFormatError(f"unexpected IonQuant filelist header: {rows[0]}")
        for flag, value in rows[1:]:
            p = _abs(value, cwd)
            if flag == "--psm":
                inputs.append((p, "filelist"))
                for f in IONQUANT_INPLACE_FILES:
                    q = os.path.join(os.path.dirname(p), f)
                    if q in files:
                        outputs.append((q, "convention:in-place rewrite"))
            elif flag == "--specdir":
                inputs += [(q, "filelist:--specdir") for q in sorted(files) if os.path.dirname(q) == p and q.endswith((".mzML", ".mzBIN"))]
            else:
                flags.append(f"ionquant_filelist_flag_not_interpreted:{flag}")
        flags.append("ionquant_rewrites_per_experiment_reports_in_place_inferred")
        outputs += [(p, "convention:ionquant outputs") for p in sorted(files) if os.path.dirname(p) == cwd and any(fnmatch.fnmatch(os.path.basename(p), g) for g in IONQUANT_OUTPUT_GLOBS)]
    step["_inputs"] = inputs
    step["_outputs"] = outputs
    step["flags"] = flags


# ------------------------------------------------------------------ build


def build_manifest(log_path: Path) -> dict:
    text = log_path.read_text(encoding="utf-8", errors="strict")
    try:
        parsed = parse_log(text)
        roots = discover_roots(parsed)
    except LogFormatError as e:
        raise LogFormatError(f"{log_path}: {e}") from e
    files = list_files(roots)
    external = {f["path"]: "{INPUT:%s}" % Path(f["path"]).name for f in parsed["lcms_files"]}
    fasta = parsed["config"]["database.db-path"]
    external[fasta] = "{INPUT:%s}" % Path(fasta).name
    norm = Normalizer({v: k for k, v in roots.items()}, external)

    state: dict = {}
    for s in parsed["steps"]:
        s["cmd_class"], s["tool"] = step_kind(s["name"])
        derive_io(s, files, parsed, roots, state)

    workdir = roots["{WORKDIR}"]
    generated = {os.path.join(workdir, k): v for k, v in FRAGPIPE_GENERATED.items()}
    # written per experiment/group dir
    generated.update({p: "CmdOPair.java:102 (group dir)" for p in files if os.path.basename(p) == "filelist_opair.txt"})
    producer: dict[str, int] = {}
    claimed: set[str] = set()
    for s in parsed["steps"]:
        ins = []
        for p, via in s["_inputs"]:
            p = os.path.normpath(p)
            if p in producer:
                origin = f"step:{producer[p]}"
            elif p in generated:
                origin = "fragpipe_generated"
            elif p in external:
                origin = "external"
            elif p.startswith(roots["{FRAGPIPE}"] + "/"):
                origin = "fragpipe_install"
            else:
                origin = "unknown"
                s["flags"].append(f"input_origin_unknown:{norm.path(p)}")
            ins.append({"path": norm.path(p), "via": via, "origin": origin, "_abs": p})
            claimed.add(p)
        outs = []
        for p, via in s["_outputs"]:
            p = os.path.normpath(p)
            outs.append({"path": norm.path(p), "via": via, "_abs": p})
            producer[p] = s["step_index"]
            claimed.add(p)
        s["inputs"], s["outputs"] = ins, outs

    # final state of each output
    by_index = {s["step_index"]: s for s in parsed["steps"]}
    last_writer: dict[str, int] = {}
    for s in parsed["steps"]:
        for o in s["outputs"]:
            prev = last_writer.get(o["_abs"])
            if prev is not None and prev != s["step_index"]:
                for po in by_index[prev]["outputs"]:
                    if po["_abs"] == o["_abs"]:
                        po["state"] = f"overwritten_by_step:{s['step_index']}"
            last_writer[o["_abs"]] = s["step_index"]
    for s in parsed["steps"]:
        for o in s["outputs"]:
            p = o.pop("_abs")
            if "state" in o:
                continue
            if p in state.get("deleted", {}):
                o["state"] = f"deleted_by_step:{state['deleted'][p]}"
            elif p in state.get("moved", {}) and state["moved"][p] > s["step_index"]:
                o["state"] = f"moved_by_step:{state['moved'][p]}"
            elif p in files:
                o["state"] = "present"
                o["size"] = files[p]["size"]
            else:
                o["state"] = "missing"
                s["flags"].append(f"output_missing:{norm.path(p)}")

    # Input availability. A golden test stages FragPipe's input files from the
    # captured run, which only holds each file's final version. An input that
    # this step or a later step rewrites in place, moves or deletes is not
    # available from this run: its golden test must take inputs from a run
    # that stops before that step (e.g. a "report" run before PTM-Shepherd/IonQuant).
    for s in parsed["steps"]:
        destroyed = []
        for i in s["inputs"]:
            p = i.pop("_abs")
            later = sorted(
                t["step_index"] for t in parsed["steps"]
                if t["step_index"] >= s["step_index"] and any(o["path"] == i["path"] for o in t["outputs"])
            )
            if i["origin"] == "fragpipe_install":
                i["golden_state"] = "fragpipe_install"
            elif later:
                i["golden_state"] = f"rewritten_by_step:{later[-1]}"
                if later[0] == s["step_index"]:
                    destroyed.append(i["path"])
            elif p in state.get("deleted", {}) and state["deleted"][p] >= s["step_index"]:
                i["golden_state"] = f"deleted_by_step:{state['deleted'][p]}"
            elif p in state.get("moved", {}) and state["moved"][p] >= s["step_index"]:
                i["golden_state"] = f"moved_by_step:{state['moved'][p]}"
            elif p in files:
                i["golden_state"] = "available"
            else:
                i["golden_state"] = "missing"
        if destroyed:
            s["flags"].append("rewrites_inputs_in_place:" + ",".join(destroyed))
        s["golden_inputs_available"] = all(i["golden_state"] in ("available", "fragpipe_install") for i in s["inputs"])
        s["flags"] = list(dict.fromkeys(s["flags"]))

    unattributed = []
    for p in sorted(files):
        if p in claimed or p in generated:
            continue
        name = os.path.basename(p)
        if name.startswith("log_") and name.endswith(".txt") and os.path.dirname(p) == workdir:
            continue  # the FragPipe log itself
        unattributed.append({"path": norm.path(p), "size": files[p]["size"], "flag": "no step or FragPipe writer identified"})

    manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "generator": "bin/fragpipe_log_to_manifest.py",
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source_log": log_path.name,
        "fragpipe_version": _version(text),
        "roots": roots,
        "lcms_files": [{**f, "path": norm.path(f["path"])} for f in parsed["lcms_files"]],
        "database": norm.path(fasta),
        "fragpipe_generated": [
            {"path": norm.path(p), "writer": w, "present": p in files} for p, w in sorted(generated.items())
        ],
        "steps": [
            {
                "step_index": s["step_index"],
                "name": s["name"],
                "cmd_name": s["tool"],
                "cmd_class": s["cmd_class"],
                "cwd": norm.path(s["cwd_raw"]) if s["cwd_raw"] is not None else None,
                "argv": norm.argv(s["argv_raw"]),
                "argv_raw": s["argv_raw"],
                "exit_code": s["exit_code"],
                "parallel_batch": s["parallel_batch"],
                "concurrent_with": s["concurrent_with"],
                "golden_inputs_available": s["golden_inputs_available"],
                "inputs": s["inputs"],
                "outputs": s["outputs"],
                "flags": (["cwd_not_logged"] if s["cwd_raw"] is None else [])
                + (
                    [f"finish_marker_ambiguous_with_steps:{s['finish_attribution_ambiguous_with']}"]
                    if s["finish_attribution_ambiguous_with"]
                    else []
                )
                + s["flags"],
            }
            for s in parsed["steps"]
        ],
        "unattributed_files": unattributed,
    }
    if contains_license_material(json.dumps(manifest)):
        raise RuntimeError("manifest still contains license material after redaction; refusing to write it")
    return manifest


def _version(text: str) -> str:
    m = re.search(r"^FragPipe(?:Plus)? version (\S+)$", text, re.MULTILINE)
    if not m:
        raise LogFormatError("no 'FragPipe version' line")
    return m.group(1)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("log", type=Path, help="FragPipe log_<timestamp>.txt (paths in it must still exist)")
    ap.add_argument("-o", "--output", type=Path, required=True, help="manifest.json to write")
    args = ap.parse_args(argv)
    manifest = build_manifest(args.log)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"{args.output}: {len(manifest['steps'])} steps, {len(manifest['unattributed_files'])} unattributed files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
