#!/usr/bin/env python3
"""Turn FragPipe's dry run (headless ``--write-plan``) into msproteomics' tool_configs.json.

A ``--write-plan`` run (a dry run that writes the plan) configures every tool as a real run does, so its work dir holds the
config files the commands read, FragPipe's saved run files (fragpipe.workflow, fragpipe-files.fp-manifest,
fragpipe.job) and the run log, which lists every command the run would start (FragpipeRun.java:412,
printProcessDescription). load_workdir reads that work dir into the plan the extractors use: each step's
argv (decode_command, told apart by the run's values; a line no single reading explains stops), its Cmd class
(step_kind), the config files its argv names, paths under the run's roots ({WORKDIR}, {LCMS_DIR},
{DB_DIR}, {FRAGPIPE}, {TMP}) and the licence as {LICENSE} (redact_lines). This script derives, per tool
key, the ``run`` / ``args`` / ``config_type`` and tool-specific fields the FragPipe subworkflows read
(the schema the retired workflow parser wrote; tests/python/test_plan_to_tool_configs.py compares with its frozen output).

``args`` holds what the msproteomics module passes through: the step's argv minus the elements the
module builds itself (binary, JVM options, threads, licence, input and output files).

Covered: the tools of FragPipe's workflows for TMT label check, LFQ-MBR, TMT quant, DIA SpecLib
Quant, glyco N-LFQ, split DB, open search + Crystal-C, diaPASEF diaTracer, DIA-Umpire, glyco O-Pair, as FragPipe 24.0 plans them):
MSFragger, MSBooster, Crystal-C, diaTracer, DIA-Umpire, Percolator, PeptideProphet, ProteinProphet,
Philosopher filter/report/abacus, IonQuant (label-free, TMT MS1 + isobaric), TMT-Integrator,
PTM-Shepherd, PairScans, O-Pair, AppendFile, database, SpecLibGen, DIA-NN, SAINTexpress,
FreeQuant, LabelQuant, PTMProphet, iProphet, MBG, FP-Meta, FPOP, FPOP coADAPTr, Write sub mzML, Export matched fragments.
A tool that has a plan step but no extractor stops the script, so nothing is silently dropped, and so do the
runs FragPipe plans but msproteomics cannot run (unsupported_runs). The legacy configs/ files the
subworkflows read (fragger*.params per data type, msbooster.config) are written next to the JSON.
"""

from __future__ import annotations

import argparse
import json
import re
import shlex
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fragpipe_log_to_manifest as fpm  # noqa: E402

# config_type per tool key, as the FragPipe subworkflows read it
CONFIG_TYPES = {
    "msfragger": "params_file", "msbooster": "params_file", "percolator": "cli", "peptideprophet": "cli",
    "proteinprophet": "cli", "filter": "cli", "ionquant": "cli", "tmtintegrator": "params_file",
    "freequant": "cli", "labelquant": "cli", "diaumpire": "params_file", "diatracer": "cli", "diann": "cli",
    "speclibgen": "bash_source", "ptmprophet": "cli", "ptmshepherd": "params_file", "crystalc": "params_file",
    "opair": "argv", "mbg": "argv", "fpop": "argv", "fpop_coadaptr": "cli", "saintexpress": "cli",
    "skyline": "cli", "metaproteomics": "argv", "transferlearning": "cli", "database": "cli",
    "tabrun": "gui_only", "exportmatchedfragments": "cli", "quantitation": "gui_only", "pairscans": "argv",
    "appendfile": "cli", "tmtintegrator_unmod": "params_file", "abacus": "cli", "iprophet": "cli",
    "writesubmzml": "cli",
}

# steps that belong to no tool key (FragPipe's own housekeeping)
# (Philosopher's database annotation is read by database_config)
INFRA_CLASSES = {"CmdCheckCentroid", "CmdPhilosopherWorkspaceCleanInit", "CmdPhilosopherWorkspaceClean", "CmdPhilosopherDbAnnotate"}


class PlanError(Exception):
    """The plan cannot be turned into tool configs."""


# CmdTransferLearning.java:175, :242 pass the Koina key; the console prints it masked (ProcessBuilderInfo.java:68-76)
MASKED_API_KEY = "********"


class Tokenizer:
    """Run paths under their roots, longest root first, at path boundaries."""

    def __init__(self, roots: dict[str, str]):
        self.roots = [(re.compile(re.escape(path) + r'(?=$|[/:"\s;,])'), name)
                      for name, path in sorted(roots.items(), key=lambda kv: len(kv[1]), reverse=True)]

    def path(self, s: str) -> str:
        for pattern, name in self.roots:
            s = pattern.sub(lambda _m, n=name: n, s)
        return s


def run_roots(lcms_files: list[str], job: list[str], argvs: list[list[str]]) -> dict[str, str]:
    """The run's path roots: {WORKDIR}, {LCMS_DIR} ({LCMS_DIR:n} per directory when there are several, in
    manifest order), {DB_DIR}, {FRAGPIPE} (the tools folder's parent; {TOOLS} when the tools folder is not
    {FRAGPIPE}/tools) and {TMP} (the directory of Philosopher's --temp workspaces)."""
    workdir, tools, fasta = (job[i].rstrip("/") for i in (3, 4, 5))
    roots = {"{WORKDIR}": workdir}
    lcms_dirs = list(dict.fromkeys(str(Path(f).parent) for f in lcms_files))
    roots.update({"{LCMS_DIR}": lcms_dirs[0]} if len(lcms_dirs) == 1 else {f"{{LCMS_DIR:{i}}}": d for i, d in enumerate(lcms_dirs, 1)})
    roots["{DB_DIR}"] = str(Path(fasta).parent)
    roots["{FRAGPIPE}"] = str(Path(tools).parent)
    if Path(tools) != Path(roots["{FRAGPIPE}"]) / "tools":
        roots["{TOOLS}"] = tools
    tmp = {str(Path(a[i + 1]).parent) for a in argvs for i in range(len(a) - 1) if a[i] == "--temp"}
    if len(tmp) > 1:
        raise PlanError(f"Philosopher --temp workspaces in several directories: {sorted(tmp)}")
    roots.update({"{TMP}": tmp.pop()} if tmp else {})
    return roots


def text_of(path: Path) -> str | None:
    """A work-dir file's text, byte-exact; None for a binary file (no tool params)."""
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        return None
    return None if "\x00" in text else text


def load_workdir(workdir: Path) -> ET.Element:
    """The plan of a ``--write-plan`` work dir: steps from its run log, configs from its files."""
    logs = sorted(workdir.glob("log_*.txt"))
    if len(logs) != 1:
        raise PlanError(f"{workdir}: expected one FragPipe run log log_*.txt, found {[p.name for p in logs]}")
    try:
        # the LC-MS paths of the manifest FragPipe saved, against which a log line split by another thread is rejoined
        manifest = frozenset(line.split("\t", 1)[0] for line in
                             (workdir / "fragpipe-files.fp-manifest").read_text(encoding="utf-8").splitlines() if line)
        parsed = fpm.parse_dry_run(logs[0].read_text(encoding="utf-8"), manifest)
    except fpm.LogFormatError as ex:
        raise PlanError(f"{logs[0].name}: {ex}") from ex
    # fragpipe.job: the run's work dir, tools folder and FASTA (FragpipeRun.java, BatchRun "runtime-save")
    jobs = [line.split("\t") for line in (workdir / "fragpipe.job").read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]
    if len(jobs) != 1:
        raise PlanError(f"{workdir}/fragpipe.job: expected one job line, found {len(jobs)}")
    # the run's free text an argv is built from: the saved workflow's values (FragPipe trims the cmd-opts
    # fields, e.g. DiannPanel.java:467-469) and the paths FragPipe wrote or was given
    files = sorted(p for p in workdir.rglob("*") if p.is_file() and p != logs[0])
    values = workflow_properties((workdir / "fragpipe.workflow").read_text(encoding="utf-8")).values()
    known = frozenset([*values, *(v.strip() for v in values), *(f["path"] for f in parsed["lcms_files"]), *jobs[0][3:6],
                       *(f"{jobs[0][3].rstrip('/')}/{p.relative_to(workdir).as_posix()}" for p in files)])
    argvs = []
    for order, s in enumerate(parsed["steps"], 1):
        try:
            a = fpm.decode_command(s["command_line"], known)
        except fpm.AmbiguousCommand as ex:
            raise PlanError(f"step {order} ({s['name']}): {ex}") from ex
        if any(a[i] == "--api-key" and a[i + 1] == MASKED_API_KEY for i in range(len(a) - 1)):
            raise PlanError(f"step {order} ({s['name']}): --api-key is masked in the log ({MASKED_API_KEY}), so its argv is unknown")
        argvs.append(a)
    roots = run_roots([f["path"] for f in parsed["lcms_files"]], jobs[0], argvs)
    tok = Tokenizer(roots)
    root = ET.Element("plan")
    roots_el = ET.SubElement(root, "roots")
    for name, path in roots.items():
        ET.SubElement(roots_el, "root", name=name, path=path)
    steps = []
    for order, (s, a) in enumerate(zip(parsed["steps"], argvs), 1):
        cls, tool = fpm.step_kind(s["name"])
        step = ET.SubElement(root, "step", order=str(order), tool=tool, **{"class": cls}, name=s["name"])
        if s["cwd_raw"] is not None:
            ET.SubElement(step, "workdir").text = tok.path(s["cwd_raw"])
        argv_el = ET.SubElement(step, "argv")
        for x in a:
            ET.SubElement(argv_el, "arg").text = tok.path(x)
        steps.append((step, a))
    # each text file of the work dir goes to the steps whose argv names it, else to the plan (FragPipe's run files)
    for path in files:
        text = text_of(path)
        if text is None:
            continue
        if fpm.contains_license_material(text):
            raise PlanError(f"{path.relative_to(workdir)} names the licence")
        run_path = f"{roots['{WORKDIR}']}/{path.relative_to(workdir).as_posix()}"
        holders = [step for step, a in steps if any(run_path in x for x in a)] or [root]
        for holder in holders:
            ET.SubElement(holder, "config", path=tok.path(run_path), encoding="text").text = text
    return root


def argv(step: ET.Element) -> list[str]:
    a = step.find("argv")
    if a is None:
        raise PlanError(f"step {step.get('order')} ({step.get('name')}) has no <argv>")
    return [x.text or "" for x in a]


def configs(step: ET.Element) -> dict[str, str]:
    """Text configs of a step, path -> content (base64 configs are never tool params)."""
    return {c.get("path"): (c.text or "") for c in step.findall("config") if c.get("encoding") == "text"}


def strip_options(args: list[str], with_value: set[str], flags: set[str] = frozenset()) -> list[str]:
    """argv without the given options (and the value after each option in ``with_value``)."""
    out, skip = [], False
    for a in args:
        if skip:
            skip = False
        elif a in with_value:
            skip = True
        elif a not in flags:
            out.append(a)
    return out


def same_args(tool: str, variants: list[list[str]]) -> list[str]:
    """All per-file steps of a tool must pass the same options; tool_configs holds one args string."""
    if any(v != variants[0] for v in variants[1:]):
        raise PlanError(f"{tool}: per-file steps pass different options: {variants}")
    return variants[0]


def by_name(steps: list[ET.Element], name: str) -> list[ET.Element]:
    return [s for s in steps if s.get("name") == name]


class Plan:
    """The parts of the plan the extractors read: steps, roots and every text config by path."""

    def __init__(self, root: ET.Element):
        self.root = root
        self.steps = root.findall("step")
        roots = root.find("roots")
        self.roots = {r.get("name"): r.get("path") for r in roots} if roots is not None else {}
        self.configs: dict[str, str] = {}
        for c in root.iter("config"):
            if c.get("encoding") == "text":
                self.configs[c.get("path")] = c.text or ""

    def config(self, path: str) -> str:
        """Content of the config at a (tokenised) plan path."""
        if path not in self.configs:
            raise PlanError(f"the plan has no config {path}")
        return self.configs[path]

    def config_named(self, name: str) -> list[tuple[str, str]]:
        return [(p, t) for p, t in self.configs.items() if Path(p).name == name]


def key_values(text: str) -> dict[str, str]:
    """``key = value  # comment`` lines of a FragPipe params file."""
    return dict(
        (line.split("=", 1)[0].strip(), line.split("=", 1)[1].split("#", 1)[0].strip())
        for line in text.splitlines()
        if "=" in line and not line.lstrip().startswith("#")
    )


# MsfraggerParams data_type -> FragPipe data type, and the params file FragPipe writes for it
# (CmdMsfragger.java:444-457: fragger.params, or fragger_dda.params when the run also has DIA, GPF-DIA or DIA-Lib)
FRAGGER_DATA_TYPES = {"0": "DDA", "1": "DIA", "2": "GPF-DIA", "3": "DDA+"}


def msfragger_config(steps: list[ET.Element], plan: Plan) -> dict:
    search = by_name(steps, "MSFragger")
    params: dict[str, set[str]] = {}
    for s in search:
        for path, text in configs(s).items():
            if Path(path).name.startswith("fragger") and path.endswith(".params"):
                params.setdefault(Path(path).name, set()).add(text)
    if not params:
        raise PlanError("msfragger: no fragger*.params in the plan's MSFragger steps")
    by_type = {}
    for name, texts in sorted(params.items()):
        if len(texts) != 1:
            raise PlanError(f"msfragger: the per-file MSFragger steps carry different {name}")
        text = texts.pop()
        values = key_values(text)
        data_type = FRAGGER_DATA_TYPES.get(values.get("data_type"))
        if data_type is None or values.get("calibrate_mass") not in ("0", "1", "2"):
            raise PlanError(f"msfragger: {name} has data_type={values.get('data_type')} calibrate_mass={values.get('calibrate_mass')}")
        by_type[data_type] = (name, text, values)
    split = [s for s in search if any(Path(a).name == "msfragger_pep_split.py" for a in argv(s))]
    if split:
        a = argv(split[0])
        num_slices = int(a[[Path(x).name for x in a].index("msfragger_pep_split.py") + 1])
    else:
        num_slices = 1
    # write_calibrated_mzml = 1 with calibration on gives <stem>_calibrated.mzML (CmdMsfragger.java:509-520)
    has_calibrated = any(v.get("write_calibrated_mzml") in ("1", "true") and int(v["calibrate_mass"]) > 0 for _n, _t, v in by_type.values())
    base = by_type.get("DDA") or next(iter(by_type.values()))
    return {
        "run": True,
        "args": base[1],
        "config_type": CONFIG_TYPES["msfragger"],
        "num_slices": num_slices,
        "calibrate_mass": {t: v["calibrate_mass"] for t, (_n, _x, v) in by_type.items()},
        "has_calibrated_mzml": has_calibrated,
        # params file per data type, as FragPipe wrote them (legacy configs/ outputs)
        "_files": {n: x for _t, (n, x, _v) in by_type.items()},
    }


# MSBooster keys CmdMSBooster.configure adds from the run (CmdMSBooster.java:160-212); msboosterParamsText re-adds them
MSBOOSTER_RUNTIME_KEYS = ("DiaNN", "fragger", "unimodObo", "mzmlDirectory", "pinPepXMLDirectory")
MSBOOSTER_LIBRARY_KEYS = ("spectraPredFile", "RTPredFile", "IMPredFile")
MSBOOSTER_LIBRARY_TOKEN = "{SPECTRAL_LIBRARY}"  # msboosterParamsText fills it with the staged library


def msbooster_config(steps: list[ET.Element], plan: Plan) -> dict:
    texts = {t for s in steps for p, t in configs(s).items() if Path(p).name == "msbooster_params.txt"}
    if len(texts) != 1:
        raise PlanError(f"msbooster: expected one msbooster_params.txt, plan has {len(texts)}")
    text = texts.pop()
    # CmdMSBooster.java:180-184 writes spectraPredFile/RTPredFile/IMPredFile = libraryPath only with a library;
    # MSBOOSTER stages the file (spectral_library) and msboosterParamsText fills the token
    library = key_values(text).get("spectraPredFile", "")
    lines = []
    for line in text.splitlines():
        key = line.split("=", 1)[0].strip()
        if key in MSBOOSTER_RUNTIME_KEYS:
            continue
        lines.append(f"{key} = {MSBOOSTER_LIBRARY_TOKEN}" if key in MSBOOSTER_LIBRARY_KEYS else line)
    text = "\n".join(lines)
    return {"run": True, "args": text, "config_type": CONFIG_TYPES["msbooster"], "spectral_library": library,
            "_files": {"msbooster.config": text + "\n"}}


def percolator_config(steps: list[ET.Element], plan: Plan) -> dict:
    variants = []
    for s in by_name(steps, "Percolator"):
        a = argv(s)[1:]  # binary
        # percolatorCommands (fragpipe_utils_b.nf:60-75) adds these, and the Write sub mzML weights: --weights
        # <stem>_percolator.weights on a first pass, --init-weights <weights> --static --override on a second
        # (CmdPercolator.java:179-201; the module's write_weights / init_weights inputs)
        if "--init-weights" in a:
            i = a.index("--init-weights")
            if a[i + 2:i + 4] != ["--static", "--override"]:
                raise PlanError(f"percolator: step {s.get('order')} has --init-weights without --static --override: {a}")
            a = a[:i] + a[i + 4:]
        a = strip_options(a, {"--num-threads", "--results-psms", "--decoy-results-psms", "--results-peptides",
                              "--decoy-results-peptides", "--protein-decoy-pattern", "--weights"})
        if not a or not a[-1].endswith(".pin"):
            raise PlanError(f"percolator: step {s.get('order')} does not end with its .pin input: {a}")
        variants.append(a[:-1])
    converts = by_name(steps, "Percolator: Convert to pepxml")
    # PercolatorOutputToPepXML <pin> <basename> <target tsv> <decoy tsv> <out> <data type> <min prob> <spectra>
    # (CmdPercolator.java:276-291)
    for s in converts:
        if argv(s)[-3] not in ("DDA", "DIA", "GPF-DIA", "DDA+", "DIA-Quant", "DIA-Lib"):
            raise PlanError(f"percolator: unexpected PercolatorOutputToPepXML argv in step {s.get('order')}: {argv(s)}")
    min_probs = {argv(s)[-2] for s in converts}
    if len(min_probs) != 1:
        raise PlanError(f"percolator: PercolatorOutputToPepXML min prob differs between files: {sorted(min_probs)}")
    return {
        "run": True,
        "args": " ".join(same_args("percolator", variants)),
        "config_type": CONFIG_TYPES["percolator"],
        "min_prob": min_probs.pop(),
        "keep_tsv_files": not by_name(steps, "Percolator delete temp"),
    }


def peptideprophet_config(steps: list[ET.Element], plan: Plan) -> dict:
    runs = by_name(steps, "PeptideProphet")
    variants = []
    for s in runs:
        a = strip_options(argv(s)[2:], {"--decoy", "--database"}, {"--combine"})  # philosopher peptideprophet
        if not a or not a[-1].lower().endswith((".pepxml", ".pep.xml")):
            raise PlanError(f"peptideprophet: step {s.get('order')} does not end with its pepXML input: {a}")
        variants.append([x for x in a if not x.lower().endswith((".pepxml", ".pep.xml"))])
    combine = {"--combine" in argv(s) for s in runs}
    if len(combine) != 1:
        raise PlanError("peptideprophet: some steps use --combine and some do not")
    return {
        "run": True,
        "args": " ".join(same_args("peptideprophet", variants)),
        "config_type": CONFIG_TYPES["peptideprophet"],
        "combine": combine.pop(),
    }


def proteinprophet_config(steps: list[ET.Element], plan: Plan) -> dict:
    a = argv(by_name(steps, "ProteinProphet")[0])
    a = strip_options(a[2:], {"--output"})  # philosopher proteinprophet
    return {"run": True, "args": " ".join(a[:-1]), "config_type": CONFIG_TYPES["proteinprophet"]}  # filelist


def filter_config(steps: list[ET.Element], plan: Plan) -> dict:
    filters = by_name(steps, "PhilosopherFilter")
    # --tag, --pepxml, --dbbin, --protxml, --probin, --razor: philosopherFilterCommands (fragpipe_utils_b.nf:130-167)
    variants = [strip_options(argv(s)[2:], {"--tag", "--pepxml", "--dbbin", "--protxml", "--probin"}, {"--razor"}) for s in filters]
    reports = [argv(s)[2:] for s in by_name(steps, "PhilosopherReport")]
    return {
        "run": True,
        "args": " ".join(same_args("filter", variants)),
        "config_type": CONFIG_TYPES["filter"],
        "dont_use_protxml": not any("--protxml" in argv(s) for s in filters),
        "report_args": " ".join(same_args("filter report", reports)),
        "report_args_single_group": "",
    }


def abacus_config(steps: list[ET.Element], plan: Plan) -> dict:
    a = argv(steps[0])[2:]  # philosopher abacus
    # the experiment group names close the argv; abacusCommand (fragpipe_utils_b.nf:95-112) adds them
    while a and not a[-1].startswith("-"):
        a = a[:-1]
    entry = {"run": True, "args": " ".join(strip_options(a, {"--plex"}, {"--labels"})), "config_type": CONFIG_TYPES["abacus"]}
    # --labels --plex <reagents of the TMT-Integrator label> with LabelQuant intensities (FragpipeRun.java:1793-1797,
    # CmdPhilosopherAbacus.java:131-135); FRAGPIPE_INFERENCE takes the label name, the workflow's tmtintegrator.channel_num
    if "--labels" in a:
        entry["labels_channel_num"] = workflow_value(plan, "tmtintegrator.channel_num")
    return entry


def abacus_off_config(props: dict[str, str]) -> dict:
    """Abacus options when FragPipe's own conditions do not run it: the pipeline still runs Abacus for
    SAINTexpress's reprint tables (docs/usage.md, SAINTexpress). CmdPhilosopherAbacus.java:106-130 options
    from the saved workflow: --picked/--uniqueonly of the filter text, --razor --reprint, --tag, and the
    protein/peptide summaries (FragPipe saves every key, so none is defaulted here)."""
    filter_text = props["phi-report.filter"]
    addons = [t for t in dict.fromkeys(filter_text.split()) if t in ("--picked", "--uniqueonly")]
    addons += [t for t in ("--razor", "--reprint") if t not in addons]
    args = addons + ["--tag", props["database.decoy-tag"]]

    def flag(key: str) -> bool:
        return props[key].strip().lower() == "true"

    if not flag("phi-report.dont-use-prot-proph-file") and flag("phi-report.prot-level-summary"):
        args.append("--protein")
    if flag("phi-report.pep-level-summary"):
        args.append("--peptide")
    return {"run": False, "args": " ".join(args), "config_type": CONFIG_TYPES["abacus"]}


METAPROTEOMICS_OPTIONS = ("qvalue", "deltaHyperscore", "minPeptCntPerProt", "minUniqPeptCntPerProt", "minUniqPeptCnt", "hostName", "iterations")


def metaproteomics_config(steps: list[ET.Element], plan: Plan) -> dict:
    """FP-Meta DbOptimizer options (CmdMetaproteomics.java:60-110); the project dir, FASTA, output dir,
    decoy tag and taxonomy files are the task's. cmd-opts close the argv after --iterations."""
    a = argv(steps[0])
    args = {key: a[a.index(f"--{key}") + 1] for key in METAPROTEOMICS_OPTIONS}
    rest = a[a.index("--iterations") + 2:]
    args["cmdLineOpts"] = " ".join(rest)
    return {"run": True, "args": args, "config_type": CONFIG_TYPES["metaproteomics"]}


def fpop_config(steps: list[ET.Element], plan: Plan) -> dict:
    a = argv(steps[0])
    a = a[[Path(x).name for x in a].index("FragPipe_FPOP_Analysis.py") + 1:]
    # the tables FragPipe resolves in its work dir (CmdFpopQuant.java:56-63, :91-98, :114-145): the first and, for TMT,
    # the last argument, relative to the work dir as the FPOP task reads them
    a = [x[len("{WORKDIR}/"):] if x.startswith("{WORKDIR}/") else x for x in a]
    return {"run": True, "args": a, "config_type": CONFIG_TYPES["fpop"]}


def fpop_coadaptr_config(steps: list[ET.Element], plan: Plan) -> dict:
    a = argv(steps[0])
    return {"run": True, "args": a[-1], "config_type": CONFIG_TYPES["fpop_coadaptr"]}  # FPOP_PSMwriter <work dir> <mass offsets>


def _ionquant_entry(step: ET.Element, plan: Plan) -> dict:
    a = argv(step)
    a = a[a.index("ionquant.IonQuant") + 1:]
    modlist = a[a.index("--modlist") + 1] if "--modlist" in a else None
    a = strip_options(a, {"--license", "--threads", "--multidir", "--filelist", "--modlist", "--annotation"})
    # FragPipe passes each value as one argv element (CmdIonquant.java:303-304); msproteomics splits the
    # IonQuant options on spaces (shellArgs), so a value with whitespace would become several
    spaced = [x for x in a if re.search(r"\s", x)]
    if spaced:
        raise PlanError(f"ionquant: IonQuant values with whitespace {spaced}: FragPipe passes each as one argument, "
                        "but msproteomics splits the IonQuant options on spaces")
    return {
        "run": True,
        "args": " ".join(a),
        "config_type": CONFIG_TYPES["ionquant"],
        "modmasses": [m for m in plan.config(modlist).splitlines() if m] if modlist else [],
    }


def ionquant_configs(steps: list[ET.Element], plan: Plan) -> dict:
    """One IonQuant step: label-free or TMT label check (key ionquant). Two: TMT quant, MS1 then isobaric
    (ionquant_ms1, ionquant_isobaric; FragpipeRun cmdTmtIonquant and cmdTmtIonquantIsobaric)."""
    runs = by_name(steps, "IonQuant")
    if len(runs) == 1:
        return {"ionquant": _ionquant_entry(runs[0], plan)}
    if len(runs) == 2:
        out = {}
        for s in runs:
            a = argv(s)
            iso = a[a.index("--perform-isoquant") + 1] == "1"
            out["ionquant_isobaric" if iso else "ionquant_ms1"] = _ionquant_entry(s, plan)
        if set(out) != {"ionquant_ms1", "ionquant_isobaric"}:
            raise PlanError(f"ionquant: two IonQuant steps but not one MS1 and one isobaric: {sorted(out)}")
        return out
    raise PlanError(f"ionquant: {len(runs)} IonQuant steps")


def tmtintegrator_config(steps: list[ET.Element], plan: Plan) -> dict:
    texts = {t for s in steps for p, t in configs(s).items() if Path(p).name == "tmt-integrator-conf.yml"}
    if len(texts) != 1:
        raise PlanError(f"tmtintegrator: expected one tmt-integrator-conf.yml, plan has {len(texts)}")
    # output: is the run's report dir; TMTINTEGRATOR fills @TMT_REPORT_DIR@
    lines = ["  output: @TMT_REPORT_DIR@" if line.startswith("  output: ") else line for line in texts.pop().splitlines()]
    text = "\n".join(lines) + "\n"
    labelquant = any(s.get("class") == "CmdLabelquant" for s in plan.steps)
    isobaric = any(s.get("class") == "CmdIonquant" and "--perform-isoquant" in argv(s)
                   and argv(s)[argv(s).index("--perform-isoquant") + 1] == "1" for s in plan.steps)
    return {
        "run": True,
        "args": text,
        "config_type": CONFIG_TYPES["tmtintegrator"],
        # TmtiPanel intensity extraction: 0 IonQuant, 1 Philosopher (LabelQuant), 2 none
        "extraction_tool": 0 if isobaric else 1 if labelquant else 2,
        # Philosopher MSstats output writes MSstatsTMT_annotation.csv (CmdTmtIntegrator.java:107)
        "philosopher_msstats": bool(plan.config_named("MSstatsTMT_annotation.csv")),
    }


PTMS_ANNOTATION_KEYWORDS = ("", "unimod", "common", "glyco")


def ptmshepherd_config(steps: list[ET.Element], plan: Plan) -> dict:
    texts = {t for s in steps for p, t in configs(s).items() if Path(p).name == "shepherd.config"}
    if len(texts) != 1:
        raise PlanError(f"ptmshepherd: expected one shepherd.config, plan has {len(texts)}")
    tools = f"{plan.roots['{FRAGPIPE}']}/lib/../tools"
    text = texts.pop()
    # PtmshepherdPanel.java:225-232 writes the keyword or the custom file's path; the keywords are those
    # loadDefaults() reads as such (PtmshepherdPanel.java:264-278). PTMSHEPHERD stages a custom file
    # (annotation_file) and fills the token.
    annotation = key_values(text).get("annotation_file", "")
    annotation = "" if annotation.lower() in PTMS_ANNOTATION_KEYWORDS else annotation
    lines = []
    for line in text.splitlines():
        key = line.split("=", 1)[0].strip()
        if not line.strip() or key in ("database", "dataset"):
            continue  # run inputs; ptmshepherdConfigText adds them
        if key == "threads":
            line = "threads = @THREADS@"
        if key == "annotation_file" and annotation:
            line = "annotation_file = @ANNOTATION_FILE@"
        lines.append(line.replace(tools, "@TOOLS_DIR@"))
    return {"run": True, "args": "\n".join(lines), "config_type": CONFIG_TYPES["ptmshepherd"], "annotation_file": annotation}


def appendfile_config(steps: list[ET.Element], plan: Plan) -> dict:
    return {"run": True, "args": "", "config_type": CONFIG_TYPES["appendfile"]}


def database_config(props: dict[str, str], plan: Plan) -> dict:
    """The decoy tag of the run (database.decoy-tag), checked against Philosopher's database annotation
    step when it runs (CmdPhilosopherDbAnnotate.java:58-74)."""
    tag = props["database.decoy-tag"]
    for s in plan.steps:
        if s.get("class") == "CmdPhilosopherDbAnnotate":
            a = argv(s)
            if a[a.index("--prefix") + 1] != tag:
                raise PlanError(f"database: PhilosopherDbAnnotate --prefix {a[a.index('--prefix') + 1]} differs from database.decoy-tag={tag}")
    return {"run": True, "args": f"--prefix {tag}", "config_type": CONFIG_TYPES["database"], "decoy_tag": tag}


# gen_con_spec_lib.py argv after the script (CmdSpecLibGen.java:137-202): 1 fasta, 2 work dir, 3 unused,
# 4 work dir, 5 overwrite, 6 unused, 7 use_easypqp, 8 "<rt cal>:<im cal>", 9 threads, 10 convert extras,
# 11 library extras, 12 keep/delete intermediate files, 13 filelist
def speclibgen_config(steps: list[ET.Element], plan: Plan) -> dict:
    a = argv(steps[0])
    i = [Path(x).name for x in a].index("gen_con_spec_lib.py")
    args = a[i + 1:]
    if len(args) != 13 or args[6] != "use_easypqp":
        raise PlanError(f"speclibgen: unexpected gen_con_spec_lib.py argv: {args}")
    rt_cal, im_cal = args[7].split(":", 1)
    parts = {"rt_cal": rt_cal, "im_cal": im_cal, "convert_args": args[9], "library_args": args[10], "intermediate": args[11]}
    rt_file = im_file = ""
    for key, value in (("rt_cal", rt_cal), ("im_cal", im_cal)):
        if "/" in value:  # a user calibration file
            raise PlanError(f"speclibgen: {key} calibration file {value} is not handled yet")
    assignments = (("SPECLIB_RT_CAL", parts["rt_cal"]), ("SPECLIB_IM_CAL", parts["im_cal"]),
                   ("SPECLIB_CONVERT_ARGS", parts["convert_args"]), ("SPECLIB_LIBRARY_ARGS", parts["library_args"]),
                   ("SPECLIB_INTERMEDIATE", parts["intermediate"]))
    return {
        "run": True,
        "args": "\n".join(f"{name}={shlex.quote(value)}" for name, value in assignments),
        "config_type": CONFIG_TYPES["speclibgen"],
        "rt_cal_file": rt_file,
        "im_cal_file": im_file,
        "argv_parts": parts,
    }


def diann_config(steps: list[ET.Element], plan: Plan) -> dict:
    """DIA-NN run, MSstats and site-report steps of CmdDiann (CmdDiann.java:190-420)."""
    run = [s for s in steps if s.get("name").endswith("run DIA-NN")]
    if len(run) != 1:
        raise PlanError(f"diann: expected one 'run DIA-NN' step, plan has {len(run)}")
    a = argv(run[0])[1:]
    if not ("--qvalue" in a and "--matrices" in a and "--smart-profiling" in a):
        raise PlanError(f"diann: unexpected DIA-NN argv: {a}")
    # --qvalue .. --matrices (generate_diann_config), then the protein-inference flags DIANN_FRAGPIPE builds,
    # then --smart-profiling .. before --cfg (generate_diann_post_args)
    args = a[a.index("--qvalue"):a.index("--matrices") + 1]
    post = a[a.index("--smart-profiling"):a.index("--cfg")]
    # diann.cmd-opts (DiannPanel.java:467-469, trimmed) is one argv element right before --cfg
    # (CmdDiann.java:385-387, :398); DIANNFRAGPIPE adds it from cmd_opts
    cmd_opts = workflow_value(plan, "diann.cmd-opts").strip()
    if cmd_opts:
        if post[-1] != cmd_opts:
            raise PlanError(f"diann: diann.cmd-opts {cmd_opts!r} is not the argv element before --cfg: {post}")
        post = post[:-1]
    lib = a[a.index("--lib") + 1]
    if lib != "library.tsv":
        raise PlanError(f"diann: a user or predicted library ({lib}) is not handled yet")
    if any(x in a for x in ("--channels", "--fixed-mod", "--lib-fixed-mod")):
        raise PlanError("diann: plexDIA is not handled yet")
    reports = [s for s in steps if s.get("name").endswith("generate site reports")]
    mod_tag, min_site_prob, levels = "", "0.75", "all"
    if reports:
        r = argv(reports[0])
        mod_tag = r[r.index("--mod-tag") + 1]
        min_site_prob = r[r.index("--min-site-prob") + 1]
        levels = r[r.index("--level") + 1]
    return {
        "run": True,
        "args": " ".join(args),
        "config_type": CONFIG_TYPES["diann"],
        "args_post": " ".join(post),
        "qvalue": args[args.index("--qvalue") + 1],
        "run_specific_protein_qvalue": "--matrix-spec-q" in args,
        "redo_protein_inference": any("protein inference" in s.get("name", "") for s in steps),
        "generate_msstats": any(s.get("name").endswith("convert DIA-NN output to MSstats.csv") for s in steps),
        "library": "",
        "run_plex": False,
        "light": "",
        "medium": "",
        "heavy": "",
        "mod_tag": mod_tag,
        "min_site_prob": min_site_prob,
        "report_levels": levels,
        "cmd_opts": cmd_opts,
    }


def saved_workflow(plan: Plan) -> dict[str, str]:
    """The run's settings FragPipe saved into its work dir (fragpipe.workflow, a top-level plan config),
    as Java properties (\\: \\= \\\\ unescaped)."""
    workflows = plan.config_named("fragpipe.workflow")
    if len(workflows) != 1:
        raise PlanError(f"the plan has {len(workflows)} saved fragpipe.workflow configs, expected 1")
    return workflow_properties(workflows[0][1])


def _java_properties_unescape(text: str) -> str:
    """Properties.loadConvert: \\uXXXX, \\t, \\n, \\r, \\f, and any other escaped character as itself
    (PropertiesUtils.saveConvert writes them, PropertiesUtils.java:111-167)."""
    out, i = [], 0
    while i < len(text):
        if text[i] == "\\" and i + 1 < len(text):
            c = text[i + 1]
            if c == "u":
                out.append(chr(int(text[i + 2:i + 6], 16)))
                i += 6
                continue
            out.append({"t": "\t", "n": "\n", "r": "\r", "f": "\f"}.get(c, c))
            i += 2
            continue
        out.append(text[i])
        i += 1
    return "".join(out)


def _java_properties_entry(line: str) -> tuple[str, str]:
    """Properties.load0 on one logical line: the key ends at the first unescaped '=', ':' or
    whitespace; the value starts after the whitespace and at most one '=' / ':' that follow."""
    i, escaped = 0, False
    while i < len(line):
        c = line[i]
        if not escaped and c in "=: \t\f":
            break
        escaped = c == "\\" and not escaped
        i += 1
    key_end, has_sep = i, i < len(line) and line[i] in "=:"
    i += 1 if i < len(line) else 0
    while i < len(line):
        c = line[i]
        if c not in " \t\f":
            if has_sep or c not in "=:":
                break
            has_sep = True
        i += 1
    return _java_properties_unescape(line[:key_end]), _java_properties_unescape(line[i:])


def workflow_properties(text: str) -> dict[str, str]:
    """java.util.Properties.load of a .workflow text (fragpipe/api/PropsFile.java:55-63), undoing everything
    Properties.store writes: natural lines joined at an odd number of trailing backslashes, leading
    whitespace dropped, blank lines and '#' / '!' comment lines skipped, the last value of a repeated key kept."""
    props = {}
    logical = None
    for raw in text.splitlines():
        line = raw.lstrip(" \t\f")
        if logical is None and (not line or line[0] in "#!"):
            continue
        logical = (logical or "") + line
        if (len(line) - len(line.rstrip("\\"))) % 2:
            logical = logical[:-1]
            continue
        key, value = _java_properties_entry(logical)
        props[key] = value
        logical = None
    if logical is not None:
        key, value = _java_properties_entry(logical)
        props[key] = value
    return props


def saintexpress_config(props: dict[str, str]) -> dict:
    """SAINTexpress runs from FragPipe's Downstream tab after the run (TabDownstream.runSaintExpress), so
    the plan has no step for it; its settings are the saved workflow's: -R max-replicates, -L virtual-controls
    (CmdSaintExpress.java:85-87, SaintexpressPanel.java:111-115); cmd-opts is one more argv element
    (CmdSaintExpress.java:88-90)."""
    return {
        "run": True,
        "args": f"-R{int(props.get('saintexpress.max-replicates', '10'))} -L{int(props.get('saintexpress.virtual-controls', '100'))}",
        "config_type": CONFIG_TYPES["saintexpress"],
        "cmd_opts": props.get("saintexpress.cmd-opts", "").strip(),
    }


# Crystal-C keys CmdCrystalc.configure sets from the run (CmdCrystalc.java:190-196); crystalcParamsText adds them
CRYSTALC_RUNTIME_KEYS = ("thread", "fasta", "raw_file_location", "raw_file_extension", "output_location")


def crystalc_config(steps: list[ET.Element], plan: Plan) -> dict:
    """The algorithm settings of the per-pepXML crystalc-<i>-<pepXML>.params, as key = value without comments."""
    variants = []
    for s in by_name(steps, "Crystal-C"):
        texts = [t for p, t in configs(s).items() if Path(p).name.startswith("crystalc-") and p.endswith(".params")]
        if len(texts) != 1:
            raise PlanError(f"crystalc: step {s.get('order')} has {len(texts)} crystalc params files")
        lines = []
        for line in texts[0].splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            key, _, value = line.partition("=")
            if key.strip() in CRYSTALC_RUNTIME_KEYS:
                continue
            lines.append(f"{key.strip()} = {value.split('#', 1)[0].strip()}")
        variants.append(lines)
    return {"run": True, "args": "\n".join(same_args("crystalc", variants)), "config_type": CONFIG_TYPES["crystalc"]}


def diatracer_config(steps: list[ET.Element], plan: Plan) -> dict:
    variants = []
    for s in by_name(steps, "diaTracer"):
        a = argv(s)
        a = a[[Path(x).name for x in a].index(next(x for x in (Path(y).name for y in a) if x.startswith("diaTracer"))) + 1:]
        # dataTracerCommand builds the license, input, work dir and thread options (CmdDiaTracer.java:90-110)
        variants.append(strip_options(a, {"--license", "--dFilePath", "--workDir", "--threadNum"}))
    return {"run": True, "args": " ".join(same_args("diatracer", variants)), "config_type": CONFIG_TYPES["diatracer"]}


def diaumpire_config(steps: list[ET.Element], plan: Plan) -> dict:
    """umpire-se.params with Thread = THREADS (the DIAUMPIRE task fills it); q_tiers from the pseudo-MS/MS
    files MSFragger searches (<stem>_Q1/_Q2/_Q3.mzML, CmdUmpireSe.java:59-75)."""
    texts = {t for s in steps for p, t in configs(s).items() if Path(p).name == "umpire-se.params"}
    if len(texts) != 1:
        raise PlanError(f"diaumpire: expected one umpire-se.params, plan has {len(texts)}")
    # the file's text as FragPipe writes it, final newline included (PropertiesUtils.java:291-321)
    text = re.sub(r"^Thread = .*$", "Thread = THREADS", texts.pop(), flags=re.MULTILINE)
    tiers = sorted({m.group(1) for s in plan.steps if s.get("class") == "CmdMsfragger" for a in argv(s)
                    for m in [re.search(r"_(Q[123])\.mzML$", a)] if m})
    return {"run": True, "args": text, "config_type": CONFIG_TYPES["diaumpire"], "q_tiers": tiers}


def pairscans_config(steps: list[ET.Element], plan: Plan) -> dict:
    variants = []
    for s in steps:
        a = argv(s)
        a = a[a.index("org.nesvilab.fragpipe.util.PairScans") + 3:]  # after the LC-MS file and the thread count
        variants.append(a)
    return {"run": True, "args": same_args("pairscans", variants), "config_type": CONFIG_TYPES["pairscans"]}


# O-Pair options the OPAIR task fills from the run (CmdOPair.java:126-160)
OPAIR_RUNTIME_OPTIONS = {"-t": "@THREADS@", "-r": "@FILELIST@", "-s": "@PSM@", "-o": "@OUTDIR@"}


def opair_config(steps: list[ET.Element], plan: Plan) -> dict:
    tools = "{FRAGPIPE}/lib/../tools"
    variants = []
    oxonium = set()
    for s in steps:
        a = argv(s)[2:]  # dotnet CMD.dll
        out = []
        for x in a:
            flag = x.split(" ", 1)[0]
            if flag in OPAIR_RUNTIME_OPTIONS:
                out.append(f"{flag} {OPAIR_RUNTIME_OPTIONS[flag]}")
            elif flag == "-f" and not x.startswith(f"-f {tools}/"):
                # CmdOPair.java:178-192: the user's oxonium file (the default is tools/opair/Glycan_Mods);
                # FRAGPIPE_GLYCO stages it (oxonium_file) and fills the token
                oxonium.add(x[len("-f "):])
                out.append("-f @OXONIUM_FILE@")
            else:
                out.append(x.replace(tools, "@TOOLS_DIR@"))
        variants.append(out)
    if len(oxonium) > 1:
        raise PlanError(f"opair: groups use different oxonium files {sorted(oxonium)}")
    lists = [t for s in steps for p, t in configs(s).items() if Path(p).name == "filelist_opair.txt"]
    if not lists:
        raise PlanError("opair: no filelist_opair.txt in the plan")
    names = [line for t in lists for line in t.splitlines() if line]
    # CmdOPair.java:100-116: _calibrated.mzML when MSFragger writes calibrated spectra, else the uncalibrated/input file
    spectra = "calibrated" if all(n.endswith("_calibrated.mzML") for n in names) else "uncalibrated_or_input"
    return {"run": True, "args": same_args("opair", variants), "config_type": CONFIG_TYPES["opair"], "spectra": spectra,
            "oxonium_file": oxonium.pop() if oxonium else ""}


def workflow_value(plan: Plan, key: str) -> str:
    """A key of the run's fragpipe.workflow in the plan, for the settings no command carries."""
    props = saved_workflow(plan)
    if key not in props:
        raise PlanError(f"the plan's fragpipe.workflow has no {key}")
    return props[key]


def freequant_config(steps: list[ET.Element], plan: Plan) -> dict:
    """philosopher freequant <options> [--raw] --dir <LC-MS dir> (CmdFreequant.java:112-125). Before LabelQuant
    it is the TMT FreeQuant (FragpipeRun.java:1932, 2049-2055): the label-free FreeQuant (key run) does not run,
    but its options are the ones the TMT FreeQuant passes."""
    variants = [strip_options(argv(s)[2:], {"--dir"}, {"--raw"}) for s in steps]  # philosopher freequant
    labelquant = [s for s in plan.steps if s.get("class") == "CmdLabelquant"]
    # label-free FreeQuant as well (FragpipeRun.java:1813) would add a second FreeQuant step per group
    if labelquant and len(steps) != len(labelquant):
        raise PlanError(f"freequant: {len(steps)} FreeQuant steps for {len(labelquant)} LabelQuant steps; "
                        "label-free FreeQuant together with LabelQuant is not handled yet")
    return {"run": not labelquant, "args": " ".join(same_args("freequant", variants)), "config_type": CONFIG_TYPES["freequant"]}


def labelquant_config(steps: list[ET.Element], plan: Plan) -> dict:
    """philosopher labelquant --tol --level --minprob --purity --removelow --plex --annot <file> --brand <type>
    [--raw] --dir <LC-MS dir> (CmdLabelquant.java:95-153); PHILOSOPHER_LABELQUANT adds --annot, --brand, --raw, --dir."""
    variants, brands = [], set()
    for s in steps:
        a = argv(s)[2:]  # philosopher labelquant
        if "--brand" not in a:
            raise PlanError(f"labelquant: step {s.get('order')} has no --brand")
        brands.add(a[a.index("--brand") + 1])
        variants.append(strip_options(a, {"--annot", "--brand", "--dir"}, {"--raw"}))
    if len(brands) != 1:
        raise PlanError(f"labelquant: the groups pass different --brand values: {sorted(brands)}")
    return {"run": True, "args": " ".join(same_args("labelquant", variants)), "config_type": CONFIG_TYPES["labelquant"],
            "brand": brands.pop()}


def ptmprophet_config(steps: list[ET.Element], plan: Plan) -> dict:
    """PTMProphetParser <options> MAXTHREADS=1 <pepXML> <mod.pep.xml> (CmdPtmProphet.java:114-122): FragPipe drops any
    MAXTHREADS= option and appends MAXTHREADS=1; the PTMPROPHET module appends it and the two file names."""
    variants = []
    for s in steps:
        a = argv(s)[1:]  # PTMProphetParser
        if len(a) < 3 or a[-3] != "MAXTHREADS=1" or not a[-1].endswith(".mod.pep.xml"):
            raise PlanError(f"ptmprophet: step {s.get('order')} does not end with MAXTHREADS=1 <pepXML> <mod.pep.xml>: {a}")
        variants.append(a[:-3])
    return {"run": True, "args": " ".join(same_args("ptmprophet", variants)), "config_type": CONFIG_TYPES["ptmprophet"]}


def iprophet_config(steps: list[ET.Element], plan: Plan) -> dict:
    """philosopher iprophet --decoy <tag> --nonsp --output combined --threads <n> <pepXMLs> (CmdIprophet.java:48-57);
    --threads and the pepXMLs are the IPROPHET task's."""
    if len(steps) != 1:
        raise PlanError(f"iprophet: expected one iProphet step, plan has {len(steps)}")
    a = strip_options(argv(steps[0])[2:], {"--threads"})  # philosopher iprophet
    return {"run": True, "args": " ".join(x for x in a if not x.lower().endswith((".pep.xml", ".pepxml"))),
            "config_type": CONFIG_TYPES["iprophet"]}


# MBG options the MBG task fills from the run (CmdMBGMatch.java:137-208), as tokens
MBG_RUNTIME_OPTIONS = {"--psm": "@PSM@", "--manifest": "@MANIFEST@", "--nopasef": "@NOPASEF@", "--numthreads": "@THREADS@"}


def mbg_config(steps: list[ET.Element], plan: Plan) -> dict:
    """MBG argv after com.mbg.MBG, one per group, the run's files and threads as tokens. im_ms is the workflow's
    IM-MS button (workflow.input.data-type.im-ms): FragPipe's --nopasef is !timsTOF, where the button decides only
    when no LC-MS file is .d or .raw (TabWorkflow.java:521-538), so the MBG task needs the button itself."""
    tools = "{FRAGPIPE}/lib/../tools"
    variants = []
    for s in steps:
        a = argv(s)
        if "com.mbg.MBG" not in a:
            raise PlanError(f"mbg: step {s.get('order')} does not run com.mbg.MBG: {a}")
        a = a[a.index("com.mbg.MBG") + 1:]
        out = []
        for i, x in enumerate(a):
            out.append(MBG_RUNTIME_OPTIONS[a[i - 1]] if i and a[i - 1] in MBG_RUNTIME_OPTIONS else x.replace(tools, "@TOOLS_DIR@"))
        variants.append(out)
    return {"run": True, "args": same_args("mbg", variants), "config_type": CONFIG_TYPES["mbg"],
            "im_ms": workflow_value(plan, "workflow.input.data-type.im-ms").strip().lower() == "true"}


def writesubmzml_config(steps: list[ET.Element], plan: Plan) -> dict:
    """WriteSubMzml <LC-MS file> <psm.tsv> <_sub.mzML> <probability threshold> 1 per LC-MS file
    (CmdWriteSubMzml.java:122-126); the threshold is the one setting, and writeSubMzmlArgv builds the rest.
    The second pass on the _sub.mzML files is the pipeline's (FragpipeRun.java:589-700)."""
    thresholds = set()
    for s in steps:
        a = argv(s)
        cls = next(i for i, x in enumerate(a) if x.endswith(".WriteSubMzml"))
        thresholds.add(a[cls + 4])
    if len(thresholds) != 1:
        raise PlanError(f"writesubmzml: the steps have different probability thresholds {sorted(thresholds)}")
    return {"run": True, "args": "", "config_type": CONFIG_TYPES["writesubmzml"], "prob_threshold": thresholds.pop()}


def exportmatchedfragments_config(steps: list[ET.Element], plan: Plan) -> dict:
    """FragmentsExportMainClass <work dir> <threads> r: no options from the workflow (CmdExportMatchedFragments.java:58-66)."""
    return {"run": True, "args": "", "config_type": CONFIG_TYPES["exportmatchedfragments"]}


EXTRACTORS = {
    "msfragger": msfragger_config, "msbooster": msbooster_config, "percolator": percolator_config,
    "peptideprophet": peptideprophet_config, "proteinprophet": proteinprophet_config, "filter": filter_config,
    "abacus": abacus_config, "tmtintegrator": tmtintegrator_config, "ptmshepherd": ptmshepherd_config,
    "appendfile": appendfile_config, "speclibgen": speclibgen_config,
    "diann": diann_config, "crystalc": crystalc_config, "diatracer": diatracer_config, "diaumpire": diaumpire_config,
    "pairscans": pairscans_config, "opair": opair_config,
    "freequant": freequant_config, "labelquant": labelquant_config, "ptmprophet": ptmprophet_config,
    "iprophet": iprophet_config, "mbg": mbg_config,
    "metaproteomics": metaproteomics_config, "fpop": fpop_config, "fpop_coadaptr": fpop_coadaptr_config,
    "writesubmzml": writesubmzml_config, "exportmatchedfragments": exportmatchedfragments_config,
}
# FragPipe Cmd class -> tool key, for every class whose steps this script understands
CLASS_TO_TOOL = {
    "CmdMsfragger": "msfragger", "CmdMSBooster": "msbooster", "CmdPercolator": "percolator",
    "CmdPeptideProphet": "peptideprophet", "CmdProteinProphet": "proteinprophet",
    "CmdPhilosopherFilter": "filter", "CmdPhilosopherReport": "filter", "CmdPhilosopherAbacus": "abacus",
    "CmdIonquant": "ionquant", "CmdTmtIntegrator": "tmtintegrator", "CmdPtmshepherd": "ptmshepherd",
    "CmdAppendFile": "appendfile", "CmdSpecLibGen": "speclibgen",
    "CmdDiann": "diann", "CmdCrystalc": "crystalc", "CmdDiaTracer": "diatracer", "CmdUmpireSe": "diaumpire",
    "CmdPairScans": "pairscans", "CmdOPair": "opair",
    "CmdFreequant": "freequant", "CmdLabelquant": "labelquant", "CmdPtmProphet": "ptmprophet", "CmdIprophet": "iprophet",
    "CmdMBGMatch": "mbg",
    "CmdMetaproteomics": "metaproteomics", "CmdFpopQuant": "fpop", "CmdFPOPcoadaptr": "fpop_coadaptr",
    "CmdWriteSubMzml": "writesubmzml", "CmdExportMatchedFragments": "exportmatchedfragments",
}
# FragPipe steps msproteomics does not run: FragPipe plans them, and the run stops here instead
UNSUPPORTED_CLASSES = {
    # CmdTransferLearning needs a credential file with the URL and API key of an external prediction server
    # (CmdTransferLearning.java:123-147); CmdMsfraggerDigest runs only with transfer learning (FragpipeRun.java:2255)
    "CmdTransferLearning": "Transfer learning (transfer-learning.run-transfer-learning=true) needs an external prediction "
                           "server, and msproteomics does not run it. Please disable Transfer learning in the workflow.",
    "CmdMsfraggerDigest": "Transfer learning's MSFragger digest (transfer-learning.run-transfer-learning=true) needs an external "
                          "prediction server, and msproteomics does not run it. Please disable Transfer learning in the workflow.",
}
# Steps that read Philosopher report's <group>/psm.tsv: CmdIonquant.java:318, CmdTmtIntegrator.java:295,
# PtmshepherdParams.java:85, CmdMBGMatch.java:139, CmdOPair.java:176, CmdSpecLibGen.java:135
PSM_TSV_READERS = ("CmdIonquant", "CmdTmtIntegrator", "CmdPtmshepherd", "CmdMBGMatch", "CmdOPair", "CmdSpecLibGen")


def unsupported_runs(root: ET.Element) -> list[str]:
    """The runs FragPipe plans but msproteomics cannot run, each with its message.

    - MSFragger off or digest only: FragPipe skips the search and hands the MSFragger outputs of its work
      directory to the later steps (FragpipeRun.java:1519-1521, outside `if isRun()`; TabMsfragger.java:1932-1934),
      i.e. it reuses a previous run's search; msproteomics starts from LC-MS files.
    - FreeQuant without the Report panel: FragPipe runs Philosopher report (FragpipeRun.java:1767,
      reportPanel.isRun() || isRunFreeQuant()) with no DbAnnotate (:1690) or Filter (:1702); msproteomics
      reports after the filter only.
    - psm.tsv readers without Philosopher report: they read the <group>/psm.tsv of a previous run.
    """
    plan = Plan(root)
    props = saved_workflow(plan)
    classes = {s.get("class") for s in plan.steps}
    errors = []
    flag = lambda key, default: props.get(key, default).strip().lower() == "true"  # noqa: E731
    if not flag("msfragger.run-msfragger", "true") or flag("msfragger.digest-only", "false"):
        errors.append("MSFragger does not search in this workflow (msfragger.run-msfragger=false or msfragger.digest-only=true): "
                      "FragPipe then reuses the search results of a previous run in its work directory, but msproteomics "
                      "starts from LC-MS files. Please enable the MSFragger search in the workflow.")
    if "CmdFreequant" in classes and "CmdPhilosopherFilter" not in classes:
        errors.append("FreeQuant runs without the Philosopher report panel (phi-report.run-report=false or run-validation-tab=false): "
                      "FragPipe then runs Philosopher report on a workspace with no database annotation and no FDR filter, "
                      "which msproteomics does not support. Please enable Generate reports or use IonQuant.")
    readers = sorted(classes & set(PSM_TSV_READERS))
    if readers and "CmdPhilosopherReport" not in classes:
        errors.append(f"Philosopher report does not run, so no psm.tsv is written, but these steps read it: {', '.join(readers)}. "
                      "Please enable Generate reports and Run validation tools, or disable these tools in the workflow.")
    errors += [UNSUPPORTED_CLASSES[c] for c in sorted(classes & set(UNSUPPORTED_CLASSES))]
    return errors


def tool_configs(root: ET.Element) -> tuple[dict, dict]:
    """tool_configs.json content and the legacy configs/ files (name -> text)."""
    plan = Plan(root)
    props = saved_workflow(plan)
    steps_by_tool: dict[str, list[ET.Element]] = {}
    for s in plan.steps:
        cls = s.get("class")
        if cls in INFRA_CLASSES:
            continue
        if cls not in CLASS_TO_TOOL:
            raise PlanError(f"step {s.get('order')} ({cls}, {s.get('name')}) has no extractor yet")
        steps_by_tool.setdefault(CLASS_TO_TOOL[cls], []).append(s)
    out = {tool: {"run": False, "args": "", "config_type": ct} for tool, ct in CONFIG_TYPES.items()}
    # msbooster.config always exists: the subworkflows read it as a value (empty when MSBooster does not run)
    files: dict[str, str] = {"msbooster.config": ""}
    out["database"] = database_config(props, plan)
    out["abacus"] = abacus_off_config(props)
    for tool, steps in steps_by_tool.items():
        if tool == "ionquant":
            out.update(ionquant_configs(steps, plan))
            continue
        entry = EXTRACTORS[tool](steps, plan)
        files.update(entry.pop("_files", {}))
        out[tool] = entry
    # FragPipe appends --no-ptm-score-in-pin after the LC-MS files (CmdMsfragger.java:586, :628-630), and
    # msfragger_pep_split.py takes it as an LC-MS file (:67-72), so FragPipe's MSFragger stops with
    # "Could not find MS file at <workdir>/--no-ptm-score-in-pin" (FragPipe 24.0 headless, 2026-10-11)
    if out["writesubmzml"]["run"] and out["msfragger"].get("num_slices", 1) > 1:
        raise PlanError("Write sub mzML (tab-run.write_sub_mzml=true) with a split database search (msfragger.misc.slice-db > 1): "
                        "FragPipe's MSFragger then stops with 'Could not find MS file at <workdir>/--no-ptm-score-in-pin'. "
                        "Please set msfragger.misc.slice-db to 1 or disable Write sub mzML in the Run tab.")
    # the GUI quant tab is "on" when label-free MS1 quant runs (IonQuant)
    out["quantitation"]["run"] = out["ionquant"]["run"]
    if props.get("saintexpress.run-saint-express", "false").strip().lower() == "true":
        out["saintexpress"] = saintexpress_config(props)
    return out, files


def main(argv_: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workdir", required=True, type=Path, help="FragPipe's work dir after headless --write-plan")
    ap.add_argument("--outdir", required=True, type=Path, help="writes tool_configs.json and configs/ (MSFragger, MSBooster params)")
    args = ap.parse_args(argv_)
    try:
        root = load_workdir(args.workdir)
        errors = unsupported_runs(root)
        if errors:
            raise PlanError("\n".join(errors))
        configs_, files = tool_configs(root)
    except PlanError as ex:
        print(f"plan_to_tool_configs.py: {ex}", file=sys.stderr)
        return 1
    (args.outdir / "configs").mkdir(parents=True, exist_ok=True)
    (args.outdir / "tool_configs.json").write_text(json.dumps(configs_, indent=2) + "\n", encoding="utf-8")
    for name, text in files.items():
        (args.outdir / "configs" / name).write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
