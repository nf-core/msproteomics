#!/usr/bin/env python3
"""
Parse FragPipe workflow files and generate tool-specific config files.

This script reads a FragPipe .workflow file and generates individual config files
for each tool. The output format depends on the tool type:

1. Tools WITHOUT subcommands (e.g., ionquant, msbooster, percolator):
   - Single line with raw CLI flags
   - Example: --mbr 1 --maxlfq 1 --mztol 10

2. Tools WITH subcommands (e.g., philosopher via phi-report, database, peptide-prophet):
   - One line per subcommand: "subcommand=flags"
   - Example: filter=--sequential --prot 0.01 --picked

3. MSFragger (special case):
   - Native MSFragger params format: "parameter = value"

4. Params file tools (e.g., diaumpire, speclibgen, opair, crystalc):
   - Java properties format: "key = value"
   - Example: SE.MS1PPM = 10

5. JSON output mode (--output-json):
   - Outputs a single JSON file with run flags and args for all tools (30+)
   - Format: {"tool_name": {"run": true/false, "args": "...", "config_type": "..."}, ...}
   - config_type: "cli", "params_file", or "gui_only"
   - Used by PARSEFRAGPIPEWORKFLOW module for dynamic workflow execution
   - Every tool prefix in the workflow file is handled; no silent parameter drops

Usage:
    python parse_fragpipe_workflow.py --workflow input.workflow --outdir configs/
    python parse_fragpipe_workflow.py --workflow input.workflow --tool ionquant
    python parse_fragpipe_workflow.py --workflow input.workflow --output-json --outdir configs/
    python parse_fragpipe_workflow.py --workflow input.workflow --output-json --outdir . --configs-outdir configs

The FragPipe templates fragger_closed.params and diaumpire_se.params are data files in the module's
assets/ directory; --fragger-template and --diaumpire-template name staged copies (a container only
sees the module's resources/usr/bin, not assets/).
"""

import argparse
import copy
import json
import math
import re
import shlex
import struct
import sys
from pathlib import Path
from collections import defaultdict
from decimal import ROUND_FLOOR, ROUND_HALF_EVEN, ROUND_HALF_UP, Context, Decimal
from typing import Any, Dict, List, Optional, Tuple


# Tool configuration: prefix -> (has_subcommands, subcommand_mapping)
# For tools with subcommands, the mapping shows which workflow keys map to which subcommands
TOOL_CONFIG = {
    # Tools WITHOUT subcommands - direct key=value to --key value
    "ionquant": {
        "has_subcommands": False,
        "cli_params": [
            "mbr",
            "maxlfq",
            "requantify",
            "mztol",
            "imtol",
            "rttol",
            "mbrmincorr",
            "mbrrttol",
            "mbrimtol",
            "mbrtoprun",
            "ionfdr",
            "proteinfdr",
            "peptidefdr",
            "normalization",
            "minisotopes",
            "intensitymode",
            "minscans",
            "writeindex",
            "minexps",
            "light",
            "medium",
            "heavy",
            "tp",
            "minfreq",
            "minions",
            "excludemods",
            "locprob",
            "uniqueness",
            "formula",
            "perform-ms1quant",
            "perform-isoquant",
            "ionmobility",
            "isotol",
            "isolevel",
            "isotype",
            "site-reports",
            "msstats",
        ],
    },
    "msbooster": {
        "has_subcommands": False,
        "cli_params": [
            "predict-rt",
            "predict-spectra",
            "predict-im",
            "rt-model",
            "spectra-model",
            "im-model",
            "fragmentation-type",
        ],
    },
    "percolator": {
        "has_subcommands": False,
        "cli_params": ["run-percolator", "cmd-opts", "min-prob"],
    },
    # Tools WITH subcommands (philosopher ecosystem)
    "phi-report": {
        "has_subcommands": True,
        "subcommand_key": "filter",  # The key that contains CLI flags
        "other_params": [
            "run-report",
            "print-decoys",
            "pep-level-summary",
            "prot-level-summary",
            "dont-use-prot-proph-file",
            "remove-contaminants",
        ],
    },
    "database": {
        "has_subcommands": True,
        "subcommand_key": None,  # Uses cmd-opts or specific params
        "other_params": ["decoy-tag", "db-path", "db-prefix"],
    },
    "peptide-prophet": {
        "has_subcommands": True,
        "subcommand_key": "cmd-opts",
        "other_params": ["run-peptide-prophet", "combine-pepxml"],
    },
    "protein-prophet": {
        "has_subcommands": True,
        "subcommand_key": "cmd-opts",
        "other_params": ["run-protein-prophet"],
    },
    # MSFragger - CLI flags format
    "msfragger": {
        "has_subcommands": False,
        "cli_params": [],  # All params converted to CLI flags
    },
    # Other tools
    "tmtintegrator": {
        "has_subcommands": False,
        "cli_params": [],  # Will collect all params
    },
    "ptmshepherd": {"has_subcommands": False, "cli_params": []},
    "freequant": {
        "has_subcommands": False,
        "cli_params": ["run-freequant", "mz-tol", "rt-tol"],
        "custom_generator": True,  # Uses generate_freequant_config()
    },
    "labelquant": {
        "has_subcommands": False,
        "cli_params": [
            "tol",
            "level",
            "minprob",
            "purity",
            "removelow",
            "plex",
            "brand",
        ],
        "custom_generator": True,  # Uses generate_labelquant_config()
    },
    # DIA tools
    "diaumpire": {
        "has_subcommands": False,
        "cli_params": [],  # All params collected; uses params file format
        "config_format": "params_file",  # Java properties format (key=value)
    },
    "diann": {
        "has_subcommands": False,
        "cli_params": [
            "q-value",
            "mbr",
            "library",
            "quantification-strategy",
            "quantification-strategy-2",
            "channel-normalization-strategy",
            "unrelated-runs",
            "redo-protein-inference",
            "run-specific-protein-q-value",
            "cmd-opts",
            "light",
            "medium",
            "heavy",
            "min-site-prob",
            "mod-tag",
            "generate-msstats",
            "gene-level-report",
            "protein-level-report",
            "peptide-level-report",
            "modified-peptide-level-report",
            "site-level-report",
        ],
    },
    "speclibgen": {
        "has_subcommands": False,
        "cli_params": [
            "convert-pepxml",
            "convert-psm",
            "keep-intermediate-files",
        ],
        "config_format": "structured",  # Returns dict merged into JSON entry
    },
    # Glycoproteomics
    "opair": {
        "has_subcommands": False,
        "cli_params": [
            "activation1",
            "activation2",
            "ms1_tol",
            "ms2_tol",
            "max_glycans",
            "min_isotope_error",
            "max_isotope_error",
            "filterOxonium",
            "oxonium_minimum_intensity",
            "reverse_scan_order",
            "single_scan_type",
            "glyco_db",
            "allowed_sites",
            "oxonium_filtering_file",
        ],
        "config_format": "params_file",  # O-Pair uses .NET CLI with short flags
    },
    # Interaction scoring
    "saintexpress": {
        "has_subcommands": False,
        "cli_params": [
            "max-replicates",
            "virtual-controls",
            "cmd-opts",
        ],
    },
    # Cross-linked peptide analysis
    "crystalc": {
        "has_subcommands": False,
        "cli_params": [],  # Crystal-C uses a params file, not CLI args
        "config_format": "params_file",
    },
    # GUI-only tab settings (not standalone tools, but stored in workflow file)
    "tab-run": {
        "has_subcommands": False,
        "cli_params": [
            "write_sub_mzml",
            "delete_temp_files",
            "export_matched_fragments",
            "sub_mzml_prob_threshold",
        ],
        "gui_only": True,  # Not a standalone tool; controls FragPipe GUI behavior
    },
    "quantitation": {
        "has_subcommands": False,
        "cli_params": [
            "run-label-free-quant",
        ],
        "gui_only": True,  # Not a standalone tool; controls quant panel behavior
    },
}

# Tool run flag mapping for all 28 FragPipe tools
# Maps tool name -> (workflow_prefix, run_flag_key)
TOOL_RUN_FLAGS = {
    # Core search/validation
    "msfragger": ("msfragger", "run-msfragger"),
    "msbooster": ("msbooster", "run-msbooster"),
    "percolator": ("percolator", "run-percolator"),
    "peptideprophet": ("peptide-prophet", "run-peptide-prophet"),
    "proteinprophet": ("protein-prophet", "run-protein-prophet"),
    "filter": ("phi-report", "run-report"),
    # Quantification
    "ionquant": ("ionquant", "run-ionquant"),
    "tmtintegrator": ("tmtintegrator", "run-tmtintegrator"),
    "freequant": ("freequant", "run-freequant"),
    "labelquant": (
        "labelquant",
        None,
    ),  # Run derived from tmtintegrator.intensity_extraction_tool==1
    # DIA
    "diaumpire": ("diaumpire", "run-diaumpire"),
    "diatracer": ("diatracer", "run-diatracer"),
    "diann": ("diann", "run-dia-nn"),
    "speclibgen": ("speclibgen", "run-speclibgen"),
    # PTM
    "ptmprophet": ("ptmprophet", "run-ptmprophet"),
    "ptmshepherd": ("ptmshepherd", "run-shepherd"),
    "crystalc": ("crystalc", "run-crystalc"),
    # Glyco
    "opair": ("opair", "run-opair"),
    "mbg": ("mbg", "run-mbg"),
    # Special
    "fpop": ("fpop", "fragpipe.fpop.run-fpop"),  # Nested key
    "fpop_coadaptr": ("fpop", "coadaptr.fpop.run-fpop-coadaptr"),  # Nested key
    "saintexpress": ("saintexpress", "run-saint-express"),
    "skyline": ("skyline", "run-skyline"),
    "metaproteomics": ("metaproteomics", "run-metaproteomics"),
    "transferlearning": ("transfer-learning", "run-transfer-learning"),
    # Database (always runs if present)
    "database": ("database", None),  # No run flag, always runs
    # GUI-only tab settings (included for completeness; not standalone tools)
    "tabrun": ("tab-run", None),  # GUI tab; run = write_sub_mzml or export_matched_fragments
    # Export matched fragments: FragpipeRun.java:2411 runs CmdExportMatchedFragments when
    # TabRun.isExportMatchedFragments() (TabRun.java:491-492), the checkbox named
    # TAB_PREFIX + "export_matched_fragments" = tab-run.export_matched_fragments (TabRun.java:92,223-225).
    "exportmatchedfragments": ("tab-run", "export_matched_fragments"),
    # Write sub mzML: FragpipeRun.java:2388-2391 runs CmdWriteSubMzml when TabRun.isWriteSubMzml()
    # (TabRun.java:487-488), the checkbox TAB_PREFIX + "write_sub_mzml" = tab-run.write_sub_mzml
    # (TabRun.java:220-221, default false).
    "writesubmzml": ("tab-run", "write_sub_mzml"),
    "quantitation": (
        "quantitation",
        "run-label-free-quant",
    ),  # GUI quant panel; controls IonQuant/FreeQuant
}


# FragPipe's decoy tag default when a .workflow omits `database.decoy-tag`:
# TabDatabase.java:149 initialises the "Decoy protein prefix" field from
# Bundle.properties:44 (`decoy.tag=rev_`), and Fragpipe.loadUi() only overwrites
# fields whose keys are present in the workflow (FragpipeCacheUtils.tabsLoad).
FRAGPIPE_DEFAULT_DECOY_TAG = "rev_"


def get_decoy_tag(all_params: Dict[str, Dict[str, str]]) -> str:
    """Return the single decoy tag (`database.decoy-tag`) used by every FragPipe tool.

    FragpipeRun.java:1442 reads tabDatabase.getDecoyTag() once and passes it to
    MSFragger, Percolator, PeptideProphet, Philosopher (database/filter/report),
    SpecLibGen and Metaproteomics. An empty tag aborts the run
    (FragpipeRun.java:1444-1447), so it raises here too.
    """
    decoy_tag = all_params.get("database", {}).get(
        "decoy-tag", FRAGPIPE_DEFAULT_DECOY_TAG
    )
    if not decoy_tag.strip():
        # FragpipeRun.java:1444-1447: SwingUtils.showErrorDialog(parent, "<b>Decoy protein prefix</b> is empty.
        # Please fix it in <b>Database</b> tab.", "Decoy tag is empty"); headless logs it without the HTML (SwingUtils.java:1110-1112)
        raise ValueError(
            "Decoy protein prefix is empty. Please fix it in Database tab. "
            "(database.decoy-tag is empty in the workflow file)"
        )
    return decoy_tag


def _java_properties_unescape(text: str) -> str:
    """Properties.loadConvert: \\uXXXX, \\t, \\n, \\r, \\f, and any other escaped character as itself."""
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


def _java_properties_entry(line: str):
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


def load_java_properties(workflow_path: str) -> Dict[str, str]:
    """java.util.Properties.load on the UTF-8 reader PropsFile.load uses (fragpipe/api/PropsFile.java:55-63):
    natural lines joined at an odd number of trailing backslashes, leading whitespace dropped, blank
    lines and '#' / '!' comment lines skipped, the last value of a repeated key kept."""
    props = {}
    logical = None
    with open(workflow_path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    for raw in lines:
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


def parse_workflow_file(workflow_path: str) -> Dict[str, Dict[str, str]]:
    """
    Parse a FragPipe workflow file into a nested dictionary.

    Returns:
        Dict mapping tool_prefix -> {param_name -> value}
    """
    params = defaultdict(dict)
    for key, value in load_java_properties(workflow_path).items():
        # Split by first dot to get tool prefix
        if "." in key:
            prefix, param = key.split(".", 1)
            params[prefix][param] = value
        else:
            # Top-level params (rare)
            params["_global"][key] = value

    # The MSFragger tab exists whatever the workflow names: without msfragger.* keys FragPipe
    # searches with the closed-search defaults (TabMsfragger.initMore loadDefaults(closed),
    # TabMsfragger.java:515-517; "Run MSFragger" UiCheck default true, :570).
    params.setdefault("msfragger", {})
    _resolve_radio_twins(params)
    return dict(params)


def _resolve_radio_twins(params: Dict[str, Dict[str, str]]) -> None:
    """Both run radios of a group true: FragPipe selects each in turn and the ButtonGroup keeps the
    last one (UiRadio.fromString setSelected(Boolean.valueOf(s)), utils/swing/UiRadio.java:40-42).
    loadUi applies the workflow in the entrySet order of a HashMap (Fragpipe.java:879 Collectors.toMap,
    remapped again by FragpipeCacheUtils.translateValuesToUi :184 / MapUtils.remapValues :57-60, then
    SwingUtils.valuesSet :861 iterates it), so the twin in the later bucket wins. The twins' buckets
    differ in their low 4 bits, so they never share one. The loser is set false, as FragPipe then
    shows and saves it."""
    capacity = _java_hashmap_capacity(len(_headless_ui_keys(params)))
    for first, second in RUN_FLAG_RADIO_PAIRS:
        twins = [TOOL_RUN_FLAGS[t] for t in (first, second)]
        if all(params.get(prefix, {}).get(key, "").lower() == "true" for prefix, key in twins):
            loser = min(twins, key=lambda twin: _java_hashmap_bucket(".".join(twin), capacity))
            params[loser[0]][loser[1]] = "false"


# FragPipe manifest data types (TabWorkflow.java:310).
MANIFEST_DATA_TYPES = ("DDA", "DDA+", "DIA", "DIA-Quant", "DIA-Lib", "GPF-DIA")


def _flag(all_params: Dict[str, Dict[str, str]], prefix: str, key: str, default: str = "false") -> bool:
    """A workflow checkbox value (true/false), with the panel default when the key is absent."""
    return all_params.get(prefix, {}).get(key, default).strip().lower() == "true"


def _disable(all_params: Dict[str, Dict[str, str]], tool: str) -> None:
    prefix, key = TOOL_RUN_FLAGS[tool]
    if prefix in all_params and key in all_params[prefix]:
        all_params[prefix][key] = "false"


def apply_manifest_data_types(all_params: Dict[str, Dict[str, str]], data_types: List[str]) -> None:
    """
    Apply FragPipe's manifest-driven tool adjustments to the loaded workflow, in place.

    FragPipe headless loads the workflow, then the manifest, then calls
    TabWorkflow.adjustToolsBasedOnDataTypes (TabWorkflow.java:1393, :1476-1513). A tool whose
    panel is disabled there does not run even when its run flag is true (isRun() is
    isEnabledAndChecked, e.g. UmpirePanel.java:151, CrystalcPanel.java:150):
      - DIA, DIA-Lib or GPF-DIA files: DIA-Umpire enabled; Crystal-C and PeptideProphet enabled
        only when DIA-Umpire or diaTracer runs; TMT-Integrator disabled; DIA-NN enabled only with DIA.
      - else DDA+ files: DIA-Umpire, Crystal-C, PeptideProphet, TMT-Integrator, DIA-NN disabled.
      - else: DIA-Umpire and DIA-NN disabled.
      - DIA-Quant files enable DIA-NN.
    PeptideProphet disabled selects Percolator (PepProphPanel.java:102-107); Percolator turning on
    sets the MSFragger output type to index 5, pepXML_pin (PercolatorPanel.java:101-106). The two
    are radio buttons, so Percolator only turns on (and the output type only changes) when the
    workflow selected PeptideProphet; FragPipe 24.0 runs with PeptideProphet off and a DIA
    manifest keep output_format and run no Percolator.

    @param all_params Parsed workflow (prefix -> key -> value), modified in place
    @param data_types Data types of the manifest files
    """
    unknown = sorted(set(data_types) - set(MANIFEST_DATA_TYPES))
    if unknown or not data_types:
        raise ValueError(f"Manifest data types must be a non-empty subset of {MANIFEST_DATA_TYPES}, got {sorted(set(data_types))}")
    has = set(data_types).__contains__

    if has("DIA") or has("DIA-Lib") or has("GPF-DIA"):
        dia_pseudo = get_run_flag(all_params, "diaumpire") or get_run_flag(all_params, "diatracer")
        peptideprophet_valid = dia_pseudo
        if not dia_pseudo:
            _disable(all_params, "crystalc")
        _disable(all_params, "tmtintegrator")
        if not has("DIA") and not has("DIA-Quant"):
            _disable(all_params, "diann")
    elif has("DDA+"):
        peptideprophet_valid = False
        for tool in ("diaumpire", "crystalc", "tmtintegrator"):
            _disable(all_params, tool)
        if not has("DIA-Quant"):
            _disable(all_params, "diann")
    else:
        peptideprophet_valid = True
        _disable(all_params, "diaumpire")
        if not has("DIA-Quant"):
            _disable(all_params, "diann")

    if not peptideprophet_valid and _flag(all_params, "peptide-prophet", "run-peptide-prophet"):
        all_params["peptide-prophet"]["run-peptide-prophet"] = "false"
        if not _flag(all_params, "percolator", "run-percolator"):
            all_params.setdefault("percolator", {})["run-percolator"] = "true"
            all_params.setdefault("msfragger", {})["output_format"] = "pepXML_pin"


def _java_float(value: str) -> float:
    """Float.parseFloat(value): the workflow number as a Java float (single precision)."""
    return struct.unpack("f", struct.pack("f", float(value)))[0]


def _mod_table_masses(raw: str) -> set:
    """Masses of a TabMsfragger mod table the way getVarModMassSet / getFixedModMassSet read it
    (TabMsfragger.java:1705-1725): enabled rows with non-blank sites and |mass| > 0.01.
    Workflow rows are "mass,sites,enabled,max-occurrences" separated by ";"."""
    masses = set()
    for entry in raw.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(",")
        if len(parts) < 3:
            raise ValueError(f"MSFragger mod table row '{entry}' needs mass,sites,enabled")
        mass = _java_float(parts[0].strip())
        if parts[2].strip().lower() == "true" and parts[1].strip() and abs(mass) > 0.01:
            masses.add(mass)
    return masses


def _parse_mod_masses(all_params: Dict[str, Dict[str, str]]):
    """modmasses_ionquant.txt lines FragPipe writes for every IonQuant run (CmdIonquant.java:340-350).

    The set is TabMsfragger.getVarModMassSet() + getFixedModMassSet() + getMassOffsetSet()
    (FragpipeRun.java:1852-1855), a TreeSet<Float> written with Float.toString one per line.
    FragPipe passes no mod mass set at all (no --modlist) when O-Pair runs
    (FragpipeRun.java:1885, 1990, 2017: oPairPanel.isRun() ? null : modMassSet).

    Returns:
        List of mass strings (possibly empty: FragPipe still writes the file and passes
        --modlist), or None when O-Pair runs.
    """
    if all_params.get("opair", {}).get("run-opair", "false").strip().lower() == "true":
        return None
    msfragger = msfragger_values(all_params)
    masses = _mod_table_masses(msfragger.get("table.var-mods", ""))
    masses |= _mod_table_masses(msfragger.get("table.fix-mods", ""))

    # getMassOffsetSet (TabMsfragger.java:1727-1750). Workflow use_detailed_offsets is the
    # "Use detailed mass offsets" checkbox (CONVERT_TO_GUI: Integer.parseInt(s) > 0, or a boolean).
    detailed = msfragger.get("use_detailed_offsets", "false").strip().lower()
    if detailed in ("true",) or (detailed.lstrip("-").isdigit() and int(detailed) > 0):
        for offset in msfragger.get("mass_offsets_detailed", "").split(";"):
            # MassOffset(String): mass is the first "([\d.-]+)\(" match, else 0
            match = re.search(r"([\d.-]+)\(", offset)
            if match:
                mass = _java_float(match.group(1))
                if abs(mass) > 0.01:
                    masses.add(mass)
    else:
        for token in re.split(r"[/\s]+", msfragger.get("mass_offsets", "").strip()):
            if not token:
                continue
            mass = _java_float(token)
            if abs(mass) > 0.01:
                masses.add(mass)
    return [_java_float_to_string(repr(m)) for m in sorted(masses)]


# IonQuant options FragPipe reads from the label-free quant panel, in CmdIonquant.java:224-283
# order (dynamicParams). Each is added only when its value is not blank; light/medium/heavy only
# when ionquant.use-labeling is true.
IONQUANT_DYNAMIC_PARAMS = (
    "mbr", "maxlfq", "requantify", "mztol", "imtol", "rttol", "mbrmincorr", "mbrrttol",
    "mbrimtol", "mbrtoprun", "ionfdr", "proteinfdr", "peptidefdr", "normalization",
    "minisotopes", "intensitymode", "minscans", "writeindex", "light", "medium", "heavy",
    "tp", "minfreq", "minions", "excludemods", "locprob", "uniqueness", "formula",
)

# Options CmdIonquant.java:226-254 uses when no label-free panel map is passed (the TMT passes).
IONQUANT_ISOBARIC_FIXED = (
    "--mbr 0 --maxlfq 0 --requantify 0 --mztol 10 --imtol 0.05 --rttol 1 --normalization 0 "
    "--minisotopes 1 --minscans 1 --writeindex 0 --tp 0 --minfreq 0 --minions 1 --locprob 0 "
    "--uniqueness 0"
)

# QuantPanelLabelfree widgets of the IonQuant options (QuantPanelLabelfree.java:274-354): kind,
# initial value (config form) and the spinner's DecimalFormat. quantPanelLabelfree.toMap() returns
# every widget's string, so a key the workflow omits gets the widget's initial value.
IONQUANT_PANEL_WIDGETS = {
    "mbr": ("check", "1"),
    "maxlfq": ("check", "1"),
    "requantify": ("check", "1"),
    "mztol": ("double", "10", "0.#"),
    "imtol": ("double", "0.05", "0.00"),
    "rttol": ("double", "0.4", "0.0"),
    "mbrmincorr": ("double", "0", "#.##"),
    "mbrrttol": ("double", "1", "#.##"),
    "mbrimtol": ("double", "0.05", "#.###"),
    "mbrtoprun": ("int", "10"),
    "ionfdr": ("double", "0.01", "#.###"),
    "proteinfdr": ("double", "1", "#.###"),
    "peptidefdr": ("double", "1", "#.###"),
    "normalization": ("check", "1"),
    "minisotopes": ("combo", "1", ("1", "2", "3"), ("1", "2", "3")),
    "intensitymode": ("combo", "2", ("0", "1", "2"), ("apex", "area", "auto")),
    "minscans": ("int", "3"),
    "writeindex": ("check", "0"),
    "light": ("text", ""),
    "medium": ("text", ""),
    "heavy": ("text", ""),
    "tp": ("int", "0"),
    "minfreq": ("double", "0", "#.##"),
    "minions": ("int", "1"),
    "excludemods": ("text", ""),
    "locprob": ("double", "0.75", "#.##"),
    "uniqueness": ("combo", "0", ("0", "1"), ("unique+razor", "unique only")),
    "formula": ("text", ""),
}




def ionquant_panel_value(key: str, raw) -> str:
    """Config-form string quantPanelLabelfree.toMap() gives for an IonQuant widget loaded with `raw`
    (None: the workflow omits the key)."""
    spec = IONQUANT_PANEL_WIDGETS[key]
    kind, default = spec[0], spec[1]
    value = default if raw is None else raw.strip()
    if kind == "check":
        if value.lower() in ("1", "true"):
            return "1"
        if value.lower() in ("0", "false"):
            return "0"
        raise ValueError(f"ionquant.{key}={value} is not a checkbox value (0/1/true/false)")
    if kind == "int":
        return str(int(value))
    if kind == "double":
        return java_decimal_format(spec[2], float(value))
    if kind == "combo":
        conf_values, ui_values = spec[2], spec[3]
        if value in ui_values:
            return conf_values[ui_values.index(value)]
        if value in conf_values:
            return value
        raise ValueError(f"ionquant.{key}={value} is not one of {ui_values}")
    return value

# TmtiPanel quant-level combo (MS2, MS3, ZOOM-HR); the workflow stores 2/3/4 (ui.translations:121-123).
QUANT_LEVEL_UI = {"2": "MS2", "3": "MS3", "4": "ZOOM-HR"}
ISOLEVEL_ARG = {"MS2": "2", "MS3": "3", "ZOOM-HR": "4"}


def _ionquant_head(
    ms1: bool, isobaric: bool, isotol: str, isolevel: str, isotype: str,
    site_reports: bool, msstats: bool, all_params: Dict[str, Dict[str, str]],
) -> list:
    """CmdIonquant.java:178-202: the options before --annotation and --minexps.

    --ionmobility is dataType == ImMsTimsTof. TabWorkflow.getInputDataType() returns
    ImMsTimsTof for any .d input, RegularMs for any .raw input, else the
    workflow.input.data-type.im-ms button (TabWorkflow.java:521-538); the IONQUANT module
    applies the file-type rules, the parser supplies the button value.
    """
    im_ms = all_params.get("workflow", {}).get("input.data-type.im-ms", "false")
    return [
        f"--perform-ms1quant {1 if ms1 else 0}",
        f"--perform-isoquant {1 if isobaric else 0}",
        f"--isotol {isotol}",
        f"--isolevel {ISOLEVEL_ARG[isolevel]}",
        f"--isotype {isotype}",
        f"--ionmobility {1 if im_ms.strip().lower() == 'true' else 0}",
        f"--site-reports {1 if site_reports else 0}",
        f"--msstats {1 if msstats else 0}",
    ]


# MassOffsetUtils.MassOffset(String) patterns (MassOffsetUtils.java:22-24, 210-219).
_MASS_OFFSET_MASS_RE = re.compile(r"([\d.-]+)\(")
_MASS_OFFSET_SITES_RE = re.compile(r"aa=([A-Z-nc\[\]^*]+)")
_MASS_OFFSET_SITE_CHAR_RE = re.compile(r"[A-Zalnc\[\]^*-]")
_ALL_RESIDUES = "ACDEFGHIKLMNPQRSTVWY"


def _float32(text) -> float:
    """Float.parseFloat(text) as a Python float (exact float32 value)."""
    return struct.unpack("f", struct.pack("f", float(text)))[0]


def _java_split_regex(text: str, pattern: str) -> list:
    """Java String.split(regex): trailing empty strings are removed."""
    parts = re.split(pattern, text)
    while len(parts) > 1 and parts[-1] == "":
        parts.pop()
    return parts


def _uses_detailed_offsets(msfragger: Dict[str, str]) -> bool:
    """TabMsfragger.uiCheckMassOffsetFile, bound to msfragger.use_detailed_offsets."""
    return msfragger.get("use_detailed_offsets", "false").strip().lower() == "true"


def _detailed_mass_offsets(msfragger: Dict[str, str]) -> list:
    """TabMsfragger.getDetailedOffsets() (TabMsfragger.java:1753-1763).

    Returns [(mass, residues), ...] in workflow order: mass is the float32 value parsed by
    MassOffset(String), residues is String.join("", allowedResidues).
    """
    text = msfragger.get("mass_offsets_detailed", "")
    offsets = []
    if not text:
        return offsets
    for entry in _java_split_regex(text, ";"):
        mass_match = _MASS_OFFSET_MASS_RE.search(entry)
        mass = _float32(mass_match.group(1)) if mass_match else 0.0
        sites_match = _MASS_OFFSET_SITES_RE.search(entry)
        residues = (
            "".join(_MASS_OFFSET_SITE_CHAR_RE.findall(sites_match.group(1)))
            if sites_match
            else ""
        )
        offsets.append((mass, residues))
    return offsets


def _mass_offset_set(msfragger: Dict[str, str]) -> list:
    """TabMsfragger.getMassOffsetSet() (TabMsfragger.java:1727-1751): sorted float32 masses, |m| > 0.01f."""
    threshold = _float32("0.01")
    if _uses_detailed_offsets(msfragger):
        masses = {m for m, _ in _detailed_mass_offsets(msfragger) if abs(m) > threshold}
    else:
        masses = set()
        for token in _java_split_regex(msfragger.get("mass_offsets", "").strip(), r"[/\s]+"):
            if token:
                value = _float32(token)
                if abs(value) > threshold:
                    masses.add(value)
    return sorted(masses)


def _table_mod_mass_set(raw: str) -> set:
    """TabMsfragger.getVarModMassSet()/getFixedModMassSet() (TabMsfragger.java:1705-1725).

    Table rows are "mass,sites,enabled,max"; a mod counts when enabled, with sites and |mass| > 0.01.
    """
    masses = set()
    for entry in raw.split(";"):
        parts = [part.strip() for part in entry.strip().split(",")]
        if len(parts) < 3:
            continue
        mass = _float32(parts[0])
        if parts[2].lower() == "true" and parts[1] and abs(mass) > 0.01:
            masses.add(mass)
    return masses




def generate_ionquant_config(
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]] = None,
) -> str:
    """Label-free IonQuant options exactly as FragPipe 24.0 builds them.

    FragpipeRun.java:1875-1889 calls the short CmdIonquant.configure with
    generateSiteReport = generateMSstatsInput = writeExperimentAnnotation = true, which is
    performMS1Quant = true, performIsobaricQuant = false, isoTol 20, "MS2", "tmt10"
    (CmdIonquant.java:79-100). Order: CmdIonquant.java:178-283; --threads, --annotation,
    --multidir, --filelist and --modlist are added by the IONQUANT module.
    """
    all_params = all_params or {}
    cli_parts = _ionquant_head(True, False, "20.0", "MS2", "tmt10", True, True, all_params)
    cli_parts.append("--minexps 1")  # CmdIonquant.java:216-219: always 1
    # CmdIonquant.java:296-301: labels are skipped only when the use-labeling radio reads "false"
    use_labeling = params.get("use-labeling", "false").strip() != "false"
    for key in IONQUANT_DYNAMIC_PARAMS:
        value = ionquant_panel_value(key, params.get(key))
        if not value:
            continue
        if key in ("light", "medium", "heavy") and not use_labeling:
            continue
        cli_parts.append(f"--{key} {value}")
    return " ".join(cli_parts)


def generate_percolator_config(params: Dict[str, str], run: bool) -> str:
    """Generate Percolator CLI flags from workflow params.

    percolator.cmd-opts is passed as is (CmdPercolator.java:118 setCmdLineParams, :152
    addFreeCommandLineParams); num-threads and the result file names are set by the module.

    @param params percolator.* workflow params
    @param run get_run_flag(all_params, "percolator"): an absent run-percolator key is FragPipe's
               radio default (selected), so the key alone does not decide
    """
    cli_parts = []

    if not run:
        return "# percolator disabled"

    # Add cmd-opts if present (contains raw CLI flags)
    if "cmd-opts" in params and params["cmd-opts"]:
        cli_parts.append(params["cmd-opts"])

    # M13: min-prob is already emitted as a separate JSON field (see generate_tool_config)
    # and passed via task.ext.min_prob in Nextflow config. Do NOT embed in CLI string
    # as "# min-prob=X" because '#' is treated as a literal argument by percolator,
    # causing it to use '#' as the decoy prefix pattern.

    return " ".join(cli_parts) if cli_parts else ""


def generate_freequant_config(params: Dict[str, str]) -> str:
    """Philosopher freequant options: QuantPanelLabelfree.getFreequantOptsAsText()
    (QuantPanelLabelfree.java:154-161) " --ptw <freequant.rt-tol> --tol <freequant.mz-tol>",
    used by both the label-free and the TMT FreeQuant (FragpipeRun.java:1817, 2052).
    --raw and --dir are added by the FREEQUANT module (CmdFreequant.java:116-125).
    toMap() holds every widget, so an absent key is the spinner's initial value, written in the
    spinner's format: rt-tol builder(0.4, ...) "0.0#", mz-tol builder(10.0, ...) "0.#"
    (QuantPanelLabelfree.java:238-243).
    """
    rt_tol = java_decimal_format("0.0#", float(params.get("rt-tol", "0.4")))
    mz_tol = java_decimal_format("0.#", float(params.get("mz-tol", "10")))
    return f"--ptw {rt_tol} --tol {mz_tol}"


# QuantLabel.LABELS (QuantLabel.java:51-70): label name -> (type, number of reagents), and
# Joiner.on(",").join(labelModMap.get(name)) (QuantLabel.java:72-87, Float.toString, HashMultimap
# order) for the TMT-Integrator label_masses; "0" for labels without masses (TmtiPanel.java:1208).
QUANT_LABELS = {
    "TMT-0": ("TMT", 1, "224.15248,295.1896"),
    "TMT-2": ("TMT", 2, "225.15584"),
    "TMT-6": ("TMT", 6, "229.16293"),
    "TMT-10": ("TMT", 10, "229.16293"),
    "TMT-11": ("TMT", 11, "229.16293"),
    "TMT-16": ("TMT", 16, "304.2071"),
    "TMT-18": ("TMT", 18, "304.2071"),
    "TMT-35": ("TMT", 35, "304.2071"),
    "iodoTMT-6": ("TMT", 6, "329.2266"),
    "iTRAQ-4": ("iTRAQ", 4, "144.10207"),
    "iTRAQ-8": ("iTRAQ", 8, "304.20535"),
    "sCLIP-6": ("sCLIP", 6, "481.2"),
    "IBT-16": ("IBT", 16, "227.1"),
    "DiLeu-12": ("DiLeu", 12, "0"),
    "DiLeu-1": ("DiLeu", 1, "0"),
    "DeAla-13": ("DeAla", 13, "0"),
    "Trp-2": ("Trp", 2, "875.3903"),
    "develop": ("xtag", 32, "0"),
    "develop2": ("xtag2", 29, "0"),
}

# ui.translations:83-101, tmtintegrator.channel_num config value -> UI label (first match wins,
# FragpipeCacheUtils.java:185-203)
CHANNEL_NUM_UI = {}
for _ui, _conf in (
    ("TMT-0", "0"), ("TMT-2", "2"), ("TMT-6", "6"), ("TMT-10", "10"), ("TMT-11", "11"),
    ("TMT-16", "16"), ("TMT-18", "18"), ("TMT-35", "35"), ("iTRAQ-4", "4"), ("iTRAQ-8", "8"),
    ("sCLIP-6", "6"), ("IBT-16", "16"), ("DiLeu-12", "12"), ("DiLeu-1", "1"), ("DeAla-13", "13"),
    ("develop", "32"), ("develop2", "29"), ("trp-2", "2"), ("iodoTMT-6", "6"),
):
    CHANNEL_NUM_UI.setdefault(_conf, _ui)


def tmt_quant_label(tmti_params: Dict[str, str]) -> str:
    """TmtiPanel.getSelectedLabel() name for tmtintegrator.channel_num (TmtiPanel.java:665-673).

    Without a value the combo keeps its first item, TMT-0 (TmtiPanel.java:323-324)."""
    value = tmti_params.get("channel_num", "").strip()
    if not value:
        return "TMT-0"
    for name in QUANT_LABELS:
        if name.lower() == value.lower():
            return name
    if value.lower() in CHANNEL_NUM_UI:
        return CHANNEL_NUM_UI[value.lower()]
    raise ValueError(f"Unknown label name in combo box: tmtintegrator.channel_num={value}")


def tmt_quant_level(tmti_params: Dict[str, str]) -> str:
    """TmtiPanel.getQuantLevel(): MS2, MS3 or ZOOM-HR (workflow 2/3/4, ui.translations:121-123)."""
    value = tmti_params.get("quant_level", "MS2").strip()
    level = QUANT_LEVEL_UI.get(value, value)
    if level not in ISOLEVEL_ARG:
        raise ValueError(f"Unknown tmtintegrator.quant_level={value}: FragPipe offers MS2, MS3, ZOOM-HR")
    return level


def _java_percent_3f(value: str) -> str:
    """String.format("%.3f", Double.parseDouble(value)) (HALF_UP on the shortest decimal)."""
    return str(Decimal(repr(float(value))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def generate_labelquant_config(
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]] = None,
) -> str:
    """Philosopher labelquant options before --annot (CmdLabelquant.java:98-116) from the
    TMT-Integrator panel (FragpipeRun.java:1962-1967):
    --tol <tolerance> --level <2|3> --minprob/--purity/--removelow <%.3f> --plex <reagents>.
    ZOOM-HR is rejected as in CmdLabelquant.java:105-107. --annot, --brand
    (labelquant_brand), --raw and --dir are added by the PHILOSOPHERLABELQUANT module.
    """
    tmti = (all_params or {}).get("tmtintegrator", {})
    level = tmt_quant_level(tmti)
    if level == "ZOOM-HR":
        raise ValueError(
            "ZOOM-HR is not supported for PhilosopherLabelquant. Please switch the "
            "`Intensity extraction tool` to `IonQuant`."
        )
    label = tmt_quant_label(tmti)
    return " ".join([
        f"--tol {int(tmti.get('tolerance', '20'))}",
        f"--level {ISOLEVEL_ARG[level]}",
        f"--minprob {_java_percent_3f(tmti.get('min_pep_prob', '0.9'))}",
        f"--purity {_java_percent_3f(tmti.get('min_purity', '0.5'))}",
        f"--removelow {_java_percent_3f(tmti.get('min_percent', '0.05'))}",
        f"--plex {QUANT_LABELS[label][1]}",
    ])


def labelquant_brand(all_params: Dict[str, Dict[str, str]]) -> str:
    """labelquant --brand: label.getType().toLowerCase() (CmdLabelquant.java:144-145)."""
    return QUANT_LABELS[tmt_quant_label(all_params.get("tmtintegrator", {}))][0].lower()


# MSBoosterPanel.modelMap: model display name (workflow / GUI value) -> MSBooster model key
MSBOOSTER_MODEL_MAP = {
    "DIA-NN": "DIA-NN",
    "DeepLC HeLa HF": "Deeplc_hela_hf",
    "AlphaPept RT Generic": "AlphaPept_rt_generic",
    "Prosit 2019 iRT": "Prosit_2019_irt",
    "Prosit 2020 iRT TMT": "Prosit_2020_irt_TMT",
    "MS2PIP 2021 HCD": "ms2pip_2021_hcd",
    "MS2PIP timsTOF 2024": "ms2pip_timsTOF2024",
    "MS2PIP TTOF 5600": "ms2pip_ttof5600",
    "MS2PIP Immuno HCD": "ms2pip_immuno_hcd",
    "MS2PIP CID TMT": "ms2pip_CID_TMT",
    "MS2PIP iTRAQ Phospho": "ms2pip_iTRAQphospho",
    "AlphaPept MS2 Generic": "AlphaPept_ms2_generic",
    "Prosit 2019 Intensity": "Prosit_2019_intensity",
    "Prosit 2023 Intensity timsTOF": "Prosit_2023_intensity_timsTOF",
    "Prosit 2020 Intensity CID": "Prosit_2020_intensity_CID",
    "Prosit 2020 Intensity TMT": "Prosit_2020_intensity_TMT",
    "Prosit 2020 Intensity HCD": "Prosit_2020_intensity_HCD",
    "AlphaPept CCS Generic": "AlphaPept_ccs_generic",
    "IM2Deep": "IM2Deep",
    "UniSpec": "UniSpec",
    "Prosit 2024 Intensity CIT": "Prosit_2024_intensity_cit",
    "Prosit 2025 Intensity MultiFrag": "Prosit_2025_intensity_MultiFrag",
    "Prosit 2024 iRT CIT": "Prosit_2024_irt_cit",
}

# ui.translations:146-170 (msbooster.<key>::<UI value> <=> <config value>): a saved workflow holds the
# config value; WorkflowTranslator turns it back into the UI value, which MSBoosterPanel.modelMap maps
# to the MSBooster key. A value without a translation is taken as the UI value (the combo box).
MSBOOSTER_CONFIG_NAMES = {
    "rt-model": {
        "Deeplc_hela_hf": "DeepLC HeLa HF",
        "AlphaPept_rt_generic": "AlphaPept RT Generic",
        "Prosit_2019_irt": "Prosit 2019 iRT",
        "Prosit_2020_irt_TMT": "Prosit 2020 iRT TMT",
        "Prosit_2024_irt_cit": "Prosit 2024 iRT CIT",
    },
    "spectra-model": {
        "ms2pip_2021_hcd": "MS2PIP 2021 HCD",
        "ms2pip_timsTOF2024": "MS2PIP timsTOF 2024",
        "ms2pip_ttof5600": "MS2PIP TTOF 5600",
        "ms2pip_Immuno_HCD": "MS2PIP Immuno HCD",
        "ms2pip_CID_TMT": "MS2PIP CID TMT",
        "ms2pip_iTRAQphospho": "MS2PIP iTRAQ Phospho",
        "AlphaPept_ms2_generic": "AlphaPept MS2 Generic",
        "Prosit_2019_intensity": "Prosit 2019 Intensity",
        "Prosit_2023_intensity_timsTOF": "Prosit 2023 Intensity timsTOF",
        "Prosit_2020_intensity_CID": "Prosit 2020 Intensity CID",
        "Prosit_2020_intensity_TMT": "Prosit 2020 Intensity TMT",
        "Prosit_2020_intensity_HCD": "Prosit 2020 Intensity HCD",
        "Prosit_2024_intensity_cit": "Prosit 2024 Intensity CIT",
        "Prosit_2025_intensity_MultiFrag": "Prosit 2025 Intensity MultiFrag",
    },
    "im-model": {
        "AlphaPept_ccs_generic": "AlphaPept CCS Generic",
    },
}


# MSBoosterPanel checkboxes: the state of a key the workflow does not set (MSBoosterPanel.java:142-148 predict-* on;
# the find-best checkboxes off, :296, :324, :362)
MSBOOSTER_FLAG_DEFAULTS = {"predict-rt": "true", "predict-spectra": "true", "predict-im": "true"}


def msbooster_flag(params: Dict[str, str], key: str) -> str:
    """"true" or "false" for an msbooster checkbox key, with the panel's state for an absent key."""
    return "true" if params.get(key, MSBOOSTER_FLAG_DEFAULTS.get(key, "false")).strip().lower() == "true" else "false"


def msbooster_model_ui_name(params: Dict[str, str], key: str) -> str:
    """The MSBoosterPanel combo value of msbooster.<key> (rt-model, spectra-model, im-model); default DIA-NN."""
    value = params.get(key, "").strip() or "DIA-NN"
    return MSBOOSTER_CONFIG_NAMES[key].get(value, value)


# MSBoosterPanel.create{Rt,Spectra,Im}ModelsPanel: the models tested by "Find best model"
MSBOOSTER_TEST_MODELS = {
    "table.rt-models": [
        "DIA-NN", "DeepLC HeLa HF", "AlphaPept RT Generic", "Prosit 2019 iRT",
        "Prosit 2020 iRT TMT", "Prosit 2024 iRT CIT",
    ],
    "table.spectra-models": [
        "DIA-NN", "MS2PIP 2021 HCD", "MS2PIP timsTOF 2024", "MS2PIP TTOF 5600",
        "MS2PIP Immuno HCD", "MS2PIP CID TMT", "MS2PIP iTRAQ Phospho", "AlphaPept MS2 Generic",
        "Prosit 2019 Intensity", "Prosit 2023 Intensity timsTOF", "Prosit 2020 Intensity CID",
        "Prosit 2020 Intensity TMT", "Prosit 2020 Intensity HCD", "UniSpec",
        "Prosit 2024 Intensity CIT", "Prosit 2025 Intensity MultiFrag",
    ],
    "table.im-models": ["DIA-NN", "AlphaPept CCS Generic", "IM2Deep"],
}

# Workflow msbooster.fragmentation-type index -> MSBoosterPanel combo item (ui.translations)
MSBOOSTER_FRAGMENTATION_TYPES = {
    "0": "auto", "1": "HCD", "2": "CID", "3": "ETD", "4": "UVPD", "5": "ECD", "6": "EID", "7": "ETCID",
}


def _msbooster_test_models(params: Dict[str, str], table_key: str) -> str:
    """MSBoosterPanel.{rt,spectra,im}TestModels(): the model keys of the table, comma-joined.

    A workflow table value is "name,enabled;name,enabled;..." (ModelsTable.asString). Loading it
    re-enables every row (ModelsTable.fromString builds the rows with MSBoosterPanel.f, which sets
    every Enabled cell to true), so every listed model is tested. Without the key the panel's
    default rows apply, all enabled.
    """
    raw = params.get(table_key, "")
    names = [chunk.split(",")[0] for chunk in raw.split(";") if len(chunk.split(",")) == 2] if raw else MSBOOSTER_TEST_MODELS[table_key]
    unknown = [name for name in names if name not in MSBOOSTER_MODEL_MAP]
    if unknown:
        raise ValueError(f"msbooster.{table_key} lists models FragPipe 24.0 does not know: {unknown}")
    return ",".join(MSBOOSTER_MODEL_MAP[name] for name in names)


# Stands for the staged spectral library in the msbooster params text (msboosterParamsText resolves it)
MSBOOSTER_LIBRARY_TOKEN = "{SPECTRAL_LIBRARY}"


def generate_msbooster_config(
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]],
    run: bool,
) -> str:
    """Workflow-derived lines of FragPipe's msbooster_params.txt (CmdMSBooster.configure).

    Keys and values as CmdMSBooster writes them. numThreads is the run's thread count
    (CmdMSBooster.java:165, TabWorkflow.getThreads() = workflow.threads, which the headless
    --threads overrides, Fragpipe.java:508-514). MSBooster's loess RT fit gives different
    delta_RT_loess_real values with another thread count, so a workflow that runs MSBooster
    without a thread count is an error (FragPipe would take the host's CPU count).
    The MSBOOSTER module adds the runtime keys
    (DiaNN, fragger, unimodObo, mzmlDirectory, pinPepXMLDirectory) and ANDs
    useIM and findBestImModel with "the run has timsTOF data" (CmdMSBooster.java:151-157,
    170, 178), which only the input files tell.

    @param run get_run_flag(all_params, "msbooster"): the run checkbox AND the Validation tab's
               "Run validation tools" (TabValidation.java:58-59, MSBoosterPanel.java:291)
    """
    def flag(key: str) -> str:
        return msbooster_flag(params, key)

    def model(key: str) -> str:
        value = msbooster_model_ui_name(params, key)
        if value not in MSBOOSTER_MODEL_MAP:
            raise ValueError(f"msbooster.{key}={value} is not an MSBooster model of FragPipe 24.0 (MSBoosterPanel.modelMap)")
        return MSBOOSTER_MODEL_MAP[value]

    threads = all_params.get("workflow", {}).get("threads", "").strip()
    if run and not (threads.isdigit() and int(threads) > 0):
        raise ValueError(
            "MSBooster needs the thread count FragPipe runs with: set workflow.threads in the workflow "
            "file or run the pipeline with --fragpipe_threads <n> (the parser's --threads), "
            f"got {threads!r}; its RT loess output depends on it"
        )
    lines = [
        "useDetect = false",
        "renamePin = 1",
        f"useRT = {flag('predict-rt')}",
        f"useSpectra = {flag('predict-spectra')}",
        f"useIM = {flag('predict-im')}",
        "deletePreds = false",
        f"rtModel = {model('rt-model')}",
        f"spectraModel = {model('spectra-model')}",
        f"imModel = {model('im-model')}",
        f"findBestRtModel = {flag('find-best-rt-model')}",
        f"findBestSpectraModel = {flag('find-best-spectra-model')}",
        f"findBestImModel = {flag('find-best-im-model')}",
        f"KoinaURL = {params.get('koina-url', '').strip()}",
    ]
    if threads:
        lines.append(f"numThreads = {threads}")
    library = params.get("spectral-library-path", "").strip()
    if library:
        # CmdMSBooster.java:180-184: "if (!libraryPath.isEmpty()) { write spectraPredFile/RTPredFile/IMPredFile = libraryPath }".
        # The file is staged by MSBOOSTER (tool_configs msbooster.spectral_library), so the text holds a token
        # that msboosterParamsText replaces with the staged path.
        lines += [f"{key} = {MSBOOSTER_LIBRARY_TOKEN}" for key in ("spectraPredFile", "RTPredFile", "IMPredFile")]
    lines += [
        f"rtSearchModelsString = {_msbooster_test_models(params, 'table.rt-models')}",
        f"ms2SearchModelsString = {_msbooster_test_models(params, 'table.spectra-models')}",
        f"imSearchModelsString = {_msbooster_test_models(params, 'table.im-models')}",
        f"FragmentationType = {MSBOOSTER_FRAGMENTATION_TYPES.get(params.get('fragmentation-type', '0').strip(), params.get('fragmentation-type', '').strip())}",
    ]
    return "\n".join(lines)


def peptideprophet_cmd_opts(all_params: Dict[str, Dict[str, str]]) -> str:
    """PeptideProphet options as CmdPeptideProphet.addFreeCommandLineParams builds them.

    CmdPeptideProphet.java:375-389: the cmd-opts tokens (dropping a literal
    "peptideprophet"), then CmdPeptideProphet.java:391-404:
    - enzyme name == "nonspecific" or a non-blank second-enzyme cut: add --nontt, --nonmc
    - else enzyme name "custom" or "nocleavage": add --nontt, --nonmc, and
      "--enzyme nonspecific" unless --enzyme is already present
    Each flag is added only when absent (addToListIfNotThere). The enzyme name and cut
    are the MSFragger tab texts (TabMsfragger.java:2143-2148), compared verbatim.
    """
    pp = all_params.get("peptide-prophet", {})
    tokens = [t for t in pp.get("cmd-opts", "").split() if t != "peptideprophet"]
    msfragger = msfragger_values(all_params)
    enzyme_name = msfragger["search_enzyme_name_1"]
    enzyme_cut_2 = msfragger["search_enzyme_cut_2"]

    def add_if_absent(opt: str) -> bool:
        if opt in tokens:
            return False
        tokens.append(opt)
        return True

    if enzyme_name == "nonspecific" or enzyme_cut_2.strip():
        add_if_absent("--nontt")
        add_if_absent("--nonmc")
    elif enzyme_name in ("custom", "nocleavage"):
        add_if_absent("--nontt")
        add_if_absent("--nonmc")
        if add_if_absent("--enzyme"):
            tokens.append("nonspecific")
    return " ".join(tokens)


def peptideprophet_combine(all_params: Dict[str, Dict[str, str]]) -> bool:
    """Per-group combined PeptideProphet (CmdPeptideProphet.java:190-202).

    The combine-pepxml checkbox, or "--combine" anywhere in the cmd-opts text
    (case-insensitive), switches to one PeptideProphet run per group.
    """
    pp = all_params.get("peptide-prophet", {})
    return (
        pp.get("combine-pepxml", "false").strip().lower() == "true"
        or "--combine" in pp.get("cmd-opts", "").lower()
    )




PERCOLATOR_DEFAULT_MIN_PROB = "0.50"  # PercolatorPanel.java:140 spinnerDouble(0.50, 0, 1, 0.01)


def percolator_min_prob(all_params: Dict[str, Dict[str, str]]) -> str:
    """The min-prob argument of PercolatorOutputToPepXML (CmdPercolator.java:298).

    FragpipeRun.java:1642 passes percolatorPanel.getMinProb() =
    Double.parseDouble(spinMinProb.asString()) (PercolatorPanel.java:77-79); the
    spinner shows the value with format "#.##" (PercolatorPanel.java:140, HALF_EVEN),
    and CmdPercolator appends it as `minProb + ""` (Double.toString).
    """
    text = all_params.get("percolator", {}).get("min-prob", PERCOLATOR_DEFAULT_MIN_PROB)
    # DecimalFormat rounds the binary double value (HALF_EVEN on its exact expansion)
    rounded = Decimal(float(text.strip())).quantize(Decimal("0.01"), rounding=ROUND_HALF_EVEN)
    return java_double_to_string(float(rounded))



def label_free_ionquant_runs(all_params: Dict[str, Dict[str, str]]) -> bool:
    """QuantPanelLabelfree.isRunIonQuant(): isRun() && IonQuant radio (QuantPanelLabelfree.java:110-116).

    isRun() is the quantitation.run-label-free-quant checkbox (QuantPanelLabelfree.java:187-188);
    the radio button starts as FragPipe's when its key is absent (get_run_flag).
    """
    return get_run_flag(all_params, "ionquant")


REPORT_DEFAULT_FILTER = "--sequential --prot 0.01"  # ReportPanel.java:114


def tmti_extraction_tool(all_params: Dict[str, Dict[str, str]]) -> int:
    """TmtiPanel.getIntensityExtractionTool(): combo index (TmtiPanel.java:259-263, :836-838).

    0 = IonQuant, 1 = Philosopher, 2 = skip extraction. The workflow stores the item text.
    """
    items = ["IonQuant", "Philosopher", "Skip extraction. Run TMT-Integrator only"]
    value = all_params.get("tmtintegrator", {}).get("extraction_tool", "IonQuant").strip()
    if value in items:
        return items.index(value)
    if value in ("0", "1", "2"):
        return int(value)
    raise ValueError(f"tmtintegrator.extraction_tool={value!r} is not one of {items}")


def philosopher_generate_msstats(all_params: Dict[str, Dict[str, str]]) -> bool:
    """tmtiPanel.isRun() && tmtiPanel.isMsstats() (FragpipeRun.java:1768), with its error check.

    FragpipeRun.java:1770-1773 refuses the run when Philosopher MSstats output is requested
    while IonQuant is the intensity extraction tool.
    """
    msstats = _flag(all_params, "tmtintegrator", "run-tmtintegrator", "false") and _flag(
        all_params, "tmtintegrator", "philosopher-msstats", "false"
    )
    if msstats and tmti_extraction_tool(all_params) == 0:
        raise ValueError(
            "Generate MSstats files (using Philosopher) was enabled but Philosopher was not "
            "selected as Intensity extraction tool (FragpipeRun.java:1770-1773)"
        )
    return msstats


def fragpipe_filter_options(text: str, run_ionquant: bool, dont_use_protxml: bool) -> str:
    """User filter options as CmdPhilosopherFilter.configure rewrites them (CmdPhilosopherFilter.java:82-97).

    --razor is removed (the module appends it after --protxml), --prot becomes 0.01 when
    label-free IonQuant runs, and without protXML --sequential and --prot are dropped.
    """
    if not text.strip():
        return ""
    text = text.replace("--razor", "")
    if run_ionquant:
        text = re.sub(r"--prot\s+[^\s]+", "--prot 0.01", text)
        if "--prot" not in text:
            text = text + " --prot 0.01"
    if dont_use_protxml:
        text = text.replace("--sequential", "")
        text = re.sub(r"--prot\s+[^\s]+", "", text)
    return " ".join(text.split())


def philosopher_filter_entry(all_params: Dict[str, Dict[str, str]]) -> Dict[str, Any]:
    """Filter and report settings of FragpipeRun.java:1702-1780.

    dont_use_protxml: true when ProteinProphet does not run (FragpipeRun.java:1716-1719; a fresh
    pipeline run has no earlier prot.xml), else the "Do not use ProteinProphet file" checkbox.
    Report flags follow CmdPhilosopherReport.java:54-62; --msstats only for a single group
    (isMultiExpReport = more than one group), so both variants are emitted.
    """
    report = all_params.get("phi-report", {})
    dont_use = (not get_run_flag(all_params, "proteinprophet")) or _flag(
        all_params, "phi-report", "dont-use-prot-proph-file", "false"
    )
    options = fragpipe_filter_options(
        report.get("filter", REPORT_DEFAULT_FILTER), label_free_ionquant_runs(all_params), dont_use
    )
    decoys = ["--decoys"] if _flag(all_params, "phi-report", "print-decoys", "false") else []
    msstats = ["--msstats"] if philosopher_generate_msstats(all_params) else []
    removecontam = ["--removecontam"] if _flag(all_params, "phi-report", "remove-contaminants", "false") else []
    return {
        "args": options,
        "dont_use_protxml": dont_use,
        "report_args": " ".join(decoys + removecontam),
        "report_args_single_group": " ".join(decoys + msstats + removecontam),
    }


def abacus_iprophet_entries(
    all_params: Dict[str, Dict[str, str]], result: Dict[str, Dict[str, Any]]
) -> Dict[str, Dict[str, Any]]:
    """Philosopher abacus and iProphet run decisions and options (FragpipeRun.java:1782-1809).

    The "more than one group" condition is known only at run time; the pipeline adds it.
    Abacus options follow CmdPhilosopherAbacus.java:106-130; "--labels --plex N" is appended
    by the pipeline from labels_channel_num (plex = number of reagents of the TMT label).
    iProphet options follow CmdIprophet.java:44-52; --threads and the pepXML list are runtime.
    """
    decoy_tag = get_decoy_tag(all_params)
    report_runs = result.get("filter", {}).get("run", False)
    lfq_ionquant = label_free_ionquant_runs(all_params)
    tmti_runs = _flag(all_params, "tmtintegrator", "run-tmtintegrator", "false")
    extraction_tool = tmti_extraction_tool(all_params) if tmti_runs else -1
    phi_msstats = philosopher_generate_msstats(all_params)
    no_protxml = _flag(all_params, "phi-report", "dont-use-prot-proph-file", "false")
    prot_summary = _flag(all_params, "phi-report", "prot-level-summary", "true")
    pep_summary = _flag(all_params, "phi-report", "pep-level-summary", "false")

    abacus_runs = (
        report_runs
        and not lfq_ionquant
        and not (tmti_runs and extraction_tool == 0)
        and (phi_msstats or (not no_protxml and prot_summary) or pep_summary)
    )
    filter_text = all_params.get("phi-report", {}).get("filter", REPORT_DEFAULT_FILTER)
    addons = []
    for token in filter_text.split():
        if token in ("--picked", "--uniqueonly") and token not in addons:
            addons.append(token)
    for token in ("--razor", "--reprint"):
        if token not in addons:
            addons.append(token)
    abacus_args = addons + ["--tag", decoy_tag]
    if not no_protxml and prot_summary:
        abacus_args.append("--protein")
    if pep_summary:
        abacus_args.append("--peptide")
    # args are emitted even when FragPipe's own conditions do not run Abacus: the pipeline also runs
    # it to produce SAINTexpress's reprint tables
    abacus = {"run": abacus_runs, "args": " ".join(abacus_args), "config_type": "cli"}
    if tmti_runs and extraction_tool == 1:
        abacus["labels_channel_num"] = tmt_quant_label(all_params["tmtintegrator"])

    iprophet_runs = abacus_runs and not lfq_ionquant and pep_summary
    iprophet = {
        "run": iprophet_runs,
        "args": f"--decoy {decoy_tag} --nonsp --output combined" if iprophet_runs else "",
        "config_type": "cli",
    }
    return {"abacus": abacus, "iprophet": iprophet}


def generate_philosopher_config(
    params: Dict[str, str],
    subcommand_prefix: str,
    all_params: Dict[str, Dict[str, str]] = None,
) -> str:
    """
    Generate Philosopher config with subcommands.

    Format: subcommand=flags

    The phi-report filter line carries the same flags as tool_configs.json
    (_get_philosopher_cli_args, CmdPhilosopherFilter.java:82-115).
    """
    lines = []

    if subcommand_prefix == "phi-report":
        # Filter subcommand: same flags as tool_configs.json (CmdPhilosopherFilter.java:82-115)
        if "filter" in params:
            lines.append(
                f"filter={_get_philosopher_cli_args(params, 'phi-report', all_params)}"
            )

        # Report options (converted to flags)
        report_flags = []
        if params.get("print-decoys", "false").lower() == "true":
            report_flags.append("--decoys")
        if params.get("remove-contaminants", "false").lower() == "true":
            report_flags.append("--removecontam")
        # NOTE: --msstats is intentionally NOT added here. FragPipe gates it behind
        # `!isMultiExpReport` (CmdPhilosopherReport.java:54-59). Our pipeline always
        # runs with multiple samples (multiple LCMS file groups), so isMultiExpReport
        # is always true and FragPipe never adds --msstats to philosopher report.
        if report_flags:
            lines.append(f"report={' '.join(report_flags)}")
        else:
            lines.append("report=")

        # Add metadata
        lines.append(f"# run-report={params.get('run-report', 'true')}")
        lines.append(f"# pep-level-summary={params.get('pep-level-summary', 'false')}")
        lines.append(f"# prot-level-summary={params.get('prot-level-summary', 'true')}")
        lines.append(
            f"# dont-use-prot-proph-file={params.get('dont-use-prot-proph-file', 'false')}"
        )

    elif subcommand_prefix == "peptide-prophet":
        if params.get("run-peptide-prophet", "false").lower() == "true":
            cmd_opts = peptideprophet_cmd_opts(all_params)
            lines.append(f"peptideprophet={cmd_opts}")
        lines.append(
            f"# run-peptide-prophet={params.get('run-peptide-prophet', 'false')}"
        )
        lines.append(f"# combine-pepxml={params.get('combine-pepxml', 'false')}")

    elif subcommand_prefix == "protein-prophet":
        if params.get("run-protein-prophet", "false").lower() == "true":
            cmd_opts = params.get("cmd-opts", "")
            lines.append(f"proteinprophet={cmd_opts}")
        lines.append(
            f"# run-protein-prophet={params.get('run-protein-prophet', 'true')}"
        )

    elif subcommand_prefix == "database":
        # Database command options
        decoy_tag = get_decoy_tag(all_params)
        lines.append(f"database=--prefix {decoy_tag}")
        lines.append(f"# decoy-tag={decoy_tag}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------------------------
# MSFragger params files, written the way FragPipe 24.0 writes them.
#
# FP = FragPipe-GUI/src/main/java/org/nesvilab/fragpipe/ at tag 24.0.
# FragPipe builds fragger.params in four stages, reproduced here in the same order:
#   1. TabMsfragger.initMore() loads the closed-search defaults into the form
#      (FP/tabs/TabMsfragger.java:513-514 -> loadDefaults -> formFrom, :1639-1650, :2088-2096).
#   2. Loading the workflow sets every msfragger.* widget the file has
#      (FP/Fragpipe.java:877-882 loadUi -> FP/api/FragpipeCacheUtils.java tabsLoad/translateValuesToUi
#      -> utils/SwingUtils.java:859-869 valuesSet). Widgets the file does not have keep stage-1 values.
#   3. TabMsfragger.getParams() reads the form back into a fresh MsfraggerParams
#      (FP/tabs/TabMsfragger.java:1652-1674 formCollect, :1780-1857 paramsFromMap).
#   4. CmdMsfragger.configure() sets database/decoy/mzBIN, derives the DIA, GPF-DIA and DDA+
#      params and saves each with Props.writeProps (FP/cmd/CmdMsfragger.java:441-468, :710-741,
#      FP/params/Props.java:240-273).
# ---------------------------------------------------------------------------------------------

# FP/tools/fragger/fragger_closed.params, copied verbatim (MsfraggerParams.DEFAULT_FILE_CLOSEDSEARCH), in the
# module's assets/; main() sets it from --fragger-template. The script's own location says nothing about
# assets/ (in a Wave container the module's resources/usr/bin is /usr/bin).
FRAGGER_TEMPLATE: Optional[Path] = None

# MSFragger params file written for each FragPipe data type (FP/cmd/CmdMsfragger.java:444-447).
# The DDA file is saved as fragger_dda.params when the run also has DIA/GPF-DIA/DIA-Lib files
# (msfragger_params_file_name); the content is the same.
MSFRAGGER_PARAMS_FILES = {
    "DDA": "fragger.params",
    "DIA": "fragger_dia.params",
    "GPF-DIA": "fragger_gpfdia.params",
    "DDA+": "fragger_dda_plus.params",
}

# MsfraggerParams.createComments() (FP/tools/fragger/MsfraggerParams.java:330-425). Props.writeProps
# appends "\t\t\t# <comment>" when the name has a non-blank comment here (FP/params/Props.java:259-264).
FRAGGER_COMMENTS = {
    "num_threads": "Number of CPU threads to use.",
    "database_name": "Path to the protein database file in FASTA format.",
    "precursor_mass_lower": "Lower bound of the precursor mass window.",
    "precursor_mass_upper": "Upper bound of the precursor mass window.",
    "precursor_mass_units": "Precursor mass tolerance units (0 for Da, 1 for ppm).",
    "precursor_true_tolerance": "True precursor mass tolerance (window is +/- this value).",
    "precursor_true_units": "True precursor mass tolerance units (0 for Da, 1 for ppm).",
    "fragment_mass_tolerance": "Fragment mass tolerance (window is +/- this value).",
    "fragment_mass_units": "Fragment mass tolerance units (0 for Da, 1 for ppm).",
    "data_type": "Data type (0 for DDA, 1 for DIA, 2 for gas-phase fractionation DIA, 3 for DDA+).",
    "calibrate_mass": "Perform mass calibration (0 for OFF, 1 for ON, 2 for ON and find optimal parameters, 4 for ON and find the optimal fragment mass tolerance).",
    "use_all_mods_in_first_search": "Use all variable modifications in first search (0 for No, 1 for Yes).",
    "write_calibrated_mzml": "Write calibrated MS2 scan to a mzML file (0 for No, 1 for Yes).",
    "write_uncalibrated_mzml": "Write uncalibrated MS2 scan to a mzML file (0 for No, 1 for Yes). Only for .raw and .d formats.",
    "decoy_prefix": "Prefix of the decoy protein entries. Used for parameter optimization only.",
    "keep_decoys": "",
    "isotope_error": "Also search for MS/MS events triggered on specified isotopic peaks.",
    "mass_offsets": "Creates multiple precursor tolerance windows with specified mass offsets.",
    "mass_offsets_detailed": "Optional detailed mass offset list. Overrides mass_offsets if use_detailed_offsets = 1.",
    "use_detailed_offsets": "Whether to use the regular (0) or detailed (1) mass offset list.",
    "restrict_deltamass_to": "Specify amino acids on which delta masses (mass offsets or search modifications) can occur. Allowed values are single letter codes (e.g. ACD) and '-', must be capitalized. Use 'all' to allow any amino acid.",
    "labile_search_mode": "type of search (nglycan, labile, or off). Off means non-labile/typical search.",
    "precursor_mass_mode": "One of isolated/selected/corrected.",
    "localize_delta_mass": "Include fragment ions mass-shifted by unknown modifications (recommended for open and mass offset searches) (0 for OFF, 1 for ON).",
    "delta_mass_exclude_ranges": "Exclude mass range for shifted ions searching.",
    "fragment_ion_series": "Ion series used in search, specify any of a,b,c,x,y,z,Y,b-18,y-18 (comma separated).",
    "labile_fragment_ion_series": "Labile ion series used in labile search, specify any of a,b,c,x,y,z (comma separated). Replaces intact mods with labile ion(s) for each specified ion type.",
    "ion_series_definitions": 'User defined ion series. Example: "b* N -17.026548;b0 N -18.010565".',
    "search_enzyme_name_1": "Name of the first enzyme.",
    "search_enzyme_name_2": "Name of the second enzyme.",
    "search_enzyme_cut_1": "First enzyme's cutting amino acid.",
    "search_enzyme_nocut_1": "First enzyme's protecting amino acid.",
    "search_enzyme_cut_2": "Second enzyme's cutting amino acid.",
    "search_enzyme_nocut_2": "Second enzyme's protecting amino acid.",
    "num_enzyme_termini": "0 for non-enzymatic, 1 for semi-enzymatic, and 2 for fully-enzymatic.",
    "allowed_missed_cleavage_1": "First enzyme's allowed number of missed cleavages per peptide. Maximum value is 5.",
    "allowed_missed_cleavage_2": "Second enzyme's allowed number of missed cleavages per peptide. Maximum value is 5.",
    "search_enzyme_sense_1": "First enzyme's cutting terminal.",
    "search_enzyme_sense_2": "Second enzyme's cutting terminal.",
    "clip_nTerm_M": "Specifies the trimming of a protein N-terminal methionine as a variable modification (0 or 1).",
    "max_variable_mods_per_peptide": "Maximum total number of variable modifications per peptide.",
    "max_variable_mods_combinations": "Maximum number of modified forms allowed for each peptide (up to 65534).",
    "mass_diff_to_variable_mod": "Put mass diff as a variable modification. 0 for no; 1 for yes and remove delta mass; 2 for yes and keep delta mass.",
    "group_variable": "Specify the variable used to decide the PSM group in the group FDR estimation. 0 = no group FDR; 1 = num_enzyme_termini; 2 = PE from protein header.",
    "output_format": "File format of output files (tsv, pin, pepxml, tsv_pin, tsv_pepxml, pepxml_pin, or tsv_pepxml_pin).",
    "output_report_topN": "Reports top N PSMs per input spectrum.",
    "output_max_expect": "Suppresses reporting of PSM if top hit has expectation value greater than this threshold.",
    "report_alternative_proteins": "Report alternative proteins for peptides that are found in multiple proteins (0 for no, 1 for yes).",
    "precursor_charge": "Assumed range of potential precursor charge states. Only relevant when override_charge is set to 1.",
    "override_charge": "Ignores precursor charge and uses charge state specified in precursor_charge range (0 or 1).",
    "digest_min_length": "Minimum length of peptides to be generated during in-silico digestion.",
    "digest_max_length": "Maximum length of peptides to be generated during in-silico digestion.",
    "digest_mass_range": "Mass range of peptides to be generated during in-silico digestion in Daltons.",
    "max_fragment_charge": "Maximum charge state for theoretical fragments to match (1-4).",
    "track_zero_topN": "Track top N unmodified peptide results separately from main results internally for boosting features.",
    "zero_bin_accept_expect": "Ranks a zero-bin hit above all non-zero-bin hit if it has expectation less than this value.",
    "zero_bin_mult_expect": "Multiplies expect value of PSMs in the zero-bin during  results ordering (set to less than 1 for boosting).",
    "minimum_peaks": "Minimum number of peaks in experimental spectrum for matching.",
    "use_topN_peaks": "Pre-process experimental spectrum to only use top N peaks.",
    "deisotope": "Perform deisotoping or not (0=no, 1=yes and assume singleton peaks single charged, 2=yes and assume singleton peaks single or double charged).",
    "deneutralloss": "Perform deneutrallossing or not (0=no, 1=yes).",
    "min_fragments_modelling": "Minimum number of matched peaks in PSM for inclusion in statistical modeling.",
    "min_matched_fragments": "Minimum number of matched peaks for PSM to be reported.",
    "minimum_ratio": "Filters out all peaks in experimental spectrum less intense than this multiple of the base peak intensity.",
    "clear_mz_range": "Removes peaks in this m/z range prior to matching.",
    "remove_precursor_peak": " Remove precursor peaks from tandem mass spectra. 0 = not remove; 1 = remove the peak with precursor charge; 2 = remove the peaks with all charge states (only for DDA mode).",
    "remove_precursor_range": "m/z range in removing precursor peaks. Only for DDA mode. Unit: Th.",
    "intensity_transform": "Transform peaks intensities with sqrt root. 0 = not transform; 1 = transform using sqrt root.",
    "check_spectral_files": "Checking spectral files before searching.",
    "reuse_dia_fragment_peaks": "Allow the same peak matches to multiple peptides. For DIA data type only. 0 = no, 1 = yes.",
    "require_precursor": "If required, PSMs with no precursor peaks will be discarded. For DIA data type only. 0 = no, 1 = yes.",
    "Y_type_masses": " [nglycan/labile search_mode only]. Specify fragments of labile mods that are commonly retained on intact peptides (e.g. Y ions for glycans). Only used if 'Y' is included in fragment_ion_series.",
    "diagnostic_fragments": "[nglycan/labile search_mode only]. Specify diagnostic fragments of labile mods that appear in the low m/z region. Only used if diagnostic_intensity_filter > 0.",
    "remainder_fragment_masses": "[labile search_mode only] List of possible remainder fragment ions to consider. Remainder masses are partial modification masses left on b/y ions after fragmentation.",
    "min_sequence_matches": "[nglycan/labile search_mode only] Minimum number of sequence-specific (not Y) ions to record a match.",
    "diagnostic_intensity_filter": "[nglycan/labile search_mode only]. Minimum relative intensity for SUM of all detected oxonium ions to achieve for spectrum to contain diagnostic fragment evidence. Calculated relative to spectrum base peak. 0 <= value.",
    "activation_types": "Filter to only search scans of provided activation type(s). Allowed: All, HCD, CID, ETD, ECD.",
    "analyzer_types": "Filter to only include scans matching the provided analyzer type(s) in search. Only support the mzML and raw format. Allowed types: all, FTMS, ITMS, ASTMS.",
}

# MsfraggerParams.PROP_NAMES (FP/tools/fragger/MsfraggerParams.java:161-246): paramsFromMap only
# copies form values whose name is in this set (FP/tabs/TabMsfragger.java:1789).
FRAGGER_PROP_NAMES = frozenset(
    "labile_search_mode restrict_deltamass_to diagnostic_intensity_filter Y_type_masses "
    "diagnostic_fragments remainder_fragment_masses min_sequence_matches database_name "
    "decoy_prefix keep_decoys fragpipe_ram num_threads precursor_mass_lower precursor_mass_upper "
    "precursor_mass_units precursor_true_tolerance precursor_true_units fragment_mass_tolerance "
    "fragment_mass_units data_type remove_precursor_peak remove_precursor_range intensity_transform "
    "reuse_dia_fragment_peaks check_spectral_files require_precursor write_calibrated_mzml "
    "write_uncalibrated_mzml write_mzbin_all mass_diff_to_variable_mod group_variable calibrate_mass "
    "use_all_mods_in_first_search isotope_error deisotope deneutralloss mass_offsets "
    "mass_offsets_detailed use_detailed_offsets precursor_mass_mode search_enzyme_name_1 "
    "search_enzyme_name_2 search_enzyme_cut_1 search_enzyme_nocut_1 search_enzyme_cut_2 "
    "search_enzyme_nocut_2 num_enzyme_termini allowed_missed_cleavage_1 allowed_missed_cleavage_2 "
    "search_enzyme_sense_1 search_enzyme_sense_2 clip_nTerm_M variable_mod "
    "max_variable_mods_per_peptide max_variable_mods_combinations output_format output_report_topN "
    "report_alternative_proteins output_max_expect precursor_charge override_charge ms_level "
    "digest_min_length digest_max_length digest_mass_range max_fragment_charge fragment_ion_series "
    "labile_fragment_ion_series ion_series_definitions track_zero_topN zero_bin_accept_expect "
    "zero_bin_mult_expect localize_delta_mass delta_mass_exclude_ranges minimum_peaks use_topN_peaks "
    "min_fragments_modelling min_matched_fragments minimum_ratio clear_mz_range add add_enabled "
    "activation_types analyzer_types".split()
)

# MsfraggerParams.ADDON_NAMES / ADDONS_HUMAN_READABLE (FP/tools/fragger/MsfraggerParams.java:253-264):
# fixed-mod table site label -> add_<name> params key, in table row order.
FRAGGER_ADDONS = (
    ("C-Term Peptide", "Cterm_peptide"), ("N-Term Peptide", "Nterm_peptide"),
    ("C-Term Protein", "Cterm_protein"), ("N-Term Protein", "Nterm_protein"),
    ("G (glycine)", "G_glycine"), ("A (alanine)", "A_alanine"), ("S (serine)", "S_serine"),
    ("P (proline)", "P_proline"), ("V (valine)", "V_valine"), ("T (threonine)", "T_threonine"),
    ("C (cysteine)", "C_cysteine"), ("L (leucine)", "L_leucine"), ("I (isoleucine)", "I_isoleucine"),
    ("N (asparagine)", "N_asparagine"), ("D (aspartic acid)", "D_aspartic_acid"),
    ("Q (glutamine)", "Q_glutamine"), ("K (lysine)", "K_lysine"), ("E (glutamic acid)", "E_glutamic_acid"),
    ("M (methionine)", "M_methionine"), ("H (histidine)", "H_histidine"),
    ("F (phenylalanine)", "F_phenylalanine"), ("R (arginine)", "R_arginine"),
    ("Y (tyrosine)", "Y_tyrosine"), ("W (tryptophan)", "W_tryptophan"),
    ("B ", "B_user_amino_acid"), ("J", "J_user_amino_acid"), ("O", "O_user_amino_acid"),
    ("U", "U_user_amino_acid"), ("X", "X_user_amino_acid"), ("Z", "Z_user_amino_acid"),
)
_ADDON_HUMAN2NAME = dict(FRAGGER_ADDONS)
_ADDON_NAME2HUMAN = {name: human for human, name in FRAGGER_ADDONS}

# msfragger.* lines of FP/ui.translations: key -> [(value in UI, value in workflow)], file order.
FRAGGER_UI_TRANSLATIONS = {
    "calibrate_mass": [("None", "0"), ("Off", "0"), ("Mass calibration", "1"),
                       ("Mass calibration, parameter optimization", "2")],
    "fragment_mass_units": [("PPM", "1"), ("Da", "0")],
    "labile_search_mode": [("off", "off"), ("off", "standard"), ("nglycan", "nglycan"), ("labile", "labile")],
    "intensity_transform": [("None", "0"), ("Square root", "1")],
    "mass_diff_to_variable_mod": [("No", "0"), ("Yes, keep delta mass", "2"), ("Yes, remove delta mass", "1")],
    "deisotope": [("No", "0"), ("Yes", "1"), ("Yes, use charge 1 and 2 for undeisotoped peaks", "2")],
    "deneutralloss": [("Yes", "1"), ("No", "0")],
    "group_variable": [("None", "0"), ("Number of enzymatic termini", "1"),
                       ("Protein evidence from FASTA file", "2")],
    "num_enzyme_termini": [("SEMI_N_TERM", "3"), ("ENZYMATIC", "2"), ("SEMI", "1"), ("NONSPECIFIC", "0")],
    "output_format": [("PEPXML", "pepXML"), ("TSV", "tsv"), ("TSV_PEPXML", "tsv_pepXML"), ("PIN", "PIN"),
                      ("TSV_PIN", "tsv_pin"), ("PEPXML_PIN", "pepXML_pin"),
                      ("TSV_PEPXML_PIN", "tsv_pepXML_pin")],
    "precursor_mass_mode": [("selected", "selected"), ("isolated", "isolated"), ("corrected", "corrected")],
    "precursor_mass_units": [("PPM", "1"), ("Da", "0")],
    "data_type": [("DDA", "0"), ("DIA", "1"), ("GPF-DIA", "2"), ("DDA+", "3")],
    "precursor_true_units": [("PPM", "1"), ("Da", "0")],
    "remove_precursor_peak": [("Do not remove", "0"), ("Only peak with precursor charge", "1"),
                              ("Peaks with all charge states", "2")],
}

# MSFragger enzymes, FP/tools/fragger/enzyme-defs-supported.txt read by EnzymeProvider (sorted by
# name): name -> (cut, nocuts, sense).
FRAGGER_ENZYMES = {
    "argc": ("R", "P", "C"), "aspn": ("D", "", "N"), "chymotrypsin": ("FLWY", "P", "C"),
    "cnbr": ("M", "", "C"), "elastase": ("AGILV", "P", "C"), "formicacid": ("D", "P", "C"),
    "gluc": ("DE", "P", "C"), "gluc_bicarb": ("E", "P", "C"), "lysc": ("K", "P", "C"),
    "lysc-p": ("K", "", "C"), "lysn": ("K", "", "N"), "lysn_promisc": ("KR", "", "N"),
    "nocleavage": ("@", "@", "C"), "nonspecific": ("@", "", "C"), "null": ("", "", "C"),
    "stricttrypsin": ("KR", "", "C"), "thermolysin": ("ALIVFM", "DE", "N"),
    "trypsin": ("KR", "P", "C"), "trypsin/chymotrypsin": ("FYWLKR", "P", "C"),
    "trypsin/cnbr": ("KR", "P", "C"), "trypsin_gluc": ("DEKR", "P", "C"),
    "trypsin_k": ("K", "P", "C"), "trypsin_r": ("R", "P", "C"),
}
_ENZYME_NAMES = sorted(FRAGGER_ENZYMES)

_CALIBRATE_LABELS = ["None", "Mass calibration", "Mass calibration, parameter optimization"]
_MASS_DIFF_TO_VAR_MOD = ["No", "Yes, keep delta mass", "Yes, remove delta mass"]
_MASS_DIFF_TO_VAR_MOD_MAP = [0, 2, 1]
_GROUP_VARIABLE = ["None", "Number of enzymatic termini", "Protein evidence from FASTA file"]
_DEISOTOPE = ["No", "Yes", "Yes, use charge 1 and 2 for undeisotoped peaks"]
_REMOVE_PRECURSOR_PEAK = ["Do not remove", "Only peak with precursor charge", "Peaks with all charge states"]
_INTENSITY_TRANSFORM = ["None", "Square root"]
_CLEAVAGE = {"ENZYMATIC": 2, "SEMI": 1, "SEMI_N_TERM": 3, "NONSPECIFIC": 0}
_OUTPUT_FORMAT = {"PEPXML": "pepXML", "TSV": "tsv", "TSV_PEPXML": "tsv_pepXML", "PIN": "pin",
                  "TSV_PIN": "tsv_pin", "PEPXML_PIN": "pepXML_pin", "TSV_PEPXML_PIN": "tsv_pepXML_pin"}
_UNITS = {"PPM": 1, "Da": 0}


def _java_trim(text: str) -> str:
    """String.trim(): strips leading/trailing chars <= U+0020."""
    start, end = 0, len(text)
    while start < end and text[start] <= " ":
        start += 1
    while end > start and text[end - 1] <= " ":
        end -= 1
    return text[start:end]


def _java_parse_boolean(text: str) -> bool:
    """Boolean.parseBoolean."""
    return text is not None and text.lower() == "true"


def _java_parse_int(text: str) -> int:
    """Integer.parseInt (no surrounding whitespace allowed)."""
    if not re.fullmatch(r"[+-]?[0-9]+", text):
        raise ValueError(f"For input string: \"{text}\" (Integer.parseInt)")
    return int(text)


def _to_float32(value: float) -> float:
    return struct.unpack("f", struct.pack("f", value))[0]


def _java_fp_to_string(value: float, is_float: bool) -> str:
    """Float.toString / Double.toString: shortest round-trip digits in Java's layout
    (plain for 1e-3 <= |v| < 1e7 with at least one fraction digit, else d.dddE<exp>)."""
    if value != value:
        return "NaN"
    if value in (float("inf"), float("-inf")):
        return "Infinity" if value > 0 else "-Infinity"
    if value == 0:
        return "-0.0" if str(value).startswith("-") else "0.0"
    # FloatingDecimal.BinaryToASCIIBuffer.dtoa: an integer value of at most 62 bits takes the
    # long path (developLongDigits), which keeps the exact integer digits, dropping only
    # insignificantDigitsForPow2(binExp - significandBits - 1) of them (half up); larger floats print
    # more digits than the shortest round trip, e.g. Float.toString(2.5e10f) = "2.49999995E10".
    significand_bits = 24 if is_float else 53
    mantissa, bin_exp = math.frexp(abs(value))
    bin_exp -= 1
    if value == int(value) and bin_exp <= 62:
        lvalue = int(abs(value))
        insignificant = 0
        if bin_exp > significand_bits:
            p2 = bin_exp - significand_bits - 1
            insignificant = _INSIGNIFICANT_DIGITS[p2] if 1 < p2 < len(_INSIGNIFICANT_DIGITS) else 0
        if insignificant:
            pow10 = 10 ** insignificant
            residue = lvalue % pow10
            lvalue //= pow10
            if residue >= pow10 >> 1:
                lvalue += 1
        text = f"{'-' if value < 0 else ''}{lvalue}e{insignificant}"
    elif is_float:
        for precision in range(1, 10):
            text = f"{value:.{precision - 1}e}"
            if _to_float32(float(text)) == value:
                break
    else:
        text = repr(value)
    sign, digits, exponent = Decimal(text).as_tuple()
    digits = list(digits)
    while len(digits) > 1 and digits[-1] == 0:
        digits.pop()
        exponent += 1
    sci = exponent + len(digits) - 1
    ds = "".join(map(str, digits))
    if -3 <= sci < 7:
        if sci >= 0:
            int_part = ds[: sci + 1].ljust(sci + 1, "0")
            frac = ds[sci + 1:] or "0"
        else:
            int_part = "0"
            frac = "0" * (-sci - 1) + ds
        out = f"{int_part}.{frac}"
    else:
        out = f"{ds[0]}.{ds[1:] or '0'}E{sci}"
    return ("-" if sign else "") + out


# FloatingDecimal.insignificantDigitsNumber
_INSIGNIFICANT_DIGITS = (0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7,
                         8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 14,
                         14, 14, 15, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19)


def java_float_to_string(value: float) -> str:
    """String.valueOf(float) of a value that is already a float32."""
    return _java_fp_to_string(_to_float32(value), True)




def java_double_to_string(value: float) -> str:
    """String.valueOf(double)."""
    return _java_fp_to_string(value, False)


def java_decimal_format(pattern: str, value: float) -> str:
    """new DecimalFormat(pattern).format(double) for the plain patterns FragPipe uses
    ("0.###", "#.##", "0.00", ...), no grouping. A pattern without '0' in the integer part
    ("#.##") still prints one integer digit (DecimalFormat.applyPattern turns "#." into "0.").

    Rounding is java.text.DigitList HALF_EVEN on the Double.toString digits: a tie in those digits
    is resolved by the exact binary value (up when the double is above the tie, down when below),
    a true tie goes to the even digit, and a tie on the first significant digit of a value below
    1e-3 (Double.toString in E notation) rounds to zero."""
    int_pattern, _, frac_pattern = pattern.partition(".")
    min_int = max(1, int_pattern.count("0"))
    min_frac = frac_pattern.count("0")
    max_frac = len(frac_pattern)
    negative = value < 0 or (value == 0 and str(value).startswith("-"))
    shortest, exact = Decimal(repr(abs(value))), Decimal(abs(value))
    unit = Decimal(1).scaleb(-max_frac)
    context = Context(prec=800)
    down = shortest.quantize(unit, rounding=ROUND_FLOOR, context=context)
    if shortest - down != unit / 2 or shortest == exact:
        rounded = shortest.quantize(unit, rounding=ROUND_HALF_EVEN, context=context)
    elif shortest > exact or (down == 0 and abs(value) < 1e-3):
        rounded = down
    else:
        rounded = down + unit
    int_digits, _, frac_digits = f"{rounded:f}".partition(".")
    frac_digits = frac_digits.rstrip("0")
    if len(frac_digits) < min_frac:
        frac_digits = frac_digits.ljust(min_frac, "0")
    int_digits = int_digits.lstrip("0")
    if len(int_digits) < min_int:
        int_digits = int_digits.rjust(min_int, "0")
    if not int_digits and not frac_digits:
        int_digits = "0"
    text = int_digits + ("." + frac_digits if frac_digits else "")
    return ("-" if negative else "") + text


def _java_format_f(value: float) -> str:
    """String.format("%f", double): 6 fraction digits, HALF_UP on Double.toString digits."""
    rounded = Decimal(java_double_to_string(value)).quantize(Decimal("0.000001"), rounding=ROUND_HALF_UP,
                                                              context=Context(prec=800))
    return f"{rounded:f}"


# Default number format of MsfraggerParams setters (FP/tools/fragger/MsfraggerParams.java:304).
_FRAGGER_DF = "0.##########"


class FraggerProps:
    """FP/params/Props.java: ordered MSFragger params, each with an enabled flag."""

    _DISABLED = re.compile(r"#\s*([^\s]+)\s*=([^#]+)(?:\s*#\s*(.+))?.*")
    _ENABLED = re.compile(r"\s*([^\s]+)\s*=([^#]*)(?:\s*#\s*(.+))?.*")
    _BLANK = "# blank line "

    def __init__(self):
        self.ordering = []
        self.values = {}  # name -> (value, enabled), insertion ordered like LinkedHashMap

    @classmethod
    def from_template(cls) -> "FraggerProps":
        """new MsfraggerParams(): Props.readProps of fragger_closed.params (Props.java:161-238)."""
        props = cls()
        blank = 0
        if FRAGGER_TEMPLATE is None:
            raise RuntimeError("fragger_closed.params is not set (--fragger-template)")
        for raw in FRAGGER_TEMPLATE.read_text(encoding="utf-8").split("\n")[:-1]:
            line = _java_trim(raw)
            if not line:
                props.ordering.append(f"{cls._BLANK}{blank}")
                blank += 1
                continue
            if line.startswith("#"):
                match = cls._DISABLED.fullmatch(line)
                if not match:
                    props.ordering.append(line)
                    continue
                name, value = match.group(1).strip(), match.group(2).strip()
                if name in props.values and props.values[name][1]:
                    continue
                props.values[name] = (value, False)
                props.ordering.append(name)
                continue
            match = cls._ENABLED.fullmatch(line)
            if match:
                name, value = match.group(1).strip(), match.group(2).strip()
                props.values[name] = (value, True)
                props.ordering.append(name)
        return props

    def copy(self) -> "FraggerProps":
        other = FraggerProps()
        other.ordering = list(self.ordering)
        other.values = dict(self.values)
        return other

    def get(self, name: str) -> str:
        return self.values[name][0]

    def set(self, name: str, value: str, enabled: bool = True) -> None:
        self.values[name] = (value, enabled)

    def remove(self, name: str) -> None:
        self.values.pop(name, None)

    def write(self) -> str:
        """Props.writeProps (FP/params/Props.java:240-273)."""
        out = []
        for name in dict.fromkeys(self.ordering + list(self.values)):
            if not name or name.startswith("#"):
                out.append(("" if name.startswith(self._BLANK) else name) + "\n")
                continue
            if name not in self.values:
                continue
            value, enabled = self.values[name]
            line = ("" if enabled else "# ") + f"{name} = {value}"
            comment = FRAGGER_COMMENTS.get(name)
            if comment is not None and comment.strip():
                line += f"\t\t\t# {comment}"
            out.append(line + "\n")
        out.append("\n")
        return "".join(out)


class _Widget:
    """One named TabMsfragger form component and its string round trip
    (utils/swing/UiCombo, UiSpinnerDouble, UiSpinnerInt, UiCheck, UiText; SwingUtils.valueGet/valueSet)."""

    def __init__(self, kind: str, init, options=None, fmt=None, filt=None):
        self.kind, self.value, self.options, self.fmt, self.filt = kind, init, options, fmt, filt
        if kind == "text":
            self.value = self._filtered(init)

    def _filtered(self, text: str) -> str:
        # DocumentFilters.getFilter(regex, replacement): applied to text set into the field.
        if self.filt is None:
            return text
        regex, replacement = self.filt
        return re.sub(regex, replacement, text)

    def set(self, text: str) -> bool:
        """valueSet(comp, text); returns True when a combo selection changed (fires ItemListeners)."""
        old = self.value
        if self.kind == "combo":
            self.value = text if text in self.options else self.options[0]
        elif self.kind == "spinD":
            self.value = float(text)
        elif self.kind == "spinI":
            self.value = _java_parse_int(text)
        elif self.kind == "check":
            self.value = _java_parse_boolean(text)
        else:
            self.value = self._filtered(text)
        return self.kind == "combo" and old != self.value

    def get(self) -> str:
        """valueGet(comp) / StringRepresentable.asString()."""
        if self.kind == "spinD":
            return java_decimal_format(self.fmt, self.value)
        if self.kind == "spinI":
            return str(self.value)
        if self.kind == "check":
            return "true" if self.value else "false"
        if self.kind == "text":
            return _java_trim(self.value)
        return self.value


def _fragger_widgets() -> Dict[str, _Widget]:
    """Named TabMsfragger components with their constructor state (FP/tabs/TabMsfragger.java:567-1600)."""
    combo, spin_d, spin_i, check, text = (
        lambda options, init=None: _Widget("combo", options[0] if init is None else init, options=options),
        lambda init, fmt: _Widget("spinD", float(init), fmt=fmt),
        lambda init: _Widget("spinI", init),
        lambda init: _Widget("check", init),
        lambda init="", filt=None: _Widget("text", init, filt=filt),
    )
    units = ["PPM", "Da"]
    enzyme_cuts = (r"[^A-Z@]", "@")
    return {
        "run-msfragger": check(True),
        "digest-only": check(False),
        "precursor_mass_units": combo(units),
        "precursor_mass_lower": spin_d(-20, "0.###"),
        "precursor_mass_upper": spin_d(20, "0.###"),
        "fragment_mass_units": combo(units),
        "fragment_mass_tolerance": spin_d(10, "0.###"),
        "calibrate_mass": combo(_CALIBRATE_LABELS),
        "isotope_error": text("-1/0/1/2", (r"[^0-9/-]+", "")),
        "misc.fragger.enzyme-dropdown-1": combo(_ENZYME_NAMES, "trypsin"),
        "misc.fragger.enzyme-dropdown-2": combo(_ENZYME_NAMES, "trypsin"),
        "search_enzyme_name_1": text(),
        "search_enzyme_cut_1": text("KR", enzyme_cuts),
        "search_enzyme_nocut_1": text("P", enzyme_cuts),
        "num_enzyme_termini": combo(list(_CLEAVAGE)),
        "allowed_missed_cleavage_1": spin_i(2),
        "search_enzyme_sense_1": combo(["N", "C"]),
        "search_enzyme_name_2": text(),
        "search_enzyme_cut_2": text("KR", enzyme_cuts),
        "search_enzyme_nocut_2": text("P", enzyme_cuts),
        "allowed_missed_cleavage_2": spin_i(2),
        "search_enzyme_sense_2": combo(["N", "C"]),
        "clip_nTerm_M": check(False),
        "digest_min_length": spin_i(7),
        "digest_max_length": spin_i(50),
        "misc.fragger.digest-mass-lo": spin_d(200, "0.#"),
        "misc.fragger.digest-mass-hi": spin_d(5000, "0.#"),
        "misc.slice-db": spin_i(1),
        "labile_search_mode": combo(["off", "nglycan", "labile"]),
        "diagnostic_intensity_filter": spin_d(0, "#.##"),
        "min_sequence_matches": spin_i(2),
        "Y_type_masses": text(),
        "diagnostic_fragments": text(),
        "remainder_fragment_masses": text(),
        "labile_fragment_ion_series": text("b,y"),
        "max_variable_mods_per_peptide": spin_i(3),
        "max_variable_mods_combinations": spin_i(5000),
        "use_all_mods_in_first_search": check(False),
        "restrict_deltamass_to": text("", (r"[^A-Zal\-nc\[\]^*]", "")),
        "mass_offsets": text(),
        "use_detailed_offsets": check(False),
        "mass_offsets_detailed": text(),
        "minimum_peaks": spin_i(15),
        "use_topN_peaks": spin_i(150),
        "minimum_ratio": spin_d(0.01, "0.00"),
        "misc.fragger.clear-mz-lo": spin_d(0, "#.###"),
        "misc.fragger.clear-mz-hi": spin_d(0, "#.###"),
        "precursor_mass_mode": combo(["selected", "isolated", "corrected"]),
        "remove_precursor_peak": combo(_REMOVE_PRECURSOR_PEAK),
        "misc.fragger.remove-precursor-range-lo": spin_d(-1.5, "0.#"),
        "misc.fragger.remove-precursor-range-hi": spin_d(1.5, "0.#"),
        "intensity_transform": combo(_INTENSITY_TRANSFORM),
        "reuse_dia_fragment_peaks": check(False),
        "check_spectral_files": check(True),
        "activation_types": combo(["all", "HCD", "ETD", "CID", "ECD"]),
        "analyzer_types": combo(["all", "FTMS", "ITMS", "ASTMS"]),
        "require_precursor": check(True),
        "min_fragments_modelling": spin_i(2),
        "min_matched_fragments": spin_i(4),
        "fragment_ion_series": text(),
        "ion_series_definitions": text(),
        "precursor_true_units": combo(units),
        "precursor_true_tolerance": spin_d(10, "0.###"),
        "misc.fragger.precursor-charge-lo": spin_i(1),
        "misc.fragger.precursor-charge-hi": spin_i(4),
        "override_charge": check(False),
        "deisotope": combo(_DEISOTOPE),
        "deneutralloss": combo(["Yes", "No"]),
        "max_fragment_charge": spin_i(2),
        "output_report_topN": spin_i(1),
        "output_report_topN_dia1": spin_i(5),
        "output_report_topN_dda_plus": spin_i(5),
        "output_max_expect": spin_d(50, "0.#"),
        "output_format": combo(list(_OUTPUT_FORMAT)),
        "report_alternative_proteins": check(False),
        "write_calibrated_mzml": check(False),
        "group_variable": combo(_GROUP_VARIABLE),
        "track_zero_topN": spin_i(0),
        "zero_bin_accept_expect": spin_d(0, "0.##########"),
        "zero_bin_mult_expect": spin_d(1, "0.##########"),
        "mass_diff_to_variable_mod": combo(_MASS_DIFF_TO_VAR_MOD),
        "delta_mass_exclude_ranges": text("(-1.5,3.5)", (r"[A-Za-z]", "")),
        "localize_delta_mass": check(False),
    }


# Widgets each ItemListener rewrites when its combo selection changes
# (FP/tabs/TabMsfragger.java:699-753 enzyme dropdowns, :1454-1462 deisotope).
_FRAGGER_LISTENER_TARGETS = {
    "misc.fragger.enzyme-dropdown-1": ("search_enzyme_name_1", "search_enzyme_cut_1",
                                       "search_enzyme_nocut_1", "search_enzyme_sense_1",
                                       "num_enzyme_termini"),
    "misc.fragger.enzyme-dropdown-2": ("search_enzyme_name_2", "search_enzyme_cut_2",
                                       "search_enzyme_nocut_2", "search_enzyme_sense_2",
                                       "num_enzyme_termini"),
    "deisotope": ("deneutralloss",),
}


class _FraggerForm:
    """TabMsfragger form state: named widgets plus the variable/fixed mod tables."""

    def __init__(self):
        self.widgets = _fragger_widgets()
        self.var_mods = []  # Mod rows: (mass float32, sites, enabled, max occurrences)
        self.fix_mods = []

    def set(self, name: str, text: str) -> None:
        widget = self.widgets.get(name)
        if widget is None:
            return
        if widget.set(text):
            self._item_listener(name)

    def _item_listener(self, name: str) -> None:
        if name.startswith("misc.fragger.enzyme-dropdown-"):
            n = name[-1]
            item = self.widgets[name].value
            cut, nocuts, sense = FRAGGER_ENZYMES[item]
            self.widgets[f"search_enzyme_name_{n}"].set("nonspecific" if item.lower() == "custom" else item)
            self.widgets[f"search_enzyme_cut_{n}"].set(cut)
            self.widgets[f"search_enzyme_nocut_{n}"].set(nocuts)
            # JComboBox.setSelectedItem ignores items that are not in the list.
            if sense in self.widgets[f"search_enzyme_sense_{n}"].options:
                self.widgets[f"search_enzyme_sense_{n}"].value = sense
            if item == "nonspecific":
                self.widgets["num_enzyme_termini"].value = "NONSPECIFIC"
        elif name == "deisotope":
            index = _DEISOTOPE.index(self.widgets["deisotope"].value)
            self.widgets["deneutralloss"].value = "Yes" if index > 0 else "No"

    def get(self, name: str) -> str:
        return self.widgets[name].get()

    @staticmethod
    def _table_rows(text: str) -> list:
        """ModsTable.fromString -> Mod.fromString per "; " chunk (FP/api/ModsTable.java,
        FP/tools/fragger/Mod.java:37-42). Chunks without exactly 4 fields are skipped."""
        rows = []
        for chunk in _java_split(text, "; "):
            parts = _java_split(chunk, ",")
            if len(parts) != 4:
                continue
            rows.append((_to_float32(float(parts[0])), parts[1], _java_parse_boolean(parts[2]),
                         _java_parse_int(parts[3])))
        return rows

    def set_var_mods_text(self, text: str) -> None:
        rows = self._table_rows(text)
        if len(rows) > 16:
            raise ValueError("loaded mods list length is longer than msfragger supports")
        self.var_mods = rows

    def set_fix_mods_text(self, text: str) -> None:
        rows = self._table_rows(text)
        if len(rows) > len(FRAGGER_ADDONS):
            raise ValueError("mod list length larger than fix mods used by fragger")
        # convertModsToFixTableData: default rows (disabled, 0.0) for every addon, then the loaded mods.
        table = [(0.0, human, False, -1) for human, _ in FRAGGER_ADDONS]
        for i, row in enumerate(rows):
            table[i] = row
        self.fix_mods = table


def _java_split(text: str, delimiter: str) -> list:
    """String.split with a literal delimiter: no match gives [text], else trailing empty
    strings are removed."""
    if delimiter not in text:
        return [text]
    parts = text.split(delimiter)
    while parts and parts[-1] == "":
        parts.pop()
    return parts


def _translate_value_to_ui(key: str, value: str) -> str:
    """FragpipeCacheUtils.translateValuesToUi for one msfragger.* key."""
    translations = FRAGGER_UI_TRANSLATIONS.get(key)
    if translations is None:
        return value
    if any(in_ui == value for in_ui, _ in translations):
        return value
    for in_ui, in_conf in translations:
        if in_conf.lower() == value.lower():
            return in_ui
    return value


def _fragger_convert_to_gui(key: str, value: str) -> str:
    """TabMsfragger.CONVERT_TO_GUI (FP/tabs/TabMsfragger.java:316-375)."""
    as_bool = lambda s: "true" if int(s) > 0 else "false"  # noqa: E731
    converters = {
        "write_calibrated_mzml": as_bool,
        "write_uncalibrated_mzml": lambda s: "true",
        "write_mzbin_all": as_bool,
        "mass_diff_to_variable_mod": lambda s: _MASS_DIFF_TO_VAR_MOD[_MASS_DIFF_TO_VAR_MOD_MAP[int(s)]],
        "group_variable": lambda s: _GROUP_VARIABLE[int(s)],
        "deisotope": lambda s: _DEISOTOPE[int(s)],
        "precursor_mass_units": lambda s: {1: "PPM", 0: "Da"}[int(s)],
        "fragment_mass_units": lambda s: {1: "PPM", 0: "Da"}[int(s)],
        "precursor_true_units": lambda s: {1: "PPM", 0: "Da"}[int(s)],
        "calibrate_mass": lambda s: {"0": "None", "1": "Mass calibration",
                                     "2": "Mass calibration, parameter optimization"}.get(s, "None"),
        "use_all_mods_in_first_search": lambda s: "true" if int(s) == 1 else "false",
        "num_enzyme_termini": lambda s: {str(v): k for k, v in _CLEAVAGE.items()}[s],
        "remove_precursor_peak": lambda s: _REMOVE_PRECURSOR_PEAK[int(s)],
        "intensity_transform": lambda s: _INTENSITY_TRANSFORM[int(s)],
        "check_spectral_files": as_bool,
        "reuse_dia_fragment_peaks": as_bool,
        "require_precursor": as_bool,
        "activation_types": lambda s: {"ALL": "all", "": "all", "HCD": "HCD", "ETD": "ETD", "CID": "CID",
                                       "ECD": "ECD"}[s.upper()],
        "analyzer_types": lambda s: {"": "all", "ALL": "all", "FTMS": "FTMS", "ITMS": "ITMS",
                                     "ASTMS": "ASTMS"}[s.upper()],
        "localize_delta_mass": as_bool,
        "clip_nTerm_M": as_bool,
        "override_charge": as_bool,
        "output_format": lambda s: next(k for k, v in _OUTPUT_FORMAT.items() if v.lower() == s.lower()),
        "report_alternative_proteins": as_bool,
        "labile_search_mode": lambda s: s,
        "mass_offsets": lambda s: "0" if not _java_trim(s) else " ".join(_java_trim(s).split("/")),
        "mass_offsets_detailed": lambda s: s if _java_trim(s) else "",
        "use_detailed_offsets": as_bool,
        "Y_type_masses": lambda s: " ".join(_java_split(s, "/")),
        "diagnostic_fragments": lambda s: " ".join(_java_split(s, "/")),
        "remainder_fragment_masses": lambda s: " ".join(_java_split(s, "/")),
        "labile_fragment_ion_series": lambda s: "b,y" if not _java_trim(s) else s,
        "deneutralloss": lambda s: "Yes" if int(s) == 1 else "No",
    }
    convert = converters.get(key)
    return value if convert is None else convert(value)


def _fragger_convert_to_file(key: str, value: str) -> str:
    """TabMsfragger.CONVERT_TO_FILE (FP/tabs/TabMsfragger.java:236-314)."""
    as_int = lambda s: "1" if _java_parse_boolean(s) else "0"  # noqa: E731

    def mass_offsets(s):
        if s is None or not _java_trim(s):
            return "0"
        parts = re.split(r"[\s/]+", _java_trim(s))
        return "/".join(sorted({java_float_to_string(float(p)) for p in parts}))

    def index_of(options, s, default):
        for i, option in enumerate(options):
            if s.lower() == option.lower():
                return str(i)
        return default

    converters = {
        "write_calibrated_mzml": as_int,
        "write_uncalibrated_mzml": lambda s: "1",
        "write_mzbin_all": as_int,
        "mass_diff_to_variable_mod": lambda s: str(_MASS_DIFF_TO_VAR_MOD_MAP[_MASS_DIFF_TO_VAR_MOD.index(s)]),
        "group_variable": lambda s: index_of(_GROUP_VARIABLE, s, "0"),
        "deisotope": lambda s: index_of(_DEISOTOPE, s, "1"),
        "precursor_mass_units": lambda s: str(_UNITS[s]),
        "fragment_mass_units": lambda s: str(_UNITS[s]),
        "precursor_true_units": lambda s: str(_UNITS[s]),
        "calibrate_mass": lambda s: str(_CALIBRATE_LABELS.index(s)) if s in _CALIBRATE_LABELS else "0",
        "use_all_mods_in_first_search": as_int,
        "num_enzyme_termini": lambda s: str(_CLEAVAGE[s]),
        "remove_precursor_peak": lambda s: str(_REMOVE_PRECURSOR_PEAK.index(s)),
        "intensity_transform": lambda s: str(_INTENSITY_TRANSFORM.index(s)),
        "check_spectral_files": as_int,
        "reuse_dia_fragment_peaks": as_int,
        "require_precursor": as_int,
        "localize_delta_mass": as_int,
        "clip_nTerm_M": as_int,
        "override_charge": as_int,
        "output_format": lambda s: _OUTPUT_FORMAT[s],
        "report_alternative_proteins": as_int,
        "fragment_ion_series": lambda s: re.sub(r"[\s,;]+", ",", _java_trim(s)),
        "labile_fragment_ion_series": lambda s: re.sub(r"[\s,;]+", ",", _java_trim(s)),
        "mass_offsets": mass_offsets,
        "mass_offsets_detailed": lambda s: "" if s is None or not _java_trim(s) else s,
        "use_detailed_offsets": as_int,
        "Y_type_masses": lambda s: re.sub(r"\s+", "/", s),
        "diagnostic_fragments": lambda s: re.sub(r"\s+", "/", s),
        "remainder_fragment_masses": lambda s: re.sub(r"\s+", "/", s),
        "deneutralloss": lambda s: "1" if s.lower() == "yes" else "0",
    }
    convert = converters.get(key)
    return value if convert is None else convert(value)


def _java_hashmap_bucket(key: str, capacity: int) -> int:
    """Bucket of `key` in a java.util.HashMap of the given table size (String.hashCode, HashMap.hash)."""
    h = 0
    for ch in key:
        h = (31 * h + ord(ch)) & 0xFFFFFFFF
    return ((h ^ (h >> 16)) & 0xFFFFFFFF) & (capacity - 1)


def _java_hashmap_capacity(size: int) -> int:
    """Table size of a HashMap built by `size` puts into new HashMap<>() (Collectors.toMap)."""
    capacity = 16
    while size > capacity * 3 // 4:
        capacity *= 2
    return capacity


def _headless_ui_keys(all_params: Dict[str, Dict[str, str]]) -> set:
    """Keys of the map FragPipe headless passes to valuesSet: the workflow keys plus what
    Fragpipe.headless() adds (FP/Fragpipe.java:500-521: workflow.ram, workflow.threads, workdir,
    fragpipe-config.tools-folder for --config-tools-folder) and the run database the pipeline always
    gives (database.db-path: FragpipeUiState, FRAGPIPEHEADLESS headlessWorkflowText), filtered by
    TabWorkflow.filterPropsForUi (FP/tabs/TabWorkflow.java:917-923)."""
    keys = {k if prefix == "_global" else f"{prefix}.{k}"
            for prefix, params in all_params.items() for k in params}
    keys |= {"workflow.ram", "workflow.threads", "workdir", "fragpipe-config.tools-folder", "database.db-path"}
    return {k for k in keys if "workflow-option" in k.lower() or "do-not-cache" not in k.lower()}


def _apply_workflow_to_form(form: _FraggerForm, params: Dict[str, str],
                            all_params: Dict[str, Dict[str, str]]) -> None:
    """Stage 2: valuesSet over the workflow map. Its iteration order is java.util.HashMap order,
    which decides whether an enzyme dropdown / deisotope ItemListener runs before or after the
    explicit values it overwrites."""
    capacity = _java_hashmap_capacity(len(_headless_ui_keys(all_params)))
    names = [k for k in params if k in form.widgets]
    bucket = {k: _java_hashmap_bucket(f"msfragger.{k}", capacity) for k in names}
    for source, targets in _FRAGGER_LISTENER_TARGETS.items():
        for target in targets:
            if source in bucket and target in bucket and bucket[source] == bucket[target]:
                raise ValueError(
                    f"msfragger.{source} and msfragger.{target} share a HashMap bucket "
                    f"(capacity {capacity}); FragPipe's load order between them is not defined by the keys"
                )
    for name in sorted(names, key=lambda k: bucket[k]):
        form.set(name, _translate_value_to_ui(name, params[name]))
    if "table.var-mods" in params:
        form.set_var_mods_text(params["table.var-mods"])
    if "table.fix-mods" in params:
        form.set_fix_mods_text(params["table.fix-mods"])


def _form_from_template(form: _FraggerForm, template: FraggerProps) -> None:
    """Stage 1: TabMsfragger.formFrom(closed-search defaults) (FP/tabs/TabMsfragger.java:1639-1650,
    confToUiMap :1859-1890)."""
    ui = {}
    for name, (value, enabled) in template.values.items():
        if enabled:
            # confToUiMap logs and skips values its converter rejects (TabMsfragger.java:1862-1871).
            try:
                converted = _fragger_convert_to_gui(name, value)
            except (KeyError, ValueError, IndexError, StopIteration):
                continue
            ui[name] = _translate_value_to_ui(name, converted)
    remove_range = [float(v) for v in re.split(r"\s*,\s*", template.get("remove_precursor_range"))]
    clear_mz = [float(v) for v in re.split(r"\s+", template.get("clear_mz_range"))]
    digest_mass = [float(v) for v in re.split(r"\s+", template.get("digest_mass_range"))]
    charge = [int(v) for v in re.split(r"\s+", template.get("precursor_charge"))]
    ui["misc.fragger.clear-mz-lo"] = java_decimal_format("0.#", clear_mz[0])
    ui["misc.fragger.clear-mz-hi"] = java_decimal_format("0.#", clear_mz[1])
    ui["misc.fragger.digest-mass-lo"] = java_decimal_format("0.#", digest_mass[0])
    ui["misc.fragger.digest-mass-hi"] = java_decimal_format("0.#", digest_mass[1])
    ui["misc.fragger.precursor-charge-lo"] = java_decimal_format("0.#", charge[0])
    ui["misc.fragger.precursor-charge-hi"] = java_decimal_format("0.#", charge[1])
    ui["misc.fragger.remove-precursor-range-lo"] = java_decimal_format("0.##", remove_range[0])
    ui["misc.fragger.remove-precursor-range-hi"] = java_decimal_format("0.##", remove_range[1])
    for name, value in ui.items():
        form.set(name, value)
    form.widgets["misc.slice-db"].value = 1
    form.var_mods = []
    for i in range(1, 17):
        name = f"variable_mod_{i:02d}"
        if name not in template.values:
            continue
        value, enabled = template.values[name]
        split = re.split(r"\s+", value)
        form.var_mods.append((_to_float32(float(split[0])), split[1], enabled, int(split[2])))
    form.fix_mods = []
    for human, addon in FRAGGER_ADDONS:
        name = f"add_{addon}"
        if name in template.values:
            value, enabled = template.values[name]
            form.fix_mods.append((_to_float32(float(value)), human, enabled, 1))


def msfragger_form(params: Dict[str, str], all_params: Dict[str, Dict[str, str]]) -> _FraggerForm:
    """TabMsfragger after a headless load: the closed-search defaults of initMore
    (loadDefaults(SearchTypeProp.closed), TabMsfragger.java:515-517, :2088-2096), then the
    workflow's msfragger.* keys (stage 1 + stage 2)."""
    form = _FraggerForm()
    _form_from_template(form, FraggerProps.from_template())
    _apply_workflow_to_form(form, params, all_params)
    return form


def _translate_value_from_ui(key: str, value: str) -> str:
    """FragpipeCacheUtils.translateValueFromUi (FP/api/FragpipeCacheUtils.java:142-157): the first
    ui.translations entry whose UI value equals `value` ignoring case gives the workflow value."""
    translations = FRAGGER_UI_TRANSLATIONS.get(key)
    if translations is None:
        return value
    for in_ui, in_conf in translations:
        if in_ui.lower() == value.lower():
            return in_conf
    raise ValueError(f"Found no ui.translation for UI element key [msfragger.{key}] with value [{value}]")


def _mods_table_text(rows: list) -> str:
    """ModsTable.asString (FP/api/ModsTable.java:50-54, Mod.asString FP/tools/fragger/Mod.java:38-40):
    "massDelta,sites,isEnabled,maxOccurrences" per row (Float.toString), joined with "; "."""
    return "; ".join(
        f"{java_float_to_string(mass)},{sites},{'true' if enabled else 'false'},{max_occurrences}"
        for mass, sites, enabled, max_occurrences in rows
    )


def msfragger_values(all_params: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    """The msfragger.* values FragPipe's MSFragger tab holds for this workflow, in workflow form
    (what a workflow saved from that tab would contain): every widget and both mod tables of
    msfragger_form, so a key the workflow omits has its closed-search default
    (fragger_closed.params). Other FragPipe panels read these widgets, not the workflow map
    (e.g. getVarModMassSet TabMsfragger.java:1705-1725, getEnzymeName :2143-2148,
    getMassDiffToVariableMod :480-482, getMassCalibration :1944)."""
    params = all_params.get("msfragger", {})
    form = msfragger_form(params, all_params)
    values = dict(params)
    for name in form.widgets:
        values[name] = _translate_value_from_ui(name, form.get(name))
    values["table.var-mods"] = _mods_table_text(form.var_mods)
    values["table.fix-mods"] = _mods_table_text(form.fix_mods)
    return values


def msfragger_writes_calibrated_mzml(all_params: Dict[str, Dict[str, str]]) -> bool:
    """tabMsf.isWriteCalMzml() && tabMsf.getMassCalibration() > 0 (FragpipeRun.java:1612, :1642, :2090; TabMsfragger.java:1944-1958)."""
    values = msfragger_values(all_params)
    return values["write_calibrated_mzml"] == "true" and int(values["calibrate_mass"]) > 0


def _fragger_params_from_form(form: _FraggerForm, num_threads: str) -> FraggerProps:
    """Stage 3: TabMsfragger.formCollect / paramsFromMap (FP/tabs/TabMsfragger.java:1652-1674,
    :1780-1857) and MsfraggerParams.setVariableMods/setFixedMods (MsfraggerParams.java:1104-1162)."""
    props = FraggerProps.from_template()
    values = {name: form.get(name) for name in form.widgets}
    values["num_threads"] = num_threads
    values["group_variable"] = form.widgets["group_variable"].value
    for name, value in values.items():
        if name not in FRAGGER_PROP_NAMES:
            continue
        converted = _fragger_convert_to_file(name, value)
        if name in ("fragment_ion_series", "labile_fragment_ion_series") and not converted.strip():
            continue
        props.set(name, converted)
    remove_range = [float(values["misc.fragger.remove-precursor-range-lo"]),
                    float(values["misc.fragger.remove-precursor-range-hi"])]
    clear_mz = [float(values["misc.fragger.clear-mz-lo"]), float(values["misc.fragger.clear-mz-hi"])]
    digest_mass = [float(values["misc.fragger.digest-mass-lo"]), float(values["misc.fragger.digest-mass-hi"])]
    charge = [_java_parse_int(values["misc.fragger.precursor-charge-lo"]),
              _java_parse_int(values["misc.fragger.precursor-charge-hi"])]
    props.set("clear_mz_range", f"{java_double_to_string(clear_mz[0])} {java_double_to_string(clear_mz[1])}")
    props.set("digest_mass_range",
              f"{java_double_to_string(digest_mass[0])} {java_double_to_string(digest_mass[1])}")
    props.set("precursor_charge", f"{charge[0]} {charge[1]}")
    props.set("remove_precursor_range", f"{_java_format_f(remove_range[0])},{_java_format_f(remove_range[1])}")

    # setVariableMods: rows with sites (ModsTableModel.getModifications), numbered from 01.
    var_names = set()
    rows = [r for r in form.var_mods if r[1] is not None and r[1].strip()]
    for i, (mass, sites, enabled, max_occurrences) in enumerate(rows, start=1):
        name = f"variable_mod_{i:02d}"
        props.set(name, f"{java_float_to_string(mass)} {sites} {max_occurrences}", enabled)
        var_names.add(name)
    for name in [k for k in props.values if k.startswith("variable_mod") and k not in var_names]:
        props.remove(name)
    # setFixedMods
    fix_names = set()
    for mass, sites, enabled, _ in form.fix_mods:
        if sites is None or not sites.strip():
            continue
        addon = _ADDON_HUMAN2NAME.get(sites)
        if addon is None:
            raise ValueError(f"Could not map fixed modification site '{sites}' to an MSFragger add_ parameter")
        props.set(f"add_{addon}", java_float_to_string(mass), enabled)
        fix_names.add(f"add_{addon}")
    for name in [k for k in props.values if k.startswith("add_") and k not in fix_names]:
        props.remove(name)
    return props


def _fragger_data_type_variant(params: FraggerProps, data_type: str, topn_dia1: int,
                               topn_dda_plus: int) -> FraggerProps:
    """CmdMsfragger.adjustMSFraggerParams (FP/cmd/CmdMsfragger.java:710-741)."""
    variant = params.copy()
    if data_type == "DDA":
        variant.set("data_type", "0")
        return variant
    variant.set("report_alternative_proteins", "1")
    variant.set("localize_delta_mass", "0")
    variant.set("labile_search_mode", "off")
    variant.set("restrict_deltamass_to", "all")
    variant.set("isotope_error", "0")
    variant.set("intensity_transform", "1")
    variant.set("remove_precursor_peak", "1")
    if data_type == "DIA":
        variant.set("data_type", "1")
        variant.set("use_topN_peaks", str(max(300, int(params.get("use_topN_peaks")))))
        variant.set("minimum_ratio", java_double_to_string(0.0))
        variant.set("output_report_topN", str(topn_dia1))
        variant.set("precursor_true_units", "1")
        variant.set("precursor_true_tolerance",
                    java_decimal_format(_FRAGGER_DF, min(10.0, float(params.get("precursor_true_tolerance")))))
        if int(params.get("precursor_mass_units")) == 1:
            variant.set("precursor_mass_lower",
                        java_decimal_format(_FRAGGER_DF, max(-10.0, float(params.get("precursor_mass_lower")))))
            variant.set("precursor_mass_upper",
                        java_decimal_format(_FRAGGER_DF, min(10.0, float(params.get("precursor_mass_upper")))))
        else:
            variant.set("precursor_mass_units", "1")
            variant.set("precursor_mass_lower", java_decimal_format(_FRAGGER_DF, -10.0))
            variant.set("precursor_mass_upper", java_decimal_format(_FRAGGER_DF, 10.0))
    elif data_type == "GPF-DIA":
        variant.set("data_type", "2")
        variant.set("output_report_topN", str(topn_dia1))
    elif data_type == "DDA+":
        variant.set("data_type", "3")
        variant.set("output_report_topN", str(topn_dda_plus))
    else:
        raise ValueError(f"Unknown data type: {data_type}")
    return variant


def msfragger_num_slices(params: Dict[str, str], all_params: Dict[str, Dict[str, str]]) -> int:
    """Split-database slice count FragPipe uses: the misc.slice-db spinner after the workflow is
    loaded (reset to 1 by loadDefaults, FP/tabs/TabMsfragger.java:2095)."""
    form = msfragger_form(params, all_params)
    slices = form.widgets["misc.slice-db"].value
    if slices < 1:
        raise ValueError(f"msfragger.misc.slice-db must be at least 1, found {slices}")
    return slices


def msfragger_params_file_name(params_type: str, data_types: List[str]) -> str:
    """File name FragPipe saves a data type's MSFragger params under (FP/cmd/CmdMsfragger.java:444-447):

        Path savedDdaParamsPath = (hasDia || hasGpfDia || hasDiaLib) ? wd.resolve("fragger_dda.params") : wd.resolve("fragger.params");

    hasDia/hasGpfDia/hasDiaLib are TabWorkflow.hasDataType over the manifest (FP/FragpipeRun.java:1474-1478).
    MSFragger is run with that file (CmdMsfragger.java:540-552) and MSBooster reads it (CmdMSBooster.java:134-139).
    """
    if params_type == "DDA" and {"DIA", "GPF-DIA", "DIA-Lib"} & set(data_types):
        return "fragger_dda.params"
    return MSFRAGGER_PARAMS_FILES[params_type]


def generate_msfragger_configs(
    params: Dict[str, str], all_params: Dict[str, Dict[str, str]]
) -> Dict[str, str]:
    """All four MSFragger params files {filename: content}, byte for byte as FragPipe 24.0 writes
    fragger[_dda].params, fragger_dia.params, fragger_gpfdia.params and fragger_dda_plus.params.

    database_name is the workflow's database.db-path and num_threads its workflow.threads (the
    MSFRAGGER module replaces both values at run time, as FragPipe does from its run inputs).
    write_mzbin_all follows tab-run.write_sub_mzml (FP/FragpipeRun.java:1481).
    """
    form = msfragger_form(params, all_params)

    template = FraggerProps.from_template()
    num_threads = all_params.get("workflow", {}).get("threads", template.get("num_threads"))
    collected = _fragger_params_from_form(form, num_threads)
    db_path = all_params.get("database", {}).get("db-path", template.get("database_name"))
    collected.set("database_name", db_path)
    collected.set("decoy_prefix", get_decoy_tag(all_params))
    write_sub_mzml = _java_parse_boolean(all_params.get("tab-run", {}).get("write_sub_mzml", "false"))
    collected.set("write_mzbin_all", "1" if write_sub_mzml else "0")

    topn_dia1 = form.widgets["output_report_topN_dia1"].value
    topn_dda_plus = form.widgets["output_report_topN_dda_plus"].value
    return {
        filename: _fragger_data_type_variant(collected, data_type, topn_dia1, topn_dda_plus).write()
        for data_type, filename in MSFRAGGER_PARAMS_FILES.items()
    }


def _params_value(text: str, key: str) -> str:
    match = re.search(r"(?m)^" + re.escape(key) + r" = ([^\t\n]*)", text)
    if match is None:
        raise ValueError(f"MSFragger params have no {key}")
    return match.group(1).strip()


def msfragger_search_settings(
    all_params: Dict[str, Dict[str, str]], data_types: List[str]
) -> Dict[str, Any]:
    """
    Per params set calibrate_mass, plus FragPipe's split database checks (CmdMsfragger.java:337-357).

    With more than one slice FragPipe stops the run when the base (DDA) params localize mass
    shifts, or when multiple rank reporting meets the split search: DIA/GPF-DIA/DIA-Lib files
    without DIA-Umpire or diaTracer and output_report_topN_dia1 > 1, or DDA+ files and
    output_report_topN_dda_plus > 1, or output_report_topN > 1 whatever the data types.
    DIA-Umpire/diaTracer pseudo spectra are searched as DDA, so only the DDA check covers them.

    @param all_params Parsed workflow after the manifest adjustments
    @param data_types Data types of the manifest files
    @return {"calibrate_mass": {params type: value}}
    @raises ValueError when FragPipe would refuse the split database search
    """
    msfragger = all_params.get("msfragger", {})
    configs = generate_msfragger_configs(msfragger, all_params)
    texts = {t: configs[f] for t, f in MSFRAGGER_PARAMS_FILES.items()}
    if msfragger_num_slices(msfragger, all_params) > 1:
        if _params_value(texts["DDA"], "localize_delta_mass") == "1":
            raise ValueError("Split database is incompatible with localize mass shift (msfragger.misc.slice-db > 1 and localize_delta_mass = 1)")
        has = set(data_types).__contains__
        dia_pseudo = get_run_flag(all_params, "diaumpire") or get_run_flag(all_params, "diatracer")
        top_n = {t: int(_params_value(texts[t], "output_report_topN")) for t in ("DDA", "DIA", "DDA+")}
        bad = (
            ((has("DIA") or has("GPF-DIA") or has("DIA-Lib")) and not dia_pseudo and top_n["DIA"] > 1)
            or (has("DDA+") and top_n["DDA+"] > 1)
            or top_n["DDA"] > 1
        )
        if bad:
            raise ValueError(
                "The split database option is not compatible with multiple rank reporting, which is typical "
                f"for DIA, GPF-DIA, DIA-Lib, or DDA+ data types (output_report_topN {top_n}); set msfragger.misc.slice-db to 1"
            )
    return {"calibrate_mass": {t: _params_value(text, "calibrate_mass") for t, text in texts.items()}}


def generate_msfragger_config(
    params: Dict[str, str], all_params: Dict[str, Dict[str, str]]
) -> str:
    """The DDA MSFragger params file (fragger.params)."""
    return generate_msfragger_configs(params, all_params)["fragger.params"]


# TMT-Integrator config keys (TmtiConfProps.PROPS) with the TmtiPanel widget kind and its initial
# value (TmtiPanel.java:320-640). formToConfig (TmtiPanel.java:1177-1224) writes each widget's
# string value for the PROPS keys, sorted (TreeMap), plus output and is_tmt_35; TmtiConfig.write
# prints "  key: value" under "tmtintegrator:" (TmtiConfig.java:332-338).
#   df2/df3: UiSpinnerDouble with DecimalFormat "#.##" / "#.###"; int: UiSpinnerInt;
#   bool: UiCheck / UiRadio; text: UiText; combo: config values of the TmtiConfProps combo.
# abn_type is not collected: its form entry (TmtiPanel.java:389) is never added to the panel, so the
# FragPipe 24.0 golden tmt-integrator-conf.yml has no abn_type line.
TMTI_WIDGETS = {
    "ref_tag": ("text", "Bridge"),
    "ref_d_tag": ("text", "Pool"),
    "groupby": ("combo", "0", ("0", "1", "2", "3", "4", "-1"),
                ("Gene level", "Protein", "Peptide sequence", "Multiple PTM sites", "Single PTM site", "All")),
    "prot_norm": ("combo", "0", ("0", "1", "2", "-1"),
                  ("None", "MD (median centering)", "GN (median centering + variance scaling)", "All")),
    "log2transformed": ("bool", "true"),
    "add_Ref": ("combo", "-1", ("-1", "1"), ("Reference sample", "Virtual")),
    "unique_gene": ("combo", "0", ("0", "1", "2"),
                    ("Keep all PSMs", "No PSMs mapping to multiple genes in dataset",
                     "No PSMs mapping to multiple genes in FASTA")),
    "min_pep_prob": ("df3", "0.9"),
    "min_purity": ("df3", "0.5"),
    "min_percent": ("df3", "0.05"),
    "min_resolution": ("int", "45000"),
    "min_snr": ("df2", "1000"),
    "min_ntt": ("int", "0"),
    "prot_exclude": ("text", "none"),
    "unique_pep": ("combo", "false", ("false", "true"), ("Unique+Razor", "Unique only")),
    "aggregation_method": ("combo", "0", ("0", "1"), ("Median", "Weighted")),
    "best_psm": ("bool", "true"),
    "psm_norm": ("bool", "false"),
    "outlier_removal": ("bool", "true"),
    "allow_overlabel": ("bool", "true"),
    "allow_unlabeled": ("bool", "true"),
    # Radio "MS1" of a ButtonGroup with an unnamed "MS2" (TmtiPanel.java:533-540): false cannot
    # clear it and no workflow key selects MS2, so FragPipe always writes ms1_int: true
    "ms1_int": ("radio", "true"),
    "print_RefInt": ("bool", "false"),
    "max_pep_prob_thres": ("df2", "0"),
    "min_site_prob": ("df2", "-1"),
    "mod_tag": ("text", "none"),
    "glyco_qval": ("df2", "-1"),
    "use_glycan_composition": ("bool", "false"),
}

# Placeholder for the TMT-Integrator report directory (formToConfig "output", the tmt-report
# directory of the FragPipe work dir, CmdTmtIntegrator.java:149-170); the TMTINTEGRATOR module
# replaces it with the task's absolute tmt-report path.
TMTI_OUTPUT_PLACEHOLDER = "@TMT_REPORT_DIR@"

# QuantLabelAnnotation.disallowedPattern (QuantLabelAnnotation.java:25, :62-64)
SAMPLE_NAME_DISALLOWED = re.compile(r"[^0-9a-zA-Z_-]")


def _tmti_widget_value(key: str, raw) -> str:
    """String a TmtiPanel widget returns after loading `raw` (None: key not in the workflow)."""
    spec = TMTI_WIDGETS[key]
    kind, default = spec[0], spec[1]
    value = default if raw is None else raw.strip()
    if kind == "bool":
        return "true" if value.lower() == "true" else "false"
    if kind == "radio":
        # UiRadio.fromString: setSelected(Boolean.valueOf(s)) (UiRadio.java:40-42); false is a no-op
        return "true" if value.lower() == "true" else default
    if kind == "int":
        return str(int(value))
    if kind in ("df2", "df3"):
        return java_decimal_format("#.##" if kind == "df2" else "#.###", float(value))
    if kind == "combo":
        conf_values, ui_values = spec[2], spec[3]
        if value in ui_values:
            return conf_values[ui_values.index(value)]
        if value in conf_values:
            return value
        raise ValueError(f"No matching {key} value: tmtintegrator.{key}={value}")
    return value


def generate_tmtintegrator_config(params: Dict[str, str], is_second_unmod_run: bool = False) -> str:
    """tmt-integrator-conf.yml as FragPipe 24.0 writes it (TmtiPanel.formToConfig, TmtiConfig.write).

    is_second_unmod_run: FPOP TMT's second run (tmt-integrator-conf-unmod.yml), which forces
    mod_tag none and min_site_prob -1 (TmtiPanel.java:1188-1199).

    channel_num is the label's reagent count with its label_masses; ref_tag is sanitised with
    QuantLabelAnnotation.unifyAnnotationSampleName; an empty mod_tag becomes none; output is
    TMTI_OUTPUT_PLACEHOLDER; is_tmt_35 is true only for TMT-35. No path/memory keys (the JAR
    and -Xmx are on the java command line, CmdTmtIntegrator.java:288-296).
    """
    label = tmt_quant_label(params)
    conf = {key: _tmti_widget_value(key, params.get(key)) for key in TMTI_WIDGETS}
    conf["ref_tag"] = SAMPLE_NAME_DISALLOWED.sub("", conf["ref_tag"])
    if not conf["mod_tag"]:
        conf["mod_tag"] = "none"
    conf["channel_num"] = str(QUANT_LABELS[label][1])
    conf["label_masses"] = QUANT_LABELS[label][2]
    conf["output"] = TMTI_OUTPUT_PLACEHOLDER
    conf["is_tmt_35"] = "true" if label == "TMT-35" else "false"
    if is_second_unmod_run:
        conf["mod_tag"] = "none"
        conf["min_site_prob"] = "-1"
    # TmtiConfig.write (TmtiConfig.java:332-338): "tmtintegrator:\n", then "  key: value\n" per entry
    return "".join(line + "\n" for line in ["tmtintegrator:"] + [f"  {key}: {conf[key]}" for key in sorted(conf)])


def _generate_tmt_ionquant_ms1_args(all_params: Dict[str, Dict[str, str]]) -> str:
    """IonQuant MS1 pass of the TMT workflow (FragpipeRun.java:1978-1995): the short configure
    with no label-free panel map and generateSiteReport = generateMSstatsInput = false, so
    CmdIonquant.java:226-254's fixed options."""
    head = _ionquant_head(True, False, "20.0", "MS2", "tmt10", False, False, all_params)
    return " ".join(head + ["--minexps 1", IONQUANT_ISOBARIC_FIXED])


def _generate_tmt_ionquant_isobaric_args(
    tmti_params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]],
) -> str:
    """IonQuant isobaric pass of the TMT workflow (FragpipeRun.java:2005-2028): performMS1Quant
    false, performIsobaricQuant true, isoTol = TMT-Integrator tolerance (int -> float),
    isoLevel = quant level, isoType = label name, writeExperimentAnnotation true. The module adds
    one --annotation per group before --minexps (CmdIonquant.java:204-214)."""
    isotol = _java_float_to_string(str(int(tmti_params.get("tolerance", "20"))))
    head = _ionquant_head(
        False, True, isotol, tmt_quant_level(tmti_params), tmt_quant_label(tmti_params),
        False, False, all_params,
    )
    return " ".join(head + ["--minexps 1", IONQUANT_ISOBARIC_FIXED])


# FragPipe-GUI/.../tools/ptmshepherd/shepherd_default_open.config (FragPipe 24.0), loaded with
# java.util.Properties as PtmshepherdParams' defaults (PtmshepherdParams.java:57-64).
# Properties keeps a value's trailing blanks ("annotation_tol = 0.01 ").
SHEPHERD_DEFAULT_OPEN_CONFIG = {
    "histo_bindivs": "5000",
    "histo_smoothbins": "2",
    "histo_normalizeTo": "psms",
    "peakpicking_promRatio": "0.3",
    "peakpicking_mass_units": "0",
    "peakpicking_width": "0.002",
    "peakpicking_topN": "500",
    "peakpicking_minPsm": "10",
    "precursor_mass_units": "0",
    "precursor_tol": "0.01",
    "annotation_tol": "0.01 ",
    "varmod_masses": "",
    "mass_offsets": "",
    "isotope_states": "",
    "localization_background": "4",
    "localization_allowed_res": "all",
    "spectra_ppmtol": "20.0",
    "compare_betweenRuns": "false",
    "output_extended": "false",
    "annotation_file": "unimod",
    "iontype_a": "0",
    "iontype_x": "0",
    "iontype_b": "1",
    "iontype_y": "1",
    "iontype_c": "0",
    "iontype_z": "0",
    "spectra_maxfragcharge": "2",
    "run_glyco_mode": "false",
    "cap_y_ions": "",
    "diag_ions": "",
    "remainder_masses": "",
    "assign_glycans": "false",
}

# PtmshepherdPanel and PTMSGlycanAssignPanel components (names without the "ptmshepherd." prefix) in the
# state a loaded workflow starts from: the constructor values, and for PtmshepherdPanel the keys of
# shepherd_default_open.config that initMore() loads over them (loadDefaults(1, SearchTypeProp.open),
# PtmshepherdPanel.java:787, 255-299: annotation_file "unimod" selects the Unimod radio and clears the text,
# CONV_TO_GUI turns the 0/1 ion types into checkbox values, PtmshepherdPanel.java:160-168).
# FragPipe writes every component (toPtmsParamsMap() and getGlycanAssignParams() are super.toMap(), i.e.
# SwingUtils.valueGet of each named component, SwingUtils.java:775-833), so an absent workflow key keeps
# this value, and a workflow key that names no component (iterloc_*, histo_bindivs, ...) is dropped.
PTMS_PANEL_WIDGETS = {
    "run-shepherd": ("check", False),  # PtmshepherdPanel.java:325 new UiCheck("Run PTM-Shepherd", null, false)
    "output_extended": ("check", False),  # :344-346 new UiCheck("Extended output", null, false)
    "run_diagextract_mode": ("check", False),  # :371
    "diag_ions": ("text", ""),  # :376 (default config "")
    "cap_y_ions": ("text", ""),  # :380 (default config "")
    "remainder_masses": ("text", ""),  # :384 (default config "")
    "run_diagmine_mode": ("check", False),  # :432
    "diagmine_minPeps": ("int", 25),  # :437-438 new UiSpinnerInt(25, ...)
    "diagmine_diagMinSpecDiff": ("double", 25.0, "00.#"),  # :446-447
    "diagmine_diagMinFoldChange": ("double", 3.0, "0.0#"),  # :452-453
    "diagmine_pepMinSpecDiff": ("double", 25.0, "00.#"),  # :469-470
    "diagmine_pepMinFoldChange": ("double", 3.0, "0.0#"),  # :475-476
    "diagmine_fragMinSpecDiff": ("double", 15.0, "00.#"),  # :492-493
    "diagmine_minIonsPerSpec": ("int", 2),  # :499-500 new UiSpinnerInt(2, ...)
    "diagmine_fragMinFoldChange": ("double", 3.0, "0.0#"),  # :503-504
    "diagmine_fragMinPropensity": ("double", 12.5, "00.#"),  # :510-511
    "histo_smoothbins": ("int", 2),  # :539-540 (default config 2)
    "peakpicking_promRatio": ("double", 0.3, "0.#"),  # :544-545 (default config 0.3)
    "peakpicking_width": ("double", 0.002, "0.####"),  # :550-551 (default config 0.002)
    "peakpicking_mass_units": ("combo", "Da", ("PPM", "Da")),  # :554 createUiCombo(MassTolUnits.values()); config 0 = Da
    "precursor_tol": ("double", 0.01, "0.###"),  # :559-560 (default config 0.01)
    "precursor_mass_units": ("combo", "Da", ("PPM", "Da")),  # :563 (default config 0 = Da)
    "peakpicking_minPsm": ("int", 10),  # :566-567 (default config 10)
    "normalization-psms": ("radio", True, "normalization"),  # :571-572 new JRadioButton("PSMs", true)
    "normalization-scans": ("radio", False, "normalization"),  # :573-574
    "annotation-unimod": ("radio", True, "annotation"),  # :579-580; loadDefaults selects it (:262-266)
    "annotation-common": ("radio", False, "annotation"),  # :581-582
    "annotation-glyco": ("radio", False, "annotation"),  # :583-584
    "annotation-custom": ("radio", False, "annotation"),  # :585-586
    "annotation_file": ("text", ""),  # :593; loadDefaults clears the "unimod" keyword (:266)
    "annotate_assigned_mods": ("check", False),  # :609-610
    "varmod_masses": ("text", ""),  # :640 (default config "")
    "annotation_tol": ("double", 0.01, "0.###"),  # :663-664 (default config "0.01 ")
    "iontype_a": ("check", False),  # :688-693, default config 0/0/1/1/0/0 through CONV_TO_GUI
    "iontype_x": ("check", False),
    "iontype_b": ("check", True),
    "iontype_y": ("check", True),
    "iontype_c": ("check", False),
    "iontype_z": ("check", False),
    "localization_allowed_res": ("text", "all"),  # :696 PROP_restrict_loc (default config "all")
    "use_msfragger_localization": ("check", False),  # :701-702
    "spectra_ppmtol": ("double", 20.0, "0.#"),  # :720-721 (default config 20.0)
    "spectra_maxfragcharge": ("int", 2),  # :726 spinnerInt(2, ...) (default config 2)
    "spectra_maxPrecursorCharge": ("int", 4),  # :732 spinnerInt(4, ...)
    "spectra_condPeaks": ("int", 150),  # :738 spinnerInt(150, ...)
    "spectra_condRatio": ("double", 0.0001, "0.#######"),  # :742-743
}
# PTMSGlycanAssignPanel.java:151-218; reach shepherd.config only through getGlycanAssignParams(), i.e. only
# when glycan assignment runs (FragpipeRun.java:2137-2138).
PTMS_GLYCAN_PANEL_WIDGETS = {
    "run_glyco_mode": ("check", False),  # :150-151 createUiCheck("Assign Glycans with FDR", false)
    "adv_params": ("check", False),  # :153-154
    "glyco_fdr": ("double", 0.01, "0.00#"),  # :157-158
    "glyco_ppm_tol": ("double", 50.0, "0.#"),  # :162-163
    "glyco_isotope_min": ("int", -1),  # :168 new UiSpinnerInt(-1, -2, 0, 1, 3)
    "glyco_isotope_max": ("int", 3),  # :171 new UiSpinnerInt(3, 0, 4, 1, 3)
    "glycodatabase": ("text", ""),  # :175
    "decoy_type": ("int", 1),  # :181 new UiSpinnerInt(1, 0, 3, 1, 1)
    "print_decoys": ("check", False),  # :188
    "remove_glycan_delta_mass": ("check", False),  # :192
    "print_full_glyco_params": ("check", False),  # :196
    "prob_mass": ("double", 0.5, "0.0#"),  # :199-200
    "n_glyco": ("check", False),  # :203
    "glyco_lda": ("check", False),  # :209
    "use_glycan_fragment_probs": ("check", False),  # :211
    "shuffle_decoy_intensities": ("check", False),  # :213
    "glyco_lda_features": ("text", ""),  # :216
}
# ui.translations:68-71: the workflow's 0/1 mass units are the UI's Da/PPM
PTMS_UI_TRANSLATIONS = {
    "peakpicking_mass_units": [("PPM", "1"), ("Da", "0")],
    "precursor_mass_units": [("PPM", "1"), ("Da", "0")],
}

# PtmshepherdPanel CONV_TO_FILE (PtmshepherdPanel.java:170-178): ion-type checkboxes 1/0, mass units 0/1
PTMS_IONTYPE_KEYS = {"iontype_a", "iontype_x", "iontype_b", "iontype_y", "iontype_c", "iontype_z"}
PTMS_MASS_UNITS_TO_FILE = {"PPM": "1", "Da": "0"}  # MassTolUnits.valueInParamsFile()

# Token PTMSHEPHERD fills with the staged copy of the custom annotation file
PTMS_ANNOTATION_FILE_TOKEN = "@ANNOTATION_FILE@"


def _ptms_panel_state(
    widgets: Dict[str, tuple], params: Dict[str, str], all_params: Dict[str, Dict[str, str]]
) -> Dict[str, str]:
    """valueGet of every component of a PTM-Shepherd panel after the workflow is loaded into it.

    FragpipeCacheUtils.tabsLoad translates the workflow values (FragpipeCacheUtils.java:182-206) and
    SwingUtils.valuesSet sets the named components in the map's HashMap order (SwingUtils.java:859-868):
    spinners parse the text (UiSpinnerDouble.fromString / UiSpinnerInt.fromString) and print it with their
    DecimalFormat / Integer.toString, checkboxes and radio buttons take Boolean.parseBoolean, a combo takes
    the first item for a value it lacks (UiCombo.java:35-41), text fields are trimmed on read (UiText.java:53).
    A radio button set to true deselects the rest of its ButtonGroup; one set to false changes nothing.
    """
    # selected button per ButtonGroup
    groups = {spec[2]: name for name, spec in widgets.items() if spec[0] == "radio" and spec[1]}
    state = {}
    for name, spec in widgets.items():
        kind = spec[0]
        if kind == "radio":
            continue
        state[name] = _Widget(
            {"double": "spinD", "int": "spinI", "check": "check", "text": "text", "combo": "combo"}[kind],
            spec[1],
            options=spec[2] if kind == "combo" else None,
            fmt=spec[2] if kind == "double" else None,
        )
    capacity = _java_hashmap_capacity(len(_headless_ui_keys(all_params)))
    names = [k for k in params if k in widgets]
    bucket = {k: _java_hashmap_bucket(f"ptmshepherd.{k}", capacity) for k in names}
    radios = [k for k in names if widgets[k][0] == "radio" and _java_parse_boolean(params[k])]
    for a in radios:
        for b in radios:
            if a < b and widgets[a][2] == widgets[b][2] and bucket[a] == bucket[b]:
                raise ValueError(
                    f"ptmshepherd.{a} and ptmshepherd.{b} are both true and share a HashMap bucket "
                    f"(capacity {capacity}); which radio button FragPipe selects is not defined by the keys"
                )
    for name in sorted(names, key=lambda k: bucket[k]):
        value = params[name]
        if widgets[name][0] == "radio":
            if _java_parse_boolean(value):
                groups[widgets[name][2]] = name
            continue
        translations = PTMS_UI_TRANSLATIONS.get(name, [])
        if not any(in_ui == value for in_ui, _ in translations):
            value = next((in_ui for in_ui, in_conf in translations if in_conf.lower() == value.lower()), value)
        state[name].set(value)
    result = {name: widget.get() for name, widget in state.items()}
    for name, spec in widgets.items():
        if spec[0] == "radio":
            result[name] = "true" if groups[spec[2]] == name else "false"
    return result


def ptms_annotation_file(params: Dict[str, str], all_params: Dict[str, Dict[str, str]]) -> str:
    """The custom annotation file path PTM-Shepherd reads, "" when annotation_file is a keyword or empty.

    PtmshepherdPanel.toPtmsParamsMap() (PtmshepherdPanel.java:229-232) copies the text field
    (ptmshepherd.annotation_file) into the config when the "annotation-glyco" or "annotation-custom"
    radio button is on; "annotation-unimod" / "annotation-common" write the keywords "unimod" /
    "common" (PtmshepherdPanel.java:225-228). loadDefaults() treats the values "unimod", "common"
    and "glyco" as keywords and any other non-blank value as the custom file
    (PtmshepherdPanel.java:264-278).
    """
    state = _ptms_panel_state(PTMS_PANEL_WIDGETS, params, all_params)
    path = state["annotation_file"]
    on = state["annotation-glyco"] == "true" or state["annotation-custom"] == "true"
    if not on or path.lower() in ("", "unimod", "common", "glyco"):
        return ""
    return path


def ptms_mass_offsets(msfragger: Dict[str, str]) -> str:
    """tabMsf.getMassOffsets() as FragpipeRun.java:2132 passes it to PTM-Shepherd (TabMsfragger.java:1972-1978).

    With detailed offsets: the detailed list's masses (MassOffset(String): first "<mass>(" match,
    parsed as a float, 0 when absent), each Double.toString'd, joined with "/".
    Otherwise uiTextMassOffsets, a text field FragPipe creates with "0" (TabMsfragger.java:1258) and
    never binds to a workflow key (msfragger.mass_offsets fills epMassOffsets), so it is "0".
    """
    if msfragger.get("use_detailed_offsets", "false").lower() == "true":
        detailed = msfragger.get("mass_offsets_detailed", "")
        masses = []
        if detailed:
            for split in detailed.split(";"):
                match = re.search(r"([\d.-]+)\(", split)
                value = struct.unpack("f", struct.pack("f", float(match.group(1))))[0] if match else 0.0
                masses.append(_java_double_to_string(repr(value)))
        return "/".join(masses)
    return "0"


def generate_ptmshepherd_config(
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]],
) -> str:
    """The property block of FragPipe's shepherd.config (PtmshepherdParams.createConfig()).

    props = shepherd_default_open.config, overlaid with the run-time map FragpipeRun.java:2127-2148
    builds (PtmshepherdPanel.toPtmsParamsMap(), threads, mass_offsets, msfragger_massdiff_to_varmod,
    getGlycanAssignParams() + glyco_only_mode, isotope_error), then written sorted by key as
    `key = value` with empty values kept and ui.* dropped (PtmshepherdParams.java:93-99).

    The module fills the run-time tokens (@THREADS@, @TOOLS_DIR@) and writes the `database =` and
    `dataset =` lines plus a blank line above this block, as createConfig() does.
    """
    # PtmshepherdPanel.toPtmsParamsMap() (PtmshepherdPanel.java:216-245): every panel component.
    map1 = _ptms_panel_state(PTMS_PANEL_WIDGETS, params, all_params)
    map2 = {}
    for key, value in map1.items():
        is_true = value.lower() == "true"
        if key.startswith("annotation_file"):
            pass  # radio-button text field, consumed by the branches below
        elif key.startswith("annotation-unimod") and is_true:
            map2["annotation_file"] = "unimod"
        elif key.startswith("annotation-common") and is_true:
            map2["annotation_file"] = "common"
        elif (key.startswith("annotation-glyco") or key.startswith("annotation-custom")) and is_true:
            # the custom-file text field (empty in every bundled workflow); a file path is staged by
            # PTMSHEPHERD, so the config names the staged copy through a token
            map2["annotation_file"] = (
                PTMS_ANNOTATION_FILE_TOKEN if ptms_annotation_file(params, all_params) else map1["annotation_file"]
            )
        else:
            map2[key] = value
        if key.startswith("normalization-scans") and is_true:
            map2["histo_normalizeTo"] = "scans"
        elif key.startswith("normalization-psms") and is_true:
            map2["histo_normalizeTo"] = "psms"
    additional = {}
    for key, value in map2.items():
        if key in PTMS_IONTYPE_KEYS:
            value = "1" if value.lower() == "true" else "0"
        elif key in PTMS_UI_TRANSLATIONS:
            value = PTMS_MASS_UNITS_TO_FILE[value]
        additional[key] = value

    # FragpipeRun.java:2128-2147
    additional["threads"] = THREADS_TOKEN
    msfragger = msfragger_values(all_params)
    mass_offsets = ptms_mass_offsets(msfragger)
    if mass_offsets.strip():
        additional["mass_offsets"] = mass_offsets
    additional["msfragger_massdiff_to_varmod"] = msfragger["mass_diff_to_variable_mod"]
    glycan = _ptms_panel_state(PTMS_GLYCAN_PANEL_WIDGETS, params, all_params)
    if glycan["run_glyco_mode"] == "true":
        # PTMSGlycanAssignPanel.getGlycanAssignParams() (PTMSGlycanAssignPanel.java:126-136)
        additional.update(glycan)
        glycan_dir = f"{TOOLS_DIR_TOKEN}/Glycan_Databases"
        additional["glyco_residue_list"] = f"{glycan_dir}/glycan_residues.txt"
        additional["glyco_mod_list"] = f"{glycan_dir}/glycan_mods.txt"
        additional["glyco_oxonium_list"] = f"{glycan_dir}/oxonium_ion_list.txt"
        if map1["run-shepherd"] != "true":
            additional["glyco_only_mode"] = "true"
    # tabMsf.getUiTextIsoErr().getNonGhostText(), put when not blank (FragpipeRun.java:2143-2145)
    isotope_error = msfragger["isotope_error"]
    if isotope_error.strip():
        additional["isotope_error"] = isotope_error

    props = dict(SHEPHERD_DEFAULT_OPEN_CONFIG)
    props.update(additional)
    return "\n".join(
        f"{key} = {value}" for key, value in sorted(props.items()) if not key.startswith("ui.")
    )


def _offsets_as_mods(msfragger: Dict[str, str]) -> list:
    """TabMsfragger.getOffsetsAsMods() (TabMsfragger.java:2252-2272) as [(mass, sites), ...].

    Detailed mode: one mod per detailed offset with its allowed residues.
    Regular mode: getMassOffsetSet() with the restrict_deltamass_to residues ("-" removed,
    "all" = the 20 amino acids).
    """
    if _uses_detailed_offsets(msfragger):
        return _detailed_mass_offsets(msfragger)
    residues = msfragger.get("restrict_deltamass_to", "").replace("-", "")
    if residues == "all":
        residues = _ALL_RESIDUES
    return [(mass, residues) for mass in _mass_offset_set(msfragger)]


def _get_ptmprophet_mods_from_msfragger(
    all_params: Dict[str, Dict[str, str]],
) -> str:
    """Extract MSFragger variable mods and mass offsets as PTMProphet mod string.

    Replicates PtmProphetPanel.getMSFraggerMods() (PtmProphetPanel.java:145-206):
    - Parses MSFragger table.var-mods (format: mass,sites,enabled,max_mods; ...)
    - Merges getOffsetsAsMods() (detailed offsets or restrict_deltamass_to sites) into var mods of
      equal mass (within 0.0001)
    - Converts MSFragger site notation to PTMProphet notation:
        n^ or [^ -> n (protein N-terminus)
        c^ or ]^ -> c (protein C-terminus)
        Plain n or [ (no ^) -> removed (not protein terminus)
        Plain c or ] (no ^) -> removed
        * or empty -> ACDEFGHIKLMNPQRSTVWY (all amino acids)
        Remaining ^ characters -> removed
    - Formats as SITES:MASS joined by commas (e.g., STY:79.96633,M:15.9949)

    Returns:
        PTMProphet mod string, or empty string if no mods found.
    """
    msfragger = msfragger_values(all_params)

    # Enabled variable mods (TabMsfragger.getVarModsTable()) as [float32 mass, sites]
    mods = []
    for entry in msfragger.get("table.var-mods", "").split(";"):
        parts = [p.strip() for p in entry.strip().split(",")]
        if len(parts) < 4:
            continue
        if parts[2].lower() == "true":
            mods.append([_float32(parts[0]), parts[1]])

    # Merge offsets into var mods of equal mass (PtmProphetPanel.java:150-175): the merged
    # mod replaces the var mod and moves to the end of the list.
    for offset_mass, offset_sites in _offsets_as_mods(msfragger):
        if offset_mass == 0:
            continue
        for mod in mods:
            if abs(_float32(mod[0] - offset_mass)) < 0.0001:
                merged = mod[1] + "".join(ch for ch in offset_sites if ch not in mod[1])
                mods.remove(mod)
                mods.append([mod[0], merged])
                break
        else:
            mods.append([offset_mass, offset_sites])

    # Convert site notation to PTMProphet format
    # Replicates PtmProphetPanel.java:177-206
    mod_strings = []
    for mass, sites in mods:
        if sites == "":
            sites = "ACDEFGHIKLMNPQRSTVWY"
        if "*" in sites:
            sites = sites.replace("*", "ACDEFGHIKLMNPQRSTVWY")

        # Handle terminal characters
        # n^ or [^ -> n (protein N-terminus)
        if "n^" in sites or "[^" in sites:
            sites = sites.replace("n^", "n")
            sites = sites.replace("[^", "n")
        else:
            sites = sites.replace("n", "")
            sites = sites.replace("[", "")

        # c^ or ]^ -> c (protein C-terminus)
        if "c^" in sites or "]^" in sites:
            sites = sites.replace("c^", "c")
            sites = sites.replace("]^", "c")
        else:
            sites = sites.replace("c", "")
            sites = sites.replace("]", "")

        # Remove remaining ^ characters
        sites = sites.replace("^", "")

        if mass != 0:
            mod_strings.append(f"{sites}:{_java_float_to_string(repr(mass))}")

    return ",".join(mod_strings)


PTMPROPHET_DEFAULT_CMDLINE = (
    "NOSTACK KEEPOLD STATIC FRAGPPMTOL=10 EM=1 NIONS=b STY:79.966331,M:15.9949 MINPROB=0.5"
)  # Bundle.properties:46 ptmprophet.cmdline.default (PtmProphetPanel.java:100-108)

PTMPROPHET_TRANSLATE = re.compile(r"--([^\s]+)\s+(((?!--).)*)")  # CmdPtmProphet.java:51


def ptmprophet_translate_cmds(cmdline: str) -> str:
    """CmdPtmProphet.translateCmds (CmdPtmProphet.java:138-151): Philosopher "--x v" style to native."""
    parts = []
    for match in PTMPROPHET_TRANSLATE.finditer(cmdline.strip()):
        name, value = match.group(1), match.group(2)
        if value is None or not value.strip():
            parts.append(name.upper())
        elif name.lower() == "mods":
            parts.append(value.strip())
        else:
            parts.append(f"{name.upper()}={value.strip()}")
    return " ".join(parts)


def generate_ptmprophet_config(
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]],
) -> str:
    """PTMProphet options before the module appends MAXTHREADS=1 and the file names.

    PtmProphetPanel.getCmdLineOpts() (PtmProphetPanel.java:241-247): with override-defaults the
    cmdline text verbatim (even empty); otherwise loadMSFraggerMods() (PtmProphetPanel.java:113-138):
    the Bundle.properties default with its first "<sites>:<mass>" token replaced by the MSFragger
    mods (FRAGPPMTOL=10 is kept). CmdPtmProphet.java:115-117 then translates "--" style options.
    """
    if params.get("override-defaults", "false").strip().lower() == "true":
        # the text field starts with the default line (PtmProphetPanel.java:86 text(getDefaults()))
        cmdline = params.get("cmdline", PTMPROPHET_DEFAULT_CMDLINE)
    else:
        splits = PTMPROPHET_DEFAULT_CMDLINE.split(" ")
        mod_pattern = re.compile(r"[^:]+:[\d+.-]+")
        new_mods = _get_ptmprophet_mods_from_msfragger(all_params)
        mod_index = next((i for i, t in enumerate(splits) if mod_pattern.search(t)), -1)
        if mod_index == -1:
            cmdline = PTMPROPHET_DEFAULT_CMDLINE + " " + new_mods
        else:
            splits[mod_index] = new_mods
            cmdline = " ".join(splits)
    if "--" in cmdline:
        cmdline = ptmprophet_translate_cmds(cmdline)
    return " ".join(cmdline.split())


def _java_number_string(shortest: str) -> str:
    """Java Float.toString / Double.toString layout of the shortest round-trip digits.

    Plain notation for 1e-3 <= |v| < 1e7 with at least one fractional digit, else
    computerized scientific notation d.dddE<n> (java.lang.Double#toString(double)).
    """
    value = Decimal(shortest)
    if value == 0:
        return "-0.0" if value.is_signed() else "0.0"
    sign = "-" if value.is_signed() else ""
    value = abs(value)
    if Decimal("1e-3") <= value < Decimal("1e7"):
        text = f"{value:f}"
        if "." in text:
            text = text.rstrip("0")
        return sign + (text + "0" if text.endswith(".") else text + ("" if "." in text else ".0"))
    _, digits, exponent = value.normalize().as_tuple()
    text = "".join(map(str, digits))
    return f"{sign}{text[0]}.{text[1:] or '0'}E{len(text) + exponent - 1}"




def _java_float_to_string(text: str) -> str:
    """String.valueOf(float) of a workflow value (Float.toString)."""
    value = struct.unpack("f", struct.pack("f", float(text)))[0]
    for digits in range(1, 10):
        candidate = f"{value:.{digits}g}"
        if struct.unpack("f", struct.pack("f", float(candidate)))[0] == value:
            return _java_number_string(candidate)
    raise ValueError(f"{text}: no float32 round-trip representation")


def _java_double_to_string(text: str) -> str:
    """String.valueOf(double) / String.format("%s", double) of a workflow value (Double.toString)."""
    return _java_number_string(repr(float(text)))


# diaTracer options in CmdDiaTracer.java:104-117 order: (workflow key, flag, kind, panel default).
# Defaults are the DiaTracerPanel widgets' initial values (DiaTracerPanel.java:121-129).
DIATRACER_OPTIONS = (
    ("write-intermediate-files", "--writeInter", "bool", "false"),
    ("delta-apex-im", "--deltaApexIM", "float", "0.01"),
    ("delta-apex-rt", "--deltaApexRT", "int", "3"),
    ("mass-defect-filter", "--massDefectFilter", "bool", "true"),
    ("mass-defect-offset", "--massDefectOffset", "float", "0.1"),
    ("corr-threshold", "--ms1MS2Corr", "float", "0.3"),
    ("rf-max", "--RFMax", "int", "500"),
)


def generate_diatracer_config(params: Dict[str, str]) -> str:
    """diaTracer options after --threadNum, as CmdDiaTracer.configure() writes them."""
    parts = []
    for key, flag, kind, default in DIATRACER_OPTIONS:
        raw = params.get(key, default).strip()
        if kind == "bool":
            value = "1" if raw.lower() == "true" else "0"
        elif kind == "float":
            value = java_float_to_string(float(raw))
        else:
            value = str(int(raw))
        parts.append(f"{flag} {value}")
    return " ".join(parts)


# Named DIA-Umpire panel params UmpirePanel.collect() overlays on the template: the ones with a
# component named "diaumpire.<name>" (UmpirePanel.java:224-244). UmpirePanel.paramNames also lists
# WindowType, WindowSize, SE.MinMSIntensity and SE.MinMSMSIntensity, but no component carries those
# names, so a workflow value never reaches them and the template value stays (UmpirePanel.java:417-424).
DIAUMPIRE_PANEL_PARAMS = (
    "RPmax", "RFmax", "CorrThreshold", "DeltaApex", "RTOverlap",
    "AdjustFragIntensity", "BoostComplementaryIon", "SE.MS1PPM", "SE.MS2PPM",
    "SE.NoMissedScan", "SE.EstimateBG",
    "SE.IsoPattern", "SE.MassDefectOffset", "SE.MassDefectFilter",
    "ExportPrecursorPeak", "SE.SN", "SE.MS2SN",
    "Q1", "Q2", "Q3",
)

# FragPipe 24.0's jar template FP/tools/umpire/diaumpire_se.params, copied verbatim, in the module's assets/;
# main() sets it from --diaumpire-template.
DIAUMPIRE_TEMPLATE: Optional[Path] = None

# Placeholder for the run's thread count; the DIAUMPIRE module puts task.cpus there.
DIAUMPIRE_THREADS = "THREADS"


def generate_diaumpire_config(params: Dict[str, str]) -> str:
    """umpire-se.params exactly as FragPipe writes it.

    UmpirePanel.collect() (UmpirePanel.java:397-428) loads the template, overlays the
    panel's named params with the workflow values and sets Thread; UmpireParams.load()
    (UmpireParams.java:180-238) and PropertiesUtils.writePropertiesContent()
    (PropertiesUtils.java:291-321) keep the template lines and write "name = value".
    """
    prop_regex = re.compile(r"^\s*([^=]+?)\s*=\s*(.+?)\s*$")
    out = []
    if DIAUMPIRE_TEMPLATE is None:
        raise RuntimeError("diaumpire_se.params is not set (--diaumpire-template)")
    for line in DIAUMPIRE_TEMPLATE.read_text(encoding="utf-8").split("\n"):
        hash_index = line.find("#")
        candidate = None
        comment = None
        if hash_index >= 0:
            if hash_index > 2 and line[:hash_index].find("=") > 2:
                candidate, comment = line[:hash_index], line[hash_index:]
        elif line.find("=") > 0:
            candidate = line
        match = prop_regex.match(candidate) if candidate is not None else None
        if match is None:
            out.append(line)
            continue
        name, value = match.group(1), match.group(2)
        if name == "Thread":
            value = DIAUMPIRE_THREADS
        elif name in DIAUMPIRE_PANEL_PARAMS and name in params:
            value = params[name]
        out.append(f"{name} = {value}" + (f"\t\t\t{comment}" if comment else ""))
    # writePropertiesContent println()s every template line (the template ends with an empty line);
    # drop the split's trailing "" after the last newline
    if out and out[-1] == "":
        out.pop()
    return "".join(f"{line}\n" for line in out)

# DIA-NN bundled with FragPipe 24.0: tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8, version
# "1.8.2 Beta 8" (Diann.java:56-58). Every NoteConfigDiann.compareVersion("1.9") and ("2.0")
# branch of CmdDiann/DiannPanel takes the "< 1.9" / "< 2.0" side for this binary.
DIANN_QUANTIFICATION_STRATEGIES = {
    # DiannPanel.getDiannQuantificationStrategy (< 1.9 branch, DiannPanel.java:482-497):
    # uiComboQuantificationStrategy index -> Collectors.toSet(). Case 3 is a HashSet whose
    # iteration order is bucket order: "--peak-center" hashes to bucket 2, "--no-ifs-removal"
    # to bucket 10 (String.hashCode ^ >>>16, & 15), so Java emits --peak-center first.
    "0": [],
    "1": ["--no-ifs-removal"],
    "2": ["--peak-center"],
    "3": ["--peak-center", "--no-ifs-removal"],
}

DIANN_LABEL_PATTERN = re.compile(r"([A-Znc*]+)([\d.+-]+)")  # CmdDiann.labelPattern (CmdDiann.java:79)


def _diann_parse_label(text: str) -> Dict[str, str]:
    """CmdDiann.parseLabel (CmdDiann.java:703-727): amino acid -> Float.toString(delta mass), TreeMap order."""
    labels: Dict[str, float] = {}
    for match in DIANN_LABEL_PATTERN.finditer(text.strip()):
        mass = struct.unpack("f", struct.pack("f", float(match.group(2))))[0]
        for aa in match.group(1):
            if aa == "*":
                for code in range(65, 91):
                    labels[chr(code)] = mass
                labels["n"] = mass
                labels["c"] = mass
            else:
                labels[aa] = mass
    return dict(sorted(labels.items()))


def _diann_plex_flags(light: str, medium: str, heavy: str) -> list:
    """CmdDiann.getPlexDiannFlags (CmdDiann.java:623-701)."""
    if not light:
        raise ValueError("diann.light must not be empty in plexDIA mode (CmdDiann.java:626-628)")
    light_map = _diann_parse_label(light)
    if not light_map:
        raise ValueError("diann.light must not be empty in plexDIA mode (CmdDiann.java:634-636)")
    other_maps = []
    for name, text in (("M", medium), ("H", heavy)):
        if not text:
            continue
        label_map = _diann_parse_label(text)
        if not label_map:
            raise ValueError(f"diann plexDIA label setting is not empty but could not parse it: {text}")
        if set(label_map) != set(light_map):
            raise ValueError("Light and medium/heavy labels must have the same amino acids (CmdDiann.java:643-655)")
        other_maps.append((name, label_map))

    def f32(value: float) -> float:
        return struct.unpack("f", struct.pack("f", value))[0]

    def channel(name: str, label_map: Dict[str, float]) -> str:
        # buildChannelString: Float subtraction, String.valueOf(float) (CmdDiann.java:689-701)
        deltas = ":".join(
            _java_float_to_string(repr(f32(mass - light_map[aa]))) for aa, mass in label_map.items()
        )
        return f"label,{name},{''.join(label_map)},{deltas}"

    flags = []
    for match in DIANN_LABEL_PATTERN.finditer(light.strip()):
        flags += ["--fixed-mod", f"label,{match.group(2)},{match.group(1)},label"]
    channels = [channel("L", light_map)] + [channel(name, label_map) for name, label_map in other_maps]
    flags += ["--channels", ";".join(channels), "--peak-translation", "--original-mods"]
    return flags


def _diann_flag(params: Dict[str, str], key: str, default: str) -> bool:
    return params.get(key, default).strip().lower() == "true"


def diann_settings(params: Dict[str, str]) -> Dict[str, Any]:
    """DIA-NN panel values as FragPipe reads them (DiannPanel.java:431-565, CmdDiann.java:100-140).

    Missing keys take the panel widgets' initial values (DiannPanel.java:170-400).
    """
    qvalue = _java_float_to_string(params.get("q-value", "0.01"))
    run_plex = _diann_flag(params, "run-dia-plex", "false")
    library = params.get("library", "").strip()
    mod_tag = params.get("mod-tag", "").strip()
    # generateMsstats(): isEnabledAndChecked, disabled while plexDIA is on (DiannPanel.java:344);
    # an external library turns it off (CmdDiann.java:137-140).
    generate_msstats = _diann_flag(params, "generate-msstats", "true") and not run_plex and not library
    # getReportLevels (DiannPanel.java:546-565)
    levels = [
        name
        for key, name, default in (
            ("gene-level-report", "gene", "false"),
            ("protein-level-report", "protein", "false"),
            ("peptide-level-report", "peptide", "true"),
            ("modified-peptide-level-report", "modified-peptide", "true"),
        )
        if _diann_flag(params, key, default)
    ]
    if _diann_flag(params, "site-level-report", "false") and mod_tag:
        levels += ["multi-site", "single-site"]
    return {
        "qvalue": qvalue,
        "run_specific_protein_qvalue": _diann_flag(params, "run-specific-protein-q-value", "false"),
        "mbr": _diann_flag(params, "mbr", "false"),
        "redo_protein_inference": _diann_flag(params, "redo-protein-inference", "false"),
        "unrelated_runs": _diann_flag(params, "unrelated-runs", "false"),
        "quantification_strategy": params.get("quantification-strategy", "3").strip(),
        "generate_msstats": generate_msstats,
        "library": "" if run_plex else library,
        "run_plex": run_plex,
        "light": params.get("light", "").strip(),
        "medium": params.get("medium", "").strip(),
        "heavy": params.get("heavy", "").strip(),
        "cmd_opts": params.get("cmd-opts", "").strip(),
        "mod_tag": mod_tag,
        "min_site_prob": _java_float_to_string(params.get("min-site-prob", "0.75")),
        "report_levels": ";".join(levels) if levels else "all",
    }


def generate_diann_config(params: Dict[str, str]) -> str:
    """DIA-NN flags of CmdDiann.configure after --qvalue and up to --matrices
    (CmdDiann.java:333-344), for the bundled DIA-NN 1.8.2 Beta 8.

    FragPipe 24.0's released jar does not pass --no-quant-files: the 24.0 tag (c2f256cb) adds it
    after --smart-profiling, but the shipped build's log does not show it and it was removed from
    the source in 91bb5c92f (2025-12-29). This follows the released 24.0 behaviour.

    The DIANNFRAGPIPE module adds the runtime arguments around them:
    --lib/--threads/--verbose 1/--out before, the protein-inference flags
    (CmdDiann.java:345-360, needs the staged FASTA) and generate_diann_post_args() after,
    then --cfg filelist_diann.txt--.
    """
    settings = diann_settings(params)
    qvalue = settings["qvalue"]
    args = ["--qvalue", qvalue]
    if settings["mbr"]:
        args.append("--reanalyse")
    if settings["run_specific_protein_qvalue"]:
        args += ["--matrix-spec-q", qvalue]
    args += ["--matrix-qvalue", qvalue, "--matrices"]
    return " ".join(args)


def generate_diann_post_args(params: Dict[str, str]) -> str:
    """DIA-NN flags of CmdDiann.configure from --smart-profiling to the extra options
    (CmdDiann.java:361-384), for the bundled DIA-NN 1.8.2 Beta 8."""
    settings = diann_settings(params)
    strategy = settings["quantification_strategy"]
    if strategy not in DIANN_QUANTIFICATION_STRATEGIES:
        raise ValueError(f"Unknown diann.quantification-strategy: {strategy!r} (DiannPanel.java:482-497)")
    args = ["--smart-profiling"] + DIANN_QUANTIFICATION_STRATEGIES[strategy]
    if settings["unrelated_runs"]:
        args += ["--individual-mass-acc", "--individual-windows"]
    if settings["generate_msstats"]:
        args.append("--report-lib-info")  # compareVersion("2.0") < 0 (CmdDiann.java:367-373)
    if settings["run_plex"]:
        # channel normalisation is only added for DIA-NN >= 1.9 (CmdDiann.java:381-383)
        args += list(_diann_plex_flags(settings["light"], settings["medium"], settings["heavy"]))
    # diann.cmd-opts is no word of this string: FragPipe adds it as ONE argv element (cmd.add(additionalCmdOpts),
    # CmdDiann.java:386-388), so the entry carries it as its own field (cmd_opts)
    return " ".join(args)


SPECLIB_RT_CAL = {
    # uiComboPqpCal items -> getEasypqpCalOption() (SpeclibPanel.java:173-175, 459-461).
    # "noiRT" is no combo item; loading it leaves the combo at index 0, which is noiRT too.
    "Automatic selection of a run as reference RT": "noiRT",
    "noiRT": "noiRT",
    "Biognosys_iRT": "Biognosys_iRT",
    "ciRT": "ciRT",
    "Pierce_iRT": "Pierce_iRT",
    "User provided RT calibration file": "a tsv file",
    # the saved form (ui.translations:79, config value of the combo item)
    "a tsv file": "a tsv file",
}
SPECLIB_IM_CAL = {
    # uiComboPqpIMCal items -> getEasypqpIMCalOption() (SpeclibPanel.java:187-188, 462-464)
    "Automatic selection of a run as reference IM": "noIM",
    "noIM": "noIM",
    "User provided IM calibration file": "a tsv file",
}
SPECLIB_LABILE_MODES = {
    # uiComboGlycoMode items -> getEasypqpGlycoOption() (SpeclibPanel.java:288, 465-467)
    "Regular (not glyco)": "",
    "O-glyco": "oglyc",
    "N-glyco": "nglyc",
    "N-glyco+HexNAc": "nglyc+",
}


def speclibgen_settings(params: Dict[str, str], all_params: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    """Arguments 8, 10, 11 and 12 of gen_con_spec_lib.py as CmdSpecLibGen.configure builds them
    (CmdSpecLibGen.java:147-187). Missing keys take the SpeclibPanel widgets' initial values
    (SpeclibPanel.java:155-300)."""
    convert_psm = params.get("convert-psm", "false").strip().lower() == "true"
    labile_mode = params.get("easypqp.labile_mode", "Regular (not glyco)")
    if labile_mode not in SPECLIB_LABILE_MODES:
        raise ValueError(f"Unknown speclibgen.easypqp.labile_mode: {labile_mode!r}")
    labile = SPECLIB_LABILE_MODES[labile_mode]
    if labile and not convert_psm:
        # SpeclibPanel.checkGlycoMode (SpeclibPanel.java:520-527), CmdSpecLibGen.java:104-112
        raise ValueError(
            "FragPipe-SpecLib glyco modes are only supported for psm.tsv conversion. "
            "Please change conversion filetype to psm.tsv and try again."
        )

    rt_cal_raw = params.get("easypqp.rt-cal", "Automatic selection of a run as reference RT")
    im_cal_raw = params.get("easypqp.im-cal", "Automatic selection of a run as reference IM")
    if rt_cal_raw not in SPECLIB_RT_CAL:
        raise ValueError(f"Unknown speclibgen.easypqp.rt-cal: {rt_cal_raw!r}")
    if im_cal_raw not in SPECLIB_IM_CAL:
        raise ValueError(f"Unknown speclibgen.easypqp.im-cal: {im_cal_raw!r}")
    rt_cal_file = ""
    im_cal_file = ""
    if SPECLIB_RT_CAL[rt_cal_raw] == "a tsv file":
        rt_cal_file = params.get("easypqp.select-file.text", "").strip()
        if not rt_cal_file:
            raise ValueError("speclibgen.easypqp.rt-cal is a user RT calibration file but speclibgen.easypqp.select-file.text is empty")
    if SPECLIB_IM_CAL[im_cal_raw] == "a tsv file":
        im_cal_file = params.get("easypqp.select-im-file.text", "").strip()
        if not im_cal_file:
            raise ValueError("speclibgen.easypqp.im-cal is a user IM calibration file but speclibgen.easypqp.select-im-file.text is empty")

    # getEasypqp_fragment_types (SpeclibPanel.java:497-512): "['b','y',]", then the quotes are
    # backslash-escaped because the convert extras are one shlex-split argument (CmdSpecLibGen.java:163)
    fragment_types = "[" + "".join(
        f"'{ion}',"
        for ion, default in (("a", "false"), ("b", "true"), ("c", "false"), ("x", "false"), ("y", "true"), ("z", "false"))
        if params.get(f"easypqp.fragment.{ion}", default).strip().lower() == "true"
    ) + "]"
    neutral_loss = params.get("easypqp.neutral_loss", "false").strip().lower() == "true"
    # String.format("--max_delta_unimod %s --max_delta_ppm %s --fragment_types %s %s%s%s%s", ...)
    # (CmdSpecLibGen.java:160-168); %s of a double is Double.toString; OsUtils.asSingleArgument trims.
    convert_args = "--max_delta_unimod {} --max_delta_ppm {} --fragment_types {} {}{}{}{}".format(
        _java_double_to_string(params.get("easypqp.extras.max_delta_unimod", "0.02")),
        _java_double_to_string(params.get("easypqp.extras.max_delta_ppm", "15")),
        fragment_types.replace("'", "\\'"),
        "--enable_unspecific_losses " if neutral_loss else "",
        f"--decoy_prefix {get_decoy_tag(all_params)} " if convert_psm else "",
        f"--labile_mods {labile} " if convert_psm and labile else "",
        "--max_glycan_q {} ".format(_java_double_to_string(params.get("easypqp.extras.max_glycan_qval", "1")))
        if convert_psm
        else "",
    ).strip()
    library_args = "--rt_lowess_fraction {}".format(
        _java_double_to_string(params.get("easypqp.extras.rt_lowess_fraction", "0.01"))
    )
    keep = params.get("keep-intermediate-files", "false").strip().lower() == "true"
    return {
        "rt_cal": "" if rt_cal_file else SPECLIB_RT_CAL[rt_cal_raw],
        "im_cal": "" if im_cal_file else SPECLIB_IM_CAL[im_cal_raw],
        "rt_cal_file": rt_cal_file,
        "im_cal_file": im_cal_file,
        "convert_args": convert_args,
        "library_args": library_args,
        "intermediate": "keep_intermediate_files" if keep else "delete_intermediate_files",
    }


def generate_speclibgen_config(params: Dict[str, str], all_params: Dict[str, Dict[str, str]]) -> str:
    """gen_con_spec_lib.py arguments from the workflow as bash-sourceable assignments.

    Variables (each value shell-quoted; CmdSpecLibGen.java:137-202 builds the argv):
        SPECLIB_RT_CAL, SPECLIB_IM_CAL: argv[8] parts (noiRT/Biognosys_iRT/ciRT/Pierce_iRT, noIM);
            empty when a user calibration file is staged instead (rt_cal_file/im_cal_file in the JSON)
        SPECLIB_CONVERT_ARGS: argv[10], FragPipe-SpecLib convert extras
        SPECLIB_LIBRARY_ARGS: argv[11], FragPipe-SpecLib library extras
        SPECLIB_INTERMEDIATE: argv[12], keep_intermediate_files or delete_intermediate_files
    gen_con_spec_lib.py picks psm.tsv conversion when SPECLIB_CONVERT_ARGS holds
    --decoy_prefix (gen_con_spec_lib.py:125-126), else pepXML conversion.
    """
    settings = speclibgen_settings(params, all_params)
    assignments = (
        ("SPECLIB_RT_CAL", settings["rt_cal"]),
        ("SPECLIB_IM_CAL", settings["im_cal"]),
        ("SPECLIB_CONVERT_ARGS", settings["convert_args"]),
        ("SPECLIB_LIBRARY_ARGS", settings["library_args"]),
        ("SPECLIB_INTERMEDIATE", settings["intermediate"]),
    )
    return "\n".join(f"{name}={shlex.quote(value)}" for name, value in assignments)




def _java_int_to_string(text: str) -> str:
    """String.valueOf(int) of a UiSpinnerInt value read from the workflow."""
    return str(int(text.strip()))


# O-Pair and MBG panel components as their constructors set them: an absent workflow key keeps this value
# in a headless run (SwingUtils.valuesSet sets only the keys in the workflow, utils/SwingUtils.java:859-868).
# Doubles as Double.toString writes them (CmdOPair / CmdMBGMatch concatenate the spinner value).
OPAIR_PANEL_DEFAULTS = {
    "ms2_tol": "20.0",  # tools/opair/OPairPanel.java:118 UiSpinnerDouble.builder(20.0, ...)
    "ms1_tol": "20.0",  # OPairPanel.java:122 builder(20.0, ...)
    "max_glycans": "2",  # OPairPanel.java:127 new UiSpinnerInt(2, 1, 10, 1, 1)
    "min_isotope_error": "0",  # OPairPanel.java:132 new UiSpinnerInt(0, -10, 0, 1, 1)
    "max_isotope_error": "2",  # OPairPanel.java:135 new UiSpinnerInt(2, 0, 10, 1, 1)
    # OPairPanel.java:139, 141 createUiCombo(ActivationTypes.values()): the first item, HCD
    # (tools/enums/ActivationTypes.java:21; UiUtils.java:138-139 maps Enum::name)
    "activation1": "HCD",
    "activation2": "HCD",
    "reverse_scan_order": "false",  # OPairPanel.java:144 createUiCheck("Reverse paired scan order", false)
    "single_scan_type": "false",  # OPairPanel.java:147 createUiCheck("Single scan type", false)
    "filterOxonium": "true",  # OPairPanel.java:150 createUiCheck("Oxonium filter", true)
    "oxonium_minimum_intensity": "0.05",  # OPairPanel.java:153 builder(0.05, ...)
}
MBG_PANEL_DEFAULTS = {
    "max_glycan_q": "0.01",  # tools/mbg/MBGPanel.java:98 UiSpinnerDouble.builder(0.01, ...)
    "fdr": "0.01",  # MBGPanel.java:104 builder(0.01, ...)
    "min_psms": "5",  # MBGPanel.java:109 new UiSpinnerInt(5, 1, 10000, 1)
    "min_glycans": "2",  # MBGPanel.java:113 new UiSpinnerInt(2, 1, 10000, 1)
    "expand_db": "2",  # MBGPanel.java:117 new UiSpinnerInt(2, 1, 10, 1)
    "max_skips": "0",  # MBGPanel.java:121 new UiSpinnerInt(0, 0, 10, 1)
    "allow_chimeric": "false",  # MBGPanel.java:125 new UiCheck("Allow inferring chimeric spectra", null, false)
}


# OPairPanel.java:139, 141 createUiCombo(ActivationTypes.values()) items (ActivationTypes.java:21)
OPAIR_ACTIVATION_TYPES = ("HCD", "ETD", "CID", "ECD")


def _opair_value(params: Dict[str, str], key: str) -> str:
    """An O-Pair panel value: the workflow's, else the panel constructor's (OPAIR_PANEL_DEFAULTS).
    An activation combo given a value it lacks selects its first item (UiCombo.fromString, UiCombo.java:35-41)."""
    value = params.get(key, OPAIR_PANEL_DEFAULTS[key])
    if key in ("activation1", "activation2") and value not in OPAIR_ACTIVATION_TYPES:
        return OPAIR_ACTIVATION_TYPES[0]
    return value


def _mbg_value(params: Dict[str, str], key: str) -> str:
    """An MBG panel value: the workflow's, else the panel constructor's (MBG_PANEL_DEFAULTS)."""
    return params.get(key, MBG_PANEL_DEFAULTS[key])


def _opair_runs_shepherd(all_params: Dict[str, Dict[str, str]]) -> bool:
    """FragpipeRun.java:2095 isRunShepherd = ptmshepherdPanel.isRun() || ptmsGlycanPanel.isRun()."""
    ptms = all_params.get("ptmshepherd", {})
    return (
        ptms.get("run-shepherd", "false").lower() == "true"
        or ptms.get("run_glyco_mode", "false").lower() == "true"
    )


def check_glyco_hard_stops(all_params: Dict[str, Dict[str, str]]) -> None:
    """FragPipe refuses these glyco workflows before running anything (exit non-zero, like FragPipe).

    - FragpipeRun.java:2099-2111: O-Pair and PTM-Shepherd (incl. glycan assignment) together.
    - CmdOPair.java:89-96: single scan type with HCD, or activation1 CID. The Java condition is
      `(single && a1 == HCD) || a1 == CID` (operator precedence), kept as is.
    - CmdMBGMatch.java:145-158: MBG without residues to add (the residue parse check,
      verifyResidues, needs the tools dir and runs in the MBG module).
    """
    if get_run_flag(all_params, "opair"):
        if _opair_runs_shepherd(all_params):
            sys.exit(
                "ERROR: O-Pair and PTM-Shepherd (including glycan composition assignment) cannot be "
                "run in the same workflow. Please disable one or the other to continue."
            )
        opair = all_params.get("opair", {})
        activation1 = _opair_value(opair, "activation1")
        single = _opair_value(opair, "single_scan_type").lower() == "true"
        if (single and activation1 == "HCD") or activation1 == "CID":
            sys.exit(
                "ERROR: Single activation type must be hybrid, not HCD/CID, for O-Pair because "
                "c/z ions are required."
            )
    if get_run_flag(all_params, "mbg"):
        if not all_params.get("mbg", {}).get("residues_to_add", "").strip():
            sys.exit(
                "ERROR: no glycan residues or mods specified for MBG. Please set residue(s) to "
                "match and try again."
            )


# Runtime values the modules fill into the argv lists below (the parser does not know them).
TOOLS_DIR_TOKEN = "@TOOLS_DIR@"  # task.ext.fragpipe_tools_dir (FragpipeLocations.getDirTools())
THREADS_TOKEN = "@THREADS@"  # task.cpus (FragPipe --threads)
PSM_TOKEN = "@PSM@"  # <group dir>/psm.tsv (group.outputDir(workdir).resolve("psm.tsv"))
OUTDIR_TOKEN = "@OUTDIR@"  # <group dir>
FILELIST_TOKEN = "@FILELIST@"  # <group dir>/filelist_opair.txt
OXONIUM_FILE_TOKEN = "@OXONIUM_FILE@"  # staged copy of the user's opair.oxonium_filtering_file
MANIFEST_TOKEN = "@MANIFEST@"  # fragpipe-files.fp-manifest
NOPASEF_TOKEN = "@NOPASEF@"  # TabWorkflow.getInputDataType() != ImMsTimsTof (from the LC-MS files)


def opair_oxonium_file(params: Dict[str, str]) -> str:
    """The user's oxonium ion definitions file ('' when unset or filterOxonium=false).

    CmdOPair.java:171-181: `if (params.isFilterOxonium()) { if (params.getOxoRulesFilePath().isEmpty()) {
    ... cmd.add("-f " + oxoPath.toAbsolutePath()); } else { cmd.add("-f " + params.getOxoRulesFilePath()); }`
    """
    if _opair_value(params, "filterOxonium").lower() != "true":
        return ""
    return params.get("oxonium_filtering_file", "").strip()


def generate_opair_args(params: Dict[str, str]) -> list:
    """O-Pair argv after `dotnet <tools>/opair/CMD.dll`, token for token as CmdOPair.java:177-209.

    FragPipe adds every flag and its value as ONE ProcessBuilder token ("-b 20.0"), and adds
    "-z " even when the allowed sites are empty. Doubles are Double.toString of the spinner value.
    activation1/2, reverse_scan_order and single_scan_type are PairScans inputs, not O-Pair flags.
    """
    a = [
        "-b " + _java_double_to_string(_opair_value(params, "ms2_tol")),
        "-c " + _java_double_to_string(_opair_value(params, "ms1_tol")),
    ]
    if _opair_value(params, "filterOxonium").lower() == "true":
        if opair_oxonium_file(params):
            # CmdOPair.java:178-180: `cmd.add("-f " + params.getOxoRulesFilePath())`; the module stages
            # the file (entry["oxonium_file"]) and fills its task path here
            a.append(f"-f {OXONIUM_FILE_TOKEN}")
        else:
            # CmdOPair.DEFAULT_OXO_FILE in tools/opair/Glycan_Mods (checked to exist by the module)
            a.append(f"-f {TOOLS_DIR_TOKEN}/opair/Glycan_Mods/oxonium_filter.txt")
        a.append(
            "-m " + _java_double_to_string(_opair_value(params, "oxonium_minimum_intensity"))
        )
    glyco_db = params.get("glyco_db", "").strip()
    if glyco_db:
        a.append("-g " + glyco_db)
    a += [
        f"-x {TOOLS_DIR_TOKEN}/Glycan_Databases/glycan_residues.txt",
        f"-y {TOOLS_DIR_TOKEN}/Glycan_Databases/glycan_mods.txt",
        "-n " + _java_int_to_string(_opair_value(params, "max_glycans")),
        "-z " + params.get("allowed_sites", "").strip(),
        f"-t {THREADS_TOKEN}",
        "-i " + _java_int_to_string(_opair_value(params, "min_isotope_error")),
        "-j " + _java_int_to_string(_opair_value(params, "max_isotope_error")),
        f"-r {FILELIST_TOKEN}",
        f"-s {PSM_TOKEN}",
        f"-o {OUTDIR_TOKEN}",
    ]
    return a


def generate_pairscans_args(params: Dict[str, str]) -> list:
    """PairScans arguments after `<spectra file> <threads>` (CmdPairScans.java:121-126)."""
    return [
        _opair_value(params, "activation1"),
        _opair_value(params, "activation2"),
        str(_opair_value(params, "reverse_scan_order").lower() == "true").lower(),
        str(_opair_value(params, "single_scan_type").lower() == "true").lower(),
    ]


def opair_spectra_kind(all_params: Dict[str, Dict[str, str]]) -> str:
    """Which spectra O-Pair's filelist names (CmdOPair.java:107-114, FragpipeRun.java:2090).

    'calibrated' when tabMsf.isWriteCalMzml() && getMassCalibration() > 0 (every file:
    <base>_calibrated.mzML); otherwise 'uncalibrated' for .raw/.d and the input file for the rest,
    decided per file by the workflow from the file extension.
    """
    return "calibrated" if msfragger_writes_calibrated_mzml(all_params) else "uncalibrated_or_input"


def generate_mbg_args(params: Dict[str, str], all_params: Dict[str, Dict[str, str]]) -> list:
    """MBG argv after `-cp <MBG jar>:<batmass-io>:<IonQuant jar> com.mbg.MBG`, token for token as
    CmdMBGMatch.java:137-208 (each flag and value a separate token).

    The tolerances come from the IonQuant panel (uiCompsRepresentation = quantPanelLabelfree.toMap());
    --runtmt is tmtiPanel.isRun().
    """
    residues = params.get("residues_to_add", "").strip()
    ionquant = all_params.get("ionquant", {})
    run_tmt = all_params.get("tmtintegrator", {}).get("run-tmtintegrator", "false").lower() == "true"
    return [
        "--match",
        "--psm", PSM_TOKEN,
        "--manifest", MANIFEST_TOKEN,
        "--toaddresiduals", residues,
        "--residuedb", f"{TOOLS_DIR_TOKEN}/Glycan_Databases/glycan_residues.txt",
        "--glycanmoddb", f"{TOOLS_DIR_TOKEN}/Glycan_Databases/glycan_mods.txt",
        "--maxq", _java_double_to_string(_mbg_value(params, "max_glycan_q")),
        "--minpsms", _java_int_to_string(_mbg_value(params, "min_psms")),
        "--minglycans", _java_int_to_string(_mbg_value(params, "min_glycans")),
        "--fdr", _java_double_to_string(_mbg_value(params, "fdr")),
        # quantPanelLabelfree.toMap() widget strings (CmdMBGMatch.java:189-195): the panel default
        # when the key is absent (QuantPanelLabelfree.java:280-285)
        "--mztol", ionquant_panel_value("mztol", ionquant.get("mztol")),
        "--rttol", ionquant_panel_value("rttol", ionquant.get("rttol")),
        "--imtol", ionquant_panel_value("imtol", ionquant.get("imtol")),
        "--nopasef", NOPASEF_TOKEN,
        "--numthreads", THREADS_TOKEN,
        "--runtmt", str(run_tmt).lower(),
        "--expanddb", _java_int_to_string(_mbg_value(params, "expand_db")),
        "--maxskips", _java_int_to_string(_mbg_value(params, "max_skips")),
        "--allowchimeric", str(_mbg_value(params, "allow_chimeric").lower() == "true").lower(),
    ]


def generate_saintexpress_config(params: Dict[str, str]) -> str:
    """SAINTexpress -R/-L flags (CmdSaintExpress.java:85-87); cmd-opts is emitted separately.

    -R is SaintexpressPanel max-replicates (default 10), -L virtual-controls (default 100)
    (SaintexpressPanel.java:111-115). CmdSaintExpress adds cmd-opts as ONE argv element
    (CmdSaintExpress.java:88-90), so it is the JSON field cmd_opts, not part of args.
    """
    max_replicates = int(params.get("max-replicates", "10"))
    virtual_controls = int(params.get("virtual-controls", "100"))
    return f"-R{max_replicates} -L{virtual_controls}"


# Crystal-C algorithm settings FragPipe uses: CrystalcPanel.toParams() loads the jar's
# crystalc.params (CrystalcParams.loadDefault) and never reads the GUI or workflow values
# (CrystalcPanel.java:139-148). Values from FP/tools/crystalc/crystalc.params:9-13.
CRYSTALC_DEFAULTS = (
    ("precursor_charge", "1 6"),
    ("isotope_number", "3"),
    ("precursor_mass", "20"),
    ("precursor_isolation_window", "0.7"),
    ("correct_isotope_error", "false"),
)


def generate_crystalc_config(params: Dict[str, str]) -> str:
    """Crystal-C algorithm params as FragPipe writes them (CrystalcPanel.toParams).

    Workflow ``crystalc.*`` values other than the run flag are ignored, as in FragPipe.
    The runtime keys (thread, fasta, raw_file_location, raw_file_extension,
    output_location) are written by the CRYSTALC module (CmdCrystalc.java:190-196).
    """
    return "\n".join(f"{key} = {value}" for key, value in CRYSTALC_DEFAULTS)


def generate_tabrun_config(params: Dict[str, str]) -> str:
    """Generate tab-run (GUI Run tab) settings as key=value pairs.

    tab-run is NOT a standalone tool. It controls FragPipe GUI behavior:
      write_sub_mzml: write subset mzML files
      delete_temp_files: delete temp files after pipeline; not replicated, the work directory
        keeps them (FragpipeRun.java:547-578: "It won't affect the result")
      export_matched_fragments: export matched fragment ions
      sub_mzml_prob_threshold: probability threshold for subset mzML

    These are stored for reference; the run decisions are the writesubmzml and
    exportmatchedfragments entries (WRITE_SUB_MZML, EXPORT_MATCHED_FRAGMENTS).

    Reference: TabRun.java, CmdWriteSubMzml.java, CmdExportMatchedFragments.java
    """
    lines = []
    for key, value in sorted(params.items()):
        if not value or value.lower() in ("null", "none"):
            continue
        lines.append(f"{key} = {value}")
    return "\n".join(lines)


def generate_quantitation_config(params: Dict[str, str]) -> str:
    """Generate quantitation (GUI quant panel) settings as key=value pairs.

    quantitation is NOT a standalone tool. It holds the GUI panel state:
      run-label-free-quant: whether label-free quantification is enabled

    This flag is read by FragPipe to decide whether to run IonQuant/FreeQuant.
    We include it in the JSON output for completeness.

    Reference: TabPtms.java, QuantPanelLabelfree.java
    """
    lines = []
    for key, value in sorted(params.items()):
        if not value or value.lower() in ("null", "none"):
            continue
        lines.append(f"{key} = {value}")
    return "\n".join(lines)


def _java_boolean(value: str) -> str:
    """String.valueOf(boolean) of a workflow checkbox value."""
    return "true" if value.strip().lower() == "true" else "false"


# TMT-Integrator prot_norm -> normalisation name in the FPOP input file names
# (TmtiConfProps.COMBO_NORM; CmdFpopQuant.getTmtiMultisiteTsv/getTmtiPeptideTsv).
FPOP_TMT_NORM_SUFFIX = {"0": "None", "-1": "None", "1": "MD", "2": "GN"}


def generate_fpop_config(params: Dict[str, str], all_params: Dict[str, Dict[str, str]]) -> List[str]:
    """FragPipe_FPOP_Analysis.py arguments of FPOP Quant (CmdFpopQuant.java:90-99), one word each.

    Widgets are fpop.fragpipe.fpop.* (FpopQuantPanel.java:42, 111-125); a key missing from the
    workflow keeps the widget default (region size 1, empty labels, unchecked boxes). Order: input
    table, region size, control label, FPOP label, subtract control, TMT flag[, unmodified peptide
    table]; empty labels stay empty words. Label-free mode reads IonQuant's
    combined_modified_peptide.tsv; TMT mode reads tmt-report/abundance_multi-site_<NORM>.tsv and
    tmt-report-unmod/abundance_peptide_<NORM>.tsv (CmdFpopQuant.java:56-63, 114-145), the names the
    FPOP module stages its inputs under.
    """
    is_tmt = _java_boolean(params.get("fragpipe.fpop.fpop-tmt", "false"))
    if is_tmt == "true":
        tmti = all_params.get("tmtintegrator", {})
        if tmti.get("run-tmtintegrator", "false").lower() != "true":
            # FragpipeRun.java:2181-2193
            raise ValueError(
                "FPOP TMT Quant was requested on Downstream tab, but TMT-Integrator is not run. "
                "Please enable TMT-Integrator"
            )
        # tmtiPanel.getNormMethod(): the Normalization combo, first item None when absent
        # (TmtiConfProps.java:78-82, TmtiPanel.java:356, :413-415)
        tmt_norm = FPOP_TMT_NORM_SUFFIX[_tmti_widget_value("prot_norm", tmti.get("prot_norm"))]
        first = [f"tmt-report/abundance_multi-site_{tmt_norm}.tsv"]
        last = [f"tmt-report-unmod/abundance_peptide_{tmt_norm}.tsv"]
    else:
        first = ["combined_modified_peptide.tsv"]
        last = []
    return first + [
        str(int(params.get("fragpipe.fpop.region_size", "1"))),
        params.get("fragpipe.fpop.label_control", ""),
        params.get("fragpipe.fpop.label_fpop", ""),
        _java_boolean(params.get("fragpipe.fpop.subtract-control", "false")),
        is_tmt,
    ] + last


def generate_fpop_coadaptr_config(params: Dict[str, str]) -> str:
    """FPOP mod masses argument of CmdFPOPcoadaptr (FpopCoadaptrPanel.getFpopMods(), java:83-100).

    The text is split on [\\s,;/]+; each number is formatted %.4f and the list comma-joined.
    Tokens that are not numbers are logged and skipped, as FragPipe does headless.
    """
    output = []
    for token in re.split(r"[\s,;/]+", params.get("coadaptr.fpop.fpop_masses", "")):
        try:
            value = Decimal(repr(float(token)))
        except ValueError:
            print(f"Could not read FPOP mod mass {token} as a number.", file=sys.stderr)
            continue
        # java.util.Formatter rounds the decimal value HALF_UP
        output.append(str(value.quantize(Decimal("0.0001"), rounding=ROUND_HALF_UP)))
    return ",".join(output)


def _append_prepend_once(base: str, quote: str) -> str:
    """org.nesvilab.utils.StringUtils.appendPrependOnce."""
    if not base.endswith(quote):
        base = base + quote
    if not base.startswith(quote):
        base = quote + base
    return base


def generate_metaproteomics_config(params: Dict[str, str]) -> Dict[str, str]:
    """FP-Meta DbOptimizer values (CmdMetaproteomics.java:70-111), keyed by their option name.

    TabDatabase widgets metaproteomics.* (TabDatabase.java:382-505); a missing key keeps the
    widget default. Floats are String.valueOf(float), the host name is wrapped in double quotes
    once (StringUtils.appendPrependOnce), and cmd-line-opts is one trimmed argv element ('' = none).
    metaproteomicsArgv (fragpipe_utils_a.nf) puts them in FragPipe's argv order.
    """
    return {
        "qvalue": _java_float_to_string(params.get("qvalue", "0.01")),
        "deltaHyperscore": _java_float_to_string(params.get("delta-hyperscore", "0.0")),
        "minPeptCntPerProt": str(int(params.get("min-pept-cnt-per-prot", "1"))),
        "minUniqPeptCntPerProt": str(int(params.get("min-uniq-pept-cnt-per-prot", "1"))),
        "minUniqPeptCnt": str(int(params.get("min-uniq-pept-cnt", "3"))),
        "hostName": _append_prepend_once(params.get("host-name", "Homo sapiens"), '"'),
        "iterations": str(int(params.get("iterations", "3"))),
        "cmdLineOpts": params.get("cmd-line-opts", "").strip(),
    }


def generate_generic_config(
    params: Dict[str, str], prefix: str, file_mode: bool = False
) -> str:
    """Generate generic config for tools not explicitly configured.

    This is the fallback handler for any tool prefix found in the workflow file
    that does not have a dedicated generator. Ensures no parameters are silently
    dropped for tools added in future FragPipe versions.

    Args:
        params: Tool parameters from the workflow file.
        prefix: Tool prefix name (for documentation).
        file_mode: If True, output key=value properties format (for file output).
                   If False, output --key value CLI format (for JSON args).
    """
    parts = []

    for key, value in sorted(params.items()):
        if value and value.lower() not in ("", "null", "none"):
            # Skip run-* flags (control flags, not CLI params)
            if key.startswith("run-"):
                continue
            if file_mode:
                parts.append(f"{key} = {value}")
            else:
                parts.append(f"--{key} {value}")

    return "\n".join(parts) if file_mode else " ".join(parts)


def _is_true(
    all_params: Dict[str, Dict[str, str]], prefix: str, key: str, default: str
) -> bool:
    """Read a .workflow checkbox ``prefix.key`` (``_global`` for top-level keys) as a boolean."""
    return all_params.get(prefix, {}).get(key, default).lower() == "true"


# FragPipe 24.0 initial state of a run checkbox / radio button, used when the .workflow omits
# its key (Fragpipe.loadUi only sets components whose key is present). Tools not listed start
# unchecked (e.g. MSBoosterPanel.java:282, CrystalcPanel.java:190, TmtiPanel.java:256).
RUN_FLAG_DEFAULTS = {
    "msfragger": "true",  # TabMsfragger.java:570 new UiCheck("Run MSFragger", null, true)
    "proteinprophet": "true",  # ProtProphPanel.java:85 createUiCheck("Run ProteinProphet", true)
    "filter": "true",  # ReportPanel.java:180 new UiCheck("Generate reports", null, true)
}

# Radio-button pairs sharing one ButtonGroup: (default-selected tool, other tool).
# ButtonGroup.add() keeps the first selected button and deselects later ones, so the first
# panel built wins: PercolatorPanel before PepProphPanel (PSMValidation.java:81-82, both
# UiRadio(..., true)), FreeQuant before IonQuant (QuantPanelLabelfree.java:222-223,
# UiRadio("FreeQuant (deprecated)", null, true) at :234, UiRadio("IonQuant", null, false) at :259).
# A workflow value reaches the button as setSelected(Boolean.valueOf(s)) (UiRadio.java:40-42), and
# in a ButtonGroup setSelected(false) does not clear the selected button: only `true` on the other
# button moves the selection, so `false` and an absent key are the same.
RUN_FLAG_RADIO_PAIRS = (("percolator", "peptideprophet"), ("freequant", "ionquant"))

# The panels of FragPipe's Validation tab (TabValidation.java:61-66); "filter" is the Report panel,
# gated separately because FreeQuant also runs Philosopher report (FragpipeRun.java:1767).
VALIDATION_TAB_TOOLS = ("crystalc", "msbooster", "percolator", "peptideprophet", "ptmprophet", "proteinprophet")


def _validation_tab_runs(all_params: Dict[str, Dict[str, str]]) -> bool:
    """TabValidation "Run validation tools" (TabValidation.java:58-59, default true)."""
    return _is_true(all_params, "_global", "run-validation-tab", "true")


def _run_checkbox(all_params: Dict[str, Dict[str, str]], tool_name: str) -> bool:
    """
    State of a tool's own run checkbox/radio, with FragPipe's initial state for an absent key.

    @param all_params Parsed workflow (parse_workflow_file)
    @param tool_name Key of TOOL_RUN_FLAGS
    @return True when the checkbox/radio is selected
    """
    def is_true(tool: str) -> bool:
        prefix, run_key = TOOL_RUN_FLAGS[tool]
        return all_params.get(prefix, {}).get(run_key, "").lower() == "true"

    for first, second in RUN_FLAG_RADIO_PAIRS:
        if tool_name in (first, second):
            # parse_workflow_file resolved both-true twins (_resolve_radio_twins)
            return is_true(second) if tool_name == second else not is_true(second)
    prefix, run_key = TOOL_RUN_FLAGS[tool_name]
    value = all_params.get(prefix, {}).get(run_key)
    if value is not None:
        return value.lower() == "true"
    return RUN_FLAG_DEFAULTS.get(tool_name, "false") == "true"


# tmtintegrator.extraction_tool is TmtiPanel's "Intensity extraction tool" combo: the workflow
# stores the item text (TmtiPanel.java:259-264) and FragPipe gates on the item index
# (getIntensityExtractionTool(), TmtiPanel.java:836-838; FragpipeRun.java:1925-1936).
EXTRACTION_TOOL_INDEX = {
    "IonQuant": 0,
    "Philosopher": 1,
    "Skip extraction. Run TMT-Integrator only": 2,
}


def get_extraction_tool_index(all_params: Dict[str, Dict[str, str]]) -> int:
    """Index of tmtintegrator.extraction_tool (0 IonQuant, 1 Philosopher, 2 skip extraction).

    A workflow without the key keeps the combo's first item (IonQuant)."""
    value = all_params.get("tmtintegrator", {}).get("extraction_tool")
    if value is None:
        return 0
    if value not in EXTRACTION_TOOL_INDEX:
        raise ValueError(
            f"Unknown tmtintegrator.extraction_tool={value}; FragPipe offers: "
            + ", ".join(EXTRACTION_TOOL_INDEX)
        )
    return EXTRACTION_TOOL_INDEX[value]


def get_run_flag(all_params: Dict[str, Dict[str, str]], tool_name: str) -> bool:
    """
    Extract the run flag for a specific tool from the workflow params.

    Returns:
        True if the tool should run, False otherwise.
    """
    if tool_name not in TOOL_RUN_FLAGS:
        return False

    prefix, run_key = TOOL_RUN_FLAGS[tool_name]

    # "Run validation tools" (TabValidation.java:58-59, :94 setEnablementUpdater(this, pContent, checkRun))
    # disables the panels of the Validation tab, and each one's isRun() is isEnabledAndChecked
    # (CrystalcPanel.java:150, MSBoosterPanel.java:291, PSMValidation.java:57, PtmProphetPanel.java:237,
    # ProtProphPanel.java:48, ReportPanel.java:194). A FragPipe 24.0 --dry-run of LFQ-MBR with
    # run-validation-tab=false plans none of their steps.
    if tool_name in VALIDATION_TAB_TOOLS and not _validation_tab_runs(all_params):
        return False

    # Labelquant: runs when TMTIntegrator is enabled AND intensity_extraction_tool==1
    # (FragpipeRun.java:1935: tmtiPanel.isRun() && tmtiPanel.getIntensityExtractionTool() == 1)
    if tool_name == "labelquant":
        tmti_params = all_params.get("tmtintegrator", {})
        tmti_runs = tmti_params.get("run-tmtintegrator", "false").lower() == "true"
        return tmti_runs and get_extraction_tool_index(all_params) == 1

    # The Run tab is no tool: its sub-mzML writer and matched-fragment export run on their own
    # checkboxes (FragpipeRun.java:2388-2391, :2411).
    if tool_name == "tabrun":
        tabrun_params = all_params.get("tab-run", {})
        return any(
            tabrun_params.get(key, "false").strip().lower() == "true"
            for key in ("write_sub_mzml", "export_matched_fragments")
        )

    # Database always runs (no run flag)
    if run_key is None:
        return prefix in all_params

    # Label-free quant (FragPipe 24.0 FP/tools/ionquant/QuantPanelLabelfree.java:110-116):
    #   isRunIonQuant()  = isRun() && uiRadioUseIonquant
    #   isRunFreeQuant() = isRun() && uiRadioUseFreequant
    # isRun() is the "Run MS1 quant" checkbox quantitation.run-label-free-quant (:187-188);
    # ionquant.run-ionquant / freequant.run-freequant are only the radio buttons.
    if tool_name in ("ionquant", "freequant"):
        return _is_true(all_params, "quantitation", "run-label-free-quant", "false") and (
            _run_checkbox(all_params, tool_name)
        )

    # PSM validation master switch (PercolatorPanel.java:64-67, PepProphPanel.java:63-66):
    #   isRun() = psmValidation.isRun() && checkRun. The top-level key run-psm-validation
    #   (PSMValidation.java:76-77) defaults to true.
    if tool_name in ("percolator", "peptideprophet"):
        return _is_true(all_params, "_global", "run-psm-validation", "true") and (
            _run_checkbox(all_params, tool_name)
        )

    # PTM-Shepherd runs for profiling or glycan assignment (FragpipeRun.java:2096):
    #   isRunShepherd = ptmshepherdPanel.isRun() || ptmsGlycanPanel.isRun()
    # The glycan-assignment checkbox is ptmshepherd.run_glyco_mode (PTMSGlycanAssignPanel.java:60,117-119,151).
    if tool_name == "ptmshepherd":
        return _run_checkbox(all_params, tool_name) or _is_true(
            all_params, prefix, "run_glyco_mode", "false"
        )

    # Philosopher report runs for the report panel or label-free FreeQuant (FragpipeRun.java:1767):
    #   new CmdPhilosopherReport(reportPanel.isRun() || quantPanelLabelfree.isRunFreeQuant(), wd)
    if tool_name == "filter":
        return (_validation_tab_runs(all_params) and _run_checkbox(all_params, tool_name)) or get_run_flag(
            all_params, "freequant"
        )

    # Handle nested keys (e.g., "fragpipe.fpop.run-fpop")
    if "." in run_key:
        parts = run_key.split(".")
        # For fpop, the prefix is "fpop" but the key structure is nested
        if prefix in all_params:
            # Try full nested key
            full_key = run_key
            if full_key in all_params.get(prefix, {}):
                return all_params[prefix][full_key].lower() == "true"
            # Try last part of key
            if parts[-1] in all_params.get(prefix, {}):
                return all_params[prefix][parts[-1]].lower() == "true"
        return False

    # Standard case: prefix.run_key
    return _run_checkbox(all_params, tool_name)


def generate_tool_args(all_params: Dict[str, Dict[str, str]], tool_name: str) -> str:
    """
    Generate the CLI arguments for a tool (used in tool_configs.json).

    Returns clean CLI argument strings that can be passed directly to tools.
    For philosopher subcommand tools, extracts just the CLI flags (not config file format).

    Returns:
        String containing CLI arguments.
    """
    if tool_name not in TOOL_RUN_FLAGS:
        return ""

    prefix, _ = TOOL_RUN_FLAGS[tool_name]

    # Labelquant params are derived from tmtintegrator, not a separate prefix
    if tool_name == "labelquant":
        params = all_params.get("labelquant", {})
        return generate_labelquant_config(params, all_params)

    # FreeQuant options come from the panel widgets, which exist without freequant.* keys
    if tool_name == "freequant":
        return generate_freequant_config(all_params.get("freequant", {}))

    if prefix not in all_params:
        return ""

    params = all_params[prefix]

    # Map tool names to their config generators
    if tool_name == "ionquant":
        return generate_ionquant_config(params, all_params)
    elif tool_name == "percolator":
        # Don't include the "disabled" comment in args
        args = generate_percolator_config(params, get_run_flag(all_params, "percolator"))
        return args if not args.startswith("#") else ""
    elif tool_name == "msbooster":
        return generate_msbooster_config(params, all_params, get_run_flag(all_params, "msbooster"))
    elif tool_name == "msfragger":
        return generate_msfragger_config(params, all_params)
    elif tool_name == "filter":
        # Return just the filter CLI flags (not config file format)
        return _get_philosopher_cli_args(
            all_params.get("phi-report", {}), "phi-report", all_params
        )
    elif tool_name == "peptideprophet":
        # Return just the CLI flags (not config file format)
        return _get_philosopher_cli_args(params, "peptide-prophet", all_params)
    elif tool_name == "proteinprophet":
        # Return just the CLI flags (not config file format)
        return _get_philosopher_cli_args(params, "protein-prophet", all_params)
    elif tool_name == "database":
        # Return just the CLI flags (not config file format)
        return _get_philosopher_cli_args(
            all_params.get("database", {}), "database", all_params
        )
    elif tool_name == "tmtintegrator":
        return generate_tmtintegrator_config(params)
    elif tool_name == "ptmprophet":
        # M14: PTMProphet auto-generates mods from MSFragger when override-defaults=false.
        # Replicates PtmProphetPanel.getCmdLineOpts() (PtmProphetPanel.java:241-247).
        return generate_ptmprophet_config(params, all_params)
    elif tool_name == "ptmshepherd":
        return generate_ptmshepherd_config(params, all_params)
    elif tool_name == "labelquant":
        # Labelquant params are derived from tmtintegrator settings, not a separate prefix
        return generate_labelquant_config(params, all_params)
    # New tool-specific generators (no silent parameter drops)
    elif tool_name == "diaumpire":
        return generate_diaumpire_config(params)
    elif tool_name == "diatracer":
        return generate_diatracer_config(params)
    elif tool_name == "diann":
        return generate_diann_config(params)
    elif tool_name == "speclibgen":
        return generate_speclibgen_config(params, all_params)
    elif tool_name == "opair":
        return generate_opair_args(params)
    elif tool_name == "mbg":
        return generate_mbg_args(params, all_params)
    elif tool_name == "saintexpress":
        return generate_saintexpress_config(params)
    elif tool_name == "crystalc":
        return generate_crystalc_config(params)
    elif tool_name == "tabrun":
        return generate_tabrun_config(params)
    elif tool_name == "exportmatchedfragments":
        # CmdExportMatchedFragments.configure() takes no options from the workflow:
        # FragmentsExportMainClass <workDir> <nThreads> r (CmdExportMatchedFragments.java:58-66)
        return ""
    elif tool_name == "writesubmzml":
        # CmdWriteSubMzml.configure() takes no options from the workflow; its argv is built by
        # writeSubMzmlArgv (fragpipe_utils_a.nf) from the probability threshold (CmdWriteSubMzml.java:97-130)
        return ""
    elif tool_name == "quantitation":
        return generate_quantitation_config(params)
    elif tool_name == "fpop":
        return generate_fpop_config(params, all_params)
    elif tool_name == "fpop_coadaptr":
        return generate_fpop_coadaptr_config(params)
    elif tool_name == "metaproteomics":
        return generate_metaproteomics_config(params)
    else:
        # Generic passthrough for any remaining tool prefix.
        # Ensures no parameters are silently dropped for tools added
        # in future FragPipe versions.
        return generate_generic_config(params, prefix)


def _get_philosopher_cli_args(
    params: Dict[str, str],
    subcommand_prefix: str,
    all_params: Dict[str, Dict[str, str]] = None,
) -> str:
    """
    Extract clean CLI argument strings for philosopher subcommand tools.

    Unlike generate_philosopher_config() which produces config file format
    (e.g., "proteinprophet=--maxppmdiff 2000000"), this returns just the
    CLI flags (e.g., "--maxppmdiff 2000000") for use in tool_configs.json.
    """
    if subcommand_prefix == "phi-report":
        # Final filter options depend on the IonQuant / ProteinProphet decisions
        # (fragpipe_filter_options, CmdPhilosopherFilter.java:82-97); --razor is appended by the
        # module after --protxml. The legacy phi-report.config carries the same options.
        if all_params is None:
            return params.get("filter", "")
        return philosopher_filter_entry(all_params)["args"]

    elif subcommand_prefix == "peptide-prophet":
        return peptideprophet_cmd_opts(all_params)

    elif subcommand_prefix == "protein-prophet":
        return params.get("cmd-opts", "")

    elif subcommand_prefix == "database":
        return f"--prefix {get_decoy_tag(all_params)}"

    return ""


def generate_json_output(
    all_params: Dict[str, Dict[str, str]],
    data_types: List[str],
) -> Dict[str, Dict[str, Any]]:
    """
    Generate JSON output with run flags and args for all tools.

    Returns:
        Dict mapping tool_name -> {"run": bool, "args": str}
    """
    result = {}

    for tool_name in TOOL_RUN_FLAGS:
        run_flag = get_run_flag(all_params, tool_name)
        # The TMT FreeQuant that precedes LabelQuant takes the label-free panel's FreeQuant
        # options too (FragpipeRun.java:2049-2055), so they are needed when either runs.
        needs_args = run_flag or (
            tool_name == "freequant" and get_run_flag(all_params, "labelquant")
        )
        args = generate_tool_args(all_params, tool_name) if needs_args else ""

        # Indicate config type: "params_file" for tools that use parameter files,
        # "cli" for tools that use CLI arguments,
        # "bash_source" for tools that output bash-sourceable KEY='value' lines,
        # "gui_only" for GUI-only tab settings that are not standalone tools
        params_file_tools = (
            "msfragger",
            "msbooster",
            "tmtintegrator",
            "ptmshepherd",
            "diaumpire",
            "crystalc",
        )
        # O-Pair, MBG and FPOP args are argv token lists, Metaproteomics args a value map, passed to the tool without a shell
        argv_tools = ("opair", "mbg", "fpop", "metaproteomics")
        bash_source_tools = ("speclibgen",)
        gui_only_tools = ("tabrun", "quantitation")
        if tool_name in argv_tools:
            config_type = "argv"
        elif tool_name in params_file_tools:
            config_type = "params_file"
        elif tool_name in bash_source_tools:
            config_type = "bash_source"
        elif tool_name in gui_only_tools:
            config_type = "gui_only"
        else:
            config_type = "cli"

        entry = {
            "run": run_flag,
            "args": args,
            "config_type": config_type,
        }

        # Database slices: TabMsfragger.getNumDbSlices() (msfragger.misc.slice-db spinner,
        # TabMsfragger.java:868-869, :1964-1966), passed to CmdMsfragger.configure (FragpipeRun.java:1470).
        if tool_name == "msfragger":
            entry["num_slices"] = msfragger_num_slices(all_params.get("msfragger", {}), all_params)
            entry.update(msfragger_search_settings(all_params, data_types))

        # modmasses_ionquant.txt lines (None: no --modlist), FragpipeRun.java:1852-1857
        if tool_name == "ionquant" and run_flag:
            entry["modmasses"] = _parse_mod_masses(all_params)

        # labelquant --brand (CmdLabelquant.java:144-145), placed after --annot by the module
        if tool_name == "labelquant" and run_flag:
            entry["brand"] = labelquant_brand(all_params)

        # Calibrated mzML written by MSFragger: tabMsf.isWriteCalMzml() && tabMsf.getMassCalibration() > 0
        # (FragpipeRun.java:1613, :1642) decides the spectra path in pepXMLs (CmdPercolator.java:219-224)
        if tool_name == "msfragger":
            entry["has_calibrated_mzml"] = msfragger_writes_calibrated_mzml(all_params)

        # Philosopher filter options and report flags (FragpipeRun.java:1702-1780)
        if tool_name == "filter" and run_flag:
            entry.update(philosopher_filter_entry(all_params))

        # ptmshepherd.annotation_file: the custom annotation file PTMSHEPHERD stages ("" for the
        # keywords unimod/common/glyco); the config line refers to the copy through @ANNOTATION_FILE@
        if tool_name == "ptmshepherd" and run_flag:
            entry["annotation_file"] = ptms_annotation_file(all_params.get("ptmshepherd", {}), all_params)

        # msbooster.spectral-library-path (CmdMSBooster.java:180-184): the library file MSBOOSTER stages, "" without one
        if tool_name == "msbooster" and run_flag:
            entry["spectral_library"] = all_params.get("msbooster", {}).get("spectral-library-path", "").strip()

        # PercolatorOutputToPepXML min-prob and intermediate TSV handling (CmdPercolator.java:233-238, :298)
        if tool_name == "percolator" and run_flag:
            entry["min_prob"] = percolator_min_prob(all_params)
            entry["keep_tsv_files"] = _flag(all_params, "percolator", "keep-tsv-files", "false")

        # Write sub mzML probability threshold: TabRun.getSubMzmlProbThreshold() (TabRun.java:495-496), the
        # spinner tab-run.sub_mzml_prob_threshold (TabRun.java:230-231, initial value 0.5). The workflow text
        # is kept; writeSubMzmlArgv prints it as FragPipe does (float, probabilityThreshold + "",
        # CmdWriteSubMzml.java:126).
        if tool_name == "writesubmzml" and run_flag:
            entry["prob_threshold"] = all_params.get("tab-run", {}).get("sub_mzml_prob_threshold", "0.5").strip()

        # PeptideProphet per-group --combine (CmdPeptideProphet.java:190-202, :302-336)
        if tool_name == "peptideprophet" and run_flag:
            entry["combine"] = peptideprophet_combine(all_params)

        # SAINTexpress: cmd-opts is ONE argv element (CmdSaintExpress.java:85-89)
        if tool_name == "saintexpress" and run_flag:
            entry["cmd_opts"] = all_params.get("saintexpress", {}).get("cmd-opts", "").strip()

        # DIA-Umpire Q tiers added to the LC-MS list (UmpirePanel.generateQ1/2/3,
        # CmdUmpireSe.java:170-184). A tier missing from the workflow keeps the
        # umpire-se template value, which is true for Q1, Q2 and Q3 (diaumpire_se.params).
        if tool_name == "diaumpire":
            umpire_params = all_params.get("diaumpire", {})
            entry["q_tiers"] = [
                tier
                for tier in ("Q1", "Q2", "Q3")
                if umpire_params.get(tier, "true").strip().lower() == "true"
            ]

        # gen_con_spec_lib.py user RT/IM calibration files (CmdSpecLibGen.java:147-152), staged
        # by FRAGPIPE_SPECLIB; "" when the RT/IM option is a keyword.
        if tool_name == "speclibgen" and run_flag:
            settings = speclibgen_settings(all_params.get("speclibgen", {}), all_params)
            entry["rt_cal_file"] = settings["rt_cal_file"]
            entry["im_cal_file"] = settings["im_cal_file"]
            # gen_con_spec_lib.py argv parts (CmdSpecLibGen.java:137-202) as values, for SPECLIBGEN's argv builder
            entry["argv_parts"] = {
                key: settings[key] for key in ("rt_cal", "im_cal", "convert_args", "library_args", "intermediate")
            }

        # DIA-NN runtime settings for DIANNFRAGPIPE (CmdDiann.java:100-547): args holds
        # --qvalue..--matrices, args_post --smart-profiling..extra options; the protein-inference
        # flags and the post-steps (PlexDiaHelper, Propagation, DiannToMsstats, FragReporter) need
        # the fields below.
        if tool_name == "diann" and run_flag:
            diann_params = all_params.get("diann", {})
            entry["args_post"] = generate_diann_post_args(diann_params)
            entry["cmd_opts"] = diann_settings(diann_params)["cmd_opts"]
            entry.update(
                {
                    key: value
                    for key, value in diann_settings(diann_params).items()
                    if key
                    in (
                        "qvalue",
                        "run_specific_protein_qvalue",
                        "redo_protein_inference",
                        "generate_msstats",
                        "library",
                        "run_plex",
                        "light",
                        "medium",
                        "heavy",
                        "mod_tag",
                        "min_site_prob",
                        "report_levels",
                    )
                }
            )


        # MBG --nopasef: TabWorkflow.getInputDataType() falls back to the "IM-MS" radio button
        # (workflow.input.data-type.im-ms) when no LC-MS file is .d or .raw (TabWorkflow.java:521-538)
        if tool_name == "mbg" and run_flag:
            entry["im_ms"] = (
                all_params.get("workflow", {}).get("input.data-type.im-ms", "false").lower() == "true"
            )

        # O-Pair spectra list (CmdOPair.java:107-114): calibrated mzML, or per file uncalibrated/input
        if tool_name == "opair" and run_flag:
            entry["spectra"] = opair_spectra_kind(all_params)
            # user oxonium file, staged by FRAGPIPE_GLYCO ('' = FragPipe's default in the tools dir)
            entry["oxonium_file"] = opair_oxonium_file(all_params.get("opair", {}))

        # Single source of the decoy tag for every downstream module
        # (read via getDecoyTag() in subworkflows/local/fragpipe_utils.nf).
        if tool_name == "database":
            entry["decoy_tag"] = get_decoy_tag(all_params)

        result[tool_name] = entry

    # PairScans: one run per LC-MS file before O-Pair when O-Pair pairs scans
    # (FragpipeRun.java:2067-2081: oPairPanel.isRun() && !isSingleScanType()).
    opair_params = all_params.get("opair", {})
    pairscans_run = (
        result["opair"]["run"]
        and _opair_value(opair_params, "single_scan_type").lower() != "true"
    )
    result["pairscans"] = {
        "run": pairscans_run,
        "args": generate_pairscans_args(opair_params) if pairscans_run else [],
        "config_type": "argv",
    }

    # AppendFile: PTM-S glycan masses appended to IonQuant's modmasses_ionquant.txt
    # (FragpipeRun.java:1822-1830: ptmsGlycanPanel.isRun() && (quantPanelLabelfree.isRunIonQuant()
    # || (tmtiPanel.isRun() && tmtiPanel.getIntensityExtractionTool() == 0))).
    lfq_ionquant = label_free_ionquant_runs(all_params)
    tmti_params_af = all_params.get("tmtintegrator", {})
    tmt_ionquant = (
        tmti_params_af.get("run-tmtintegrator", "false").lower() == "true"
        and tmti_params_af.get("extraction_tool", "0") in ("0", "", "IonQuant")
    )
    glycan_assign = all_params.get("ptmshepherd", {}).get("run_glyco_mode", "false").lower() == "true"
    result["appendfile"] = {
        "run": glycan_assign and (lfq_ionquant or tmt_ionquant),
        "args": "",
        "config_type": "cli",
    }

    # When TMTIntegrator is enabled, add TMT-specific IonQuant two-pass entries.
    # This replicates FragPipe's two-pass IonQuant pattern (FragpipeRun.java:1978-2031):
    #   Pass 1 (MS1): --perform-ms1quant 1, adds precursor intensity to psm.tsv
    #   Pass 2 (Isobaric): --perform-isoquant 1, adds TMT reporter ion intensities
    # Both passes use the same base IonQuant params from the workflow file.
    #
    # IMP-6: Check intensityExtractionTool before generating TMT IonQuant entries.
    # FragPipe checks tmtiPanel.getIntensityExtractionTool() (FragpipeRun.java:1972, 2032):
    #   0 = IonQuant (default), 1 = Philosopher Freequant+Labelquant, 2 = skip
    tmti_run = result.get("tmtintegrator", {}).get("run", False)
    if tmti_run:
        tmti_params = all_params.get("tmtintegrator", {})
        extraction_tool = get_extraction_tool_index(all_params)
        result["tmtintegrator"]["extraction_tool"] = extraction_tool

        # FragpipeRun.java:1913-1916: TMT-Integrator only refuses "Skip extraction" when Report
        # runs; without Report the group psm.tsv files TMT-Integrator reads are never written.
        if extraction_tool == 2:
            if result.get("filter", {}).get("run", False):
                raise ValueError(
                    "'Intensity extraction tool' in 'Quant (Isobaric)' tab was set to 'Skip "
                    "extraction. Run TMT-Integrator only'. Please change it to 'IonQuant'."
                )
            raise ValueError(
                "'Skip extraction. Run TMT-Integrator only' reads psm.tsv files of an earlier run, "
                "which this pipeline does not take; set the intensity extraction tool to IonQuant "
                "or Philosopher."
            )

        # Philosopher MSstats (tmtintegrator.philosopher-msstats, TmtiPanel.java:369-370):
        # FragpipeRun.java:1768-1773 refuses it with IonQuant extraction; with Philosopher
        # extraction TMT-Integrator writes MSstatsTMT_annotation.csv (CmdTmtIntegrator.java:105-140).
        philosopher_msstats = (
            tmti_params.get("philosopher-msstats", "false").strip().lower() == "true"
        )
        if philosopher_msstats and extraction_tool == 0:
            raise ValueError(
                "Generate MSstats files (using Philosopher) was enabled but Philosopher was not "
                "selected as Intensity extraction tool."
            )
        result["tmtintegrator"]["philosopher_msstats"] = philosopher_msstats

        # CmdTmtIntegrator.java:245-253: a min site probability needs a mod tag
        # (getModTag() is the raw text, so only a blank tag is refused)
        if (
            float(_tmti_widget_value("min_site_prob", tmti_params.get("min_site_prob"))) >= 0
            and not tmti_params.get("mod_tag", "none").strip()
        ):
            raise ValueError(
                "TMT Integrator Min mod site probability is specified, but no mods are provided. "
                "Please set Min site prob to -1 or specify mods."
            )

        # FragpipeRun.java:2033-2047: TMT with Philosopher extraction refuses label-free FreeQuant
        if extraction_tool == 1 and result.get("freequant", {}).get("run", False):
            raise ValueError(
                "Both FreeQuant and TMT-Integrator were enabled. If you want to perform TMT "
                "analysis, please disable FreeQuant. If you want to perform LFQ analysis, "
                "please disable TMT-Integrator."
            )

        if extraction_tool == 0:
            modmasses = _parse_mod_masses(all_params)
            result["ionquant_ms1"] = {
                "run": True,
                "args": _generate_tmt_ionquant_ms1_args(all_params),
                "config_type": "cli",
                "modmasses": modmasses,
            }
            result["ionquant_isobaric"] = {
                "run": True,
                "args": _generate_tmt_ionquant_isobaric_args(tmti_params, all_params),
                "config_type": "cli",
                "modmasses": modmasses,
            }

    # FPOP TMT: TMT-Integrator runs a second time for the unmodified peptides
    # (FragpipeRun.java:2195-2209, output tmt-report-unmod). The DIA data type gate is checked
    # by the pipeline, which knows the input files (checkFpopTools in fragpipe_utils.nf).
    fpop_tmt = result.get("fpop", {}).get("run", False) and _java_boolean(
        all_params.get("fpop", {}).get("fragpipe.fpop.fpop-tmt", "false")
    ) == "true"
    result["tmtintegrator_unmod"] = {
        "run": tmti_run and fpop_tmt,
        # tmt-integrator-conf-unmod.yml text (CmdTmtIntegrator.CONFIG_FN_2), as tmtintegrator.args holds the first run's
        "args": generate_tmtintegrator_config(all_params.get("tmtintegrator", {}), True) if tmti_run and fpop_tmt else "",
        "config_type": "params_file",
    }

    # Multi-experiment report and iProphet (FragpipeRun.java:1782-1809)
    result.update(abacus_iprophet_entries(all_params, result))

    return result


def generate_tool_config(
    tool_prefix: str,
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]] = None,
) -> Tuple[str, str]:
    """
    Generate config file content for a specific tool.

    Returns:
        Tuple of (filename, content)
    """
    filename = f"{tool_prefix}.config"

    if tool_prefix == "ionquant":
        content = generate_ionquant_config(params, all_params)
    elif tool_prefix == "percolator":
        content = generate_percolator_config(params, get_run_flag(all_params, "percolator"))
    elif tool_prefix == "msbooster":
        content = generate_msbooster_config(params, all_params, get_run_flag(all_params, "msbooster"))
    elif tool_prefix == "msfragger":
        content = generate_msfragger_config(params, all_params)
        filename = "fragger.params"  # Native MSFragger params file format
    elif tool_prefix in (
        "phi-report",
        "peptide-prophet",
        "protein-prophet",
        "database",
    ):
        content = generate_philosopher_config(params, tool_prefix, all_params)
    elif tool_prefix == "tmtintegrator":
        content = generate_tmtintegrator_config(params)
        filename = "tmtintegrator.yml"  # YAML format required by TMTIntegrator
    elif tool_prefix == "ptmprophet":
        content = generate_ptmprophet_config(params, all_params)
    elif tool_prefix == "ptmshepherd":
        content = generate_ptmshepherd_config(params, all_params)
        filename = "shepherd.config"  # Java properties format required by PTMShepherd
    elif tool_prefix == "freequant":
        # C2: Use dedicated freequant config generator with correct CLI flag names
        content = generate_freequant_config(params)
    elif tool_prefix == "labelquant":
        content = generate_labelquant_config(params, all_params)
    # New tool-specific generators
    elif tool_prefix == "diann":
        content = generate_diann_config(params)
    elif tool_prefix == "speclibgen":
        content = generate_speclibgen_config(params, all_params)
        filename = "speclibgen.params"
    elif tool_prefix == "saintexpress":
        content = generate_saintexpress_config(params)
    elif tool_prefix == "crystalc":
        content = generate_crystalc_config(params)
        filename = "crystalc.params"  # Java properties format
    elif tool_prefix == "tab-run":
        content = generate_tabrun_config(params)
    elif tool_prefix == "quantitation":
        content = generate_quantitation_config(params)
    else:
        content = generate_generic_config(params, tool_prefix, file_mode=True)

    return filename, content


def write_tool_configs(
    outdir: Path,
    tool_prefix: str,
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]],
    data_types: List[str],
) -> None:
    """Write a tool's config file(s); MSFragger gets one params file per data type, written
    byte for byte and named as FragPipe writes them (msfragger_params_file_name)."""
    if tool_prefix == "msfragger":
        configs = generate_msfragger_configs(params, all_params)
        files = {msfragger_params_file_name(t, data_types): configs[f] for t, f in MSFRAGGER_PARAMS_FILES.items()}
    elif tool_prefix == "tmtintegrator":
        # tmtintegrator_unmod.yml: FPOP TMT's second run (CmdTmtIntegrator.CONFIG_FN_2)
        files = {
            "tmtintegrator.yml": generate_tmtintegrator_config(params),
            "tmtintegrator_unmod.yml": generate_tmtintegrator_config(params, True),
        }
    elif tool_prefix == "diaumpire":
        # the umpire-se.params text as FragPipe writes it, final newline included
        files = {"diaumpire.params": generate_diaumpire_config(params)}
    else:
        filename, content = generate_tool_config(tool_prefix, params, all_params)
        files = {filename: content + "\n"}
    for filename, content in files.items():
        output_path = outdir / filename
        with open(output_path, "w", encoding="utf-8", newline="\n") as f:
            f.write(content)
        print(f"Generated: {output_path}")


# FragPipe 24.0 settings that msproteomics does not support. FragPipe aborts the
# run when a stage's configure() returns false (FragpipeRun.java:2529-2536); msproteomics stops
# the same way instead of finishing without the stage's outputs.
#   (prefix, key, unsupported value, message)
UNSUPPORTED_STAGES = (
    (
        # Digest only: TabMsfragger.isRun() is isEnabledAndChecked(checkRun) && !isDigestOnly()
        # (TabMsfragger.java:1932-1934), so FragPipe skips the search exactly as with
        # run-msfragger=false (a FragPipe 24.0 --dry-run plans no MSFragger step) and hands the
        # expected MSFragger outputs of the work directory to later stages (FragpipeRun.java:1519-1521).
        "msfragger",
        "digest-only",
        "true",
        "MSFragger digest only was enabled (msfragger.digest-only=true): FragPipe then skips the "
        "MSFragger search and reuses the search results of a previous run in its work directory, "
        "but msproteomics starts from LC-MS files and cannot take a previous work directory as input. "
        "Please disable digest only in the workflow.",
    ),
    (
        # With MSFragger off FragPipe still hands the expected MSFragger outputs of the work
        # directory to every later stage (FragpipeRun.java:1519-1521, outside `if isRun()`), i.e.
        # it reuses a previous run's search results. msproteomics starts from LC-MS files only.
        # "Run MSFragger" defaults to true (TabMsfragger.java:570-571).
        "msfragger",
        "run-msfragger",
        "false",
        "MSFragger is disabled (msfragger.run-msfragger=false). FragPipe then reuses the "
        "search results of a previous run in its work directory, but msproteomics starts "
        "from LC-MS files and cannot take a previous work directory as input. "
        "Please enable MSFragger in the workflow.",
    ),
    (
        # FragPipe on Linux refuses the run: SkylinePanel.getSkylinePath() returns null when
        # !isWindows() (SkylinePanel.java:332-345), and CmdSkyline.configure() logs this
        # message and returns false (CmdSkyline.java:80-89).
        "skyline",
        "run-skyline",
        "true",
        "Skyline was requested but no installed version could be found! "
        "Please verify that it is installed or located at the specified path. "
        "(skyline.run-skyline=true: FragPipe runs Skyline on Windows only, "
        "and msproteomics does not run it. Please disable Skyline in the workflow.)",
    ),
    (
        # CmdTransferLearning needs a credential file with the URL and API key of an external
        # prediction server (CmdTransferLearning.java:123-147); CmdMsfraggerDigest runs only
        # with transfer learning (FragpipeRun.java:2255).
        "transfer-learning",
        "run-transfer-learning",
        "true",
        "Transfer learning was enabled (transfer-learning.run-transfer-learning=true), "
        "but msproteomics does not run Transfer learning or its MSFragger digest: it needs an "
        "external prediction server (credential file with URL and API key). "
        "Please disable Transfer learning in the workflow.",
    ),
)


# Tools that read Philosopher report's <group>/psm.tsv: CmdIonquant.java:318, CmdTmtIntegrator.java:295,
# PtmshepherdParams.java:85, CmdMBGMatch.java:139, CmdOPair.java:176, CmdSpecLibGen.java:135.
PSM_TSV_READERS = ("ionquant", "tmtintegrator", "ptmshepherd", "mbg", "opair", "speclibgen")


def apply_run_tab_listeners(all_params: Dict[str, Dict[str, str]]) -> None:
    """
    The Run tab's UI side effects on a loaded workflow, in place.

    Checking "Write sub mzML" checks the report panel's "Print decoys" (TabRun.java:240-243):

        uiCheckWriteSubMzml.addItemListener(e -> {
          uiSpinnerProbThreshold.setEnabled(uiCheckWriteSubMzml.isSelected());
          reportPanel.setPrintDecoys(uiCheckWriteSubMzml.isSelected());
        });

    FragPipe loads the Run tab after the report panel, so a workflow with tab-run.write_sub_mzml=true
    runs Philosopher Report with --decoys (reportPanel.isPrintDecoys(), FragpipeRun.java:1778) whatever
    phi-report.print-decoys says; FragPipe's own load and save (assets/fragpipe_ui_state) writes
    phi-report.print-decoys=true for assets/LFQ-MBR.workflow with write_sub_mzml=true.

    @param all_params Parsed workflow (prefix -> key -> value), modified in place
    """
    if _flag(all_params, "tab-run", "write_sub_mzml"):
        all_params.setdefault("phi-report", {})["print-decoys"] = "true"


def check_write_sub_mzml(all_params: Dict[str, Dict[str, str]], data_types: List[str]) -> None:
    """
    Fail like CmdWriteSubMzml.configure (CmdWriteSubMzml.java:72-80) when Write sub mzML is on:

        if (!isRunMSFragger) { ... "Write sub mzML was enabled but MSFragger was not." ... return false; }
        if (hasDia || hasGpfDia || hasDiaLib || hasDdaPlus) { ... "there are non-DDA data types." ... return false; }

    (hasDataType("DIA"), ("GPF-DIA"), ("DIA-Lib"), ("DDA+") over the manifest, FragpipeRun.java:2400-2403).
    A split database search fails too: FragPipe appends --no-ptm-score-in-pin after the LC-MS files
    (CmdMsfragger.java:586, :628-630), and msfragger_pep_split.py takes every argument after the params
    file as an LC-MS file (msfragger_pep_split.py:67, :72, :111) and passes it to MSFragger's calibration
    search (:502-514). FragPipe 24.0 headless (licensed image, LFQ-MBR with msfragger.misc.slice-db=2 and
    tab-run.write_sub_mzml=true, 2026-10-11) stops there: MSFragger prints
    "Could not find MS file at <workdir>/--no-ptm-score-in-pin" and exits 1, msfragger_pep_split.py raises
    AssertionError (line 514) and FragPipe stops ("Process returned non-zero exit code, stopping").

    @param all_params Parsed workflow after the manifest adjustments
    @param data_types Data types of the manifest files
    @raises SystemExit with FragPipe's message
    """
    if not _flag(all_params, "tab-run", "write_sub_mzml"):
        return
    if not get_run_flag(all_params, "msfragger"):
        sys.exit(
            "ERROR: Write sub mzML was enabled but MSFragger was not.\n"
            "Please either disable Write sub mzML in the Run tab or enable MSFragger in the MSFragger tab."
        )
    if {"DIA", "GPF-DIA", "DIA-Lib", "DDA+"} & set(data_types):
        sys.exit(
            "ERROR: Write sub mzML was enabled but there are non-DDA data types.\n"
            "Please disable Write sub mzML in the Run tab."
        )
    if msfragger_num_slices(all_params.get("msfragger", {}), all_params) > 1:
        sys.exit(
            "ERROR: Write sub mzML was enabled with a split database search (msfragger.misc.slice-db > 1). "
            "FragPipe appends --no-ptm-score-in-pin after the LC-MS files (CmdMsfragger.java:628-630) and "
            "msfragger_pep_split.py reads it as an LC-MS file (msfragger_pep_split.py:67-72), so FragPipe's "
            "MSFragger stops with 'Could not find MS file at <workdir>/--no-ptm-score-in-pin'. "
            "Please set msfragger.misc.slice-db to 1 or disable Write sub mzML in the Run tab."
        )


def check_supported_stages(all_params: Dict[str, Dict[str, str]]) -> None:
    """
    Fail when the workflow sets a FragPipe stage switch that msproteomics does not support.

    @param all_params Parsed workflow (parse_workflow_file)
    @raises SystemExit with every enabled unsupported stage's message
    """
    errors = [
        f"ERROR: {message}"
        for prefix, key, value, message in UNSUPPORTED_STAGES
        if all_params.get(prefix, {}).get(key, "").lower() == value
    ]
    report_panel_runs = _validation_tab_runs(all_params) and _run_checkbox(all_params, "filter")
    # FreeQuant without the Report panel: FragPipe runs Philosopher report (FragpipeRun.java:1767
    # reportPanel.isRun() || quantPanelLabelfree.isRunFreeQuant()) but no DbAnnotate (:1690) and no
    # Filter (:1701), i.e. a report of an unfiltered workspace; a FragPipe 24.0 --dry-run of
    # TMT-labelcheck-304 with FreeQuant and run-report=false plans FreeQuant, PhilosopherReport and
    # PhilosopherAbacus but no DbAnnotate or Filter. msproteomics runs report after filter only.
    freequant_without_report = get_run_flag(all_params, "freequant") and not report_panel_runs
    if freequant_without_report:
        errors.append(
            "ERROR: FreeQuant runs without the Philosopher report panel (phi-report.run-report=false "
            "or run-validation-tab=false): FragPipe then runs Philosopher report on a workspace that "
            "has no database annotation and no FDR filter, which msproteomics does not support. "
            "Please enable Generate reports (phi-report.run-report=true) or use IonQuant."
        )
    # Report panel off: FragPipe still plans the tools that read Philosopher report's <group>/psm.tsv,
    # which then never exists (a FragPipe 24.0 --dry-run of LFQ-MBR with run-validation-tab=false plans
    # MSFragger and IonQuant only), so its run fails there. (With FreeQuant Report does run; refused above.)
    if not report_panel_runs and not freequant_without_report:
        readers = [tool for tool in PSM_TSV_READERS if get_run_flag(all_params, tool)]
        if readers:
            errors.append(
                "ERROR: Philosopher report does not run (phi-report.run-report=false or run-validation-tab=false), "
                f"so no psm.tsv is written, but these tools read it: {', '.join(readers)}. "
                "Please enable Generate reports and Run validation tools, or disable these tools in the workflow."
            )
    if errors:
        sys.exit("\n".join(errors))


# CmdMsfragger.pattern (CmdMsfragger.java:65): the allowed sites of an enabled variable modification
VAR_MOD_SITES = re.compile(r"[A-Znc\[\]\^\*]+")
# PtmshepherdPanel.validateForm (PtmshepherdPanel.java:816-829): "<name>:<mass>, <name>:<mass>, ..."
PTMSHEPHERD_VARMODS = re.compile(r"[^\s]+:-?\d+(?:\.\d+)?(?:\s*,\s*[^\s]+:-?\d+(?:\.\d+)?)*")
# IonQuant values FragPipe adds as ONE argv element (CmdIonquant.java:303-304); the pipeline splits the options on spaces
IONQUANT_SINGLE_ARGUMENTS = ("light", "medium", "heavy", "excludemods", "formula")


def workflow_refusals(all_params: Dict[str, Dict[str, str]]) -> list:
    """
    The checks after which FragPipe aborts a run (or its tool's configure() returns false) on values that
    msproteomics would otherwise run. Each message is FragPipe's text verbatim; a second line names the
    workflow key and value. A key that is absent keeps its panel default.

    - MSBooster: open search, neither prediction on, Koina URL missing (CmdMSBooster.java:106-131)
    - Percolator / PeptideProphet need the MSFragger output formats (FragpipeRun.java:252-266)
    - labile mode with remainder masses needs localize_delta_mass (CmdMsfragger.java:418-427)
    - an enabled variable modification with invalid sites (CmdMsfragger.java:331-335)
    - PTM-Shepherd variable mass list (PtmshepherdPanel.java:816-829; FragpipeRun.java:2113-2122 checks it with
      PTM-Shepherd off too)
    - IonQuant text values with whitespace (CmdIonquant.java:303-304 adds each as one argv element; this is
      msproteomics' own limit: it splits the options on spaces)

    @param all_params Parsed workflow (parse_workflow_file)
    @return the error messages
    """
    errors = []
    msfragger = all_params.get("msfragger", {})

    # FragpipeRun.java:1462 configures MSFragger only when it runs, and :1558 runs MSBooster only on MSFragger's pepXML
    runs_msfragger = get_run_flag(all_params, "msfragger")

    if runs_msfragger and get_run_flag(all_params, "msbooster"):
        msbooster = all_params.get("msbooster", {})
        # TabMsfragger.isOpenSearch() reads the precursor spinners and combo (TabMsfragger.java:477): the
        # panel defaults when the workflow has no such keys
        form = msfragger_values(all_params)
        lower = float(form["precursor_mass_lower"])
        upper = float(form["precursor_mass_upper"])
        if form["precursor_mass_units"] == "0" and lower < -3 and upper > 3:
            errors.append(
                "MSBooster is incompatible with open search. Please disable MSBooster."
                f"\n  (msfragger.precursor_mass_units=0, precursor_mass_lower={lower:g}, precursor_mass_upper={upper:g})"
            )
        if msbooster_flag(msbooster, "predict-rt") == "false" and msbooster_flag(msbooster, "predict-spectra") == "false":
            errors.append(
                "At least one of <b>Predict RT</b> or <b>Predict spectra</b> must be enabled if you want to run MSBooster."
                "\n  (msbooster.predict-rt=false, msbooster.predict-spectra=false)"
            )
        koina = msbooster.get("koina-url", "").strip()
        models = [msbooster_model_ui_name(msbooster, key) for key in ("rt-model", "spectra-model", "im-model")]
        if not koina and any(model != "DIA-NN" for model in models):
            errors.append(
                "Koina URL is required for non DIA-NN models.\nPlease go to <b>Validation</b> tab and adjust the settings."
                f"\n  (msbooster.koina-url is empty, models {models})"
            )
        if not koina and any(
            msbooster_flag(msbooster, key) == "true"
            for key in ("find-best-rt-model", "find-best-spectra-model", "find-best-im-model")
        ):
            errors.append(
                "Koina URL is required for <b>Find best RT/spectral/IM model</b>.\nPlease go to <b>Validation</b> tab and adjust the settings."
                "\n  (msbooster.koina-url is empty)"
            )

    output_format = msfragger.get("output_format", "").strip().lower()
    if output_format:
        detail = f"\n  (msfragger.output_format={msfragger['output_format'].strip()})"
        if get_run_flag(all_params, "percolator") and ("pin" not in output_format or "pepxml" not in output_format):
            errors.append("Percolator was enabled but MSFragger's output formats did not contain pin or pepXML." + detail)
        if get_run_flag(all_params, "peptideprophet") and "pepxml" not in output_format:
            errors.append("PeptideProphet was enabled but MSFragger's output formats did not contain pepXML." + detail)

    # ui.translations:23-24: the saved "standard" is the combo's "off"; localize_delta_mass is
    # Boolean.parseBoolean (TabMsfragger.java:284), and an absent key keeps the panel's 0 (fragger_closed.params:38)
    if (
        runs_msfragger
        and msfragger.get("labile_search_mode", "off").strip() not in ("off", "standard")
        and msfragger.get("remainder_fragment_masses", "") != ""
        and msfragger.get("localize_delta_mass", "").strip().lower() != "true"
    ):
        errors.append(
            "Fragment remainder masses are specified for MSFragger, but localize_delta_mass is disabled. "
            "Please enable it to continue (it is required for fragment remainder masses)."
            f"\n  (msfragger.labile_search_mode={msfragger['labile_search_mode'].strip()}, "
            f"remainder_fragment_masses={msfragger['remainder_fragment_masses']})"
        )

    for entry in msfragger.get("table.var-mods", "").split(";") if runs_msfragger else []:
        fields = [field.strip() for field in entry.split(",")]
        if len(fields) >= 3 and fields[2].lower() == "true" and not VAR_MOD_SITES.fullmatch(fields[1]):
            errors.append(
                f"Invalid variable modification sites: {fields[1]}.\nPlease check the tooltip or documentations."
                f"\n  (msfragger.table.var-mods entry {entry.strip()!r})"
            )

    varmods = all_params.get("ptmshepherd", {}).get("varmod_masses", "").strip()
    if varmods and not PTMSHEPHERD_VARMODS.fullmatch(varmods):
        errors.append(
            "There are errors in PTM-Shepherd configuration panel on Report tab."
            f"\n  (ptmshepherd.varmod_masses={varmods!r} does not match <name>:<mass>[,<name>:<mass>...])"
        )

    if get_run_flag(all_params, "ionquant"):
        ionquant = all_params.get("ionquant", {})
        # CmdIonquant.java:296-301: light, medium and heavy are skipped unless the use-labeling radio is on
        labeling = ionquant.get("use-labeling", "false").strip() != "false"
        for key in IONQUANT_SINGLE_ARGUMENTS:
            if key in ("light", "medium", "heavy") and not labeling:
                continue
            if re.search(r"\s", ionquant.get(key, "").strip()):
                errors.append(
                    f"ionquant.{key}={ionquant[key].strip()!r} contains whitespace: FragPipe passes it to IonQuant "
                    "as one argument, but msproteomics splits the IonQuant options on spaces."
                )
    return errors


def check_workflow_refusals(all_params: Dict[str, Dict[str, str]]) -> None:
    """Fail with every message of workflow_refusals (FragPipe aborts on each of them)."""
    errors = workflow_refusals(all_params)
    if errors:
        sys.exit("\n".join(f"ERROR: {message}" for message in errors))


def main():
    global FRAGGER_TEMPLATE, DIAUMPIRE_TEMPLATE
    parser = argparse.ArgumentParser(
        description="Parse FragPipe workflow files and generate tool-specific configs"
    )
    parser.add_argument(
        "--workflow", "-w", required=True, help="Path to FragPipe .workflow file"
    )
    parser.add_argument(
        "--outdir", "-o", default=".", help="Output directory for config files"
    )
    parser.add_argument("--tool", "-t", help="Generate config for specific tool only")
    parser.add_argument(
        "--list-tools",
        action="store_true",
        help="List all tools found in workflow file",
    )
    parser.add_argument(
        "--data-types",
        required=True,
        help="Comma-separated FragPipe data types of the manifest files (DDA, DDA+, DIA, DIA-Quant, DIA-Lib, GPF-DIA)",
    )
    parser.add_argument(
        "--output-json",
        action="store_true",
        help="Output JSON with run flags and args for all tools",
    )

    parser.add_argument(
        "--configs-outdir",
        help="With --output-json: also write the tool config files to this directory in the same run",
    )
    parser.add_argument(
        "--fragger-template",
        type=Path,
        required=True,
        help="FragPipe 24.0's fragger_closed.params (the module's assets/)",
    )
    parser.add_argument(
        "--diaumpire-template",
        type=Path,
        required=True,
        help="FragPipe 24.0's diaumpire_se.params (the module's assets/)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Thread count of the run; overrides workflow.threads like FragPipe headless --threads (Fragpipe.java:508-514)",
    )

    args = parser.parse_args()
    if args.configs_outdir and not args.output_json:
        parser.error("--configs-outdir needs --output-json")
    FRAGGER_TEMPLATE = args.fragger_template
    DIAUMPIRE_TEMPLATE = args.diaumpire_template

    # Parse workflow file
    all_params = parse_workflow_file(args.workflow)
    if args.threads is not None:
        all_params.setdefault("workflow", {})["threads"] = str(args.threads)
    apply_manifest_data_types(all_params, args.data_types.split(","))
    apply_run_tab_listeners(all_params)
    check_supported_stages(all_params)
    check_write_sub_mzml(all_params, args.data_types.split(","))
    check_workflow_refusals(all_params)

    if args.list_tools:
        print("Tools found in workflow file:")
        for tool in sorted(all_params.keys()):
            if tool != "_global":
                print(f"  {tool}: {len(all_params[tool])} parameters")
        return 0

    # JSON output mode
    if args.output_json:
        # the config files start from the parsed workflow as a run without --output-json does
        config_params = copy.deepcopy(all_params) if args.configs_outdir else None
        check_glyco_hard_stops(all_params)
        json_output = generate_json_output(all_params, args.data_types.split(","))
        # Create output directory
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        # Write JSON file
        output_path = outdir / "tool_configs.json"
        with open(output_path, "w") as f:
            json.dump(json_output, f, indent=2)
        print(f"Generated: {output_path}")
        if not config_params:
            return 0
        all_params = config_params
        args.outdir = args.configs_outdir

    # Create output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Generate configs
    if args.tool:
        # Single tool
        if args.tool not in all_params:
            print(
                f"Error: Tool '{args.tool}' not found in workflow file", file=sys.stderr
            )
            return 1

        write_tool_configs(outdir, args.tool, all_params[args.tool], all_params, args.data_types.split(","))
    else:
        # All tools
        for tool_prefix, params in all_params.items():
            if tool_prefix == "_global":
                continue
            # FragPipe builds PTM-Shepherd's config only when PTM-Shepherd runs (FragpipeRun.java:2124-2126,
            # cmdPtmshepherd.isRun()); it reads MSFragger values a workflow without PTM-Shepherd need not have
            if tool_prefix == "ptmshepherd" and not get_run_flag(all_params, "ptmshepherd"):
                continue

            write_tool_configs(outdir, tool_prefix, params, all_params, args.data_types.split(","))

    return 0


if __name__ == "__main__":
    sys.exit(main())
