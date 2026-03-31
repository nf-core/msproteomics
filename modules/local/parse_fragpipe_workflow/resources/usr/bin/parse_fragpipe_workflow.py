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
   - Used by PARSE_FRAGPIPE_WORKFLOW module for dynamic workflow execution
   - Every tool prefix in the workflow file is handled; no silent parameter drops

Usage:
    python parse_fragpipe_workflow.py --workflow input.workflow --outdir configs/
    python parse_fragpipe_workflow.py --workflow input.workflow --tool ionquant
    python parse_fragpipe_workflow.py --workflow input.workflow --output-json --outdir configs/
"""

import argparse
import json
import re
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, Tuple, Any


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
            "use-corr",
            "use-spectra",
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
    "saintexpress": ("saintexpress", "run-saint-express"),
    "skyline": ("skyline", "run-skyline"),
    "metaproteomics": ("metaproteomics", "run-metaproteomics"),
    "transferlearning": ("transfer-learning", "run-transfer-learning"),
    # Database (always runs if present)
    "database": ("database", None),  # No run flag, always runs
    # GUI-only tab settings (included for completeness; not standalone tools)
    "tabrun": ("tab-run", None),  # GUI tab; run derived from write_sub_mzml
    "quantitation": (
        "quantitation",
        "run-label-free-quant",
    ),  # GUI quant panel; controls IonQuant/FreeQuant
}


def _unescape_java_properties(value: str) -> str:
    """Unescape Java Properties file escape sequences in values.

    FragPipe .workflow files use Java Properties format where:
      \\= → =    \\: → :    \\\\ → \\    \\n → newline    \\t → tab
    """
    result = []
    i = 0
    while i < len(value):
        if value[i] == "\\" and i + 1 < len(value):
            next_char = value[i + 1]
            if next_char == "=":
                result.append("=")
                i += 2
            elif next_char == ":":
                result.append(":")
                i += 2
            elif next_char == "\\":
                result.append("\\")
                i += 2
            elif next_char == "n":
                result.append("\n")
                i += 2
            elif next_char == "t":
                result.append("\t")
                i += 2
            else:
                result.append(value[i])
                i += 1
        else:
            result.append(value[i])
            i += 1
    return "".join(result)


def parse_workflow_file(workflow_path: str) -> Dict[str, Dict[str, str]]:
    """
    Parse a FragPipe workflow file into a nested dictionary.

    Returns:
        Dict mapping tool_prefix -> {param_name -> value}
    """
    params = defaultdict(dict)

    with open(workflow_path, "r") as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith("#"):
                continue

            # Parse key=value
            if "=" not in line:
                continue

            key, value = line.split("=", 1)
            key = key.strip()
            value = _unescape_java_properties(value.strip())

            # Split by first dot to get tool prefix
            if "." in key:
                prefix, param = key.split(".", 1)
                params[prefix][param] = value
            else:
                # Top-level params (rare)
                params["_global"][key] = value

    return dict(params)


def _parse_mod_masses(all_params: Dict[str, Dict[str, str]]) -> str:
    """Extract enabled modification masses from MSFragger configuration.

    Replicates FragPipe's modmasses_ionquant.txt generation (FragpipeRun.java:1852-1855).
    Combines enabled variable mods, fixed mods, and mass offsets with |mass| > 0.01.

    Returns:
        Comma-separated string of unique modification masses (sorted), or empty string.
    """
    masses = set()
    msfragger = all_params.get("msfragger", {})

    # Parse variable mods: "15.9949,M,true,2; 42.0106,[^,true,1; ..."
    for table_key in ("table.var-mods", "table.fix-mods"):
        raw = msfragger.get(table_key, "")
        if not raw:
            continue
        for entry in raw.split(";"):
            entry = entry.strip()
            if not entry:
                continue
            parts = entry.split(",")
            if len(parts) < 3:
                continue
            try:
                mass = float(parts[0].strip())
            except ValueError:
                print(
                    f"WARNING: Skipping unparseable mass value '{parts[0].strip()}' in variable mod entry: {entry}",
                    file=sys.stderr,
                )
                continue
            enabled = parts[2].strip().lower() == "true"
            if enabled and abs(mass) > 0.01:
                masses.add(mass)

    # Parse mass offsets: "0" or "0/79.96633/..."
    raw_offsets = msfragger.get("mass_offsets", "")
    if raw_offsets:
        for token in raw_offsets.replace("/", " ").split():
            try:
                mass = float(token.strip())
                if abs(mass) > 0.01:
                    masses.add(mass)
            except ValueError:
                print(
                    f"WARNING: Skipping unparseable mass offset token '{token.strip()}' in mass_offsets: {raw_offsets}",
                    file=sys.stderr,
                )
                continue

    # IMP-2: For glycoproteomics workflows, FragPipe appends glycan masses from
    # O-Pair output to the modlist (FragpipeRun.java after glyco, CmdAppendFile).
    # We replicate this by reading opair.glycan-masses if the opair tool is enabled.
    # The glycan masses are stored in the workflow file as a space-separated list.
    opair_params = all_params.get("opair", {}) if isinstance(all_params, dict) else {}
    if opair_params.get("run-opair", "false").lower() == "true":
        glycan_masses_str = opair_params.get("glycan-masses", "")
        if glycan_masses_str:
            for token in glycan_masses_str.split():
                try:
                    mass = float(token.strip())
                    if abs(mass) > 0.01:
                        masses.add(mass)
                except ValueError:
                    print(
                        f"WARNING: Skipping unparseable glycan mass token '{token.strip()}' in glycan-masses",
                        file=sys.stderr,
                    )
                    continue

    if not masses:
        return ""
    return ",".join(f"{m:.10g}" for m in sorted(masses))


def generate_ionquant_config(
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]] = None,
) -> str:
    """Generate IonQuant CLI flags from workflow params.

    Note: Runtime parameters (threads, multidir, specdir, psm, filelist) are
    excluded - they are set by the module at runtime.

    Replicates FragPipe's programmatic parameter logic from CmdIonquant.java:
    - Default flags: perform-ms1quant, perform-isoquant, site-reports, msstats
      are always set explicitly (CmdIonquant.java:181-202), not from workflow file.
    - use-labeling: when false, suppresses --light/--medium/--heavy flags
      (CmdIonquant.java:296-301).
    - ionmobility: auto-detected from workflow.input.data-type.im-ms
      (CmdIonquant.java:197-198).
    """
    cli_parts = []

    # Check use-labeling flag. When false, suppress light/medium/heavy even if
    # they have values (CmdIonquant.java:296-301).
    use_labeling = params.get("use-labeling", "false").lower() == "true"

    # Map workflow params to CLI flags
    param_mapping = {
        "mbr": "mbr",
        "maxlfq": "maxlfq",
        "requantify": "requantify",
        "mztol": "mztol",
        "imtol": "imtol",
        "rttol": "rttol",
        "mbrmincorr": "mbrmincorr",
        "mbrrttol": "mbrrttol",
        "mbrimtol": "mbrimtol",
        "mbrtoprun": "mbrtoprun",
        "ionfdr": "ionfdr",
        "proteinfdr": "proteinfdr",
        "peptidefdr": "peptidefdr",
        "normalization": "normalization",
        "minisotopes": "minisotopes",
        "intensitymode": "intensitymode",
        "minscans": "minscans",
        "tp": "tp",
        "minfreq": "minfreq",
        "minions": "minions",
        "excludemods": "excludemods",
        "locprob": "locprob",
        "uniqueness": "uniqueness",
        "formula": "formula",
        "ionmobility": "ionmobility",
        "minexps": "minexps",
        "writeindex": "writeindex",
    }

    # Label params: only included when use-labeling=true (CmdIonquant.java:296-301)
    label_params = {
        "light": "light",
        "medium": "medium",
        "heavy": "heavy",
    }
    if use_labeling:
        param_mapping.update(label_params)

    # Note: threads, multidir excluded - set at runtime
    # Note: perform-ms1quant, perform-isoquant, isotol, isolevel, isotype,
    #   site-reports, msstats are set as defaults below (not from workflow file)

    # isolevel requires special translation: MS2→2, MS3→3, ZOOM-HR→4
    # (CmdIonquant.java:187-194). Numeric values pass through.
    isolevel_map = {"MS2": "2", "MS3": "3", "ZOOM-HR": "4"}

    for workflow_key, cli_flag in param_mapping.items():
        if workflow_key in params and params[workflow_key]:
            value = params[workflow_key]
            if value and value.lower() not in ("", "null", "none"):
                cli_parts.append(f"--{cli_flag} {value}")

    # Handle isolevel separately due to string-to-numeric translation
    if "isolevel" in params and params["isolevel"]:
        raw = params["isolevel"]
        if raw and raw.lower() not in ("", "null", "none"):
            numeric = isolevel_map.get(raw, raw)  # Pass through if already numeric
            cli_parts.append(f"--isolevel {numeric}")

    # FragPipe always sets these flags programmatically (CmdIonquant.java:181-202).
    # They are NOT in the workflow file — CmdIonquant's configure() method sets them
    # from its parameters: performMS1Quant=true, performIsobaricQuant=false for LFQ.
    # We add them as defaults if not already set by the workflow params above.
    #
    # site-reports: Derived from workflow context (FragpipeRun.java:1887 vs 1995/2028):
    #   LFQ path (ionquant.run-ionquant=true, no tmtintegrator) → site-reports=1
    #   TMT path (tmtintegrator.run-tmtintegrator=true)          → site-reports=0
    is_tmt = get_run_flag(all_params, "tmtintegrator") if all_params else False
    site_reports_val = "0" if is_tmt else "1"

    # msstats: FragPipe reads diann.generate-msstats from workflow file
    # (FragpipeRun.java:2324 → diannPanel.generateMsstats() → CmdIonquant.java:201)
    # Default is "1" (true) matching DiannPanel.java:239 UI default
    msstats_val = "1"
    if all_params:
        diann_params = all_params.get("diann", {})
        if diann_params.get("generate-msstats", "true").lower() == "false":
            msstats_val = "0"

    programmatic_defaults = {
        "perform-ms1quant": "1",  # CmdIonquant.java:181 (true for LFQ)
        "perform-isoquant": "0",  # CmdIonquant.java:183 (false for LFQ)
        "isotol": "20",  # MIN-1: CmdIonquant.java:98 (dummy, isoquant disabled)
        "isolevel": "2",  # MIN-1: CmdIonquant.java:99 (dummy, isoquant disabled)
        "isotype": "tmt10",  # MIN-1: CmdIonquant.java:100 (dummy, isoquant disabled)
        "intensitymode": "2",  # QuantPanelLabelfree.java:296 GUI default "auto"=2
        "site-reports": site_reports_val,  # CmdIonquant.java:199; context-dependent
        "msstats": msstats_val,  # CmdIonquant.java:201; from diann.generate-msstats
        "minexps": "1",  # CmdIonquant.java:180 (default minimum experiments)
        "writeindex": "0",  # CmdIonquant.java (index writing disabled by default)
    }
    for flag, default_val in programmatic_defaults.items():
        has_flag = any(f"--{flag}" in p for p in cli_parts)
        if not has_flag:
            cli_parts.append(f"--{flag} {default_val}")

    # Auto-detect ionmobility from workflow data type if not explicitly set.
    # FragPipe auto-detects from file extension (.d → timsTOF → ionmobility=1)
    # and writes workflow.input.data-type.im-ms=true (CmdIonquant.java:197-198).
    has_ionmobility = any("--ionmobility" in p for p in cli_parts)
    if not has_ionmobility and all_params:
        workflow_params = all_params.get("workflow", {})
        im_ms = workflow_params.get("input.data-type.im-ms", "false")
        ionmobility_val = "1" if im_ms.lower() == "true" else "0"
        cli_parts.append(f"--ionmobility {ionmobility_val}")

    return " ".join(cli_parts)


def generate_percolator_config(params: Dict[str, str]) -> str:
    """Generate Percolator CLI flags from workflow params.

    Note: Runtime parameters (num-threads, results-psms, decoy-results-psms,
    only-psms, no-terminate, post-processing-tdc) are set by the module.
    """
    cli_parts = []

    # Check if percolator should run
    if params.get("run-percolator", "false").lower() != "true":
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
    """Generate Philosopher FreeQuant CLI flags from workflow params.

    FreeQuant CLI flag names differ from the workflow file keys:
    - mz-tol -> --tol (CmdFreequant.java / QuantPanelLabelfree.java:154-161)
    - rt-tol -> --ptw (peak time window)

    Note: run-freequant is a run flag, not a CLI param.
    """
    cli_parts = []

    # Mapping from workflow file keys to Philosopher freequant CLI flags
    flag_mapping = {
        "mz-tol": "--tol",
        "rt-tol": "--ptw",
    }

    for workflow_key, cli_flag in flag_mapping.items():
        if workflow_key in params and params[workflow_key]:
            value = params[workflow_key]
            if value and value.lower() not in ("", "null", "none"):
                cli_parts.append(f"{cli_flag} {value}")

    return " ".join(cli_parts)


def generate_labelquant_config(
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]] = None,
) -> str:
    """Generate Philosopher labelquant CLI flags from workflow/TMTIntegrator params.

    Labelquant extracts TMT/isobaric reporter ion intensities from mzML spectra
    using Philosopher. It runs when tmtintegrator.intensity_extraction_tool=1
    (Philosopher mode, as opposed to IonQuant mode=0).

    CLI flags are derived from TMTIntegrator panel settings (CmdLabelquant.java:98-145):
    - --tol: mass tolerance in ppm (from tmtintegrator.tolerance)
    - --level: MS level for quantification (2=MS2, 3=MS3, from tmtintegrator.quant_level)
    - --minprob: minimum PSM probability (from tmtintegrator.min_pep_prob)
    - --purity: minimum purity (from tmtintegrator.min_purity)
    - --removelow: minimum intensity percent (from tmtintegrator.min_percent)
    - --plex: number of channels (from label type)
    - --brand: label type name lowercase (from label type)

    Note: --annot and --dir are runtime parameters set by the module.
    Note: --raw flag is auto-detected by the module (CmdLabelquant.java:147-152).
    """
    cli_parts = []

    # Get TMTIntegrator params for label info
    tmti_params = all_params.get("tmtintegrator", {}) if all_params else {}

    # Tolerance (CmdLabelquant.java:98-99): from tmtintegrator.tolerance
    tol = tmti_params.get("tolerance", params.get("tol", "20"))
    if tol:
        cli_parts.append(f"--tol {tol}")

    # Level (CmdLabelquant.java:100-108): MS level from tmtintegrator.quant_level
    quant_level = tmti_params.get("quant_level", params.get("level", "2"))
    # Handle string formats: "MS2" -> "2", "MS3" -> "3"
    isolevel_map = {"MS2": "2", "MS3": "3"}
    level = isolevel_map.get(str(quant_level), str(quant_level))
    if level:
        cli_parts.append(f"--level {level}")

    # Minprob (CmdLabelquant.java:109-110): from tmtintegrator.min_pep_prob
    minprob = tmti_params.get("min_pep_prob", params.get("minprob", "0.7"))
    if minprob:
        cli_parts.append(f"--minprob {minprob}")

    # Purity (CmdLabelquant.java:111-112): from tmtintegrator.min_purity
    purity = tmti_params.get("min_purity", params.get("purity", "0.5"))
    if purity:
        cli_parts.append(f"--purity {purity}")

    # Removelow (CmdLabelquant.java:113-114): from tmtintegrator.min_percent
    removelow = tmti_params.get("min_percent", params.get("removelow", "0.0"))
    if removelow:
        cli_parts.append(f"--removelow {removelow}")

    # Plex and brand from label type (CmdLabelquant.java:115-145)
    # Label info: label_type -> (num_channels, brand_name)
    label_info = {
        "TMT-0": (1, "tmt"),
        "TMT-2": (2, "tmt"),
        "TMT-6": (6, "tmt"),
        "TMT-10": (10, "tmt"),
        "TMT-11": (11, "tmt"),
        "TMT-16": (16, "tmt"),
        "TMT-18": (18, "tmt"),
        "TMT-35": (35, "tmt"),
        "iodoTMT-6": (6, "iodotmt"),
        "iTRAQ-4": (4, "itraq"),
        "iTRAQ-8": (8, "itraq"),
        "sCLIP-6": (6, "sclip"),
        "IBT-16": (16, "ibt"),
        "Trp-2": (2, "trp"),
    }

    # Determine label type from tmtintegrator settings
    raw_channel = tmti_params.get("channel_num", "TMT-16")
    label_type = tmti_params.get("label-type", "")

    if label_type and label_type in label_info:
        num_channels, brand = label_info[label_type]
    elif raw_channel in label_info:
        num_channels, brand = label_info[raw_channel]
    elif raw_channel.isdigit():
        # Numeric channel count - map to TMT label type
        _channel_to_label = {
            "0": "TMT-0",
            "1": "TMT-0",
            "2": "TMT-2",
            "6": "TMT-6",
            "10": "TMT-10",
            "11": "TMT-11",
            "16": "TMT-16",
            "18": "TMT-18",
            "35": "TMT-35",
        }
        mapped_label = _channel_to_label.get(raw_channel, f"TMT-{raw_channel}")
        if mapped_label in label_info:
            num_channels, brand = label_info[mapped_label]
        else:
            num_channels, brand = int(raw_channel), "tmt"
    else:
        num_channels, brand = 16, "tmt"

    cli_parts.append(f"--plex {num_channels}")
    # Brand is label type lowercase (CmdLabelquant.java:144-145: label.getType().toLowerCase())
    cli_parts.append(f"--brand {brand}")

    return " ".join(cli_parts)


def generate_msbooster_config(
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]] = None,
) -> str:
    """Generate MSBooster params file content matching FragPipe format.

    Key names match CmdMSBooster.java in FragPipe source.
    Runtime parameters (numThreads, DiaNN, unimodObo, mzmlDirectory,
    pinPepXMLDirectory) are set by the module at runtime.

    Replicates FragPipe's programmatic parameter logic:
    - FragmentationType: numeric→string translation via ui.translations (lines 171-178)
    - use-correlated-features: maps to useMultipleCorrelatedFeatures (CmdMSBooster.java)
    - useIM: requires BOTH predict-im=true AND timsTOF data (CmdMSBooster.java:170)
    - Koina params: find-best-*-model, koina-url (CmdMSBooster.java:176-184)
    """
    lines = []

    # FragmentationType numeric-to-string translation (FragPipe ui.translations:171-178).
    # Workflow files store numeric values (e.g., msbooster.fragmentation-type=0).
    # MSBoosterPanel.fragmentationType() returns the display string, which
    # CmdMSBooster writes to the params file (line 188).
    frag_type_map = {
        "0": "auto",
        "1": "HCD",
        "2": "CID",
        "3": "ETD",
        "4": "UVPD",
        "5": "ECD",
        "6": "EID",
        "7": "ETCID",
    }

    # MSBooster params file key names (from FragPipe CmdMSBooster.java)
    param_mapping = {
        "predict-rt": "useRT",
        "predict-spectra": "useSpectra",
        # Note: predict-im is handled separately (AND with timsTOF data type)
        "use-corr": "useCorr",
        "use-correlated-features": "useMultipleCorrelatedFeatures",
        "rt-model": "rtModel",
        "spectra-model": "spectraModel",
        "im-model": "imModel",
        "fragmentation-type": "FragmentationType",  # CmdMSBooster.java:188
        # M11: Koina model selection params (CmdMSBooster.java:176-184)
        "find-best-rt-model": "findBestRtModel",
        "find-best-spectra-model": "findBestSpectraModel",
        "find-best-im-model": "findBestImModel",
        "koina-url": "KoinaURL",
    }

    # Default values matching FragPipe
    defaults = {
        "useDetect": "false",
        "renamePin": "1",
        "deletePreds": "false",
    }

    # Write defaults first
    for key, value in defaults.items():
        lines.append(f"{key} = {value}")

    # Map workflow params to MSBooster params file keys
    for workflow_key, params_key in param_mapping.items():
        if workflow_key in params and params[workflow_key]:
            value = params[workflow_key]
            if value and value.lower() not in ("", "null", "none"):
                # Translate FragmentationType numeric to string
                if workflow_key == "fragmentation-type":
                    value = frag_type_map.get(value, value)
                lines.append(f"{params_key} = {value}")

    # M12: useIM requires BOTH predict-im=true AND timsTOF data (CmdMSBooster.java:170).
    # FragPipe sets useIM = hasTimsTof AND predictIm. We check workflow.input.data-type.im-ms.
    predict_im = params.get("predict-im", "true").lower() == "true"
    has_tims_tof = False
    if all_params:
        workflow_params = all_params.get("workflow", {})
        has_tims_tof = (
            workflow_params.get("input.data-type.im-ms", "false").lower() == "true"
        )
    use_im = predict_im and has_tims_tof
    # Remove any existing useIM line (from the param_mapping loop - predict-im was removed)
    # and add the correctly conditioned value
    lines = [line for line in lines if not line.startswith("useIM =")]
    lines.append(f"useIM = {str(use_im).lower()}")

    # Ensure default model names if not specified
    has_rt_model = any(line.startswith("rtModel") for line in lines)
    has_spectra_model = any(line.startswith("spectraModel") for line in lines)
    has_im_model = any(line.startswith("imModel") for line in lines)
    if not has_rt_model:
        lines.append("rtModel = DIA-NN")
    if not has_spectra_model:
        lines.append("spectraModel = DIA-NN")
    if not has_im_model:
        lines.append("imModel = DIA-NN")

    return "\n".join(lines)


def _get_peptideprophet_enzyme_flags(
    all_params: Dict[str, Dict[str, str]] = None,
) -> str:
    """Get PeptideProphet enzyme-conditional flags based on MSFragger enzyme settings.

    IMP-10: Replicates CmdPeptideProphet.java:391-404 addFreeCommandLineParams():
    - Nonspecific enzyme OR dual-enzyme (enzyme_cut_2 set): appends --nontt --nonmc
    - Custom or nocleavage enzyme: appends --nontt --nonmc --enzyme nonspecific
    """
    if not all_params:
        return ""

    msfragger = all_params.get("msfragger", {})
    enzyme_name_1 = msfragger.get(
        "search_enzyme_name_1", msfragger.get("search-enzyme-name-1", "")
    )
    _enzyme_name_2 = msfragger.get(  # noqa: F841 -- reserved for future dual-enzyme logic
        "search_enzyme_name_2", msfragger.get("search-enzyme-name-2", "")
    )
    enzyme_cut_2 = msfragger.get(
        "search_enzyme_cut_2", msfragger.get("search-enzyme-cut-2", "")
    )
    num_enzyme_termini = msfragger.get(
        "num_enzyme_termini", msfragger.get("num-enzyme-termini", "")
    )

    # Nonspecific: num_enzyme_termini=0 or enzyme_name contains "nonspecific"
    is_nonspecific = num_enzyme_termini == "0" or "nonspecific" in enzyme_name_1.lower()

    # Dual enzyme: enzyme_cut_2 is set and non-empty (and not '@' which is the
    # converted form of '-' meaning no cut)
    has_dual_enzyme = bool(enzyme_cut_2) and enzyme_cut_2 not in ("", "@")

    # Custom or nocleavage enzyme (CmdPeptideProphet.java:397-399)
    is_nocleavage = enzyme_name_1.lower() in ("nocleavage", "null")

    if is_nocleavage:
        return "--nontt --nonmc --enzyme nonspecific"
    elif is_nonspecific or has_dual_enzyme:
        return "--nontt --nonmc"

    return ""


def generate_philosopher_config(
    params: Dict[str, str],
    subcommand_prefix: str,
    all_params: Dict[str, Dict[str, str]] = None,
) -> str:
    """
    Generate Philosopher config with subcommands.

    Format: subcommand=flags

    For the phi-report filter subcommand, this function replicates FragPipe's
    programmatic behavior from CmdPhilosopherFilter.java:
    - --razor is always added when protxml is used (line 84-85: removed from user flags,
      line 114: re-added when protxml provided). Since our pipeline always uses protxml
      (ProteinProphet output), we always add --razor.
    - --picked is the default picked FDR algorithm flag. FragPipe workflows include it
      by default, but some bundled workflow files may omit it.
    - M9: When dont-use-prot-proph-file=true, --razor, --sequential, --prot are stripped
      from filter flags (CmdPhilosopherFilter.java:86-115).
    """
    lines = []

    if subcommand_prefix == "phi-report":
        # M9: Check dont-use-prot-proph-file flag (CmdPhilosopherFilter.java:86-115)
        dont_use_protxml = (
            params.get("dont-use-prot-proph-file", "false").lower() == "true"
        )

        # Filter subcommand
        if "filter" in params:
            filter_flags = params["filter"]

            if dont_use_protxml:
                # M9: When protxml is not used, strip --razor, --sequential, --prot
                # (CmdPhilosopherFilter.java:86-115). These flags require protxml.
                filter_flags = re.sub(r"\s*--razor\b", "", filter_flags)
                filter_flags = re.sub(r"\s*--sequential\b", "", filter_flags)
                filter_flags = re.sub(r"\s*--prot\s+\S+", "", filter_flags)
            else:
                # FragPipe always adds --razor when protxml is used (CmdPhilosopherFilter.java:84-114)
                # Our pipeline always uses ProteinProphet, so always ensure --razor is present
                if "--razor" not in filter_flags:
                    filter_flags = filter_flags.rstrip() + " --razor"

            # Note: --picked is NOT added by FragPipe automatically.
            # It should only be present if explicitly set in the workflow file.
            # Adding --picked for nonspecific enzyme searches can significantly
            # reduce identifications (picked FDR is more conservative).

            lines.append(f"filter={filter_flags}")

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
            cmd_opts = params.get("cmd-opts", "")
            # IMP-10: Append enzyme-conditional flags (CmdPeptideProphet.java:391-404)
            enzyme_flags = _get_peptideprophet_enzyme_flags(all_params)
            if enzyme_flags and enzyme_flags not in cmd_opts:
                cmd_opts = f"{cmd_opts} {enzyme_flags}".strip()
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
        decoy_tag = params.get("decoy-tag", "rev_")
        lines.append(f"database=--prefix {decoy_tag}")
        lines.append(f"# decoy-tag={decoy_tag}")

    return "\n".join(lines)


def _set_param(lines: list, key: str, value: str) -> None:
    """Set or override a parameter in the fragger.params lines list.

    If the parameter already exists, replace it. Otherwise, append it.
    Used by adjustMSFraggerParams logic (CmdMsfragger.java:710-741).
    """
    prefix = f"{key} ="
    for i, line in enumerate(lines):
        if line.startswith(prefix):
            lines[i] = f"{key} = {value}"
            return
    lines.append(f"{key} = {value}")


def generate_msfragger_config(
    params: Dict[str, str], all_params: Dict[str, Dict[str, str]] = None
) -> str:
    """
    Generate MSFragger parameter file (fragger.params) from workflow params.

    Format: native MSFragger params file (key = value, one per line).
    This matches how FragPipe invokes MSFragger - with a parameter file,
    NOT CLI flags.

    Handles:
    - Direct params: msfragger.param_name → param_name = value
    - Variable mods: msfragger.table.var-mods → variable_mod_NN = mass site maxocc
    - Fixed mods: msfragger.table.fix-mods → add_X_name = mass
    - Range params: msfragger.misc.fragger.* → param = lo hi (or lo,hi)
    - Boolean conversion: true/false → 1/0

    Note: Runtime parameters (database_name, num_threads) are excluded -
    they are set by the module at runtime. calibrate_mass and write_calibrated_mzml
    are included in the params file; the MSFRAGGER module only overrides them
    when task.ext values are explicitly set.
    """
    # Parameters that are set at runtime by the module
    runtime_params = {
        "database_name",
        "database-name",
        "num_threads",
        "num-threads",
    }

    # Parameters handled specially (table, misc)
    special_prefixes = {"table.", "misc."}
    special_params = {
        "delta_mass_exclude_ranges",
        "fragment_ion_series",
        "isotope_error",
    }

    # FragPipe GUI-only parameters that are NOT valid MSFragger params.
    gui_only_params = {
        "run_msfragger",
        "run-msfragger",
        "output_report_topN_dda_plus",
        "output_report_topN_dia1",
        "output_report_topN_dia2",  # IMP-9: GUI-only DIA2 topN (not a valid MSFragger param)
        "write_uncalibrated_mgf",
        "ext-thermo",  # GUI path to ThermoRawFileParser binary
        "ext_thermo",  # Same key with underscore variant
        "check_spectral_files",  # GUI validation, not needed in pipeline
        "fragpipe_ram",  # IMP-8: JVM -Xmx setting (CmdMsfragger.java:515), not a search param
        "search_enzyme_sense_1",  # MSFragger defaults to C, not needed in params file
        "search_enzyme_sense_2",  # Orphaned without enzyme_2 name/cut/nocut, causes "Unknown parameters" error
    }

    # Build params file lines: key = value
    lines = []
    for key, value in sorted(params.items()):
        if not value or value.lower() in ("null", "none", ""):
            continue
        if key in runtime_params:
            continue
        if key in gui_only_params:
            continue
        if any(key.startswith(p) for p in special_prefixes):
            continue
        if key in special_params:
            continue

        # Convert hyphens to underscores for MSFragger compatibility
        param_name = key.replace("-", "_")

        # Convert boolean strings to integers
        if value.lower() == "true":
            value = "1"
        elif value.lower() == "false":
            value = "0"

        # C3: output_format enum translation (FraggerOutputType.java).
        # Workflow files may store enum names; MSFragger expects params-file form.
        if param_name == "output_format":
            OUTPUT_FORMAT_MAP = {
                "PEPXML": "pepXML",
                "TSV": "tsv",
                "TSV_PEPXML": "tsv_pepXML",
                "PIN": "PIN",
                "TSV_PIN": "tsv_pin",
                "PEPXML_PIN": "pepXML_pin",
                "TSV_PEPXML_PIN": "tsv_pepXML_pin",
            }
            value = OUTPUT_FORMAT_MAP.get(value, value)

        # M1: Additional enum translations (TabMsfragger.java:220-290 CONVERT_TO_FILE).
        # These translate display-name/enum values to numeric forms expected by MSFragger.
        # Only applied when the value matches an enum name (numeric values pass through).
        ENUM_TRANSLATIONS = {
            "precursor_mass_units": {"PPM": "1", "Da": "0"},
            "fragment_mass_units": {"PPM": "1", "Da": "0"},
            "precursor_true_units": {"PPM": "1", "Da": "0"},
            "num_enzyme_termini": {"ENZYMATIC": "2", "SEMI": "1", "NONSPECIFIC": "0"},
        }
        if param_name in ENUM_TRANSLATIONS:
            enum_map = ENUM_TRANSLATIONS[param_name]
            if value in enum_map:
                value = enum_map[value]

        # FragPipe uses '-' for nonspecific enzyme cuts in .workflow files,
        # but MSFragger requires '@' (see FragPipe Version.java, enzyme-defs-supported.txt)
        if param_name in (
            "search_enzyme_cut_1",
            "search_enzyme_cut_2",
            "search_enzyme_nocut_1",
            "search_enzyme_nocut_2",
        ):
            value = value.replace("-", "@")

        lines.append(f"{param_name} = {value}")

    # Detect DDA+ data type from workflow.input.data-type.im-ms
    # FragPipe classifies timsTOF PASEF data as DDA+ (InputLcmsFile.java:54-101).
    is_dda_plus = False
    if all_params:
        workflow_params = all_params.get("workflow", {})
        im_ms = workflow_params.get("input.data-type.im-ms", "false")
        if im_ms.lower() == "true":
            is_dda_plus = True

    # adjustMSFraggerParams (CmdMsfragger.java:710-741): parameter overrides by data type.
    # For DDA+, FragPipe applies specific parameter adjustments that MSFragger expects.
    if is_dda_plus:
        # DDA+ overrides (CmdMsfragger.java:734-739)
        _set_param(lines, "data_type", "3")
        _set_param(lines, "isotope_error", "0")
        _set_param(lines, "intensity_transform", "1")
        _set_param(lines, "remove_precursor_peak", "1")
        _set_param(lines, "localize_delta_mass", "0")
        _set_param(lines, "labile_search_mode", "off")
        _set_param(lines, "deltamass_allowed_residues", "all")
        _set_param(lines, "report_alternative_proteins", "1")
        # M2: FragPipe always sets shifted_ions=false for DDA+ (CmdMsfragger.java:737)
        _set_param(lines, "shifted_ions", "0")
        # output_report_topN: FragPipe uses outputReportTopNDdaPlus=5 (TabMsfragger.java:1089)
        # GUI-only param output_report_topN_dda_plus overrides this if set.
        topn_dda_plus = params.get(
            "output_report_topN_dda_plus",
            params.get("output-report-topN-dda-plus", "5"),
        )
        _set_param(lines, "output_report_topN", topn_dda_plus)
    else:
        # MIN-3: Explicit data_type = 0 for DDA (CmdMsfragger.java:454).
        has_data_type = any(line.startswith("data_type =") for line in lines)
        if not has_data_type:
            lines.append("data_type = 0")

    # Decoy prefix from database config
    decoy_tag = "rev_"
    if all_params and "database" in all_params:
        decoy_tag = all_params["database"].get("decoy-tag", "rev_")
    lines.append(f"decoy_prefix = {decoy_tag}")

    # Range parameters from misc.fragger.* (FragPipe GUI properties)
    misc = {k: v for k, v in params.items() if k.startswith("misc.fragger.")}
    _add_range_param_file(
        lines,
        misc,
        "misc.fragger.precursor-charge-lo",
        "misc.fragger.precursor-charge-hi",
        "precursor_charge",
        sep=" ",
    )
    _add_range_param_file(
        lines,
        misc,
        "misc.fragger.remove-precursor-range-lo",
        "misc.fragger.remove-precursor-range-hi",
        "remove_precursor_range",
        sep=",",
    )
    _add_range_param_file(
        lines,
        misc,
        "misc.fragger.digest-mass-lo",
        "misc.fragger.digest-mass-hi",
        "digest_mass_range",
        sep=" ",
    )
    _add_range_param_file(
        lines,
        misc,
        "misc.fragger.clear-mz-lo",
        "misc.fragger.clear-mz-hi",
        "clear_mz_range",
        sep=" ",
    )

    # Special string parameters (no quoting needed in params file format)
    delta_excl = params.get("delta_mass_exclude_ranges", "")
    if delta_excl:
        lines.append(f"delta_mass_exclude_ranges = {delta_excl}")

    ion_series = params.get("fragment_ion_series", "")
    if ion_series:
        lines.append(f"fragment_ion_series = {ion_series}")

    # isotope_error: for DDA+ data, adjustMSFraggerParams forces "0" (CmdMsfragger.java:712).
    # The workflow file value (e.g. "0/1/2") must NOT override the DDA+ forced value.
    isotope_error = params.get("isotope_error", "")
    if isotope_error and not is_dda_plus:
        lines.append(f"isotope_error = {isotope_error}")

    # Variable modifications from table.var-mods
    var_mods = params.get("table.var-mods", "")
    if var_mods:
        lines.extend(_parse_variable_mods_file(var_mods))

    # Fixed modifications from table.fix-mods
    fix_mods = params.get("table.fix-mods", "")
    if fix_mods:
        lines.extend(_parse_fixed_mods_file(fix_mods))

    # L5: FragPipe always forces write_uncalibrated_mzml=1 (TabMsfragger.java:237).
    # This ensures uncalibrated mzML files are always written, which is needed
    # for downstream tools (MSBooster, etc.) that require mzML input.
    _set_param(lines, "write_uncalibrated_mzml", "1")

    return "\n".join(lines)


def _add_range_param(
    cli_parts: list,
    misc: Dict[str, str],
    lo_key: str,
    hi_key: str,
    cli_param: str,
    sep: str = " ",
) -> None:
    """Add a quoted range parameter (e.g., --precursor_charge '1 4')."""
    lo = misc.get(lo_key, "")
    hi = misc.get(hi_key, "")
    if lo and hi:
        cli_parts.append(f"--{cli_param} '{lo}{sep}{hi}'")


def _add_range_param_file(
    lines: list,
    misc: Dict[str, str],
    lo_key: str,
    hi_key: str,
    param_name: str,
    sep: str = " ",
) -> None:
    """Add a range parameter in params-file format (e.g., precursor_charge = 1 4)."""
    lo = misc.get(lo_key, "")
    hi = misc.get(hi_key, "")
    if lo and hi:
        lines.append(f"{param_name} = {lo}{sep}{hi}")


def _parse_variable_mods(var_mods_str: str) -> list:
    """
    Parse variable modifications from FragPipe table format.

    Format: "mass,site,isEnabled,maxOccurrences; ..."
    Only enabled mods with non-zero mass are included.
    Output: --variable_mod_NN 'mass site maxocc'
    """
    args = []
    mod_idx = 1
    entries = [e.strip() for e in var_mods_str.split(";") if e.strip()]

    for entry in entries:
        parts = [p.strip() for p in entry.split(",")]
        if len(parts) >= 4:
            mass, site, is_enabled, max_occ = parts[0], parts[1], parts[2], parts[3]
            if is_enabled.lower() == "true" and float(mass) != 0.0:
                idx = f"{mod_idx:02d}"
                args.append(f"--variable_mod_{idx} '{mass} {site} {max_occ}'")
                mod_idx += 1

    return args


def _parse_variable_mods_file(var_mods_str: str) -> list:
    """
    Parse variable modifications for params-file format.

    Format: "mass,site,isEnabled,maxOccurrences; ..."
    Output: variable_mod_NN = mass site maxocc
    """
    lines = []
    mod_idx = 1
    entries = [e.strip() for e in var_mods_str.split(";") if e.strip()]

    for entry in entries:
        parts = [p.strip() for p in entry.split(",")]
        if len(parts) >= 4:
            mass, site, is_enabled, max_occ = parts[0], parts[1], parts[2], parts[3]
            if is_enabled.lower() == "true" and float(mass) != 0.0:
                idx = f"{mod_idx:02d}"
                lines.append(f"variable_mod_{idx} = {mass} {site} {max_occ}")
                mod_idx += 1

    return lines


def _parse_fixed_mods(fix_mods_str: str) -> list:
    """
    Parse fixed modifications from FragPipe table format.

    Format: "mass,label,isEnabled,maxAllowed; ..."
    Maps to MSFragger --add_X_name flags.
    Only includes modifications where isEnabled is true.
    """
    args = []
    entries = [e.strip() for e in fix_mods_str.split(";") if e.strip()]

    for entry in entries:
        parts = [p.strip() for p in entry.split(",")]
        if len(parts) >= 3:
            mass = float(parts[0])
            label = parts[1]
            is_enabled = parts[2]

            if is_enabled.lower() == "true" and mass != 0.0:
                term_key = next((k for k in _TERM_MAP if label == k), None)
                if term_key:
                    args.append(f"--add_{_TERM_MAP[term_key]} {mass}")
                else:
                    aa = label[0]
                    param_suffix = _AA_MAP.get(aa)
                    if param_suffix:
                        args.append(f"--add_{param_suffix} {mass}")

    return args


def _parse_fixed_mods_file(fix_mods_str: str) -> list:
    """
    Parse fixed modifications for params-file format.

    Format: "mass,label,isEnabled,maxAllowed; ..."
    Output: add_X_name = mass (only for enabled mods)
    Only includes modifications where isEnabled is true.
    """
    lines = []
    entries = [e.strip() for e in fix_mods_str.split(";") if e.strip()]

    for entry in entries:
        parts = [p.strip() for p in entry.split(",")]
        if len(parts) >= 3:
            mass = float(parts[0])
            label = parts[1]
            is_enabled = parts[2]

            if is_enabled.lower() == "true" and mass != 0.0:
                term_key = next((k for k in _TERM_MAP if label == k), None)
                if term_key:
                    lines.append(f"add_{_TERM_MAP[term_key]} = {mass}")
                else:
                    aa = label[0]
                    param_suffix = _AA_MAP.get(aa)
                    if param_suffix:
                        lines.append(f"add_{param_suffix} = {mass}")

    return lines


# Shared amino acid and terminal modification maps
_AA_MAP = {
    "G": "G_glycine",
    "A": "A_alanine",
    "S": "S_serine",
    "P": "P_proline",
    "V": "V_valine",
    "T": "T_threonine",
    "C": "C_cysteine",
    "L": "L_leucine",
    "I": "I_isoleucine",
    "N": "N_asparagine",
    "D": "D_aspartic_acid",
    "Q": "Q_glutamine",
    "K": "K_lysine",
    "E": "E_glutamic_acid",
    "M": "M_methionine",
    "H": "H_histidine",
    "F": "F_phenylalanine",
    "R": "R_arginine",
    "Y": "Y_tyrosine",
    "W": "W_tryptophan",
    # MIN-5: User amino acids (MsfraggerParams.java:256-257, ADDON_NAMES)
    "B": "B_user_amino_acid",
    "J": "J_user_amino_acid",
    "O": "O_user_amino_acid",
    "U": "U_user_amino_acid",
    "X": "X_user_amino_acid",
    "Z": "Z_user_amino_acid",
}

_TERM_MAP = {
    "C-Term Peptide": "Cterm_peptide",
    "N-Term Peptide": "Nterm_peptide",
    "C-Term Protein": "Cterm_protein",
    "N-Term Protein": "Nterm_protein",
}


def generate_tmtintegrator_config(params: Dict[str, str]) -> str:
    """
    Generate TMTIntegrator YAML config file.

    TMTIntegrator requires a YAML file with the following structure:
    tmtintegrator:
      path: /path/to/jar
      memory: 30
      output: /path/to/output
      channel_num: 16
      ref_tag: Bridge
      ...

    Note: path, memory, and output are appended at runtime by the module.
    This function only generates the user-configurable parameters.
    """
    lines = ["tmtintegrator:"]

    # TMTIntegrator parameter mapping from workflow keys to YAML keys.
    # FragPipe .workflow files use underscores for TMTIntegrator params
    # (e.g., channel_num, ref_tag), which map directly to YAML keys.
    param_mapping = {
        "channel_num": "channel_num",
        "ref_tag": "ref_tag",
        "ref_d_tag": "ref_d_tag",
        "min_pep_prob": "min_pep_prob",
        "min_purity": "min_purity",
        "min_percent": "min_percent",
        "min_ntt": "min_ntt",
        "min_site_prob": "min_site_prob",
        "unique_gene": "unique_gene",
        "prot_exclude": "prot_exclude",
        "mod_tag": "mod_tag",
        "groupby": "groupby",
        "prot_norm": "prot_norm",
        "add_Ref": "add_Ref",
        "psm_norm": "psm_norm",
        "unique_pep": "unique_pep",
        "aggregation_method": "aggregation_method",
        "outlier_removal": "outlier_removal",
        "best_psm": "best_psm",
        "allow_overlabel": "allow_overlabel",
        "allow_unlabeled": "allow_unlabeled",
        "ms1_int": "ms1_int",
        "print_RefInt": "print_RefInt",
        "max_pep_prob_thres": "max_pep_prob_thres",
        "log2transformed": "log2transformed",
        "tolerance": "tolerance",
        "min_resolution": "min_resolution",
        "min_snr": "min_snr",
        "use_glycan_composition": "use_glycan_composition",
        "glyco_qval": "glyco_qval",
        "prefix": "prefix",
        "top3_pep": "top3_pep",
        "quant_level": "quant_level",
        # M6: TMT-35 flag (TmtiPanel.java:1217-1221)
        "is_tmt_35": "is_tmt_35",
        # M7: Abundance type (TmtiConfProps.PROP_abundance_type)
        "abn_type": "abn_type",
        # Intensity extraction tool (FragpipeRun.java:1972,2032):
        # Workflow key is "extraction_tool" (e.g., "IonQuant" or numeric 0/1/2)
        "extraction_tool": "intensity_extraction_tool",
        # Label type string (e.g., "TMT-16") - used for IonQuant isobaric pass
        "label-type": "label_type",
    }

    # Default values for parameters in TmtiConfProps.PROPS (the official TMTIntegrator
    # config properties). Only properties listed here are written to the YAML config,
    # matching FragPipe's formToConfig() behavior which skips non-PROPS keys.
    # Note: path, memory, output are appended at runtime by the module.
    defaults = {
        "channel_num": "16",
        "ref_tag": "Bridge",
        "ref_d_tag": "Pool",
        "min_pep_prob": "0.9",
        "min_purity": "0.5",
        "min_percent": "0.05",
        "min_ntt": "0",
        "min_site_prob": "-1",
        "unique_gene": "0",
        "prot_exclude": "none",
        "mod_tag": "none",
        "groupby": "-1",
        "prot_norm": "1",
        "add_Ref": "1",
        "psm_norm": "false",
        "unique_pep": "false",
        "aggregation_method": "0",
        "outlier_removal": "true",
        "best_psm": "true",
        "allow_overlabel": "true",
        "allow_unlabeled": "false",
        "ms1_int": "true",
        "print_RefInt": "false",
        "max_pep_prob_thres": "0.9",
        "log2transformed": "true",
        "min_resolution": "45000",
        "min_snr": "1000",
        "glyco_qval": "-1",
        "use_glycan_composition": "false",
    }

    # Label mass mapping: label_type -> label_mass (from QuantLabel.java:73-87)
    # TMTIntegrator needs numeric channel_num and label_masses derived from label type
    label_info = {
        "TMT-0": (1, "295.1896,224.152478"),
        "TMT-2": (2, "225.155833"),
        "TMT-6": (6, "229.162932"),
        "TMT-10": (10, "229.162932"),
        "TMT-11": (11, "229.162932"),
        "TMT-16": (16, "304.2071"),
        "TMT-18": (18, "304.2071"),
        "TMT-35": (35, "304.2071"),
        "iodoTMT-6": (6, "329.226595"),
        "iTRAQ-4": (4, "144.102063"),
        "iTRAQ-8": (8, "304.20536"),
        "sCLIP-6": (6, "481.2"),
        "IBT-16": (16, "227.1"),
        "Trp-2": (2, "875.3903"),
    }

    # Collect values from params, using defaults if not present
    yaml_params = {}
    for workflow_key, yaml_key in param_mapping.items():
        if workflow_key in params and params[workflow_key]:
            value = params[workflow_key]
            if value and value.lower() not in ("null", "none", ""):
                yaml_params[yaml_key] = value

    # Convert label type string (e.g., "TMT-16") to numeric channel_num and label_masses.
    # FragPipe stores label type in channel_num (TmtiPanel.java:1205-1209), but
    # TMTIntegrator expects numeric channel_num + separate label_masses.
    channel_value = yaml_params.get("channel_num", defaults.get("channel_num", "16"))
    if channel_value in label_info:
        num_channels, label_masses = label_info[channel_value]
        yaml_params["channel_num"] = str(num_channels)
        yaml_params["label_masses"] = label_masses
    elif not channel_value.isdigit():
        # Unknown label type - try to extract number
        import re

        match = re.search(r"(\d+)", channel_value)
        if match:
            yaml_params["channel_num"] = match.group(1)

    # Derive is_tmt_35 from label type (TmtiPanel.java:1217-1221).
    # FragPipe always writes this to the YAML even though it's not in TmtiConfProps.PROPS.
    if "is_tmt_35" not in yaml_params:
        yaml_params["is_tmt_35"] = (
            "true" if channel_value in ("TMT-35", "35") else "false"
        )

    # Output all parameters with defaults
    for yaml_key, default in defaults.items():
        value = yaml_params.get(yaml_key, default)
        lines.append(f"  {yaml_key}: {value}")

    # Append label_masses if derived (not in defaults, but needed by TMTIntegrator)
    if "label_masses" in yaml_params:
        lines.append(f"  label_masses: {yaml_params['label_masses']}")

    # Append extra params that FragPipe writes outside of the TmtiConfProps.PROPS check.
    # - is_tmt_35: explicitly added by formToConfig() (TmtiPanel.java:1217-1221)
    # - abn_type: in TmtiConfProps.PROPS as PROP_abundance_type
    #
    # NOT included (not in TmtiConfProps.PROPS, FragPipe's formToConfig skips them):
    # - tolerance, top3_pep, quant_level, intensity_extraction_tool, label_type
    # These remain in all_params["tmtintegrator"] for use by other functions
    # (e.g., _generate_tmt_ionquant_isobaric_args reads quant_level for --isolevel).
    extra_params = (
        "is_tmt_35",
        "abn_type",
    )
    for extra_key in extra_params:
        if extra_key in yaml_params and extra_key not in defaults:
            lines.append(f"  {extra_key}: {yaml_params[extra_key]}")

    return "\n".join(lines)


def generate_ptmshepherd_config(
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]] = None,
) -> str:
    """Generate PTMShepherd config file in Java properties format.

    PTMShepherd expects a `key = value` properties file (shepherd.config),
    NOT CLI flags. This replicates FragPipe's PtmshepherdParams.createConfig()
    and PtmshepherdPanel.toPtmsParamsMap() logic:

    - Format: `key = value` (one per line)
    - database and dataset lines are NOT included (set at runtime by the module)
    - run-shepherd is a run flag, not a tool parameter
    - GUI radio buttons (annotation-*, normalization-*) are converted to their
      config file equivalents (annotation_file, histo_normalizeTo)
    - iontype_* booleans are converted from true/false to 1/0
    - diagmine_mode/diagextract_mode map to run_diagmine_mode/run_diagextract_mode
    - adv_params is a GUI-only toggle, skipped

    Reference: CmdPtmshepherd.java, PtmshepherdParams.java, PtmshepherdPanel.java
    """
    lines = []

    # GUI-only params that should not appear in the config file.
    # These control the FragPipe GUI state, not PTMShepherd behavior.
    gui_only_params = {
        "run-shepherd",  # Run flag, not a tool param
        "adv_params",  # GUI toggle for advanced params panel
        "annotation-unimod",  # GUI radio button (handled below)
        "annotation-common",  # GUI radio button (handled below)
        "annotation-custom",  # GUI radio button (handled below)
        "annotation-glyco",  # GUI radio button (handled below)
        "normalization-psms",  # GUI radio button (handled below)
        "normalization-scans",  # GUI radio button (handled below)
    }

    # Convert GUI radio buttons to config file values
    # (PtmshepherdPanel.toPtmsParamsMap() lines 219-241)
    annotation_file = None
    if params.get("annotation-unimod", "").lower() == "true":
        annotation_file = "unimod"
    elif params.get("annotation-common", "").lower() == "true":
        annotation_file = "common"
    elif params.get("annotation-glyco", "").lower() == "true":
        # M15-Bug1: "glyco" is the special keyword PTMShepherd recognizes
        # (PtmshepherdPanel.java uses this, not the annotation_file param)
        annotation_file = "glyco"
    elif params.get("annotation-custom", "").lower() == "true":
        annotation_file = params.get("annotation_file", "")

    histo_normalize_to = None
    if params.get("normalization-scans", "").lower() == "true":
        histo_normalize_to = "scans"
    elif params.get("normalization-psms", "").lower() == "true":
        histo_normalize_to = "psms"

    # iontype_* conversion: true/false -> 1/0 (CONV_TO_FILE in PtmshepherdPanel)
    iontype_keys = {
        "iontype_a",
        "iontype_b",
        "iontype_c",
        "iontype_x",
        "iontype_y",
        "iontype_z",
    }

    # diagmine_mode/diagextract_mode -> run_diagmine_mode/run_diagextract_mode
    # (PtmshepherdPanel uses these as GUI checkboxes that map to config keys
    # with "run_" prefix in the actual config file)
    mode_remap = {
        "diagmine_mode": "run_diagmine_mode",
        "diagextract_mode": "run_diagextract_mode",
    }

    for key, value in sorted(params.items()):
        if not value or value.lower() in ("", "null", "none"):
            continue
        if key in gui_only_params:
            continue
        # MIN-17: Filter ui.* properties (PtmshepherdParams.java:96-97)
        if key.startswith("ui."):
            continue

        # Remap mode keys
        config_key = mode_remap.get(key, key)

        # Convert hyphens to underscores for consistency with PTMShepherd config
        config_key = config_key.replace("-", "_")

        # Convert iontype booleans: true/false -> 1/0
        if config_key in iontype_keys:
            value = "1" if value.lower() == "true" else "0"

        lines.append(f"{config_key} = {value}")

    # Add converted radio button values
    if annotation_file is not None:
        # Only add if not already present from direct annotation_file param
        has_annotation = any(line.startswith("annotation_file =") for line in lines)
        if has_annotation:
            # Replace the existing annotation_file line
            lines = [
                f"annotation_file = {annotation_file}"
                if line.startswith("annotation_file =")
                else line
                for line in lines
            ]
        else:
            lines.append(f"annotation_file = {annotation_file}")

    if histo_normalize_to is not None:
        lines.append(f"histo_normalizeTo = {histo_normalize_to}")

    # Cross-tool parameter injection from MSFragger (FragpipeRun.java:2132-2145).
    # PTMShepherd needs several MSFragger parameters that are stored under the
    # msfragger.* prefix in the workflow file, not ptmshepherd.*.
    if all_params:
        msfragger = all_params.get("msfragger", {})

        # IMP-3: Copy mass_offsets (FragpipeRun.java:2132-2135)
        mass_offsets = msfragger.get("mass_offsets", "")
        if mass_offsets:
            has_mass_offsets = any(line.startswith("mass_offsets =") for line in lines)
            if not has_mass_offsets:
                lines.append(f"mass_offsets = {mass_offsets}")

        # IMP-4: Copy mass_diff_to_variable_mod -> msfragger_massdiff_to_varmod
        # (FragpipeRun.java:2136). Note the key name translation:
        # MSFragger param: mass_diff_to_variable_mod
        # PTMShepherd config: msfragger_massdiff_to_varmod
        massdiff = msfragger.get("mass_diff_to_variable_mod", "")
        if massdiff:
            has_massdiff = any(
                line.startswith("msfragger_massdiff_to_varmod =") for line in lines
            )
            if not has_massdiff:
                lines.append(f"msfragger_massdiff_to_varmod = {massdiff}")

        # IMP-5: Copy isotope_error (FragpipeRun.java:2143-2145)
        isotope_error = msfragger.get("isotope_error", "")
        if isotope_error:
            has_isotope_error = any(
                line.startswith("isotope_error =") for line in lines
            )
            if not has_isotope_error:
                lines.append(f"isotope_error = {isotope_error}")

    # M15-Bug2: Inject glycan database file paths when glyco mode is enabled.
    # Replicates PTMSGlycanAssignPanel.getGlycanAssignParams() (lines 130-134).
    # FragPipe injects these paths from the tools directory; in Docker the path
    # is /opt/fragpipe/tools/Glycan_Databases/.
    run_glyco_mode = params.get("run_glyco_mode", "false").lower() == "true"
    if run_glyco_mode:
        glycan_db_dir = "/opt/fragpipe/tools/Glycan_Databases"
        glyco_file_params = {
            "glyco_residue_list": f"{glycan_db_dir}/glycan_residues.txt",
            "glyco_mod_list": f"{glycan_db_dir}/glycan_mods.txt",
            "glyco_oxonium_list": f"{glycan_db_dir}/oxonium_ion_list.txt",
        }
        for key, value in glyco_file_params.items():
            has_key = any(line.startswith(f"{key} =") for line in lines)
            if not has_key:
                lines.append(f"{key} = {value}")

        # M15-Bug3: Set glyco_only_mode when glyco is enabled but PTMShepherd
        # profiling is disabled. Replicates FragpipeRun.java:2139-2141:
        #   if (ptmsGlycanPanel.isRun() && !ptmshepherdPanel.isRun())
        #       additionalShepherdParams.put("glyco_only_mode", "true")
        run_shepherd = params.get("run-shepherd", "false").lower() == "true"
        if not run_shepherd:
            has_glyco_only = any(line.startswith("glyco_only_mode =") for line in lines)
            if not has_glyco_only:
                lines.append("glyco_only_mode = true")

    return "\n".join(lines)


def _get_ptmprophet_mods_from_msfragger(
    all_params: Dict[str, Dict[str, str]],
) -> str:
    """Extract MSFragger variable mods and mass offsets as PTMProphet mod string.

    Replicates PtmProphetPanel.getMSFraggerMods() (PtmProphetPanel.java:145-206):
    - Parses MSFragger table.var-mods (format: mass,sites,enabled,max_mods; ...)
    - Parses mass_offsets and merges with var mods (dedup by mass within 0.0001)
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
    msfragger = all_params.get("msfragger", {})

    # Parse variable mods into list of (mass, sites, enabled) tuples
    mods = []  # list of [mass_float, sites_str]
    var_mods_raw = msfragger.get("table.var-mods", "")
    if var_mods_raw:
        for entry in var_mods_raw.split(";"):
            entry = entry.strip()
            if not entry:
                continue
            parts = [p.strip() for p in entry.split(",")]
            if len(parts) < 4:
                continue
            try:
                mass = float(parts[0])
            except ValueError:
                print(
                    f"WARNING: Skipping unparseable mass value '{parts[0]}' in var-mods entry: {entry}",
                    file=sys.stderr,
                )
                continue
            sites = parts[1]
            is_enabled = parts[2].strip().lower() == "true"
            if is_enabled and mass != 0.0:
                mods.append([mass, sites])

    # Parse mass offsets and merge with var mods (dedup by mass within 0.0001)
    # Replicates PtmProphetPanel.java:150-175
    raw_offsets = msfragger.get("mass_offsets", "")
    if raw_offsets:
        # Mass offsets use "/" separator: "0/79.96633/15.9949"
        # Sites for offsets default to all amino acids
        offset_sites = "ACDEFGHIKLMNPQRSTVWY"
        for token in raw_offsets.replace("/", " ").split():
            try:
                offset_mass = float(token.strip())
            except ValueError:
                print(
                    f"WARNING: Skipping unparseable mass offset token '{token.strip()}' in mass_offsets: {raw_offsets}",
                    file=sys.stderr,
                )
                continue
            if offset_mass == 0.0:
                continue
            # Check if this mass already exists in var mods (within 0.0001)
            found_match = False
            for mod in mods:
                if abs(mod[0] - offset_mass) < 0.0001:
                    # Merge sites: add any offset sites not already in var mod sites
                    for ch in offset_sites:
                        if ch not in mod[1]:
                            mod[1] += ch
                    found_match = True
                    break
            if not found_match:
                mods.append([offset_mass, offset_sites])

    # Convert site notation to PTMProphet format
    # Replicates PtmProphetPanel.java:177-206
    mod_strings = []
    for mass, sites in mods:
        if sites == "" or sites == "*":
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

        if sites:
            mod_strings.append(f"{sites}:{mass}")

    return ",".join(mod_strings)


def generate_ptmprophet_config(
    params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]],
) -> str:
    """Generate PTMProphet command line config.

    Replicates PtmProphetPanel.getCmdLineOpts() (PtmProphetPanel.java:241-247):
    - If override-defaults=true AND cmdline is non-empty: return cmdline verbatim
    - Otherwise: auto-generate mods from MSFragger variable modifications and
      build the command from the default template

    The default template base tokens are:
        NOSTACK KEEPOLD STATIC EM=1 NIONS=b <mods> MINPROB=0.5

    When MSFragger uses PPM fragment mass units (fragment_mass_units=1),
    FRAGPPMTOL=<fragment_mass_tolerance> is inserted after STATIC.
    FragPipe hardcodes FRAGPPMTOL=10 (PtmProphetPanel.java:104), but we use
    the actual msfragger.fragment_mass_tolerance value for correctness.

    Reference: PtmProphetPanel.java:100-108, 113-138, 241-247
    """
    override = params.get("override-defaults", "false").lower() == "true"
    cmdline = params.get("cmdline", "").strip()

    if override and cmdline:
        return cmdline

    # Build FRAGPPMTOL token conditionally based on MSFragger fragment mass units.
    # fragment_mass_units: 0 = Daltons, 1 = PPM
    msfragger = all_params.get("msfragger", {})
    frag_units = msfragger.get("fragment_mass_units", "1")
    frag_tol = msfragger.get("fragment_mass_tolerance", "20")
    fragppmtol_token = ""
    if frag_units == "1":
        # PPM mode: include FRAGPPMTOL with the actual tolerance value
        # Strip trailing zeros for clean output (e.g., "20.0" -> "20", "10" -> "10")
        try:
            tol_val = float(frag_tol)
            fragppmtol_token = f"FRAGPPMTOL={tol_val:g}"
        except ValueError:
            print(
                f"WARNING: Cannot parse fragment_mass_tolerance '{frag_tol}' as float, using raw value",
                file=sys.stderr,
            )
            fragppmtol_token = f"FRAGPPMTOL={frag_tol}"

    # Build the template dynamically
    tokens = ["NOSTACK", "KEEPOLD", "STATIC"]
    if fragppmtol_token:
        tokens.append(fragppmtol_token)
    tokens.extend(["EM=1", "NIONS=b"])
    default_mods = "STY:79.966331,M:15.9949"
    tokens.append(default_mods)
    tokens.append("MINPROB=0.5")

    # Auto-generate: replicate loadMSFraggerMods() (PtmProphetPanel.java:113-138)
    new_mods = _get_ptmprophet_mods_from_msfragger(all_params)

    if not new_mods:
        # No mods found; return template as-is (keep defaults)
        return " ".join(tokens)

    # Find and replace the mod token in the template
    # Mod token pattern: one or more comma-separated SITES:MASS entries
    # e.g., STY:79.966331,M:15.9949 or just M:15.9949
    # Replicates PtmProphetPanel.java:117-126: Pattern.compile("[^:]+:[\\d+.-]+")
    mod_pattern = re.compile(r"[A-Za-z]+:-?[\d.]+(?:,[A-Za-z]+:-?[\d.]+)*")
    mod_index = -1
    for i, token in enumerate(tokens):
        if mod_pattern.fullmatch(token):
            mod_index = i
            break

    if mod_index == -1:
        # No previous mod argument found, append
        tokens.append(new_mods)
    else:
        tokens[mod_index] = new_mods
    return " ".join(tokens)


def generate_diaumpire_config(params: Dict[str, str]) -> str:
    """Generate DIA-Umpire params file in Java properties format.

    DIA-Umpire SE expects a properties file (umpire-se.params) with key=value pairs.
    The parameters are written by UmpireParams.save() using PropertiesUtils.
    The workflow file stores them under the "diaumpire" prefix with nested SE.* keys.

    Reference: CmdUmpireSe.java, UmpireParams.java
    """
    lines = []
    # GUI-only params that should not appear in the params file
    gui_only = {"run-diaumpire"}

    for key, value in sorted(params.items()):
        if not value or value.lower() in ("null", "none"):
            continue
        if key in gui_only:
            continue
        lines.append(f"{key} = {value}")

    return "\n".join(lines)


def generate_diann_config(params: Dict[str, str]) -> str:
    """Generate DIA-NN CLI flags from workflow params.

    DIA-NN uses standard --key value CLI flags. Some parameters are boolean
    flags (mbr), some are numeric (q-value), and cmd-opts contains additional
    raw CLI options passed verbatim.

    Reference: CmdDiann.java
    """
    cli_parts = []

    # GUI/internal-only params that should not become CLI flags
    gui_only = {
        "run-dia-nn",
        "run-dia-plex",
        "generate-msstats",  # Controls downstream MSstats generation, not DIA-NN flag
        "gene-level-report",  # FragPipe post-processing flag
        "protein-level-report",  # FragPipe post-processing flag
        "peptide-level-report",  # FragPipe post-processing flag
        "modified-peptide-level-report",  # FragPipe post-processing flag
        "site-level-report",  # FragPipe post-processing flag
    }

    # Boolean params that map to DIA-NN flags (only passed when true)
    boolean_flags = {
        "mbr",
        "unrelated-runs",
        "redo-protein-inference",
        "run-specific-protein-q-value",
    }

    # cmd-opts is passed verbatim
    cmd_opts = params.get("cmd-opts", "").strip()

    for key, value in sorted(params.items()):
        if not value or value.lower() in ("", "null", "none"):
            continue
        if key in gui_only:
            continue
        if key == "cmd-opts":
            continue  # Handled separately below

        if key in boolean_flags:
            if value.lower() == "true":
                cli_parts.append(f"--{key}")
        else:
            cli_parts.append(f"--{key} {value}")

    if cmd_opts:
        cli_parts.append(cmd_opts)

    return " ".join(cli_parts)


def generate_speclibgen_config(params: Dict[str, str]) -> str:
    """Generate SpecLibGen config as a bash-sourceable string.

    Returns a single string of KEY='value' lines that the module sources
    directly in bash. The module takes this as ONE val(config_cli) input
    and does zero Groovy parsing.

    Bash variables set:
        CONVERT_ARGS: CLI flags for EasyPQP convertpsm
        FRAGMENT_TYPES: Python list format for --fragment_types (needs quoting)
        LIBRARY_ARGS: CLI flags for EasyPQP library
        RT_CAL: RT calibration mode (noiRT, ciRT, Pierce_iRT, etc.)
        IM_CAL: IM calibration mode (noIM or file path)
        KEEP_INTERMEDIATE: Whether to keep intermediate files (true/false)

    Runtime params (set by module, NOT here):
        --threads, file paths, --decoy_prefix

    Reference: CmdSpecLibGen.java, SpeclibPanel.java
    """
    # --- Pre-compute fragment_types from easypqp.fragment.{a,b,c,x,y,z} booleans ---
    enabled = [
        ion
        for ion in ("a", "b", "c", "x", "y", "z")
        if params.get(f"easypqp.fragment.{ion}", "false") == "true"
    ]
    # EasyPQP expects Python list format: ["b","y"] (see SpeclibPanel.getEasypqp_fragment_types())
    if not enabled:
        enabled = ["b", "y"]
    fragment_types = "[" + ",".join(f'"{ion}"' for ion in enabled) + "]"

    # --- Build convert_args (CmdSpecLibGen.java:160-168) ---
    max_delta_unimod = params.get("easypqp.extras.max_delta_unimod", "0.02")
    max_delta_ppm = params.get("easypqp.extras.max_delta_ppm", "15")
    convert_parts = [
        f"--max_delta_unimod {max_delta_unimod}",
        f"--max_delta_ppm {max_delta_ppm}",
    ]
    neutral_loss = params.get("easypqp.neutral_loss", "false")
    if neutral_loss == "true":
        convert_parts.append("--enable_unspecific_losses")

    # --- Build library_args (CmdSpecLibGen.java:170) ---
    rt_lowess_fraction = params.get("easypqp.extras.rt_lowess_fraction", "0")

    # --- Map RT calibration display text to internal values (SpeclibPanel.java:459-463) ---
    rt_cal_raw = params.get("easypqp.rt-cal", "noiRT")
    rt_cal_map = {
        "Automatic selection of a run as reference RT": "noiRT",
        "noiRT": "noiRT",
        "Biognosys_iRT": "Biognosys_iRT",
        "ciRT": "ciRT",
        "Pierce_iRT": "Pierce_iRT",
    }
    rt_cal = rt_cal_map.get(rt_cal_raw, rt_cal_raw)

    # --- Map IM calibration display text to internal values ---
    im_cal_raw = params.get("easypqp.im-cal", "noIM")
    im_cal_map = {
        "Automatic selection of a run as reference IM": "noIM",
        "noIM": "noIM",
    }
    im_cal = im_cal_map.get(im_cal_raw, im_cal_raw)

    keep_intermediate = params.get("keep-intermediate-files", "false")

    # Output as bash-sourceable KEY='value' assignments separated by semicolons.
    # Single line — no newlines, no quoting headaches in Groovy/Nextflow.
    parts = [
        f"CONVERT_ARGS='{' '.join(convert_parts)}'",
        f"FRAGMENT_TYPES='{fragment_types}'",
        f"LIBRARY_ARGS='--rt_lowess_fraction {rt_lowess_fraction}'",
        f"RT_CAL='{rt_cal}'",
        f"IM_CAL='{im_cal}'",
        f"KEEP_INTERMEDIATE='{keep_intermediate}'",
    ]
    return "; ".join(parts)


def generate_opair_config(params: Dict[str, str]) -> str:
    """Generate O-Pair params in key=value format.

    O-Pair is a .NET tool that uses short CLI flags (-b, -c, -n, etc.).
    The mapping from workflow params to CLI flags is done by CmdOPair.java:
      ms2_tol -> -b (product PPM tolerance)
      ms1_tol -> -c (precursor PPM tolerance)
      max_glycans -> -n
      min_isotope_error -> -i
      max_isotope_error -> -j
      allowed_sites -> -z
      filterOxonium + oxonium file -> -f
      oxonium_minimum_intensity -> -m

    We output all params as key=value for the Nextflow module to map.

    Reference: CmdOPair.java, OPairParams.java
    """
    lines = []
    gui_only = {"run-opair"}

    for key, value in sorted(params.items()):
        if not value or value.lower() in ("null", "none"):
            continue
        if key in gui_only:
            continue
        lines.append(f"{key} = {value}")

    return "\n".join(lines)


def generate_saintexpress_config(params: Dict[str, str]) -> str:
    """Generate SAINTexpress CLI flags from workflow params.

    SAINTexpress uses short CLI flags:
      max-replicates -> -R{value}
      virtual-controls -> -L{value}
      cmd-opts -> passed verbatim

    Reference: CmdSaintExpress.java:85-89
    """
    cli_parts = []

    max_replicates = params.get("max-replicates", "")
    if max_replicates:
        cli_parts.append(f"-R{max_replicates}")

    virtual_controls = params.get("virtual-controls", "")
    if virtual_controls:
        cli_parts.append(f"-L{virtual_controls}")

    cmd_opts = params.get("cmd-opts", "").strip()
    if cmd_opts:
        cli_parts.append(cmd_opts)

    return " ".join(cli_parts)


def generate_crystalc_config(params: Dict[str, str]) -> str:
    """Generate Crystal-C params in key=value format.

    Crystal-C uses a Java properties params file (crystalc-*.params).
    The only workflow-level parameter is run-crystalc; actual params
    (raw file location, fasta, output, etc.) are set at runtime.

    Reference: CmdCrystalc.java, CrystalcParams.java
    """
    lines = []
    gui_only = {"run-crystalc"}

    for key, value in sorted(params.items()):
        if not value or value.lower() in ("null", "none"):
            continue
        if key in gui_only:
            continue
        lines.append(f"{key} = {value}")

    return "\n".join(lines)


def generate_tabrun_config(params: Dict[str, str]) -> str:
    """Generate tab-run (GUI Run tab) settings as key=value pairs.

    tab-run is NOT a standalone tool. It controls FragPipe GUI behavior:
      write_sub_mzml: write subset mzML files
      delete_temp_files: delete temp files after pipeline
      export_matched_fragments: export matched fragment ions
      sub_mzml_prob_threshold: probability threshold for subset mzML

    These are stored for reference and may be used by Nextflow modules
    that replicate CmdWriteSubMzml and CmdExportMatchedFragments behavior.

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


def generate_fpop_config(params: Dict[str, str]) -> str:
    """Generate FPOP config as a bash-sourceable string.

    Returns KEY='value' lines that the module sources in bash.
    The module takes this as ONE val(config_cli) and does zero Groovy parsing.

    Bash variables set:
        REGION_SIZE, CONTROL_LABEL, FPOP_LABEL, SUBTRACT_CONTROL, IS_TMT

    Runtime params (file paths, threads) are set by the module.
    Reference: FragPipe_FPOP_Analysis.py positional args
    """
    lines = [
        f"REGION_SIZE='{params.get('region-size', '50')}'",
        f"CONTROL_LABEL='{params.get('control-label', 'control')}'",
        f"FPOP_LABEL='{params.get('fpop-label', 'fpop')}'",
        f"SUBTRACT_CONTROL='{params.get('subtract-control', 'true')}'",
        f"IS_TMT='{params.get('is-tmt', 'false')}'",
    ]
    return "; ".join(lines)


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


def get_run_flag(all_params: Dict[str, Dict[str, str]], tool_name: str) -> bool:
    """
    Extract the run flag for a specific tool from the workflow params.

    Returns:
        True if the tool should run, False otherwise.
    """
    if tool_name not in TOOL_RUN_FLAGS:
        return False

    prefix, run_key = TOOL_RUN_FLAGS[tool_name]

    # Labelquant: runs when TMTIntegrator is enabled AND intensity_extraction_tool==1
    # (FragpipeRun.java:1935: tmtiPanel.isRun() && tmtiPanel.getIntensityExtractionTool() == 1)
    if tool_name == "labelquant":
        tmti_params = all_params.get("tmtintegrator", {})
        tmti_runs = tmti_params.get("run-tmtintegrator", "false").lower() == "true"
        intensity_tool = tmti_params.get("extraction_tool", "0")
        return tmti_runs and intensity_tool == "1"

    # Database always runs (no run flag)
    if run_key is None:
        return prefix in all_params

    # MSFragger special case: FragPipe's run-msfragger flag is a GUI checkbox state,
    # not an actual run control. In TMT Label Check workflows, FragPipe sets
    # run-msfragger=false but still runs MSFragger as the search engine.
    # If MSFragger has search parameters (table mods, misc params), it should run.
    if tool_name == "msfragger" and prefix in all_params:
        msfragger_params = all_params[prefix]
        has_search_params = any(
            k.startswith("table.") or k.startswith("misc.") for k in msfragger_params
        )
        if has_search_params:
            return True

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
    if prefix in all_params:
        value = all_params[prefix].get(run_key, "false")
        if value.lower() == "true":
            return True

    # L20: Filter/report should also run when FreeQuant is enabled
    # (FragpipeRun.java:1767: reportPanel.isRun() || quantPanelLabelfree.isRunFreeQuant())
    if tool_name == "filter":
        freequant_params = all_params.get("freequant", {})
        if freequant_params.get("run-freequant", "false").lower() == "true":
            return True

    return False


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

    if prefix not in all_params:
        return ""

    params = all_params[prefix]

    # Map tool names to their config generators
    if tool_name == "ionquant":
        return generate_ionquant_config(params, all_params)
    elif tool_name == "percolator":
        # Don't include the "disabled" comment in args
        args = generate_percolator_config(params)
        return args if not args.startswith("#") else ""
    elif tool_name == "msbooster":
        return generate_msbooster_config(params, all_params)
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
    elif tool_name == "freequant":
        # C2: Use dedicated freequant config generator with correct CLI flag names
        return generate_freequant_config(params)
    elif tool_name == "labelquant":
        # Labelquant params are derived from tmtintegrator settings, not a separate prefix
        return generate_labelquant_config(params, all_params)
    # New tool-specific generators (no silent parameter drops)
    elif tool_name == "diaumpire":
        return generate_diaumpire_config(params)
    elif tool_name == "diann":
        return generate_diann_config(params)
    elif tool_name == "speclibgen":
        # Returns dict (not string) — handled specially in generate_json_output()
        return generate_speclibgen_config(params)
    elif tool_name == "opair":
        return generate_opair_config(params)
    elif tool_name == "saintexpress":
        return generate_saintexpress_config(params)
    elif tool_name == "crystalc":
        return generate_crystalc_config(params)
    elif tool_name == "tabrun":
        return generate_tabrun_config(params)
    elif tool_name == "quantitation":
        return generate_quantitation_config(params)
    elif tool_name == "fpop":
        return generate_fpop_config(params)
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
        filter_flags = params.get("filter", "")
        # M9: Check dont-use-prot-proph-file flag
        dont_use_protxml = (
            params.get("dont-use-prot-proph-file", "false").lower() == "true"
        )
        if dont_use_protxml:
            # Strip --razor, --sequential, --prot when protxml is not used
            filter_flags = re.sub(r"\s*--razor\b", "", filter_flags)
            filter_flags = re.sub(r"\s*--sequential\b", "", filter_flags)
            filter_flags = re.sub(r"\s*--prot\s+\S+", "", filter_flags)
        else:
            # FragPipe always adds --razor when protxml is used
            if filter_flags and "--razor" not in filter_flags:
                filter_flags = filter_flags.rstrip() + " --razor"
        # Note: --picked is NOT added by FragPipe automatically.
        # Only include it if already present in workflow file.
        return filter_flags

    elif subcommand_prefix == "peptide-prophet":
        cmd_opts = params.get("cmd-opts", "")
        # IMP-10: Append enzyme-conditional flags (CmdPeptideProphet.java:391-404)
        enzyme_flags = _get_peptideprophet_enzyme_flags(all_params)
        if enzyme_flags and enzyme_flags not in cmd_opts:
            cmd_opts = f"{cmd_opts} {enzyme_flags}".strip()
        return cmd_opts

    elif subcommand_prefix == "protein-prophet":
        return params.get("cmd-opts", "")

    elif subcommand_prefix == "database":
        decoy_tag = params.get("decoy-tag", "rev_")
        return f"--prefix {decoy_tag}"

    return ""


def _override_cli_flags(base_cli: str, overrides: Dict[str, str]) -> str:
    """Override specific flags in a CLI argument string.

    For each key in overrides, replaces the existing --key value pair in base_cli,
    or appends it if not present.

    Args:
        base_cli: Base CLI flags string (e.g., "--mbr 1 --maxlfq 1")
        overrides: Dict of flag_name -> value to override (e.g., {"mbr": "0"})

    Returns:
        Updated CLI flags string with overrides applied
    """
    result = base_cli
    for flag, value in overrides.items():
        # Replace existing flag or append
        pattern = rf"--{re.escape(flag)}\s+\S+"
        replacement = f"--{flag} {value}"
        if re.search(pattern, result):
            result = re.sub(pattern, replacement, result)
        else:
            result = result.rstrip() + f" {replacement}"
    return result


def _generate_tmt_ionquant_ms1_args(
    all_params: Dict[str, Dict[str, str]],
) -> str:
    """Generate IonQuant MS1 pass args for TMT workflows.

    FragPipe hardcodes all IonQuant parameters for TMT two-pass mode
    (CmdIonquant.java:224-254). When uiCompsRepresentation is null/empty (TMT mode),
    FragPipe does NOT read from the workflow file — it uses fixed values optimized
    for maximum sensitivity in isobaric quantification.

    MS1 pass adds precursor intensity columns to psm.tsv.

    Args:
        all_params: All parsed workflow parameters keyed by tool prefix

    Returns:
        CLI flags string for IonQuant MS1 pass
    """
    ionquant_params = all_params.get("ionquant", {})

    # Auto-detect ionmobility from workflow data type (CmdIonquant.java:197-198)
    im_ms = "0"
    if all_params:
        workflow_input = all_params.get("workflow.input", all_params.get("input", {}))
        if workflow_input.get("data-type.im-ms", "false").lower() == "true":
            im_ms = "1"
        elif ionquant_params.get("ionmobility", "") == "1":
            im_ms = "1"

    # Hardcoded TMT two-pass params (CmdIonquant.java:224-254).
    # These are NOT read from the workflow file — FragPipe ignores workflow values
    # when running in TMT mode (uiCompsRepresentation is null/empty).
    cli_parts = [
        "--perform-ms1quant 1",
        "--perform-isoquant 0",
        "--isotol 20",
        "--isolevel 2",
        "--isotype tmt10",
        f"--ionmobility {im_ms}",
        "--site-reports 0",
        "--msstats 0",
        "--mbr 0",
        "--maxlfq 0",
        "--requantify 0",
        "--mztol 10",
        "--imtol 0.05",
        "--rttol 1",
        "--normalization 0",
        "--minisotopes 1",
        "--minscans 1",
        "--tp 0",
        "--minfreq 0",
        "--minions 1",
        "--locprob 0",
        "--uniqueness 0",
    ]
    return " ".join(cli_parts)


def _generate_tmt_ionquant_isobaric_args(
    tmti_params: Dict[str, str],
    all_params: Dict[str, Dict[str, str]],
) -> str:
    """Generate IonQuant isobaric pass args for TMT workflows.

    FragPipe hardcodes all IonQuant parameters for TMT two-pass mode
    (CmdIonquant.java:224-254). The only values read from the workflow/TMTIntegrator
    settings are isotype (label name), isolevel (MS2/MS3), and isotol (ppm tolerance),
    which are set earlier in CmdIonquant.java:185-196.

    Isobaric pass extracts TMT reporter ion intensities into psm.tsv.

    Args:
        tmti_params: TMTIntegrator workflow parameters
        all_params: All parsed workflow parameters keyed by tool prefix

    Returns:
        CLI flags string for IonQuant isobaric pass
    """
    ionquant_params = all_params.get("ionquant", {})

    # IMP-1: IonQuant --isotype expects the label type name (e.g., "TMT-16"), not the
    # numeric channel count. FragPipe passes label.getName() (FragpipeRun.java:2026,
    # CmdIonquant.java:195-196). The workflow file may store either "TMT-16" (label
    # type string) or "16" (numeric). Also check tmtintegrator.label-type which stores
    # the label type name directly (TmtiPanel.java).
    raw_channel = tmti_params.get("channel_num", "TMT-16")

    # Prefer label-type if available (it's the label name directly)
    label_type = tmti_params.get("label-type", "")
    if label_type:
        isotype = label_type
    elif raw_channel and not raw_channel[0].isdigit():
        # Already a label type string (e.g., "TMT-16")
        isotype = raw_channel
    else:
        # Numeric channel count -- map back to label type name.
        # IonQuant expects label names (CmdIonquant.java:195-196).
        _channel_to_label = {
            "0": "TMT-0",
            "1": "TMT-0",
            "2": "TMT-2",
            "6": "TMT-6",
            "10": "TMT-10",
            "11": "TMT-11",
            "16": "TMT-16",
            "18": "TMT-18",
            "35": "TMT-35",
        }
        isotype = _channel_to_label.get(raw_channel, f"TMT-{raw_channel}")

    # Quant level: workflow stores as integer (2=MS2, 3=MS3)
    # Also handle string formats: "MS2"->"2", "MS3"->"3"
    quant_level = str(tmti_params.get("quant_level", "2"))
    isolevel_map = {"MS2": "2", "MS3": "3", "ZOOM-HR": "4"}
    isolevel = isolevel_map.get(quant_level, quant_level)

    # Tolerance for isobaric extraction (ppm) - from TMTIntegrator settings
    # (CmdIonquant.java:185-186, reads from panel isoTol)
    tolerance = tmti_params.get("tolerance", "20")

    # Auto-detect ionmobility from workflow data type (CmdIonquant.java:197-198)
    im_ms = "0"
    if all_params:
        workflow_input = all_params.get("workflow.input", all_params.get("input", {}))
        if workflow_input.get("data-type.im-ms", "false").lower() == "true":
            im_ms = "1"
        elif ionquant_params.get("ionmobility", "") == "1":
            im_ms = "1"

    # Hardcoded TMT two-pass params (CmdIonquant.java:224-254).
    # These are NOT read from the workflow file — FragPipe ignores workflow values
    # when running in TMT mode (uiCompsRepresentation is null/empty).
    # Only isotype, isolevel, isotol are from TMTIntegrator settings.
    cli_parts = [
        "--perform-ms1quant 0",
        "--perform-isoquant 1",
        f"--isotol {tolerance}",
        f"--isolevel {isolevel}",
        f"--isotype {isotype}",
        f"--ionmobility {im_ms}",
        "--site-reports 0",
        "--msstats 0",
        "--mbr 0",
        "--maxlfq 0",
        "--requantify 0",
        "--mztol 10",
        "--imtol 0.05",
        "--rttol 1",
        "--normalization 0",
        "--minisotopes 1",
        "--minscans 1",
        "--tp 0",
        "--minfreq 0",
        "--minions 1",
        "--locprob 0",
        "--uniqueness 0",
    ]
    return " ".join(cli_parts)


def generate_json_output(
    all_params: Dict[str, Dict[str, str]],
) -> Dict[str, Dict[str, Any]]:
    """
    Generate JSON output with run flags and args for all tools.

    Returns:
        Dict mapping tool_name -> {"run": bool, "args": str}
    """
    result = {}

    for tool_name in TOOL_RUN_FLAGS:
        run_flag = get_run_flag(all_params, tool_name)
        args = generate_tool_args(all_params, tool_name) if run_flag else ""

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
            "opair",
            "crystalc",
        )
        bash_source_tools = ("speclibgen", "fpop")
        gui_only_tools = ("tabrun", "quantitation")
        if tool_name in params_file_tools:
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

        # Add modification masses for IonQuant MBR feature matching
        # (replicates FragPipe's modmasses_ionquant.txt, FragpipeRun.java:1852-1855)
        if tool_name == "ionquant" and run_flag:
            modmasses = _parse_mod_masses(all_params)
            if modmasses:
                entry["modmasses"] = modmasses

        # Add report flags for philosopher filter (CmdPhilosopherFilter.java report command)
        if tool_name == "filter" and run_flag:
            phi_report_params = all_params.get("phi-report", {})
            report_flags = []
            if phi_report_params.get("print-decoys", "false").lower() == "true":
                report_flags.append("--decoys")
            if phi_report_params.get("remove-contaminants", "false").lower() == "true":
                report_flags.append("--removecontam")
            # NOTE: --msstats is intentionally NOT added here. FragPipe gates it behind
            # `!isMultiExpReport` (CmdPhilosopherReport.java:54-59). Our pipeline always
            # runs with multiple samples (multiple LCMS file groups), so isMultiExpReport
            # is always true and FragPipe never adds --msstats to philosopher report.
            entry["report_args"] = " ".join(report_flags)

            # IMP-7: Signal multi-experiment mode for downstream --dbbin/--probin
            # injection (CmdPhilosopherFilter.java:103-113). The actual paths are
            # runtime-determined by the Nextflow module; we only emit the flag here.
            # The dont-use-prot-proph-file flag also affects filter behavior.
            entry["dont_use_protxml"] = (
                phi_report_params.get("dont-use-prot-proph-file", "false").lower()
                == "true"
            )

        # M13: Include percolator min-prob in JSON output for PercolatorOutputToPepXML.
        # CmdPercolator.java:278,298 passes minProb from percolator.min-prob.
        if tool_name == "percolator" and run_flag:
            perc_params = all_params.get("percolator", {})
            min_prob = perc_params.get("min-prob", "")
            if min_prob:
                entry["min_prob"] = min_prob

        result[tool_name] = entry

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
        modmasses = result.get("ionquant", {}).get("modmasses", "")

        # IMP-6: Only generate IonQuant two-pass entries when intensity extraction
        # tool is IonQuant (0, "IonQuant", or absent/default). Skip for Philosopher (1) or none (2).
        intensity_tool = tmti_params.get("extraction_tool", "0")
        use_ionquant_extraction = intensity_tool in ("0", "", "IonQuant")

        if use_ionquant_extraction:
            result["ionquant_ms1"] = {
                "run": True,
                "args": _generate_tmt_ionquant_ms1_args(all_params),
                "config_type": "cli",
            }
            if modmasses:
                result["ionquant_ms1"]["modmasses"] = modmasses

            result["ionquant_isobaric"] = {
                "run": True,
                "args": _generate_tmt_ionquant_isobaric_args(tmti_params, all_params),
                "config_type": "cli",
            }
            if modmasses:
                result["ionquant_isobaric"]["modmasses"] = modmasses

        # Disable standalone LFQ IonQuant when TMT two-pass is active
        if "ionquant" in result:
            result["ionquant"]["run"] = False

    # Enforce --prot 0.01 in philosopher filter when any IonQuant variant is enabled.
    # Replicates CmdPhilosopherFilter.java:86-91: when isRunIonQuant=true,
    # any existing --prot value is replaced with 0.01, or --prot 0.01 is appended.
    # This ensures consistent protein-level FDR for quantification.
    ionquant_enabled = result.get("ionquant", {}).get("run", False) or result.get(
        "ionquant_ms1", {}
    ).get("run", False)
    if ionquant_enabled and "filter" in result and result["filter"].get("run", False):
        filter_args = result["filter"].get("args", "")
        filter_args = re.sub(r"--prot\s+\S+", "--prot 0.01", filter_args)
        if "--prot" not in filter_args:
            filter_args = filter_args.rstrip() + " --prot 0.01"
        result["filter"]["args"] = filter_args

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
        content = generate_percolator_config(params)
    elif tool_prefix == "msbooster":
        content = generate_msbooster_config(params, all_params)
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
    elif tool_prefix == "diaumpire":
        content = generate_diaumpire_config(params)
        filename = "diaumpire.params"  # Java properties format
    elif tool_prefix == "diann":
        content = generate_diann_config(params)
    elif tool_prefix == "speclibgen":
        content = generate_speclibgen_config(params)
        filename = "speclibgen.params"
    elif tool_prefix == "opair":
        content = generate_opair_config(params)
        filename = "opair.params"  # Key=value properties format
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


def main():
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
        "--output-json",
        action="store_true",
        help="Output JSON with run flags and args for all tools",
    )

    args = parser.parse_args()

    # Parse workflow file
    all_params = parse_workflow_file(args.workflow)

    if args.list_tools:
        print("Tools found in workflow file:")
        for tool in sorted(all_params.keys()):
            if tool != "_global":
                print(f"  {tool}: {len(all_params[tool])} parameters")
        return 0

    # JSON output mode
    if args.output_json:
        json_output = generate_json_output(all_params)
        # Create output directory
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        # Write JSON file
        output_path = outdir / "tool_configs.json"
        with open(output_path, "w") as f:
            json.dump(json_output, f, indent=2)
        print(f"Generated: {output_path}")
        return 0

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

        filename, content = generate_tool_config(
            args.tool, all_params[args.tool], all_params
        )
        output_path = outdir / filename
        with open(output_path, "w") as f:
            f.write(content + "\n")
        print(f"Generated: {output_path}")
    else:
        # All tools
        for tool_prefix, params in all_params.items():
            if tool_prefix == "_global":
                continue

            filename, content = generate_tool_config(tool_prefix, params, all_params)
            output_path = outdir / filename
            with open(output_path, "w") as f:
                f.write(content + "\n")
            print(f"Generated: {output_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
