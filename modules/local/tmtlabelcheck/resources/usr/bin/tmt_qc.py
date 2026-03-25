#!/usr/bin/env python3
"""
TMT QC Analysis Tool

A Python CLI tool for calculating TMT labeling efficiency and channel mixing ratios
from FragPipe/MSFragger proteomics search results.

Usage:
    tmt_qc.py analyze --psm-file psm.tsv --tmt-type TMT16 --output-dir results/
    tmt_qc.py analyze --psm-file psm.tsv --abundance-file abundance.tsv --tmt-type TMT10

Author: nf-fragpipe pipeline
"""

import argparse
import sys
import re
import base64
import io
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

# matplotlib is required for plot generation.
# Import eagerly so missing dependency fails loudly at startup.
import matplotlib

matplotlib.use("Agg")  # Non-interactive backend for server use
import matplotlib.pyplot as plt

# TMT mass configurations based on FragPipe workflow specifications
TMT_MASSES: Dict[str, float] = {
    "TMT0": 224.1525,
    "TMT2": 225.1558,
    "TMT6": 229.1629,
    "TMT10": 229.1629,
    "TMT11": 229.1629,
    "TMT16": 304.2071,
    "TMT18": 304.2071,
    "TMT35": 304.2071,
    "TMTPRO": 304.2071,
}

# TMT channel definitions for different reagent types
TMT_CHANNELS: Dict[str, List[str]] = {
    "TMT6": ["126", "127", "128", "129", "130", "131"],
    "TMT10": [
        "126",
        "127N",
        "127C",
        "128N",
        "128C",
        "129N",
        "129C",
        "130N",
        "130C",
        "131",
    ],
    "TMT11": [
        "126",
        "127N",
        "127C",
        "128N",
        "128C",
        "129N",
        "129C",
        "130N",
        "130C",
        "131N",
        "131C",
    ],
    "TMT16": [
        "126",
        "127N",
        "127C",
        "128N",
        "128C",
        "129N",
        "129C",
        "130N",
        "130C",
        "131N",
        "131C",
        "132N",
        "132C",
        "133N",
        "133C",
        "134N",
    ],
    "TMT18": [
        "126",
        "127N",
        "127C",
        "128N",
        "128C",
        "129N",
        "129C",
        "130N",
        "130C",
        "131N",
        "131C",
        "132N",
        "132C",
        "133N",
        "133C",
        "134N",
        "134C",
        "135N",
    ],
    "TMTPRO": [
        "126",
        "127N",
        "127C",
        "128N",
        "128C",
        "129N",
        "129C",
        "130N",
        "130C",
        "131N",
        "131C",
        "132N",
        "132C",
        "133N",
        "133C",
        "134N",
    ],
}


def parse_annotation_file(file_path: Path) -> Dict[str, str]:
    """
    Parse TMT annotation file mapping channels to sample names.

    Expected format (tab-separated):
        126     Sample_01
        127N    Sample_02
        ...

    Args:
        file_path: Path to annotation file

    Returns:
        Dictionary mapping sample name -> channel
    """
    channel_to_sample: Dict[str, str] = {}
    sample_to_channel: Dict[str, str] = {}

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                channel = parts[0].strip()
                sample = parts[1].strip()
                if channel and sample:
                    channel_to_sample[channel] = sample
                    sample_to_channel[sample] = channel

    return sample_to_channel


def parse_tsv_file(file_path: Path) -> List[Dict[str, Any]]:
    """
    Parse a TSV file into a list of dictionaries.

    Args:
        file_path: Path to the TSV file

    Returns:
        List of dictionaries, one per row
    """
    rows = []
    with open(file_path, "r") as f:
        lines = f.readlines()
        if not lines:
            return rows

        headers = lines[0].strip().split("\t")
        for line in lines[1:]:
            if line.strip():
                values = line.strip().split("\t")
                row = {}
                for i, header in enumerate(headers):
                    row[header] = values[i] if i < len(values) else ""
                rows.append(row)
    return rows


def validate_workflow(workflow_path: Path, tmt_mass: float) -> bool:
    """
    Validate that the workflow used variable TMT modifications (not fixed).

    This is a critical safety check - labeling efficiency CANNOT be calculated
    if TMT was set as a fixed modification (100% of peptides would appear labeled).

    Args:
        workflow_path: Path to fragpipe.workflow file
        tmt_mass: Expected TMT mass value

    Returns:
        True if workflow is valid for labeling efficiency calculation

    Raises:
        SystemExit if TMT is set as fixed modification
    """
    if not workflow_path.exists():
        print(
            f"Warning: Workflow file not found at {workflow_path}, skipping validation"
        )
        return True

    content = workflow_path.read_text()
    tmt_mass_prefix = f"{tmt_mass:.4f}"[:6]  # e.g., "229.16" or "304.20"

    # Find the fixed mods line
    for line in content.split("\n"):
        if "table.fix-mods=" in line.lower() or "table.fix-mods=" in line:
            # Check if TMT mass is present on K (lysine) as fixed
            if tmt_mass_prefix in line and (
                "K (lysine)" in line or f"{tmt_mass_prefix}" in line
            ):
                # Verify it's actually set (mass > 0)
                # Pattern: 229.16293,K (lysine),true,-1
                pattern = rf"{tmt_mass_prefix}\d*,K \(lysine\),true"
                if re.search(pattern, line):
                    print("\033[91m" + "=" * 70 + "\033[0m")
                    print(
                        "\033[91mERROR: TMT is set as FIXED modification on Lysine!\033[0m"
                    )
                    print("\033[91m" + "=" * 70 + "\033[0m")
                    print()
                    print("Labeling efficiency CANNOT be calculated from this search.")
                    print(
                        "When TMT is fixed, 100% of peptides appear labeled by definition."
                    )
                    print()
                    print("To calculate labeling efficiency, re-run the search with:")
                    print(
                        f"  - TMT ({tmt_mass_prefix}...) as VARIABLE modification on K and n^"
                    )
                    print("  - Use the TMT Label Check workflow in FragPipe")
                    print()
                    sys.exit(1)

    return True


def count_tmt_modifications(assigned_mods: str, tmt_mass: float) -> int:
    """
    Count the number of TMT modifications in an assigned modifications string.

    Handles multiple formats:
    - Philosopher format: "N-term(229.1629), 5K(229.1629)"
    - Mass-only format: "229.1629, 229.1629"
    - Percolator format: "n[229], K[229]"

    Args:
        assigned_mods: The assigned modifications string from PSM file
        tmt_mass: The TMT mass to look for

    Returns:
        Number of TMT modifications found
    """
    if not assigned_mods or assigned_mods == "nan":
        return 0

    assigned_mods = str(assigned_mods)
    tmt_mass_str = f"{tmt_mass:.4f}"
    tmt_mass_prefix = str(tmt_mass)[:6]  # "229.16" or "304.20"

    count = 0

    # Count occurrences of TMT mass in various formats
    # Format 1: exact mass match (229.1629)
    count += assigned_mods.count(tmt_mass_str)

    # Format 2: shorter mass (229.16)
    if count == 0:
        count += assigned_mods.count(tmt_mass_prefix)

    # Format 3: Percolator bracket format n[229] or K[229]
    if count == 0:
        bracket_pattern = rf"\[{int(tmt_mass)}\]"
        count += len(re.findall(bracket_pattern, assigned_mods))

    return count


def has_nterm_acetylation(assigned_mods: str) -> bool:
    """
    Check if the peptide has N-terminal acetylation (biological blocking).

    N-terminal acetylation (42.0106 Da) is a biological modification that
    blocks the N-terminus from TMT labeling. This is NOT a labeling failure
    and should be excluded from efficiency calculations.

    Args:
        assigned_mods: The assigned modifications string

    Returns:
        True if N-terminal acetylation is present
    """
    if not assigned_mods or assigned_mods == "nan":
        return False

    assigned_mods = str(assigned_mods)

    # Check for acetylation patterns
    acetyl_patterns = [
        "42.0106",  # Mass value
        "42.01",  # Short mass
        "Acetyl",  # Named mod
        "acetyl",  # Lowercase
        "n[42]",  # Percolator format
        "N-term(42",  # Philosopher format
    ]

    return any(pattern in assigned_mods for pattern in acetyl_patterns)


def extract_sample_name(row: Dict[str, Any]) -> str:
    """
    Extract sample name from a PSM row.

    Tries multiple columns in order of preference:
    1. Spectrum File (Philosopher format)
    2. PSMId (Percolator format - extract from spectrum ID)
    3. Spectrum (general format)

    Args:
        row: PSM dictionary

    Returns:
        Sample name string, or 'Unknown' if not found
    """
    # Try Spectrum File column (Philosopher format)
    # Format: "interact-20250101_INST_NEO1_TMT18_User_SAMPLE001_B01s002.pep.xml"
    spectrum_file = row.get("Spectrum File", "")
    if spectrum_file:
        # Remove prefix and suffix to get sample name
        name = spectrum_file
        if name.startswith("interact-"):
            name = name[9:]  # Remove 'interact-'
        if ".pep.xml" in name:
            name = name.replace(".pep.xml", "")
        if ".mzML" in name:
            name = name.replace(".mzML", "")
        return name

    # Try Spectrum column (format: "sample.scan.scan.charge")
    spectrum = row.get("Spectrum", "")
    if spectrum:
        parts = spectrum.split(".")
        if len(parts) >= 1:
            return parts[0]

    # Try PSMId (Percolator format)
    psm_id = row.get("PSMId", "")
    if psm_id:
        # Format varies, try to extract sample portion
        parts = psm_id.split("_")
        if len(parts) >= 2:
            # Take all parts except the last few (scan numbers)
            return "_".join(parts[:-2]) if len(parts) > 2 else parts[0]

    return "Unknown"


def find_common_prefix(names: List[str]) -> str:
    """
    Find the common prefix among a list of sample names.

    Args:
        names: List of sample names

    Returns:
        Common prefix string (may be empty)
    """
    if not names or len(names) < 2:
        return ""

    # Find common prefix
    prefix = names[0]
    for name in names[1:]:
        while not name.startswith(prefix) and prefix:
            prefix = prefix[:-1]

    # Don't strip the entire name - keep at least the unique part
    # Also, try to break at a natural boundary (underscore, dash, etc.)
    if prefix:
        # Find last separator in prefix
        for sep in ["_", "-", "."]:
            last_sep = prefix.rfind(sep)
            if last_sep > 0:
                prefix = prefix[: last_sep + 1]
                break

    return prefix


def strip_common_prefix(names: List[str]) -> Tuple[List[str], str]:
    """
    Strip common prefix from sample names for cleaner display.

    Args:
        names: List of sample names

    Returns:
        Tuple of (stripped_names, common_prefix)
    """
    prefix = find_common_prefix(names)
    if prefix:
        stripped = [
            name[len(prefix) :] if name.startswith(prefix) else name for name in names
        ]
        return stripped, prefix
    return names, ""


def generate_labeling_plots(
    efficiency_data: Dict[str, Any], oxidation_data: Dict[str, Any]
) -> Optional[str]:
    """
    Generate a single compact figure with TMT labeling and oxidation plots.

    Creates a 2x2 grid similar to the R Shiny app:
    - Top left: TMT Labeling Efficiency [%] (diverging bar chart)
    - Top right: TMT Labeling Efficiency [count] (grouped bar chart)
    - Bottom left: Methionine Oxidation [%] (diverging bar chart)
    - Bottom right: Methionine Oxidation [count] (grouped bar chart)

    Args:
        efficiency_data: Dictionary from calculate_labeling_efficiency()
        oxidation_data: Dictionary from calculate_oxidation_efficiency()

    Returns:
        Base64-encoded PNG string or None if no data
    """

    # Colors: cyan (right side), yellow (partial), red/orange (left side)
    COLORS = {
        "fully_labeled": "#5AC2F1",  # Cyan/blue (right)
        "partially_labeled": "#FBCF35",  # Yellow
        "not_labeled": "#ED4C1C",  # Red/orange (left)
        "fully_oxidized": "#5AC2F1",  # Cyan/blue (right)
        "partially_oxidized": "#FBCF35",  # Yellow
        "not_oxidized": "#ED4C1C",  # Red/orange (left)
    }

    per_sample = efficiency_data.get("per_sample", [])
    ox_per_sample = oxidation_data.get("per_sample", [])

    if not per_sample:
        return None

    # Strip common prefix for cleaner labels
    sample_names = [s["sample"] for s in per_sample]
    if len(sample_names) == 1:
        # For single sample, extract short identifier from the end
        # e.g., "20250101_INST_NEO1_TMT18_User_SAMPLE001_A01s001" -> "A01s001"
        full_name = sample_names[0]
        parts = full_name.split("_")
        # Use last part as the short name
        short_name = parts[-1] if parts else full_name
        stripped_names = [short_name]
        prefix = (
            full_name[: len(full_name) - len(short_name)].rstrip("_") + "_"
            if short_name != full_name
            else ""
        )
    else:
        stripped_names, prefix = strip_common_prefix(sample_names)
    name_map = dict(zip(sample_names, stripped_names))

    # Build data for labeling plot
    labeling_data = []
    for s in per_sample:
        display_name = name_map.get(s["sample"], s["sample"])
        labeling_data.append(
            {
                "sample": display_name,
                "fully_labeled": s.get("fully_labeled_pct", 0),
                "partially_labeled": s.get("partially_labeled_pct", 0),
                "not_labeled": s.get("not_labeled_pct", 0),
                "fully_labeled_count": s.get("fully_labeled", 0),
                "partially_labeled_count": s.get("partially_labeled", 0),
                "not_labeled_count": s.get("not_labeled", 0),
            }
        )

    # Build data for oxidation plot
    ox_sample_dict = {s["sample"]: s for s in ox_per_sample}
    oxidation_data_list = []
    for s in per_sample:
        display_name = name_map.get(s["sample"], s["sample"])
        ox = ox_sample_dict.get(s["sample"], {})
        oxidation_data_list.append(
            {
                "sample": display_name,
                "not_oxidized": ox.get("not_oxidized_pct", 0),
                "partially_oxidized": ox.get("partially_oxidized_pct", 0),
                "fully_oxidized": ox.get("fully_oxidized_pct", 0),
                "not_oxidized_count": ox.get("not_oxidized", 0),
                "partially_oxidized_count": ox.get("partially_oxidized", 0),
                "fully_oxidized_count": ox.get("fully_oxidized", 0),
            }
        )

    # Generate combined 2x2 figure
    return _create_combined_figure(labeling_data, oxidation_data_list, COLORS)


def _create_combined_figure(
    labeling_data: List[Dict], oxidation_data: List[Dict], colors: Dict[str, str]
) -> Optional[str]:
    """
    Create two stacked plots (2x1) for TMT labeling and oxidation percentages.

    Layout matches R Shiny app style with diverging horizontal bars:
    - Top: TMT Labeling [%] - diverging bar with values on bars
    - Bottom: Met Oxidation [%] - diverging bar with values on bars

    Each sample is ONE horizontal row with:
    - LHS (negative direction): not labeled / not oxidized
    - RHS (positive direction): partially (yellow) + fully stacked

    All samples share the same symmetric x-axis (-100 to 100 for percentages).
    """
    if not labeling_data:
        return None

    # Calculate figure height based on number of samples
    n_samples = len(labeling_data)
    bar_height = 0.6
    row_height = 0.5  # Height per sample in inches
    plot_height = max(2.5, n_samples * row_height + 1)
    fig_height = plot_height * 2 + 0.5  # Two plots stacked

    fig, (ax_label_pct, ax_ox_pct) = plt.subplots(2, 1, figsize=(10, fig_height))

    samples = [d["sample"] for d in labeling_data]
    y_pos = list(range(len(samples)))

    # Helper function to add value labels on bars
    def add_bar_labels(
        ax, y_positions, values, x_offsets, color="white", fontsize=10, min_width=3
    ):
        """Add value labels centered on horizontal bars if they're wide enough."""
        for y, val, x_off in zip(y_positions, values, x_offsets):
            if abs(val) >= min_width:  # Only show label if bar is wide enough
                label_x = x_off + val / 2
                ax.text(
                    label_x,
                    y,
                    f"{abs(val):.1f}",
                    ha="center",
                    va="center",
                    color=color,
                    fontsize=fontsize,
                    fontweight="bold",
                )

    # ==================== TOP: TMT Labeling [%] ====================
    # Layout: Not Labeled (LEFT, red) | Partially Labeled (yellow) | Fully Labeled (cyan, RIGHT)
    not_labeled_vals = [d["not_labeled"] for d in labeling_data]
    not_labeled_pct = [-v for v in not_labeled_vals]  # Negative = left
    partial_pct = [d["partially_labeled"] for d in labeling_data]
    fully_pct = [d["fully_labeled"] for d in labeling_data]

    # Draw bars: not labeled (far left), then partially (middle), then fully (right)
    ax_label_pct.barh(
        y_pos,
        not_labeled_pct,
        height=bar_height,
        color=colors["not_labeled"],
        alpha=0.9,
        label="Not Labeled",
    )
    ax_label_pct.barh(
        y_pos,
        partial_pct,
        height=bar_height,
        color=colors["partially_labeled"],
        alpha=0.9,
        label="Partially Labeled",
    )
    ax_label_pct.barh(
        y_pos,
        fully_pct,
        left=partial_pct,
        height=bar_height,
        color=colors["fully_labeled"],
        alpha=0.9,
        label="Fully Labeled",
    )

    # Add value labels on bars (all black for visibility)
    add_bar_labels(
        ax_label_pct,
        y_pos,
        not_labeled_vals,
        not_labeled_pct,
        color="black",
        min_width=2,
    )
    add_bar_labels(
        ax_label_pct, y_pos, partial_pct, [0] * len(y_pos), color="black", min_width=2
    )
    add_bar_labels(
        ax_label_pct, y_pos, fully_pct, partial_pct, color="black", min_width=5
    )

    ax_label_pct.set_yticks(y_pos)
    # Hide y-axis labels for single sample
    if n_samples == 1:
        ax_label_pct.set_yticklabels([""])
    else:
        ax_label_pct.set_yticklabels(samples, fontsize=11)
    ax_label_pct.set_xlabel("Percent", fontsize=11)
    ax_label_pct.set_title("TMT Labeling [%]", fontweight="bold", fontsize=13)
    ax_label_pct.axvline(x=0, color="gray", linewidth=0.8)
    ax_label_pct.legend(loc="lower right", fontsize=9)
    ax_label_pct.set_xlim(-105, 105)  # Symmetric axis for all samples
    ax_label_pct.invert_yaxis()

    # ==================== BOTTOM: Met Oxidation [%] ====================
    has_ox_data = any(
        d.get("not_oxidized", 0)
        + d.get("partially_oxidized", 0)
        + d.get("fully_oxidized", 0)
        > 0
        for d in oxidation_data
    )

    if has_ox_data:
        # Layout: Not Oxidized (LEFT, cyan - good) | Partially Oxidized (yellow) | Fully Oxidized (red, RIGHT - bad)
        not_ox_vals = [d.get("not_oxidized", 0) for d in oxidation_data]
        not_ox_pct = [-v for v in not_ox_vals]  # Negative = left (good)
        partial_ox_pct = [d.get("partially_oxidized", 0) for d in oxidation_data]
        fully_ox_pct = [d.get("fully_oxidized", 0) for d in oxidation_data]

        # Draw bars: not oxidized (far left), then partially (middle), then fully (right)
        ax_ox_pct.barh(
            y_pos,
            not_ox_pct,
            height=bar_height,
            color=colors["not_oxidized"],
            alpha=0.9,
            label="Not Oxidized",
        )
        ax_ox_pct.barh(
            y_pos,
            partial_ox_pct,
            height=bar_height,
            color=colors["partially_oxidized"],
            alpha=0.9,
            label="Partially Oxidized",
        )
        ax_ox_pct.barh(
            y_pos,
            fully_ox_pct,
            left=partial_ox_pct,
            height=bar_height,
            color=colors["fully_oxidized"],
            alpha=0.9,
            label="Fully Oxidized",
        )

        # Add value labels on bars (all black for visibility)
        add_bar_labels(
            ax_ox_pct, y_pos, not_ox_vals, not_ox_pct, color="black", min_width=5
        )
        add_bar_labels(
            ax_ox_pct,
            y_pos,
            partial_ox_pct,
            [0] * len(y_pos),
            color="black",
            min_width=2,
        )
        add_bar_labels(
            ax_ox_pct, y_pos, fully_ox_pct, partial_ox_pct, color="black", min_width=2
        )

        ax_ox_pct.set_yticks(y_pos)
        # Hide y-axis labels for single sample
        if n_samples == 1:
            ax_ox_pct.set_yticklabels([""])
        else:
            ax_ox_pct.set_yticklabels(samples, fontsize=11)
        ax_ox_pct.set_xlabel("Percent", fontsize=11)
        ax_ox_pct.set_title("Met Oxidation [%]", fontweight="bold", fontsize=13)
        ax_ox_pct.axvline(x=0, color="gray", linewidth=0.8)
        ax_ox_pct.legend(loc="lower right", fontsize=9)
        ax_ox_pct.set_xlim(-105, 105)  # Symmetric axis for all samples
        ax_ox_pct.invert_yaxis()
    else:
        # No oxidation data - show placeholder
        ax_ox_pct.text(
            0.5,
            0.5,
            "No Met-containing PSMs",
            ha="center",
            va="center",
            transform=ax_ox_pct.transAxes,
            fontsize=12,
            color="gray",
        )
        ax_ox_pct.set_title("Met Oxidation [%]", fontweight="bold", fontsize=13)
        ax_ox_pct.set_xlim(-105, 105)
        ax_ox_pct.invert_yaxis()

    plt.tight_layout()

    # Convert to base64
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def calculate_oxidation_efficiency(
    psm_data: List[Dict[str, Any]],
    prob_threshold: float = 0.95,
    prob_column: str = "Probability",
) -> Dict[str, Any]:
    """
    Calculate methionine oxidation rate per sample.

    Low oxidation is desirable - high oxidation indicates sample degradation.

    Formula:
        For peptides containing M:
        - not oxidized: no M(15.9949) modifications
        - partially oxidized: some but not all M residues oxidized
        - fully oxidized: all M residues oxidized

        Non-oxidation rate = (not oxidized / total M-containing PSMs) * 100

    Args:
        psm_data: List of PSM dictionaries
        prob_threshold: Minimum probability for PSM inclusion
        prob_column: Column name for probability scores

    Returns:
        Dictionary with oxidation metrics and per-sample breakdown
    """
    # Detect available columns
    sample_row = psm_data[0] if psm_data else {}
    peptide_col = "Peptide" if "Peptide" in sample_row else "peptide"

    mods_col = None
    for col in [
        "Assigned Modifications",
        "Assigned_Modifications",
        "modifications",
        "Modifications",
    ]:
        if col in sample_row:
            mods_col = col
            break

    modified_peptide_col = None
    for col in ["Modified Peptide", "Modified_Peptide", "modified_peptide"]:
        if col in sample_row:
            modified_peptide_col = col
            break

    # Per-sample tracking
    sample_stats: Dict[str, Dict[str, int]] = {}

    for row in psm_data:
        # Filter by probability
        if prob_column in row:
            prob = float(row[prob_column])
            if prob < prob_threshold:
                continue

        peptide = row.get(peptide_col, "")
        if not peptide:
            continue

        # Count methionines in sequence
        m_count = peptide.upper().count("M")
        if m_count == 0:
            continue  # Skip peptides without methionine

        sample_name = extract_sample_name(row)

        if sample_name not in sample_stats:
            sample_stats[sample_name] = {
                "total_m_psms": 0,
                "not_oxidized": 0,
                "partially_oxidized": 0,
                "fully_oxidized": 0,
                "total_m_sites": 0,
                "oxidized_m_sites": 0,
            }

        # Count oxidized methionines
        assigned_mods = row.get(mods_col, "") if mods_col else ""
        modified_peptide = (
            row.get(modified_peptide_col, "") if modified_peptide_col else ""
        )

        # Count M(15.9949) in assigned modifications
        ox_count = 0
        if assigned_mods:
            ox_count = len(re.findall(r"M\s*\(\s*15\.99", assigned_mods))

        # Also check modified peptide for M[147] pattern
        if ox_count == 0 and modified_peptide:
            ox_count = len(re.findall(r"M\[147\]", modified_peptide))

        # Classify oxidation status
        sample_stats[sample_name]["total_m_psms"] += 1
        sample_stats[sample_name]["total_m_sites"] += m_count
        sample_stats[sample_name]["oxidized_m_sites"] += ox_count

        if ox_count == 0:
            sample_stats[sample_name]["not_oxidized"] += 1
        elif ox_count >= m_count:
            sample_stats[sample_name]["fully_oxidized"] += 1
        else:
            sample_stats[sample_name]["partially_oxidized"] += 1

    # Calculate global metrics
    total_m_psms = sum(s["total_m_psms"] for s in sample_stats.values())
    total_not_oxidized = sum(s["not_oxidized"] for s in sample_stats.values())
    total_m_sites = sum(s["total_m_sites"] for s in sample_stats.values())
    total_oxidized_sites = sum(s["oxidized_m_sites"] for s in sample_stats.values())

    # Non-oxidation rate (higher is better)
    non_oxidation_rate = (
        (total_not_oxidized / total_m_psms * 100) if total_m_psms > 0 else 100.0
    )
    # Site-level oxidation rate
    site_oxidation_rate = (
        (total_oxidized_sites / total_m_sites * 100) if total_m_sites > 0 else 0.0
    )

    # Calculate per-sample metrics
    per_sample = []
    for sample_name in sorted(sample_stats.keys()):
        stats = sample_stats[sample_name]
        total_psms = stats["total_m_psms"]
        sample_non_ox = (
            (stats["not_oxidized"] / total_psms * 100) if total_psms > 0 else 100.0
        )
        sample_site_ox = (
            (stats["oxidized_m_sites"] / stats["total_m_sites"] * 100)
            if stats["total_m_sites"] > 0
            else 0.0
        )
        per_sample.append(
            {
                "sample": sample_name,
                "total_m_psms": total_psms,
                "not_oxidized": stats["not_oxidized"],
                "partially_oxidized": stats["partially_oxidized"],
                "fully_oxidized": stats["fully_oxidized"],
                "not_oxidized_pct": (stats["not_oxidized"] / total_psms * 100)
                if total_psms > 0
                else 0.0,
                "partially_oxidized_pct": (
                    stats["partially_oxidized"] / total_psms * 100
                )
                if total_psms > 0
                else 0.0,
                "fully_oxidized_pct": (stats["fully_oxidized"] / total_psms * 100)
                if total_psms > 0
                else 0.0,
                "non_oxidation_rate": sample_non_ox,
                "site_oxidation_rate": sample_site_ox,
            }
        )

    return {
        "total_m_psms": total_m_psms,
        "total_not_oxidized": total_not_oxidized,
        "total_m_sites": total_m_sites,
        "total_oxidized_sites": total_oxidized_sites,
        "non_oxidation_rate": non_oxidation_rate,
        "site_oxidation_rate": site_oxidation_rate,
        "per_sample": per_sample,
    }


def calculate_labeling_efficiency(
    psm_data: List[Dict[str, Any]],
    tmt_mass: float,
    prob_threshold: float = 0.95,
    prob_column: str = "Probability",
) -> Dict[str, Any]:
    """
    Calculate per-peptide and per-sample TMT labeling efficiency.

    Formula:
        Total_Sites = peptide.count('K') + 1  (lysines + N-terminus)
        Labeled_Sites = count(TMT modifications)
        Efficiency = sum(Labeled_Sites) / sum(Total_Sites) * 100

    N-terminal acetylation is handled specially:
    - Acetylated N-termini are biologically blocked (not labeling failures)
    - They are excluded from the total labelable sites

    Args:
        psm_data: List of PSM dictionaries
        tmt_mass: TMT mass value
        prob_threshold: Minimum probability for PSM inclusion
        prob_column: Column name for probability scores

    Returns:
        Dictionary with efficiency metrics, per-peptide data, and per-sample breakdown
    """
    results = []
    filtered_count = 0

    # Detect available columns
    sample_row = psm_data[0] if psm_data else {}
    peptide_col = "Peptide" if "Peptide" in sample_row else "peptide"
    mods_col = None
    for col in [
        "Assigned Modifications",
        "Assigned_Modifications",
        "modifications",
        "Modifications",
    ]:
        if col in sample_row:
            mods_col = col
            break

    # Also check for Modified Peptide column (contains inline modifications)
    modified_peptide_col = None
    for col in ["Modified Peptide", "Modified_Peptide", "modified_peptide"]:
        if col in sample_row:
            modified_peptide_col = col
            break

    # Per-sample tracking
    sample_stats: Dict[str, Dict[str, int]] = {}

    for row in psm_data:
        # Filter by probability if column exists
        if prob_column in row:
            prob = float(row[prob_column])
            if prob < prob_threshold:
                filtered_count += 1
                continue

        # Get peptide sequence
        peptide = row.get(peptide_col, "")
        if not peptide:
            continue

        # Extract sample name
        sample_name = extract_sample_name(row)

        # Initialize sample stats if needed
        if sample_name not in sample_stats:
            sample_stats[sample_name] = {
                "psm_count": 0,
                "total_sites": 0,
                "labeled_sites": 0,
                "nterm_blocked": 0,
                "fully_labeled": 0,
                "partially_labeled": 0,
                "not_labeled": 0,
            }

        # Count lysines in peptide sequence
        k_count = peptide.upper().count("K")

        # Total labelable sites = lysines + 1 (N-terminus)
        total_sites = k_count + 1

        # Get modifications
        assigned_mods = row.get(mods_col, "") if mods_col else ""
        modified_peptide = (
            row.get(modified_peptide_col, "") if modified_peptide_col else ""
        )

        # Count labeled sites
        labeled_sites = count_tmt_modifications(assigned_mods, tmt_mass)

        # If no mods found in assigned mods, check modified peptide
        if labeled_sites == 0 and modified_peptide:
            labeled_sites = count_tmt_modifications(modified_peptide, tmt_mass)

        # Handle N-terminal acetylation
        nterm_blocked = has_nterm_acetylation(assigned_mods)
        if not nterm_blocked and modified_peptide:
            nterm_blocked = has_nterm_acetylation(modified_peptide)

        if nterm_blocked:
            # N-terminus was biologically blocked, not a labeling failure
            total_sites -= 1
            sample_stats[sample_name]["nterm_blocked"] += 1

        # Ensure we don't have negative sites
        total_sites = max(0, total_sites)
        labeled_sites = min(
            labeled_sites, total_sites
        )  # Can't label more than available

        # Classify PSM labeling status (like R app)
        # fully labeled: all sites labeled (N-term or acetylated + all K)
        # partially labeled: some but not all sites labeled
        # not labeled: no TMT modifications at all
        if labeled_sites == 0:
            label_status = "not_labeled"
        elif labeled_sites >= total_sites:
            label_status = "fully_labeled"
        else:
            label_status = "partially_labeled"

        # Update sample stats
        sample_stats[sample_name]["psm_count"] += 1
        sample_stats[sample_name]["total_sites"] += total_sites
        sample_stats[sample_name]["labeled_sites"] += labeled_sites
        sample_stats[sample_name][label_status] += 1

        results.append(
            {
                "peptide": peptide,
                "sample": sample_name,
                "k_count": k_count,
                "total_sites": total_sites,
                "labeled_sites": labeled_sites,
                "nterm_blocked": nterm_blocked,
                "unlabeled_sites": total_sites - labeled_sites,
                "label_status": label_status,
            }
        )

    # Aggregate global metrics
    total_labelable = sum(r["total_sites"] for r in results)
    total_labeled = sum(r["labeled_sites"] for r in results)
    total_unlabeled = sum(r["unlabeled_sites"] for r in results)
    nterm_blocked_count = sum(1 for r in results if r["nterm_blocked"])

    # Count PSMs by label status
    total_fully_labeled = sum(
        1 for r in results if r["label_status"] == "fully_labeled"
    )
    total_partially_labeled = sum(
        1 for r in results if r["label_status"] == "partially_labeled"
    )
    total_not_labeled = sum(1 for r in results if r["label_status"] == "not_labeled")

    efficiency = (total_labeled / total_labelable * 100) if total_labelable > 0 else 0.0

    # Calculate per-sample efficiency
    per_sample = []
    for sample_name in sorted(sample_stats.keys()):
        stats = sample_stats[sample_name]
        sample_efficiency = (
            (stats["labeled_sites"] / stats["total_sites"] * 100)
            if stats["total_sites"] > 0
            else 0.0
        )
        # Calculate label status percentages
        psm_total = stats["psm_count"]
        per_sample.append(
            {
                "sample": sample_name,
                "psm_count": stats["psm_count"],
                "total_sites": stats["total_sites"],
                "labeled_sites": stats["labeled_sites"],
                "unlabeled_sites": stats["total_sites"] - stats["labeled_sites"],
                "nterm_blocked": stats["nterm_blocked"],
                "efficiency": sample_efficiency,
                "fully_labeled": stats["fully_labeled"],
                "partially_labeled": stats["partially_labeled"],
                "not_labeled": stats["not_labeled"],
                "fully_labeled_pct": (stats["fully_labeled"] / psm_total * 100)
                if psm_total > 0
                else 0.0,
                "partially_labeled_pct": (stats["partially_labeled"] / psm_total * 100)
                if psm_total > 0
                else 0.0,
                "not_labeled_pct": (stats["not_labeled"] / psm_total * 100)
                if psm_total > 0
                else 0.0,
            }
        )

    total_psms = len(results)
    return {
        "total_psms": total_psms,
        "filtered_psms": filtered_count,
        "total_labelable_sites": total_labelable,
        "total_labeled_sites": total_labeled,
        "total_unlabeled_sites": total_unlabeled,
        "nterm_blocked_count": nterm_blocked_count,
        "labeling_efficiency": efficiency,
        "fully_labeled": total_fully_labeled,
        "partially_labeled": total_partially_labeled,
        "not_labeled": total_not_labeled,
        "fully_labeled_pct": (total_fully_labeled / total_psms * 100)
        if total_psms > 0
        else 0.0,
        "partially_labeled_pct": (total_partially_labeled / total_psms * 100)
        if total_psms > 0
        else 0.0,
        "not_labeled_pct": (total_not_labeled / total_psms * 100)
        if total_psms > 0
        else 0.0,
        "per_peptide": results,
        "per_sample": per_sample,
    }


def calculate_mixing_ratios(
    abundance_data: List[Dict[str, Any]], tmt_type: str
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """
    Calculate channel intensities and correction factors for TMT mixing ratios.

    Correction factors normalize channel intensities to account for unequal
    sample mixing. A correction factor > 1 means the channel was under-loaded.

    Formula:
        Correction_Factor = Mean_Intensity / Channel_Intensity

    Args:
        abundance_data: List of peptide abundance dictionaries
        tmt_type: TMT reagent type (TMT6, TMT10, TMT16, etc.)

    Returns:
        Tuple of (mixing_table, skipped_channels):
        - mixing_table: List of dictionaries with channel metrics (only channels with signal)
        - skipped_channels: List of channel names with zero signal
    """
    if not abundance_data:
        return [], []

    # Get expected channels for this TMT type
    expected_channels = TMT_CHANNELS.get(
        tmt_type.upper(), TMT_CHANNELS.get("TMT16", [])
    )

    # Find intensity columns in the data
    sample_row = abundance_data[0]
    intensity_cols = {}

    for col in sample_row.keys():
        col_lower = col.lower()
        # Match patterns like "126 Intensity", "127N Intensity", "channel_126", etc.
        for channel in expected_channels:
            if channel.lower() in col_lower and (
                "intensity" in col_lower or "abundance" in col_lower
            ):
                intensity_cols[channel] = col
                break

    if not intensity_cols:
        # Try to find any numeric columns that might be channel data
        for col in sample_row.keys():
            for channel in expected_channels:
                if channel in col:
                    intensity_cols[channel] = col
                    break

    if not intensity_cols:
        print("Warning: No TMT channel intensity columns found in abundance file")
        return [], []

    # Sum intensities for each channel
    channel_sums = {channel: 0.0 for channel in intensity_cols.keys()}

    for row in abundance_data:
        for channel, col in intensity_cols.items():
            raw_value = row.get(col, 0) or 0
            channel_sums[channel] += float(raw_value)

    # Separate channels with signal from those without (skipped channels)
    used_channels = {ch: total for ch, total in channel_sums.items() if total > 0}
    skipped_channels = [
        ch for ch in expected_channels if ch in channel_sums and channel_sums[ch] == 0
    ]

    # Also include channels that weren't found in the data at all
    missing_channels = [ch for ch in expected_channels if ch not in channel_sums]
    skipped_channels.extend(missing_channels)

    # Calculate correction factors using only channels with signal
    valid_sums = list(used_channels.values())
    mean_sum = sum(valid_sums) / len(valid_sums) if valid_sums else 0

    mixing_table = []
    for channel in expected_channels:
        if channel in used_channels:
            total_intensity = used_channels[channel]
            correction_factor = mean_sum / total_intensity if total_intensity > 0 else 0
            percent_of_mean = (total_intensity / mean_sum * 100) if mean_sum > 0 else 0

            mixing_table.append(
                {
                    "channel": channel,
                    "total_intensity": total_intensity,
                    "correction_factor": correction_factor,
                    "percent_of_mean": percent_of_mean,
                }
            )

    if skipped_channels:
        print(
            f"  Note: {len(skipped_channels)} channels with no signal (skipped): {', '.join(skipped_channels)}"
        )

    return mixing_table, skipped_channels


def determine_qc_status(
    efficiency: float, pass_threshold: float, warn_threshold: float
) -> Tuple[str, str]:
    """
    Determine QC status based on labeling efficiency.

    Args:
        efficiency: Labeling efficiency percentage
        pass_threshold: Threshold for PASS status
        warn_threshold: Threshold for WARNING status

    Returns:
        Tuple of (status, css_class)
    """
    if efficiency >= pass_threshold * 100:
        return "PASS", "pass"
    elif efficiency >= warn_threshold * 100:
        return "WARNING", "warning"
    else:
        return "FAIL", "fail"


def generate_html_report(
    efficiency_data: Dict[str, Any],
    oxidation_data: Dict[str, Any],
    mixing_table: List[Dict[str, Any]],
    skipped_channels: List[str],
    tmt_type: str,
    pass_threshold: float,
    warn_threshold: float,
    output_path: Path,
) -> None:
    """Generate comprehensive HTML report with methodology and recommendations."""

    status, status_class = determine_qc_status(
        efficiency_data["labeling_efficiency"], pass_threshold, warn_threshold
    )

    status_emoji = {"PASS": "✅", "WARNING": "⚠️", "FAIL": "❌"}.get(status, "❓")

    # Get TMT mass for this type
    tmt_mass = TMT_MASSES.get(tmt_type.upper(), 229.1629)

    # Extract experiment ID from common prefix of sample names
    # For single sample, use the full sample name as experiment ID
    per_sample_data = efficiency_data.get("per_sample", [])
    sample_names = [s["sample"] for s in per_sample_data]
    if len(sample_names) == 1:
        experiment_id = sample_names[0]
    else:
        _, experiment_id = strip_common_prefix(sample_names)
        # Clean up experiment ID (remove trailing underscore if present)
        experiment_id = experiment_id.rstrip("_") if experiment_id else None

    # Generate plots if matplotlib is available
    combined_plot_b64 = generate_labeling_plots(efficiency_data, oxidation_data)
    oxidation_per_sample = {
        s["sample"]: s for s in oxidation_data.get("per_sample", [])
    }

    problematic_samples = []
    for sample in per_sample_data:
        sample_status, _ = determine_qc_status(
            sample["efficiency"], pass_threshold, warn_threshold
        )
        ox_data = oxidation_per_sample.get(sample["sample"], {})
        if sample_status != "PASS":
            problematic_samples.append(
                {
                    "sample": sample["sample"],
                    "labeling_status": sample_status,
                    "labeling_efficiency": sample["efficiency"],
                    "oxidation_rate": ox_data.get("site_oxidation_rate", 0),
                }
            )

    # Generate recommendation based on status
    if status == "PASS":
        recommendation = """
        <div class="recommendation pass-bg">
            <h3>✅ Recommendation: Proceed with Analysis</h3>
            <p>The labeling efficiency meets the quality threshold. Your samples are ready for full TMT quantitative analysis.</p>
            <ul>
                <li>Proceed with the complete FragPipe TMT workflow</li>
                <li>Use TMT as a <strong>fixed modification</strong> in the main analysis</li>
                <li>Expected quantification accuracy: High</li>
            </ul>
        </div>
        """
    elif status == "WARNING":
        recommendation = f"""
        <div class="recommendation warning-bg">
            <h3>⚠️ Recommendation: Review Before Proceeding</h3>
            <p>The labeling efficiency is below optimal ({efficiency_data["labeling_efficiency"]:.1f}% vs {pass_threshold * 100:.0f}% threshold).
            This may affect quantification accuracy.</p>
            <ul>
                <li><strong>Option 1:</strong> Proceed with caution - results may have higher variability</li>
                <li><strong>Option 2:</strong> Consider re-labeling samples if critical experiment</li>
                <li><strong>Option 3:</strong> Check reagent quality and labeling protocol</li>
            </ul>
            <p><strong>Potential causes:</strong> Old/degraded TMT reagent, insufficient reagent, suboptimal pH, competing amines in sample</p>
        </div>
        """
    else:
        recommendation = f"""
        <div class="recommendation fail-bg">
            <h3>❌ Recommendation: Do Not Proceed</h3>
            <p>The labeling efficiency ({efficiency_data["labeling_efficiency"]:.1f}%) is significantly below the acceptable threshold ({warn_threshold * 100:.0f}%).
            Proceeding with analysis will produce unreliable quantification results.</p>
            <ul>
                <li><strong>Action Required:</strong> Re-label samples before proceeding</li>
                <li>Check TMT reagent expiration date and storage conditions</li>
                <li>Verify labeling protocol: pH 8.5, sufficient reagent (0.8mg per 100μg protein)</li>
                <li>Ensure sample is free of primary amines (Tris, ammonium salts)</li>
            </ul>
            <p><strong>Common issues:</strong> Expired reagent, incorrect pH, amine contamination, insufficient quenching time</p>
        </div>
        """

    # Generate plots HTML if matplotlib is available
    plots_html = ""
    if combined_plot_b64:
        plots_html = f"""
        <h2 id="plots">📊 Visual Summary</h2>
        <p>Bar charts showing labeling and oxidation status across samples.
        In the percentage plots, values extend left for undesirable outcomes (not labeled/oxidized)
        and right for desirable outcomes (fully labeled/not oxidized).</p>
        <div class="plot-container">
            <img src="data:image/png;base64,{combined_plot_b64}" alt="TMT Labeling and Oxidation Summary" style="max-width: 100%; height: auto;">
        </div>
        """

    # Generate mixing table HTML if available
    mixing_html = ""
    if mixing_table:
        mixing_rows = "\n".join(
            [
                f"<tr><td>{row['channel']}</td><td>{row['total_intensity']:.2e}</td>"
                f"<td>{row['correction_factor']:.3f}</td><td>{row['percent_of_mean']:.1f}%</td></tr>"
                for row in mixing_table
            ]
        )

        # Add skipped channels section if any
        skipped_html = ""
        if skipped_channels:
            skipped_html = f"""
        <div class="info-card" style="margin-top: 15px; background: #fef9e7; border-left: 4px solid #f39c12;">
            <div class="label" style="color: #9a7b4f; font-weight: bold;">⚠️ Skipped Channels (No Signal)</div>
            <div class="value" style="font-size: 16px; color: #7d6608;">{", ".join(skipped_channels)}</div>
            <p style="font-size: 12px; color: #9a7b4f; margin-top: 8px;">
                These {len(skipped_channels)} channel(s) had zero intensity and are excluded from mixing ratio calculations.
                This is normal when using fewer samples than the TMT reagent capacity.
            </p>
        </div>
            """

        mixing_html = f"""
        <h2>📊 Channel Mixing Ratios</h2>
        <p>The table below shows the signal distribution across TMT channels. Correction factors normalize
        for unequal sample loading. Values close to 1.0 indicate balanced mixing.</p>
        <table>
            <tr><th>Channel</th><th>Total Intensity</th><th>Correction Factor</th><th>% of Mean</th></tr>
            {mixing_rows}
        </table>
        {skipped_html}
        """

    # Generate per-sample efficiency tables (separate for labeling and oxidation)
    per_sample_html = ""
    per_sample_data = efficiency_data.get("per_sample", [])
    if per_sample_data and len(per_sample_data) >= 1:
        # Generate rows with color-coded efficiency
        def get_efficiency_class(eff: float) -> str:
            if eff >= pass_threshold * 100:
                return "pass"
            elif eff >= warn_threshold * 100:
                return "warning"
            return "fail"

        # Strip common prefix from sample names
        sample_names = [row["sample"] for row in per_sample_data]
        stripped_names, common_prefix = strip_common_prefix(sample_names)
        name_map = dict(zip(sample_names, stripped_names))

        # Build labeling table rows
        labeling_rows_list = []
        for row in per_sample_data:
            display_name = name_map.get(row["sample"], row["sample"])
            labeling_rows_list.append(
                f"<tr><td>{display_name}</td><td>{row['psm_count']:,}</td>"
                f'<td class="{get_efficiency_class(row["efficiency"])}">{row["efficiency"]:.1f}%</td>'
                f'<td class="pass">{row.get("fully_labeled_pct", 0):.1f}%</td>'
                f'<td class="warning">{row.get("partially_labeled_pct", 0):.1f}%</td>'
                f'<td class="fail">{row.get("not_labeled_pct", 0):.1f}%</td></tr>'
            )
        labeling_rows = "\n".join(labeling_rows_list)

        # Build oxidation table rows
        oxidation_rows_list = []
        for row in per_sample_data:
            display_name = name_map.get(row["sample"], row["sample"])
            ox_data = oxidation_per_sample.get(row["sample"], {})
            ox_rate = ox_data.get("site_oxidation_rate", 0)
            not_ox_pct = ox_data.get("not_oxidized_pct", 0)
            partial_ox_pct = ox_data.get("partially_oxidized_pct", 0)
            fully_ox_pct = ox_data.get("fully_oxidized_pct", 0)
            m_psms = ox_data.get("total_m_psms", 0)

            oxidation_rows_list.append(
                f"<tr><td>{display_name}</td><td>{m_psms:,}</td>"
                f"<td>{ox_rate:.1f}%</td>"
                f'<td class="fail">{fully_ox_pct:.1f}%</td>'
                f'<td class="warning">{partial_ox_pct:.1f}%</td>'
                f'<td class="pass">{not_ox_pct:.1f}%</td></tr>'
            )
        oxidation_rows = "\n".join(oxidation_rows_list)

        # Add prefix note if stripped
        prefix_note = ""
        if common_prefix:
            prefix_note = f'<p style="font-size: 12px; color: #7f8c8d;">Common prefix removed: <code>{common_prefix}</code></p>'

        per_sample_html = f"""
        <h2>🧫 Per-Sample Quality</h2>
        {prefix_note}

        <h3>TMT Labeling</h3>
        <p>PSM breakdown by labeling status: fully labeled (all sites), partially labeled (some sites), or not labeled (no TMT).</p>
        <table class="centered-table">
            <tr>
                <th>Sample</th>
                <th>PSMs</th>
                <th>Site Eff.</th>
                <th>Fully</th>
                <th>Partial</th>
                <th>None</th>
            </tr>
            {labeling_rows}
        </table>

        <h3>Methionine Oxidation</h3>
        <p>Oxidation status for Met-containing PSMs. Lower oxidation indicates better sample quality.</p>
        <table class="centered-table">
            <tr>
                <th>Sample</th>
                <th>Met PSMs</th>
                <th>Site Ox</th>
                <th>Fully</th>
                <th>Partial</th>
                <th>Not Ox</th>
            </tr>
            {oxidation_rows}
        </table>
        """

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>TMT Labeling QC Report - {tmt_type}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            margin: 0;
            padding: 40px;
            background: #f5f7fa;
            color: #2c3e50;
            line-height: 1.6;
        }}
        .container {{
            max-width: 1000px;
            margin: 0 auto;
            background: white;
            padding: 50px;
            border-radius: 12px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 15px;
            margin-top: 0;
        }}
        h2 {{
            color: #34495e;
            margin-top: 40px;
            border-bottom: 1px solid #ecf0f1;
            padding-bottom: 10px;
        }}
        h3 {{
            color: #2c3e50;
            margin-top: 0;
            margin-bottom: 10px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 12px 15px;
            text-align: center;
        }}
        th {{
            background-color: #3498db;
            color: white;
            font-weight: 600;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .pass {{ color: #27ae60; }}
        .warning {{ color: #f39c12; }}
        .fail {{ color: #e74c3c; }}
        .metric-box {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 30px;
            border-radius: 12px;
            margin: 25px 0;
            color: white;
            text-align: center;
        }}
        .big-number {{
            font-size: 72px;
            font-weight: bold;
            margin: 10px 0;
        }}
        .status-badge {{
            display: flex;
            align-items: center;
            justify-content: center;
            width: 60px;
            height: 60px;
            border-radius: 30px;
            font-size: 32px;
            margin: 15px auto 0 auto;
        }}
        .status-badge.pass {{ background: #27ae60; }}
        .status-badge.warning {{ background: #f39c12; }}
        .status-badge.fail {{ background: #e74c3c; }}
        .status-text {{
            font-weight: bold;
            font-size: 18px;
            margin-top: 8px;
        }}
        .info-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 25px 0;
        }}
        .info-card {{
            background: #ecf0f1;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .info-card .value {{
            font-size: 28px;
            font-weight: bold;
            color: #2c3e50;
        }}
        .info-card .label {{
            font-size: 14px;
            color: #7f8c8d;
            margin-top: 5px;
        }}
        .recommendation {{
            padding: 25px;
            border-radius: 8px;
            margin: 25px 0;
        }}
        .recommendation ul {{
            margin: 15px 0;
            padding-left: 20px;
        }}
        .recommendation li {{
            margin: 8px 0;
        }}
        .pass-bg {{
            background: #d5f4e6;
            border-left: 5px solid #27ae60;
        }}
        .warning-bg {{
            background: #fef9e7;
            border-left: 5px solid #f39c12;
        }}
        .fail-bg {{
            background: #fdedec;
            border-left: 5px solid #e74c3c;
        }}
        .method-box {{
            background: #eaf2f8;
            padding: 25px;
            border-radius: 8px;
            margin: 20px 0;
            border-left: 5px solid #3498db;
        }}
        .method-box h3 {{
            color: #2980b9;
            margin-top: 0;
        }}
        code {{
            background: #f4f4f4;
            padding: 2px 6px;
            border-radius: 4px;
            font-family: 'Monaco', 'Menlo', monospace;
            font-size: 0.9em;
        }}
        .formula {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            font-family: 'Monaco', 'Menlo', monospace;
            margin: 15px 0;
            text-align: center;
            font-size: 1.1em;
        }}
        .plot-container {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
            border: 1px solid #e0e0e0;
            text-align: center;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}
        .plot-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 4px;
        }}
        footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ecf0f1;
            color: #7f8c8d;
            font-size: 14px;
            text-align: center;
        }}
        .centered-table th, .centered-table td {{
            text-align: center;
        }}
        .metric-cards {{
            display: flex;
            gap: 20px;
            margin: 25px 0;
            flex-wrap: wrap;
        }}
        .metric-card {{
            flex: 1;
            min-width: 200px;
            padding: 25px;
            border-radius: 12px;
            text-align: center;
        }}
        .metric-card.labeling {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
        }}
        .metric-card.oxidation {{
            background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
            color: white;
        }}
        .metric-card .metric-value {{
            font-size: 48px;
            font-weight: bold;
            margin: 10px 0;
        }}
        .metric-card .metric-label {{
            font-size: 14px;
            opacity: 0.9;
        }}
        .metric-card .metric-status {{
            font-size: 24px;
            margin-top: 10px;
        }}
        .alert-box {{
            padding: 15px 20px;
            border-radius: 8px;
            margin: 20px 0;
        }}
        .alert-box.warning {{
            background: #fef9e7;
            border-left: 5px solid #f39c12;
        }}
        .alert-box.danger {{
            background: #fdedec;
            border-left: 5px solid #e74c3c;
        }}
        .alert-box h4 {{
            margin: 0 0 10px 0;
        }}
        .alert-box ul {{
            margin: 0;
            padding-left: 20px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧪 TMT Labeling Efficiency QC Report</h1>
        <p>{
        f"<strong>Experiment:</strong> {experiment_id} &nbsp;|&nbsp; "
        if experiment_id
        else ""
    }<strong>TMT Reagent:</strong> {tmt_type} ({tmt_mass:.4f} Da) &nbsp;|&nbsp;
           <strong>Generated:</strong> {datetime.now().strftime("%Y-%m-%d")}</p>

        <h2 id="results">📈 Results Summary</h2>

        {
        "".join(
            [
                f'''
        <div class="alert-box {'danger' if any(s['labeling_status'] == 'FAIL' for s in [p for p in problematic_samples if p['labeling_status'] != 'PASS']) else 'warning'}">
            <h4>⚠️ Samples with Labeling Efficiency Issues</h4>
            <ul>
            {"".join([f"<li><strong>{s['sample'].split('_')[-1] if '_' in s['sample'] else s['sample']}</strong>: {s['labeling_efficiency']:.1f}% ({s['labeling_status']})</li>" for s in problematic_samples if s['labeling_status'] != 'PASS'][:5])}
            {f"<li>...and {len([s for s in problematic_samples if s['labeling_status'] != 'PASS']) - 5} more</li>" if len([s for s in problematic_samples if s['labeling_status'] != 'PASS']) > 5 else ""}
            </ul>
        </div>
        '''
            ]
        )
        if any(s["labeling_status"] != "PASS" for s in problematic_samples)
        else ""
    }

        <div class="metric-cards">
            <div class="metric-card labeling">
                <div class="metric-label">Site Labeling Efficiency</div>
                <div class="metric-value">{
        efficiency_data["labeling_efficiency"]:.1f}%</div>
                <div class="metric-status">{status_emoji} {status}</div>
            </div>
            <div class="metric-card oxidation">
                <div class="metric-label">Met Oxidation Rate</div>
                <div class="metric-value">{
        oxidation_data["site_oxidation_rate"]:.1f}%</div>
            </div>
        </div>

        <div class="info-grid">
            <div class="info-card">
                <div class="value">{len(per_sample_data)}</div>
                <div class="label">Samples</div>
            </div>
            <div class="info-card">
                <div class="value pass">{
        efficiency_data.get("fully_labeled_pct", 0):.1f}%</div>
                <div class="label">Fully Labeled PSMs</div>
            </div>
            <div class="info-card">
                <div class="value warning">{
        efficiency_data.get("partially_labeled_pct", 0):.1f}%</div>
                <div class="label">Partially Labeled PSMs</div>
            </div>
            <div class="info-card">
                <div class="value fail">{
        efficiency_data.get("not_labeled_pct", 0):.1f}%</div>
                <div class="label">Not Labeled PSMs</div>
            </div>
        </div>

        <h2 id="recommendation">🎯 Recommendation</h2>
        {recommendation}

        {plots_html}

        <h2 id="statistics">📊 Detailed Statistics</h2>
        <table>
            <tr><th>Metric</th><th>Value</th><th>Description</th></tr>
            <tr><td><strong>Site Labeling Efficiency</strong></td><td><strong class="{
        status_class
    }">{
        efficiency_data[
            "labeling_efficiency"
        ]:.2f}%</strong></td><td>Primary QC metric (labeled sites / total sites)</td></tr>
            <tr><td>Total PSMs Analyzed</td><td>{
        efficiency_data[
            "total_psms"
        ]:,}</td><td>Peptide-spectrum matches passing probability threshold (&gt;0.95)</td></tr>
            <tr><td>PSMs Filtered Out</td><td>{
        efficiency_data[
            "filtered_psms"
        ]:,}</td><td>Low-confidence identifications excluded from analysis</td></tr>
            <tr><td>Total Labelable Sites</td><td>{
        efficiency_data[
            "total_labelable_sites"
        ]:,}</td><td>Sum of (K residues + N-termini) across all PSMs</td></tr>
            <tr><td>Sites Labeled (TMT)</td><td>{
        efficiency_data[
            "total_labeled_sites"
        ]:,}</td><td>Sites with detected TMT modification</td></tr>
            <tr><td>Sites Unlabeled</td><td>{
        efficiency_data[
            "total_unlabeled_sites"
        ]:,}</td><td>Sites without TMT modification (labeling failures)</td></tr>
            <tr><td>N-term Acetylated</td><td>{
        efficiency_data[
            "nterm_blocked_count"
        ]:,}</td><td>Biologically blocked N-termini (not counted as failures)</td></tr>
        </table>

        {mixing_html}

        {per_sample_html}

        <h2 id="thresholds">⚖️ QC Thresholds</h2>
        <h3>TMT Labeling Efficiency</h3>
        <table>
            <tr><th>Status</th><th>Labeling Threshold</th><th>Interpretation</th></tr>
            <tr>
                <td class="pass">✅ PASS</td>
                <td>≥ {pass_threshold * 100:.0f}%</td>
                <td>Excellent labeling - proceed with full analysis</td>
            </tr>
            <tr>
                <td class="warning">⚠️ WARNING</td>
                <td>{warn_threshold * 100:.0f}% – {pass_threshold * 100:.0f}%</td>
                <td>Suboptimal labeling - review before proceeding</td>
            </tr>
            <tr>
                <td class="fail">❌ FAIL</td>
                <td>&lt; {warn_threshold * 100:.0f}%</td>
                <td>Poor labeling - re-label samples recommended</td>
            </tr>
        </table>

        <h2 id="methodology">🔬 Methodology</h2>

        <div class="method-box">
            <h3>What is TMT Labeling Efficiency?</h3>
            <p>TMT (Tandem Mass Tag) labeling efficiency measures the percentage of available labeling sites
            (lysine residues and peptide N-termini) that successfully received a TMT tag. High efficiency
            is critical for accurate multiplexed quantification.</p>

            <h3>How is it Calculated?</h3>
            <p>For each identified peptide, we count:</p>
            <ul>
                <li><strong>Labelable sites</strong> = Number of lysine (K) residues + 1 (N-terminus)</li>
                <li><strong>Labeled sites</strong> = Number of detected TMT modifications ({
        tmt_mass:.4f} Da)</li>
            </ul>

            <div class="formula">
                Labeling Efficiency = (Total Labeled Sites / Total Labelable Sites) × 100%
            </div>

            <h3>Special Handling: N-terminal Acetylation</h3>
            <p>Peptides with N-terminal acetylation (42.01 Da) are handled specially. This is a natural
            biological modification that blocks the N-terminus from labeling - it is <strong>not</strong>
            counted as a labeling failure. These peptides have their labelable site count reduced by 1.</p>

            <h3>Search Configuration</h3>
            <p>This analysis uses a <strong>variable modification</strong> search strategy where TMT is
            set as a variable (not fixed) modification. This allows detection of both labeled and unlabeled
            peptides, enabling accurate efficiency calculation.</p>
        </div>

        <h2 id="interpretation">📖 Interpretation Guide</h2>

        <div class="method-box">
            <h3>Why Does Labeling Efficiency Matter?</h3>
            <p>In TMT experiments, quantification relies on comparing reporter ion intensities across channels.
            If labeling is incomplete:</p>
            <ul>
                <li><strong>Unlabeled peptides</strong> contribute no signal to any channel</li>
                <li><strong>Partially labeled peptides</strong> may have skewed ratios</li>
                <li><strong>Channel cross-talk</strong> effects are amplified</li>
            </ul>

            <h3>What Causes Low Labeling Efficiency?</h3>
            <table>
                <tr><th>Cause</th><th>Solution</th></tr>
                <tr><td>Expired or degraded TMT reagent</td><td>Use fresh reagent, store at -20°C desiccated</td></tr>
                <tr><td>Insufficient reagent amount</td><td>Use 0.8mg TMT per 100μg protein minimum</td></tr>
                <tr><td>Incorrect pH</td><td>Ensure pH 8.5 with HEPES or TEAB buffer</td></tr>
                <tr><td>Amine contamination</td><td>Avoid Tris, ammonium salts; use TEAB buffer</td></tr>
                <tr><td>Incomplete quenching</td><td>Quench with hydroxylamine for 15+ minutes</td></tr>
                <tr><td>Sample aggregation</td><td>Add SDS or use higher urea concentration</td></tr>
            </table>

            <h3>Expected Efficiency Values</h3>
            <ul>
                <li><strong>Excellent:</strong> ≥97% (PASS - well-optimized protocol)</li>
                <li><strong>Acceptable:</strong> 90-97% (WARNING - standard conditions, review recommended)</li>
                <li><strong>Poor:</strong> &lt;90% (FAIL - re-labeling recommended)</li>
            </ul>

            <h3>Site Efficiency vs Fully Labeled PSMs</h3>
            <p>These two metrics measure different aspects of labeling:</p>
            <table>
                <tr><th>Metric</th><th>What it Counts</th><th>Formula</th></tr>
                <tr><td><strong>Site Labeling Efficiency</strong></td><td>Individual sites (K + N-term)</td><td>labeled_sites / total_sites</td></tr>
                <tr><td><strong>Fully Labeled PSMs</strong></td><td>Whole peptides with ALL sites labeled</td><td>fully_labeled_PSMs / total_PSMs</td></tr>
            </table>
            <p><strong>Why they differ:</strong> Site efficiency is typically higher because even partially labeled PSMs
            contribute their labeled sites to the numerator. A peptide with 3 sites where 2 are labeled counts as
            0% for "Fully Labeled" but contributes 67% (2/3) to Site Efficiency.</p>
            <p><strong>Which matters more?</strong> Site Efficiency is the standard QC metric for overall labeling completeness.
            Fully Labeled % provides additional insight into how many peptides achieved perfect labeling.</p>
        </div>

        <footer>
            <p><strong>nf-fragpipe TMT QC Tool v1.0.0</strong></p>
            <p>FragPipe-based TMT labeling efficiency analysis pipeline</p>
            <p style="font-size: 12px; color: #95a5a6;">
                Analysis pipeline: ThermoRawFileParser → MSFragger → MSBooster → Percolator → ProteinProphet → Philosopher
            </p>
        </footer>
    </div>
</body>
</html>"""

    output_path.write_text(html_content)


def generate_markdown_report(
    efficiency_data: Dict[str, Any],
    oxidation_data: Dict[str, Any],
    mixing_table: List[Dict[str, Any]],
    skipped_channels: List[str],
    tmt_type: str,
    pass_threshold: float,
    warn_threshold: float,
    output_path: Path,
) -> None:
    """Generate Markdown report."""

    status, _ = determine_qc_status(
        efficiency_data["labeling_efficiency"], pass_threshold, warn_threshold
    )

    status_emoji = {"PASS": "✅", "WARNING": "⚠️", "FAIL": "❌"}.get(status, "❓")

    # Generate mixing table if available
    mixing_md = ""
    if mixing_table:
        mixing_rows = "\n".join(
            [
                f"| {row['channel']} | {row['total_intensity']:.2e} | {row['correction_factor']:.3f} | {row['percent_of_mean']:.1f}% |"
                for row in mixing_table
            ]
        )

        # Add skipped channels note if any
        skipped_note = ""
        if skipped_channels:
            skipped_note = f"""
> **⚠️ Skipped Channels:** {", ".join(skipped_channels)}
>
> These {len(skipped_channels)} channel(s) had zero intensity and are excluded from mixing ratio calculations.
"""

        mixing_md = f"""
## Channel Mixing Ratios

| Channel | Total Intensity | Correction Factor | % of Mean |
|---------|-----------------|-------------------|-----------|
{mixing_rows}
{skipped_note}"""

    # Generate per-sample tables for markdown (separate for labeling and oxidation)
    per_sample_md = ""
    per_sample_data = efficiency_data.get("per_sample", [])
    oxidation_per_sample = {
        s["sample"]: s for s in oxidation_data.get("per_sample", [])
    }

    # Extract experiment ID from common prefix
    # For single sample, use the full sample name as experiment ID
    sample_names = [row["sample"] for row in per_sample_data]
    stripped_names, common_prefix = strip_common_prefix(sample_names)
    if len(sample_names) == 1:
        experiment_id = sample_names[0]
    else:
        experiment_id = common_prefix.rstrip("_") if common_prefix else None

    if per_sample_data and len(per_sample_data) >= 1:
        name_map = dict(zip(sample_names, stripped_names))

        # Labeling table rows
        labeling_rows_list = []
        for row in per_sample_data:
            display_name = name_map.get(row["sample"], row["sample"])
            labeling_rows_list.append(
                f"| {display_name} | {row['psm_count']:,} | {row['efficiency']:.1f}% | {row.get('fully_labeled_pct', 0):.1f}% | {row.get('partially_labeled_pct', 0):.1f}% | {row.get('not_labeled_pct', 0):.1f}% |"
            )
        labeling_rows = "\n".join(labeling_rows_list)

        # Oxidation table rows
        oxidation_rows_list = []
        for row in per_sample_data:
            display_name = name_map.get(row["sample"], row["sample"])
            ox_data = oxidation_per_sample.get(row["sample"], {})
            ox_rate = ox_data.get("site_oxidation_rate", 0)
            not_ox_pct = ox_data.get("not_oxidized_pct", 0)
            partial_ox_pct = ox_data.get("partially_oxidized_pct", 0)
            fully_ox_pct = ox_data.get("fully_oxidized_pct", 0)
            m_psms = ox_data.get("total_m_psms", 0)
            oxidation_rows_list.append(
                f"| {display_name} | {m_psms:,} | {ox_rate:.1f}% | {fully_ox_pct:.1f}% | {partial_ox_pct:.1f}% | {not_ox_pct:.1f}% |"
            )
        oxidation_rows = "\n".join(oxidation_rows_list)

        prefix_note = (
            f"*Common prefix removed: `{common_prefix}`*\n\n" if common_prefix else ""
        )

        per_sample_md = f"""
## Per-Sample Quality

{prefix_note}### TMT Labeling

| Sample | PSMs | Site Eff. | Fully | Partial | None |
|--------|------|-----------|-------|---------|------|
{labeling_rows}

### Methionine Oxidation

| Sample | Met PSMs | Site Ox | Fully | Partial | Not Ox |
|--------|----------|---------|------|---------|--------|
{oxidation_rows}
"""

    md_content = f"""# TMT Labeling QC Report

{f"**Experiment:** {experiment_id}" + chr(10) if experiment_id else ""}**TMT Type:** {tmt_type}
**Generated:** {datetime.now().strftime("%Y-%m-%d")}

---

## Summary

| Metric | Value |
|--------|-------|
| **Site Labeling Efficiency** | **{efficiency_data["labeling_efficiency"]:.1f}%** {status_emoji} {status} |
| **Met Oxidation Rate** | **{oxidation_data["site_oxidation_rate"]:.1f}%** |
| Total PSMs | {efficiency_data["total_psms"]:,} |
| Samples | {len(per_sample_data)} |
| Sites Labeled (TMT) | {efficiency_data["total_labeled_sites"]:,} |
| Sites Unlabeled | {efficiency_data["total_unlabeled_sites"]:,} |

## QC Status

{status_emoji} **{status}** - Labeling efficiency {efficiency_data["labeling_efficiency"]:.1f}% {"exceeds" if status == "PASS" else "is below"} threshold of {pass_threshold * 100:.0f}%

## QC Thresholds

### TMT Labeling Efficiency

| Status | Labeling Threshold |
|--------|-----------|
| ✅ PASS | ≥ {pass_threshold * 100:.0f}% |
| ⚠️ WARNING | {warn_threshold * 100:.0f}% – {pass_threshold * 100:.0f}% |
| ❌ FAIL | < {warn_threshold * 100:.0f}% |
{mixing_md}
{per_sample_md}
---

*Generated by nf-fragpipe TMT QC Tool v1.0.0*
"""

    output_path.write_text(md_content)


def generate_summary_tsv(
    efficiency_data: Dict[str, Any], status: str, output_path: Path
) -> None:
    """Generate machine-readable TSV summary."""

    lines = [
        "metric\tvalue\tpercentage",
        f"total_psms\t{efficiency_data['total_psms']}\t100.0",
        f"filtered_psms\t{efficiency_data['filtered_psms']}\t-",
        f"total_labelable_sites\t{efficiency_data['total_labelable_sites']}\t-",
        f"labeled_sites\t{efficiency_data['total_labeled_sites']}\t{efficiency_data['total_labeled_sites'] / efficiency_data['total_labelable_sites'] * 100:.1f}"
        if efficiency_data["total_labelable_sites"] > 0
        else f"labeled_sites\t{efficiency_data['total_labeled_sites']}\t0.0",
        f"unlabeled_sites\t{efficiency_data['total_unlabeled_sites']}\t{efficiency_data['total_unlabeled_sites'] / efficiency_data['total_labelable_sites'] * 100:.1f}"
        if efficiency_data["total_labelable_sites"] > 0
        else f"unlabeled_sites\t{efficiency_data['total_unlabeled_sites']}\t0.0",
        f"nterm_acetylated\t{efficiency_data['nterm_blocked_count']}\t-",
        f"labeling_efficiency\t{efficiency_data['labeling_efficiency']:.2f}\t-",
        f"qc_status\t{status}\t-",
    ]

    output_path.write_text("\n".join(lines) + "\n")


def generate_mixing_csv(mixing_table: List[Dict[str, Any]], output_path: Path) -> None:
    """Generate mixing ratios CSV."""

    lines = ["channel,total_intensity,correction_factor,percent_of_mean"]
    for row in mixing_table:
        lines.append(
            f"{row['channel']},{row['total_intensity']:.6e},{row['correction_factor']:.4f},{row['percent_of_mean']:.2f}"
        )

    output_path.write_text("\n".join(lines) + "\n")


def generate_per_sample_csv(
    per_sample_data: List[Dict[str, Any]], output_path: Path
) -> None:
    """Generate per-sample labeling efficiency CSV."""

    lines = [
        "sample,psm_count,total_sites,labeled_sites,unlabeled_sites,nterm_blocked,efficiency"
    ]
    for row in per_sample_data:
        lines.append(
            f"{row['sample']},{row['psm_count']},{row['total_sites']},"
            f"{row['labeled_sites']},{row['unlabeled_sites']},{row['nterm_blocked']},{row['efficiency']:.2f}"
        )

    output_path.write_text("\n".join(lines) + "\n")


def cmd_analyze(args: argparse.Namespace) -> int:
    """Execute the analyze command."""

    # Validate TMT type
    tmt_type = args.tmt_type.upper()
    if tmt_type not in TMT_MASSES:
        print(f"Error: Unknown TMT type '{args.tmt_type}'")
        print(f"Valid types: {', '.join(TMT_MASSES.keys())}")
        return 1

    tmt_mass = TMT_MASSES[tmt_type]

    # Check PSM file exists
    psm_path = Path(args.psm_file)
    if not psm_path.exists():
        print(f"Error: PSM file not found: {psm_path}")
        return 1

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate workflow if present
    if args.workflow_file:
        workflow_path = Path(args.workflow_file)
        validate_workflow(workflow_path, tmt_mass)

    print("TMT QC Analysis")
    print("=" * 50)
    print(f"TMT Type: {tmt_type} (mass: {tmt_mass})")
    print(f"PSM File: {psm_path}")
    print(f"Probability Threshold: {args.prob_threshold}")
    print()

    # Load and analyze PSM data
    print("Loading PSM data...")
    psm_data = parse_tsv_file(psm_path)
    print(f"  Loaded {len(psm_data)} PSMs")

    print("Calculating labeling efficiency...")
    efficiency_data = calculate_labeling_efficiency(
        psm_data, tmt_mass, prob_threshold=args.prob_threshold
    )

    print("Calculating oxidation metrics...")
    oxidation_data = calculate_oxidation_efficiency(
        psm_data, prob_threshold=args.prob_threshold
    )

    # Calculate mixing ratios if abundance file provided
    mixing_table = []
    skipped_channels = []
    if args.abundance_file:
        abundance_path = Path(args.abundance_file)
        if abundance_path.exists():
            print("Loading abundance data...")
            abundance_data = parse_tsv_file(abundance_path)
            print(f"  Loaded {len(abundance_data)} peptides")

            print("Calculating mixing ratios...")
            mixing_table, skipped_channels = calculate_mixing_ratios(
                abundance_data, tmt_type
            )

    # Determine QC status
    status, status_class = determine_qc_status(
        efficiency_data["labeling_efficiency"], args.pass_threshold, args.warn_threshold
    )

    # Generate outputs
    print("Generating reports...")

    generate_html_report(
        efficiency_data,
        oxidation_data,
        mixing_table,
        skipped_channels,
        tmt_type,
        args.pass_threshold,
        args.warn_threshold,
        output_dir / "report.html",
    )

    generate_markdown_report(
        efficiency_data,
        oxidation_data,
        mixing_table,
        skipped_channels,
        tmt_type,
        args.pass_threshold,
        args.warn_threshold,
        output_dir / "report.md",
    )

    generate_summary_tsv(efficiency_data, status, output_dir / "labeling_summary.tsv")

    if mixing_table:
        generate_mixing_csv(mixing_table, output_dir / "mixing_table.csv")

    # Generate per-sample CSV
    per_sample_data = efficiency_data.get("per_sample", [])
    if per_sample_data and len(per_sample_data) >= 1:
        generate_per_sample_csv(
            per_sample_data, output_dir / "per_sample_efficiency.csv"
        )

    # Print summary to console
    print()
    print("=" * 50)
    print("RESULTS")
    print("=" * 50)
    print(f"Labeling Efficiency: {efficiency_data['labeling_efficiency']:.1f}%")
    print(f"Met Oxidation Rate: {oxidation_data['site_oxidation_rate']:.1f}%")
    print(f"QC Status: {status}")
    print(f"  - Total PSMs: {efficiency_data['total_psms']:,}")
    print(f"  - Samples: {len(per_sample_data)}")
    print(f"  - Labelable Sites: {efficiency_data['total_labelable_sites']:,}")
    print(f"  - Labeled Sites: {efficiency_data['total_labeled_sites']:,}")
    print(f"  - Met-containing PSMs: {oxidation_data['total_m_psms']:,}")
    print(f"  - N-term Acetylated: {efficiency_data['nterm_blocked_count']:,}")
    if mixing_table:
        print(f"  - Active Channels: {len(mixing_table)}")
    if skipped_channels:
        print(
            f"  - Skipped Channels: {len(skipped_channels)} ({', '.join(skipped_channels)})"
        )

    # Print per-sample summary with oxidation
    oxidation_per_sample = {
        s["sample"]: s for s in oxidation_data.get("per_sample", [])
    }
    if per_sample_data and len(per_sample_data) >= 1:
        print()
        print("Per-Sample Quality:")
        for sample in per_sample_data:
            label_status = (
                "✅"
                if sample["efficiency"] >= args.pass_threshold * 100
                else (
                    "⚠️" if sample["efficiency"] >= args.warn_threshold * 100 else "❌"
                )
            )
            ox_data = oxidation_per_sample.get(sample["sample"], {})
            ox_rate = ox_data.get("site_oxidation_rate", 0)
            print(
                f"  {label_status} {sample['sample']}: Label {sample['efficiency']:.1f}%, Ox {ox_rate:.1f}%"
            )

    print()
    print(f"Reports saved to: {output_dir}")
    print("  - report.html")
    print("  - report.md")
    print("  - labeling_summary.tsv")
    if mixing_table:
        print("  - mixing_table.csv")
    if per_sample_data and len(per_sample_data) >= 1:
        print("  - per_sample_efficiency.csv")

    return 0


def main() -> int:
    """Main entry point."""

    parser = argparse.ArgumentParser(
        description="TMT QC Analysis Tool - Calculate labeling efficiency and mixing ratios",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  tmt_qc.py analyze --psm-file psm.tsv --tmt-type TMT16 --output-dir results/

  # With mixing ratio analysis
  tmt_qc.py analyze --psm-file psm.tsv --abundance-file abundance.tsv --tmt-type TMT10 --output-dir results/

  # Custom thresholds
  tmt_qc.py analyze --psm-file psm.tsv --tmt-type TMT16 --pass-threshold 0.90 --warn-threshold 0.75 --output-dir results/

TMT Types: TMT0, TMT2, TMT6, TMT10, TMT11, TMT16, TMT18, TMT35, TMTPRO
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Analyze subcommand
    analyze_parser = subparsers.add_parser(
        "analyze", help="Analyze TMT labeling efficiency and mixing ratios"
    )
    analyze_parser.add_argument(
        "--psm-file",
        "-p",
        required=True,
        help="Path to PSM TSV file (Philosopher psm.tsv or Percolator target_psms.tsv)",
    )
    analyze_parser.add_argument(
        "--abundance-file",
        "-a",
        help="Path to peptide abundance TSV file (optional, for mixing ratios)",
    )
    analyze_parser.add_argument(
        "--workflow-file",
        "-w",
        help="Path to fragpipe.workflow file (optional, for validation)",
    )
    analyze_parser.add_argument(
        "--tmt-type",
        "-t",
        required=True,
        help="TMT reagent type (TMT6, TMT10, TMT11, TMT16, TMT18, TMTPRO)",
    )
    analyze_parser.add_argument(
        "--prob-threshold",
        type=float,
        default=0.95,
        help="Minimum probability for PSM inclusion (default: 0.95)",
    )
    analyze_parser.add_argument(
        "--pass-threshold",
        type=float,
        default=0.97,
        help="Efficiency threshold for PASS status (default: 0.97)",
    )
    analyze_parser.add_argument(
        "--warn-threshold",
        type=float,
        default=0.90,
        help="Efficiency threshold for WARNING status (default: 0.90)",
    )
    analyze_parser.add_argument(
        "--output-dir",
        "-o",
        default=".",
        help="Output directory for reports (default: current directory)",
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "analyze":
        return cmd_analyze(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
