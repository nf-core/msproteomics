#!/usr/bin/env python3
"""TMT QC Analysis Tool (IonQuant version) - Calculate TMT labeling efficiency from combined_modified_peptide.tsv."""

# =============================================================================
# Imports
# =============================================================================
import argparse
import base64
import io
import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# TMT mass configurations
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


def parse_tsv_file(file_path: Path) -> Tuple[List[str], List[Dict[str, Any]]]:
    """Parse a TSV file into headers and list of dictionaries."""
    rows = []
    headers = []
    with open(file_path, "r") as f:
        lines = f.readlines()
        if not lines:
            return headers, rows

        headers = lines[0].strip().split("\t")
        for line in lines[1:]:
            if line.strip():
                values = line.strip().split("\t")
                row = {}
                for i, header in enumerate(headers):
                    row[header] = values[i] if i < len(values) else ""
                rows.append(row)
    return headers, rows


def extract_sample_columns(headers: List[str]) -> List[str]:
    """Extract sample names from column headers (e.g., 'Sample1 Spectral Count')."""
    samples = []
    for header in headers:
        if header.endswith(" Spectral Count"):
            sample_name = header.replace(" Spectral Count", "")
            samples.append(sample_name)
    return samples


def count_tmt_modifications_ionquant(modified_sequence: str, tmt_mass: float) -> int:
    """Count TMT modifications from IonQuant's Modified Sequence format.

    IonQuant format: n[304.2071]PEPTIDEK[304.2071]R
    """
    if not modified_sequence or modified_sequence == "nan":
        return 0

    modified_sequence = str(modified_sequence)

    # Try exact mass match first
    tmt_mass_str = f"{tmt_mass:.4f}"
    count = modified_sequence.count(f"[{tmt_mass_str}]")

    # Try shorter precision
    if count == 0:
        tmt_mass_short = f"{tmt_mass:.2f}"
        count = modified_sequence.count(f"[{tmt_mass_short}]")

    # Try integer mass
    if count == 0:
        mass_int = str(int(tmt_mass))
        count = modified_sequence.count(f"[{mass_int}]")

    return count


def has_nterm_acetylation_ionquant(modified_sequence: str) -> bool:
    """Check if peptide has N-terminal acetylation from Modified Sequence."""
    if not modified_sequence or modified_sequence == "nan":
        return False

    modified_sequence = str(modified_sequence)

    # Check for acetylation at N-terminus (first ~20 chars)
    nterm_region = modified_sequence[:30]

    acetyl_patterns = [
        "[42.0106]",
        "[42.01]",
        "[42]",
        "n[42",
    ]

    return any(pattern in nterm_region for pattern in acetyl_patterns)


def count_oxidized_m(modified_sequence: str) -> int:
    """Count oxidized methionines from Modified Sequence."""
    if not modified_sequence or modified_sequence == "nan":
        return 0

    modified_sequence = str(modified_sequence)

    # Look for M[15.9949] or similar patterns
    count = 0
    count += modified_sequence.count("M[15.99")
    count += modified_sequence.count("M[16]")
    count += modified_sequence.count("M[147]")  # Met + oxidation combined mass

    return count


def find_common_prefix(names: List[str]) -> str:
    """Find common prefix among sample names."""
    if not names or len(names) < 2:
        return ""

    prefix = names[0]
    for name in names[1:]:
        while not name.startswith(prefix) and prefix:
            prefix = prefix[:-1]

    if prefix:
        for sep in ["_", "-", "."]:
            last_sep = prefix.rfind(sep)
            if last_sep > 0:
                prefix = prefix[: last_sep + 1]
                break

    return prefix


def strip_common_prefix(names: List[str]) -> Tuple[List[str], str]:
    """Strip common prefix from sample names for cleaner display."""
    prefix = find_common_prefix(names)
    if prefix:
        stripped = [
            name[len(prefix) :] if name.startswith(prefix) else name for name in names
        ]
        return stripped, prefix
    return names, ""


def calculate_labeling_efficiency_ionquant(
    peptide_data: List[Dict[str, Any]],
    headers: List[str],
    tmt_mass: float,
) -> Dict[str, Any]:
    """Calculate per-peptide and per-sample TMT labeling efficiency from IonQuant data."""
    results = []

    # Detect column names
    sample_row = peptide_data[0] if peptide_data else {}
    peptide_col = "Peptide Sequence" if "Peptide Sequence" in sample_row else "Peptide"
    modified_col = (
        "Modified Sequence" if "Modified Sequence" in sample_row else "Modified Peptide"
    )
    assigned_mods_col = (
        "Assigned Modifications" if "Assigned Modifications" in sample_row else None
    )

    # Extract sample names from headers
    sample_names = extract_sample_columns(headers)

    # Initialize per-sample stats
    sample_stats: Dict[str, Dict[str, int]] = {}
    for sample_name in sample_names:
        sample_stats[sample_name] = {
            "peptide_count": 0,
            "total_sites": 0,
            "labeled_sites": 0,
            "nterm_blocked": 0,
            "fully_labeled": 0,
            "partially_labeled": 0,
            "not_labeled": 0,
        }

    for row in peptide_data:
        peptide = row.get(peptide_col, "")
        if not peptide:
            continue

        modified_sequence = row.get(modified_col, "")

        # Calculate labelable sites: K residues + N-terminus
        k_count = peptide.upper().count("K")
        total_sites = k_count + 1

        # Count TMT modifications
        labeled_sites = count_tmt_modifications_ionquant(modified_sequence, tmt_mass)

        # Also try from Assigned Modifications if available
        if labeled_sites == 0 and assigned_mods_col and assigned_mods_col in row:
            assigned_mods = row.get(assigned_mods_col, "")
            # Format: 8K(304.2071), n(304.2071)
            if assigned_mods:
                tmt_mass_str = f"{tmt_mass:.4f}"
                labeled_sites = assigned_mods.count(tmt_mass_str)

        # Check for N-terminal acetylation
        nterm_blocked = has_nterm_acetylation_ionquant(modified_sequence)
        if nterm_blocked:
            total_sites -= 1

        total_sites = max(0, total_sites)
        labeled_sites = min(labeled_sites, total_sites)

        # Determine labeling status
        if labeled_sites == 0:
            label_status = "not_labeled"
        elif labeled_sites >= total_sites:
            label_status = "fully_labeled"
        else:
            label_status = "partially_labeled"

        # Aggregate by sample using spectral counts
        for sample_name in sample_names:
            spec_count_col = f"{sample_name} Spectral Count"
            if spec_count_col in row:
                spec_count = int(float(row[spec_count_col]))

                if spec_count > 0:
                    sample_stats[sample_name]["peptide_count"] += spec_count
                    sample_stats[sample_name]["total_sites"] += total_sites * spec_count
                    sample_stats[sample_name]["labeled_sites"] += (
                        labeled_sites * spec_count
                    )
                    if nterm_blocked:
                        sample_stats[sample_name]["nterm_blocked"] += spec_count
                    sample_stats[sample_name][label_status] += spec_count

        results.append(
            {
                "peptide": peptide,
                "k_count": k_count,
                "total_sites": total_sites,
                "labeled_sites": labeled_sites,
                "nterm_blocked": nterm_blocked,
                "unlabeled_sites": total_sites - labeled_sites,
                "label_status": label_status,
            }
        )

    # Calculate totals
    total_labelable = sum(s["total_sites"] for s in sample_stats.values())
    total_labeled = sum(s["labeled_sites"] for s in sample_stats.values())
    total_unlabeled = total_labelable - total_labeled
    nterm_blocked_count = sum(s["nterm_blocked"] for s in sample_stats.values())
    total_psms = sum(s["peptide_count"] for s in sample_stats.values())

    total_fully_labeled = sum(s["fully_labeled"] for s in sample_stats.values())
    total_partially_labeled = sum(s["partially_labeled"] for s in sample_stats.values())
    total_not_labeled = sum(s["not_labeled"] for s in sample_stats.values())

    efficiency = (total_labeled / total_labelable * 100) if total_labelable > 0 else 0.0

    # Per-sample summary
    per_sample = []
    for sample_name in sorted(sample_stats.keys()):
        stats = sample_stats[sample_name]
        if stats["peptide_count"] == 0:
            continue

        sample_efficiency = (
            (stats["labeled_sites"] / stats["total_sites"] * 100)
            if stats["total_sites"] > 0
            else 0.0
        )
        psm_total = stats["peptide_count"]
        per_sample.append(
            {
                "sample": sample_name,
                "psm_count": stats["peptide_count"],
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

    return {
        "total_psms": total_psms,
        "filtered_psms": 0,  # No filtering at peptide level
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


def calculate_oxidation_efficiency_ionquant(
    peptide_data: List[Dict[str, Any]],
    headers: List[str],
) -> Dict[str, Any]:
    """Calculate methionine oxidation rate from IonQuant data."""
    sample_row = peptide_data[0] if peptide_data else {}
    peptide_col = "Peptide Sequence" if "Peptide Sequence" in sample_row else "Peptide"
    modified_col = (
        "Modified Sequence" if "Modified Sequence" in sample_row else "Modified Peptide"
    )

    sample_names = extract_sample_columns(headers)

    sample_stats: Dict[str, Dict[str, int]] = {}
    for sample_name in sample_names:
        sample_stats[sample_name] = {
            "total_m_psms": 0,
            "not_oxidized": 0,
            "partially_oxidized": 0,
            "fully_oxidized": 0,
            "total_m_sites": 0,
            "oxidized_m_sites": 0,
        }

    for row in peptide_data:
        peptide = row.get(peptide_col, "")
        if not peptide:
            continue

        m_count = peptide.upper().count("M")
        if m_count == 0:
            continue

        modified_sequence = row.get(modified_col, "")
        ox_count = count_oxidized_m(modified_sequence)

        for sample_name in sample_names:
            spec_count_col = f"{sample_name} Spectral Count"
            if spec_count_col in row:
                spec_count = int(float(row[spec_count_col]))

                if spec_count > 0:
                    sample_stats[sample_name]["total_m_psms"] += spec_count
                    sample_stats[sample_name]["total_m_sites"] += m_count * spec_count
                    sample_stats[sample_name]["oxidized_m_sites"] += (
                        ox_count * spec_count
                    )

                    if ox_count == 0:
                        sample_stats[sample_name]["not_oxidized"] += spec_count
                    elif ox_count >= m_count:
                        sample_stats[sample_name]["fully_oxidized"] += spec_count
                    else:
                        sample_stats[sample_name]["partially_oxidized"] += spec_count

    total_m_psms = sum(s["total_m_psms"] for s in sample_stats.values())
    total_not_oxidized = sum(s["not_oxidized"] for s in sample_stats.values())
    total_m_sites = sum(s["total_m_sites"] for s in sample_stats.values())
    total_oxidized_sites = sum(s["oxidized_m_sites"] for s in sample_stats.values())

    non_oxidation_rate = (
        (total_not_oxidized / total_m_psms * 100) if total_m_psms > 0 else 100.0
    )
    site_oxidation_rate = (
        (total_oxidized_sites / total_m_sites * 100) if total_m_sites > 0 else 0.0
    )

    per_sample = []
    for sample_name in sorted(sample_stats.keys()):
        stats = sample_stats[sample_name]
        if stats["total_m_psms"] == 0:
            continue

        total_psms = stats["total_m_psms"]
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


def determine_qc_status(
    efficiency: float, pass_threshold: float, warn_threshold: float
) -> Tuple[str, str]:
    """Determine QC status based on labeling efficiency."""
    if efficiency >= pass_threshold * 100:
        return "PASS", "pass"
    elif efficiency >= warn_threshold * 100:
        return "WARNING", "warning"
    else:
        return "FAIL", "fail"


def generate_labeling_plots(
    efficiency_data: Dict[str, Any], oxidation_data: Dict[str, Any]
) -> Optional[str]:
    """Generate TMT labeling and oxidation plots. Returns base64 PNG or None."""
    COLORS = {
        "fully_labeled": "#5AC2F1",
        "partially_labeled": "#FBCF35",
        "not_labeled": "#ED4C1C",
        "fully_oxidized": "#5AC2F1",
        "partially_oxidized": "#FBCF35",
        "not_oxidized": "#ED4C1C",
    }

    per_sample = efficiency_data.get("per_sample", [])
    ox_per_sample = oxidation_data.get("per_sample", [])

    if not per_sample:
        return None

    sample_names = [s["sample"] for s in per_sample]
    if len(sample_names) == 1:
        full_name = sample_names[0]
        parts = full_name.split("_")
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

    labeling_data = []
    for s in per_sample:
        display_name = name_map.get(s["sample"], s["sample"])
        labeling_data.append(
            {
                "sample": display_name,
                "fully_labeled": s.get("fully_labeled_pct", 0),
                "partially_labeled": s.get("partially_labeled_pct", 0),
                "not_labeled": s.get("not_labeled_pct", 0),
            }
        )

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
            }
        )

    # Create figure
    n_samples = len(labeling_data)
    bar_height = 0.6
    row_height = 0.5
    plot_height = max(2.5, n_samples * row_height + 1)
    fig_height = plot_height * 2 + 0.5

    fig, (ax_label_pct, ax_ox_pct) = plt.subplots(2, 1, figsize=(10, fig_height))

    samples = [d["sample"] for d in labeling_data]
    y_pos = list(range(len(samples)))

    def add_bar_labels(
        ax, y_positions, values, x_offsets, color="white", fontsize=10, min_width=3
    ):
        for y, val, x_off in zip(y_positions, values, x_offsets):
            if abs(val) >= min_width:
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

    # TMT Labeling plot
    not_labeled_vals = [d["not_labeled"] for d in labeling_data]
    not_labeled_pct = [-v for v in not_labeled_vals]
    partial_pct = [d["partially_labeled"] for d in labeling_data]
    fully_pct = [d["fully_labeled"] for d in labeling_data]

    ax_label_pct.barh(
        y_pos,
        not_labeled_pct,
        height=bar_height,
        color=COLORS["not_labeled"],
        alpha=0.9,
        label="Not Labeled",
    )
    ax_label_pct.barh(
        y_pos,
        partial_pct,
        height=bar_height,
        color=COLORS["partially_labeled"],
        alpha=0.9,
        label="Partially Labeled",
    )
    ax_label_pct.barh(
        y_pos,
        fully_pct,
        left=partial_pct,
        height=bar_height,
        color=COLORS["fully_labeled"],
        alpha=0.9,
        label="Fully Labeled",
    )

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
    ax_label_pct.set_yticklabels(samples if n_samples > 1 else [""], fontsize=11)
    ax_label_pct.set_xlabel("Percent", fontsize=11)
    ax_label_pct.set_title("TMT Labeling [%]", fontweight="bold", fontsize=13)
    ax_label_pct.axvline(x=0, color="gray", linewidth=0.8)
    ax_label_pct.legend(loc="lower right", fontsize=9)
    ax_label_pct.set_xlim(-105, 105)
    ax_label_pct.invert_yaxis()

    # Oxidation plot
    has_ox_data = any(
        d.get("not_oxidized", 0)
        + d.get("partially_oxidized", 0)
        + d.get("fully_oxidized", 0)
        > 0
        for d in oxidation_data_list
    )

    if has_ox_data:
        not_ox_vals = [d.get("not_oxidized", 0) for d in oxidation_data_list]
        not_ox_pct = [-v for v in not_ox_vals]
        partial_ox_pct = [d.get("partially_oxidized", 0) for d in oxidation_data_list]
        fully_ox_pct = [d.get("fully_oxidized", 0) for d in oxidation_data_list]

        ax_ox_pct.barh(
            y_pos,
            not_ox_pct,
            height=bar_height,
            color=COLORS["not_oxidized"],
            alpha=0.9,
            label="Not Oxidized",
        )
        ax_ox_pct.barh(
            y_pos,
            partial_ox_pct,
            height=bar_height,
            color=COLORS["partially_oxidized"],
            alpha=0.9,
            label="Partially Oxidized",
        )
        ax_ox_pct.barh(
            y_pos,
            fully_ox_pct,
            left=partial_ox_pct,
            height=bar_height,
            color=COLORS["fully_oxidized"],
            alpha=0.9,
            label="Fully Oxidized",
        )

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
        ax_ox_pct.set_yticklabels(samples if n_samples > 1 else [""], fontsize=11)
        ax_ox_pct.set_xlabel("Percent", fontsize=11)
        ax_ox_pct.set_title("Met Oxidation [%]", fontweight="bold", fontsize=13)
        ax_ox_pct.axvline(x=0, color="gray", linewidth=0.8)
        ax_ox_pct.legend(loc="lower right", fontsize=9)
        ax_ox_pct.set_xlim(-105, 105)
        ax_ox_pct.invert_yaxis()
    else:
        ax_ox_pct.text(
            0.5,
            0.5,
            "No Met-containing peptides",
            ha="center",
            va="center",
            transform=ax_ox_pct.transAxes,
            fontsize=12,
            color="gray",
        )
        ax_ox_pct.set_title("Met Oxidation [%]", fontweight="bold", fontsize=13)

    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def generate_html_report(
    efficiency_data: Dict[str, Any],
    oxidation_data: Dict[str, Any],
    tmt_type: str,
    pass_threshold: float,
    warn_threshold: float,
    output_path: Path,
) -> None:
    """Generate comprehensive HTML report."""

    status, status_class = determine_qc_status(
        efficiency_data["labeling_efficiency"], pass_threshold, warn_threshold
    )

    status_emoji = {
        "PASS": "&#x2705;",
        "WARNING": "&#x26A0;&#xFE0F;",
        "FAIL": "&#x274C;",
    }.get(status, "&#x2753;")

    tmt_mass = TMT_MASSES.get(tmt_type.upper(), 229.1629)
    per_sample_data = efficiency_data.get("per_sample", [])
    sample_names = [s["sample"] for s in per_sample_data]
    if len(sample_names) == 1:
        experiment_id = sample_names[0]
    else:
        _, experiment_id = strip_common_prefix(sample_names)
        experiment_id = experiment_id.rstrip("_") if experiment_id else None

    combined_plot_b64 = generate_labeling_plots(efficiency_data, oxidation_data)
    oxidation_per_sample = {
        s["sample"]: s for s in oxidation_data.get("per_sample", [])
    }

    # Build list of problematic samples (below pass threshold)
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

    if status == "PASS":
        recommendation = '<div class="recommendation pass-bg"><h3>&#x2705; Proceed with Analysis</h3><p>Labeling meets threshold. Use TMT as fixed modification.</p></div>'
    elif status == "WARNING":
        recommendation = f'<div class="recommendation warning-bg"><h3>&#x26A0;&#xFE0F; Review Before Proceeding</h3><p>Efficiency {efficiency_data["labeling_efficiency"]:.1f}% below {pass_threshold * 100:.0f}% threshold.</p></div>'
    else:
        recommendation = f'<div class="recommendation fail-bg"><h3>&#x274C; Do Not Proceed</h3><p>Efficiency {efficiency_data["labeling_efficiency"]:.1f}% below {warn_threshold * 100:.0f}%. Re-label samples.</p></div>'

    plots_html = ""
    if combined_plot_b64:
        plots_html = f'<h2>&#x1F4CA; Visual Summary</h2><div class="plot-container"><img src="data:image/png;base64,{combined_plot_b64}" alt="TMT QC"></div>'

    # Per-sample tables
    per_sample_html = ""
    if per_sample_data:

        def get_efficiency_class(eff: float) -> str:
            if eff >= pass_threshold * 100:
                return "pass"
            elif eff >= warn_threshold * 100:
                return "warning"
            return "fail"

        stripped_names, common_prefix = strip_common_prefix(sample_names)
        name_map = dict(zip(sample_names, stripped_names))
        labeling_rows = []
        for row in per_sample_data:
            display_name = name_map.get(row["sample"], row["sample"])
            labeling_rows.append(
                f"<tr><td>{display_name}</td><td>{row['psm_count']:,}</td>"
                f'<td class="{get_efficiency_class(row["efficiency"])}">{row["efficiency"]:.1f}%</td>'
                f'<td class="pass">{row.get("fully_labeled_pct", 0):.1f}%</td>'
                f'<td class="warning">{row.get("partially_labeled_pct", 0):.1f}%</td>'
                f'<td class="fail">{row.get("not_labeled_pct", 0):.1f}%</td></tr>'
            )

        oxidation_rows = []
        for row in per_sample_data:
            display_name = name_map.get(row["sample"], row["sample"])
            ox_data = oxidation_per_sample.get(row["sample"], {})
            oxidation_rows.append(
                f"<tr><td>{display_name}</td><td>{ox_data.get('total_m_psms', 0):,}</td>"
                f"<td>{ox_data.get('site_oxidation_rate', 0):.1f}%</td>"
                f'<td class="fail">{ox_data.get("fully_oxidized_pct", 0):.1f}%</td>'
                f'<td class="warning">{ox_data.get("partially_oxidized_pct", 0):.1f}%</td>'
                f'<td class="pass">{ox_data.get("not_oxidized_pct", 0):.1f}%</td></tr>'
            )

        prefix_note = (
            f'<p style="font-size: 12px; color: #7f8c8d;">Common prefix removed: <code>{common_prefix}</code></p>'
            if common_prefix
            else ""
        )

        per_sample_html = f"""
        <h2>&#x1F9EB; Per-Sample Quality</h2>
        {prefix_note}
        <h3>TMT Labeling</h3>
        <table class="centered-table">
            <tr><th>Sample</th><th>PSMs</th><th>Site Eff.</th><th>Fully</th><th>Partial</th><th>None</th></tr>
            {chr(10).join(labeling_rows)}
        </table>
        <h3>Methionine Oxidation</h3>
        <table class="centered-table">
            <tr><th>Sample</th><th>Met PSMs</th><th>Site Ox</th><th>Fully</th><th>Partial</th><th>Not Ox</th></tr>
            {chr(10).join(oxidation_rows)}
        </table>
        """

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>TMT QC (Combined) - {tmt_type}</title>
    <style>
        body{{font-family:system-ui,sans-serif;margin:0;padding:40px;background:#f5f7fa;color:#2c3e50;line-height:1.6}}
        .container{{max-width:1000px;margin:0 auto;background:#fff;padding:50px;border-radius:12px;box-shadow:0 2px 10px rgba(0,0,0,.1)}}
        h1{{border-bottom:3px solid #3498db;padding-bottom:15px;margin-top:0}}
        h2{{color:#34495e;margin-top:40px;border-bottom:1px solid #ecf0f1;padding-bottom:10px}}
        table{{border-collapse:collapse;width:100%;margin:20px 0}}
        th,td{{border:1px solid #ddd;padding:12px;text-align:center}}
        th{{background:#3498db;color:#fff}}
        tr:nth-child(even){{background:#f9f9f9}}
        .pass{{color:#27ae60}}.warning{{color:#f39c12}}.fail{{color:#e74c3c}}
        .info-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:20px;margin:25px 0}}
        .info-card{{background:#ecf0f1;padding:20px;border-radius:8px;text-align:center}}
        .info-card .value{{font-size:28px;font-weight:bold}}.info-card .label{{font-size:14px;color:#7f8c8d}}
        .recommendation{{padding:25px;border-radius:8px;margin:25px 0}}
        .pass-bg{{background:#d5f4e6;border-left:5px solid #27ae60}}
        .warning-bg{{background:#fef9e7;border-left:5px solid #f39c12}}
        .fail-bg{{background:#fdedec;border-left:5px solid #e74c3c}}
        .plot-container{{padding:20px;border-radius:8px;margin:20px 0;border:1px solid #e0e0e0;text-align:center}}
        .plot-container img{{max-width:100%;height:auto}}
        .metric-cards{{display:flex;gap:20px;margin:25px 0;flex-wrap:wrap}}
        .metric-card{{flex:1;min-width:200px;padding:25px;border-radius:12px;text-align:center;color:#fff}}
        .metric-card.labeling{{background:linear-gradient(135deg,#667eea,#764ba2)}}
        .metric-card.oxidation{{background:linear-gradient(135deg,#11998e,#38ef7d)}}
        .metric-card .metric-value{{font-size:48px;font-weight:bold;margin:10px 0}}
        .metric-card .metric-label{{font-size:14px;opacity:.9}}.metric-card .metric-status{{font-size:24px;margin-top:10px}}
        .alert-box{{background:#fff3cd;border:1px solid #ffc107;border-radius:8px;padding:15px;margin:20px 0}}
        .alert-box.warning{{background:#fff3cd;border-color:#ffc107}}
        .alert-box.danger{{background:#f8d7da;border-color:#dc3545}}
        .alert-box h4{{margin:0 0 10px 0;color:#856404}}
        .alert-box.danger h4{{color:#721c24}}
        .alert-box ul{{margin:0;padding-left:20px}}
        footer{{margin-top:40px;padding-top:20px;border-top:1px solid #ecf0f1;color:#7f8c8d;font-size:14px;text-align:center}}
    </style>
</head>
<body>
    <div class="container">
        <h1>&#x1F9EA; TMT QC Report (Combined - IonQuant)</h1>
        <p>{
        f"<strong>Experiment:</strong> {experiment_id} &nbsp;|&nbsp; "
        if experiment_id
        else ""
    }<strong>TMT:</strong> {tmt_type} ({tmt_mass:.4f} Da) &nbsp;|&nbsp;
           <strong>Generated:</strong> {datetime.now().strftime("%Y-%m-%d")}</p>

        <h2>&#x1F4C8; Results Summary</h2>

        {
        "".join(
            [
                f'''
        <div class="alert-box {'danger' if any(s['labeling_status'] == 'FAIL' for s in problematic_samples) else 'warning'}">
            <h4>&#x26A0;&#xFE0F; Samples with Labeling Efficiency Issues</h4>
            <ul>
            {"".join([f"<li><strong>{s['sample'].split('_')[-1] if '_' in s['sample'] else s['sample']}</strong>: {s['labeling_efficiency']:.1f}% ({s['labeling_status']})</li>" for s in problematic_samples][:5])}
            {f"<li>...and {len(problematic_samples) - 5} more</li>" if len(problematic_samples) > 5 else ""}
            </ul>
        </div>
        '''
            ]
        )
        if problematic_samples
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
            <div class="info-card"><div class="value">{
        len(per_sample_data)
    }</div><div class="label">Samples</div></div>
            <div class="info-card"><div class="value pass">{
        efficiency_data.get(
            "fully_labeled_pct", 0
        ):.1f}%</div><div class="label">Fully Labeled</div></div>
            <div class="info-card"><div class="value warning">{
        efficiency_data.get(
            "partially_labeled_pct", 0
        ):.1f}%</div><div class="label">Partially Labeled</div></div>
            <div class="info-card"><div class="value fail">{
        efficiency_data.get(
            "not_labeled_pct", 0
        ):.1f}%</div><div class="label">Not Labeled</div></div>
        </div>

        <h2>&#x1F3AF; Recommendation</h2>
        {recommendation}

        {plots_html}

        <h2>&#x1F4CA; Detailed Statistics</h2>
        <table>
            <tr><th>Metric</th><th>Value</th><th>Description</th></tr>
            <tr><td><strong>Site Labeling Efficiency</strong></td><td><strong class="{
        status_class
    }">{
        efficiency_data[
            "labeling_efficiency"
        ]:.2f}%</strong></td><td>Labeled sites / total sites</td></tr>
            <tr><td>Total PSMs (weighted)</td><td>{
        efficiency_data[
            "total_psms"
        ]:,}</td><td>Sum of spectral counts across all peptides</td></tr>
            <tr><td>Total Labelable Sites</td><td>{
        efficiency_data[
            "total_labelable_sites"
        ]:,}</td><td>Sum of (K + N-term) x spectral count</td></tr>
            <tr><td>Sites Labeled</td><td>{
        efficiency_data[
            "total_labeled_sites"
        ]:,}</td><td>Sites with TMT modification</td></tr>
            <tr><td>Sites Unlabeled</td><td>{
        efficiency_data["total_unlabeled_sites"]:,}</td><td>Sites without TMT</td></tr>
        </table>

        {per_sample_html}

        <h2>&#x2696;&#xFE0F; QC Thresholds</h2>
        <table>
            <tr><th>Status</th><th>Threshold</th><th>Interpretation</th></tr>
            <tr><td class="pass">&#x2705; PASS</td><td>&ge; {
        pass_threshold * 100:.0f}%</td><td>Excellent - proceed with analysis</td></tr>
            <tr><td class="warning">&#x26A0;&#xFE0F; WARNING</td><td>{
        warn_threshold * 100:.0f}% &ndash; {
        pass_threshold * 100:.0f}%</td><td>Review before proceeding</td></tr>
            <tr><td class="fail">&#x274C; FAIL</td><td>&lt; {
        warn_threshold * 100:.0f}%</td><td>Re-label samples recommended</td></tr>
        </table>

        <footer><p><strong>nf-fragpipe TMT QC Tool (IonQuant)</strong></p></footer>
    </div>
</body>
</html>"""

    output_path.write_text(html_content)


def generate_markdown_report(
    efficiency_data: Dict[str, Any],
    oxidation_data: Dict[str, Any],
    tmt_type: str,
    pass_threshold: float,
    warn_threshold: float,
    output_path: Path,
) -> None:
    """Generate Markdown report."""

    status, _ = determine_qc_status(
        efficiency_data["labeling_efficiency"], pass_threshold, warn_threshold
    )
    status_emoji = {"PASS": "PASS", "WARNING": "WARNING", "FAIL": "FAIL"}.get(
        status, "?"
    )

    per_sample_data = efficiency_data.get("per_sample", [])
    oxidation_per_sample = {
        s["sample"]: s for s in oxidation_data.get("per_sample", [])
    }

    sample_names = [row["sample"] for row in per_sample_data]
    stripped_names, common_prefix = strip_common_prefix(sample_names)
    experiment_id = common_prefix.rstrip("_") if common_prefix else None

    per_sample_md = ""
    if per_sample_data:
        name_map = dict(zip(sample_names, stripped_names))
        labeling_rows = []
        for row in per_sample_data:
            display_name = name_map.get(row["sample"], row["sample"])
            labeling_rows.append(
                f"| {display_name} | {row['psm_count']:,} | {row['efficiency']:.1f}% | {row.get('fully_labeled_pct', 0):.1f}% | {row.get('partially_labeled_pct', 0):.1f}% | {row.get('not_labeled_pct', 0):.1f}% |"
            )

        per_sample_md = f"""
## Per-Sample Quality

### TMT Labeling

| Sample | PSMs | Site Eff. | Fully | Partial | None |
|--------|------|-----------|-------|---------|------|
{chr(10).join(labeling_rows)}
"""

    md_content = f"""# TMT QC Report (Combined - IonQuant)

{f"**Experiment:** {experiment_id}" + chr(10) if experiment_id else ""}**TMT Type:** {tmt_type}
**Generated:** {datetime.now().strftime("%Y-%m-%d")}

---

## Summary

| Metric | Value |
|--------|-------|
| **Site Labeling Efficiency** | **{efficiency_data["labeling_efficiency"]:.1f}%** {status_emoji} |
| **Met Oxidation Rate** | **{oxidation_data["site_oxidation_rate"]:.1f}%** |
| Total PSMs | {efficiency_data["total_psms"]:,} |
| Samples | {len(per_sample_data)} |

## QC Status

**{status}** - Labeling efficiency {efficiency_data["labeling_efficiency"]:.1f}%
{per_sample_md}
---

*Generated by nf-fragpipe TMT QC Tool (IonQuant) v1.0.0*
"""

    output_path.write_text(md_content)


def generate_summary_tsv(
    efficiency_data: Dict[str, Any], status: str, output_path: Path
) -> None:
    """Generate machine-readable TSV summary."""

    lines = [
        "metric\tvalue\tpercentage",
        f"total_psms\t{efficiency_data['total_psms']}\t100.0",
        f"total_labelable_sites\t{efficiency_data['total_labelable_sites']}\t-",
        f"labeled_sites\t{efficiency_data['total_labeled_sites']}\t{efficiency_data['total_labeled_sites'] / efficiency_data['total_labelable_sites'] * 100:.1f}"
        if efficiency_data["total_labelable_sites"] > 0
        else f"labeled_sites\t{efficiency_data['total_labeled_sites']}\t0.0",
        f"unlabeled_sites\t{efficiency_data['total_unlabeled_sites']}\t{efficiency_data['total_unlabeled_sites'] / efficiency_data['total_labelable_sites'] * 100:.1f}"
        if efficiency_data["total_labelable_sites"] > 0
        else f"unlabeled_sites\t{efficiency_data['total_unlabeled_sites']}\t0.0",
        f"labeling_efficiency\t{efficiency_data['labeling_efficiency']:.2f}\t-",
        f"qc_status\t{status}\t-",
    ]

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


def write_versions_yml(task_process: str, output_path: Path) -> None:
    """Write versions.yml file."""
    python_version = platform.python_version()

    content = (
        f'"{task_process}":\n    python: {python_version}\n    tmt_qc_ionquant: 1.0.0\n'
    )
    output_path.write_text(content)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="TMT QC Analysis Tool (IonQuant version) - Calculate TMT labeling efficiency from combined_modified_peptide.tsv."
    )
    parser.add_argument(
        "--combined-modified-peptide",
        required=True,
        help="Path to combined_modified_peptide.tsv file from IonQuant",
    )
    parser.add_argument(
        "--tmt-type",
        required=True,
        help="TMT reagent type (e.g., TMT10, TMT16, TMTPRO)",
    )
    parser.add_argument(
        "--pass-threshold",
        type=float,
        default=0.95,
        help="Labeling efficiency threshold for PASS status",
    )
    parser.add_argument(
        "--warn-threshold",
        type=float,
        default=0.85,
        help="Labeling efficiency threshold for WARNING status",
    )
    parser.add_argument(
        "--task-process",
        default="TMT_LABELCHECK_ANALYZE_IONQUANT",
        help="Nextflow task process name for versions.yml",
    )
    return parser.parse_args()


def main() -> int:
    """Main entry point for TMT QC analysis (IonQuant version)."""
    args = parse_args()
    tmt_type = args.tmt_type
    pass_threshold = args.pass_threshold
    warn_threshold = args.warn_threshold
    task_process = args.task_process

    tmt_type_upper = tmt_type.upper()
    if tmt_type_upper not in TMT_MASSES:
        print(
            f"Error: Unknown TMT type '{tmt_type}'. Valid types: {', '.join(TMT_MASSES.keys())}"
        )
        return 1

    tmt_mass = TMT_MASSES[tmt_type_upper]
    peptide_path = Path(args.combined_modified_peptide)
    if not peptide_path.exists():
        print(f"Error: Peptide file not found: {peptide_path}")
        return 1

    output_dir = Path(".")

    print("TMT QC Analysis (IonQuant)")
    print("=" * 50)
    print(f"TMT Type: {tmt_type_upper} (mass: {tmt_mass})")
    print(f"Input: {peptide_path}")
    print()

    print("Loading peptide data...")
    headers, peptide_data = parse_tsv_file(peptide_path)
    print(f"  Loaded {len(peptide_data)} peptides")

    sample_names = extract_sample_columns(headers)
    print(f"  Found {len(sample_names)} samples")

    if len(peptide_data) == 0:
        print("No peptides found - generating empty reports")
        (output_dir / "report.html").write_text(
            "<html><body>No peptides found</body></html>"
        )
        (output_dir / "report.md").write_text("No peptides found")
        (output_dir / "labeling_summary.tsv").write_text(
            "metric\tvalue\nqc_status\tNO_DATA"
        )
        (output_dir / "per_sample_efficiency.csv").write_text("sample,psm_count")
        write_versions_yml(task_process, output_dir / "versions.yml")
        return 0

    print("Calculating labeling efficiency...")
    efficiency_data = calculate_labeling_efficiency_ionquant(
        peptide_data, headers, tmt_mass
    )

    print("Calculating oxidation metrics...")
    oxidation_data = calculate_oxidation_efficiency_ionquant(peptide_data, headers)

    status, _ = determine_qc_status(
        efficiency_data["labeling_efficiency"], pass_threshold, warn_threshold
    )

    print("Generating reports...")

    generate_html_report(
        efficiency_data,
        oxidation_data,
        tmt_type_upper,
        pass_threshold,
        warn_threshold,
        output_dir / "report.html",
    )

    generate_markdown_report(
        efficiency_data,
        oxidation_data,
        tmt_type_upper,
        pass_threshold,
        warn_threshold,
        output_dir / "report.md",
    )

    generate_summary_tsv(efficiency_data, status, output_dir / "labeling_summary.tsv")

    per_sample_data = efficiency_data.get("per_sample", [])
    if per_sample_data:
        generate_per_sample_csv(
            per_sample_data, output_dir / "per_sample_efficiency.csv"
        )

    write_versions_yml(task_process, output_dir / "versions.yml")

    print("\n" + "=" * 50 + "\nRESULTS\n" + "=" * 50)
    print(
        f"Labeling: {efficiency_data['labeling_efficiency']:.1f}% | Oxidation: {oxidation_data['site_oxidation_rate']:.1f}% | Status: {status}"
    )
    print(f"Peptides: {len(peptide_data):,} | Samples: {len(per_sample_data)}")

    if per_sample_data:
        print("\nPer-Sample:")
        ox_map = {s["sample"]: s for s in oxidation_data.get("per_sample", [])}
        for s in per_sample_data:
            icon = (
                "PASS"
                if s["efficiency"] >= pass_threshold * 100
                else ("WARN" if s["efficiency"] >= warn_threshold * 100 else "FAIL")
            )
            ox = ox_map.get(s["sample"], {}).get("site_oxidation_rate", 0)
            print(
                f"  [{icon}] {s['sample']}: {s['efficiency']:.1f}% label, {ox:.1f}% ox"
            )

    print(f"\nReports saved to: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
