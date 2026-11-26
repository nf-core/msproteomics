"""TMT label check QC report (tmt_qc.py): the "Expected Efficiency Values" list of the HTML report shows the
--pass-threshold / --warn-threshold values, and the CLI defaults are the pipeline's thresholds
(TMTLABELCHECKANALYZE ext.pass_threshold / ext.warn_threshold in conf/modules.config)."""

import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_tmt_qc_acetyl import load

REPO = Path(__file__).resolve().parents[2]

# Efficiency far from every threshold pair below, so that no other text of the report has the threshold numbers.
EFFICIENCY = 73.0


def per_sample(name, efficiency):
    return {
        "sample": name,
        "psm_count": 100,
        "total_sites": 200,
        "labeled_sites": int(2 * efficiency),
        "unlabeled_sites": 200 - int(2 * efficiency),
        "nterm_blocked": 0,
        "efficiency": efficiency,
        "fully_labeled": int(efficiency),
        "partially_labeled": 0,
        "not_labeled": 100 - int(efficiency),
        "fully_labeled_pct": efficiency,
        "partially_labeled_pct": 0.0,
        "not_labeled_pct": 100.0 - efficiency,
    }


def efficiency_data():
    return {
        "total_psms": 100,
        "filtered_psms": 0,
        "total_labelable_sites": 200,
        "total_labeled_sites": int(2 * EFFICIENCY),
        "total_unlabeled_sites": 200 - int(2 * EFFICIENCY),
        "nterm_blocked_count": 0,
        "labeling_efficiency": EFFICIENCY,
        "fully_labeled": int(EFFICIENCY),
        "partially_labeled": 0,
        "not_labeled": 100 - int(EFFICIENCY),
        "fully_labeled_pct": EFFICIENCY,
        "partially_labeled_pct": 0.0,
        "not_labeled_pct": 100.0 - EFFICIENCY,
        "per_peptide": [],
        "per_sample": [per_sample("S1", EFFICIENCY)],
    }


def oxidation_data():
    return {
        "total_m_psms": 0,
        "total_not_oxidized": 0,
        "total_m_sites": 0,
        "total_oxidized_sites": 0,
        "non_oxidation_rate": 100.0,
        "site_oxidation_rate": 0.0,
        "per_sample": [],
    }


def pipeline_thresholds(name):
    """Values of ext.<name> in conf/modules.config (set for TMTLABELCHECKANALYZE and _IONQUANT only)."""
    text = (REPO / "conf/modules.config").read_text()
    return [float(v) for v in re.findall(rf"ext\.{name}\s*=\s*([0-9.]+)", text)]


class TestExpectedEfficiencyList(unittest.TestCase):
    qc = load("tmt_qc")

    def report(self, pass_threshold, warn_threshold):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "report.html"
            self.qc.generate_html_report(
                efficiency_data(), oxidation_data(), [], [], "TMT16", pass_threshold, warn_threshold, out
            )
            return out.read_text(encoding="utf-8")

    def check(self, pass_threshold, warn_threshold, pass_pct, warn_pct, absent):
        html = self.report(pass_threshold, warn_threshold)
        self.assertIn(f"<strong>Excellent:</strong> ≥{pass_pct}% (PASS", html)
        self.assertIn(f"<strong>Acceptable:</strong> {warn_pct}-{pass_pct}% (WARNING", html)
        self.assertIn(f"<strong>Poor:</strong> &lt;{warn_pct}% (FAIL", html)
        for pct in absent:
            self.assertNotIn(f"{pct}%", html)

    def test_pipeline_thresholds(self):
        self.check(0.95, 0.85, 95, 85, absent=(97, 90))

    def test_other_thresholds(self):
        # a second pair: fixed "95%" / "85%" text would pass the first test but not this one
        self.check(0.80, 0.60, 80, 60, absent=(97, 95, 90, 85))


class TestCliDefaults(unittest.TestCase):
    qc = load("tmt_qc")

    def parsed(self):
        argv = ["tmt_qc.py", "analyze", "--psm-file", "psm.tsv", "--tmt-type", "TMT16"]
        with mock.patch.object(sys, "argv", argv), mock.patch.object(self.qc, "cmd_analyze", return_value=0) as cmd:
            self.assertEqual(self.qc.main(), 0)
        return cmd.call_args.args[0]

    def test_defaults_are_pipeline_thresholds(self):
        # TMTLABELCHECKANALYZE and TMTLABELCHECKANALYZE_IONQUANT
        self.assertEqual(pipeline_thresholds("pass_threshold"), [0.95, 0.95])
        self.assertEqual(pipeline_thresholds("warn_threshold"), [0.85, 0.85])
        args = self.parsed()
        self.assertEqual((args.pass_threshold, args.warn_threshold), (0.95, 0.85))


if __name__ == "__main__":
    unittest.main()
