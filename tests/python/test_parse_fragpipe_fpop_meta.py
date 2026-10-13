"""FPOP Quant and Metaproteomics outputs of parse_fragpipe_workflow.py (FragPipe 24.0).

generate_fpop_config gives FragPipe_FPOP_Analysis.py's positional argv (CmdFpopQuant.java:90-99);
generate_metaproteomics_config gives the DbOptimizer values the module puts in CmdMetaproteomics order.
"""

import importlib.util
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "modules/local/parsefragpipeworkflow/resources/usr/bin/parse_fragpipe_workflow.py"

spec = importlib.util.spec_from_file_location("parse_fragpipe_workflow", SCRIPT)
pfw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pfw)


class TestFpopArgv(unittest.TestCase):
    def test_label_free_defaults_keep_empty_labels(self):
        self.assertEqual(pfw.generate_fpop_config({}, {}), ["combined_modified_peptide.tsv", "1", "", "", "false", "false"])

    def test_label_free_values(self):
        params = {"fragpipe.fpop.region_size": "2", "fragpipe.fpop.label_control": "CTRL", "fragpipe.fpop.label_fpop": "OX 1",
                  "fragpipe.fpop.subtract-control": "TRUE"}
        self.assertEqual(pfw.generate_fpop_config(params, {}), ["combined_modified_peptide.tsv", "2", "CTRL", "OX 1", "true", "false"])

    def test_tmt_reads_the_normalisation_tables(self):
        params = {"fragpipe.fpop.fpop-tmt": "true"}
        tmti = {"tmtintegrator": {"run-tmtintegrator": "true", "prot_norm": "2"}}
        self.assertEqual(
            pfw.generate_fpop_config(params, tmti),
            ["tmt-report/abundance_multi-site_GN.tsv", "1", "", "", "false", "true", "tmt-report-unmod/abundance_peptide_GN.tsv"],
        )

    def test_tmt_without_tmt_integrator_fails_like_fragpipe(self):
        with self.assertRaisesRegex(ValueError, "TMT-Integrator is not run"):
            pfw.generate_fpop_config({"fragpipe.fpop.fpop-tmt": "true"}, {"tmtintegrator": {"run-tmtintegrator": "false"}})


class TestMetaproteomicsValues(unittest.TestCase):
    def test_defaults(self):
        self.assertEqual(
            pfw.generate_metaproteomics_config({}),
            {"qvalue": "0.01", "deltaHyperscore": "0.0", "minPeptCntPerProt": "1", "minUniqPeptCntPerProt": "1",
             "minUniqPeptCnt": "3", "hostName": '"Homo sapiens"', "iterations": "3", "cmdLineOpts": ""},
        )

    def test_host_name_quoted_once_and_opts_trimmed(self):
        got = pfw.generate_metaproteomics_config({"host-name": '"Mus musculus"', "cmd-line-opts": "  --foo 1 "})
        self.assertEqual((got["hostName"], got["cmdLineOpts"]), ('"Mus musculus"', "--foo 1"))


if __name__ == "__main__":
    unittest.main()
