"""Unit tests for tests/tools/golden_coverage.py."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

import golden_coverage as cov  # noqa: E402

MANIFEST = {
    "steps": [
        {"step_index": 1, "name": "Percolator", "cmd_name": "Percolator", "cwd": "{WORKDIR}/C02"},
        {"step_index": 2, "name": "Percolator", "cmd_name": "Percolator", "cwd": "{WORKDIR}/D02"},
        {"step_index": 3, "name": "ProteinProphet", "cmd_name": "ProteinProphet", "cwd": "{WORKDIR}"},
    ]
}


class TestCoverage(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.repo = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def write_test(self, rel: str, text: str) -> None:
        p = self.repo / rel
        p.parent.mkdir(parents=True)
        p.write_text(text)

    def test_cwd_selector_covers_only_that_instance(self):
        self.write_test("modules/local/a/tests/golden.nf.test", "tag \"golden\"\nit.cmd_name == 'Percolator' && it.cwd == '{WORKDIR}/C02'\n")
        self.write_test("modules/local/b/tests/golden.nf.test", "tag \"golden\"\nit.cmd_name == 'ProteinProphet'\n")
        rows = cov.coverage(MANIFEST, cov.find_tests(self.repo))
        self.assertEqual([bool(r["tests"]) for r in rows], [True, False, True])

    def test_golden_test_without_selector_is_an_error(self):
        self.write_test("modules/local/a/tests/golden.nf.test", 'tag "golden"\n')
        with self.assertRaisesRegex(ValueError, "selects no manifest step"):
            cov.find_tests(self.repo)


if __name__ == "__main__":
    unittest.main()
