"""merge_split_search.py fails at once when --msfragger_cmd names no existing MSFragger jar.

nf-core/msproteomics 5b43675 added the check; FragPipe's msfragger_cmd puts JVM options between
-jar and the jar (msfragger_pep_split.py:73), so the jar is the first non-option token after -jar.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "modules/local/mergesplitsearch/resources/usr/bin"))
import merge_split_search as mss  # noqa: E402


class TestRequireMsfraggerJar(unittest.TestCase):
    def test_existing_jar_after_jvm_options_passes(self):
        with TemporaryDirectory() as d:
            jar = Path(d) / "MSFragger-4.4.1.jar"
            jar.write_text("")
            mss.require_msfragger_jar(["java", "-jar", "-Dfile.encoding=UTF-8", "-Xmx4G", str(jar), "--license", "l"])

    def test_missing_jar_exits_with_error(self):
        with self.assertRaisesRegex(SystemExit, "MSFragger jar missing"):
            mss.require_msfragger_jar(["java", "-jar", "-Dfile.encoding=UTF-8", "/nonexistent/MSFragger-4.4.1.jar"])

    def test_no_jar_after_options_exits_with_error(self):
        with self.assertRaisesRegex(SystemExit, "got ''"):
            mss.require_msfragger_jar(["java", "-jar", "-Dfile.encoding=UTF-8"])

    def test_command_without_jar_flag_is_not_checked(self):
        mss.require_msfragger_jar(["msfragger", "--version"])


if __name__ == "__main__":
    unittest.main()
