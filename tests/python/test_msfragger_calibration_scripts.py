"""sort_fasta.py and parse_calibration_output.py vs FragPipe 24.0 tools/msfragger_pep_split.py.

sort_fasta.py is sample_fasta() (lines 614-620) with n = 1; parse_calibration_output.py is the
post-search half of calibrate() (lines 524-611).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "modules/local/sortfasta/resources/usr/bin"))
sys.path.insert(0, str(REPO / "modules/local/msfraggerapplycalibration/resources/usr/bin"))
import parse_calibration_output as pco  # noqa: E402
import sort_fasta  # noqa: E402


def fragpipe_sample_fasta(data: bytes) -> bytes:
    """msfragger_pep_split.py:615-620 verbatim, n = 1."""
    fasta_prots = [b">" + e + b"\n" for e in data[1:].split(b"\n>")]
    return b"".join(sorted(fasta_prots)[::1])


class TestSortFasta:
    @pytest.mark.parametrize(
        "data",
        [
            b">sp|B|b\nMKB\n>rev_sp|A|a\nMKA\n>sp|A|a\nMKA\n",
            b">only\nMK",
            b">z\nAA\nBB\n>a\nCC\n",
        ],
    )
    def test_matches_fragpipe_sample_fasta(self, data):
        assert sort_fasta.sorted_fasta(data) == fragpipe_sample_fasta(data)

    def test_cli_writes_the_sorted_file(self, tmp_path):
        src, dst = tmp_path / "db.fasta", tmp_path / "sorted_db.fasta"
        src.write_bytes(b">b\nB\n>a\nA\n")
        assert sort_fasta.main([str(src), str(dst)]) == 0
        assert dst.read_bytes() == b">a\nA\n\n>b\nB\n"


class TestCalibratedParams:
    PARAMS = "calibrate_mass = 2\nfragment_mass_tolerance = 20\nfragment_mass_units = 0\nuse_topN_peaks = 300\ncheck_spectral_files = 1\t\t\t# d\n"

    def test_calibrate_mass_2_applies_new_values_and_disables_spectral_check(self):
        out = "New fragment_mass_tolerance = 12.000000 PPM\nNew use_topN_peaks = 150\n"
        got = pco.calibrated_params(self.PARAMS, out)
        assert got == "calibrate_mass = 2\nfragment_mass_tolerance = 12.000000\nfragment_mass_units = 1\nuse_topN_peaks = 150\ncheck_spectral_files = 0\t\t\t# d\n"

    def test_calibrate_mass_1_only_disables_spectral_check(self):
        params = self.PARAMS.replace("calibrate_mass = 2", "calibrate_mass = 1")
        assert pco.calibrated_params(params, "") == params.replace("check_spectral_files = 1", "check_spectral_files = 0")

    def test_missing_fragment_tolerance_fails_like_fragpipe(self):
        with pytest.raises(ValueError):
            pco.calibrated_params(self.PARAMS, "")

    def test_spectra_replaced_by_mzbin_calibrated_or_copied(self, tmp_path):
        (tmp_path / "a.mzML").write_text("a")
        (tmp_path / "a.mzBIN_calibrated").write_text("cal")
        (tmp_path / "b.mzML").write_text("b")
        out = tmp_path / "out"
        pco.collect_calibrated_spectra([str(tmp_path / "a.mzML"), str(tmp_path / "b.mzML")], out)
        assert sorted(p.name for p in out.iterdir()) == ["a.mzBIN_calibrated", "b.mzML"]
