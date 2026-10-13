"""reprint_to_saint_input.py vs FragPipe 24.0 CmdSaintExpress.reprintTsvToSaintInput (CmdSaintExpress.java:98-167).

Expected bytes were produced by FragPipe's own method (FragPipe 24.0 lib/*) on the same table.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "modules/local/reprinttosaintinput/resources/usr/bin"))
import reprint_to_saint_input as rts  # noqa: E402

TABLE = "PROTID\tGENEID\tPROTLEN\tB1_SPC\tC1_SPC\nna\tna\tna\tBAIT\tcontrol\nP1\tG1\t100\t3\t0\nP2\tG2\t200\t0.0\t2\n"


@pytest.fixture
def converted(tmp_path):
    src = tmp_path / "reprint.spc.tsv"
    src.write_text(TABLE)
    rts.convert(src, tmp_path / "out")
    return tmp_path / "out"


class TestReprintToSaintInput:
    def test_inter_skips_zero_and_ends_with_blank_line(self, converted):
        # CmdSaintExpress.java:136-141, :152-158
        assert (converted / "inter.dat").read_text() == "B1\tBAIT\tP1\t3\nC1\tC1\tP2\t2\n\n"

    def test_bait_control_uses_ip_name(self, converted):
        # CmdSaintExpress.java:121-129: CONTROL -> IP name, type C when the name starts with CONTROL
        assert (converted / "bait.dat").read_text() == "B1\tBAIT\tT\nC1\tC1\tT\n"

    def test_prey_with_protlen_and_gene(self, converted):
        # CmdSaintExpress.java:133-134
        assert (converted / "prey.dat").read_text() == "P1\t100\tG1\nP2\t200\tG2\n"

    def test_columns_must_end_in_int_or_spc(self, tmp_path):
        src = tmp_path / "reprint.spc.tsv"
        src.write_text("PROTID\tGENEID\tX\nna\tna\tBAIT\n")
        with pytest.raises(AssertionError):
            rts.convert(src, tmp_path / "out")


class TestMissingTable:
    """CmdSaintExpress.java:72-79: an inaccessible table logs FragPipe's error and skips the mode."""

    def test_missing_table_writes_fragpipe_message_and_exits_zero(self, tmp_path, monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        assert rts.main(["reprint.int.tsv", "saintexpress-int"]) == 0
        message = f"Input file {tmp_path / 'reprint.int.tsv'} is not accessible."
        assert (tmp_path / "saintexpress-int" / "skipped_mode.txt").read_text() == message + "\n"
        assert capsys.readouterr().err == message + "\n"

    def test_missing_table_writes_no_saint_inputs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        rts.main(["reprint.spc.tsv", "saintexpress-spc"])
        assert sorted(p.name for p in (tmp_path / "saintexpress-spc").iterdir()) == ["skipped_mode.txt"]

    def test_directory_is_not_a_table(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "reprint.spc.tsv").mkdir()
        rts.main(["reprint.spc.tsv", "saintexpress-spc"])
        assert (tmp_path / "saintexpress-spc" / "skipped_mode.txt").exists()

    def test_present_table_writes_the_three_inputs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "reprint.spc.tsv").write_text(TABLE)
        assert rts.main(["reprint.spc.tsv", "saintexpress-spc"]) == 0
        assert sorted(p.name for p in (tmp_path / "saintexpress-spc").iterdir()) == ["bait.dat", "inter.dat", "prey.dat"]
