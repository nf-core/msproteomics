"""merge_split_search.py names the full FASTA in the merged pepXML header and needs it in the task.

msfragger_pep_split.py:92-94 takes fasta_path from the params' database_name (the path the tools run
with) and :302-303 replaces the chunk FASTA name with it. nf-core/msproteomics gave the merge a
head-node S3 string (/<bucket>/<key>), which no task can open, so ProteinProphet wrote prot_length 0.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "modules/local/mergesplitsearch/resources/usr/bin"))
import merge_split_search as mss  # noqa: E402

CHUNK_HEADER = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<msms_pipeline_analysis date="2026-10-13T10:00:00" summary_xml="{chunk}/s.pepXML">\n'
    '<msms_run_summary base_name="{chunk}/s" raw_data_type="mzML" raw_data="mzML">\n'
    '<search_summary base_name="{chunk}/s" precursor_mass_type="monoisotopic">\n'
    '<search_database local_path="{chunk_fasta}" type="AA"/>\n'
    '<parameter name="database_name" value="{chunk_fasta}"/>\n'
    "</search_summary>\n"
    '<spectrum_query spectrum="s.1.1.2" start_scan="1">\n'
)


class TestPepxmlHeaderFasta(unittest.TestCase):
    def header(self, chunk: str, chunk_fasta: str, fasta_path: str) -> bytes:
        with TemporaryDirectory() as d:
            p = Path(d) / "s.pepXML"
            p.write_text(CHUNK_HEADER.format(chunk=chunk, chunk_fasta=chunk_fasta))
            return mss.get_pepxml_header(p, fasta_path)

    def test_fusion_chunk_fasta_becomes_the_merge_task_fasta(self):
        chunk = "/fusion/s3/bucket/work/f8/19983c"
        merged = "/fusion/s3/bucket/work/03/9b6640/db.fasta"
        header = self.header(chunk, f"{chunk}/db.fasta", merged)
        self.assertIn(f'<search_database local_path="{merged}" type="AA"/>'.encode(), header)
        self.assertIn(f'<parameter name="database_name" value="{merged}"/>'.encode(), header)
        self.assertNotIn(f"{chunk}/db.fasta".encode(), header)

    def test_relative_chunk_fasta_becomes_the_full_path(self):
        # FragPipe's chunk params name the FASTA without a directory (msfragger_pep_split.py:141-143)
        header = self.header("/w/c", "db.fasta", "/tmp/run/db.fasta")
        self.assertIn(b'<search_database local_path="/tmp/run/db.fasta" type="AA"/>', header)
        self.assertNotIn(b'"db.fasta"', header)

    def test_date_and_summary_xml_are_placeholders(self):
        header = self.header("/w/c", "/w/c/db.fasta", "/w/m/db.fasta")
        self.assertIn(b'date="%b"', header)
        self.assertIn(b'summary_xml="%b"', header)


class TestRequireFasta(unittest.TestCase):
    def test_existing_fasta_passes(self):
        with TemporaryDirectory() as d:
            fasta = Path(d) / "db.fasta"
            fasta.write_text(">p\nPEPTIDEK\n")
            mss.require_fasta(str(fasta))

    def test_head_node_s3_string_exits_with_error(self):
        with self.assertRaisesRegex(SystemExit, "is not a file in this task"):
            mss.require_fasta("/bucket/work/27/525735/db.fasta")


if __name__ == "__main__":
    unittest.main()
