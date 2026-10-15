"""Unit tests for tests/tools/tower_compare.py (path prefixes of Tower/Fusion runs, pairing) and the ragged-row read of
tests/tools/philosopher_tables.py (FragPipe's own psm.tsv has rows with fewer fields than the header)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

import philosopher_tables as pt  # noqa: E402
import tower_compare as tc  # noqa: E402

TASK = "/fusion/s3/bucket-a/{}/ab/0123456789abcdef0123456789abcd"


def found(text):
    return [m.decode() for m in tc.PATTERN.findall(text.encode())]


def test_scratch_task_dir():
    d = TASK.format("scratch/4XkpD3D4pOSiNq")
    assert found(f"x {d}/file.pepXML y") == [d]


def test_any_work_dir_prefix():
    for sub in ("project/run1/PXD014415/work", "a/b/c/work", "work"):
        d = TASK.format(sub)
        assert found(f"database_name = {d}/db.fasta") == [d], sub


def test_s3_url_task_dir():
    line = b"database_name = s3://my-bucket/p/work/9b/75c3927bd47ac2481526276b588793/db.fasta\t\t\t# x"
    assert tc.PATTERN.sub(b"<RUN>", line) == b"database_name = <RUN>/db.fasta\t\t\t# x"


def test_project_resolved_maps_like_the_workdir():
    d = TASK.format("x/work")
    assert found(f"projectDir = {d}/project_resolved") == [f"{d}/project_resolved"]
    assert found(f"projectDir = {d}/project") == [f"{d}/project"]
    assert found(f"fastaFile = {d}/projects.fasta") == [d]


def test_nextflow_scratch_tmp_dirs():
    assert found("a /tmp/nxf.04u5PJ9AKZ/results/x b /tmp/nxf.04u5PJ9AKZ/raw_files/y") == [
        "/tmp/nxf.04u5PJ9AKZ/results", "/tmp/nxf.04u5PJ9AKZ/raw_files"]
    assert found("/tmp/other/x") == []


def test_input_dirs_of_s3_urls_and_fusion_paths():
    assert found("s3://bucket-b/SDI/p/PXD1/inputs/raw/a.raw") == ["/bucket-b/SDI/p/PXD1/inputs/raw"]
    assert found("/fusion/s3/bucket-b/p/PXD1/inputs/d/a.d") == ["/fusion/s3/bucket-b/p/PXD1/inputs/d"]
    assert found("/usr/share/x/inputs/y") == []


def test_is_task_dir():
    assert tc.is_task_dir(TASK.format("project/run1/EXP/work"))
    assert not tc.is_task_dir("/fusion/s3/bucket-b/x/inputs")


def test_choose_candidate_prefers_the_last_rewriting_stage(tmp_path):
    names = ["philosopher/C02/psm.tsv", "ionquant_ms1/fragpipe/C02/psm.tsv", "ionquant_isobaric/fragpipe/C02/psm.tsv"]
    cands = [tmp_path / n for n in names]
    assert tc.choose_candidate(cands, tmp_path / "gold/C02/psm.tsv") == [cands[2]]
    glyco = [tmp_path / n for n in ("philosopher/stce/psm.tsv", "opair/stce/psm.tsv")]
    assert tc.choose_candidate(glyco, tmp_path / "gold/stce/psm.tsv") == [glyco[1]]
    ptms = [tmp_path / n for n in ("philosopher/A01/psm.tsv", "ptmshepherd/groups/A01/psm.tsv")]
    assert tc.choose_candidate(ptms, tmp_path / "gold/A01/psm.tsv") == [ptms[1]]


def test_choose_candidate_parent_dir_name_first(tmp_path):
    cands = [tmp_path / "x/A01/ion.tsv", tmp_path / "x/B01/ion.tsv"]
    assert tc.choose_candidate(cands, tmp_path / "gold/B01/ion.tsv") == [cands[1]]


def test_choose_candidate_pairs_fragger_params_with_the_run_level_copy(tmp_path):
    for name in ("fragger.params", "fragger_dia.params"):
        cands = [tmp_path / f"out/msfragger/A01/{name}", tmp_path / f"out/msfragger/{name}", tmp_path / f"out/msfragger/B01/{name}"]
        assert tc.choose_candidate(cands, tmp_path / f"gold/results/{name}") == [cands[1]]


def test_choose_candidate_leaves_per_file_fragger_params_unpaired_without_a_run_level_copy(tmp_path):
    cands = [tmp_path / "out/msfragger/A01/fragger_dia.params", tmp_path / "out/msfragger/B01/fragger_dia.params"]
    assert len(tc.choose_candidate(cands, tmp_path / "gold/results/fragger_dia.params")) == 2


def test_table_pads_missing_trailing_fields(tmp_path):
    p = tmp_path / "psm.tsv"
    p.write_text("a\tb\tc\n1\t2\t3\n4\t5\n")
    t = pt.Table(p)
    assert t.rows == [["1", "2", "3"], ["4", "5", ""]]


def test_table_still_rejects_extra_fields(tmp_path):
    p = tmp_path / "psm.tsv"
    p.write_text("a\tb\n1\t2\t3\n")
    try:
        pt.Table(p)
    except ValueError as e:
        assert "has 3 fields, header has 2" in str(e)
    else:
        raise AssertionError("a row longer than the header must still be an error")


def test_headless_roots_from_fragpipe_job(tmp_path):
    (tmp_path / "fragpipe.job").write_text("runtime-save\t/tmp/nxf.AbC123xyz9/results/fragpipe-workflow.workflow\t/tmp/nxf.AbC123xyz9/results\n")
    roots = tc.headless_roots(tmp_path)
    assert roots["{WORKDIR}"] == "/tmp/nxf.AbC123xyz9/results"
    assert roots["{LCMS_DIR}"] == "/tmp/nxf.AbC123xyz9/raw_files"
    assert roots["{DB_DIR}"] == "/tmp/nxf.AbC123xyz9"
    (tmp_path / "fragpipe.job").write_text("no task dir here\n")
    try:
        tc.headless_roots(tmp_path)
    except ValueError as e:
        assert "no /tmp/nxf" in str(e)
    else:
        raise AssertionError("a fragpipe.job without a task dir must be an error")
