"""Tests for bin/fragpipe_sdrf.py (port of FragPipe 24.0 SDRF writer)."""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "fragpipe_sdrf.py"
FRAGGER_FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures" / "fragger_params"
_spec = importlib.util.spec_from_file_location("fragpipe_sdrf", SCRIPT)
fs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fs)

OBO = """format-version: 1.4\r
default-namespace: UNIMOD\r
\r
[Term]\r
id: UNIMOD:0\r
name: unimod root node\r
def: "The root node of the unimod modifications ontology." []\r
\r
[Term]\r
id: UNIMOD:1\r
name: Acetyl\r
xref: delta_mono_mass "42.010565"\r
xref: spec_1_site "K"\r
xref: spec_2_site "N-term"\r
xref: spec_3_site "S"\r
\r
[Term]\r
id: UNIMOD:4\r
name: Carbamidomethyl\r
xref: delta_mono_mass "57.021464"\r
xref: spec_1_site "C"\r
xref: spec_2_site "N-term"\r
\r
[Term]\r
id: UNIMOD:35\r
name: Oxidation\r
xref: delta_mono_mass "15.994915"\r
xref: spec_1_site "M"\r
xref: spec_2_site "W"\r
\r
[Term]\r
id: UNIMOD:2016\r
name: TMTpro\r
xref: delta_mono_mass "304.207146"\r
xref: spec_1_site "K"\r
xref: spec_2_site "N-term"\r
\r
[Term]\r
id: UNIMOD:21\r
name: Phospho\r
xref: delta_mono_mass "79.966331"\r
xref: spec_1_site "S"\r
xref: spec_2_site "T"\r
xref: spec_3_site "Y"\r
\r
[Term]\r
id: UNIMOD:9001\r
name: Label:2H(4)+Phospho\r
xref: delta_mono_mass "79.966331"\r
xref: spec_1_site "S"\r
"""

VAR_MODS = (
    "15.9949,M,true,3; 42.0106,[^,true,1; 79.96633,STY,false,3; 0.0,site_10,false,1"
)
FIX_MODS = "0.0,C-Term Peptide,true,-1; 0.0,N-Term Peptide,true,-1; 57.02146,C (cysteine),true,-1; 0.0,K (lysine),true,-1"


def workflow_text(**over):
    props = {
        "msfragger.search_enzyme_name_1": "stricttrypsin",
        "msfragger.search_enzyme_cut_1": "KR",
        "msfragger.search_enzyme_nocut_1": "",
        "msfragger.search_enzyme_sense_1": "C",
        "msfragger.search_enzyme_name_2": "null",
        "msfragger.search_enzyme_cut_2": "",
        "msfragger.search_enzyme_nocut_2": "",
        "msfragger.search_enzyme_sense_2": "C",
        "msfragger.table.var-mods": VAR_MODS,
        "msfragger.table.fix-mods": FIX_MODS,
        "msfragger.mass_offsets": "0",
        "msfragger.mass_offsets_detailed": "",
        "msfragger.use_detailed_offsets": "false",
        "msfragger.restrict_deltamass_to": "all",
        "msfragger.labile_search_mode": "off",
        "msfragger.precursor_mass_units": "1",
        "msfragger.precursor_mass_lower": "-20",
        "msfragger.precursor_mass_upper": "20",
        "msfragger.fragment_mass_units": "1",
        "msfragger.fragment_mass_tolerance": "20",
        "msfragger.calibrate_mass": "2",
        "workflow.misc.save-sdrf": "true",
        "workflow.misc.sdrf-type": "Default",
        "database.db-path": "C\\:\\\\db\\\\x.fasta",
    }
    props.update(over)
    return "# FragPipe workflow\n" + "\n".join("%s=%s" % kv for kv in props.items()) + "\n"


# exact line formats FragPipe parses (MSFragger 4.x log)
LOG_A = (
    "MSFragger version MSFragger-4.3\n"
    "/data/run/a.mzML: Scans = 451; MS2 ITMS = false; MS2 FTMS = true; MS2 ASTMS = false; "
    "MS1 ITMS = false; Isolation sizes = [2.0]; Instrument = Q Exactive\n"
    "fragment_mass_tolerance = 20.0\n"
    "New fragment_mass_tolerance = 5.000000 PPM\n"
)
LOG_B = (
    "b.mzML: Scans = 12; MS2 ITMS = true; MS2 FTMS = false; MS2 ASTMS = false; "
    "MS1 ITMS = false; Isolation sizes = [0.7]; Instrument = Orbitrap Ascend\r\n"
    "New fragment_mass_tolerance = 7.000000 PPM\r\n"
)

FASTA = (
    ">sp|P1|A_HUMAN Alpha OS=Homo sapiens OX=9606 GN=A PE=1 SV=1\nMKR\n"
    ">sp|P2|B_HUMAN Beta OS=Homo sapiens OX=9606 GN=B PE=1 SV=1\nMKR\n"
    ">sp|P3|C_MOUSE Gamma OS=Mus musculus OX=10090 GN=C PE=1 SV=1\nMKR\n"
    ">rev_sp|P3|C_MOUSE no organism here\nMKR\n"
)


@pytest.fixture
def inputs(tmp_path):
    def make(workflow=None, manifest=None, logs=None, fasta=FASTA):
        wf = tmp_path / "fragpipe.workflow"
        wf.write_text(workflow if workflow is not None else workflow_text())
        mf = tmp_path / "fragpipe-files.fp-manifest"
        mf.write_text(manifest if manifest is not None else "/data/run/a.mzML\texp1\t1\tDDA\n/data/run/b.mzML\texp2\t\tDDA")
        fa = tmp_path / "db.fasta"
        fa.write_text(fasta)
        ob = tmp_path / "unimod.obo"
        ob.write_bytes(OBO.encode())
        log_paths = []
        for i, text in enumerate(logs if logs is not None else [LOG_A, LOG_B]):
            p = tmp_path / ("msfragger_%d.log" % i)
            p.write_bytes(text.encode())
            log_paths.append(str(p))
        return dict(workflow=str(wf), manifest=str(mf), fasta=str(fa), unimod=str(ob), logs=log_paths)

    return make


def build(i, tmt_label=""):
    return fs.build_sdrf(i["workflow"], i["manifest"], i["fasta"], i["unimod"], tmt_label, i["logs"])


def rows(text):
    lines = text.split("\n")
    return lines[0].split("\t"), [line.split("\t") for line in lines[1:]]


class TestLfqTable:
    def test_lfq_header_layout_default_type(self, inputs):
        header, _ = rows(build(inputs()))
        assert header == [
            "source name",
            "characteristics[organism]",
            "characteristics[disease]",
            "characteristics[organism part]",
            "characteristics[cell type]",
            "characteristics[biological replicate]",
            "technology type",
            "assay name",
            "comment[data file]",
            "comment[technical replicate]",
            "comment[fraction identifier]",
            "comment[label]",
            "comment[instrument]",
            "comment[precursor mass tolerance]",
            "comment[fragment mass tolerance]",
            "comment[cleavage agent details]",
            "comment[modification parameters]",
            "comment[modification parameters]",
            "comment[modification parameters]",
        ]

    def test_lfq_rows_match_fragpipe_format(self, inputs):
        _, body = rows(build(inputs()))
        assert body[0] == [
            "not available",
            "Homo sapiens",
            "not available",
            "not available",
            "not available",
            "not available",
            "not available",
            "not available",
            "a.mzML",
            "1",
            "not available",
            "not available",
            "NT=Q Exactive",
            "40.0 ppm",
            "5.0 ppm",
            "NT=Trypsin/P;CS=(?<=[KR])",
            "NT=Oxidation;MT=Variable;PP=Anywhere;TA=M;MM=15.99490",
            "NT=Acetyl;MT=Variable;PP=Protein N-term;MM=42.01060",
            "NT=Carbamidomethyl;MT=Fixed;PP=Anywhere;TA=C;MM=57.02146",
        ]

    def test_blank_replicate_is_empty_not_not_available(self, inputs):
        _, body = rows(build(inputs()))
        assert body[1][8:10] == ["b.mzML", ""]

    def test_instrument_from_second_log_with_crlf(self, inputs):
        _, body = rows(build(inputs()))
        assert body[1][12] == "NT=Orbitrap Ascend"

    def test_converted_raw_takes_instrument_of_searched_file(self, inputs, tmp_path):
        # Manifest lists a.raw; MSFragger searched (and logged) the converted a.mzML.
        i = inputs(manifest="/data/run/a.raw\texp1\t1\tDDA\n/data/run/b.mzML\texp2\t\tDDA")
        sn = tmp_path / "searched_files.tsv"
        sn.write_text("a.raw\ta.mzML\nb.mzML\tb.mzML\n")
        text = fs.build_sdrf(i["workflow"], i["manifest"], i["fasta"], i["unimod"], "", i["logs"], str(sn))
        _, body = rows(text)
        assert [body[0][8], body[0][12]] == ["a.raw", "NT=Q Exactive"]

    def test_missing_instrument_is_not_available(self, inputs):
        _, body = rows(build(inputs(logs=[LOG_A])))
        assert body[1][12] == "not available"

    def test_no_trailing_newline(self, inputs, tmp_path):
        i = inputs()
        out = tmp_path / "sdrf.tsv"
        args = ["--workflow", i["workflow"], "--manifest", i["manifest"], "--fasta", i["fasta"]]
        args += ["--unimod", i["unimod"], "--tmt-label", "", "--output", str(out)]
        for log in i["logs"]:
            args += ["--log", log]
        subprocess.run([sys.executable, str(SCRIPT)] + args, check=True)
        data = out.read_bytes()
        assert not data.endswith(b"\n") and data.count(b"\n") == 2

    def test_manifest_windows_paths_and_comments(self, inputs):
        manifest = "# comment\n// other\n\nC:\\raw\\x.raw\te\t02\tDDA\n"
        _, body = rows(build(inputs(manifest=manifest)))
        assert [r[8:10] for r in body] == [["x.raw", "2"]]


class TestTmtTable:
    def test_tmt_rows_one_per_channel(self, inputs):
        _, body = rows(build(inputs(), tmt_label="tmt-6"))
        assert [r[11] for r in body[:6]] == ["TMT126", "TMT127N", "TMT128C", "TMT129N", "TMT130C", "TMT131N"]
        assert len(body) == 12 and body[6][8] == "b.mzML"

    def test_unknown_label_raises(self, inputs):
        with pytest.raises(ValueError):
            build(inputs(), tmt_label="TMT-7")


class TestEnzymes:
    @pytest.mark.parametrize(
        "name,expected",
        [
            ("stricttrypsin", "Trypsin/P"),
            ("trypsin", "Trypsin"),
            ("lysc", "Lys-C"),
            ("nonspecific", "unspecific cleavage"),
            ("null", "no cleavage"),
            ("nocleavage", "nocleavage"),  # FragPipe key is "nocleavage " (trailing space)
            ("mycustom", "mycustom"),
        ],
    )
    def test_enzyme_name_mapping(self, name, expected):
        assert fs.map_enzyme_to_sdrf(name) == expected

    def test_two_enzymes_with_nocut_and_n_sense(self, inputs):
        wf = workflow_text(
            **{
                "msfragger.search_enzyme_name_1": "trypsin",
                "msfragger.search_enzyme_nocut_1": "P",
                "msfragger.search_enzyme_name_2": "aspn",
                "msfragger.search_enzyme_cut_2": "D",
                "msfragger.search_enzyme_sense_2": "N",
            }
        )
        header, body = rows(build(inputs(workflow=wf)))
        assert header.count("comment[cleavage agent details]") == 2
        assert body[0][15:17] == ["NT=Trypsin;CS=(?<=[KR])(?![P])", "NT=Asp-N;CS=(?=[D])"]


class TestMods:
    def test_ambiguous_unimod_falls_back_to_generic_name(self, inputs):
        wf = workflow_text(**{"msfragger.table.var-mods": "79.96633,S,true,3"})
        _, body = rows(build(inputs(workflow=wf)))
        assert body[0][16] == "NT=S:80;MT=Variable;PP=Anywhere;TA=S;MM=79.96633"

    def test_terminal_fixed_mod_and_nterm_var_mod(self, inputs):
        wf = workflow_text(
            **{
                "msfragger.table.var-mods": "304.2071,n^,true,1",
                "msfragger.table.fix-mods": "304.2071,N-Term Peptide,true,-1; 10.0,C-Term Protein,true,-1",
            }
        )
        _, body = rows(build(inputs(workflow=wf)))
        assert body[0][16:19] == [
            "NT=TMTpro;MT=Variable;PP=Any N-term;MM=304.20709",
            "NT=TMTpro;MT=Fixed;PP=Any N-term;MM=304.20709",
            "NT=C-Term Protein:10;MT=Fixed;PP=Protein C-term;MM=10.00000",
        ]

    def test_mass_offsets_become_variable_mods(self, inputs):
        wf = workflow_text(
            **{
                "msfragger.table.var-mods": "",
                "msfragger.table.fix-mods": "",
                "msfragger.mass_offsets": "0/79.966331 15.5",
                "msfragger.restrict_deltamass_to": "STY",
            }
        )
        _, body = rows(build(inputs(workflow=wf)))
        assert body[0][16:18] == [
            "NT=offset:15.50000;MT=Variable;PP=Anywhere;TA=STY;MM=15.50000",
            "NT=Phospho;MT=Variable;PP=Anywhere;TA=STY;MM=79.96633",
        ]


class TestOrganism:
    def test_majority_organism(self):
        headers = [">x OS=Mus musculus OX=1", ">y OS=Homo sapiens OX=2", ">z OS=Homo sapiens GN=Q"]
        assert fs.guess_organism(headers) == "Homo sapiens"

    def test_no_uniprot_headers_is_not_available(self):
        assert fs.guess_organism([">contaminant_1", ">rev_xyz"]) == "not available"

    def test_tie_resolved_by_java_hashmap_order(self):
        # java.util.HashMap iterates "Mus musculus" before "Homo sapiens" (verified with JDK 17)
        headers = [">a OS=Homo sapiens OX=1", ">b OS=Mus musculus OX=2"]
        assert fs.guess_organism(headers) == "Mus musculus"

    def test_header_without_space_suffix_not_matched(self):
        assert fs.guess_organism([">a OS=Homo sapiens"]) == "not available"


class TestTolerances:
    def test_log_tolerance_used_only_with_parameter_optimization(self, inputs):
        wf = workflow_text(**{"msfragger.calibrate_mass": "1"})
        _, body = rows(build(inputs(workflow=wf)))
        assert body[0][14] == "20.0 ppm"

    def test_first_log_tolerance_wins(self, inputs):
        _, body = rows(build(inputs(logs=[LOG_B, LOG_A])))
        assert body[0][14] == "7.0 ppm"

    @pytest.mark.parametrize(
        "manifest, logs, expected",
        [
            # glyco-N-LFQ PXD030670: the manifest lists SALIVA-1 (10 ppm) first, MSFragger ran CANCER-3
            # (group cancer < saliva, 15 ppm) first
            (
                "/r/SALIVA-1.raw\tsaliva\t\tDDA\n/r/CANCER-3.raw\tcancer\t\tDDA",
                ["New fragment_mass_tolerance = 15.000000 PPM\n", "New fragment_mass_tolerance = 10.000000 PPM\n"],
                [("SALIVA-1.raw", "15.0 ppm"), ("CANCER-3.raw", "15.0 ppm")],
            ),
            # FPOP-coADAPTr-quant: the manifest lists A01 (15 ppm) first, MSFragger ran B01
            # (group Control_B01 < FPOP_A01, 10 ppm) first
            (
                "/r/run_A01.mzML\tFPOP_A01\t\tDDA\n/r/run_B01.mzML\tControl_B01\t\tDDA",
                ["New fragment_mass_tolerance = 10.000000 PPM\n", "New fragment_mass_tolerance = 15.000000 PPM\n"],
                [("run_A01.mzML", "10.0 ppm"), ("run_B01.mzML", "10.0 ppm")],
            ),
        ],
        ids=["glyco-N-LFQ", "FPOP-coADAPTr-quant"],
    )
    def test_two_files_take_the_first_logs_tolerance_not_the_first_manifest_files(self, inputs, manifest, logs, expected):
        # FragPipe's single matcher.find() (TabMsfragger.java:2359-2363) puts the first logged file's value on
        # every row; the workflow passes the logs in MSFragger's run order (msfraggerConsoleLogs)
        header, body = rows(build(inputs(manifest=manifest, logs=logs)))
        name, tol = header.index("comment[data file]"), header.index("comment[fragment mass tolerance]")
        assert [(r[name], r[tol]) for r in body] == expected

    def test_dalton_units_use_two_decimals(self, inputs):
        wf = workflow_text(
            **{
                "msfragger.precursor_mass_units": "0",
                "msfragger.precursor_mass_lower": "-0.5",
                "msfragger.precursor_mass_upper": "0.505",
                "msfragger.fragment_mass_units": "0",
                "msfragger.fragment_mass_tolerance": "0.025",
                "msfragger.calibrate_mass": "0",
            }
        )
        _, body = rows(build(inputs(workflow=wf)))
        assert body[0][13:15] == ["1.01 Da", "0.03 Da"]

    def test_instrument_map_basename_and_last_wins(self):
        text = LOG_A + "/other/a.mzML: Scans = 1; x = y; Instrument = Astral\n"
        assert fs.get_instrument_map(text) == {"a.mzML": "Astral"}


class TestHeaderTypes:
    def test_nonvertebrates_falls_through_like_fragpipe(self):
        h = fs.header_characteristics("NonVertebrates")
        assert h[2:12] == [
            "characteristics[strain/breed]",
            "characteristics[age]",
            "characteristics[developmental stage]",
            "characteristics[sex]",
            "characteristics[individual]",
            "characteristics[ecotype/cultivar]",
            "characteristics[age]",
            "characteristics[developmental stage]",
            "characteristics[individual]",
            "characteristics[cell line]",
        ]

    def test_human_type_from_workflow(self, inputs):
        header, body = rows(build(inputs(workflow=workflow_text(**{"workflow.misc.sdrf-type": "Human"}))))
        assert header[2] == "characteristics[ancestry category]" and len(body[0]) == len(header)

    def test_unknown_type_falls_back_to_default(self, inputs):
        header, _ = rows(build(inputs(workflow=workflow_text(**{"workflow.misc.sdrf-type": "Bogus"}))))
        assert header[2] == "characteristics[disease]"


class TestJavaEmulation:
    @pytest.mark.parametrize(
        "value,digits,expected",
        [(0.15, 1, "0.2"), (2.675, 2, "2.68"), (0.5, 0, "1"), (40.0, 1, "40.0")],
    )
    def test_java_fixed_half_up(self, value, digits, expected):
        assert fs.java_fixed(value, digits) == expected

    def test_float_mass_formatting(self):
        assert fs.java_fixed(fs.java_parse_float("57.02146"), 5) == "57.02146"

    @pytest.mark.parametrize("value,expected", [(5.0, "5.0"), (1e-4, "1.0E-4"), (1e7, "1.0E7"), (12.5, "12.5")])
    def test_double_to_string(self, value, expected):
        assert fs.java_double_to_string(value) == expected

    def test_properties_unescape(self, tmp_path):
        p = tmp_path / "w.workflow"
        p.write_text("a.b=C\\:\\\\x\\\\y\nkey value\\\n  continued\n! comment\n")
        assert fs.read_properties(str(p)) == {"a.b": "C:\\x\\y", "key": "valuecontinued"}

    def test_unimod_reader_names_and_sites(self, tmp_path):
        p = tmp_path / "u.obo"
        p.write_bytes(OBO.encode())
        db = {name: sites for _, name, sites in fs.read_unimod(str(p))}
        assert db["Acetyl"] == "KnS" and db["Label:2H(4)+Phospho"] == "S" and "unimod root node" in db


def test_missing_workflow_key_raises(inputs):
    # the input is FragPipe's saved UI state, which names every MSFragger key
    wf = "\n".join(line for line in workflow_text().splitlines() if "calibrate_mass" not in line)
    with pytest.raises(KeyError):
        build(inputs(workflow=wf))


def test_fragpipe_saved_state_of_a_workflow_without_msfragger_keys(inputs):
    # The msfragger.* keys FragPipe 24.0 saved for a workflow without any (closed-search state,
    # TabMsfragger.initMore loadDefaults(closed) :515-517): trypsin KR, 15.9949 M and 42.0106 [^
    # variable, 57.02146 C fixed, 20 ppm precursor, the values workflow_text() spells out
    saved = (FRAGGER_FIXTURES / "closed_defaults_saved" / "msfragger_keys.workflow").read_text()
    wf = "\n".join(line for line in workflow_text().splitlines() if not line.startswith("msfragger.")) + "\n" + saved
    assert build(inputs(workflow=wf)) == build(inputs())


def test_version_flag():
    out = subprocess.run([sys.executable, str(SCRIPT), "--version"], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "fragpipe_sdrf.py " + fs.__version__
