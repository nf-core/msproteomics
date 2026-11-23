"""Static checks: every FragPipe route converts Thermo .raw to mzML with FRAGPIPE_RAW_TO_MZML.

FragPipe never gets a .raw: Thermo RawFileReader under Mono misreads a staged symlink
("No open raw file", 0 scans), and one conversion gives every route the same mzML.
Every workflow that calls a FragPipe search (FRAGPIPE, FRAGPIPE_HEADLESS_WF, MSFRAGGER) must give it the
`spectra` output of FRAGPIPE_RAW_TO_MZML, and build the manifest from those files, not from the
samplesheet files. workflows/tests/fragpipe_raw_routes.nf.test checks the run-time behaviour.
"""

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
WORKFLOWS = sorted((REPO / "workflows").glob("*.nf"))
SEARCHES = ("FRAGPIPE", "FRAGPIPE_HEADLESS_WF", "MSFRAGGER")
_INCLUDE = re.compile(r"\binclude\s*\{([^}]*)\}\s*from\s*'([^']+)'")


def included(text):
    """Return {name or alias: (process or workflow name, source)} for each include of a Nextflow file."""
    names = {}
    for block, source in _INCLUDE.findall(text):
        for item in block.split(";"):
            parts = item.split()
            if parts:
                names[parts[-1]] = (parts[0], source)
    return names


def first_arguments(text, name):
    """Return the first argument of each call `name(` (a channel variable name)."""
    return re.findall(rf"(?<![\w.]){name}\(\s*([\w.]+)", text)


def search_calls():
    """Yield (workflow file, alias) for each FragPipe search a workflow includes."""
    for path in WORKFLOWS:
        for alias, (name, source) in included(path.read_text()).items():
            if name in SEARCHES and source.startswith(("../subworkflows/local/", "../modules/local/")):
                yield path, alias


CALLS = list(search_calls())


def test_routes_found():
    found = {(p.name, alias) for p, alias in CALLS}
    assert found == {("fragpipe.nf", "FRAGPIPE_PIPELINE"), ("fragpipe.nf", "FRAGPIPE_HEADLESS_WF")}, (
        f"FragPipe search callers changed; every new route must call FRAGPIPE_RAW_TO_MZML first: {found}"
    )


@pytest.mark.parametrize("path,alias", CALLS, ids=[f"{p.name}:{a}" for p, a in CALLS])
def test_search_gets_the_converted_spectra(path, alias):
    text = path.read_text()
    assert included(text).get("FRAGPIPE_RAW_TO_MZML", ("", ""))[1] == "../subworkflows/local/fragpipe_raw_to_mzml/main", (
        f"{path.name} calls {alias} but does not include FRAGPIPE_RAW_TO_MZML"
    )
    converted = set(re.findall(r"(\w+)\s*=\s*FRAGPIPE_RAW_TO_MZML\.out\.spectra\b", text))
    assert converted, f"{path.name}: no variable takes FRAGPIPE_RAW_TO_MZML.out.spectra"
    args = first_arguments(text, alias)
    assert args and all(arg in converted for arg in args), f"{path.name}: {alias} gets {args}, not {sorted(converted)}"


def test_manifests_name_the_converted_files():
    """The samplesheet list (.raw paths) feeds the conversion only, never a manifest or a staged file list."""
    text = (REPO / "workflows" / "fragpipe.nf").read_text()
    assert "fragpipeManifestEntries(spectra_files" not in text
    assert re.findall(r"channel\.fromList\(spectra_files\)", text) == ["channel.fromList(spectra_files)"]
    assert "FRAGPIPE_RAW_TO_MZML(channel.fromList(spectra_files))" in text
    assert "msfraggerConsoleLogs(items, spectra_files)" not in text


def test_tmt_labelcheck_runs_the_fragpipe_workflow():
    """The TMT label check is MSPROTEOMICS_FRAGPIPE (so it converts .raw too), not its own FragPipe search."""
    text = (REPO / "workflows" / "tmt_labelcheck.nf").read_text()
    assert "MSPROTEOMICS_FRAGPIPE()" in text
    assert not [alias for alias, (name, _src) in included(text).items() if name in SEARCHES]
