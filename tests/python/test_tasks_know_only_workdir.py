"""A task knows only its work dir: modules get files and folders as `path` inputs, never host paths.

Under Fusion an input is a link to the input bucket, so resolving it (readlink, realpath, canonical
paths) or writing next to or through it writes into the bucket. These checks scan the local module sources,
the Groovy helpers that render module scripts, the scripts modules run (modules/**/resources, bin/)
and the configs that could set a staging mode.
"""

import io
import re
import tokenize
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# Local modules only: nf-core modules are maintained upstream (nf-core lint and their own tests)
MODULES = sorted((ROOT / "modules/local").rglob("main.nf"))


def included_helpers():
    """The .nf files modules include (their functions render module script lines)."""
    found = set()
    for f in MODULES:
        for rel in re.findall(r"include\s*\{[^}]*\}\s*from\s*'([^']+)'", f.read_text()):
            p = (f.parent / rel).resolve()
            p = p if p.suffix == ".nf" else p.with_suffix(".nf")
            if p.is_file() and "modules" not in p.relative_to(ROOT).parts[:1]:
                found.add(p)
    return sorted(found)


HELPERS = included_helpers()
# Scripts a task runs: module resources and the bin/ scripts named in a module or helper
_TASK_TEXT = "\n".join(f.read_text() for f in MODULES + HELPERS)
SCRIPTS = sorted([p for p in (ROOT / "modules").rglob("resources/usr/bin/*") if p.suffix in (".py", ".sh", ".R", ".java")]
                 + [p for p in (ROOT / "bin").glob("*") if p.is_file() and p.name in _TASK_TEXT])


def strip_groovy_comments(text):
    """`text` without // and /* */ comments; string literals (', ", ''', \"\"\") are kept whole."""
    out, i, n = [], 0, len(text)
    while i < n:
        if text.startswith("//", i):
            j = text.find("\n", i)
            i = n if j < 0 else j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            i = n if j < 0 else j + 2
        elif text[i] in "'\"":
            q = text[i:i + 3] if text[i:i + 3] in ("'''", '"""') else text[i]
            j = i + len(q)
            while j < n and not text.startswith(q, j):
                j += 2 if text[j] == "\\" else 1
            out.append(text[i:j + len(q)])
            i = j + len(q)
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def strip_script_comments(path):
    text = path.read_text()
    if path.suffix == ".py":
        toks = [t for t in tokenize.generate_tokens(io.StringIO(text).readline) if t.type != tokenize.COMMENT]
        return tokenize.untokenize(toks)
    if path.suffix == ".java":
        return strip_groovy_comments(text)
    return "\n".join(re.sub(r"(^|\s)#.*$", "", line) for line in text.splitlines())


def hits(files, pattern, strip):
    rx = re.compile(pattern)
    return [f"{f.relative_to(ROOT)}: {m.group(0)}" for f in files for m in rx.finditer(strip(f))]


# Link resolution and locations outside the task, in module code and the helpers it includes
NF_FORBIDDEN = r"\b(readlink|realpath|toRealPath|getCanonicalPath|getCanonicalFile|projectDir|launchDir)\b(?<!--projectDir)|\bnew\s+File\s*\("
# The same in the scripts tasks run. Path(__file__).resolve() locates the script's own data files.
SCRIPT_FORBIDDEN = (r"\b(readlink|realpath|toRealPath|getCanonicalPath|getCanonicalFile)\b"
                    r"|(?<!__file__\))\.resolve\(\s*\)")

# `val` inputs whose name looks like a path. Each is text the module writes into a tool file or
# argv, never a location the task opens.
VAL_PATH_NAMES = {
    "root_group_dir": "name of the staged group dir holding FragPipe's blank group",
    "pin_dirs": "FragPipe group names of the PINs",
    "filelist": "spectra file names for filelist_opair.txt",
    "pepxml_paths": "pepXML locations relative to the FragPipe work dir",
    "filelist_paths": "pepXML locations relative to the FragPipe work dir",
    "ptmprophet_cli": "PTMProphet options",
    "workdir": "FragPipe's workdir value recorded in the saved fragpipe.workflow; FragpipeUiState only writes it",
    "database_source": "FRAGPIPEPLAN: the database's original path (s3:// kept), text for the portable copies' source map",
}
PATH_LIKE = re.compile(r"path|dir|file|fasta|database|db|folder|home|tmp|temp", re.I)

# LC-MS path inputs (names) of local modules; each such module replaces staged .raw links by copies
# (rawLinkCopies: batmass-io under Mono reads a symlinked .raw as an empty file), unless listed here.
LCMS_INPUT = re.compile(r"path\s*\(?\s*(\w*(?:spectra|lcms|mzml|spec_files|dia_file|raw|ms_file)\w*)", re.I)
NO_RAW_READ = {
    "fragpipeheadless": "headlessStagingLines calls rawLinkCopies",
    "fragpipeplan": "a dry run reads no spectra: FragPipe only checks that the manifest files exist",
    "fragpipesaveworkflow": "FragpipeUiState only checks that the manifest files exist",
    "freequant": "FreeQuant takes only mzML (CmdFreequant.java:44)",
    "philosopherlabelquant": "LabelQuant takes only mzML (CmdLabelquant.java:51)",
    "msfraggerapplycalibration": "rewrites params; the spectra are mzML/mzBIN names",
    "mzmlindexing": "OpenMS on mzML",
    "openmspeakpicker": "OpenMS on mzML",
    "tdf2mzml": "Bruker .d",
    "percolator": "PercolatorOutputToPepXML gets _uncalibrated.mzML for .raw (CmdPercolator)",
    "philosopherpeptideprophet": "RewritePepxml only names the LC-MS files",
    "ptmprophet": "PTMProphet does not read the LC-MS files",
}
# Tools that write through a staged input: the module replaces exactly those links by copies (linkCopies)
WRITE_THROUGH = {
    "msbooster": ("msboosterPinStaging", "writes <pin>_edited.pin next to the canonical PIN path"),
    "philosopherpeptideprophet": ("linkCopiesWhenPathContains(pepxml_files", "rewrites in place an input pepXML whose path contains 'interact'"),
}
# stageInMode settings allowed in modules and configs (DECOMPRESS: retry staging modes)
STAGE_IN_MODE = {"modules/local/decompress/main.nf"}
# test configs too: a golden or integration config that copy-stages hides a write through a staged link
CONFIGS = sorted([f for d in ("conf", "modules", "subworkflows", "tests") for f in (ROOT / d).rglob("*.config")]
                 + [ROOT / "nextflow.config"])


SECTION = re.compile(r"^\s*(input|output|when|script|shell|exec|stub):\s*$", re.M)


def sections(text):
    """{label: body} of a process: each body runs to the next section label (or the end of the file)."""
    marks = list(SECTION.finditer(text))
    return {m.group(1): text[m.end():marks[i + 1].start() if i + 1 < len(marks) else len(text)]
            for i, m in enumerate(marks)}


def module_sections(f):
    found = sections(strip_groovy_comments(f.read_text()))
    if "input" not in found or "script" not in found:
        raise AssertionError(f"{f.relative_to(ROOT)}: no input: or script: section found; the checks would skip it")
    return found


def module_key(f):
    return str(f.parent.relative_to(ROOT / "modules/local"))


class TestTasksKnowOnlyWorkdir(unittest.TestCase):
    def test_scope_covers_helpers_and_task_scripts(self):
        names = {p.name for p in HELPERS + SCRIPTS}
        for expected in ("fragpipe_utils_a.nf", "fragpipe_utils_b.nf", "philosopher_relocate_workspace.py",
                         "merge_split_search.py"):
            self.assertIn(expected, names)

    def test_every_local_module_parses(self):
        for f in MODULES:
            module_sections(f)

    def test_sections_end_at_the_next_label(self):
        found = sections("process P {\n    input:\n    path x\n\n    script:\n    \"\"\"\n    cat x\n    \"\"\"\n}\n")
        self.assertEqual(found["input"].split(), ["path", "x"])
        self.assertIn("cat x", found["script"])

    def test_comment_strip_keeps_strings(self):
        self.assertEqual(strip_groovy_comments("a // realpath\nb"), "a \nb")
        self.assertEqual(strip_groovy_comments("x = 'http://h/realpath' /* c */"), "x = 'http://h/realpath' ")

    def test_no_link_resolution_in_modules_or_script_helpers(self):
        found = hits(MODULES + HELPERS, NF_FORBIDDEN, lambda f: strip_groovy_comments(f.read_text()))
        self.assertEqual(found, [], "a task must not resolve its staged links or name locations outside it")

    def test_no_link_resolution_in_task_scripts(self):
        found = hits(SCRIPTS, SCRIPT_FORBIDDEN, strip_script_comments)
        self.assertEqual(found, [], "scripts a task runs must keep paths as given")

    def test_no_params_in_modules(self):
        found = hits(MODULES, r"(^|[^A-Za-z_])params\.[a-z]", lambda f: strip_groovy_comments(f.read_text()))
        self.assertEqual(found, [], "modules take paths as inputs, not params")

    def test_path_like_val_inputs_are_listed(self):
        found = {}
        for f in MODULES:
            for name in re.findall(r"\bval\s*\(?\s*([A-Za-z_]\w*)", module_sections(f)["input"]):
                if not name.startswith("meta") and PATH_LIKE.search(name):
                    found.setdefault(name, []).append(str(f.relative_to(ROOT)))
        unlisted = {n: fs for n, fs in found.items() if n not in VAL_PATH_NAMES}
        self.assertEqual(unlisted, {}, "files and folders are path inputs; a path-like val needs a reason in VAL_PATH_NAMES")
        self.assertEqual(sorted(set(VAL_PATH_NAMES) - set(found)), [], "stale VAL_PATH_NAMES entries")

    def test_lcms_inputs_get_raw_copies(self):
        missing, seen = [], set()
        for f in MODULES:
            found = module_sections(f)
            if LCMS_INPUT.search(found["input"]):
                seen.add(module_key(f))
                # the copy lines must be in the task script (directly or through a helper it calls there)
                if module_key(f) not in NO_RAW_READ and "rawLinkCopies(" not in found["script"]:
                    missing.append(module_key(f))
        self.assertEqual(missing, [], "a module with LC-MS inputs replaces staged .raw links by copies (rawLinkCopies)")
        self.assertEqual(sorted(set(NO_RAW_READ) - seen), [], "stale NO_RAW_READ entries")

    def test_write_through_inputs_are_copies(self):
        missing = [m for m, (call, _) in WRITE_THROUGH.items()
                   if call not in module_sections(ROOT / "modules/local" / m / "main.nf")["script"]]
        self.assertEqual(missing, [], "inputs a tool writes through must be copies")

    def test_stage_in_mode_only_where_listed(self):
        found = [str(f.relative_to(ROOT)) for f in MODULES + CONFIGS
                 if re.search(r"\bstageInMode\b", strip_groovy_comments(f.read_text()))]
        self.assertEqual(sorted(set(found) - STAGE_IN_MODE), [],
                         "copy only the inputs a tool misreads or writes through (linkCopies), not every input")


if __name__ == "__main__":
    unittest.main()
