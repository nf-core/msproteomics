"""Every subworkflows/local/*/meta.yml lists exactly the take: inputs and emit: outputs of its workflow.

The workflow is the one named like meta.yml `name` (upper-cased); main.nf may hold helper workflows too.
"""

from __future__ import annotations

import re
import unittest
from pathlib import Path

import yaml

SUBWORKFLOWS = Path(__file__).resolve().parents[2] / "subworkflows" / "local"


def workflow_take_emit(main_nf: str, name: str) -> tuple[list[str], list[str]]:
    """take: and emit: names of `workflow <name> { ... }` in main.nf text (the block ends at the next workflow)."""
    start = re.search(rf"^workflow\s+{name}\s*\{{", main_nf, re.M)
    if start is None:
        raise ValueError(f"no workflow {name}")
    nxt = re.search(r"^workflow\s+\w+\s*\{", main_nf[start.end():], re.M)
    body = main_nf[start.end(): start.end() + nxt.start()] if nxt else main_nf[start.end():]
    take = re.search(r"^\s*take:\s*$(.*?)^\s*main:\s*$", body, re.M | re.S)
    emit = re.search(r"^\s*emit:\s*$(.*)", body, re.M | re.S)
    take_names = re.findall(r"^\s*(\w+)\b(?!\s*=)", take.group(1), re.M) if take else []
    emit_names = re.findall(r"^\s*(\w+)\s*=", emit.group(1), re.M) if emit else []
    return take_names, emit_names


class TestSubworkflowMeta(unittest.TestCase):
    def test_meta_yml_matches_take_and_emit(self):
        dirs = sorted(d for d in SUBWORKFLOWS.iterdir() if (d / "meta.yml").exists())
        self.assertTrue(dirs)
        for d in dirs:
            with self.subTest(subworkflow=d.name):
                meta = yaml.safe_load((d / "meta.yml").read_text())
                take, emit = workflow_take_emit((d / "main.nf").read_text(), meta["name"].upper())
                self.assertEqual([next(iter(e)) for e in meta.get("input") or []], take, "meta.yml input vs take:")
                self.assertEqual([next(iter(e)) for e in meta.get("output") or []], emit, "meta.yml output vs emit:")

    def test_helper_workflow_is_not_read(self):
        text = "workflow A {\n    take:\n    x // c\n    main:\n    emit:\n    y = x\n}\nworkflow B {\n    take:\n    z\n    main:\n    emit:\n    w = z\n}\n"
        self.assertEqual(workflow_take_emit(text, "A"), (["x"], ["y"]))
        self.assertEqual(workflow_take_emit(text, "B"), (["z"], ["w"]))


if __name__ == "__main__":
    unittest.main()
