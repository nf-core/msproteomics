"""parse_fragpipe_workflow.load_java_properties reads a .workflow as java.util.Properties.load does
(FragPipe's PropsFile.load, FragPipe-GUI/src/main/java/org/nesvilab/fragpipe/api/PropsFile.java:55-63).

The expected maps were printed by Properties.load on OpenJDK 17 for the same text (2026-10-12).
"""

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "resources" / "usr" / "bin" / "parse_fragpipe_workflow.py"
spec = importlib.util.spec_from_file_location("parse_fragpipe_workflow", SCRIPT)
pfw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pfw)


def load(tmp_path, text):
    path = tmp_path / "t.workflow"
    path.write_text(text, encoding="utf-8")
    return pfw.load_java_properties(str(path))


def test_escapes_separators_continuations_and_comments_follow_java(tmp_path):
    text = (
        "a.esc=x\\=y\\:z\\\\w\\u0041\\q\n"
        "a.trail=v  \n"
        "  a.lead = v2\n"
        "a.multi=line1\\\n"
        "    line2\\\\\n"
        "a.next=after\n"
        "!bang=1\n"
        "# hash=1\n"
        "a.empty\n"
    )
    assert load(tmp_path, text) == {
        "a.empty": "",
        "a.esc": "x=y:z\\wAq",
        "a.lead": "v2",
        "a.multi": "line1line2\\",
        "a.next": "after",
        "a.trail": "v  ",
    }


def test_colon_and_whitespace_separate_key_and_value(tmp_path):
    assert load(tmp_path, "filler.colon:1\nfiller.space 2\nfiller.both = : 3\n") == {
        "filler.colon": "1",
        "filler.space": "2",
        "filler.both": ": 3",
    }


def test_repeated_key_keeps_the_last_value(tmp_path):
    assert load(tmp_path, "k=1\nk=2\n") == {"k": "2"}
