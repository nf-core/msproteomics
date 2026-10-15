"""bin/philosopher_relocate_workspace.py rewrites the workspace paths in a msgpack met.Data record."""

import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin/philosopher_relocate_workspace.py"
spec = importlib.util.spec_from_file_location("relocate", SCRIPT)
relocate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(relocate)


def s(text):
    raw = text.encode()
    return relocate._str_header(len(raw)) + raw


def record(home, temp="/tmp/abc"):
    # map of 5: Home, MetaDir (str16-long path), Temp, UUID, Nested{List:[int, float]}
    long_dir = home + "/.meta/" + "x" * 300
    return (b"\x85" + s("Home") + s(home) + s("MetaDir") + s(long_dir) + s("Temp") + s(temp)
            + s("UUID") + s("u-1") + s("Nested") + b"\x81" + s("List") + b"\x92\x05\xcb" + b"\x40\x09\x21\xfb\x54\x44\x2d\x18")


class TestRelocate(unittest.TestCase):
    def test_paths_under_home_are_rewritten_others_kept(self):
        old, new = "/old/task/dir", "/new/place/g1"
        out = relocate.Relocator(record(old), old, new)
        self.assertEqual(out.element(0), len(record(old)))
        self.assertEqual(bytes(out.out), record(new))
        self.assertEqual(out.changed, 2)

    def test_home_and_temp_are_read_and_file_rewritten(self):
        # Filter's workspace: Temp is <filter task>/phi_temp/<uuid>, a sibling of Home, not under it
        old = record("/upstream/task/g1", temp="/upstream/task/phi_temp/u-1")
        with tempfile.TemporaryDirectory() as tmp:
            meta = Path(tmp, "g1", ".meta")
            meta.mkdir(parents=True)
            (meta / "meta.bin").write_bytes(old)
            self.assertEqual(relocate.field_of(old, "Home"), "/upstream/task/g1")
            self.assertEqual(relocate.field_of(old, "Temp"), "/upstream/task/phi_temp/u-1")
            import sys
            argv = sys.argv
            sys.argv = ["x", str(Path(tmp, "g1")), str(Path(tmp, "phi_temp", "u-2"))]
            try:
                relocate.main()
            finally:
                sys.argv = argv
            here = Path(tmp)  # as given: a symlinked work dir keeps its logical path
            self.assertEqual((meta / "meta.bin").read_bytes(), record(str(here / "g1"), temp=str(here / "phi_temp" / "u-2")))

    def test_relative_paths_fail(self):
        import sys
        argv = sys.argv
        sys.argv = ["x", "g1", "/abs/temp"]
        try:
            with self.assertRaises(SystemExit):
                relocate.main()
        finally:
            sys.argv = argv

    def test_empty_temp_fails(self):
        with self.assertRaises(ValueError):
            relocate.relocate(record("/old", temp=""), "", "/new/temp")

    def test_sibling_prefix_is_not_rewritten(self):
        out = relocate.Relocator(s("/old/task/dir2"), "/old/task/dir", "/new")
        out.element(0)
        self.assertEqual(out.changed, 0)


if __name__ == "__main__":
    unittest.main()
