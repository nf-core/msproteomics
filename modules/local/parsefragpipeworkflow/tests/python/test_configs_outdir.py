"""Unit tests for --configs-outdir of parse_fragpipe_workflow.py.

One run with --output-json --outdir A --configs-outdir A/configs must write what two runs wrote before:
tool_configs.json (--output-json) and the config files (a plain run).
"""

import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "resources" / "usr" / "bin" / "parse_fragpipe_workflow.py"
# FragPipe 24.0 templates the module stages (the script takes them as arguments)
TEMPLATES = Path(__file__).resolve().parents[2] / "assets"
TEMPLATE_ARGS = ["--fragger-template", str(TEMPLATES / "fragger_closed.params"), "--diaumpire-template", str(TEMPLATES / "diaumpire_se.params")]
LFQ_MBR = Path(__file__).resolve().parents[5] / "assets" / "LFQ-MBR.workflow"


def run(*extra):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *TEMPLATE_ARGS, "--workflow", str(LFQ_MBR), "--data-types", "DDA", "--threads", "4", *extra],
        capture_output=True,
        text=True,
    )


def tree(root):
    return {str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()}


def test_one_run_writes_json_and_configs_of_two_runs(tmp_path):
    one = tmp_path / "one"
    result = run("--output-json", "--outdir", str(one), "--configs-outdir", str(one / "configs"))
    assert result.returncode == 0, result.stderr

    json_only = tmp_path / "json_only"
    assert run("--output-json", "--outdir", str(json_only)).returncode == 0
    configs_only = tmp_path / "configs_only"
    assert run("--outdir", str(configs_only)).returncode == 0

    written = tree(one)
    assert written["tool_configs.json"] == (json_only / "tool_configs.json").read_bytes()
    assert {k[len("configs/") :]: v for k, v in written.items() if k.startswith("configs/")} == tree(configs_only)
    assert len(tree(configs_only)) > 0


def test_configs_outdir_without_output_json_fails(tmp_path):
    result = run("--outdir", str(tmp_path / "out"), "--configs-outdir", str(tmp_path / "out" / "configs"))
    assert result.returncode != 0
    assert "--configs-outdir needs --output-json" in result.stderr


def test_script_runs_outside_the_module_tree(tmp_path):
    # In a Wave container the module's resources/usr/bin is /usr/bin: the script must not look for its
    # data next to itself, only at the templates it is given
    shallow = tmp_path / "parse_fragpipe_workflow.py"
    shallow.write_bytes(SCRIPT.read_bytes())
    out = tmp_path / "out"
    result = subprocess.run(
        [sys.executable, str(shallow), *TEMPLATE_ARGS, "--workflow", str(LFQ_MBR), "--data-types", "DDA", "--threads", "4",
         "--output-json", "--outdir", str(out), "--configs-outdir", str(out / "configs")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert (out / "tool_configs.json").is_file() and (out / "configs" / "fragger.params").is_file()
