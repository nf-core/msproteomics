#!/usr/bin/env python3
"""Copy a FRAGPIPEPLAN work dir (fp_work) to a fixture: the task dir as /plan, the run log licence-redacted
(tests/lib/record_fragpipe_plans.nf.test records the stub tests' dry runs).

Usage: normalise_fragpipe_plan.py <fp_work> <fixture dir>
The task dir is the parent of the output directory in FragPipe's fragpipe.job (column 4). The fixture
must name no licence (fragpipe_log_to_manifest.contains_license_material), as plan_to_tool_configs.py requires.
"""

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
import fragpipe_log_to_manifest as fpm  # noqa: E402


def main(src: Path, dst: Path) -> None:
    job = [line for line in (src / "fragpipe.job").read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]
    task_dir = job[0].split("\t")[3].rstrip("/").rsplit("/", 1)[0]
    if dst.exists():
        shutil.rmtree(dst)
    for path in sorted(p for p in src.rglob("*") if p.is_file()):
        out = dst / path.relative_to(src)
        out.parent.mkdir(parents=True, exist_ok=True)
        text = path.read_text(encoding="utf-8").replace(task_dir, "/plan")
        if path.parent == src and path.name.startswith("log_") and path.suffix == ".txt":
            text = "\n".join(fpm.redact_lines(text.split("\n")))
        if fpm.contains_license_material(text):
            sys.exit(f"{out} names the licence")
        out.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
