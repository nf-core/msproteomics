#!/usr/bin/env python3
"""FRAGPIPEPLAN inputs of every case in cases.json: plans/<case>.plan_inputs.json for tests/lib/record_fragpipe_plans.nf.test.

The same run as capture.sh: two DDA files in groups exp1_1 and exp2_1, a TMT workflow with each group's
annotation (the first channel_num channels, the 126 channel the pool), FragPipe --threads 4.

    python3 -I tests/fixtures/ptmshepherd_parity/generate_specs.py
"""

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from capture import case_workflow, cases  # noqa: E402

CHANNELS = "126 127N 127C 128N 128C 129N 129C 130N 130C 131N 131C 132N 132C 133N 133C 134N 134C 135N".split()


def spec(workflow: str) -> dict:
    rows = [["exp1.mzML", "exp1", 1, "DDA", "exp1_1"], ["exp2.mzML", "exp2", 1, "DDA", "exp2_1"]]
    tmt = re.search(r"^tmtintegrator.run-tmtintegrator=true$", workflow, re.M)
    plex = int(re.search(r"^tmtintegrator.channel_num=TMT-(\d+)$", workflow, re.M)[1]) if tmt else 0
    annotation = [[group, c, f"pool_{exp}" if c == "126" else f"{exp}_{c}"]
                  for exp, group in (("exp1", "exp1_1"), ("exp2", "exp2_1")) for c in CHANNELS[:plex]]
    return {"workflow_text": workflow, "manifest_rows": rows, "annotation_rows": annotation,
            "sources": {r[0]: r[0] for r in rows}, "database_name": "db.fasta", "database_source": "db.fasta",
            "ext_args": "--threads 4"}


def main() -> None:
    out = HERE / "plans"
    out.mkdir(exist_ok=True)
    for name, case in cases().items():
        (out / f"{name}.plan_inputs.json").write_text(json.dumps(spec(case_workflow(case)), indent=2) + "\n")


if __name__ == "__main__":
    main()
