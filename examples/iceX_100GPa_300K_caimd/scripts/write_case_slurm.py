#!/usr/bin/env python3
"""Write manual 24 h deformation Slurm scripts from the GPU QE template."""

from __future__ import annotations

import re
from pathlib import Path


EXAMPLE_ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = EXAMPLE_ROOT / "scripts" / "run_qe_gpu.slurm"

JOBS = {
    "deformation/tension": "icex_tension_xx",
    "deformation/compression": "icex_compression_xx",
    "deformation/shear": "icex_shear_xy",
}


def specialize(template: str, job_name: str) -> str:
    text = re.sub(r"#SBATCH --time=.*", "#SBATCH --time=24:00:00", template)
    text = re.sub(r"#SBATCH --job-name=.*", f"#SBATCH --job-name={job_name}", text)
    text = text.replace("JOB DONE: Ice X QE case completed", f"JOB DONE: {job_name} completed")
    return text


def main() -> int:
    template = TEMPLATE.read_text(encoding="utf-8")
    for relative, job_name in JOBS.items():
        target = EXAMPLE_ROOT / relative / "run.slurm"
        target.write_text(specialize(template, job_name), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
