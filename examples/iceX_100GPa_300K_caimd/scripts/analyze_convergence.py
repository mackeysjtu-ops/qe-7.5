#!/usr/bin/env python3
"""Summarize bounded SCF/MD outputs for the Ice X example."""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path


EXAMPLE_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = EXAMPLE_ROOT / "convergence" / "analysis"


def parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def parse_qe_output(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"status": "missing", "job_done": False}
    text = path.read_text(encoding="utf-8", errors="replace")
    energies = [parse_float(match) for match in re.findall(r"!\s+total energy\s+=\s+([-+0-9.eEdD]+)\s+Ry", text)]
    pressures = [parse_float(match) for match in re.findall(r"P=\s*([-+0-9.eEdD]+)", text)]
    total_forces = [parse_float(match) for match in re.findall(r"Total force\s+=\s*([-+0-9.eEdD]+)", text)]
    temperatures = [parse_float(match) for match in re.findall(r"temperature\s+=\s*([-+0-9.eEdD]+)\s+K", text, flags=re.IGNORECASE)]
    scf_not_converged = bool(re.search(r"convergence\s+NOT\s+achieved", text, flags=re.IGNORECASE))
    wall_times = re.findall(r"PWSCF\s+:\s+(.+?)CPU\s+(.+?)WALL", text)
    return {
        "status": "ok" if "JOB DONE." in text and not scf_not_converged else "incomplete",
        "job_done": "JOB DONE." in text,
        "scf_not_converged": scf_not_converged,
        "final_total_energy_Ry": energies[-1] if energies else None,
        "final_pressure_kbar": pressures[-1] if pressures else None,
        "final_pressure_GPa": pressures[-1] / 10.0 if pressures else None,
        "final_total_force_Ry_bohr": total_forces[-1] if total_forces else None,
        "min_temperature_K": min(temperatures) if temperatures else None,
        "max_temperature_K": max(temperatures) if temperatures else None,
        "last_temperature_K": temperatures[-1] if temperatures else None,
        "wall_time_text": wall_times[-1][1].strip() if wall_times else None,
    }


def read_case_metadata(path: Path) -> dict[str, object]:
    metadata = json.loads(path.read_text(encoding="utf-8"))
    metadata["case_dir"] = str(path.parent.relative_to(EXAMPLE_ROOT))
    return metadata


def collect_cases() -> list[dict[str, object]]:
    cases: list[dict[str, object]] = []
    for metadata_path in sorted(EXAMPLE_ROOT.glob("**/CASE_METADATA.json")):
        metadata = read_case_metadata(metadata_path)
        parsed = parse_qe_output(metadata_path.parent / "qe.out")
        case_dir = str(metadata["case_dir"])
        if parsed["status"] == "missing" and case_dir.startswith("deformation/"):
            parsed["status"] = "not_submitted"
        if parsed["status"] == "missing" and case_dir == "reference/optional_vc_relax":
            parsed["status"] = "optional_not_run"
        cases.append({**metadata, **parsed})
    return cases


def energy_delta_mev_per_atom(reference: float | None, value: float | None, nat: int = 144) -> float | None:
    if reference is None or value is None:
        return None
    return (value - reference) * 13605.693122994 / nat


def annotate_convergence(cases: list[dict[str, object]]) -> dict[str, object]:
    ecut_cases = [case for case in cases if str(case["case_dir"]).startswith("convergence/ecut/")]
    k_cases = [case for case in cases if str(case["case_dir"]).startswith("convergence/kpoints/")]
    completed_ecut = sorted(
        [case for case in ecut_cases if case.get("job_done")],
        key=lambda case: float(case["ecutwfc_Ry"]),
    )
    completed_k = sorted(
        [case for case in k_cases if case.get("job_done")],
        key=lambda case: int(case["kgrid"][0]) * int(case["kgrid"][1]) * int(case["kgrid"][2]),
    )
    ecut_reference = completed_ecut[-1]["final_total_energy_Ry"] if completed_ecut else None
    k_reference = completed_k[-1]["final_total_energy_Ry"] if completed_k else None
    for case in ecut_cases:
        case["delta_to_highest_ecut_meV_per_atom"] = energy_delta_mev_per_atom(
            ecut_reference if isinstance(ecut_reference, float) else None,
            case.get("final_total_energy_Ry") if isinstance(case.get("final_total_energy_Ry"), float) else None,
        )
    for case in k_cases:
        case["delta_to_highest_k_meV_per_atom"] = energy_delta_mev_per_atom(
            k_reference if isinstance(k_reference, float) else None,
            case.get("final_total_energy_Ry") if isinstance(case.get("final_total_energy_Ry"), float) else None,
        )
    return {
        "ecut_reference_energy_Ry": ecut_reference,
        "kpoint_reference_energy_Ry": k_reference,
        "completed_cases": sum(1 for case in cases if case.get("job_done")),
        "pending_validation_cases": [
            case["case_dir"]
            for case in cases
            if case.get("status") in {"missing", "incomplete"}
        ],
        "final_not_submitted_cases": [
            case["case_dir"] for case in cases if case.get("status") == "not_submitted"
        ],
        "optional_not_run_cases": [
            case["case_dir"] for case in cases if case.get("status") == "optional_not_run"
        ],
    }


def write_outputs(cases: list[dict[str, object]], summary: dict[str, object]) -> None:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    json_path = ANALYSIS_DIR / "convergence_summary.json"
    csv_path = ANALYSIS_DIR / "convergence_summary.csv"
    md_path = ANALYSIS_DIR / "convergence_summary.md"
    payload = {"summary": summary, "cases": cases}
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    fieldnames = sorted({key for case in cases for key in case.keys()})
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for case in cases:
            writer.writerow(case)

    with md_path.open("w", encoding="utf-8") as handle:
        handle.write("# Ice X bounded convergence summary\n\n")
        handle.write(f"- Completed cases: {summary['completed_cases']}\n")
        handle.write(f"- Pending validation cases: {summary['pending_validation_cases']}\n")
        handle.write(f"- Final not-submitted cases: {summary['final_not_submitted_cases']}\n")
        handle.write(f"- Optional not-run cases: {summary['optional_not_run_cases']}\n\n")
        handle.write("| Case | Status | Ecut Ry | K-grid | Energy Ry | Pressure GPa | Force Ry/bohr |\n")
        handle.write("| --- | --- | ---: | --- | ---: | ---: | ---: |\n")
        for case in cases:
            energy = case.get("final_total_energy_Ry")
            pressure = case.get("final_pressure_GPa")
            force = case.get("final_total_force_Ry_bohr")
            kgrid = case.get("kgrid")
            handle.write(
                f"| {case['case_dir']} | {case.get('status')} | {case.get('ecutwfc_Ry')} | {kgrid} | "
                f"{energy if energy is not None else ''} | {pressure if pressure is not None else ''} | "
                f"{force if force is not None else ''} |\n"
            )
        if math.isnan(float("nan")):
            handle.write("\n")


def main() -> int:
    cases = collect_cases()
    summary = annotate_convergence(cases)
    write_outputs(cases, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
