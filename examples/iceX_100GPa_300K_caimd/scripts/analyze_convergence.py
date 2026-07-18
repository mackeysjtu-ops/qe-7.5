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
RY_BOHR_TO_EV_ANG = 25.71104309541616


def parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def parse_qe_output(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"status": "missing", "job_done": False}
    text = path.read_text(encoding="utf-8", errors="replace")
    energies = [parse_float(match) for match in re.findall(r"!\s+total energy\s+=\s+([-+0-9.eEdD]+)\s+Ry", text)]
    pressures = [parse_float(match) for match in re.findall(r"P=\s*([-+0-9.eEdD]+)", text)]
    total_forces = [parse_float(match) for match in re.findall(r"Total force\s+=\s*([-+0-9.eEdD]+)", text)]
    atomic_forces = [
        math.sqrt(sum(parse_float(value) ** 2 for value in match))
        for match in re.findall(
            r"force\s+=\s*([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)",
            text,
            flags=re.IGNORECASE,
        )
    ]
    temperatures = [parse_float(match) for match in re.findall(r"temperature\s+=\s*([-+0-9.eEdD]+)\s+K", text, flags=re.IGNORECASE)]
    temperatures.extend(parse_float(match) for match in re.findall(r"\bT\s*=\s*([-+0-9.eEdD]+)\s+K", text))
    scf_iterations = [int(match) for match in re.findall(r"convergence has been achieved in\s+([0-9]+)\s+iterations", text)]
    scf_not_converged = bool(re.search(r"convergence\s+NOT\s+achieved", text, flags=re.IGNORECASE))
    wall_times = re.findall(r"PWSCF\s+:\s+(.+?)CPU\s+(.+?)WALL", text)
    stress_tensors_kbar = parse_stress_tensors_kbar(text)
    return {
        "status": "ok" if "JOB DONE." in text and not scf_not_converged else "incomplete",
        "job_done": "JOB DONE." in text,
        "scf_not_converged": scf_not_converged,
        "final_total_energy_Ry": energies[-1] if energies else None,
        "final_pressure_kbar": pressures[-1] if pressures else None,
        "final_pressure_GPa": pressures[-1] / 10.0 if pressures else None,
        "final_total_force_Ry_bohr": total_forces[-1] if total_forces else None,
        "final_total_force_eV_A": total_forces[-1] * RY_BOHR_TO_EV_ANG if total_forces else None,
        "max_atomic_force_Ry_bohr": max(atomic_forces) if atomic_forces else None,
        "max_atomic_force_eV_A": max(atomic_forces) * RY_BOHR_TO_EV_ANG if atomic_forces else None,
        "final_stress_kbar": stress_tensors_kbar[-1] if stress_tensors_kbar else None,
        "final_stress_GPa": [[value / 10.0 for value in row] for row in stress_tensors_kbar[-1]] if stress_tensors_kbar else None,
        "min_temperature_K": min(temperatures) if temperatures else None,
        "max_temperature_K": max(temperatures) if temperatures else None,
        "last_temperature_K": temperatures[-1] if temperatures else None,
        "scf_iterations_last": scf_iterations[-1] if scf_iterations else None,
        "wall_time_text": wall_times[-1][1].strip() if wall_times else None,
    }


def parse_stress_tensors_kbar(text: str) -> list[list[list[float]]]:
    lines = text.splitlines()
    tensors: list[list[list[float]]] = []
    for index, line in enumerate(lines):
        if "total   stress" not in line or "(kbar)" not in line:
            continue
        rows: list[list[float]] = []
        for candidate in lines[index + 1 : index + 4]:
            values = re.findall(r"[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eEdD][-+]?[0-9]+)?", candidate)
            if len(values) >= 6:
                rows.append([parse_float(value) for value in values[-3:]])
        if len(rows) == 3:
            tensors.append(rows)
    return tensors


def parse_gpu_log(case_dir: Path) -> dict[str, object]:
    logs = sorted(case_dir.glob("nvidia-smi-*.log"))
    gpu_utils: list[float] = []
    mem_used: list[float] = []
    power: list[float] = []
    for log in logs:
        for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
            fields = [field.strip() for field in line.split(",")]
            if len(fields) < 7:
                continue
            try:
                gpu_utils.append(float(fields[3]))
                mem_used.append(float(fields[5]))
                power.append(float(fields[6]))
            except ValueError:
                continue
    return {
        "gpu_log_count": len(logs),
        "max_gpu_utilization_pct": max(gpu_utils) if gpu_utils else None,
        "max_gpu_memory_MiB": max(mem_used) if mem_used else None,
        "max_gpu_power_W": max(power) if power else None,
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
        gpu = parse_gpu_log(metadata_path.parent)
        case_dir = str(metadata["case_dir"])
        if parsed["status"] == "missing" and case_dir.startswith("deformation/"):
            parsed["status"] = "not_submitted"
        if parsed["status"] == "missing" and case_dir == "reference/optional_vc_relax":
            parsed["status"] = "optional_not_run"
        cases.append({**metadata, **parsed, **gpu})
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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
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
