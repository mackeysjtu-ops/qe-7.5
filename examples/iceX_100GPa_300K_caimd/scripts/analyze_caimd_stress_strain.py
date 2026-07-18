#!/usr/bin/env python3
"""Analyze CAIMD stress-strain diagnostic files."""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path


EXAMPLE_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = EXAMPLE_ROOT / "convergence" / "analysis"
COLUMNS = [
    "ionic_step",
    "time_ps",
    "imposed_events_done",
    "following_move_is_event",
    "nominal_xx",
    "nominal_yy",
    "nominal_zz",
    "nominal_yz",
    "nominal_xz",
    "nominal_xy",
    "actual_xx",
    "actual_yy",
    "actual_zz",
    "actual_yz",
    "actual_xz",
    "actual_xy",
    "stress_xx_GPa",
    "stress_yy_GPa",
    "stress_zz_GPa",
    "stress_yz_GPa",
    "stress_xz_GPa",
    "stress_xy_GPa",
    "temperature_K",
    "volume_bohr3",
    "velocity_jump_max",
    "velocity_jump_rms",
]


def parse_caimd(path: Path) -> list[dict[str, float | int | bool]]:
    rows: list[dict[str, float | int | bool]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) < len(COLUMNS):
            continue
        row: dict[str, float | int | bool] = {}
        for name, value in zip(COLUMNS, fields):
            if name in {"ionic_step", "imposed_events_done"}:
                row[name] = int(value)
            elif name == "following_move_is_event":
                row[name] = value.upper().startswith("T")
            else:
                row[name] = float(value.replace("D", "E").replace("d", "e"))
        rows.append(row)
    return rows


def summarize(path: Path, rows: list[dict[str, float | int | bool]]) -> dict[str, object]:
    if not rows:
        return {"path": str(path.relative_to(EXAMPLE_ROOT)), "rows": 0, "status": "empty"}
    temps = [float(row["temperature_K"]) for row in rows]
    jumps = [float(row["velocity_jump_rms"]) for row in rows]
    final = rows[-1]
    first = rows[0]
    return {
        "path": str(path.relative_to(EXAMPLE_ROOT)),
        "rows": len(rows),
        "status": "ok",
        "first_step": first["ionic_step"],
        "last_step": final["ionic_step"],
        "events_done": final["imposed_events_done"],
        "initial_temperature_K": temps[0],
        "final_temperature_K": temps[-1],
        "min_temperature_K": min(temps),
        "max_temperature_K": max(temps),
        "mean_temperature_K": statistics.fmean(temps),
        "max_velocity_jump_rms": max(jumps),
        "final_nominal_xx": final["nominal_xx"],
        "final_nominal_yy": final["nominal_yy"],
        "final_nominal_zz": final["nominal_zz"],
        "final_nominal_yz": final["nominal_yz"],
        "final_nominal_xz": final["nominal_xz"],
        "final_nominal_xy": final["nominal_xy"],
        "final_actual_xx": final["actual_xx"],
        "final_actual_yy": final["actual_yy"],
        "final_actual_zz": final["actual_zz"],
        "final_actual_yz": final["actual_yz"],
        "final_actual_xz": final["actual_xz"],
        "final_actual_xy": final["actual_xy"],
        "final_stress_xx_GPa": final["stress_xx_GPa"],
        "final_stress_yy_GPa": final["stress_yy_GPa"],
        "final_stress_zz_GPa": final["stress_zz_GPa"],
        "final_stress_yz_GPa": final["stress_yz_GPa"],
        "final_stress_xz_GPa": final["stress_xz_GPa"],
        "final_stress_xy_GPa": final["stress_xy_GPa"],
    }


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    caimd_files = sorted(EXAMPLE_ROOT.glob("**/*.caimd"))
    summaries = [summarize(path, parse_caimd(path)) for path in caimd_files]
    json_path = OUTPUT_DIR / "caimd_stress_strain_summary.json"
    csv_path = OUTPUT_DIR / "caimd_stress_strain_summary.csv"
    md_path = OUTPUT_DIR / "caimd_stress_strain_summary.md"
    payload = {"files": summaries}
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    fieldnames = sorted({key for row in summaries for key in row.keys()}) if summaries else ["path", "rows", "status"]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in summaries:
            writer.writerow(row)
    with md_path.open("w", encoding="utf-8") as handle:
        handle.write("# CAIMD stress-strain summary\n\n")
        if not summaries:
            handle.write("No `.caimd` files found yet.\n")
        else:
            handle.write("| File | Rows | Events | Final nominal xx | Final nominal xy | Final T K |\n")
            handle.write("| --- | ---: | ---: | ---: | ---: | ---: |\n")
            for row in summaries:
                handle.write(
                    f"| {row['path']} | {row['rows']} | {row.get('events_done', '')} | "
                    f"{row.get('final_nominal_xx', '')} | {row.get('final_nominal_xy', '')} | "
                    f"{row.get('final_temperature_K', '')} |\n"
                )
    print(json.dumps({"caimd_files": len(caimd_files)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
