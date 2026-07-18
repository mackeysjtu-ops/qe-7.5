#!/usr/bin/env python3
"""Generate deterministic 300 K ATOMIC_VELOCITIES for the Ice X example."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import random
from pathlib import Path


RYDBERG_SI = 2.1798723611035e-18
K_BOLTZMANN_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
ELECTRONMASS_SI = 9.1093837015e-31
AMU_RY = (AMU_SI / ELECTRONMASS_SI) / 2.0
RY_TO_KELVIN = RYDBERG_SI / K_BOLTZMANN_SI

MASS_AMU = {
    "H": 1.0080,
    "O": 15.9994,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_qe_atom_labels(path: Path) -> list[str]:
    labels: list[str] = []
    in_positions = False
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        upper = stripped.upper()
        if upper.startswith("ATOMIC_POSITIONS"):
            in_positions = True
            continue
        if in_positions and (upper.startswith("K_POINTS") or upper.startswith("CELL_PARAMETERS") or upper.startswith("ATOMIC_VELOCITIES")):
            break
        if in_positions:
            fields = stripped.split()
            if fields and fields[0] in MASS_AMU:
                labels.append(fields[0])
    if not labels:
        raise ValueError(f"No ATOMIC_POSITIONS labels found in {path}")
    return labels


def kinetic_temperature(velocities: list[list[float]], masses_ry: list[float], ndof: int) -> tuple[float, float]:
    ekin = 0.0
    for velocity, mass in zip(velocities, masses_ry):
        ekin += 0.5 * mass * sum(component * component for component in velocity)
    temp = 2.0 * ekin / float(ndof) * RY_TO_KELVIN
    return ekin, temp


def center_of_mass_velocity(velocities: list[list[float]], masses_ry: list[float]) -> list[float]:
    total_mass = sum(masses_ry)
    return [
        sum(mass * velocity[axis] for velocity, mass in zip(velocities, masses_ry)) / total_mass
        for axis in range(3)
    ]


def generate(labels: list[str], target_temperature: float, seed: int) -> dict[str, object]:
    rng = random.Random(seed)
    masses_ry = [MASS_AMU[label] * AMU_RY for label in labels]
    velocities: list[list[float]] = []
    kt_ry = target_temperature / RY_TO_KELVIN
    for mass in masses_ry:
        sigma = math.sqrt(kt_ry / mass)
        velocities.append([rng.gauss(0.0, sigma) for _ in range(3)])

    initial_ekin, initial_temp = kinetic_temperature(velocities, masses_ry, 3 * len(labels))
    v_com = center_of_mass_velocity(velocities, masses_ry)
    for velocity in velocities:
        for axis in range(3):
            velocity[axis] -= v_com[axis]

    ndof = 3 * len(labels) - 3
    removed_ekin, removed_temp = kinetic_temperature(velocities, masses_ry, ndof)
    if removed_temp <= 0.0:
        raise ValueError("Generated zero kinetic temperature after COM removal")
    scale = math.sqrt(target_temperature / removed_temp)
    for velocity in velocities:
        for axis in range(3):
            velocity[axis] *= scale
    final_ekin, final_temp = kinetic_temperature(velocities, masses_ry, ndof)
    final_com = center_of_mass_velocity(velocities, masses_ry)
    speeds = [math.sqrt(sum(component * component for component in velocity)) for velocity in velocities]

    return {
        "labels": labels,
        "velocities_bohr_per_rydberg_au": velocities,
        "metadata": {
            "target_temperature_K": target_temperature,
            "seed": seed,
            "nat": len(labels),
            "ndof": ndof,
            "temperature_formula": "T = 2 * E_kin / ndof * ry_to_kelvin, with masses in amu_ry and velocities in Bohr/Ry atomic-time units",
            "initial_temperature_before_com_removal_K": initial_temp,
            "temperature_after_com_removal_before_rescale_K": removed_temp,
            "final_temperature_K": final_temp,
            "initial_ekin_ry": initial_ekin,
            "final_ekin_ry": final_ekin,
            "scale_after_com_removal": scale,
            "center_of_mass_velocity_before_removal_bohr_per_rydberg_au": v_com,
            "center_of_mass_velocity_after_rescale_bohr_per_rydberg_au": final_com,
            "max_abs_com_velocity_after_rescale": max(abs(component) for component in final_com),
            "rms_speed_bohr_per_rydberg_au": math.sqrt(sum(speed * speed for speed in speeds) / len(speeds)),
            "max_speed_bohr_per_rydberg_au": max(speeds),
            "masses_amu": MASS_AMU,
            "amu_ry": AMU_RY,
            "ry_to_kelvin": RY_TO_KELVIN,
            "python_version": platform.python_version(),
        },
    }


def write_velocity_include(path: Path, labels: list[str], velocities: list[list[float]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("ATOMIC_VELOCITIES { a.u }\n")
        for label, velocity in zip(labels, velocities):
            handle.write(
                f"{label:<2s} {velocity[0]:18.12e} {velocity[1]:18.12e} {velocity[2]:18.12e}\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structure",
        type=Path,
        default=Path("examples/iceX_100GPa_300K_caimd/structure/iceX_qe_structure.inc"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("examples/iceX_100GPa_300K_caimd/structure/iceX_velocities_300K.inc"),
    )
    parser.add_argument("--metadata", type=Path, default=None)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--seed", type=int, default=20260718)
    args = parser.parse_args()

    labels = read_qe_atom_labels(args.structure)
    generated = generate(labels, args.temperature, args.seed)
    velocities = generated["velocities_bohr_per_rydberg_au"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_velocity_include(args.output, labels, velocities)  # type: ignore[arg-type]
    metadata_path = args.metadata or args.output.with_suffix(".json")
    metadata = dict(generated["metadata"])  # type: ignore[arg-type]
    metadata["structure_include_sha256"] = sha256(args.structure)
    metadata["velocity_include_sha256"] = sha256(args.output)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"nat": len(labels), "ndof": metadata["ndof"], "final_temperature_K": metadata["final_temperature_K"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
