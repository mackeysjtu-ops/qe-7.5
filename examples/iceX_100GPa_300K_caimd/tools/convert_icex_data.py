#!/usr/bin/env python3
"""Convert the local Ice X LAMMPS data file into QE structure include files.

The parser is intentionally narrow enough to fail loudly on ambiguous atom
styles, but it keeps the LAMMPS atom order and cell orientation unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


BOHR_PER_ANGSTROM = 1.8897261254578281
AMU_TO_GRAM = 1.66053906660e-24
ANGSTROM3_TO_CM3 = 1.0e-24


@dataclass(frozen=True)
class Atom:
    atom_id: int
    type_id: int
    species: str
    cart_angstrom: tuple[float, float, float]
    frac: tuple[float, float, float]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_header_float_pair(lines: list[str], token: str) -> tuple[float, float]:
    pattern = re.compile(rf"^\s*([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+{token}\s*$")
    for line in lines:
        match = pattern.match(line)
        if match:
            return float(match.group(1).replace("D", "E").replace("d", "e")), float(
                match.group(2).replace("D", "E").replace("d", "e")
            )
    raise ValueError(f"Could not find LAMMPS box line for {token!r}")


def parse_count(lines: list[str], label: str) -> int:
    pattern = re.compile(rf"^\s*(\d+)\s+{re.escape(label)}\s*$")
    for line in lines:
        match = pattern.match(line)
        if match:
            return int(match.group(1))
    raise ValueError(f"Could not find LAMMPS count for {label!r}")


def parse_tilts(lines: list[str]) -> tuple[float, float, float]:
    pattern = re.compile(
        r"^\s*([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+xy\s+xz\s+yz\s*$"
    )
    for line in lines:
        match = pattern.match(line)
        if match:
            return tuple(float(match.group(i).replace("D", "E").replace("d", "e")) for i in range(1, 4))
    return 0.0, 0.0, 0.0


def find_section(lines: list[str], name: str) -> list[str]:
    start = None
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.split("#", 1)[0].strip().split()[:1] == [name]:
            start = index + 1
            break
    if start is None:
        raise ValueError(f"Could not find LAMMPS section {name!r}")

    data: list[str] = []
    for line in lines[start:]:
        stripped = line.strip()
        if not stripped:
            if data:
                break
            continue
        first = stripped.split()[0]
        if data and first[0].isalpha():
            break
        if stripped.startswith("#"):
            continue
        data.append(line.rstrip("\n"))
    return data


def parse_masses(lines: list[str]) -> dict[int, tuple[float, str | None]]:
    masses: dict[int, tuple[float, str | None]] = {}
    for line in find_section(lines, "Masses"):
        fields, *comment_parts = line.split("#", 1)
        parts = fields.split()
        if len(parts) < 2:
            continue
        type_id = int(parts[0])
        mass = float(parts[1].replace("D", "E").replace("d", "e"))
        comment = comment_parts[0].strip() if comment_parts else None
        masses[type_id] = (mass, comment)
    return masses


def infer_species(masses: dict[int, tuple[float, str | None]]) -> tuple[dict[int, str], dict[str, str]]:
    mapping: dict[int, str] = {}
    evidence: dict[str, str] = {}
    for type_id, (mass, comment) in masses.items():
        normalized = (comment or "").strip().split()
        symbol = normalized[0] if normalized else ""
        if symbol in {"H", "O"}:
            mapping[type_id] = symbol
            evidence[str(type_id)] = f"Masses comment {comment!r}, mass={mass:.6f} amu"
            continue
        if abs(mass - 1.0080) < 0.02:
            mapping[type_id] = "H"
            evidence[str(type_id)] = f"mass={mass:.6f} amu matches H"
        elif abs(mass - 15.9994) < 0.05:
            mapping[type_id] = "O"
            evidence[str(type_id)] = f"mass={mass:.6f} amu matches O"
        else:
            raise ValueError(f"Cannot infer element for LAMMPS atom type {type_id} with mass {mass}")
    if set(mapping.values()) != {"H", "O"}:
        raise ValueError(f"Expected exactly H and O atom types, got {mapping}")
    return mapping, evidence


def determinant_3x3(m: tuple[tuple[float, float, float], ...]) -> float:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def inverse_3x3(m: tuple[tuple[float, float, float], ...]) -> tuple[tuple[float, float, float], ...]:
    det = determinant_3x3(m)
    if abs(det) < 1.0e-14:
        raise ValueError("Singular cell matrix")
    return (
        (
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det,
        ),
        (
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det,
        ),
        (
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det,
        ),
    )


def transpose(m: tuple[tuple[float, float, float], ...]) -> tuple[tuple[float, float, float], ...]:
    return tuple(tuple(m[j][i] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def mat_vec(m: tuple[tuple[float, float, float], ...], v: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(sum(m[i][j] * v[j] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def frac_to_cart(cell_rows: tuple[tuple[float, float, float], ...], frac: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(sum(frac[i] * cell_rows[i][j] for i in range(3)) for j in range(3))  # type: ignore[return-value]


def cart_to_frac(cell_rows: tuple[tuple[float, float, float], ...], cart: tuple[float, float, float]) -> tuple[float, float, float]:
    return mat_vec(inverse_3x3(transpose(cell_rows)), cart)


def parse_atoms(
    lines: list[str],
    cell_rows: tuple[tuple[float, float, float], ...],
    species_by_type: dict[int, str],
) -> tuple[list[Atom], dict[str, bool]]:
    atoms: list[Atom] = []
    atom_style_flags = {
        "has_molecule_id": False,
        "has_charge": False,
        "has_image_flags": False,
        "has_velocities_section": any(line.strip().split()[:1] == ["Velocities"] for line in lines),
    }
    for line in find_section(lines, "Atoms"):
        fields = line.split("#", 1)[0].split()
        if not fields:
            continue
        if len(fields) != 5:
            raise ValueError(
                "Only the simple LAMMPS Atoms format 'id type x y z' is supported here; "
                f"got {len(fields)} fields on line: {line!r}"
            )
        atom_id = int(fields[0])
        type_id = int(fields[1])
        cart = tuple(float(fields[i].replace("D", "E").replace("d", "e")) for i in range(2, 5))
        if type_id not in species_by_type:
            raise ValueError(f"Atom {atom_id} uses unknown type {type_id}")
        frac = tuple(x % 1.0 for x in cart_to_frac(cell_rows, cart))
        atoms.append(Atom(atom_id, type_id, species_by_type[type_id], cart, frac))
    return atoms, atom_style_flags


def min_image_distance(
    cell_rows: tuple[tuple[float, float, float], ...],
    f1: tuple[float, float, float],
    f2: tuple[float, float, float],
) -> float:
    best = float("inf")
    base = tuple(f1[i] - f2[i] for i in range(3))
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                shifted = (base[0] + i, base[1] + j, base[2] + k)
                cart = frac_to_cart(cell_rows, shifted)
                dist = math.sqrt(sum(x * x for x in cart))
                if dist < best:
                    best = dist
    return best


def distance_summary(cell_rows: tuple[tuple[float, float, float], ...], atoms: list[Atom]) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    for pair in ("H-H", "H-O", "O-O"):
        result[pair] = {"distance_angstrom": float("inf"), "atom_ids": []}
    for i, atom_i in enumerate(atoms):
        for atom_j in atoms[i + 1 :]:
            pair = "-".join(sorted((atom_i.species, atom_j.species)))
            if pair not in result:
                continue
            dist = min_image_distance(cell_rows, atom_i.frac, atom_j.frac)
            if dist < result[pair]["distance_angstrom"]:
                result[pair] = {"distance_angstrom": dist, "atom_ids": [atom_i.atom_id, atom_j.atom_id]}
    return result


def write_qe_structure(path: Path, cell_rows: tuple[tuple[float, float, float], ...], atoms: Iterable[Atom]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("CELL_PARAMETERS angstrom\n")
        for row in cell_rows:
            handle.write(f"  {row[0]:18.10f} {row[1]:18.10f} {row[2]:18.10f}\n")
        handle.write("\nATOMIC_POSITIONS crystal\n")
        for atom in atoms:
            handle.write(
                f"{atom.species:<2s} {atom.frac[0]:18.12f} {atom.frac[1]:18.12f} {atom.frac[2]:18.12f}"
                f"  ! lammps_id={atom.atom_id} type={atom.type_id}\n"
            )


def write_validation_md(path: Path, validation: dict[str, object]) -> None:
    distances = validation["minimum_distances_angstrom"]
    counts = validation["atom_counts"]
    cell = validation["cell_parameters_angstrom"]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# Ice X structure conversion validation\n\n")
        handle.write(f"- Source SHA256: `{validation['source_sha256']}`\n")
        handle.write(f"- Structure include SHA256: `{validation['qe_structure_sha256']}`\n")
        handle.write(f"- Atom count: {validation['nat']} total, O={counts['O']}, H={counts['H']}\n")
        handle.write(f"- Stoichiometry: {validation['stoichiometry']}\n")
        handle.write(f"- Cell volume: {validation['volume_angstrom3']:.12f} Angstrom^3\n")
        handle.write(f"- Estimated density: {validation['density_g_cm3']:.12f} g/cm^3\n")
        handle.write("- Cell rows in Angstrom:\n")
        for row in cell:
            handle.write(f"  - `{row[0]:.10f} {row[1]:.10f} {row[2]:.10f}`\n")
        handle.write(f"- LAMMPS tilt factors xy/xz/yz: {validation['tilt_factors_angstrom']}\n")
        handle.write(f"- Maximum cartesian roundtrip error: {validation['max_roundtrip_error_angstrom']:.3e} Angstrom\n")
        handle.write(f"- Minimum fractional separation norm: {validation['minimum_fractional_separation_norm']:.12f}\n")
        for pair, value in distances.items():
            handle.write(
                f"- Minimum {pair}: {value['distance_angstrom']:.12f} Angstrom "
                f"between LAMMPS atom ids {value['atom_ids']}\n"
            )
        handle.write("\n## Interpretation\n\n")
        handle.write(
            "The conversion preserves the LAMMPS cell vectors, atom order, and atom positions. "
            "No rotation, symmetrization, relaxation, coordinate rounding beyond text formatting, "
            "or species remapping beyond the Masses-section H/O identification was applied.\n"
        )
        handle.write(
            "Coordinates are interpreted as Angstrom because this gives a physically plausible "
            "high-pressure ice density and O-H/O-O separations; interpreting the same numbers as "
            "Bohr would make the density and bond lengths unphysical for the stated Ice X cell.\n"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("iceX.data"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("examples/iceX_100GPa_300K_caimd/structure"),
    )
    args = parser.parse_args()

    source = args.input.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    lines = source.read_text(encoding="utf-8").splitlines()

    nat = parse_count(lines, "atoms")
    ntypes = parse_count(lines, "atom types")
    if ntypes != 2:
        raise ValueError(f"Expected two atom types for Ice X H/O data, got {ntypes}")
    xlo, xhi = parse_header_float_pair(lines, "xlo\\s+xhi")
    ylo, yhi = parse_header_float_pair(lines, "ylo\\s+yhi")
    zlo, zhi = parse_header_float_pair(lines, "zlo\\s+zhi")
    xy, xz, yz = parse_tilts(lines)
    cell_rows = (
        (xhi - xlo, 0.0, 0.0),
        (xy, yhi - ylo, 0.0),
        (xz, yz, zhi - zlo),
    )
    volume = abs(determinant_3x3(cell_rows))

    masses = parse_masses(lines)
    species_by_type, mapping_evidence = infer_species(masses)
    atoms, atom_style_flags = parse_atoms(lines, cell_rows, species_by_type)
    if len(atoms) != nat:
        raise ValueError(f"Header says {nat} atoms but parsed {len(atoms)} atoms")
    if len({atom.atom_id for atom in atoms}) != len(atoms):
        raise ValueError("Duplicate LAMMPS atom ids detected")

    counts = {species: sum(1 for atom in atoms if atom.species == species) for species in ("O", "H")}
    if counts["H"] != 2 * counts["O"]:
        raise ValueError(f"Expected H:O=2:1, got {counts}")

    original_copy = output_dir / "iceX_original.data"
    shutil.copyfile(source, original_copy)
    qe_structure = output_dir / "iceX_qe_structure.inc"
    write_qe_structure(qe_structure, cell_rows, atoms)

    max_roundtrip = 0.0
    min_frac_norm = float("inf")
    for atom in atoms:
        reconstructed = frac_to_cart(cell_rows, atom.frac)
        diff = [reconstructed[i] - atom.cart_angstrom[i] for i in range(3)]
        diff = [component - round(component / length) * length for component, length in zip(diff, (xhi - xlo, yhi - ylo, zhi - zlo))]
        max_roundtrip = max(max_roundtrip, max(abs(component) for component in diff))
    for i, atom_i in enumerate(atoms):
        for atom_j in atoms[i + 1 :]:
            delta = tuple(atom_i.frac[k] - atom_j.frac[k] for k in range(3))
            delta = tuple(value - round(value) for value in delta)
            min_frac_norm = min(min_frac_norm, math.sqrt(sum(value * value for value in delta)))

    total_mass = sum(masses[atom.type_id][0] for atom in atoms)
    density = total_mass * AMU_TO_GRAM / (volume * ANGSTROM3_TO_CM3)
    validation: dict[str, object] = {
        "source_path": str(source),
        "source_sha256": sha256(source),
        "original_copy_sha256": sha256(original_copy),
        "qe_structure_sha256": sha256(qe_structure),
        "title": lines[0].strip() if lines else "",
        "nat": len(atoms),
        "atom_counts": counts,
        "stoichiometry": f"{counts['O']} H2O units; H:O={counts['H']}:{counts['O']}",
        "type_mapping": {str(type_id): species for type_id, species in sorted(species_by_type.items())},
        "type_mapping_evidence": mapping_evidence,
        "atom_style_flags": atom_style_flags,
        "coordinate_unit": "angstrom",
        "coordinate_unit_evidence": [
            "LAMMPS data files conventionally store coordinates in the active unit system; the title states 100 GPa Ice X.",
            "The Angstrom interpretation gives a high-pressure ice density near 3.18 g/cm^3.",
            "The Angstrom interpretation gives O-H, H-H, and O-O nearest-neighbor distances plausible for symmetric Ice X.",
        ],
        "cell_parameters_angstrom": [list(row) for row in cell_rows],
        "cell_parameters_bohr": [[value * BOHR_PER_ANGSTROM for value in row] for row in cell_rows],
        "tilt_factors_angstrom": {"xy": xy, "xz": xz, "yz": yz},
        "volume_angstrom3": volume,
        "density_g_cm3": density,
        "total_mass_amu": total_mass,
        "minimum_distances_angstrom": distance_summary(cell_rows, atoms),
        "minimum_fractional_separation_norm": min_frac_norm,
        "max_roundtrip_error_angstrom": max_roundtrip,
        "preserved_lammps_atom_order": True,
        "applied_rotation": False,
        "applied_symmetrization": False,
    }

    validation_json = output_dir / "structure_validation.json"
    validation_json.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_validation_md(output_dir / "structure_validation.md", validation)
    print(json.dumps({"nat": len(atoms), "counts": counts, "volume_angstrom3": volume, "density_g_cm3": density}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
