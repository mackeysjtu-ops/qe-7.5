#!/usr/bin/env python3
"""Prepare bounded QE input decks for the Ice X CAIMD example."""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path


EXAMPLE_ROOT = Path(__file__).resolve().parents[1]
STRUCTURE = EXAMPLE_ROOT / "structure" / "iceX_qe_structure.inc"
VELOCITIES = EXAMPLE_ROOT / "structure" / "iceX_velocities_300K.inc"
VALIDATION = EXAMPLE_ROOT / "structure" / "structure_validation.json"

PRESS_GPA = 100.0
PRESS_KBAR = 1000.0
NAT = 144
NTYP = 2
ECUT_RHO_FACTOR = 8.0
TENTATIVE_ECUT = 65.0
TENTATIVE_KGRID = (2, 2, 2)
TENTATIVE_DT = 5.0
TENTATIVE_STRIDE = 5
TENTATIVE_STEP_STRAIN = 5.0e-5
FINAL_NSTEP = 100
PILOT_NSTEP = 20


@dataclass(frozen=True)
class Case:
    directory: Path
    calculation: str
    prefix: str
    ecutwfc: float
    kgrid: tuple[int, int, int]
    nstep: int | None = None
    dt: float | None = None
    wmass: float | None = None
    caimd: bool = False
    caimd_stride: int = TENTATIVE_STRIDE
    caimd_step_strain: tuple[float, float, float, float, float, float] = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    caimd_mask: tuple[bool, bool, bool, bool, bool, bool] = (False, False, False, False, False, False)
    include_velocities: bool = False
    max_seconds_hint: int = 7200


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8").strip() + "\n"


def qe_bool(value: bool) -> str:
    return ".true." if value else ".false."


def d_float(value: float) -> str:
    return f"{value:.10e}".replace("e", "d")


def relative_pseudo_dir(case_dir: Path) -> str:
    return os.path.relpath(EXAMPLE_ROOT / "pseudo", case_dir)


def default_wmass_amu() -> float:
    validation = json.loads(VALIDATION.read_text(encoding="utf-8"))
    return 0.75 * float(validation["total_mass_amu"]) / math.pi**2


def electronic_block() -> str:
    return """&ELECTRONS
  conv_thr = 1.0d-8,
  electron_maxstep = 200,
  mixing_beta = 0.30,
/
"""


def species_block() -> str:
    return """ATOMIC_SPECIES
O  15.9994  O.upf
H   1.0080  H.upf
"""


def kpoints_block(kgrid: tuple[int, int, int]) -> str:
    if kgrid == (1, 1, 1):
        return """K_POINTS gamma
"""
    return f"""K_POINTS automatic
  {kgrid[0]} {kgrid[1]} {kgrid[2]} 0 0 0
"""


def cell_block(case: Case) -> str:
    if case.calculation not in {"vc-md", "vc-relax"}:
        return ""
    lines = [
        "&CELL",
        "  cell_dynamics = 'pr',",
        "  cell_dofree = 'all',",
        f"  press = {PRESS_KBAR:.1f},",
    ]
    if case.wmass is not None:
        lines.append(f"  wmass = {case.wmass:.10f},")
    if case.calculation == "vc-md":
        lines.append("  cell_temperature = 'not_controlled',")
    if case.caimd:
        lines.append("  caimd = .true.,")
        lines.append(f"  caimd_stride = {case.caimd_stride},")
        for i, value in enumerate(case.caimd_step_strain, start=1):
            if value != 0.0:
                lines.append(f"  caimd_step_strain({i}) = {d_float(value)},")
        for i, value in enumerate(case.caimd_mask, start=1):
            if value:
                lines.append(f"  caimd_mask({i}) = {qe_bool(value)},")
    else:
        lines.append("  caimd = .false.,")
    lines.append("/")
    return "\n".join(lines) + "\n"


def ions_block(case: Case) -> str:
    if case.calculation == "vc-md":
        return """&IONS
  ion_dynamics = 'beeman',
  ion_temperature = 'not_controlled',
/
"""
    if case.calculation == "vc-relax":
        return """&IONS
  ion_dynamics = 'bfgs',
/
"""
    return ""


def control_block(case: Case) -> str:
    case_dir = case.directory
    pseudo_dir = relative_pseudo_dir(case_dir)
    lines = [
        "&CONTROL",
        f"  calculation = '{case.calculation}',",
        "  restart_mode = 'from_scratch',",
        f"  prefix = '{case.prefix}',",
        f"  pseudo_dir = '{pseudo_dir}',",
        "  outdir = './tmp',",
        "  tstress = .true.,",
        "  tprnfor = .true.,",
        "  verbosity = 'high',",
        "  disk_io = 'low',",
    ]
    if case.nstep is not None:
        lines.append(f"  nstep = {case.nstep},")
    if case.dt is not None:
        lines.append(f"  dt = {case.dt:.4f},")
    lines.append("/")
    return "\n".join(lines) + "\n"


def system_block(case: Case) -> str:
    ecutrho = case.ecutwfc * ECUT_RHO_FACTOR
    return f"""&SYSTEM
  ibrav = 0,
  nat = {NAT},
  ntyp = {NTYP},
  ecutwfc = {case.ecutwfc:.1f},
  ecutrho = {ecutrho:.1f},
  occupations = 'fixed',
  nosym = .true.,
  noinv = .true.,
/
"""


def write_case(case: Case) -> None:
    case.directory.mkdir(parents=True, exist_ok=True)
    parts = [
        control_block(case),
        system_block(case),
        electronic_block(),
        ions_block(case),
        cell_block(case),
        species_block(),
        kpoints_block(case.kgrid),
        read(STRUCTURE),
    ]
    if case.include_velocities:
        parts.append(read(VELOCITIES))
    (case.directory / "qe.in").write_text("\n".join(part for part in parts if part), encoding="utf-8")
    (case.directory / "CASE_METADATA.json").write_text(
        json.dumps(
            {
                "prefix": case.prefix,
                "calculation": case.calculation,
                "ecutwfc_Ry": case.ecutwfc,
                "ecutrho_Ry": case.ecutwfc * ECUT_RHO_FACTOR,
                "kgrid": case.kgrid,
                "press_input_kbar": PRESS_KBAR,
                "press_target_GPa": PRESS_GPA,
                "dt_rydberg_au": case.dt,
                "wmass_amu": case.wmass,
                "nstep": case.nstep,
                "caimd": case.caimd,
                "caimd_stride": case.caimd_stride,
                "caimd_step_strain_order_xx_yy_zz_yz_xz_xy": case.caimd_step_strain,
                "caimd_mask_order_xx_yy_zz_yz_xz_xy": case.caimd_mask,
                "include_300K_velocities": case.include_velocities,
                "max_seconds_hint": case.max_seconds_hint,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def md_case(
    directory: Path,
    prefix: str,
    *,
    dt: float = TENTATIVE_DT,
    wmass: float | None = None,
    nstep: int = PILOT_NSTEP,
    caimd: bool = False,
    step: tuple[float, float, float, float, float, float] = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    mask: tuple[bool, bool, bool, bool, bool, bool] = (False, False, False, False, False, False),
    ecutwfc: float = TENTATIVE_ECUT,
    kgrid: tuple[int, int, int] = TENTATIVE_KGRID,
    max_seconds_hint: int = 7200,
) -> Case:
    return Case(
        directory=directory,
        calculation="vc-md",
        prefix=prefix,
        ecutwfc=ecutwfc,
        kgrid=kgrid,
        nstep=nstep,
        dt=dt,
        wmass=wmass,
        caimd=caimd,
        caimd_step_strain=step,
        caimd_mask=mask,
        include_velocities=True,
        max_seconds_hint=max_seconds_hint,
    )


def all_cases() -> list[Case]:
    wmass = default_wmass_amu()
    half = 0.5 * wmass
    double = 2.0 * wmass
    root = EXAMPLE_ROOT
    cases: list[Case] = []

    for ecut in (35.0, 45.0, 55.0, 65.0):
        cases.append(
            Case(
                directory=root / "convergence" / "ecut" / f"ecut_{int(ecut)}Ry",
                calculation="scf",
                prefix=f"icex_ecut_{int(ecut)}Ry",
                ecutwfc=ecut,
                kgrid=(1, 1, 1),
                max_seconds_hint=5400,
            )
        )
    for kgrid in ((1, 1, 1), (2, 2, 2), (3, 3, 3)):
        label = "gamma" if kgrid == (1, 1, 1) else f"{kgrid[0]}x{kgrid[1]}x{kgrid[2]}"
        cases.append(
            Case(
                directory=root / "convergence" / "kpoints" / f"k_{label}",
                calculation="scf",
                prefix=f"icex_k_{label}",
                ecutwfc=TENTATIVE_ECUT,
                kgrid=kgrid,
                max_seconds_hint=7200,
            )
        )

    cases.extend(
        [
            Case(
                directory=root / "reference" / "scf",
                calculation="scf",
                prefix="icex_reference_scf",
                ecutwfc=TENTATIVE_ECUT,
                kgrid=TENTATIVE_KGRID,
                max_seconds_hint=7200,
            ),
            Case(
                directory=root / "reference" / "optional_vc_relax",
                calculation="vc-relax",
                prefix="icex_optional_vc_relax",
                ecutwfc=TENTATIVE_ECUT,
                kgrid=TENTATIVE_KGRID,
                nstep=20,
                wmass=wmass,
                max_seconds_hint=7200,
            ),
            md_case(root / "reference" / "zero_strain", "icex_zero_strain", wmass=wmass, nstep=20),
            md_case(root / "pilot" / "parameter_screen" / "dt5_wmass_default", "icex_dt5_wmass_default", dt=5.0, wmass=wmass),
            md_case(root / "pilot" / "parameter_screen" / "dt7p5_wmass_default", "icex_dt7p5_wmass_default", dt=7.5, wmass=wmass),
            md_case(root / "pilot" / "parameter_screen" / "dt10_wmass_default", "icex_dt10_wmass_default", dt=10.0, wmass=wmass),
            md_case(root / "pilot" / "parameter_screen" / "dt5_wmass_half", "icex_dt5_wmass_half", dt=5.0, wmass=half),
            md_case(root / "pilot" / "parameter_screen" / "dt5_wmass_double", "icex_dt5_wmass_double", dt=5.0, wmass=double),
            md_case(
                root / "pilot" / "tension" / "xx_step1e-5_stride5",
                "icex_pilot_tension_xx_1e5",
                wmass=wmass,
                caimd=True,
                step=(1.0e-5, 0.0, 0.0, 0.0, 0.0, 0.0),
                mask=(True, False, False, False, False, False),
            ),
            md_case(
                root / "pilot" / "tension" / "xx_step5e-5_stride5",
                "icex_pilot_tension_xx_5e5",
                wmass=wmass,
                caimd=True,
                step=(TENTATIVE_STEP_STRAIN, 0.0, 0.0, 0.0, 0.0, 0.0),
                mask=(True, False, False, False, False, False),
            ),
            md_case(
                root / "pilot" / "compression" / "xx_stepm5e-5_stride5",
                "icex_pilot_compression_xx_m5e5",
                wmass=wmass,
                caimd=True,
                step=(-TENTATIVE_STEP_STRAIN, 0.0, 0.0, 0.0, 0.0, 0.0),
                mask=(True, False, False, False, False, False),
            ),
            md_case(
                root / "pilot" / "shear" / "xy_step5e-5_stride5",
                "icex_pilot_shear_xy_5e5",
                wmass=wmass,
                caimd=True,
                step=(0.0, 0.0, 0.0, 0.0, 0.0, TENTATIVE_STEP_STRAIN),
                mask=(False, False, False, False, False, True),
            ),
            md_case(
                root / "deformation" / "tension",
                "icex_caimd_tension_xx",
                wmass=wmass,
                nstep=FINAL_NSTEP,
                caimd=True,
                step=(TENTATIVE_STEP_STRAIN, 0.0, 0.0, 0.0, 0.0, 0.0),
                mask=(True, False, False, False, False, False),
                max_seconds_hint=86400,
            ),
            md_case(
                root / "deformation" / "compression",
                "icex_caimd_compression_xx",
                wmass=wmass,
                nstep=FINAL_NSTEP,
                caimd=True,
                step=(-TENTATIVE_STEP_STRAIN, 0.0, 0.0, 0.0, 0.0, 0.0),
                mask=(True, False, False, False, False, False),
                max_seconds_hint=86400,
            ),
            md_case(
                root / "deformation" / "shear",
                "icex_caimd_shear_xy",
                wmass=wmass,
                nstep=FINAL_NSTEP,
                caimd=True,
                step=(0.0, 0.0, 0.0, 0.0, 0.0, TENTATIVE_STEP_STRAIN),
                mask=(False, False, False, False, False, True),
                max_seconds_hint=86400,
            ),
        ]
    )
    return cases


def write_case_lists(cases: list[Case]) -> None:
    convergence = [case for case in cases if "convergence" in case.directory.parts]
    pilots = [case for case in cases if "pilot" in case.directory.parts or case.directory.parts[-2:] == ("reference", "zero_strain")]
    (EXAMPLE_ROOT / "convergence" / "cases.txt").write_text(
        "\n".join(os.path.relpath(case.directory, EXAMPLE_ROOT) for case in convergence) + "\n",
        encoding="utf-8",
    )
    (EXAMPLE_ROOT / "pilot" / "pilot_cases.txt").write_text(
        "\n".join(os.path.relpath(case.directory, EXAMPLE_ROOT) for case in pilots) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    cases = all_cases()
    for case in cases:
        write_case(case)
    write_case_lists(cases)
    print(json.dumps({"cases_written": len(cases), "default_wmass_amu": default_wmass_amu()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
