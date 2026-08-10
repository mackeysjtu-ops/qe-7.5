#!/usr/bin/env python3
"""Generate the CONTROL and TEST inputs for the single-factor thermostat isolation.

The ONLY intended difference is the ionic thermostat block:

  CONTROL : ion_temperature='initial'                       (no continuing thermostat)
  TEST    : ion_temperature='nose', fnosep(1)=50.0, nhpcl=4 (continuing Nose-Hoover)

Everything else - executable, structure, cell, positions, ordering, masses,
pseudopotential, cutoffs, SCF settings, timestep, nstep, PR settings, caimd=.false.
- is byte-identical apart from prefix/outdir.

Structure and electronic settings are taken verbatim from the completed
reference workflow input
  examples/diamond8_no_velocity_qe74_vs_qe75_caimd_false_completed_20260724T113644Z/
  inputs/qe75_caimd_false_no_velocity.in
(sha256 a7ad97b873f1fda24892695802164e7b6e86e2ceb06e6f5193d2238980f4fd39),
which produced the bounded control trajectory (mean 170.014 K, no step above 600 K).

Note: ecutwfc=25 Ry is deliberately retained from that reference deck. This is a
controlled thermostat diagnosis, NOT a converged pressure calculation; absolute
pressures from this input are dominated by Pulay stress and must not be quoted
as physical.

`fnosep(1) = 50.0` is the QE-validated indexed form of `fnosep = 50.0`
(a Fortran namelist assignment `fnosep = 50.0` sets exactly the same element).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

CELL_ROWS_BOHR = [
    [3.89171787574137,  4.76636150921721,  2.75186010040163],
    [3.89171787574137, -4.76636150921721,  2.75186010040163],
    [3.89171787574137,  0.00000000000000, -5.50372020080326],
]
FRAC = [
    [0.00, 0.00, 0.00], [0.00, 0.50, 0.50], [0.50, 0.00, 0.50], [0.50, 0.50, 0.00],
    [0.25, 0.25, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75], [0.75, 0.75, 0.25],
]

NSTEP = 1000
DT_QE = 41.34137333665
TEMPW = 300.0
FNOSEP = 50.0
NHPCL = 4


def build(prefix, ions_lines, header):
    L = []
    L.append(f"! {header}")
    L.append("! Single-factor thermostat isolation, QE 7.5 af7ebfc, CAIMD disabled.")
    L.append("! Explicit velocity input is deliberately omitted so that QE generates")
    L.append("! its own Maxwell distribution at tempw via the native vcinit path.")
    L.append("&CONTROL")
    L.append("  calculation = 'vc-md',")
    L.append("  restart_mode = 'from_scratch',")
    L.append(f"  prefix = '{prefix}',")
    L.append("  pseudo_dir = '../pseudo',")
    L.append("  outdir = './tmp',")
    L.append(f"  nstep = {NSTEP},")
    L.append(f"  dt = {DT_QE:.11f},")
    L.append("  iprint = 1,")
    L.append("  tstress = .true.,")
    L.append("  tprnfor = .true.,")
    L.append("  disk_io = 'low',")
    L.append("  verbosity = 'high',")
    L.append("/")
    L.append("&SYSTEM")
    L.append("  ibrav = 0,")
    L.append("  nat = 8,")
    L.append("  ntyp = 1,")
    L.append("  ecutwfc = 25.0,")
    L.append("  ecutrho = 120.0,")
    L.append("  occupations = 'smearing',")
    L.append("  smearing = 'gaussian',")
    L.append("  degauss = 0.01,")
    L.append("  nosym = .true.,")
    L.append("  noinv = .true.,")
    L.append("/")
    L.append("&ELECTRONS")
    L.append("  conv_thr = 1.0000000000000000e-04,")
    L.append("  electron_maxstep = 40,")
    L.append("  mixing_beta = 0.300,")
    L.append("/")
    L.append("&IONS")
    L.append("  ion_dynamics = 'beeman',")
    L.extend(ions_lines)
    L.append("/")
    L.append("&CELL")
    L.append("  cell_dynamics = 'pr',")
    L.append("  cell_dofree = 'all',")
    L.append("  press = 0.0,")
    L.append("  wmass = 2000000.00000000,")
    L.append("  cell_temperature = 'not_controlled',")
    L.append("  caimd = .false.")
    L.append("/")
    L.append("ATOMIC_SPECIES")
    L.append("C 12.01070000 C.UPF")
    L.append("CELL_PARAMETERS bohr")
    for r in CELL_ROWS_BOHR:
        L.append("  %.14f  %.14f  %.14f" % tuple(r))
    L.append("ATOMIC_POSITIONS crystal")
    for p in FRAC:
        L.append("C  %.14f  %.14f  %.14f" % tuple(p))
    L.append("K_POINTS gamma")
    L.append("")
    return "\n".join(L)


def main():
    idir = os.path.join(ROOT, "inputs")
    os.makedirs(idir, exist_ok=True)

    control = build(
        prefix="diamond8_control_initial",
        ions_lines=[f"  ion_temperature = 'initial',", f"  tempw = {TEMPW:.1f},"],
        header="CONTROL: native 300 K initialization only, no continuing ionic thermostat.",
    )
    nose = build(
        prefix="diamond8_test_nose",
        ions_lines=[f"  ion_temperature = 'nose',", f"  tempw = {TEMPW:.1f},",
                    f"  fnosep(1) = {FNOSEP:.1f},", f"  nhpcl = {NHPCL},"],
        header="TEST: continuing Nose-Hoover ionic thermostat (fnosep=50 THz, nhpcl=4).",
    )

    with open(os.path.join(idir, "control_initial_no_velocity.in"), "w") as fh:
        fh.write(control)
    with open(os.path.join(idir, "test_nose_no_velocity.in"), "w") as fh:
        fh.write(nose)
    print("wrote inputs/control_initial_no_velocity.in")
    print("wrote inputs/test_nose_no_velocity.in")


if __name__ == "__main__":
    main()
