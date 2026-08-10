#!/usr/bin/env python3
"""Generate the two QE inputs for the no-explicit-velocity 300 K comparison.

Both inputs describe exactly the same physical system and dynamics.
The only permitted differences are:
  - prefix
  - outdir
  - the QE 7.5-only line `caimd = .false.`
  - the leading comment identifying the target executable

The structure is the IDEAL, UNSTRAINED 8-atom diamond conventional cell
recovered from the validated Diamond8 setup metadata (a = 3.567 Angstrom
= 6.740653089508073 bohr), NOT a terminal state from a failed trajectory.

Neither input contains ATOMIC_VELOCITIES or ion_velocities: both rely on
QE's native from-scratch Maxwell initialization via vcinit -> ranv(temp),
driven by ion_temperature='initial' with tempw=300.0.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# ---- ideal unstrained Diamond8 structure (validated setup metadata) --------
# Source: examples/caimd_diamond8_local_nose_four_dof_tests_0p03_.../
#         00_metadata/DIAMOND8_SETUP.json  (lattice_constant_angstrom = 3.567)
CELL_ROWS_BOHR = [
    [3.891717875741369,  4.7663615092172105,  2.7518601004016277],
    [3.891717875741369, -4.7663615092172105,  2.7518601004016277],
    [3.891717875741369,  0.0,                -5.503720200803255],
]
FRAC_POSITIONS = [
    [0.00, 0.00, 0.00],
    [0.00, 0.50, 0.50],
    [0.50, 0.00, 0.50],
    [0.50, 0.50, 0.00],
    [0.25, 0.25, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
]

# ---- common physical / numerical parameters -------------------------------
NSTEP = 1000
DT_QE = 41.34137333665          # 2 fs in QE Rydberg atomic time units
ECUTWFC = 25.0
ECUTRHO = 120.0
DEGAUSS = 0.01
CONV_THR = 1.0e-4
ELECTRON_MAXSTEP = 40
MIXING_BETA = 0.300
TEMPW = 300.0
PRESS = 0.0
WMASS = 2000000.0
MASS_C = 12.01070000


def build(prefix, outdir, caimd_line, header):
    L = []
    L.append("! %s" % header)
    L.append("! Native from-scratch 300 K initialization: ion_temperature='initial'.")
    L.append("! Explicit velocity input is deliberately omitted so that QE generates")
    L.append("! its own Maxwell distribution at tempw via the native vcinit path.")
    L.append("&CONTROL")
    L.append("  calculation = 'vc-md',")
    L.append("  restart_mode = 'from_scratch',")
    L.append("  prefix = '%s'," % prefix)
    L.append("  pseudo_dir = '../pseudo',")
    L.append("  outdir = '%s'," % outdir)
    L.append("  nstep = %d," % NSTEP)
    L.append("  dt = %.11f," % DT_QE)
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
    L.append("  ecutwfc = %.1f," % ECUTWFC)
    L.append("  ecutrho = %.1f," % ECUTRHO)
    L.append("  occupations = 'smearing',")
    L.append("  smearing = 'gaussian',")
    L.append("  degauss = %.2f," % DEGAUSS)
    L.append("  nosym = .true.,")
    L.append("  noinv = .true.,")
    L.append("/")
    L.append("&ELECTRONS")
    L.append("  conv_thr = %.16e," % CONV_THR)
    L.append("  electron_maxstep = %d," % ELECTRON_MAXSTEP)
    L.append("  mixing_beta = %.3f," % MIXING_BETA)
    L.append("/")
    L.append("&IONS")
    L.append("  ion_dynamics = 'beeman',")
    L.append("  ion_temperature = 'initial',")
    L.append("  tempw = %.1f," % TEMPW)
    L.append("/")
    L.append("&CELL")
    L.append("  cell_dynamics = 'pr',")
    L.append("  cell_dofree = 'all',")
    L.append("  press = %.1f," % PRESS)
    L.append("  wmass = %.8f," % WMASS)
    L.append("  cell_temperature = 'not_controlled',")
    if caimd_line:
        L.append("  %s" % caimd_line)
    L.append("/")
    L.append("ATOMIC_SPECIES")
    L.append("C %.8f C.UPF" % MASS_C)
    L.append("CELL_PARAMETERS bohr")
    for r in CELL_ROWS_BOHR:
        L.append("  %.14f  %.14f  %.14f" % tuple(r))
    L.append("ATOMIC_POSITIONS crystal")
    for p in FRAC_POSITIONS:
        L.append("C  %.14f  %.14f  %.14f" % tuple(p))
    L.append("K_POINTS gamma")
    L.append("")
    return "\n".join(L)


def main():
    idir = os.path.join(ROOT, "inputs")
    os.makedirs(idir, exist_ok=True)

    qe74 = build(
        prefix="diamond8_qe74_novel",
        outdir="./tmp",
        caimd_line=None,
        header="Official unmodified Quantum ESPRESSO 7.4.1 from the Leonardo module stack.",
    )
    qe75 = build(
        prefix="diamond8_qe75_novel",
        outdir="./tmp",
        caimd_line="caimd = .false.",
        header="QE 7.5 CAIMD-development build at commit af7ebfc, CAIMD explicitly disabled.",
    )

    with open(os.path.join(idir, "qe74_no_velocity.in"), "w") as fh:
        fh.write(qe74)
    with open(os.path.join(idir, "qe75_caimd_false_no_velocity.in"), "w") as fh:
        fh.write(qe75)
    print("wrote inputs/qe74_no_velocity.in")
    print("wrote inputs/qe75_caimd_false_no_velocity.in")


if __name__ == "__main__":
    main()
