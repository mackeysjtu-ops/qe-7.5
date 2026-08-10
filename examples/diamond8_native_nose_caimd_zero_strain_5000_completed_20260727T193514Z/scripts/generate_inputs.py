#!/usr/bin/env python3
"""Generate the three 5000-step (10 ps) QE 7.5 inputs for the CAIMD zero-strain long test.

All three share: same structure/electronics/PR, ion_temperature='nose' (tempw=300,
fnosep=50, nhpcl=4), 2 fs timestep, nstep=5000, no explicit velocities, from_scratch.

The ONLY intended differences are the CAIMD lines in &CELL:
  CASE A: caimd = .false.                                    (ordinary control)
  CASE B: caimd = .true. + stride=1 + all six masks .false. + all six step_strain 0.0  (no-op)
  CASE C: caimd = .true. + stride=1 + all six masks .true.  + all six step_strain 0.0  (zero-strain hold)

Structure/electronics are taken verbatim from the validated reference input
examples/diamond8_af7ebfc_no_velocity_initial_vs_nose_completed_20260727T081614Z/
inputs/test_nose_no_velocity.in (sha256 95e4ac6a...), which produced the bounded
1000-step Nose trajectory (mean ~170 K, no step above 600 K).

ecutwfc=25 Ry is retained deliberately for strict comparison with prior results;
this is a stability/temperature-control diagnostic, NOT a converged pressure run.

CAIMD namelist syntax (per-element indexed) matches the proven form used by
test-suite/caimd/caimd_validation.py (e.g. caimd_mask(1)=.true., caimd_step_strain(1)=0.0).
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
NSTEP = 5000
DT_QE = 41.34137333665
TEMPW = 300.0
FNOSEP = 50.0
NHPCL = 4


def build(prefix, caimd_lines, header):
    L = []
    L.append(f"! {header}")
    L.append("! 5000-step (10 ps) QE 7.5 af7ebfc; ion_temperature='nose', tempw=300, fnosep=50, nhpcl=4.")
    L.append("! No explicit velocities: QE vcinit generates the native 300 K Maxwell distribution.")
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
    L.append("  ion_temperature = 'nose',")
    L.append(f"  tempw = {TEMPW:.1f},")
    L.append(f"  fnosep(1) = {FNOSEP:.1f},")
    L.append(f"  nhpcl = {NHPCL},")
    L.append("/")
    L.append("&CELL")
    L.append("  cell_dynamics = 'pr',")
    L.append("  cell_dofree = 'all',")
    L.append("  press = 0.0,")
    L.append("  wmass = 2000000.00000000,")
    L.append("  cell_temperature = 'not_controlled',")
    for cl in caimd_lines:
        L.append(f"  {cl}")
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


def mask_false_lines():
    lines = ["caimd = .true.", "caimd_stride = 1"]
    lines += [f"caimd_step_strain({i}) = 0.0," for i in range(1, 7)]
    lines += [f"caimd_mask({i}) = .false.," for i in range(1, 7)]
    # strip trailing comma on the last line for cleanliness (namelist tolerates it anyway)
    lines[-1] = lines[-1].rstrip(",")
    return lines


def mask_true_lines():
    lines = ["caimd = .true.", "caimd_stride = 1"]
    lines += [f"caimd_step_strain({i}) = 0.0," for i in range(1, 7)]
    lines += [f"caimd_mask({i}) = .true.," for i in range(1, 7)]
    lines[-1] = lines[-1].rstrip(",")
    return lines


def main():
    idir = os.path.join(ROOT, "inputs")
    os.makedirs(idir, exist_ok=True)
    cases = {
        "case_A_caimd_false_nose_5000.in":
            (build("diamond8_case_A_caimd_false", ["caimd = .false."],
                   "CASE A: ordinary PR vc-md + Nose, CAIMD disabled (control).")),
        "case_B_caimd_true_all_false_nose_5000.in":
            (build("diamond8_case_B_caimd_noop", mask_false_lines(),
                   "CASE B: CAIMD enabled, all six masks .false., zero step_strain (no-op boundary).")),
        "case_C_caimd_true_all_fixed_zero_strain_nose_5000.in":
            (build("diamond8_case_C_caimd_zero", mask_true_lines(),
                   "CASE C: active CAIMD, all six masks .true., zero step_strain (zero-strain hold).")),
    }
    for fname, txt in cases.items():
        with open(os.path.join(idir, fname), "w") as fh:
            fh.write(txt)
        print("wrote inputs/" + fname)


if __name__ == "__main__":
    main()
