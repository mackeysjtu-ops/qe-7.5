#!/usr/bin/env python3
"""Robust final-structure extraction + downstream input generation for the chain.

Modes:
  --mode vcrelax_gate     : parse JOB1 vc-relax qe.out; write final_structure.json / fragment /
                            pressure_history.csv / volume_history.csv; apply the vc-relax pass gate.
  --mode gen_nvt          : read JOB1 final_structure.json -> write JOB2 fixed-cell md input (ibrav=1).
  --mode gen_npt          : read JOB1 final_structure.json + JOB2 final .cel/.pos -> JOB3 vc-md input (ibrav=1).
  --mode gen_caimd_zero   : read JOB3 final cell/pos -> JOB4 stage-4A input (ibrav=0, cell_dofree=all).
  --mode gen_caimd_nonzero: read stage-4A final cell/pos -> JOB4 stage-4B input.

Structure handoff carries cell + positions + ordering only. NO velocities are ever written.
For diamond volume-only relaxation the fractional positions are symmetry-fixed, but this script
reads the ACTUAL relaxed positions from output and reuses them, validating 8 atoms / positive
determinant / finite values / no duplicates.
"""
import argparse
import json
import math
import os
import re
import sys

NUM = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")
BOHR_TO_ANG = 0.529177210903


def text(p):
    return open(p, errors="replace").read() if os.path.exists(p) else ""


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def parse_last_cell_positions_from_out(out_txt):
    """Extract the LAST 'CELL_PARAMETERS' + 'ATOMIC_POSITIONS' blocks from a QE vc-relax stdout.
    Returns (cell_bohr_rows, alat_bohr, species_positions[(sym,[x,y,z],units)]) or raises."""
    lines = out_txt.splitlines()
    # find all CELL_PARAMETERS blocks
    cell_idx = [i for i, l in enumerate(lines) if l.strip().startswith("CELL_PARAMETERS")]
    pos_idx = [i for i, l in enumerate(lines) if l.strip().startswith("ATOMIC_POSITIONS")]
    if not cell_idx or not pos_idx:
        raise ValueError("no CELL_PARAMETERS/ATOMIC_POSITIONS blocks in output")
    ci = cell_idx[-1]
    header = lines[ci].strip()
    m = re.search(r"alat=\s*([0-9.EeDd+-]+)", header)
    unit = "alat" if "alat" in header else ("bohr" if "bohr" in header else "angstrom")
    cell = []
    for r in range(1, 4):
        v = [float(x.replace("D", "E")) for x in NUM.findall(lines[ci + r])]
        cell.append(v[:3])
    alat = float(m.group(1).replace("D", "E")) if m else None
    if unit == "alat":
        if alat is None:
            raise ValueError("alat cell without alat value")
        cell = [[c * alat for c in row] for row in cell]     # -> bohr
    elif unit == "angstrom":
        cell = [[c / BOHR_TO_ANG for c in row] for row in cell]
    # positions: last block
    pi = pos_idx[-1]
    punit = "alat"
    hp = lines[pi].strip()
    if "crystal" in hp:
        punit = "crystal"
    elif "bohr" in hp:
        punit = "bohr"
    elif "angstrom" in hp:
        punit = "angstrom"
    sp = []
    j = pi + 1
    while j < len(lines):
        f = lines[j].split()
        if len(f) >= 4 and re.match(r"^[A-Za-z]{1,2}$", f[0]):
            sp.append((f[0], [float(x.replace("D", "E")) for x in f[1:4]], punit))
            j += 1
        else:
            break
    return cell, alat, sp, punit


def to_crystal(sp, punit, cell_bohr, alat):
    """Return fractional (crystal) positions from whatever unit QE printed."""
    import numpy as np
    H = np.array(cell_bohr)          # rows = lattice vectors (bohr)
    Hinv = np.linalg.inv(H.T)        # cartesian(bohr) = frac @ H  -> frac = cart @ Hinv? careful
    out = []
    for sym, xyz, _ in sp:
        if punit == "crystal":
            out.append((sym, list(xyz)))
        else:
            xyz = list(xyz)
            if punit == "alat":
                cart = [c * alat for c in xyz]           # bohr
            elif punit == "angstrom":
                cart = [c / BOHR_TO_ANG for c in xyz]
            else:
                cart = xyz                                # bohr
            frac = (np.array(cart) @ np.linalg.inv(H)).tolist()  # cart = frac@H (rows) -> frac = cart @ H^-1
            out.append((sym, frac))
    return out


def validate_struct(cell_bohr, frac):
    errs = []
    if len(frac) != 8:
        errs.append(f"atom count {len(frac)} != 8")
    d = det3(cell_bohr)
    if not (d > 0):
        errs.append("cell determinant not positive")
    for row in cell_bohr:
        if any((not math.isfinite(x)) for x in row):
            errs.append("nonfinite cell value")
    for sym, xyz in frac:
        if any((not math.isfinite(x)) for x in xyz):
            errs.append("nonfinite position")
    # duplicates in fractional (mod 1)
    import itertools
    for a, b in itertools.combinations(range(len(frac)), 2):
        if all(abs(((frac[a][1][k] - frac[b][1][k]) + 0.5) % 1.0 - 0.5) < 1e-6 for k in range(3)):
            errs.append(f"duplicate atoms {a},{b}")
    return errs, d


DFT_HEAD = """&CONTROL
  calculation = '{calc}',
  restart_mode = 'from_scratch',
  prefix = '{prefix}',
  pseudo_dir = '{pseudo_dir}',
  outdir = './tmp',
{control_extra}/
&SYSTEM
{system}/
&ELECTRONS
  conv_thr = 1.0d-8,
  electron_maxstep = 120,
  mixing_beta = 0.30,
  diagonalization = 'david',
/
"""


def system_ibrav1(celldm1, extra=""):
    return (f"  ibrav = 1,\n  celldm(1) = {celldm1:.10f},\n  nat = 8,\n  ntyp = 1,\n"
            f"  ecutwfc = 80.0,\n  ecutrho = 320.0,\n  occupations = 'fixed',\n"
            f"  nosym = .true.,\n  noinv = .true.,\n{extra}")


def system_ibrav0(extra=""):
    return ("  ibrav = 0,\n  nat = 8,\n  ntyp = 1,\n  ecutwfc = 80.0,\n  ecutrho = 320.0,\n"
            "  occupations = 'fixed',\n  nosym = .true.,\n  noinv = .true.,\n" + extra)


def species_cards(frac, cell_bohr=None, ibrav=1):
    L = ["ATOMIC_SPECIES", "C 12.01070000 C.UPF"]
    if ibrav == 0:
        L.append("CELL_PARAMETERS bohr")
        for row in cell_bohr:
            L.append("  %.12f  %.12f  %.12f" % tuple(row))
    L.append("ATOMIC_POSITIONS crystal")
    for sym, xyz in frac:
        L.append("C  %.12f  %.12f  %.12f" % tuple(xyz))
    L.append("K_POINTS automatic")
    L.append("2 2 2 0 0 0")
    L.append("")
    return "\n".join(L)


def write_json_struct(path, cell_bohr, frac, celldm1, source):
    json.dump({"cell_bohr_rows": cell_bohr, "fractional_positions": [[s, list(x)] for s, x in frac],
               "celldm1_bohr": celldm1, "determinant_bohr3": det3(cell_bohr), "n_atoms": len(frac),
               "source": source}, open(path, "w"), indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True)
    ap.add_argument("--job1-out")
    ap.add_argument("--prev-cel")
    ap.add_argument("--prev-pos")
    ap.add_argument("--struct-json")
    ap.add_argument("--out-input")
    ap.add_argument("--out-struct")
    ap.add_argument("--out-report")
    ap.add_argument("--out-fragment")
    ap.add_argument("--out-pressure")
    ap.add_argument("--out-volume")
    ap.add_argument("--out-handoff")
    ap.add_argument("--prefix", default="d8")
    ap.add_argument("--pseudo-dir", default="../pseudo")
    a = ap.parse_args()
    import numpy as np

    if a.mode == "vcrelax_gate":
        out = text(a.job1_out)
        fails = []
        job_done = "JOB DONE" in out
        converged = "End of BFGS Geometry Optimization" in out or "bfgs converged" in out.lower()
        err = re.search(r"Error in routine\s+(\S+)", out)
        # pressure & volume history
        pmatches = re.findall(r"P=\s*([-+0-9.EeDd]+)", out)
        pressures = [float(x.replace("D", "E")) for x in pmatches]
        # total energy / volume history (unit-cell volume)
        vmatches = re.findall(r"unit-cell volume\s*=\s*([0-9.EeDd+-]+)", out)
        volumes = [float(x.replace("D", "E")) for x in vmatches]
        # final force
        fmax = None
        for m in re.finditer(r"Total force =\s*([0-9.EeDd+-]+)", out):
            fmax = float(m.group(1).replace("D", "E"))
        try:
            cell, alat, sp, punit = parse_last_cell_positions_from_out(out)
            frac = to_crystal(sp, punit, cell, alat)
            errs, d = validate_struct(cell, frac)
        except Exception as e:
            errs, d, cell, frac, alat = [f"parse error: {e}"], None, None, None, None
        # celldm1 from relaxed cubic cell: a = det^(1/3) for cubic
        celldm1 = (abs(d)) ** (1.0 / 3.0) if d else None
        final_P = pressures[-1] if pressures else None
        if not job_done:
            fails.append("no JOB DONE")
        if err:
            fails.append(f"QE error {err.group(1)}")
        if not converged:
            fails.append("BFGS did not converge (hit step limit or incomplete)")
        fails += errs
        if final_P is not None and abs(final_P) > 5.0:
            fails.append(f"final |P| {abs(final_P):.2f} kbar > 5")
        if final_P is None:
            fails.append("no pressure parsed")
        summary = {"stage": "vcrelax", "job_done": job_done, "converged": converged,
                   "final_pressure_kbar": final_P, "final_max_force": fmax,
                   "n_atoms": len(frac) if frac else None, "determinant_bohr3": d,
                   "celldm1_bohr": celldm1, "n_pressure_samples": len(pressures), "n_volume_samples": len(volumes),
                   "gate_failures": fails, "gate_pass": not fails}
        os.makedirs(os.path.dirname(a.out_struct), exist_ok=True)
        if frac and d:
            write_json_struct(a.out_struct, cell, frac, celldm1, "JOB1 vc-relax final")
            with open(a.out_fragment, "w") as fh:
                fh.write("# relaxed structure fragment (ibrav=1)\n")
                fh.write(f"celldm(1) = {celldm1:.10f}\n")
                fh.write("ATOMIC_POSITIONS crystal\n")
                for sym, xyz in frac:
                    fh.write("C  %.12f  %.12f  %.12f\n" % tuple(xyz))
        with open(a.out_pressure, "w") as fh:
            fh.write("index,P_kbar\n")
            for i, p in enumerate(pressures):
                fh.write(f"{i},{p}\n")
        with open(a.out_volume, "w") as fh:
            fh.write("index,volume_bohr3\n")
            for i, v in enumerate(volumes):
                fh.write(f"{i},{v}\n")
        json.dump(summary, open(os.path.join(os.path.dirname(a.out_struct), "analysis_summary.json"), "w"), indent=2)
        with open(a.out_report, "w") as fh:
            fh.write("# Final vc-relax report\n\n")
            fh.write(f"Status: {'PASS' if not fails else 'FAIL'}\n\n")
            fh.write(f"- JOB DONE: {job_done}; BFGS converged: {converged}\n")
            fh.write(f"- final pressure: {final_P} kbar; final max force: {fmax}\n")
            fh.write(f"- relaxed celldm(1): {celldm1} bohr; determinant: {d} bohr^3; atoms: {len(frac) if frac else None}\n")
            fh.write(f"- gate failures: {fails or 'none'}\n")
        print(json.dumps({"gate_pass": not fails, "failures": fails, "final_P": final_P, "celldm1": celldm1}, indent=2))
        sys.exit(0 if not fails else 1)

    # ---- input generation modes ----
    def load_struct(pth):
        j = json.load(open(pth))
        cell = j["cell_bohr_rows"]; frac = [(s, x) for s, x in j["fractional_positions"]]
        return cell, frac, j.get("celldm1_bohr")

    def load_final_cel_pos(cel, pos):
        # last frame
        from analyze_stage import parse_cel, parse_pos
        cs, cm = parse_cel(cel); ps, pp = parse_pos(pos)
        if not len(cm):
            raise ValueError("no cel frames")
        H = cm[-1]  # columns = lattice vectors
        Hrows = np.array(H).T.tolist()
        cart = pp[-1] if len(pp) else None
        if cart is not None:
            frac = (np.array(cart) @ np.linalg.inv(np.array(Hrows))).tolist()
        else:
            frac = None
        return Hrows, frac

    if a.mode == "gen_nvt":
        cell, frac, celldm1 = load_struct(a.struct_json)
        errs, d = validate_struct(cell, frac)
        if errs:
            print("handoff fail", errs); sys.exit(1)
        control = "  nstep = 5000,\n  dt = 20.6706866683,\n  iprint = 1,\n  tstress = .true.,\n  tprnfor = .true.,\n  disk_io = 'low',\n  verbosity = 'high',\n"
        ions = "&IONS\n  ion_dynamics = 'verlet',\n  ion_temperature = 'nose',\n  tempw = 300.0,\n  fnosep(1) = 5.0,\n  nhpcl = 3,\n/\n"
        inp = DFT_HEAD.format(calc="md", prefix="diamond8_nvt", pseudo_dir=a.pseudo_dir,
                              control_extra=control, system=system_ibrav1(celldm1))
        inp += ions + species_cards([(s, x) for s, x in frac], ibrav=1)
        open(a.out_input, "w").write(inp)
        json.dump({"transition": "01_to_02", "n_atoms": len(frac), "determinant_bohr3": d,
                   "celldm1_bohr": celldm1, "velocities_transferred": False, "ok": True},
                  open(a.out_handoff, "w"), indent=2)
        print("wrote NVT input")

    elif a.mode == "gen_npt":
        cell, frac0, celldm1 = load_struct(a.struct_json)   # relaxed cubic cell (fixed by NVT)
        # positions from JOB2 final .pos; cell stays the vc-relaxed cubic (NVT did not change it)
        _, frac2 = load_final_cel_pos(a.prev_cel, a.prev_pos)
        frac = [("C", list(x)) for x in frac2] if frac2 is not None else frac0
        errs, d = validate_struct(cell, frac)
        if errs:
            print("handoff fail", errs); sys.exit(1)
        control = "  nstep = 5000,\n  dt = 20.6706866683,\n  iprint = 1,\n  tstress = .true.,\n  tprnfor = .true.,\n  disk_io = 'low',\n  verbosity = 'high',\n"
        ions = "&IONS\n  ion_dynamics = 'beeman',\n  ion_temperature = 'nose',\n  tempw = 300.0,\n  fnosep(1) = 5.0,\n  nhpcl = 3,\n/\n"
        cellnml = ("&CELL\n  cell_dynamics = 'pr',\n  cell_temperature = 'nose',\n  temph = 300.0,\n"
                   "  fnoseh = 1.0,\n  press = 0.0,\n  cell_dofree = 'volume',\n  caimd = .false.\n/\n")
        inp = DFT_HEAD.format(calc="vc-md", prefix="diamond8_npt", pseudo_dir=a.pseudo_dir,
                              control_extra=control, system=system_ibrav1(celldm1))
        inp += ions + cellnml + species_cards(frac, ibrav=1)
        open(a.out_input, "w").write(inp)
        json.dump({"transition": "02_to_03", "n_atoms": len(frac), "determinant_bohr3": d,
                   "celldm1_bohr": celldm1, "velocities_transferred": False, "ok": True},
                  open(a.out_handoff, "w"), indent=2)
        print("wrote NPT input")

    elif a.mode in ("gen_caimd_zero", "gen_caimd_nonzero"):
        Hrows, frac2 = load_final_cel_pos(a.prev_cel, a.prev_pos)
        frac = [("C", list(x)) for x in frac2]
        errs, d = validate_struct(Hrows, frac)
        if errs:
            print("handoff fail", errs); sys.exit(1)
        if a.mode == "gen_caimd_zero":
            nstep, caimd = 3000, ("  caimd = .true.,\n  caimd_stride = 1,\n"
                                  + "".join(f"  caimd_mask({i}) = .true.,\n" for i in range(1, 7))
                                  + "".join(f"  caimd_step_strain({i}) = 0.0,\n" for i in range(1, 7)))
            prefix = "diamond8_caimd_zero"; trans = "03_to_04"
        else:
            nstep = 2000
            caimd = ("  caimd = .true.,\n  caimd_stride = 10,\n"
                     + "".join(f"  caimd_mask({i}) = .true.,\n" for i in range(1, 7))
                     + "  caimd_step_strain(1) = 5.0d-6,\n"
                     + "".join(f"  caimd_step_strain({i}) = 0.0,\n" for i in range(2, 7)))
            prefix = "diamond8_caimd_nonzero"; trans = "zero_to_nonzero"
        control = f"  nstep = {nstep},\n  dt = 20.6706866683,\n  iprint = 1,\n  tstress = .true.,\n  tprnfor = .true.,\n  disk_io = 'low',\n  verbosity = 'high',\n"
        ions = "&IONS\n  ion_dynamics = 'beeman',\n  ion_temperature = 'nose',\n  tempw = 300.0,\n  fnosep(1) = 5.0,\n  nhpcl = 3,\n/\n"
        cellnml = ("&CELL\n  cell_dynamics = 'pr',\n  cell_temperature = 'not_controlled',\n  press = 0.0,\n"
                   "  cell_dofree = 'all',\n" + caimd + "/\n")
        inp = DFT_HEAD.format(calc="vc-md", prefix=prefix, pseudo_dir=a.pseudo_dir,
                              control_extra=control, system=system_ibrav0())
        inp += ions + cellnml + species_cards(frac, cell_bohr=Hrows, ibrav=0)
        open(a.out_input, "w").write(inp)
        json.dump({"transition": trans, "n_atoms": len(frac), "determinant_bohr3": d,
                   "velocities_transferred": False, "ok": True}, open(a.out_handoff, "w"), indent=2)
        print(f"wrote {a.mode} input")
    else:
        print("unknown mode"); sys.exit(3)


if __name__ == "__main__":
    main()
