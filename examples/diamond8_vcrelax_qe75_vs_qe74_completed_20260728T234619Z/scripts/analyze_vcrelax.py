#!/usr/bin/env python3
"""Analyze a Diamond8 vc-relax stdout and apply the scientific pass gate.

Writes, under --outdir:
  analysis_summary.json, final_structure.json, final_structure.in.fragment,
  pressure_history.csv, volume_history.csv, force_history.csv
and a case report at --report.

Exit code 0 iff the scientific gate passes; nonzero otherwise (so the Slurm
job exits 0 only when the vc-relax gate passes).

Gate (all must hold):
  QE return code == 0; 'JOB DONE' present; BFGS explicitly converged;
  exactly 8 C atoms; final cell finite; final cell determinant > 0;
  final positions finite; final max atomic force <= forc_conv_thr;
  |final scalar pressure| <= 5 kbar; no NaN/Inf, singular cell, SCF failure,
  or fatal QE error.

Pure stdlib (no numpy needed); 3x3 inverse implemented locally.
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


def fnum(s):
    return float(s.replace("D", "E").replace("d", "e"))


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def inv3(m):
    d = det3(m)
    if d == 0:
        raise ValueError("singular matrix")
    c = [[(m[1][1] * m[2][2] - m[1][2] * m[2][1]), -(m[0][1] * m[2][2] - m[0][2] * m[2][1]),
          (m[0][1] * m[1][2] - m[0][2] * m[1][1])],
         [-(m[1][0] * m[2][2] - m[1][2] * m[2][0]), (m[0][0] * m[2][2] - m[0][2] * m[2][0]),
          -(m[0][0] * m[1][2] - m[0][2] * m[1][0])],
         [(m[1][0] * m[2][1] - m[1][1] * m[2][0]), -(m[0][0] * m[2][1] - m[0][1] * m[2][0]),
          (m[0][0] * m[1][1] - m[0][1] * m[1][0])]]
    return [[c[i][j] / d for j in range(3)] for i in range(3)]


def matvec(M, v):
    return [sum(M[i][k] * v[k] for k in range(3)) for i in range(3)]


def parse_last_cell_positions(out_txt):
    """Return (cell_bohr_rows, alat_bohr, [(sym,[x,y,z]),...] fractional)."""
    lines = out_txt.splitlines()
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
        v = [fnum(x) for x in NUM.findall(lines[ci + r])]
        cell.append(v[:3])
    alat = fnum(m.group(1)) if m else None
    if unit == "alat":
        if alat is None:
            raise ValueError("alat cell without alat value")
        cell = [[c * alat for c in row] for row in cell]
    elif unit == "angstrom":
        cell = [[c / BOHR_TO_ANG for c in row] for row in cell]
    pi = pos_idx[-1]
    hp = lines[pi].strip()
    punit = ("crystal" if "crystal" in hp else "bohr" if "bohr" in hp
             else "angstrom" if "angstrom" in hp else "alat")
    raw = []
    j = pi + 1
    while j < len(lines):
        f = lines[j].split()
        if len(f) >= 4 and re.match(r"^[A-Za-z]{1,2}$", f[0]):
            raw.append((f[0], [fnum(x) for x in f[1:4]]))
            j += 1
        else:
            break
    Hinv = inv3(cell)  # cart(row) = frac @ H  => frac = cart @ H^-1 ; here use column form
    frac = []
    for sym, xyz in raw:
        if punit == "crystal":
            frac.append((sym, list(xyz)))
        else:
            if punit == "alat":
                cart = [c * alat for c in xyz]
            elif punit == "angstrom":
                cart = [c / BOHR_TO_ANG for c in xyz]
            else:
                cart = list(xyz)
            # frac_k = sum_j cart_j * Hinv[j][k]  (rows of H are lattice vectors)
            fr = [sum(cart[j] * Hinv[j][k] for j in range(3)) for k in range(3)]
            frac.append((sym, fr))
    return cell, alat, frac


def parse_last_forces(out_txt):
    """Return max atomic force magnitude (Ry/au) from the last force block, or None."""
    blocks = [m.start() for m in re.finditer(r"Forces acting on atoms", out_txt)]
    if not blocks:
        return None
    tail = out_txt[blocks[-1]:]
    fmax = None
    for m in re.finditer(r"force\s*=\s*([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)", tail):
        fx, fy, fz = fnum(m.group(1)), fnum(m.group(2)), fnum(m.group(3))
        mag = math.sqrt(fx * fx + fy * fy + fz * fz)
        fmax = mag if fmax is None else max(fmax, mag)
        # stop at end of this block (Total force line follows)
    return fmax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qe-out", required=True)
    ap.add_argument("--qe-rc", type=int, required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--version-label", default="")
    ap.add_argument("--forc-conv-thr", type=float, default=1.0e-4)
    ap.add_argument("--press-tol-kbar", type=float, default=5.0)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    out = text(a.qe_out)
    fails = []

    job_done = "JOB DONE" in out
    converged = ("End of BFGS Geometry Optimization" in out) or ("bfgs converged" in out.lower())
    err = re.search(r"Error in routine\s+(\S+)", out)
    crash = "%%%%%%%%%%%%%%" in out
    scf_fail = "convergence NOT achieved" in out
    nan_inf = bool(re.search(r"\bNaN\b", out)) or bool(re.search(r"\bInfinity\b", out)) \
        or ("*****" in out.replace("*****", "*****", 1) and re.search(r"[-+]?\*{4,}", out) is not None)
    banner = re.search(r"Program PWSCF\s+v\.?\s*([0-9.]+)", out)
    version = banner.group(1) if banner else None

    pressures = [fnum(x) for x in re.findall(r"P=\s*([-+0-9.EeDd]+)", out)]
    volumes = [fnum(x) for x in re.findall(r"unit-cell volume\s*=\s*([0-9.EeDd+-]+)", out)]
    total_forces = [fnum(x) for x in re.findall(r"Total force\s*=\s*([0-9.EeDd+-]+)", out)]
    energies = [fnum(x) for x in re.findall(r"!\s+total energy\s*=\s*([-0-9.EeDd]+)", out)]
    bfgs_steps = None
    mbf = re.search(r"bfgs converged in\s+(\d+)\s+scf cycles and\s+(\d+)\s+bfgs steps", out)
    if mbf:
        bfgs_steps = int(mbf.group(2))
    n_cell_blocks = len(re.findall(r"^\s*CELL_PARAMETERS", out, re.M))

    fmax = parse_last_forces(out)
    final_P = pressures[-1] if pressures else None
    final_E = energies[-1] if energies else None
    final_V = volumes[-1] if volumes else None

    cell = frac = alat = None
    det = None
    n_atoms = None
    try:
        cell, alat, frac = parse_last_cell_positions(out)
        det = det3(cell)
        n_atoms = len(frac)
    except Exception as e:
        fails.append(f"structure parse error: {e}")

    celldm1 = (abs(det)) ** (1.0 / 3.0) if det else None  # cubic: a = V^(1/3)

    # ---- gate ----
    if a.qe_rc != 0:
        fails.append(f"QE return code {a.qe_rc} != 0")
    if not job_done:
        fails.append("no JOB DONE")
    if not converged:
        fails.append("BFGS did not explicitly converge")
    if err:
        fails.append(f"fatal QE error in routine {err.group(1)}")
    if crash:
        fails.append("QE CRASH block present")
    if scf_fail:
        fails.append("SCF convergence NOT achieved")
    if nan_inf:
        fails.append("NaN/Inf detected")
    if n_atoms is not None and n_atoms != 8:
        fails.append(f"atom count {n_atoms} != 8")
    if cell is not None:
        if any((not math.isfinite(x)) for row in cell for x in row):
            fails.append("nonfinite cell value")
        if not (det is not None and det > 0):
            fails.append("cell determinant not positive")
    if frac is not None and any((not math.isfinite(x)) for _, xyz in frac for x in xyz):
        fails.append("nonfinite atomic position")
    if fmax is None:
        fails.append("no force block parsed")
    elif fmax > a.forc_conv_thr:
        fails.append(f"final max force {fmax:.3e} > forc_conv_thr {a.forc_conv_thr:.3e}")
    if final_P is None:
        fails.append("no pressure parsed")
    elif abs(final_P) > a.press_tol_kbar:
        fails.append(f"final |P| {abs(final_P):.3f} kbar > {a.press_tol_kbar} kbar")

    gate_pass = not fails

    summary = {
        "case": a.case, "version_label": a.version_label, "qe_version_banner": version,
        "qe_return_code": a.qe_rc, "job_done": job_done, "bfgs_converged": converged,
        "bfgs_steps": bfgs_steps, "n_ionic_cell_blocks": n_cell_blocks,
        "n_atoms": n_atoms, "determinant_bohr3": det, "final_celldm1_bohr": celldm1,
        "final_lattice_const_ang": (celldm1 * BOHR_TO_ANG) if celldm1 else None,
        "final_volume_bohr3": final_V, "final_volume_ang3": (final_V * BOHR_TO_ANG ** 3) if final_V else None,
        "final_pressure_kbar": final_P, "final_max_atomic_force_Ry_au": fmax,
        "final_total_energy_Ry": final_E,
        "final_total_energy_Ry_per_atom": (final_E / n_atoms) if (final_E is not None and n_atoms) else None,
        "n_pressure_samples": len(pressures), "n_volume_samples": len(volumes),
        "n_scf_energy_samples": len(energies),
        "forc_conv_thr": a.forc_conv_thr, "press_tol_kbar": a.press_tol_kbar,
        "fatal_error": err.group(1) if err else None, "scf_failure": scf_fail,
        "nan_inf": nan_inf, "crash_block": crash,
        "gate_failures": fails, "gate_pass": gate_pass,
    }
    json.dump(summary, open(os.path.join(a.outdir, "analysis_summary.json"), "w"), indent=2)

    # structure json + fragment
    if frac and det:
        json.dump({"cell_bohr_rows": cell,
                   "fractional_positions": [[s, list(x)] for s, x in frac],
                   "celldm1_bohr": celldm1, "determinant_bohr3": det, "n_atoms": len(frac),
                   "source": f"{a.case} vc-relax final (cell_dofree='ibrav')"},
                  open(os.path.join(a.outdir, "final_structure.json"), "w"), indent=2)
        with open(os.path.join(a.outdir, "final_structure.in.fragment"), "w") as fh:
            fh.write(f"! {a.case} relaxed structure fragment (ibrav=1, cubic)\n")
            fh.write(f"  celldm(1) = {celldm1:.10f},\n")
            fh.write("ATOMIC_POSITIONS crystal\n")
            for sym, xyz in frac:
                fh.write("C  %.12f  %.12f  %.12f\n" % tuple(xyz))
    else:
        # still emit stubs so required outputs always exist
        json.dump({"error": "structure not parsed", "gate_failures": fails},
                  open(os.path.join(a.outdir, "final_structure.json"), "w"), indent=2)
        open(os.path.join(a.outdir, "final_structure.in.fragment"), "w").write(
            "! structure not parsed (see analysis_summary.json)\n")

    with open(os.path.join(a.outdir, "pressure_history.csv"), "w") as fh:
        fh.write("index,P_kbar\n")
        for i, p in enumerate(pressures):
            fh.write(f"{i},{p}\n")
    with open(os.path.join(a.outdir, "volume_history.csv"), "w") as fh:
        fh.write("index,volume_bohr3\n")
        for i, v in enumerate(volumes):
            fh.write(f"{i},{v}\n")
    with open(os.path.join(a.outdir, "force_history.csv"), "w") as fh:
        fh.write("index,total_force_Ry_au\n")
        for i, f in enumerate(total_forces):
            fh.write(f"{i},{f}\n")

    with open(a.report, "w") as fh:
        fh.write(f"# vc-relax report: {a.case}\n\n")
        fh.write(f"Version: {a.version_label} (banner v.{version})\n\n")
        fh.write(f"**Gate: {'PASS' if gate_pass else 'FAIL'}**\n\n")
        fh.write(f"- QE return code: {a.qe_rc}; JOB DONE: {job_done}; BFGS converged: {converged} "
                 f"(bfgs steps: {bfgs_steps})\n")
        fh.write(f"- atoms: {n_atoms}; determinant: {det} bohr^3; det>0: {bool(det and det>0)}\n")
        fh.write(f"- final celldm(1): {celldm1} bohr ({(celldm1*BOHR_TO_ANG) if celldm1 else None} Ang)\n")
        fh.write(f"- final volume: {final_V} bohr^3; final pressure: {final_P} kbar "
                 f"(|P|<= {a.press_tol_kbar})\n")
        fh.write(f"- final max atomic force: {fmax} Ry/au (thr {a.forc_conv_thr})\n")
        fh.write(f"- final total energy: {final_E} Ry "
                 f"({(final_E/n_atoms) if (final_E is not None and n_atoms) else None} Ry/atom)\n")
        fh.write(f"- gate failures: {fails or 'none'}\n")

    print(json.dumps({"case": a.case, "gate_pass": gate_pass, "failures": fails,
                      "final_P_kbar": final_P, "final_celldm1_bohr": celldm1,
                      "final_max_force": fmax}, indent=2))
    sys.exit(0 if gate_pass else 1)


if __name__ == "__main__":
    main()
