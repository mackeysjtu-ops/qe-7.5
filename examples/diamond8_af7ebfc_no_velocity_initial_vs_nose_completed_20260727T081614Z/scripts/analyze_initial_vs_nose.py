#!/usr/bin/env python3
"""Post-processing for the single-factor thermostat isolation (CONTROL initial vs TEST nose).

Runs inside the Slurm job after both QE calculations finish or terminate.
Standard library + numpy/matplotlib only; no network access.

Column semantics verified against PW/src/vcsmd.f90 (af7ebfc, lines 410-419):
  .e    : ut, ekint, edyn, pv, nst
  .eal  : uta, eka, eta, utl, ekla, etl, nst
  .ave  : avu, avk, nst
  .p    : press, p, avp, nst           [Ry/bohr^3]
  .avec : |a|,|b|,|c|, theta12, theta23, theta31, nst
  .tv   : vcell, tnew, nst
  .md   : TERMINAL state only (seqopn rewrite): rati, ratd, rat2d, rat2di, tauold /
          aveci, avecd, avec2d, avec2di / ...

Temperature source rules enforced here:
  * instantaneous ionic T  = 2*eka/((3N-3) k_B)      <- primary metric
  * direct T               = from .md ratd x cell     <- physical confirmation
  * .tv / stdout T         = 2*avk/(3(N+1) k_B)       <- RUNNING AVERAGE, never used
                                                         as the instantaneous value
Energy caveat: ions_nose_energy / cell_nose_energy are computed in vcsmd but are
NOT written to any output file, so the full Nose+PR conserved Hamiltonian is
unavailable. `edyn` is therefore reported as a PARTIAL extended-energy-like sum;
for the Nose run it is expected to change because the chain exchanges energy.
"""
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ANA = os.path.join(ROOT, "analysis")
EVD = os.path.join(ROOT, "evidence")
os.makedirs(ANA, exist_ok=True)

RY_TO_K = 157887.51240116
RY_BOHR3_TO_KBAR = 147105.078
BOHR_TO_ANG = 0.529177210903
NAT = 8
NDOF = 3 * NAT - 3
DT_FS = 2.0
MASS_C_AMU = 12.0107
AMU_RY = 911.4442421
TARGET_T = 300.0

RUNS = {
    "control": {"dir": "control_initial", "prefix": "diamond8_control_initial",
                "label": "CONTROL  ion_temperature='initial'"},
    "nose": {"dir": "test_nose", "prefix": "diamond8_test_nose",
             "label": "TEST  ion_temperature='nose' (fnosep=50, nhpcl=4)"},
}

# reference values verified from the audited completed workflow before submission
REFERENCE = {
    "control_previous_no_velocity_initial": {
        "mean_K": 170.01436653663114, "final200_K": 174.82712827452696,
        "volume_change_percent": -4.622395517125142, "fraction_above_600K": 0.0},
    "nose_previous_with_explicit_velocities_job50127909": {
        "max_K": 8123.838804792437, "final_K": 4231.686070493259,
        "mean_K": 4201.45228096483, "first_step_above_600K": 26,
        "first_step_above_2000K": 130},
}

NUM = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")


def rd(p):
    out = []
    if not os.path.exists(p):
        return np.zeros((0, 0))
    for line in open(p, errors="replace"):
        v = [float(t.replace("D", "E").replace("d", "e")) for t in NUM.findall(line)]
        if v:
            out.append(v)
    if not out:
        return np.zeros((0, 0))
    n = min(len(r) for r in out)
    return np.array([r[:n] for r in out])


def text(p):
    return open(p, errors="replace").read() if os.path.exists(p) else ""


def roll(a, w):
    if len(a) < w:
        return np.full(len(a), np.nan)
    return np.concatenate([np.full(w - 1, np.nan), np.convolve(a, np.ones(w) / w, "valid")])


def blocks(a, w=100):
    return np.array([a[i:i + w].mean() for i in range(0, len(a), w) if len(a[i:i + w])])


def slope(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    return float(np.polyfit(x[m], y[m], 1)[0]) if m.sum() >= 3 else None


def parse_cel(p):
    steps, mats = [], []
    toks = [l.split() for l in text(p).splitlines() if l.strip()]
    i = 0
    while i + 3 <= len(toks):
        if len(toks[i]) == 2:
            try:
                s = int(float(toks[i][0]))
                m = [[float(x) for x in toks[i + r]] for r in (1, 2, 3)]
            except (ValueError, IndexError):
                i += 1
                continue
            steps.append(s); mats.append(m); i += 4
        else:
            i += 1
    return np.array(steps), np.array(mats)


def parse_pos(p, nat=NAT):
    steps, pos = [], []
    toks = [l.split() for l in text(p).splitlines() if l.strip()]
    i = 0
    while i < len(toks):
        if len(toks[i]) == 2:
            try:
                s = int(float(toks[i][0]))
                blk = [[float(x) for x in toks[i + 1 + k]] for k in range(nat)]
            except (ValueError, IndexError):
                i += 1
                continue
            steps.append(s); pos.append(blk); i += 1 + nat
        else:
            i += 1
    return np.array(steps), np.array(pos)


def parse_stress(p):
    lines = text(p).splitlines()
    P, comp = [], []
    for i, line in enumerate(lines):
        if "stress" in line and "P=" in line:
            m = re.search(r"P\s*=\s*([-+0-9.EeDd]+)", line)
            P.append(float(m.group(1).replace("D", "E")) if m else np.nan)
            mat, ok = [], True
            for j in range(1, 4):
                if i + j >= len(lines):
                    ok = False; break
                v = [float(t.replace("D", "E")) for t in NUM.findall(lines[i + j])]
                if len(v) < 6:
                    ok = False; break
                mat.append(v[3:6])
            comp.append(mat if ok else [[np.nan] * 3] * 3)
    return np.array(P), np.array(comp)


def parse_md(p, nat=NAT):
    if not os.path.exists(p):
        return None
    allnum = [float(t.replace("D", "E")) for t in NUM.findall(text(p))]
    n3 = 3 * nat
    need = n3 * 4 + 3 * nat * 3
    if len(allnum) < need + 9:
        return None
    return {"rati": np.array(allnum[0:n3]).reshape(nat, 3),
            "ratd": np.array(allnum[n3:2 * n3]).reshape(nat, 3),
            "aveci": np.array(allnum[need:need + 9]).reshape(3, 3)}


def det3(m):
    return float(np.linalg.det(m))


def cond_sv(m):
    s = np.linalg.svd(m, compute_uv=False)
    return float(s[0] / s[-1]), s


def strain6(H, H0):
    """.cel stores lattice vectors as COLUMNS -> F = H H0^-1 acts on Cartesian columns."""
    F = H @ np.linalg.inv(H0)
    E = 0.5 * (F + F.T) - np.eye(3)
    return [E[0, 0], E[1, 1], E[2, 2], E[1, 2], E[0, 2], E[0, 1]]


def cc_distances(Hcol, cart):
    """min / mean-nearest-neighbour periodic C-C distance (bohr). Hcol columns = lattice vecs."""
    a = Hcol.T  # rows = lattice vectors
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]) @ a
    dmin, nn = np.inf, []
    for i in range(len(cart)):
        best = np.inf
        for j in range(len(cart)):
            d = cart[j] - cart[i] + shifts
            r = np.sqrt((d * d).sum(axis=1))
            r = r[r > 1e-8]
            if len(r):
                dmin = min(dmin, r.min())
                if j != i:
                    best = min(best, r.min())
        nn.append(best)
    return float(dmin), float(np.mean(nn))


def stats(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    if not len(a):
        return {"n": 0}
    return {"n": int(len(a)), "initial": float(a[0]), "final": float(a[-1]),
            "mean": float(a.mean()), "std": float(a.std()), "median": float(np.median(a)),
            "min": float(a.min()), "max": float(a.max()),
            "first100_mean": float(a[:100].mean()),
            "final100_mean": float(a[-100:].mean()),
            "final200_mean": float(a[-200:].mean())}


def first_cross(t, y, thr):
    idx = np.nonzero(y > thr)[0]
    if not len(idx):
        return None
    i = int(idx[0])
    return {"nst": i, "time_ps": float(t[i]), "value": float(y[i])}


def longest_run(y, thr):
    best = cur = 0
    for v in y:
        cur = cur + 1 if v > thr else 0
        best = max(best, cur)
    return int(best)


def wcsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join("" if v is None or (isinstance(v, float) and not np.isfinite(v))
                              else (f"{v:.10g}" if isinstance(v, float) else str(v))
                              for v in r) + "\n")


# ------------------------------------------------------------------ load ----
rc_txt = text(os.path.join(EVD, "run_return_codes.txt"))
RCS = dict(re.findall(r"(\w+)=(-?\d+)", rc_txt))

D = {}
for tag, meta in RUNS.items():
    base = os.path.join(ROOT, meta["dir"])
    pref = meta["prefix"]
    out = os.path.join(base, "qe.out")
    eal = rd(os.path.join(base, pref + ".eal"))
    pf = rd(os.path.join(base, pref + ".p"))
    tv = rd(os.path.join(base, pref + ".tv"))
    ef = rd(os.path.join(base, pref + ".e"))
    av = rd(os.path.join(base, pref + ".ave"))
    cs, cm = parse_cel(os.path.join(base, "tmp", pref + ".cel"))
    ps, pp = parse_pos(os.path.join(base, "tmp", pref + ".pos"))
    P, S = parse_stress(out)
    md = parse_md(os.path.join(base, "tmp", pref + ".md"))
    txt = text(out)
    n = len(eal)
    d = {"label": meta["label"], "n": n, "job_done": "JOB DONE" in txt,
         "version": (re.search(r"Program PWSCF v\.(\S+)", txt) or [None, None])[1]
         if re.search(r"Program PWSCF v\.(\S+)", txt) else None,
         "error": (re.search(r"Error in routine\s+(\S+)", txt) or [None, None])[1]
         if re.search(r"Error in routine\s+(\S+)", txt) else None,
         "scf_fail": txt.count("convergence NOT achieved"),
         "rc": int(RCS.get("RC_CONTROL" if tag == "control" else "RC_NOSE", -1)),
         "cel_steps": cs, "cel": cm, "pos_steps": ps, "pos": pp,
         "P_kbar_stdout": P, "S": S, "md": md}
    if n:
        d["nst"] = eal[:, 6].astype(int)
        d["t_ps"] = d["nst"] * DT_FS / 1000.0
        d["eka"] = eal[:, 1]
        d["ekla"] = eal[:, 4]
        d["T"] = 2 * d["eka"] * RY_TO_K / NDOF
        d["T_3Np1"] = 2 * d["eka"] * RY_TO_K / (3 * (NAT + 1))
        d["tv_T"] = tv[:, 1] if len(tv) else np.full(n, np.nan)
        d["vcell"] = tv[:, 0] if len(tv) else np.full(n, np.nan)
        d["avk"] = av[:, 1] if len(av) else np.full(n, np.nan)
        d["T_runavg"] = 2 * d["avk"] * RY_TO_K / (3 * (NAT + 1))
        d["p_kbar"] = pf[:, 1] * RY_BOHR3_TO_KBAR if len(pf) else np.full(n, np.nan)
        d["edyn"] = ef[:, 2] if len(ef) else np.full(n, np.nan)
        d["ut"] = ef[:, 0] if len(ef) else np.full(n, np.nan)
        d["ekint"] = ef[:, 1] if len(ef) else np.full(n, np.nan)
    D[tag] = d

# ------------------------------------------------- initial-state comparison --
c, s = D["control"], D["nose"]
init = {"note": ("Both runs omit explicit velocities and use the same af7ebfc executable, so "
                 "vcinit->ranv(temp=tempw) should generate the SAME deterministic Maxwell "
                 "realization (QE randy is seeded with idum=0). ranv also zeroes the total "
                 "linear momentum, so the initial COM velocity is zero by construction. "
                 "tmp/*.md holds only the TERMINAL state, so the first committed ratd cannot "
                 "be read back; step-1 kinetic energy and the step-2 positions are used as the "
                 "equivalent evidence.")}
if c["n"] and s["n"]:
    init["first_step_eka_Ry"] = {"control": float(c["eka"][0]), "nose": float(s["eka"][0])}
    init["first_step_eka_abs_diff_Ry"] = abs(float(c["eka"][0] - s["eka"][0]))
    init["first_step_T_K"] = {"control": float(c["T"][0]), "nose": float(s["T"][0])}
    init["first_step_T_abs_diff_K"] = abs(float(c["T"][0] - s["T"][0]))
    init["identical_initial_kinetic_energy"] = bool(c["eka"][0] == s["eka"][0])
    # step-2 positions reflect v*dt from the initial velocities
    if len(c["pos"]) > 1 and len(s["pos"]) > 1:
        dp = np.abs(c["pos"][1] - s["pos"][1]).max()
        init["step2_position_max_abs_diff_bohr"] = float(dp)
        init["identical_step2_positions"] = bool(dp == 0.0)
    for tag, dd in (("control", c), ("nose", s)):
        if dd["md"] is not None and len(dd["cel"]):
            v = dd["md"]["ratd"] @ dd["cel"][-1].T
            init.setdefault("terminal_COM_velocity_au", {})[tag] = [float(x) for x in v.mean(axis=0)]
    init["strict_single_factor"] = bool(init.get("identical_initial_kinetic_energy")
                                        and init.get("identical_step2_positions", False))
    init["interpretation"] = ("Identical initial realization: this is an especially strict "
                              "single-factor comparison (thermostat mode is the only variable)."
                              if init["strict_single_factor"] else
                              "Initial realizations differ: the random realization is an "
                              "additional uncontrolled difference; compare statistically.")
json.dump(init, open(os.path.join(ANA, "initial_state_comparison.json"), "w"), indent=2)

# --------------------------------------------------------- per-run analysis --
tcross, blk_rows, energy_rows, runstat_rows = [], [], [], []
RES = {}
for tag, d in D.items():
    if not d["n"]:
        RES[tag] = {"label": d["label"], "classification": "INCOMPLETE", "n": 0,
                    "rc": d["rc"], "job_done": d["job_done"], "error": d["error"]}
        continue
    t, T, P = d["t_ps"], d["T"], d["p_kbar"]

    wcsv(os.path.join(ANA, f"{tag}_temperature.csv"),
         ["nst", "time_ps", "eka_Ry", "T_instant_3Nm3_K", "T_instant_3Np1_K",
          "T_running_avg_reported_K", "roll25_K", "roll100_K"],
         [[int(d["nst"][i]), t[i], d["eka"][i], T[i], d["T_3Np1"][i], d["tv_T"][i],
           roll(T, 25)[i], roll(T, 100)[i]] for i in range(d["n"])])

    wcsv(os.path.join(ANA, f"{tag}_pressure_stress.csv"),
         ["nst", "time_ps", "p_kbar", "p_GPa", "stdout_P_kbar",
          "sxx_kbar", "syy_kbar", "szz_kbar", "syz_kbar", "sxz_kbar", "sxy_kbar"],
         [[int(d["nst"][i]), t[i], P[i], P[i] / 10.0,
           d["P_kbar_stdout"][i] if i < len(d["P_kbar_stdout"]) else None]
          + ([d["S"][i][0][0], d["S"][i][1][1], d["S"][i][2][2],
              d["S"][i][1][2], d["S"][i][0][2], d["S"][i][0][1]]
             if i < len(d["S"]) else [None] * 6) for i in range(d["n"])])

    # cell
    cel_rows = []
    vol = np.array([abs(det3(m)) for m in d["cel"]]) if len(d["cel"]) else np.array([])
    conds = np.array([cond_sv(m)[0] for m in d["cel"]]) if len(d["cel"]) else np.array([])
    dmins, nns = [], []
    if len(d["cel"]):
        H0 = d["cel"][0]
        sample = range(0, len(d["cel"]), max(1, len(d["cel"]) // 200))
        for i in range(len(d["cel"])):
            e6 = strain6(d["cel"][i], H0)
            dmin = nn = None
            if i in sample and i < len(d["pos"]):
                dmin, nn = cc_distances(d["cel"][i], d["pos"][i])
                dmins.append(dmin); nns.append(nn)
            sv = cond_sv(d["cel"][i])[1]
            cel_rows.append([int(d["cel_steps"][i]), (d["cel_steps"][i] - d["cel_steps"][0]) * DT_FS / 1000.0,
                             vol[i], 100 * (vol[i] - vol[0]) / vol[0], conds[i],
                             float(sv[0]), float(sv[1]), float(sv[2])] + [float(x) for x in e6]
                            + [dmin, nn])
    wcsv(os.path.join(ANA, f"{tag}_cell.csv"),
         ["step", "time_ps", "volume_bohr3", "volume_change_percent", "condition_number",
          "sv1", "sv2", "sv3", "strain_xx", "strain_yy", "strain_zz", "strain_yz",
          "strain_xz", "strain_xy", "min_CC_bohr", "mean_NN_CC_bohr"], cel_rows)

    # temperature source crosscheck
    tcross.append([tag, "T_instant_3Np1 vs reported .tv (lag 0)", 0,
                   float(np.nanmax(np.abs(d["T_3Np1"] - d["tv_T"])))])
    tcross.append([tag, "T_running_avg(2*avk/3(N+1)) vs reported .tv (lag 0)", 0,
                   float(np.nanmax(np.abs(d["T_runavg"] - d["tv_T"])))])
    if d["md"] is not None and len(d["cel"]):
        v = d["md"]["ratd"] @ d["cel"][-1].T
        vcom = v.mean(axis=0)
        K = 0.5 * MASS_C_AMU * AMU_RY * ((v - vcom) ** 2).sum()
        Tdir = 2 * K * RY_TO_K / NDOF
        tcross.append([tag, "direct .md ratd x terminal .cel (3N-3)", 0, float(Tdir)])
        tcross.append([tag, "eal instantaneous terminal (3N-3)", 0, float(T[-1])])
        tcross.append([tag, "direct-vs-eal terminal abs diff", 0, float(abs(Tdir - T[-1]))])
        direct_conf = bool(abs(Tdir - T[-1]) < max(5.0, 0.05 * abs(T[-1])))
    else:
        Tdir, direct_conf = None, False

    st, sp = stats(T), stats(P)
    r25, r100 = roll(T, 25), roll(T, 100)
    bT, bP = blocks(T), blocks(P)
    for i, (bt, bp) in enumerate(zip(bT, bP)):
        blk_rows.append([tag, i + 1, i * 100, min((i + 1) * 100 - 1, d["n"] - 1), bt, bp])

    # energy (PARTIAL for nose: Nose chain terms are not written by QE)
    ed = d["edyn"]
    drift_total = float(ed[-1] - ed[0])
    drift_ps = drift_total / max(t[-1], 1e-12)
    drift_at = drift_ps * 13.605693122994 / NAT
    bounded = bool(np.all(np.isfinite(ed)) and np.abs(ed - ed[0]).max() < 5.0)
    partial = (tag == "nose")
    energy_rows.append([tag, "edyn (.e col3 = ekint + ut)",
                        "PARTIAL - excludes Nose chain energy (not written by QE)"
                        if partial else "extended PR total at press=0",
                        float(ed[0]), float(ed[-1]), drift_total, drift_ps, drift_at, bounded,
                        float(np.corrcoef(ed, T)[0, 1]) if np.all(np.isfinite(ed)) else None])

    # safety condition check (post-hoc; no external killer was armed)
    safety = {
        "nonpositive_determinant": bool(len(vol) and vol.min() <= 0),
        "condition_number_gt_1e8": bool(len(conds) and conds.max() > 1e8),
        "nonfinite_values": bool(not np.all(np.isfinite(T)) or not np.all(np.isfinite(P))
                                 or (len(vol) and not np.all(np.isfinite(vol)))),
        "qe_singular_cell": bool(d["error"] and "singular" in str(d["error"]).lower()),
        "T_above_20000K_for_10_consecutive": bool(longest_run(T, 20000.0) >= 10),
    }
    safety["any_triggered"] = any(safety.values())

    # ---- classification ----
    fin200 = st["final200_mean"]
    bslope = slope(np.arange(len(bT), dtype=float), bT) if len(bT) >= 3 else None
    late_slope = (slope(np.arange(len(bT[len(bT) // 2:]), dtype=float), bT[len(bT) // 2:])
                  if len(bT) >= 6 else None)
    if tag == "control":
        if d["n"] < 50:
            cls = "CONTROL_INCOMPLETE"
        elif safety["nonpositive_determinant"] or safety["qe_singular_cell"] or safety["nonfinite_values"]:
            cls = "CONTROL_CELL_FAILURE"
        elif st["max"] > 2000 and fin200 > 1000:
            cls = "CONTROL_HEATING"
        elif fin200 < 30:
            cls = "CONTROL_COOLING"
        else:
            cls = "CONTROL_BOUNDED"
    else:
        early_max = float(T[:200].max()) if d["n"] >= 200 else st["max"]
        overshoot = early_max > 1.5 * TARGET_T
        near_target = 0.6 * TARGET_T <= fin200 <= 1.4 * TARGET_T
        if d["n"] < 50:
            cls = "INCOMPLETE"
        elif safety["nonfinite_values"] or safety["qe_singular_cell"] or safety["nonpositive_determinant"]:
            cls = "NUMERICAL_FAILURE"
        elif fin200 > 1000 and st["max"] > 2000 and (late_slope is None or late_slope > -0.05 * TARGET_T):
            cls = "SUSTAINED_PHYSICAL_HEATING"
        elif overshoot and near_target and (late_slope is not None and late_slope < 0):
            cls = "INITIAL_TRANSIENT_THEN_CONVERGENCE"
        elif 0.75 * TARGET_T <= fin200 <= 1.25 * TARGET_T:
            cls = "BOUNDED_OSCILLATION_AROUND_TARGET"
        else:
            cls = "BOUNDED_BUT_OFFSET"

    RES[tag] = {
        "label": d["label"], "rc": d["rc"], "job_done": d["job_done"],
        "version": d["version"], "error": d["error"], "scf_failures": d["scf_fail"],
        "completed_steps": d["n"], "simulated_time_ps": float(t[-1]),
        "classification": cls,
        "temperature": st,
        "temperature_initialized_K": float(T[0]),
        "temperature_drift_K_per_ps": slope(t, T),
        "block100_means_K": [float(x) for x in bT],
        "block_slope_K_per_block": bslope,
        "late_half_block_slope_K_per_block": late_slope,
        "roll25_final": float(r25[-1]) if np.isfinite(r25[-1]) else None,
        "roll100_final": float(r100[-1]) if np.isfinite(r100[-1]) else None,
        "crossings": {f"first_above_{k}K": first_cross(t, T, k) for k in (600, 1000, 2000)},
        "fraction_above": {f"{k}K": float(np.mean(T > k)) for k in (600, 1000, 2000)},
        "longest_run_above": {f"{k}K": longest_run(T, k) for k in (600, 1000, 2000)},
        "direct_terminal_T_K": Tdir,
        "direct_confirms_reported": direct_conf,
        "pressure_kbar": sp,
        "pressure_GPa": {k: (v / 10.0 if isinstance(v, float) and k != "n" else v)
                         for k, v in sp.items()},
        "stress_mean_kbar": ({n_: float(np.nanmean(d["S"][:, i, j])) for n_, (i, j) in
                              zip(("sxx", "syy", "szz", "syz", "sxz", "sxy"),
                                  ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)))}
                             if len(d["S"]) else {}),
        "cell": ({"volume_initial_bohr3": float(vol[0]), "volume_final_bohr3": float(vol[-1]),
                  "volume_change_percent": float(100 * (vol[-1] - vol[0]) / vol[0]),
                  "determinant_min": float(vol.min()),
                  "condition_number_max": float(conds.max()),
                  "strain6_final": [float(x) for x in strain6(d["cel"][-1], d["cel"][0])],
                  "min_CC_bohr": float(min(dmins)) if dmins else None,
                  "mean_NN_CC_bohr": float(np.mean(nns)) if nns else None}
                 if len(vol) else {}),
        "energy_partial": {"quantity": "edyn = ekint + ut",
                           "is_full_conserved_hamiltonian": False,
                           "note": ("PARTIAL: Nose chain energy is computed in vcsmd but never "
                                    "written to an output file, so a full Nose+PR conserved "
                                    "quantity cannot be formed. For the Nose run a change in "
                                    "edyn is EXPECTED (thermostat exchanges energy)."),
                           "initial_Ry": float(ed[0]), "final_Ry": float(ed[-1]),
                           "drift_total_Ry": drift_total, "drift_Ry_per_ps": drift_ps,
                           "drift_eV_per_atom_per_ps": drift_at, "bounded": bounded},
        "safety_conditions": safety,
    }
    runstat_rows.append([tag, d["label"], cls, d["rc"], d["job_done"], d["n"],
                         st["initial"], st["mean"], st["std"], st["min"], st["max"],
                         st["final"], st["first100_mean"], st["final100_mean"],
                         st["final200_mean"], RES[tag]["temperature_drift_K_per_ps"],
                         sp["mean"], sp["final200_mean"],
                         RES[tag]["cell"].get("volume_change_percent"),
                         RES[tag]["cell"].get("condition_number_max"), drift_at])

wcsv(os.path.join(ANA, "temperature_source_crosscheck.csv"),
     ["run", "comparison", "lag_steps", "value_or_max_abs_diff_K"], tcross)
wcsv(os.path.join(ANA, "block_statistics.csv"),
     ["run", "block", "nst_start", "nst_end", "temperature_mean_K", "pressure_mean_kbar"], blk_rows)
wcsv(os.path.join(ANA, "energy_diagnostics.csv"),
     ["run", "quantity", "status", "initial_Ry", "final_Ry", "drift_total_Ry",
      "drift_Ry_per_ps", "drift_eV_per_atom_per_ps", "bounded", "corr_with_temperature"],
     energy_rows)
wcsv(os.path.join(ANA, "run_statistics.csv"),
     ["run", "label", "classification", "return_code", "job_done", "steps", "initial_K",
      "mean_K", "std_K", "min_K", "max_K", "final_K", "first100_mean_K", "final100_mean_K",
      "final200_mean_K", "drift_K_per_ps", "mean_P_kbar", "final200_P_kbar",
      "volume_change_percent", "condition_number_max", "energy_drift_eV_atom_ps"], runstat_rows)

# ----------------------------------------------------------- overall status --
cc_, nn_ = RES["control"], RES["nose"]
ccls, ncls = cc_.get("classification"), nn_.get("classification")
if cc_.get("completed_steps", 0) < 50 or nn_.get("completed_steps", 0) < 50:
    overall = "PARTIAL"
elif ncls == "SUSTAINED_PHYSICAL_HEATING" and not nn_.get("direct_confirms_reported", False):
    overall = "NOSE_DIAGNOSTIC_MISMATCH"
elif ccls != "CONTROL_BOUNDED" and ncls in ("SUSTAINED_PHYSICAL_HEATING", "NUMERICAL_FAILURE"):
    overall = "BOTH_RUNS_UNSTABLE"
elif ccls != "CONTROL_BOUNDED":
    overall = "CONTROL_UNSTABLE_NOSE_INTERPRETATION_BLOCKED"
elif ncls == "SUSTAINED_PHYSICAL_HEATING":
    overall = "NOSE_ALONE_CAUSES_HEATING"
elif ncls == "INITIAL_TRANSIENT_THEN_CONVERGENCE":
    overall = "NOSE_TRANSIENT_THEN_CONVERGES"
elif ncls in ("BOUNDED_OSCILLATION_AROUND_TARGET",):
    overall = "NOSE_STABLE_WITH_NATIVE_INITIALIZATION"
elif ncls == "BOUNDED_BUT_OFFSET":
    overall = "NOSE_STABLE_WITH_NATIVE_INITIALIZATION"
else:
    overall = "PARTIAL"

# heating-vs-cell-deformation ordering
order = {}
if nn_.get("completed_steps"):
    fc = nn_["crossings"]["first_above_600K"]
    if fc and D["nose"]["n"] and len(D["nose"]["cel"]):
        vol = np.array([abs(det3(m)) for m in D["nose"]["cel"]])
        rel = np.abs(vol - vol[0]) / vol[0]
        idx = np.nonzero(rel > 0.02)[0]
        first_def = float(idx[0] * DT_FS / 1000.0) if len(idx) else None
        order = {"first_600K_time_ps": fc["time_ps"],
                 "first_2percent_volume_change_time_ps": first_def,
                 "heating_precedes_cell_deformation":
                     bool(first_def is None or fc["time_ps"] < first_def)}
    elif not fc:
        order = {"first_600K_time_ps": None, "note": "Nose never crossed 600 K"}

summary = {
    "workflow": "diamond8_af7ebfc_no_velocity_initial_vs_nose_300K",
    "purpose": ("strict single-factor isolation of the ionic thermostat mode; same af7ebfc "
                "executable, same structure, no explicit velocities, caimd=.false. in both"),
    "changing_variable": "ion_temperature 'initial' -> 'nose' (+ fnosep=50, nhpcl=4)",
    "caveat_pressure": ("ecutwfc=25 Ry retained from the reference deck; absolute pressures are "
                        "dominated by Pulay stress and are NOT converged diamond pressures. "
                        "Pressure/volume are used only to compare the two runs and to diagnose "
                        "cell stability."),
    "temperature_convention": {
        "instantaneous": "2*eka/((3N-3) k_B), eka from .eal col2",
        "reported_tv_stdout": "2*avk/(3(N+1) k_B) - RUNNING AVERAGE, not used as instantaneous",
        "direct": "from terminal .md ratd x terminal .cel, COM removed"},
    "initial_state": init,
    "runs": RES,
    "overall_status": overall,
    "heating_vs_cell_order": order,
    "reference_values_verified_before_submission": REFERENCE,
    "comparison_with_reference": {
        "control_vs_previous_no_velocity_initial": {
            "previous_mean_K": REFERENCE["control_previous_no_velocity_initial"]["mean_K"],
            "new_mean_K": cc_.get("temperature", {}).get("mean"),
            "previous_final200_K": REFERENCE["control_previous_no_velocity_initial"]["final200_K"],
            "new_final200_K": cc_.get("temperature", {}).get("final200_mean"),
        },
        "nose_vs_previous_explicit_velocity_nose": {
            "previous_max_K": REFERENCE["nose_previous_with_explicit_velocities_job50127909"]["max_K"],
            "new_max_K": nn_.get("temperature", {}).get("max"),
            "previous_final_K": REFERENCE["nose_previous_with_explicit_velocities_job50127909"]["final_K"],
            "new_final_K": nn_.get("temperature", {}).get("final"),
        },
    },
}
json.dump(summary, open(os.path.join(ANA, "comparison_summary.json"), "w"), indent=2)

# ------------------------------------------------------------------ plots ----
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.color": "#e1e0d9",
                         "grid.linewidth": .8, "axes.axisbelow": True,
                         "axes.unicode_minus": False, "savefig.dpi": 150,
                         "figure.facecolor": "white", "axes.facecolor": "white"})
    COL = {"control": "#2a78d6", "nose": "#eb6834"}
    SRCN = ("Source: this workflow; T from .eal eka (3N-3=21 dof); P from .p col2; cell from .cel. "
            "ecutwfc=25 Ry: pressures are not converged.")

    def fig1(fname, xf, yf, title, xl, yl, hline=None, marker=""):
        fig, ax = plt.subplots(figsize=(10, 5))
        drew = False
        for tag in ("control", "nose"):
            d = D[tag]
            if not d["n"]:
                continue
            x, y = xf(tag), yf(tag)
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum():
                ax.plot(x[m], y[m], marker + "-" if marker else "-", lw=1.2 if marker else .9,
                        ms=4, color=COL[tag], label=RES[tag]["label"], alpha=.9)
                drew = True
        if hline is not None:
            ax.axhline(hline, color="#898781", ls="--", lw=1.2)
        ax.set_title(title, fontsize=11.5); ax.set_xlabel(xl); ax.set_ylabel(yl)
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
        if drew:
            ax.legend(frameon=False, fontsize=9)
        fig.text(0.01, -0.02, SRCN, fontsize=7.4, color="#52514e")
        fig.savefig(os.path.join(ANA, fname), bbox_inches="tight"); plt.close(fig)

    fig1("temperature_initial_vs_nose.png", lambda k: D[k]["t_ps"], lambda k: D[k]["T"],
         "Instantaneous ionic temperature: initial vs Nose (no explicit velocities)",
         "time (ps)", "T (K, 3N-3)", hline=TARGET_T)
    fig1("temperature_rolling100_initial_vs_nose.png", lambda k: D[k]["t_ps"],
         lambda k: roll(D[k]["T"], 100),
         "Temperature, 100-step rolling mean", "time (ps)", "rolling mean T (K)", hline=TARGET_T)
    fig1("temperature_block100_initial_vs_nose.png",
         lambda k: np.arange(1, len(blocks(D[k]["T"])) + 1, dtype=float),
         lambda k: blocks(D[k]["T"]),
         "Temperature, 100-step block means", "100-step block", "block-mean T (K)",
         hline=TARGET_T, marker="o")
    fig1("ionic_kinetic_energy_initial_vs_nose.png", lambda k: D[k]["t_ps"],
         lambda k: D[k]["eka"], "Ionic kinetic energy", "time (ps)", "eka (Ry)")
    fig1("pressure_initial_vs_nose.png", lambda k: D[k]["t_ps"], lambda k: D[k]["p_kbar"],
         "MD pressure (target press = 0; NOT a converged pressure)", "time (ps)", "P (kbar)",
         hline=0.0)
    fig1("volume_initial_vs_nose.png",
         lambda k: (D[k]["cel_steps"] - D[k]["cel_steps"][0]) * DT_FS / 1000.0 if len(D[k]["cel"]) else np.array([]),
         lambda k: np.array([abs(det3(m)) for m in D[k]["cel"]]) if len(D[k]["cel"]) else np.array([]),
         "Parrinello-Rahman cell volume", "time (ps)", "volume (bohr^3)")
    fig1("cell_condition_initial_vs_nose.png",
         lambda k: (D[k]["cel_steps"] - D[k]["cel_steps"][0]) * DT_FS / 1000.0 if len(D[k]["cel"]) else np.array([]),
         lambda k: np.array([cond_sv(m)[0] for m in D[k]["cel"]]) if len(D[k]["cel"]) else np.array([]),
         "Cell condition number", "time (ps)", "cond(H)")
    fig1("energy_drift_initial_vs_nose.png", lambda k: D[k]["t_ps"],
         lambda k: (D[k]["edyn"] - D[k]["edyn"][0]) * 13.605693122994 / NAT,
         "Partial extended-energy drift (edyn = ekint + ut; excludes Nose chain)",
         "time (ps)", "drift per atom (eV/atom)", hline=0.0)
    summary["plots_generated"] = True
except Exception as exc:
    summary["plots_generated"] = False
    summary["plot_error"] = repr(exc)
json.dump(summary, open(os.path.join(ANA, "comparison_summary.json"), "w"), indent=2)


# ----------------------------------------------------------------- report ----
def f(x, n=3):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{n}f}"


ct, nt = cc_.get("temperature", {}), nn_.get("temperature", {})
cp, npr = cc_.get("pressure_kbar", {}), nn_.get("pressure_kbar", {})
L = []
L.append("# Final Thermostat Isolation Report")
L.append("")
L.append("Diamond8, QE 7.5 `af7ebfc`, **no explicit velocities**, `caimd=.false.` in both runs.")
L.append("Single changing variable: `ion_temperature='initial'` (CONTROL) vs `'nose'` "
         "with `fnosep=50`, `nhpcl=4` (TEST).")
L.append("")
L.append(f"**Overall status: `{overall}`**")
L.append("")
L.append(f"- CONTROL classification: `{ccls}`")
L.append(f"- NOSE classification: `{ncls}`")
L.append("")
L.append("## Required answers")
L.append("")
L.append("1. **Explicit velocities absent from both inputs?** Yes — see "
         "`evidence/no_explicit_velocity_audit.json`: no `ATOMIC_VELOCITIES` header, no velocity "
         "rows, no explicit-velocity keyword anywhere in the complete text, no restart-state "
         "input, both `restart_mode='from_scratch'`.")
L.append("2. **`caimd=.false.` in both?** Yes — `evidence/caimd_disabled_audit.json`: exactly one "
         "`caimd = .false.` line each, no `caimd_stride`/`caimd_step_strain`/`caimd_mask`, no "
         "`caimd=.true.`.")
L.append("3. **All non-thermostat settings identical?** Yes — `evidence/input_identity_audit.json` "
         "reports maximum difference 0 for cell, positions, ordering, masses, cutoffs, timestep, "
         "`nstep`, `press`, `wmass`, `tempw` and every SCF parameter. Only `prefix`, "
         "`ion_temperature`, `fnosep(1)` and `nhpcl` differ.")
L.append(f"4. **Same native initial velocity realization?** "
         f"{'Yes' if init.get('strict_single_factor') else 'Not identical'} — first-step ionic "
         f"kinetic energy {f(init.get('first_step_eka_Ry',{}).get('control'),8)} Ry (control) vs "
         f"{f(init.get('first_step_eka_Ry',{}).get('nose'),8)} Ry (nose), |diff| = "
         f"{init.get('first_step_eka_abs_diff_Ry')}; step-2 position max |diff| = "
         f"{init.get('step2_position_max_abs_diff_bohr')} bohr. {init.get('interpretation','')}")
L.append(f"5. **Did the control remain bounded?** `{ccls}` — mean {f(ct.get('mean'))} K, "
         f"final-200 {f(ct.get('final200_mean'))} K, max {f(ct.get('max'))} K, "
         f"min {f(ct.get('min'))} K, drift {f(cc_.get('temperature_drift_K_per_ps'))} K/ps.")
L.append(f"6. **Did Nose produce an early overshoot?** max {f(nt.get('max'))} K; first >600 K at "
         f"{nn_.get('crossings',{}).get('first_above_600K')}.")
L.append(f"7. **Did the rolling/block mean later return toward 300 K?** final-200 mean "
         f"{f(nt.get('final200_mean'))} K; 100-step rolling mean at the end "
         f"{f(nn_.get('roll100_final'))} K; late-half block slope "
         f"{f(nn_.get('late_half_block_slope_K_per_block'))} K per 100-step block "
         f"(negative = declining).")
L.append(f"8. **Did Nose cross 600 / 1000 / 2000 K?** "
         f"600 K: {nn_.get('crossings',{}).get('first_above_600K')}; "
         f"1000 K: {nn_.get('crossings',{}).get('first_above_1000K')}; "
         f"2000 K: {nn_.get('crossings',{}).get('first_above_2000K')}.")
L.append(f"9. **How long did high-temperature intervals persist?** longest consecutive run above "
         f"600 K = {nn_.get('longest_run_above',{}).get('600K')} steps, above 1000 K = "
         f"{nn_.get('longest_run_above',{}).get('1000K')}, above 2000 K = "
         f"{nn_.get('longest_run_above',{}).get('2000K')}; fractions "
         f"{nn_.get('fraction_above')}.")
L.append(f"10. **Did direct `.md` velocities confirm the reported temperature?** "
         f"control: direct {f(cc_.get('direct_terminal_T_K'))} K vs `.eal` "
         f"{f(ct.get('final'))} K; nose: direct {f(nn_.get('direct_terminal_T_K'))} K vs "
         f"{f(nt.get('final'))} K. Confirmed: control="
         f"{cc_.get('direct_confirms_reported')}, nose={nn_.get('direct_confirms_reported')}.")
L.append(f"11. **Did heating precede significant cell deformation?** {order}")
L.append(f"12. **Was the extended-energy behavior bounded?** control edyn drift "
         f"{cc_.get('energy_partial',{}).get('drift_eV_per_atom_per_ps'):.3e} eV/atom/ps "
         f"(bounded={cc_.get('energy_partial',{}).get('bounded')}); nose "
         f"{nn_.get('energy_partial',{}).get('drift_eV_per_atom_per_ps'):.3e} eV/atom/ps "
         f"(bounded={nn_.get('energy_partial',{}).get('bounded')}). "
         f"**This is a PARTIAL quantity**: `ions_nose_energy` is computed in `vcsmd` but never "
         f"written to an output file, so the full Nose+PR conserved Hamiltonian cannot be formed; "
         f"a change in `edyn` for the Nose run is expected by construction.")
if overall == "NOSE_ALONE_CAUSES_HEATING":
    L.append("13. **Are explicit velocities required to reproduce the earlier heating?** **No.** "
             "The Nose run heats with native initialization and no velocity card.")
    L.append("14. **Is active CAIMD required?** **No.** `caimd=.false.` in both runs.")
elif overall in ("NOSE_TRANSIENT_THEN_CONVERGES", "NOSE_STABLE_WITH_NATIVE_INITIALIZATION"):
    L.append("13. **Are explicit velocities required to reproduce the earlier heating?** On this "
             "evidence, **plausibly yes** — Nose alone with native initialization did not "
             "reproduce sustained heating, so the explicit-velocity path becomes the prime "
             "remaining suspect. This is an inference from two runs, not yet a proof.")
    L.append("14. **Is active CAIMD required?** No — CAIMD was disabled in both runs here, and in "
             "the earlier heating run as well.")
else:
    L.append("13. **Are explicit velocities required?** Not determined by this comparison.")
    L.append("14. **Is active CAIMD required?** No — CAIMD was disabled in both runs.")
L.append("15. **Narrowest justified next test:** see the recommendation below.")
L.append("")
L.append("## Run comparison")
L.append("")
L.append("| quantity | CONTROL (`initial`) | TEST (`nose`) |")
L.append("|---|---:|---:|")
for lab, ka, kb in (
        ("classification", ccls, ncls),
        ("return code", cc_.get("rc"), nn_.get("rc")),
        ("JOB DONE", cc_.get("job_done"), nn_.get("job_done")),
        ("steps", cc_.get("completed_steps"), nn_.get("completed_steps")),
        ("initialized T (K)", f(cc_.get("temperature_initialized_K"), 4),
         f(nn_.get("temperature_initialized_K"), 4)),
        ("mean T (K)", f(ct.get("mean")), f(nt.get("mean"))),
        ("final-200 mean T (K)", f(ct.get("final200_mean")), f(nt.get("final200_mean"))),
        ("max T (K)", f(ct.get("max")), f(nt.get("max"))),
        ("min T (K)", f(ct.get("min")), f(nt.get("min"))),
        ("T drift (K/ps)", f(cc_.get("temperature_drift_K_per_ps")),
         f(nn_.get("temperature_drift_K_per_ps"))),
        ("mean P (kbar)", f(cp.get("mean"), 2), f(npr.get("mean"), 2)),
        ("final-200 P (kbar)", f(cp.get("final200_mean"), 2), f(npr.get("final200_mean"), 2)),
        ("volume change (%)", f(cc_.get("cell", {}).get("volume_change_percent")),
         f(nn_.get("cell", {}).get("volume_change_percent"))),
        ("max cond(H)", f(cc_.get("cell", {}).get("condition_number_max"), 5),
         f(nn_.get("cell", {}).get("condition_number_max"), 5)),
        ("min C-C (bohr)", f(cc_.get("cell", {}).get("min_CC_bohr")),
         f(nn_.get("cell", {}).get("min_CC_bohr"))),
):
    L.append(f"| {lab} | {ka} | {kb} |")
L.append("")
L.append("## Comparison with previously verified results")
L.append("")
L.append(f"- Previous `initial`/no-velocity control (audited): mean "
         f"{REFERENCE['control_previous_no_velocity_initial']['mean_K']:.3f} K, final-200 "
         f"{REFERENCE['control_previous_no_velocity_initial']['final200_K']:.3f} K, volume change "
         f"{REFERENCE['control_previous_no_velocity_initial']['volume_change_percent']:.3f}%. "
         f"New control: mean {f(ct.get('mean'))} K, final-200 {f(ct.get('final200_mean'))} K, "
         f"volume {f(cc_.get('cell',{}).get('volume_change_percent'))}%.")
L.append(f"- Previous Nose **with explicit velocities** (job 50127909): max "
         f"{REFERENCE['nose_previous_with_explicit_velocities_job50127909']['max_K']:.3f} K, final "
         f"{REFERENCE['nose_previous_with_explicit_velocities_job50127909']['final_K']:.3f} K, "
         f"first >600 K at step "
         f"{REFERENCE['nose_previous_with_explicit_velocities_job50127909']['first_step_above_600K']}. "
         f"New Nose **without** explicit velocities: max {f(nt.get('max'))} K, final "
         f"{f(nt.get('final'))} K.")
L.append("")
L.append("## Safety")
L.append("")
L.append("No external runtime killer was armed: each run is ~2 minutes of wall time and QE's own "
         "singular-cell / non-finite guards are the effective safety mechanism. The extreme "
         "conditions were evaluated post-hoc:")
L.append("")
L.append(f"- CONTROL: {cc_.get('safety_conditions')}")
L.append(f"- TEST: {nn_.get('safety_conditions')}")
L.append("")
L.append("## Recommended minimum next test (not submitted)")
L.append("")
if overall == "NOSE_ALONE_CAUSES_HEATING":
    L.append("Keep no explicit velocities and `caimd=.false.`, and vary **one** Nose parameter at "
             "a time, starting from a conservative point: `fnosep=1.0`, `nhpcl=1`, `dt` = 1 fs. "
             "Then isolate timestep, thermostat frequency and chain length separately. "
             "Optionally compare fixed-cell `calculation='md'` against PR `vc-md` at the same "
             "settings to separate the thermostat from the barostat coupling.")
elif overall in ("NOSE_TRANSIENT_THEN_CONVERGES", "NOSE_STABLE_WITH_NATIVE_INITIALIZATION"):
    L.append("Nose is not sufficient on its own here, so move the variable to the velocity path: "
             "(a) `ion_temperature='initial'` **with** an explicit `ATOMIC_VELOCITIES` card and no "
             "Nose; then (b) Nose **plus** explicit velocities. Inspect how the Nose chain is "
             "initialized when velocities are supplied (`vcinit` velocity branch vs `ranv`), "
             "since that is the remaining difference against job 50127909.")
elif overall in ("BOTH_RUNS_UNSTABLE", "CONTROL_UNSTABLE_NOSE_INTERPRETATION_BLOCKED"):
    L.append("The control is not a valid bounded reference in this run, so fix the common "
             "baseline first: re-establish the PR/electronic/cell settings that previously gave "
             "the bounded control (mean 170 K), then repeat this isolation.")
else:
    L.append("Evidence is partial; diagnose the incomplete run before drawing a thermostat "
             "conclusion.")
L.append("")
L.append("## Caveats")
L.append("")
L.append("- `ecutwfc=25 Ry` is retained deliberately from the reference deck; absolute pressures "
         "are dominated by Pulay stress and are **not** converged diamond pressures. Pressure and "
         "volume are used only to compare the two runs and to judge cell stability.")
L.append("- The reported `.tv`/stdout temperature is a **running average**; all instantaneous "
         "statistics above use `2*eka/((3N-3) k_B)` from `.eal`.")
L.append("- `tmp/*.vel` is not used: it is identically zero in this PR path.")
L.append("- `tmp/*.md` holds only the terminal state, so direct velocity reconstruction is "
         "available at the final frame only.")
L.append("- `sbatch` acceptance is not scientific validation; this report reflects parsed output.")
L.append("")
with open(os.path.join(ROOT, "FINAL_THERMOSTAT_ISOLATION_REPORT.md"), "w") as fh:
    fh.write("\n".join(L) + "\n")

print(json.dumps({"overall_status": overall, "control": ccls, "nose": ncls,
                  "control_mean_K": ct.get("mean"), "nose_mean_K": nt.get("mean"),
                  "nose_max_K": nt.get("max"), "nose_final200_K": nt.get("final200_mean"),
                  "strict_single_factor": init.get("strict_single_factor")}, indent=2))
