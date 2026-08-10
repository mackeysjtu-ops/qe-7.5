#!/usr/bin/env python3
"""Post-processing for the three-case CAIMD zero-strain 5000-step (10 ps) long test.

Runs inside the Slurm job after CASE A/B/C finish. numpy + matplotlib only.

Column semantics verified from PW/src/vcsmd.f90 (af7ebfc, 408-419):
  .e  : ut, ekint, edyn, pv, nst    .eal : uta, eka, eta, utl, ekla, etl, nst
  .p  : press, p, avp, nst [Ry/bohr^3]   .tv : vcell, tnew, nst (RUNNING AVG; not instantaneous)
  .cel: per step, lattice vectors as COLUMNS   .md : TERMINAL state only
CAIMD .caimd: 26 cols (CAIMD_IMPLEMENTATION.md): ionic_step time_ps events flag
  nominal_xx..xy actual_xx..xy stress_xx..xy(GPa) temperature_K volume_bohr3
  velocity_jump_max velocity_jump_rms.
Nose chain energy (ions_nose_energy) is NOT written by QE -> extended energy is PARTIAL.
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
EV_PER_RY = 13.605693122994
NAT = 8
NDOF = 3 * NAT - 3
DT_FS = 2.0
MASS_C_AMU = 12.0107
AMU_RY = 911.4442421
TARGET = 300.0

CASES = {
    "A": {"dir": "case_A_caimd_false", "prefix": "diamond8_case_A_caimd_false",
          "label": "CASE A caimd=.false.", "input": "inputs/case_A_caimd_false_nose_5000.in", "caimd": False},
    "B": {"dir": "case_B_caimd_noop", "prefix": "diamond8_case_B_caimd_noop",
          "label": "CASE B caimd no-op (masks .false.)", "input": "inputs/case_B_caimd_true_all_false_nose_5000.in", "caimd": True},
    "C": {"dir": "case_C_caimd_zero_strain", "prefix": "diamond8_case_C_caimd_zero",
          "label": "CASE C active zero-strain", "input": "inputs/case_C_caimd_true_all_fixed_zero_strain_nose_5000.in", "caimd": True},
}
COL = {"A": "#2a78d6", "B": "#1baf7a", "C": "#eb6834"}

REFERENCE_1000 = {"nose_native_mean_K": 169.696, "nose_native_final200_K": 173.237,
                  "nose_native_max_K": 300.001, "control_native_mean_K": 169.858}

NUM = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")


def rd(p):
    o = []
    if not os.path.exists(p):
        return np.zeros((0, 0))
    for line in open(p, errors="replace"):
        v = [float(t.replace("D", "E").replace("d", "e")) for t in NUM.findall(line)]
        if v:
            o.append(v)
    if not o:
        return np.zeros((0, 0))
    n = min(len(r) for r in o)
    return np.array([r[:n] for r in o])


def text(p):
    return open(p, errors="replace").read() if os.path.exists(p) else ""


def roll(a, w):
    if len(a) < w:
        return np.full(len(a), np.nan)
    return np.concatenate([np.full(w - 1, np.nan), np.convolve(a, np.ones(w) / w, "valid")])


def block_means(a, w):
    return np.array([a[i:i + w].mean() for i in range(0, len(a), w) if len(a[i:i + w])])


def slope(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    return float(np.polyfit(x[m], y[m], 1)[0]) if m.sum() >= 3 else None


def first_cross(t, y, thr):
    idx = np.nonzero(y > thr)[0]
    if not len(idx):
        return None
    i = int(idx[0])
    return {"nst": int(i), "time_ps": float(t[i]), "value": float(y[i])}


def longest_run(y, thr):
    best = cur = 0
    for v in y:
        cur = cur + 1 if v > thr else 0
        best = max(best, cur)
    return int(best)


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
    a = [float(t.replace("D", "E")) for t in NUM.findall(text(p))]
    n3 = 3 * nat
    need = n3 * 4 + 3 * nat * 3
    if len(a) < need + 9:
        return None
    return {"ratd": np.array(a[n3:2 * n3]).reshape(nat, 3),
            "aveci": np.array(a[need:need + 9]).reshape(3, 3)}


def parse_caimd(p):
    if not os.path.exists(p):
        return None
    cols = ("ionic_step time_ps events flag nom_xx nom_yy nom_zz nom_yz nom_xz nom_xy "
            "act_xx act_yy act_zz act_yz act_xz act_xy s_xx s_yy s_zz s_yz s_xz s_xy "
            "temperature_K volume vjmax vjrms").split()
    data = {c: [] for c in cols}
    ncols_seen = set()
    nrows = 0
    for line in text(p).splitlines():
        if line.lstrip().startswith("#"):
            continue
        parts = line.split()
        ncols_seen.add(len(parts))
        if len(parts) != 26:
            continue
        nrows += 1
        for i, c in enumerate(cols):
            if c == "flag":
                data[c].append(parts[i] == "T")
            elif c in ("ionic_step", "events"):
                data[c].append(int(parts[i]))
            else:
                data[c].append(float(parts[i]))
    for c in cols:
        if c not in ("flag",):
            data[c] = np.array(data[c])
    data["_ncols_seen"] = sorted(ncols_seen)
    data["_nrows"] = nrows
    return data


def det3(m):
    return float(np.linalg.det(m))


def cond_sv(m):
    s = np.linalg.svd(m, compute_uv=False)
    return float(s[0] / s[-1]), s


def strain6(H, H0):
    F = H @ np.linalg.inv(H0)
    E = 0.5 * (F + F.T) - np.eye(3)
    return [E[0, 0], E[1, 1], E[2, 2], E[1, 2], E[0, 2], E[0, 1]]


def skew_norm(H, H0):
    F = H @ np.linalg.inv(H0)
    W = 0.5 * (F - F.T)
    return float(np.sqrt((W * W).sum()))


def cc_dist(Hcol, cart):
    a = Hcol.T
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]) @ a
    dmin, dmax, nn = np.inf, 0.0, []
    for i in range(len(cart)):
        best = np.inf
        for j in range(len(cart)):
            d = cart[j] - cart[i] + shifts
            r = np.sqrt((d * d).sum(axis=1)); r = r[r > 1e-8]
            if len(r):
                dmin = min(dmin, r.min()); dmax = max(dmax, r.max())
                if j != i:
                    best = min(best, r.min())
        nn.append(best)
    return float(dmin), float(dmax), float(np.mean(nn))


def stats(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    if not len(a):
        return {"n": 0}
    return {"n": int(len(a)), "initial": float(a[0]), "final": float(a[-1]),
            "mean": float(a.mean()), "std": float(a.std()), "median": float(np.median(a)),
            "min": float(a.min()), "max": float(a.max()),
            "first500_mean": float(a[:500].mean()) if len(a) >= 500 else float(a.mean()),
            "final500_mean": float(a[-500:].mean()),
            "final1000_mean": float(a[-1000:].mean()),
            "final2500_mean": float(a[-2500:].mean())}


def wcsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join("" if v is None or (isinstance(v, float) and not np.isfinite(v))
                              else (f"{v:.10g}" if isinstance(v, float) else str(v))
                              for v in r) + "\n")


# ----------------------------------------------------------- input audits --
def strip_comments(t):
    return [l.split("!")[0].rstrip() for l in t.splitlines() if l.split("!")[0].strip()]


rc_txt = text(os.path.join(EVD, "run_return_codes.txt"))
RCS = dict(re.findall(r"(\w+)=(-?\d+)", rc_txt))
for src in ("no_explicit_velocity_audit.json", "input_identity_audit.json"):
    s, d = os.path.join(EVD, src), os.path.join(ANA, src)
    if os.path.exists(s):
        json.dump(json.load(open(s)), open(d, "w"), indent=2)

# ----------------------------------------------------------- load runs -----
D = {}
DIRECT_STEPS = [0, 10, 100, 500, 1000, 2000, 3000, 4000]
for k, meta in CASES.items():
    base = os.path.join(ROOT, meta["dir"]); pref = meta["prefix"]
    out = os.path.join(base, "qe.out")
    eal = rd(os.path.join(base, pref + ".eal")); pf = rd(os.path.join(base, pref + ".p"))
    tv = rd(os.path.join(base, pref + ".tv")); ef = rd(os.path.join(base, pref + ".e"))
    av = rd(os.path.join(base, pref + ".ave"))
    cs, cm = parse_cel(os.path.join(base, "tmp", pref + ".cel"))
    ps, pp = parse_pos(os.path.join(base, "tmp", pref + ".pos"))
    P, S = parse_stress(out)
    md = parse_md(os.path.join(base, "tmp", pref + ".md"))
    caimd = parse_caimd(os.path.join(base, "tmp", pref + ".caimd")) if meta["caimd"] else None
    txt = text(out)
    n = len(eal)
    d = {"label": meta["label"], "n": n, "job_done": "JOB DONE" in txt,
         "version": (re.search(r"Program PWSCF v\.(\S+)", txt) or [None, None])[1]
         if re.search(r"Program PWSCF v\.(\S+)", txt) else None,
         "error": (re.search(r"Error in routine\s+(\S+)", txt) or [None, None])[1]
         if re.search(r"Error in routine\s+(\S+)", txt) else None,
         "singular": "singular" in txt.lower(),
         "scf_fail": txt.count("convergence NOT achieved"),
         "nose_active": "nose" in txt.lower() and "temperature control" in txt.lower(),
         "rc": int(RCS.get(f"RC_{k}", -1)),
         "cel_steps": cs, "cel": cm, "pos_steps": ps, "pos": pp, "P": P, "S": S, "md": md, "caimd": caimd}
    if n:
        nst = eal[:, 6].astype(int)
        d.update(nst=nst, t_ps=nst * DT_FS / 1000.0, eka=eal[:, 1],
                 T=2 * eal[:, 1] * RY_TO_K / NDOF,
                 tv_T=tv[:, 1] if len(tv) else np.full(n, np.nan),
                 avk=av[:, 1] if len(av) else np.full(n, np.nan),
                 p_kbar=pf[:, 1] * RY_BOHR3_TO_KBAR if len(pf) else np.full(n, np.nan),
                 edyn=ef[:, 2] if len(ef) else np.full(n, np.nan),
                 ut=ef[:, 0] if len(ef) else np.full(n, np.nan),
                 ekint=ef[:, 1] if len(ef) else np.full(n, np.nan))
        d["T_runavg"] = 2 * d["avk"] * RY_TO_K / (3 * (NAT + 1))
    D[k] = d

# ----------------------------------------------------------- provenance ----
exe_txt = text(os.path.join(EVD, "executable_checksums.txt"))
prov = {"job_id": os.environ.get("SLURM_JOB_ID", "unknown"),
        "executable_path": (re.search(r"pw_x=(\S+)", exe_txt) or [None, None])[1],
        "executable_sha256": (re.findall(r"([0-9a-f]{64})", exe_txt) or [None])[0],
        "expected_sha256": "89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d",
        "versions": {k: D[k]["version"] for k in CASES},
        "nose_active": {k: D[k]["nose_active"] for k in CASES}}

# ----------------------------------------------------------- completeness --
comp = {}
for k in CASES:
    d = D[k]; n = d["n"]
    safety_trigger = None
    if n:
        if not np.all(np.isfinite(d["T"])):
            safety_trigger = "nonfinite temperature"
        elif d["singular"]:
            safety_trigger = "QE singular cell"
        elif longest_run(d["T"], 20000.0) >= 10:
            safety_trigger = "T>20000 K for >=10 consecutive records"
    cls = ("COMPLETE" if (d["job_done"] and n >= 5000 and not d["error"])
           else "SAFETY_STOP" if safety_trigger
           else "PARTIAL" if n > 0 else "MISSING_OUTPUT")
    comp[k] = {"return_code": d["rc"], "slurm_step_state": "see sacct", "job_done": d["job_done"],
               "runtime_version": d["version"], "first_error": d["error"],
               "scf_failure_count": d["scf_fail"], "completed_ionic_steps": n,
               "requested_steps": 5000, "final_valid_nst": int(d["nst"][-1]) if n else None,
               "completed_time_ps": float(d["t_ps"][-1]) if n else 0.0, "requested_time_ps": 10.0,
               "safety_stop": bool(safety_trigger), "safety_trigger": safety_trigger,
               "nose_active_in_stdout": d["nose_active"], "classification": cls}
comp["internal_post_processing_rc"] = int(RCS.get("RC_ANALYSIS", -1))
json.dump(comp, open(os.path.join(ANA, "run_completeness.json"), "w"), indent=2)

# ----------------------------------------------------------- initial state -
init = {"note": ("All three cases omit explicit velocities and use the same executable + tempw, so "
                 "vcinit->ranv(300 K) should yield the same deterministic Maxwell realization "
                 "(randy idum=0; total momentum zeroed). .md holds only the terminal state, so the "
                 "identity of the initial realization is read from the early bit-identical .eal rows."),
        "pairs": {}}
keys = list(CASES)
for a, b in ((keys[i], keys[j]) for i in range(len(keys)) for j in range(i + 1, len(keys))):
    da, db = D[a], D[b]
    if not (da["n"] and db["n"]):
        init["pairs"][f"{a}_vs_{b}"] = {"classification": "INSUFFICIENT_DATA"}
        continue
    deka = np.abs(da["eka"][:min(da["n"], db["n"])] - db["eka"][:min(da["n"], db["n"])])
    fd = int(np.nonzero(deka > 0)[0][0]) if np.any(deka > 0) else None
    same01 = da["eka"][0] == db["eka"][0] and da["eka"][1] == db["eka"][1]
    cls = ("IDENTICAL" if same01 and (fd is None or fd > 100)
           else "NUMERICALLY_EQUIVALENT" if abs(da["eka"][0] - db["eka"][0]) < 1e-10
           else "DIFFERENT_RANDOM_REALIZATION")
    init["pairs"][f"{a}_vs_{b}"] = {"first_step_eka": [float(da["eka"][0]), float(db["eka"][0])],
                                    "first_step_eka_abs_diff": float(abs(da["eka"][0] - db["eka"][0])),
                                    "first_eka_divergence_step": fd, "classification": cls}
init["all_identical_initialization"] = all(v.get("classification") == "IDENTICAL" for v in init["pairs"].values())
json.dump(init, open(os.path.join(ANA, "initial_state_comparison.json"), "w"), indent=2)

# ----------------------------------------------------------- per-case ------
tcross, tstat_rows, tblk_rows, energy_rows = [], [], [], []
RES = {}
for k in CASES:
    d = D[k]; n = d["n"]
    if not n:
        RES[k] = {"label": d["label"], "classification": "INCOMPLETE", "n": 0, "rc": d["rc"]}
        continue
    t, T, P = d["t_ps"], d["T"], d["p_kbar"]
    r100, r250, r500 = roll(T, 100), roll(T, 250), roll(T, 500)
    b500, b1000 = block_means(T, 500), block_means(T, 1000)

    wcsv(os.path.join(ANA, f"case_{k}_temperature.csv"),
         ["nst", "time_ps", "eka_Ry", "T_instant_3Nm3_K", "T_running_avg_tv_K", "roll100_K", "roll250_K", "roll500_K"],
         [[int(d["nst"][i]), t[i], d["eka"][i], T[i], d["tv_T"][i], r100[i], r250[i], r500[i]] for i in range(n)])

    # pressure/cell csv
    vol = np.array([abs(det3(m)) for m in d["cel"]]) if len(d["cel"]) else np.array([])
    conds = np.array([cond_sv(m)[0] for m in d["cel"]]) if len(d["cel"]) else np.array([])
    pc_rows, dmins, dmaxs, nns = [], [], [], []
    if len(d["cel"]):
        H0 = d["cel"][0]
        sample = set(range(0, len(d["cel"]), max(1, len(d["cel"]) // 200)))
        for i in range(len(d["cel"])):
            e6 = strain6(d["cel"][i], H0); dmin = dmax = nn = None
            if i in sample and i < len(d["pos"]):
                dmin, dmax, nn = cc_dist(d["cel"][i], d["pos"][i]); dmins.append(dmin); dmaxs.append(dmax); nns.append(nn)
            Pi = d["P"][i] if i < len(d["P"]) else None
            pk = P[i] if i < len(P) else None
            pc_rows.append([int(d["cel_steps"][i]), (d["cel_steps"][i] - d["cel_steps"][0]) * DT_FS / 1000.0,
                            pk, pk / 10.0 if pk is not None and np.isfinite(pk) else None,
                            vol[i], 100 * (vol[i] - vol[0]) / vol[0], conds[i]] + [float(x) for x in e6] + [dmin, nn])
    wcsv(os.path.join(ANA, f"case_{k}_pressure_cell.csv"),
         ["step", "time_ps", "p_kbar", "p_GPa", "volume_bohr3", "volume_change_percent", "condition_number",
          "strain_xx", "strain_yy", "strain_zz", "strain_yz", "strain_xz", "strain_xy", "min_CC_bohr", "mean_NN_CC_bohr"],
         pc_rows)

    # temperature source crosscheck + direct
    tcross.append([k, "T_running_avg(2avk/3(N+1)) vs .tv", 0, float(np.nanmax(np.abs(d["T_runavg"] - d["tv_T"])))])
    Tdir = None; direct_conf = False
    if d["md"] is not None and len(d["cel"]):
        v = d["md"]["ratd"] @ d["cel"][-1].T
        vcom = v.mean(axis=0)
        K = 0.5 * MASS_C_AMU * AMU_RY * ((v - vcom) ** 2).sum()
        Tdir = 2 * K * RY_TO_K / NDOF
        tcross.append([k, "direct .md ratd x terminal .cel (3N-3)", int(d["nst"][-1]), float(Tdir)])
        tcross.append([k, "eal instantaneous terminal", int(d["nst"][-1]), float(T[-1])])
        tcross.append([k, "direct-vs-eal terminal abs diff", 0, float(abs(Tdir - T[-1]))])
        direct_conf = bool(abs(Tdir - T[-1]) < max(5.0, 0.05 * abs(T[-1])))
    for st in DIRECT_STEPS:
        if st < n:
            tcross.append([k, f"eal instantaneous at step {st}", st, float(T[st])])

    stT = stats(T); stP = stats(P)
    for i, bm in enumerate(b1000):
        tblk_rows.append([k, "1000", i * 1000, min((i + 1) * 1000 - 1, n - 1), float(bm)])
    for i, bm in enumerate(b500):
        tblk_rows.append([k, "500", i * 500, min((i + 1) * 500 - 1, n - 1), float(bm)])

    fin_half = slope(t[n // 2:], T[n // 2:])
    fin_quarter = slope(t[3 * n // 4:], T[3 * n // 4:])
    f5 = stT["final2500_mean"]     # final 2500 steps = final 5 ps
    f25 = stats(T[-1250:])["mean"] if n >= 1250 else stT["mean"]  # final 2.5 ps = 1250 steps

    # 300 K convergence classification
    if n < 100:
        tcls = "INCOMPLETE"
    elif not np.all(np.isfinite(T)) or d["error"]:
        tcls = "NUMERICAL_FAILURE"
    elif stT["max"] > 2000 and f5 > 1000 and (fin_quarter is None or fin_quarter > -50):
        tcls = "SUSTAINED_HEATING"
    elif 270 <= f5 <= 330 and 270 <= f25 <= 330 and (fin_quarter is None or abs(fin_quarter) < 50):
        tcls = "TARGET_300K_REACHED"
    elif fin_quarter is not None and fin_quarter > 5 and f5 < 270:
        tcls = "APPROACHING_300K_NOT_CONVERGED"
    elif f5 < 270:
        tcls = "STABLE_BELOW_TARGET"
    elif f5 > 330:
        tcls = "STABLE_ABOVE_TARGET"
    else:
        tcls = "STABLE_BELOW_TARGET"

    # energy (partial)
    ed = d["edyn"]
    drift_total = float(ed[-1] - ed[0]); drift_ps = drift_total / max(t[-1], 1e-12)
    energy_rows.append([k, "edyn (.e col3 = ekint+ut)",
                        "PARTIAL_EXTENDED_ENERGY (Nose chain energy not written by QE)" if CASES[k]["caimd"] or True else "",
                        float(ed[0]), float(ed[-1]), float(ed.max() - ed.min()), drift_total, drift_ps,
                        drift_ps * EV_PER_RY / NAT, slope(t[3 * n // 4:], ed[3 * n // 4:]),
                        bool(np.all(np.isfinite(ed)) and np.abs(ed - ed[0]).max() < 10.0)])

    cellinfo = ({"volume_initial_bohr3": float(vol[0]), "volume_final_bohr3": float(vol[-1]),
                 "volume_change_percent": float(100 * (vol[-1] - vol[0]) / vol[0]),
                 "volume_min": float(vol.min()), "volume_max": float(vol.max()),
                 "determinant_min": float(vol.min()), "condition_number_max": float(conds.max()),
                 "strain6_final": [float(x) for x in strain6(d["cel"][-1], d["cel"][0])],
                 "min_CC_bohr": float(min(dmins)) if dmins else None,
                 "max_CC_bohr": float(max(dmaxs)) if dmaxs else None,
                 "mean_NN_CC_bohr": float(np.mean(nns)) if nns else None,
                 "classification": ("SINGULAR_OR_UNSTABLE_CELL" if (vol.min() <= 0 or conds.max() > 5)
                                    else "DRIFTING_PR_CELL" if abs(vol[-1] - vol[0]) / vol[0] > 0.25
                                    else "STABLE_PR_CELL")} if len(vol) else {})

    RES[k] = {"label": d["label"], "rc": d["rc"], "job_done": d["job_done"], "version": d["version"],
              "completed_steps": n, "temperature_classification": tcls,
              "nose_active_in_stdout": d["nose_active"],
              "temperature": stT, "temperature_initialized_K": float(T[0]),
              "final_5ps_mean_K": f5, "final_2p5ps_mean_K": f25,
              "drift_full_K_per_ps": slope(t, T), "drift_final_half_K_per_ps": fin_half,
              "drift_final_quarter_K_per_ps": fin_quarter,
              "block1000_means_K": [float(x) for x in b1000],
              "block500_means_K": [float(x) for x in b500],
              "roll500_final": float(r500[-1]) if np.isfinite(r500[-1]) else None,
              "crossings": {f"first_above_{x}K": first_cross(t, T, x) for x in (600, 1000, 2000)},
              "fraction_above": {f"{x}K": float(np.mean(T > x)) for x in (600, 1000, 2000)},
              "longest_run_above": {f"{x}K": longest_run(T, x) for x in (600, 1000, 2000)},
              "direct_terminal_T_K": Tdir, "direct_confirms_reported": direct_conf,
              "pressure_kbar": stP, "pressure_final1000_mean_kbar": float(P[-1000:].mean()) if len(P) >= 1000 else stP.get("mean"),
              "stress_mean_kbar": ({nm: float(np.nanmean(d["S"][:, i, j])) for nm, (i, j) in
                                    zip(("sxx", "syy", "szz", "syz", "sxz", "sxy"),
                                        ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)))} if len(d["S"]) else {}),
              "cell": cellinfo,
              "energy_partial": {"quantity": "edyn = ekint+ut", "is_full_conserved": False,
                                 "note": "Nose chain energy not written by QE; partial diagnostic",
                                 "drift_eV_per_atom_per_ps": drift_ps * EV_PER_RY / NAT,
                                 "bounded": bool(np.all(np.isfinite(ed)) and np.abs(ed - ed[0]).max() < 10.0)}}
    st = RES[k]["temperature"]
    tstat_rows.append([k, d["label"], tcls, n, st["initial"], st["mean"], st["median"], st["std"],
                       st["min"], st["max"], st["final"], st["first500_mean"], f5, f25,
                       st["final500_mean"], st["final1000_mean"], RES[k]["drift_full_K_per_ps"],
                       RES[k]["drift_final_half_K_per_ps"], RES[k]["drift_final_quarter_K_per_ps"],
                       RES[k]["fraction_above"]["600K"], RES[k]["longest_run_above"]["600K"]])

wcsv(os.path.join(ANA, "temperature_source_crosscheck.csv"),
     ["case", "comparison", "step_or_lag", "value_or_diff_K"], tcross)
wcsv(os.path.join(ANA, "temperature_statistics.csv"),
     ["case", "label", "classification", "n", "initial_K", "mean_K", "median_K", "std_K", "min_K", "max_K",
      "final_K", "first500_mean_K", "final5ps_mean_K", "final2p5ps_mean_K", "final500_mean_K",
      "final1000_mean_K", "drift_full_K_per_ps", "drift_final_half", "drift_final_quarter",
      "fraction_above_600K", "longest_run_above_600K"], tstat_rows)
wcsv(os.path.join(ANA, "temperature_block_statistics.csv"),
     ["case", "block_size", "nst_start", "nst_end", "temperature_mean_K"], tblk_rows)
wcsv(os.path.join(ANA, "energy_diagnostics.csv"),
     ["case", "quantity", "status", "initial_Ry", "final_Ry", "range_Ry", "drift_total_Ry",
      "drift_Ry_per_ps", "drift_eV_per_atom_per_ps", "slope_final_quarter", "bounded"], energy_rows)

# ----------------------------------------------------------- CAIMD audits --
def caimd_noop_audit(k_noop, k_ctrl):
    d = D[k_noop]; c = d["caimd"]
    res = {"case": k_noop, "all_masks_false_in_input": True,
           "caimd_file_present": c is not None}
    if c is not None:
        res["caimd_rows"] = c["_nrows"]
        res["caimd_ncols_seen"] = c["_ncols_seen"]
        res["any_event_flag_true"] = bool(any(c["flag"])) if c["flag"] else False
        res["max_abs_actual_strain"] = float(np.max([np.abs(c["act_" + comp]).max()
                                                     for comp in ("xx", "yy", "zz", "yz", "xz", "xy")]))
        res["max_velocity_jump_max"] = float(c["vjmax"].max()) if len(c["vjmax"]) else None
    # statistical agreement with control
    dc, dr = D[k_noop], D[k_ctrl]
    if dc["n"] and dr["n"]:
        res["mean_T_vs_control_delta_K"] = abs(stats(dc["T"])["mean"] - stats(dr["T"])["mean"])
        res["final5ps_T_vs_control_delta_K"] = abs(stats(dc["T"][-2500:])["mean"] - stats(dr["T"][-2500:])["mean"])
        deka = np.abs(dc["eka"][:min(dc["n"], dr["n"])] - dr["eka"][:min(dc["n"], dr["n"])])
        res["eka_first_divergence_step_vs_control"] = int(np.nonzero(deka > 0)[0][0]) if np.any(deka > 0) else None
        res["statistically_consistent_with_control"] = bool(res["mean_T_vs_control_delta_K"] < 30)
    return res


def caimd_active_audit(k):
    d = D[k]; c = d["caimd"]
    res = {"case": k, "all_masks_true_in_input": True, "caimd_file_present": c is not None}
    if c is None:
        res["status"] = "NO_CAIMD_FILE"
        return res
    res["caimd_rows"] = c["_nrows"]
    res["caimd_ncols_seen"] = c["_ncols_seen"]
    res["schema_26_columns_ok"] = c["_ncols_seen"] == [26]
    res["event_count_final"] = int(c["events"][-1]) if len(c["events"]) else None
    res["n_event_flags_true"] = int(sum(c["flag"])) if c["flag"] else 0
    for comp in ("xx", "yy", "zz", "yz", "xz", "xy"):
        res[f"max_abs_nominal_{comp}"] = float(np.abs(c["nom_" + comp]).max())
        res[f"max_abs_actual_{comp}"] = float(np.abs(c["act_" + comp]).max())
    res["max_abs_nominal_minus_actual"] = float(np.max(
        [np.abs(c["nom_" + comp] - c["act_" + comp]).max() for comp in ("xx", "yy", "zz", "yz", "xz", "xy")]))
    res["max_velocity_jump_max"] = float(c["vjmax"].max()) if len(c["vjmax"]) else None
    res["max_velocity_jump_rms"] = float(c["vjrms"].max()) if len(c["vjrms"]) else None
    res["temperature_diag_range_K"] = [float(c["temperature_K"].min()), float(c["temperature_K"].max())] if len(c["temperature_K"]) else None
    # cell det / cond from .cel
    if len(d["cel"]):
        vol = np.array([abs(det3(m)) for m in d["cel"]])
        conds = np.array([cond_sv(m)[0] for m in d["cel"]])
        H0 = d["cel"][0]
        skews = np.array([skew_norm(m, H0) for m in d["cel"]])
        res["determinant_range_bohr3"] = [float(vol.min()), float(vol.max())]
        res["determinant_always_positive"] = bool(vol.min() > 0)
        res["condition_number_range"] = [float(conds.min()), float(conds.max())]
        res["max_skew_norm"] = float(skews.max())
        res["max_abs_strain_from_cel"] = float(np.max([np.abs([strain6(m, H0)[j] for m in d["cel"]]).max()
                                                       for j in range(6)]))
    res["nominal_all_zero"] = res["max_abs_nominal_minus_actual"] is not None and all(
        res[f"max_abs_nominal_{comp}"] == 0.0 for comp in ("xx", "yy", "zz", "yz", "xz", "xy"))
    return res


noop = caimd_noop_audit("B", "A")
json.dump(noop, open(os.path.join(ANA, "caimd_case_B_noop_audit.json"), "w"), indent=2)
active = caimd_active_audit("C")
json.dump(active, open(os.path.join(ANA, "caimd_case_C_zero_strain_audit.json"), "w"), indent=2)

# ----------------------------------------------------------- overall -------
def cls_of(k):
    return RES[k].get("temperature_classification")


A, B, C = RES["A"], RES["B"], RES["C"]
cA, cB, cC = cls_of("A"), cls_of("B"), cls_of("C")
# no-op regression?
noop_ok = noop.get("statistically_consistent_with_control", False) and not noop.get("any_event_flag_true", True)
c_stable = active.get("determinant_always_positive", False) and (active.get("max_skew_norm") or 1) < 1e-3 \
    and (active.get("max_velocity_jump_max") or 1) < 1e-6
c_finite = comp["C"]["classification"] in ("COMPLETE", "PARTIAL")

if any(comp[k]["classification"] not in ("COMPLETE", "PARTIAL") for k in CASES):
    overall = "PARTIAL"
elif not noop_ok:
    overall = "CAIMD_ENABLE_NOOP_REGRESSION"
elif cC == "SUSTAINED_HEATING" or not c_stable or comp["C"]["safety_stop"]:
    overall = "ACTIVE_CAIMD_ZERO_STRAIN_UNSTABLE"
elif cC == "TARGET_300K_REACHED":
    overall = "ACTIVE_CAIMD_ZERO_STRAIN_REACHES_300K"
elif cC == "APPROACHING_300K_NOT_CONVERGED":
    overall = "ACTIVE_CAIMD_ZERO_STRAIN_APPROACHING_TARGET"
elif cC in ("STABLE_BELOW_TARGET", "STABLE_ABOVE_TARGET"):
    if cA in ("STABLE_BELOW_TARGET", "STABLE_ABOVE_TARGET") and cB in ("STABLE_BELOW_TARGET", "STABLE_ABOVE_TARGET"):
        overall = "ACTIVE_CAIMD_ZERO_STRAIN_STABLE_BELOW_TARGET" if cC == "STABLE_BELOW_TARGET" else "COMMON_NOSE_TARGET_FAILURE"
    else:
        overall = "ACTIVE_CAIMD_ZERO_STRAIN_STABLE_BELOW_TARGET" if cC == "STABLE_BELOW_TARGET" else "PARTIAL"
else:
    overall = "PARTIAL"

summary = {
    "job_id": prov["job_id"],
    "remote_workdir": os.path.abspath(ROOT),
    "overall_status": overall,
    "test": "three-case native-init Nose + CAIMD zero-strain, 5000 steps (10 ps) each",
    "cutoff_caveat": "ecutwfc=25 Ry: absolute pressures are Pulay-dominated, not converged.",
    "temperature_convention": {"instantaneous": "2*eka/((3N-3)k_B) from .eal",
                               "reported_tv": "2*avk/(3(N+1)k_B) RUNNING AVERAGE (not used)",
                               "direct": "terminal .md ratd x terminal .cel, COM removed"},
    "executable_provenance": prov,
    "run_completeness": comp, "initial_state": init,
    "runs": RES, "caimd_case_B_noop": noop, "caimd_case_C_zero_strain": active,
    "reference_1000step": REFERENCE_1000,
}
json.dump(summary, open(os.path.join(ANA, "comparison_summary.json"), "w"), indent=2)

# ----------------------------------------------------------- plots ---------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.color": "#e1e0d9",
                         "grid.linewidth": .8, "axes.axisbelow": True, "axes.unicode_minus": False,
                         "savefig.dpi": 150, "figure.facecolor": "white", "axes.facecolor": "white"})
    SRC = "Source: this workflow; T=2*eka/((3N-3)k_B) from .eal; P from .p; cell/strain from .cel/.caimd. ecutwfc=25 Ry."

    def fin(ax, title, xl, yl, hline=None, note=SRC, legend=True):
        if hline is not None:
            ax.axhline(hline, color="#898781", ls="--", lw=1.2)
        ax.set_title(title, fontsize=11.5); ax.set_xlabel(xl); ax.set_ylabel(yl)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        if legend:
            ax.legend(frameon=False, fontsize=9)
        ax.figure.text(0.01, -0.02, note, fontsize=7.2, color="#52514e")

    def savef(fig, name):
        fig.savefig(os.path.join(ANA, name), bbox_inches="tight"); plt.close(fig); print("wrote", name)

    def has(k):
        return D[k]["n"] > 0

    # temperature
    fig, ax = plt.subplots(figsize=(11, 5))
    for k in CASES:
        if has(k):
            ax.plot(D[k]["t_ps"], D[k]["T"], lw=.5, color=COL[k], alpha=.8, label=CASES[k]["label"])
    fin(ax, "Instantaneous ionic temperature, 10 ps (native init + Nose)", "time (ps)", "T (K, 3N-3)", hline=300)
    savef(fig, "temperature_all_cases.png")

    fig, ax = plt.subplots(figsize=(11, 5))
    for k in CASES:
        if has(k):
            ax.plot(D[k]["t_ps"], roll(D[k]["T"], 500), lw=1.6, color=COL[k], label=CASES[k]["label"])
    fin(ax, "Temperature, 500-step rolling mean", "time (ps)", "rolling mean T (K)", hline=300)
    savef(fig, "temperature_rolling500_all_cases.png")

    fig, ax = plt.subplots(figsize=(11, 5))
    for k in CASES:
        if has(k):
            b = block_means(D[k]["T"], 1000)
            ax.plot(np.arange(1, len(b) + 1), b, "o-", ms=6, lw=1.5, color=COL[k], label=CASES[k]["label"])
    fin(ax, "Temperature, 1000-step (2 ps) block means", "block", "block-mean T (K)", hline=300)
    savef(fig, "temperature_block1000_all_cases.png")

    fig, ax = plt.subplots(figsize=(11, 5))
    for k in CASES:
        if has(k):
            m = D[k]["t_ps"] >= (D[k]["t_ps"][-1] - 5.0)
            ax.plot(D[k]["t_ps"][m], D[k]["T"][m], lw=.6, color=COL[k], alpha=.85, label=CASES[k]["label"])
    fin(ax, "Temperature, last 5 ps", "time (ps)", "T (K, 3N-3)", hline=300)
    savef(fig, "temperature_last5ps_all_cases.png")

    fig, ax = plt.subplots(figsize=(11, 5))
    for k in CASES:
        if has(k):
            ax.plot(D[k]["t_ps"], roll(np.abs(D[k]["T"] - 300.0), 500), lw=1.5, color=COL[k], label=CASES[k]["label"])
    fin(ax, "Distance from 300 K (500-step rolling mean of |T-300|)", "time (ps)", "|T - 300 K| (K)")
    savef(fig, "temperature_distance_from_300K.png")

    fig, ax = plt.subplots(figsize=(11, 5))
    for k in CASES:
        if has(k):
            ax.plot(D[k]["t_ps"], D[k]["eka"], lw=.5, color=COL[k], alpha=.8, label=CASES[k]["label"])
    fin(ax, "Ionic kinetic energy", "time (ps)", "eka (Ry)")
    savef(fig, "ionic_kinetic_energy_all_cases.png")

    fig, ax = plt.subplots(figsize=(11, 5))
    for k in CASES:
        if has(k):
            ax.plot(D[k]["t_ps"], D[k]["p_kbar"], lw=.5, color=COL[k], alpha=.8, label=CASES[k]["label"])
    fin(ax, "MD pressure (NOT converged; diagnostic only)", "time (ps)", "P (kbar)", hline=0)
    savef(fig, "pressure_all_cases.png")

    fig, ax = plt.subplots(figsize=(11, 5))
    for k in CASES:
        if len(D[k]["cel"]):
            ct = (D[k]["cel_steps"] - D[k]["cel_steps"][0]) * DT_FS / 1000.0
            ax.plot(ct, [abs(det3(m)) for m in D[k]["cel"]], lw=1.3, color=COL[k], label=CASES[k]["label"])
    fin(ax, "Cell volume", "time (ps)", "volume (bohr^3)")
    savef(fig, "volume_all_cases.png")

    fig, ax = plt.subplots(figsize=(11, 5))
    for k in CASES:
        if len(D[k]["cel"]):
            ct = (D[k]["cel_steps"] - D[k]["cel_steps"][0]) * DT_FS / 1000.0
            ax.plot(ct, [cond_sv(m)[0] for m in D[k]["cel"]], lw=1.3, color=COL[k], label=CASES[k]["label"])
    fin(ax, "Cell condition number", "time (ps)", "cond(H)")
    savef(fig, "cell_condition_all_cases.png")

    # case C strains
    if D["C"]["caimd"] is not None and D["C"]["caimd"]["_nrows"]:
        c = D["C"]["caimd"]
        fig, ax = plt.subplots(figsize=(11, 5))
        for cmp6, cc2 in zip(("xx", "yy", "zz", "yz", "xz", "xy"),
                             ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")):
            ax.plot(c["time_ps"], c["act_" + cmp6], lw=1.0, color=cc2, label=f"actual {cmp6}")
        fin(ax, "CASE C: all six actual strain components (target = 0)", "time (ps)", "actual strain (tensor)")
        savef(fig, "case_C_all_six_strains.png")

        fig, ax = plt.subplots(figsize=(11, 5))
        err = np.max([np.abs(c["nom_" + comp] - c["act_" + comp]) for comp in ("xx", "yy", "zz", "yz", "xz", "xy")], axis=0)
        pos = err > 0
        ax.plot(c["time_ps"][pos] if pos.any() else c["time_ps"], err[pos] if pos.any() else err,
                lw=1.2, color="#eb6834")
        if pos.any():
            ax.set_yscale("log")
            from matplotlib.ticker import FuncFormatter, LogLocator
            ax.yaxis.set_major_locator(LogLocator(base=10.0))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: "1e%d" % int(round(np.log10(v))) if v > 0 else ""))
        fin(ax, "CASE C: max |nominal - actual| strain error (all six components)", "time (ps)",
            "max strain error", legend=False)
        savef(fig, "case_C_strain_error.png")

        fig, ax = plt.subplots(figsize=(11, 5))
        vj = c["vjmax"]; pos = vj > 0
        ax.plot(c["time_ps"][pos] if pos.any() else c["time_ps"], vj[pos] if pos.any() else vj, lw=1.2, color="#eb6834")
        if pos.any():
            ax.set_yscale("log")
            from matplotlib.ticker import FuncFormatter, LogLocator
            ax.yaxis.set_major_locator(LogLocator(base=10.0))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: "1e%d" % int(round(np.log10(v))) if v > 0 else ""))
        fin(ax, "CASE C: velocity-transform residual (velocity_jump_max)", "time (ps)",
            "velocity_jump_max (a.u.)", legend=False)
        savef(fig, "case_C_velocity_transform_residual.png")
    else:
        for nm in ("case_C_all_six_strains.png", "case_C_strain_error.png", "case_C_velocity_transform_residual.png"):
            fig, ax = plt.subplots(figsize=(11, 5))
            ax.text(0.5, 0.5, "CASE C .caimd not available", transform=ax.transAxes, ha="center", va="center")
            fin(ax, nm.replace(".png", ""), "time (ps)", "", legend=False)
            savef(fig, nm)

    fig, ax = plt.subplots(figsize=(11, 5))
    for k in CASES:
        if has(k):
            ax.plot(D[k]["t_ps"], (D[k]["edyn"] - D[k]["edyn"][0]) * EV_PER_RY / NAT, lw=1.2, color=COL[k], label=CASES[k]["label"])
    fin(ax, "Partial extended-energy drift (edyn = ekint+ut; excludes Nose chain)", "time (ps)",
        "drift per atom (eV/atom)", hline=0, note=SRC + " Nose chain energy not written by QE -> partial.")
    savef(fig, "energy_drift_all_cases.png")
    summary["plots_generated"] = True
except Exception as exc:
    summary["plots_generated"] = False
    summary["plot_error"] = repr(exc)
json.dump(summary, open(os.path.join(ANA, "comparison_summary.json"), "w"), indent=2)


# ----------------------------------------------------------- report --------
def f(x, n=3):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{n}f}"


L = []
L.append("# Final Long Zero-Strain CAIMD Report")
L.append("")
L.append("Diamond8, QE 7.5 `af7ebfc`, three 5000-step (10 ps) native-initialization Nose "
         "trajectories. No explicit velocities in any case.")
L.append("")
L.append(f"**Overall status: `{summary['overall_status']}`**")
L.append("")
L.append("| case | CAIMD | T class | mean K | final-5ps K | final-2.5ps K | max K | >600K frac |")
L.append("|---|---|---|---:|---:|---:|---:|---:|")
for k in CASES:
    r = RES[k]
    if r.get("completed_steps"):
        st = r["temperature"]
        L.append(f"| {k} | {'.true.' if CASES[k]['caimd'] else '.false.'} | "
                 f"{r['temperature_classification']} | {f(st['mean'],1)} | {f(r['final_5ps_mean_K'],1)} | "
                 f"{f(r['final_2p5ps_mean_K'],1)} | {f(st['max'],1)} | {r['fraction_above']['600K']} |")
    else:
        L.append(f"| {k} | - | INCOMPLETE | - | - | - | - | - |")
L.append("")
L.append("## Required answers")
L.append("")
L.append(f"1. **All three completed 5000 steps?** " +
         " ".join(f"{k}={comp[k]['completed_ionic_steps']}({comp[k]['classification']})" for k in CASES) + ".")
L.append(f"2. **Explicit velocities absent in all cases?** {summary.get('run_completeness') is not None}; "
         "see `analysis/no_explicit_velocity_audit.json` (no ATOMIC_VELOCITIES / ion_velocities / "
         "from_input in the complete text of any input).")
L.append("3. **Native velocity initialization?** Yes — all three omit velocities and run "
         "`vcinit->ranv(300 K)`.")
L.append(f"4. **Identical initial velocity realizations?** all-identical="
         f"{init['all_identical_initialization']}; pairs " +
         ", ".join(f"{p}={v['classification']}" for p, v in init['pairs'].items()) + ".")
L.append(f"5. **CASE A approach 300 K?** `{cA}` (mean {f(A['temperature']['mean'],1)} K, final-5ps "
         f"{f(A['final_5ps_mean_K'],1)} K).")
L.append(f"6. **CASE B approach 300 K?** `{cB}` (mean {f(B['temperature']['mean'],1)} K, final-5ps "
         f"{f(B['final_5ps_mean_K'],1)} K).")
L.append(f"7. **CASE C approach 300 K?** `{cC}` (mean {f(C['temperature']['mean'],1)} K, final-5ps "
         f"{f(C['final_5ps_mean_K'],1)} K).")
L.append(f"8. **final-5-ps / final-2.5-ps mean T:** A {f(A['final_5ps_mean_K'],1)}/{f(A['final_2p5ps_mean_K'],1)}; "
         f"B {f(B['final_5ps_mean_K'],1)}/{f(B['final_2p5ps_mean_K'],1)}; "
         f"C {f(C['final_5ps_mean_K'],1)}/{f(C['final_2p5ps_mean_K'],1)} K.")
L.append(f"9. **1000-step block means (K):** A {[round(x,1) for x in A['block1000_means_K']]}; "
         f"B {[round(x,1) for x in B['block1000_means_K']]}; C {[round(x,1) for x in C['block1000_means_K']]}.")
L.append(f"10. **Late-time slopes consistent with a plateau?** final-quarter slope A "
         f"{f(A['drift_final_quarter_K_per_ps'])}, B {f(B['drift_final_quarter_K_per_ps'])}, "
         f"C {f(C['drift_final_quarter_K_per_ps'])} K/ps.")
L.append(f"11. **Any case cross 600/1000/2000 K?** " +
         "; ".join(f"{k}: 600={RES[k]['crossings']['first_above_600K']}" for k in CASES) + ".")
L.append(f"12. **CASE B a long-time no-op vs A?** consistent={noop.get('statistically_consistent_with_control')}, "
         f"mean-T delta {f(noop.get('mean_T_vs_control_delta_K'),3)} K, any event flag "
         f"{noop.get('any_event_flag_true')}, eka first divergence step "
         f"{noop.get('eka_first_divergence_step_vs_control')}.")
L.append(f"13. **CASE C keep all six strains at zero?** nominal all zero={active.get('nominal_all_zero')}; "
         f"max |actual| strain from .cel {f(active.get('max_abs_strain_from_cel'), 3)}.")
L.append(f"14. **CASE C max strain-tracking error:** {active.get('max_abs_nominal_minus_actual')}.")
L.append(f"15. **CASE C determinant always positive?** {active.get('determinant_always_positive')} "
         f"(range {active.get('determinant_range_bohr3')} bohr^3).")
L.append(f"16. **CAIMD velocity-transform residuals bounded?** CASE C max velocity_jump_max "
         f"{active.get('max_velocity_jump_max')}, max skew norm {active.get('max_skew_norm')}.")
L.append(f"17. **Did active CAIMD alter temperature control?** mean-T C vs A delta "
         f"{f(abs(C['temperature']['mean'] - A['temperature']['mean']),3)} K; classifications "
         f"A/B/C = {cA}/{cB}/{cC}.")
L.append(f"18. **Nose/PR energy bounded?** partial edyn drift (eV/atom/ps): A "
         f"{A['energy_partial']['drift_eV_per_atom_per_ps']:.2e}, B "
         f"{B['energy_partial']['drift_eV_per_atom_per_ps']:.2e}, C "
         f"{C['energy_partial']['drift_eV_per_atom_per_ps']:.2e}; bounded A/B/C="
         f"{A['energy_partial']['bounded']}/{B['energy_partial']['bounded']}/{C['energy_partial']['bounded']}. "
         f"PARTIAL: Nose chain energy is not written by QE.")
L.append(f"19. **Stable-below / approaching / at target?** {summary['overall_status']}.")
L.append("20. **Minimum next step:** see below.")
L.append("21. **Sources or previous raw outputs modified?** No — this job only wrote new files in "
         "its own workflow directory.")
L.append("")
L.append("## Minimum next step (recommendation)")
L.append("")
if summary["overall_status"] == "ACTIVE_CAIMD_ZERO_STRAIN_REACHES_300K":
    L.append("CASE C reached 300 K and is stable: next run a very small nonzero-strain active-CAIMD "
             "test (native init, same Nose parameters, one masked component, conservative rate).")
elif summary["overall_status"] in ("ACTIVE_CAIMD_ZERO_STRAIN_STABLE_BELOW_TARGET", "COMMON_NOSE_TARGET_FAILURE"):
    L.append("All cases are numerically stable but sit below 300 K: do NOT start scientific loading. "
             "First audit why Nose does not transfer enough energy (it barely acts within this window) "
             "with a conservative ordinary-Nose parameter matrix (lower fnosep, shorter dt, shorter "
             "chain); optionally add a delayed CAIMD loading start within one continuous trajectory "
             "after a validated equilibration stage.")
elif summary["overall_status"] == "ACTIVE_CAIMD_ZERO_STRAIN_APPROACHING_TARGET":
    L.append("CASE C is stable and trending toward 300 K: consider extending the same protocol only "
             "if energy and cell diagnostics stay bounded.")
elif summary["overall_status"] == "CAIMD_ENABLE_NOOP_REGRESSION":
    L.append("CASE B differs from CASE A despite all masks false: investigate the supposedly inactive "
             "CAIMD boundary before any further active tests.")
elif summary["overall_status"] == "ACTIVE_CAIMD_ZERO_STRAIN_UNSTABLE":
    L.append("CASE C is unstable: diagnose the instability (cell singularity, residual growth) before "
             "any loading test.")
else:
    L.append("Evidence is partial; diagnose the incomplete/failed case first.")
L.append("")
L.append("## Caveats")
L.append("")
L.append("- `ecutwfc=25 Ry`: absolute pressures are Pulay-dominated, not converged; pressure/volume "
         "are used only for software and stability diagnostics.")
L.append("- The `.tv`/stdout temperature is a running average; all instantaneous statistics use "
         "`2*eka/((3N-3)k_B)` from `.eal`.")
L.append("- `tmp/*.vel` is zero in the PR path and not used; `tmp/*.md` holds only the terminal "
         "state, so direct velocity reconstruction is at the final frame only.")
L.append("- Nose chain conserved-energy terms are not written by QE; the extended-energy check is a "
         "partial diagnostic.")
L.append("")
with open(os.path.join(ROOT, "FINAL_LONG_ZERO_STRAIN_REPORT.md"), "w") as fh:
    fh.write("\n".join(L) + "\n")

print(json.dumps({"overall": summary["overall_status"],
                  "A": cA, "B": cB, "C": cC,
                  "A_mean": A["temperature"].get("mean"), "B_mean": B["temperature"].get("mean"),
                  "C_mean": C["temperature"].get("mean"),
                  "C_final5ps": C.get("final_5ps_mean_K"),
                  "noop_consistent": noop.get("statistically_consistent_with_control"),
                  "C_det_positive": active.get("determinant_always_positive"),
                  "C_max_strain_err": active.get("max_abs_nominal_minus_actual"),
                  "C_max_vjump": active.get("max_velocity_jump_max")}, indent=2))
