#!/usr/bin/env python3
"""Shared analysis + pass/fail gate for every stage of the diamond8 vc-relax/NVT/NPT/CAIMD chain.

CLI:
  analyze_stage.py --stage {vcrelax|nvt|npt|caimd_zero|caimd_nonzero} --dir <casedir>
                   --prefix <prefix> --out <outdir> [--ref-volume <bohr3>] [--rc <qe_return_code>]

Writes stage CSV/JSON/report and exits 0 iff the stage pass gate is satisfied, nonzero otherwise.
Standard library + numpy only (no network, no seaborn, matplotlib optional and not required here).

Column semantics verified from PW/src/vcsmd.f90 (af7ebfc) and dynamics_module:
  .eal : uta, eka, eta, utl, ekla, etl, nst   (eka = ionic kinetic energy, Ry)
  .e   : ut, ekint, edyn, pv, nst
  .p   : press, p, avp, nst   (p in Ry/bohr^3)
  .tv  : vcell, tnew, nst     (tnew is a RUNNING AVERAGE -> not used as instantaneous)
  .cel : per step, lattice vectors as COLUMNS (bohr)
  .caimd : 26 cols: ionic_step time_ps events flag nom(6) act(6) stress(6,GPa) T volume vjmax vjrms
Instantaneous ionic temperature: T = 2*eka/((3N-3)*k_B), 3N-3 = 21 (= ndega=get_ndof()).
Fixed-cell md temperature also uses 3N-3 (get_ndof convention).
Nose chain energy is NOT written by QE -> extended energy is PARTIAL.
"""
import argparse
import json
import math
import os
import re
import sys

import numpy as np

RY_TO_K = 157887.51240116
RY_BOHR3_TO_KBAR = 147105.078
BOHR_TO_ANG = 0.529177210903
EV_PER_RY = 13.605693122994
NAT = 8
NDOF = 3 * NAT - 3
MASS_C_AMU = 12.0107
AMU_RY = 911.4442421
NUM = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")


def text(p):
    return open(p, errors="replace").read() if os.path.exists(p) else ""


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


def roll(a, w):
    if len(a) < w:
        return np.full(len(a), np.nan)
    return np.concatenate([np.full(w - 1, np.nan), np.convolve(a, np.ones(w) / w, "valid")])


def blocks(a, w):
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
    tk = [l.split() for l in text(p).splitlines() if l.strip()]
    i = 0
    while i + 3 <= len(tk):
        if len(tk[i]) == 2:
            try:
                s = int(float(tk[i][0])); m = [[float(x) for x in tk[i + r]] for r in (1, 2, 3)]
            except (ValueError, IndexError):
                i += 1; continue
            steps.append(s); mats.append(m); i += 4
        else:
            i += 1
    return np.array(steps), np.array(mats)


def parse_pos(p, nat=NAT):
    steps, pos = [], []
    tk = [l.split() for l in text(p).splitlines() if l.strip()]
    i = 0
    while i < len(tk):
        if len(tk[i]) == 2:
            try:
                s = int(float(tk[i][0])); blk = [[float(x) for x in tk[i + 1 + k]] for k in range(nat)]
            except (ValueError, IndexError):
                i += 1; continue
            steps.append(s); pos.append(blk); i += 1 + nat
        else:
            i += 1
    return np.array(steps), np.array(pos)


def parse_stress_stdout(p):
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
    n3 = 3 * nat; need = n3 * 4 + 3 * nat * 3
    if len(a) < need + 9:
        return None
    return {"ratd": np.array(a[n3:2 * n3]).reshape(nat, 3), "aveci": np.array(a[need:need + 9]).reshape(3, 3)}


def parse_caimd(p):
    if not os.path.exists(p):
        return None
    cols = ("step time events flag nom_xx nom_yy nom_zz nom_yz nom_xz nom_xy "
            "act_xx act_yy act_zz act_yz act_xz act_xy s_xx s_yy s_zz s_yz s_xz s_xy "
            "temperature volume vjmax vjrms").split()
    data = {c: [] for c in cols}
    ncols, nrows = set(), 0
    for line in text(p).splitlines():
        if line.lstrip().startswith("#"):
            continue
        parts = line.split()
        ncols.add(len(parts))
        if len(parts) != 26:
            continue
        nrows += 1
        for i, c in enumerate(cols):
            if c == "flag":
                data[c].append(parts[i] == "T")
            elif c in ("step", "events"):
                data[c].append(int(parts[i]))
            else:
                data[c].append(float(parts[i].replace("D", "E")))
    for c in cols:
        if c != "flag":
            data[c] = np.array(data[c])
    data["_ncols"] = sorted(ncols); data["_nrows"] = nrows
    return data


def det3(m):
    return float(np.linalg.det(m))


def cond2(m):
    s = np.linalg.svd(np.array(m), compute_uv=False)
    return float(s[0] / s[-1])


def strain6(H, H0):
    F = np.array(H) @ np.linalg.inv(np.array(H0))
    E = 0.5 * (F + F.T) - np.eye(3)
    return [E[0, 0], E[1, 1], E[2, 2], E[1, 2], E[0, 2], E[0, 1]]


def skew_norm(H, H0):
    F = np.array(H) @ np.linalg.inv(np.array(H0))
    W = 0.5 * (F - F.T)
    return float(np.sqrt((W * W).sum()))


def cc_dist(Hcol, cart):
    a = np.array(Hcol).T
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]) @ a
    dmin, dmax, nn = np.inf, 0.0, []
    for i in range(len(cart)):
        best = np.inf
        for j in range(len(cart)):
            d = np.array(cart[j]) - np.array(cart[i]) + shifts
            r = np.sqrt((d * d).sum(axis=1)); r = r[r > 1e-8]
            if len(r):
                dmin = min(dmin, r.min()); dmax = max(dmax, r.max())
                if j != i:
                    best = min(best, r.min())
        nn.append(best)
    return float(dmin), float(dmax), float(np.mean(nn))


def tstats(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    if not len(a):
        return {"n": 0}
    return {"n": int(len(a)), "mean": float(a.mean()), "std": float(a.std()), "min": float(a.min()),
            "max": float(a.max()), "final": float(a[-1]), "median": float(np.median(a))}


def block_unc(a, w):
    half = a[len(a) // 2:]
    bm = blocks(half, w)
    if len(bm) < 2:
        return {"n_blocks": int(len(bm)), "block_mean": float(bm.mean()) if len(bm) else None, "sem": None}
    return {"n_blocks": int(len(bm)), "block_mean": float(bm.mean()),
            "block_std": float(bm.std(ddof=1)), "sem": float(bm.std(ddof=1) / math.sqrt(len(bm)))}


def wcsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join("" if v is None or (isinstance(v, float) and not np.isfinite(v))
                              else (f"{v:.10g}" if isinstance(v, float) else str(v)) for v in r) + "\n")


def load_md(casedir, prefix, dt_fs):
    base = casedir
    out = os.path.join(base, "qe.out")
    eal = rd(os.path.join(base, prefix + ".eal"))
    pf = rd(os.path.join(base, prefix + ".p"))
    ef = rd(os.path.join(base, prefix + ".e"))
    cs, cm = parse_cel(os.path.join(base, "tmp", prefix + ".cel"))
    ps, pp = parse_pos(os.path.join(base, "tmp", prefix + ".pos"))
    P, S = parse_stress_stdout(out)
    md = parse_md(os.path.join(base, "tmp", prefix + ".md"))
    caimd = parse_caimd(os.path.join(base, "tmp", prefix + ".caimd"))
    txt = text(out)
    d = {"txt": txt, "job_done": "JOB DONE" in txt,
         "error": (re.search(r"Error in routine\s+(\S+)", txt) or [None, None])[1]
         if re.search(r"Error in routine\s+(\S+)", txt) else None,
         "version": (re.search(r"Program PWSCF v\.(\S+)", txt) or [None, None])[1]
         if re.search(r"Program PWSCF v\.(\S+)", txt) else None,
         "scf_fail": txt.count("convergence NOT achieved"),
         "nose_ion": "nose" in txt.lower() and "temperature control" in txt.lower(),
         "P_stdout": P, "S": S, "cel_steps": cs, "cel": cm, "pos_steps": ps, "pos": pp,
         "md": md, "caimd": caimd}
    n = len(eal)
    d["n"] = n
    if n:
        nst = eal[:, 6].astype(int)
        d.update(nst=nst, t_ps=nst * dt_fs / 1000.0, eka=eal[:, 1], T=2 * eal[:, 1] * RY_TO_K / NDOF,
                 p_kbar=pf[:, 1] * RY_BOHR3_TO_KBAR if len(pf) else np.full(n, np.nan),
                 edyn=ef[:, 2] if len(ef) else np.full(n, np.nan))
    return d


def temp_report(d, dt_fs):
    T, t = d["T"], d["t_ps"]; n = d["n"]
    st = tstats(T)
    b250, b500 = blocks(T, 250), blocks(T, 500)
    final1ps = T[-int(1000 / dt_fs):].mean() if n >= int(1000 / dt_fs) else T.mean()
    final2ps = T[-int(2000 / dt_fs):].mean() if n >= int(2000 / dt_fs) else T.mean()
    bw = int(2000 / dt_fs)   # 2 ps block
    unc = block_unc(T, min(bw, max(1, n // 2)))
    out = {**st, "final_1ps_mean": float(final1ps), "final_2ps_mean": float(final2ps),
           "block250_means": [float(x) for x in b250], "block500_means": [float(x) for x in b500],
           "slope_full_K_per_ps": slope(t, T), "slope_final_half": slope(t[n // 2:], T[n // 2:]),
           "slope_final_quarter": slope(t[3 * n // 4:], T[3 * n // 4:]),
           "block_uncertainty": unc,
           "crossings": {f"first_above_{x}K": first_cross(t, T, x) for x in (600, 1000, 2000)},
           "fraction_above": {f"{x}K": float(np.mean(T > x)) for x in (600, 1000, 2000)},
           "longest_run_above": {f"{x}K": longest_run(T, x) for x in (600, 1000, 2000)}}
    # direct velocity crosscheck at terminal state
    if d["md"] is not None and len(d["cel"]):
        v = d["md"]["ratd"] @ np.array(d["cel"][-1]).T
        vcom = v.mean(axis=0)
        K = 0.5 * MASS_C_AMU * AMU_RY * ((v - vcom) ** 2).sum()
        out["direct_terminal_T_K"] = float(2 * K * RY_TO_K / NDOF)
        out["direct_confirms"] = bool(abs(out["direct_terminal_T_K"] - T[-1]) < max(8.0, 0.06 * abs(T[-1])))
    return out


def pressure_report(d):
    P = d["p_kbar"]; st = tstats(P)
    bw = min(2000, max(1, len(P) // 2))
    unc = block_unc(P, bw)
    st.update({"final_block_mean_kbar": unc.get("block_mean"), "block_sem_kbar": unc.get("sem"),
               "slope_kbar_per_step": slope(np.arange(len(P), dtype=float), P),
               "mean_GPa": st.get("mean", 0) / 10.0})
    if len(d["S"]):
        st["stress_mean_kbar"] = {nm: float(np.nanmean(d["S"][:, i, j])) for nm, (i, j) in
                                  zip(("sxx", "syy", "szz", "syz", "sxz", "sxy"),
                                      ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)))}
    return st


def cell_report(d, ref_vol):
    if not len(d["cel"]):
        return {}
    vol = np.array([abs(det3(m)) for m in d["cel"]])
    conds = np.array([cond2(m) for m in d["cel"]])
    H0 = d["cel"][0]
    dmins, nns = [], []
    for i in range(0, len(d["cel"]), max(1, len(d["cel"]) // 150)):
        if i < len(d["pos"]):
            a1, _, a3 = cc_dist(d["cel"][i], d["pos"][i]); dmins.append(a1); nns.append(a3)
    return {"volume_initial_bohr3": float(vol[0]), "volume_final_bohr3": float(vol[-1]),
            "volume_min": float(vol.min()), "volume_max": float(vol.max()),
            "volume_change_vs_initial_pct": float(100 * (vol[-1] - vol[0]) / vol[0]),
            "volume_change_vs_refvol_pct": float(100 * (vol[-1] - ref_vol) / ref_vol) if ref_vol else None,
            "max_abs_volume_dev_vs_refvol_pct": float(100 * np.max(np.abs(vol - ref_vol)) / ref_vol) if ref_vol else None,
            "condition_number_max": float(conds.max()),
            "determinant_min": float(vol.min()), "determinant_always_positive": bool(vol.min() > 0),
            "strain6_final": [float(x) for x in strain6(d["cel"][-1], H0)],
            "min_CC_bohr": float(min(dmins)) if dmins else None,
            "mean_NN_CC_bohr": float(np.mean(nns)) if nns else None}


def energy_report(d, dt_fs):
    ed = d["edyn"]; t = d["t_ps"]
    drift_ps = float(ed[-1] - ed[0]) / max(t[-1], 1e-12)
    return {"quantity": "edyn (.e col3 = ekint+ut)", "status": "PARTIAL_EXTENDED_ENERGY_DIAGNOSTIC",
            "note": "Nose chain energy not written by QE", "initial_Ry": float(ed[0]), "final_Ry": float(ed[-1]),
            "drift_total_Ry": float(ed[-1] - ed[0]), "drift_eV_per_atom_per_ps": drift_ps * EV_PER_RY / NAT,
            "bounded": bool(np.all(np.isfinite(ed)))}


def caimd_report(d):
    c = d["caimd"]
    if c is None:
        return {"caimd_present": False}
    r = {"caimd_present": True, "rows": c["_nrows"], "ncols_seen": c["_ncols"], "schema_26_ok": c["_ncols"] == [26],
         "event_count_final": int(c["events"][-1]) if len(c["events"]) else None,
         "n_event_flags_true": int(sum(c["flag"]))}
    comps = ("xx", "yy", "zz", "yz", "xz", "xy")
    for cp in comps:
        r[f"max_abs_nominal_{cp}"] = float(np.abs(c["nom_" + cp]).max())
        r[f"max_abs_actual_{cp}"] = float(np.abs(c["act_" + cp]).max())
        r[f"final_nominal_{cp}"] = float(c["nom_" + cp][-1])
        r[f"final_actual_{cp}"] = float(c["act_" + cp][-1])
    r["max_abs_nominal_minus_actual"] = float(np.max([np.abs(c["nom_" + cp] - c["act_" + cp]).max() for cp in comps]))
    r["max_velocity_jump_max"] = float(c["vjmax"].max())
    r["max_velocity_jump_rms"] = float(c["vjrms"].max())
    if len(d["cel"]):
        H0 = d["cel"][0]
        skews = [skew_norm(m, H0) for m in d["cel"]]
        vol = np.array([abs(det3(m)) for m in d["cel"]])
        r["max_skew_norm"] = float(max(skews))
        r["determinant_range"] = [float(vol.min()), float(vol.max())]
        r["condition_number_max"] = float(max(cond2(m) for m in d["cel"]))
    return r


def dump_temp_csv(path, d, dt_fs):
    T = d["T"]
    wcsv(path, ["nst", "time_ps", "eka_Ry", "T_instant_3Nm3_K", "roll250_K", "roll500_K"],
         [[int(d["nst"][i]), d["t_ps"][i], d["eka"][i], T[i], roll(T, 250)[i], roll(T, 500)[i]] for i in range(d["n"])])


def dump_blocks_csv(path, d, w):
    T = d["T"]; b = blocks(T, w)
    wcsv(path, ["block", "nst_start", "nst_end", "block_mean_T_K"],
         [[i + 1, i * w, min((i + 1) * w - 1, d["n"] - 1), float(x)] for i, x in enumerate(b)])


def dump_press_csv(path, d):
    P = d["p_kbar"]
    rows = []
    for i in range(d["n"]):
        s = d["S"][i] if i < len(d["S"]) else [[None] * 3] * 3
        rows.append([int(d["nst"][i]), d["t_ps"][i], P[i], P[i] / 10.0,
                     s[0][0], s[1][1], s[2][2], s[1][2], s[0][2], s[0][1]])
    wcsv(path, ["nst", "time_ps", "p_kbar", "p_GPa", "sxx", "syy", "szz", "syz", "sxz", "sxy"], rows)


# ============================================================ stage runners
def gate_nvt(d, dt_fs):
    tr = temp_report(d, dt_fs)
    fails = []
    if not (d["job_done"] and d["n"] >= 5000):
        fails.append("did not complete 5000 steps / JOB DONE")
    if d["error"]:
        fails.append(f"QE error {d['error']}")
    if not np.all(np.isfinite(d["T"])):
        fails.append("nonfinite temperature")
    if not (250 <= tr["final_2ps_mean"] <= 350):
        fails.append(f"final-2ps mean {tr['final_2ps_mean']:.1f} K outside 250-350")
    unc = tr["block_uncertainty"]
    if not (unc.get("sem") and abs((unc.get("block_mean") or 0) - 300) < 2 * unc["sem"] + 25):
        # compatible within 2 SEM or within 25 K of 300 (permissive for 8-atom)
        if not (275 <= (unc.get("block_mean") or 0) <= 325):
            fails.append(f"late-time block mean {unc.get('block_mean')} not compatible with 300 K")
    if tr["slope_final_quarter"] is not None and abs(tr["slope_final_quarter"]) > 200:
        fails.append(f"secular final-quarter slope {tr['slope_final_quarter']:.1f} K/ps")
    if tr["longest_run_above"]["1000K"] >= 100:
        fails.append("interval above 1000 K persists >=100 steps")
    return tr, fails


def gate_npt(d, dt_fs, ref_vol):
    tr = temp_report(d, dt_fs); pr = pressure_report(d); cr = cell_report(d, ref_vol)
    fails = []
    if not (d["job_done"] and d["n"] >= 5000):
        fails.append("did not complete 5000 steps / JOB DONE")
    if d["error"]:
        fails.append(f"QE error {d['error']}")
    if d["caimd"] is not None:
        fails.append("CAIMD .caimd present but must be disabled")
    if not d["nose_ion"]:
        fails.append("ionic Nose not active in stdout")
    if not cr.get("determinant_always_positive", False):
        fails.append("cell determinant non-positive")
    if cr.get("condition_number_max", 99) > 5:
        fails.append(f"cell condition number {cr.get('condition_number_max')} too large")
    if not (250 <= tr["final_2ps_mean"] <= 350):
        fails.append(f"final-2ps T mean {tr['final_2ps_mean']:.1f} K outside 250-350")
    # pressure: within 2 block SEM OR 20 kbar
    pbm, psem = pr.get("final_block_mean_kbar"), pr.get("block_sem_kbar")
    press_ok = (pbm is not None and (abs(pbm) < 20.0 or (psem and abs(pbm) < 2 * psem)))
    if not press_ok:
        fails.append(f"final-block pressure {pbm} kbar not compatible with 0 (>20 kbar and >2 SEM)")
    if cr.get("max_abs_volume_dev_vs_refvol_pct") is not None and cr["max_abs_volume_dev_vs_refvol_pct"] > 5.0:
        fails.append(f"volume deviates {cr['max_abs_volume_dev_vs_refvol_pct']:.1f}% from vc-relaxed (>5%)")
    if tr["temperature"] if False else tr["max"] > 3000 and tr["final_2ps_mean"] > 1000:
        fails.append("sustained several-thousand-K heating")
    return tr, pr, cr, fails


def gate_caimd_zero(d, dt_fs):
    tr = temp_report(d, dt_fs); cr = caimd_report(d)
    fails = []
    if not (d["job_done"] and d["n"] >= 3000):
        fails.append("did not complete 3000 steps / JOB DONE")
    if d["error"]:
        fails.append(f"QE error {d['error']}")
    if not cr.get("caimd_present"):
        fails.append("no .caimd file")
    else:
        if cr.get("n_event_flags_true", 0) < 1:
            fails.append("no CAIMD events fired")
        if any(cr.get(f"max_abs_nominal_{c}", 1) != 0.0 for c in ("xx", "yy", "zz", "yz", "xz", "xy")):
            fails.append("nominal strain not all zero")
        if cr.get("max_abs_nominal_minus_actual", 1) > 1e-10:
            fails.append(f"strain tracking error {cr.get('max_abs_nominal_minus_actual'):.2e} > 1e-10")
        if not cr.get("determinant_range", [1, 1]) or cr.get("determinant_range", [-1])[0] <= 0:
            fails.append("cell determinant non-positive")
        if cr.get("max_velocity_jump_max", 1) > 1e-6:
            fails.append(f"velocity-transform residual {cr.get('max_velocity_jump_max'):.2e} > 1e-6")
    if not np.all(np.isfinite(d["T"])):
        fails.append("nonfinite temperature")
    f15 = d["T"][-int(1500 / dt_fs):].mean() if d["n"] >= int(1500 / dt_fs) else d["T"].mean()
    if not (250 <= f15 <= 350):
        fails.append(f"final-1.5ps T mean {f15:.1f} K outside 250-350")
    if tr["max"] > 3000 and tr["final_1ps_mean"] > 1000:
        fails.append("sustained heating")
    tr["final_1p5ps_mean"] = float(f15)
    return tr, cr, fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dt-fs", type=float, default=1.0)
    ap.add_argument("--ref-volume", type=float, default=0.0)
    ap.add_argument("--rc", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    d = load_md(a.dir, a.prefix, a.dt_fs)
    if d["n"] == 0:
        json.dump({"stage": a.stage, "classification": "MISSING_OUTPUT", "gate_pass": False},
                  open(os.path.join(a.out, "analysis_summary.json"), "w"), indent=2)
        print("MISSING_OUTPUT"); sys.exit(2)

    summary = {"stage": a.stage, "qe_return_code": a.rc, "job_done": d["job_done"], "version": d["version"],
               "error": d["error"], "completed_steps": d["n"], "nose_ion_active": d["nose_ion"]}

    if a.stage == "nvt":
        tr, fails = gate_nvt(d, a.dt_fs)
        dump_temp_csv(os.path.join(a.out, "temperature.csv"), d, a.dt_fs)
        dump_blocks_csv(os.path.join(a.out, "temperature_blocks.csv"), d, 500)
        dump_press_csv(os.path.join(a.out, "pressure.csv"), d)
        wcsv(os.path.join(a.out, "energy_diagnostics.csv"),
             ["quantity", "status", "initial_Ry", "final_Ry", "drift_eV_per_atom_per_ps", "bounded"],
             [[er["quantity"], er["status"], er["initial_Ry"], er["final_Ry"], er["drift_eV_per_atom_per_ps"], er["bounded"]]
              for er in [energy_report(d, a.dt_fs)]])
        summary.update({"temperature": tr, "pressure": pressure_report(d), "energy": energy_report(d, a.dt_fs)})
        gate = not fails
    elif a.stage == "npt":
        tr, pr, cr, fails = gate_npt(d, a.dt_fs, a.ref_volume)
        dump_temp_csv(os.path.join(a.out, "temperature.csv"), d, a.dt_fs)
        dump_blocks_csv(os.path.join(a.out, "temperature_blocks.csv"), d, 500)
        dump_press_csv(os.path.join(a.out, "pressure.csv"), d)
        P = d["p_kbar"]; wcsv(os.path.join(a.out, "pressure_blocks.csv"),
             ["block", "nst_start", "nst_end", "block_mean_p_kbar"],
             [[i + 1, i * 500, min((i + 1) * 500 - 1, d["n"] - 1), float(x)] for i, x in enumerate(blocks(P, 500))])
        vol = np.array([abs(det3(m)) for m in d["cel"]]) if len(d["cel"]) else np.array([])
        ct = (d["cel_steps"] - d["cel_steps"][0]) * a.dt_fs / 1000.0 if len(d["cel"]) else np.array([])
        wcsv(os.path.join(a.out, "volume.csv"), ["step", "time_ps", "volume_bohr3"],
             [[int(d["cel_steps"][i]), float(ct[i]), float(vol[i])] for i in range(len(vol))])
        wcsv(os.path.join(a.out, "cell_diagnostics.csv"), ["step", "time_ps", "volume_bohr3", "condition_number"],
             [[int(d["cel_steps"][i]), float(ct[i]), float(vol[i]), float(cond2(d["cel"][i]))] for i in range(len(vol))])
        wcsv(os.path.join(a.out, "energy_diagnostics.csv"),
             ["quantity", "status", "drift_eV_per_atom_per_ps", "bounded"],
             [[er["quantity"], er["status"], er["drift_eV_per_atom_per_ps"], er["bounded"]] for er in [energy_report(d, a.dt_fs)]])
        summary.update({"temperature": tr, "pressure": pr, "cell": cr, "caimd_disabled": d["caimd"] is None,
                        "energy": energy_report(d, a.dt_fs)})
        gate = not fails
    elif a.stage == "caimd_zero":
        tr, cr, fails = gate_caimd_zero(d, a.dt_fs)
        dump_temp_csv(os.path.join(a.out, "temperature.csv"), d, a.dt_fs)
        c = d["caimd"]
        if c is not None:
            wcsv(os.path.join(a.out, "strain_diagnostics.csv"),
                 ["step", "time_ps", "events", "nom_xx", "act_xx", "max_abs_actual6", "vjmax", "vjrms", "temperature_K", "volume"],
                 [[int(c["step"][i]), float(c["time"][i]), int(c["events"][i]), float(c["nom_xx"][i]), float(c["act_xx"][i]),
                   float(max(abs(c["act_" + q][i]) for q in ("xx", "yy", "zz", "yz", "xz", "xy"))),
                   float(c["vjmax"][i]), float(c["vjrms"][i]), float(c["temperature"][i]), float(c["volume"][i])]
                  for i in range(c["_nrows"])])
        summary.update({"temperature": tr, "caimd": cr})
        gate = not fails
    elif a.stage == "caimd_nonzero":
        tr, cr, fails = gate_caimd_zero(d, a.dt_fs)   # reuse numerical checks; relax the "nominal all zero"
        # for nonzero: nominal xx SHOULD be nonzero and tracked
        fails = [f for f in fails if "nominal strain not all zero" not in f and "3000" not in f]
        if not (d["job_done"] and d["n"] >= 2000):
            fails.append("did not complete 2000 steps / JOB DONE")
        c = d["caimd"]
        exp_events = 200
        exp_cum = 0.01   # 2000 steps * 5e-6/step (implementation: n_events * stride * step_strain = 200*10*5e-6)
        if c is not None:
            tracking = float(np.max(np.abs(c["nom_xx"] - c["act_xx"]))) if len(c["nom_xx"]) else None
            cr["stage4B_event_count"] = int(c["events"][-1]) if len(c["events"]) else None
            cr["stage4B_expected_events"] = exp_events
            cr["stage4B_final_nominal_xx"] = float(c["nom_xx"][-1]) if len(c["nom_xx"]) else None
            cr["stage4B_final_actual_xx"] = float(c["act_xx"][-1]) if len(c["act_xx"]) else None
            cr["stage4B_expected_cumulative_xx"] = exp_cum
            cr["stage4B_xx_tracking_error"] = tracking
            cr["prompt_estimate_note"] = ("prompt stated ~0.1%; implementation event increment is "
                                          "stride*step_strain=5e-5, so 200 events -> ~1.0% cumulative "
                                          "(= per-step rate 5e-6 * 2000 steps). Keywords used exactly as specified.")
            if tracking is not None and tracking > 1e-9:
                fails.append(f"xx strain tracking error {tracking:.2e} > 1e-9")
            dump_temp_csv(os.path.join(a.out, "temperature.csv"), d, a.dt_fs)
            wcsv(os.path.join(a.out, "strain_diagnostics.csv"),
                 ["step", "time_ps", "events", "flag", "nom_xx", "act_xx", "max_abs_offdiag_actual", "vjmax", "vjrms", "temperature_K", "volume"],
                 [[int(c["step"][i]), float(c["time"][i]), int(c["events"][i]), 1 if c["flag"][i] else 0,
                   float(c["nom_xx"][i]), float(c["act_xx"][i]),
                   float(max(abs(c["act_" + q][i]) for q in ("yy", "zz", "yz", "xz", "xy"))),
                   float(c["vjmax"][i]), float(c["vjrms"][i]), float(c["temperature"][i]), float(c["volume"][i])]
                  for i in range(c["_nrows"])])
            dump_press_csv(os.path.join(a.out, "stress.csv"), d)
        summary.update({"temperature": tr, "caimd": cr})
        gate = not fails
    else:
        print("unknown stage"); sys.exit(3)

    summary["gate_failures"] = fails
    summary["gate_pass"] = bool(gate)
    json.dump(summary, open(os.path.join(a.out, "analysis_summary.json"), "w"), indent=2, default=float)
    print(json.dumps({"stage": a.stage, "gate_pass": gate, "failures": fails,
                      "final_2ps_T": summary.get("temperature", {}).get("final_2ps_mean")}, indent=2, default=float))
    sys.exit(0 if gate else 1)


if __name__ == "__main__":
    main()
