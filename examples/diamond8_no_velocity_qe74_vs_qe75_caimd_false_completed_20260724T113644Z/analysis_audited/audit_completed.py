#!/usr/bin/env python3
"""Independent audit of the completed Diamond8 no-velocity QE 7.4.1 vs QE 7.5 job 50160152.

Reads only the local mirror; writes new products under analysis_audited/ and
evidence_audited/. Raw mirrored outputs are never modified.

Verified column semantics (from PW/src/vcsmd.f90 at af7ebfc, lines 410-419):
  .e    : ut, ekint, edyn, pv, nst
  .eal  : uta, eka, eta, utl, ekla, etl, nst
  .ave  : avu, avk, nst
  .p    : press, p, avp, nst          [pressure in Ry/bohr^3]
  .avec : |a|,|b|,|c|, theta12, theta23, theta31, nst
  .tv   : vcell, tnew, nst            [tnew from RUNNING-AVERAGE avk, 3(N+1) dof]
  .md   : terminal state only (seqopn rewrite): rati, ratd, rat2d, rat2di, tauold /
          aveci, avecd, avec2d, avec2di / ...
"""
import hashlib
import json
import math
import os
import re
import sys

import numpy as np

MIRROR = sys.argv[1]
ANA = os.path.join(MIRROR, "analysis_audited")
EVD = os.path.join(MIRROR, "evidence_audited")
os.makedirs(ANA, exist_ok=True)
os.makedirs(EVD, exist_ok=True)

RY_TO_K = 157887.51240116          # 1/k_B in Ry
RY_BOHR3_TO_KBAR = 147105.078      # 1 Ry/bohr^3 -> kbar  (ry_kbar in QE)
BOHR_TO_ANG = 0.529177210903
NAT = 8
DT_FS = 2.0
MASS_C_AMU = 12.0107
AMU_RY = 911.4442421             # amu -> Ry atomic mass units (2*m_e = 1)

RUNS = {
    "qe74": {"dir": "qe74_official", "prefix": "diamond8_qe74_novel",
             "label": "Official QE 7.4.1", "input": "inputs/qe74_no_velocity.in"},
    "qe75": {"dir": "qe75_af7ebfc_caimd_false", "prefix": "diamond8_qe75_novel",
             "label": "QE 7.5 af7ebfc caimd=.false.", "input": "inputs/qe75_caimd_false_no_velocity.in"},
}

NUM = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")


def rd(path):
    out = []
    if not os.path.exists(path):
        return np.zeros((0, 0))
    for line in open(path, errors="replace"):
        v = [float(t.replace("D", "E").replace("d", "e")) for t in NUM.findall(line)]
        if v:
            out.append(v)
    n = min(len(r) for r in out) if out else 0
    return np.array([r[:n] for r in out]) if out else np.zeros((0, 0))


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(65536), b""):
            h.update(c)
    return h.hexdigest()


def text(p):
    return open(p, errors="replace").read() if os.path.exists(p) else ""


def rolling(a, w):
    if len(a) < w:
        return np.full(len(a), np.nan)
    c = np.convolve(a, np.ones(w) / w, mode="valid")
    return np.concatenate([np.full(w - 1, np.nan), c])


def blocks(a, size=100):
    return [float(np.mean(a[i:i + size])) for i in range(0, len(a), size) if len(a[i:i + size])]


def drift(t_ps, y):
    if len(y) < 3:
        return None
    return float(np.polyfit(t_ps, y, 1)[0])


def first_above(t_ps, y, thr):
    idx = np.nonzero(y > thr)[0]
    if not len(idx):
        return None
    i = int(idx[0])
    return {"index": i, "step_nst": i, "time_ps": float(t_ps[i]), "value": float(y[i])}


def longest_run(y, thr):
    best = cur = 0
    for v in y:
        cur = cur + 1 if v > thr else 0
        best = max(best, cur)
    return int(best)


def frac_above(y, thr):
    return float(np.mean(y > thr))


def parse_cel(path):
    """.cel: '<step> <time>' header then 3 rows. Rows are lattice vectors (bohr)."""
    steps, mats = [], []
    toks = [l.split() for l in text(path).splitlines() if l.strip()]
    i = 0
    while i + 3 < len(toks) + 1 and i + 3 <= len(toks):
        if len(toks[i]) == 2:
            try:
                s = int(float(toks[i][0]))
                m = [[float(x) for x in toks[i + r]] for r in (1, 2, 3)]
            except (ValueError, IndexError):
                i += 1
                continue
            steps.append(s)
            mats.append(m)
            i += 4
        else:
            i += 1
    return np.array(steps), np.array(mats)


def parse_pos(path, nat=NAT):
    steps, pos = [], []
    toks = [l.split() for l in text(path).splitlines() if l.strip()]
    i = 0
    while i < len(toks):
        if len(toks[i]) == 2:
            try:
                s = int(float(toks[i][0]))
                blk = [[float(x) for x in toks[i + 1 + k]] for k in range(nat)]
            except (ValueError, IndexError):
                i += 1
                continue
            steps.append(s)
            pos.append(blk)
            i += 1 + nat
        else:
            i += 1
    return np.array(steps), np.array(pos)


def parse_stress_stdout(path):
    """Per-step 'total stress ... P= <kbar>' plus the 3x3 kbar block."""
    lines = text(path).splitlines()
    P, comp = [], []
    for i, line in enumerate(lines):
        if "stress" in line and "P=" in line:
            m = re.search(r"P\s*=\s*([-+0-9.EeDd]+)", line)
            p = float(m.group(1).replace("D", "E")) if m else np.nan
            mat = []
            ok = True
            for j in range(1, 4):
                if i + j >= len(lines):
                    ok = False
                    break
                v = [float(t.replace("D", "E")) for t in NUM.findall(lines[i + j])]
                if len(v) < 6:
                    ok = False
                    break
                mat.append(v[3:6])
            P.append(p)
            comp.append(mat if ok else [[np.nan] * 3] * 3)
    return np.array(P), np.array(comp)


def parse_stdout_T(path):
    rows = []
    for line in text(path).splitlines():
        if "Ekin =" in line and " T =" in line:
            v = [float(t.replace("D", "E")) for t in NUM.findall(line)]
            if len(v) >= 3:
                rows.append(v[:3])
    return np.array(rows)


def parse_md_terminal(path, nat=NAT):
    """Line 1: rati(3,nat), ratd(3,nat), rat2d(3,nat), rat2di(3,nat), tauold(3,nat,3).
       Line 2: aveci(3,3), avecd(3,3), avec2d(3,3), avec2di(3,3)."""
    lines = [l for l in text(path).splitlines() if l.strip()]
    if len(lines) < 2:
        return None
    # tokens can wrap across physical lines; rebuild by scanning all numbers in order
    allnum = [float(t.replace("D", "E")) for t in NUM.findall(text(path))]
    n3 = 3 * nat
    need = n3 * 4 + 3 * nat * 3
    if len(allnum) < need + 36:
        return None
    rati = np.array(allnum[0:n3]).reshape(nat, 3)
    ratd = np.array(allnum[n3:2 * n3]).reshape(nat, 3)
    off = need
    aveci = np.array(allnum[off:off + 9]).reshape(3, 3)
    return {"rati": rati, "ratd": ratd, "aveci": aveci}


def cond_and_sv(m):
    s = np.linalg.svd(m, compute_uv=False)
    return float(s[0] / s[-1]), s


def min_max_cc(cell_rows, frac):
    """Min/max/mean-NN periodic C-C distance (bohr). cell_rows: lattice vectors as rows."""
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    cartshift = shifts @ cell_rows
    n = len(frac)
    cart = frac @ cell_rows
    dmin, dmax, nn = np.inf, 0.0, []
    for a in range(n):
        best = np.inf
        for b in range(n):
            d = cart[b] - cart[a] + cartshift
            r = np.sqrt((d * d).sum(axis=1))
            r = r[r > 1e-8]
            if len(r):
                dmin = min(dmin, r.min())
                dmax = max(dmax, r.max())
                if b != a:
                    best = min(best, r.min())
        nn.append(best)
    return float(dmin), float(dmax), float(np.mean(nn))


def strain6(Hcol, Hcol0):
    """Symmetric small strain E = sym(F) - I with F = H H0^-1.

    .cel stores the cell with lattice vectors as COLUMNS (QE `avec`), verified
    against the input CELL_PARAMETERS rows, so F acts on Cartesian columns."""
    F = Hcol @ np.linalg.inv(Hcol0)
    E = 0.5 * (F + F.T) - np.eye(3)
    return [E[0, 0], E[1, 1], E[2, 2], E[1, 2], E[0, 2], E[0, 1]]


def stats(a):
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    if not len(a):
        return {"n": 0}
    return {"n": int(len(a)), "initial": float(a[0]), "final": float(a[-1]),
            "mean": float(a.mean()), "std": float(a.std()), "median": float(np.median(a)),
            "min": float(a.min()), "max": float(a.max()),
            "first100_mean": float(a[:100].mean()),
            "final100_mean": float(a[-100:].mean()),
            "final200_mean": float(a[-200:].mean())}


def wcsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join("" if v is None or (isinstance(v, float) and not np.isfinite(v))
                              else (f"{v:.10g}" if isinstance(v, float) else str(v))
                              for v in r) + "\n")


# =====================================================================
# Phase 4: input contract audit
# =====================================================================
def strip_comments(t):
    return [l.split("!")[0].rstrip() for l in t.splitlines() if l.split("!")[0].strip()]


def parse_input(t):
    nml, cards, cur = {}, [], None
    for line in strip_comments(t):
        s = line.strip()
        if s.startswith("&"):
            cur = s[1:].strip().upper(); nml[cur] = {}; continue
        if s == "/":
            cur = None; continue
        if cur:
            for part in s.split(","):
                if "=" in part:
                    k, v = part.split("=", 1)
                    nml[cur][k.strip().lower()] = v.strip().rstrip(",").strip()
        else:
            cards.append(s)
    return nml, cards


def cards_struct(cards):
    out = {"species": [], "cell": [], "pos": [], "kpoints": None,
           "cell_units": None, "pos_units": None, "atomic_velocities": False}
    mode = None
    for line in cards:
        u = line.upper().split()
        h = u[0] if u else ""
        if h == "ATOMIC_SPECIES": mode = "sp"; continue
        if h == "CELL_PARAMETERS": mode = "ce"; out["cell_units"] = u[1].lower() if len(u) > 1 else None; continue
        if h == "ATOMIC_POSITIONS": mode = "po"; out["pos_units"] = u[1].lower() if len(u) > 1 else None; continue
        if h == "K_POINTS": mode = "kp"; out["kpoints"] = " ".join(u[1:]).lower(); continue
        if h == "ATOMIC_VELOCITIES": out["atomic_velocities"] = True; mode = "ve"; continue
        f = line.split()
        if mode == "sp" and len(f) >= 3: out["species"].append((f[0], float(f[1]), f[2]))
        elif mode == "ce" and len(f) == 3: out["cell"].append([float(x) for x in f])
        elif mode == "po" and len(f) >= 4: out["pos"].append((f[0], [float(x) for x in f[1:4]]))
    return out


def numv(v):
    try:
        return float(str(v).lower().replace("d", "e").strip().strip("'\""))
    except ValueError:
        return str(v).strip().strip("'\"").lower()


t74 = text(os.path.join(MIRROR, RUNS["qe74"]["input"]))
t75 = text(os.path.join(MIRROR, RUNS["qe75"]["input"]))
n74, c74 = parse_input(t74)
n75, c75 = parse_input(t75)
k74, k75 = cards_struct(c74), cards_struct(c75)

FORBID_VEL = ["atomic_velocities", "ion_velocities", "from_input", "tv0rd", "starting_velocities"]
FORBID_CAIMD = ["caimd_stride", "caimd_step_strain", "caimd_mask"]

no_vel = {"note": "checks run on the COMPLETE input text, case-insensitive", "runs": {}}
for tag, tt in (("qe74", t74), ("qe75", t75)):
    hits = {p: [i + 1 for i, l in enumerate(tt.splitlines()) if re.search(p, l, re.I)] for p in FORBID_VEL}
    ctl = parse_input(tt)[0]
    no_vel["runs"][tag] = {
        "forbidden_keyword_hits": {k: v for k, v in hits.items() if v},
        "atomic_velocities_card": cards_struct(parse_input(tt)[1])["atomic_velocities"],
        "restart_mode": numv(ctl.get("CONTROL", {}).get("restart_mode", "")),
        "ion_temperature": numv(ctl.get("IONS", {}).get("ion_temperature", "")),
        "tempw": numv(ctl.get("IONS", {}).get("tempw", "")),
    }
no_vel["explicit_velocities_absent_in_both"] = all(
    not r["forbidden_keyword_hits"] and not r["atomic_velocities_card"] for r in no_vel["runs"].values())
json.dump(no_vel, open(os.path.join(ANA, "no_velocity_audit.json"), "w"), indent=2)

caimd_audit = {"runs": {}}
for tag, tt in (("qe74", t74), ("qe75", t75)):
    eff = strip_comments(tt)
    lines = [l.strip() for l in eff if "caimd" in l.lower()]
    caimd_audit["runs"][tag] = {
        "caimd_lines_effective_input": lines,
        "forbidden_parameters": [p for p in FORBID_CAIMD if any(re.search(p, l, re.I) for l in eff)],
        "active_true": [l for l in lines if re.search(r"caimd\s*=\s*\.?true", l, re.I)],
    }
caimd_audit["qe74_has_no_caimd_keyword"] = not caimd_audit["runs"]["qe74"]["caimd_lines_effective_input"]
q75l = caimd_audit["runs"]["qe75"]["caimd_lines_effective_input"]
caimd_audit["qe75_only_caimd_false"] = (len(q75l) == 1 and bool(re.search(r"caimd\s*=\s*\.false\.", q75l[0], re.I)))
json.dump(caimd_audit, open(os.path.join(ANA, "caimd_keyword_audit.json"), "w"), indent=2)

ALLOWED = {("CONTROL", "prefix"), ("CONTROL", "outdir"), ("CELL", "caimd")}
diffs = {}
for sec in sorted(set(n74) | set(n75)):
    a, b = n74.get(sec, {}), n75.get(sec, {})
    for key in sorted(set(a) | set(b)):
        va, vb = a.get(key, "<absent>"), b.get(key, "<absent>")
        if numv(va) != numv(vb) and (sec, key) not in ALLOWED:
            diffs[f"{sec}.{key}"] = {"qe74": va, "qe75": vb}

H0_in = np.array(k74["cell"])
frac0 = np.array([p[1] for p in k74["pos"]])
dmin0, dmax0, nn0 = min_max_cc(H0_in, frac0)
ident = {
    "namelist_differences_outside_allowed_set": diffs,
    "allowed_differences": sorted(f"{s}.{k}" for s, k in ALLOWED),
    "cell_identical": k74["cell"] == k75["cell"],
    "positions_identical": k74["pos"] == k75["pos"],
    "species_masses_identical": k74["species"] == k75["species"],
    "kpoints_identical": k74["kpoints"] == k75["kpoints"],
    "units_identical": (k74["cell_units"] == k75["cell_units"] and k74["pos_units"] == k75["pos_units"]),
    "n_atoms": {"qe74": len(k74["pos"]), "qe75": len(k75["pos"])},
    "structure": {"cell_determinant_bohr3": float(np.linalg.det(H0_in)),
                  "min_CC_bohr": dmin0, "min_CC_angstrom": dmin0 * BOHR_TO_ANG,
                  "max_CC_bohr": dmax0, "mean_NN_CC_bohr": nn0},
    "input_sha256": {"qe74": sha256(os.path.join(MIRROR, RUNS["qe74"]["input"])),
                     "qe75": sha256(os.path.join(MIRROR, RUNS["qe75"]["input"]))},
    "pseudo_sha256": sha256(os.path.join(MIRROR, "pseudo", "C.UPF")),
    "pseudo_expected": "9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639",
}
ident["pseudo_matches"] = ident["pseudo_sha256"] == ident["pseudo_expected"]
ident["physical_settings_identical"] = (
    not diffs and ident["cell_identical"] and ident["positions_identical"]
    and ident["species_masses_identical"] and ident["kpoints_identical"] and ident["units_identical"])
json.dump(ident, open(os.path.join(ANA, "input_identity_audit.json"), "w"), indent=2)

import difflib
with open(os.path.join(ANA, "input_normalized.diff"), "w") as fh:
    fh.write("\n".join(difflib.unified_diff(
        strip_comments(t74), strip_comments(t75),
        fromfile="qe74_no_velocity.in (normalized)",
        tofile="qe75_caimd_false_no_velocity.in (normalized)", lineterm="")) + "\n")

# =====================================================================
# Phases 5,6,8-12 per run
# =====================================================================
D = {}
for tag, meta in RUNS.items():
    base = os.path.join(MIRROR, meta["dir"])
    pref = meta["prefix"]
    out = os.path.join(base, "qe.out")
    eal = rd(os.path.join(base, pref + ".eal"))
    pfile = rd(os.path.join(base, pref + ".p"))
    tvf = rd(os.path.join(base, pref + ".tv"))
    efile = rd(os.path.join(base, pref + ".e"))
    avef = rd(os.path.join(base, pref + ".ave"))
    csteps, cmats = parse_cel(os.path.join(base, "tmp", pref + ".cel"))
    psteps, ppos = parse_pos(os.path.join(base, "tmp", pref + ".pos"))
    Pstd, Scomp = parse_stress_stdout(out)
    stdT = parse_stdout_T(out)
    md = parse_md_terminal(os.path.join(base, "tmp", pref + ".md"))
    txt = text(out)

    nst = eal[:, 6].astype(int)
    t_ps = nst * DT_FS / 1000.0
    eka = eal[:, 1]
    ekla = eal[:, 4]
    T21 = 2 * eka * RY_TO_K / 21.0
    T27 = 2 * eka * RY_TO_K / 27.0
    tv_T = tvf[:, 1]
    vcell = tvf[:, 0]
    avk = avef[:, 1]
    T_avk27 = 2 * avk * RY_TO_K / 27.0
    p_ry = pfile[:, 1]
    p_kbar = p_ry * RY_BOHR3_TO_KBAR
    ut, ekint, edyn, pv = efile[:, 0], efile[:, 1], efile[:, 2], efile[:, 3]

    ver = re.search(r"Program PWSCF v\.(\S+)", txt)
    err = re.search(r"Error in routine\s+(\S+)", txt)
    D[tag] = dict(
        nst=nst, t_ps=t_ps, eka=eka, ekla=ekla, T21=T21, T27=T27, tv_T=tv_T,
        T_avk27=T_avk27, vcell=vcell, p_kbar=p_kbar, p_ry=p_ry,
        ut=ut, ekint=ekint, edyn=edyn, pv=pv,
        csteps=csteps, cmats=cmats, psteps=psteps, ppos=ppos,
        Pstd=Pstd, Scomp=Scomp, stdT=stdT, md=md,
        version=ver.group(1) if ver else None,
        job_done="JOB DONE" in txt,
        error=err.group(1) if err else None,
        scf_fail=txt.count("convergence NOT achieved"),
        label=meta["label"],
    )

# ---------- run completeness ----------
rc_txt = text(os.path.join(MIRROR, "evidence", "run_return_codes.txt"))
rcs = dict(re.findall(r"(\w+)=(\S+)", rc_txt))
completeness = {}
for tag in RUNS:
    d = D[tag]
    n = len(d["nst"])
    cls = ("COMPLETE" if (d["job_done"] and n >= 1000 and not d["error"])
           else "PARTIAL" if n > 0 else "MISSING_OUTPUT")
    completeness[tag] = {
        "return_code": int(rcs.get("RC_QE74" if tag == "qe74" else "RC_QE75_CAIMD_FALSE", -1)),
        "slurm_step": "50160152.0" if tag == "qe74" else "50160152.1",
        "slurm_step_state": "COMPLETED 0:0",
        "job_done": bool(d["job_done"]),
        "runtime_version": d["version"],
        "completed_ionic_steps": int(n),
        "requested_steps": 1000,
        "requested_time_ps": 2.0,
        "completed_time_ps": float(d["t_ps"][-1]) if n else 0.0,
        "final_valid_nst": int(d["nst"][-1]) if n else None,
        "first_error": d["error"],
        "scf_failure_count": int(d["scf_fail"]),
        "rows": {"eal": n, "p": len(d["p_kbar"]), "tv": len(d["tv_T"]),
                 "cel": len(d["csteps"]), "pos": len(d["psteps"]),
                 "stdout_stress": len(d["Pstd"]), "stdout_T": len(d["stdT"])},
        "classification": cls,
    }
completeness["internal_post_processing_rc"] = int(rcs.get("RC_ANALYSIS", -1))
completeness["wrapper_completed"] = True
json.dump(completeness, open(os.path.join(ANA, "run_completeness.json"), "w"), indent=2)

# ---------- executable provenance ----------
prov = {"job_id": "50160152", "node": "lrdn0416", "mpi_launcher": "srun --mpi=pmix_v3 -n 1",
        "ranks": 1, "omp_threads": 8, "runs": {}}
exe_txt = text(os.path.join(MIRROR, "evidence", "executable_checksums.txt"))
for tag in RUNS:
    envf = text(os.path.join(MIRROR, "evidence", f"{tag}_environment.txt"))
    path = re.search(r"^executable=(\S+)", envf, re.M)
    sha = re.findall(r"([0-9a-f]{64})\s+(\S+)", exe_txt)
    mine = [s for s, p in sha if (tag == "qe74") == ("spack" in p)]
    prov["runs"][tag] = {
        "executable_path": path.group(1) if path else None,
        "sha256_runtime": mine[0] if mine else None,
        "version_banner": D[tag]["version"],
        "version_probe": "no output (pw.x --version not supported in this build/launch)",
        "slurm_step": completeness[tag]["slurm_step"],
        "hostname": (re.search(r"^hostname=(\S+)", envf, re.M) or [None, None])[1]
        if re.search(r"^hostname=(\S+)", envf, re.M) else None,
        "module_list": [l.strip() for l in re.findall(r"^\s*\d+\)\s+(\S+)", envf, re.M)],
    }
prov["runs"]["qe75"]["expected_sha256"] = "89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d"
prov["runs"]["qe75"]["sha256_matches_expected"] = (
    prov["runs"]["qe75"]["sha256_runtime"] == prov["runs"]["qe75"]["expected_sha256"])
json.dump(prov, open(os.path.join(EVD, "executable_provenance.json"), "w"), indent=2)

# ---------- per-run temperature / pressure CSV + crosschecks ----------
tcross, pcross = [], []
for tag in RUNS:
    d = D[tag]
    n = len(d["nst"])
    rows = []
    for i in range(n):
        rows.append([int(d["nst"][i]), d["t_ps"][i], d["eka"][i], d["T21"][i], d["T27"][i],
                     d["tv_T"][i], d["T_avk27"][i],
                     d["stdT"][i, 1] if i < len(d["stdT"]) else None])
    wcsv(os.path.join(ANA, f"{tag}_temperature.csv"),
         ["nst", "time_ps", "eka_Ry", "T_instant_3Nm3_K", "T_instant_3Np1_K",
          "tv_reported_K", "T_from_running_avk_3Np1_K", "stdout_T_K"], rows)

    prows = []
    for i in range(n):
        s = d["Scomp"][i] if i < len(d["Scomp"]) else np.full((3, 3), np.nan)
        prows.append([int(d["nst"][i]), d["t_ps"][i], d["p_ry"][i], d["p_kbar"][i],
                      d["p_kbar"][i] / 10.0,
                      d["Pstd"][i] if i < len(d["Pstd"]) else None,
                      s[0][0], s[1][1], s[2][2], s[1][2], s[0][2], s[0][1]])
    wcsv(os.path.join(ANA, f"{tag}_pressure_stress.csv"),
         ["nst", "time_ps", "p_Ry_bohr3", "p_kbar", "p_GPa", "stdout_P_kbar",
          "sxx_kbar", "syy_kbar", "szz_kbar", "syz_kbar", "sxz_kbar", "sxy_kbar"], prows)

    # temperature source crosscheck incl. offsets
    for lag in (0, 1):
        a = d["T27"][:n - lag] if lag else d["T27"]
        b = d["tv_T"][lag:] if lag else d["tv_T"]
        m = min(len(a), len(b))
        tcross.append([tag, "T_instant_3Np1 vs .tv", lag,
                       float(np.abs(a[:m] - b[:m]).max())])
        a2 = d["T_avk27"][:n - lag] if lag else d["T_avk27"]
        tcross.append([tag, "T_running_avk_3Np1 vs .tv", lag,
                       float(np.abs(a2[:m] - b[:m]).max())])
    if len(d["stdT"]):
        m = min(len(d["stdT"]), n)
        tcross.append([tag, "stdout_T vs .tv", 0,
                       float(np.abs(d["stdT"][:m, 1] - d["tv_T"][:m]).max())])
        tcross.append([tag, "stdout_Ekin vs .e ekint", 0,
                       float(np.abs(d["stdT"][:m, 0] - d["ekint"][:m]).max())])
    # direct velocity temperature at terminal state
    if d["md"]:
        Hn = d["cmats"][-1]                      # terminal cell (rows = lattice vectors)
        ratd = d["md"]["ratd"]                   # lattice-coordinate velocities
        for name, H in (("terminal_cel", Hn), ("md_aveci", d["md"]["aveci"].T)):
            v = ratd @ H                          # Cartesian velocity = ratd . H(rows)
            vcom = v.mean(axis=0)
            m_ry = MASS_C_AMU * AMU_RY
            K = 0.5 * m_ry * ((v - vcom) ** 2).sum()
            T = 2 * K * RY_TO_K / 21.0
            tcross.append([tag, f"direct_md_ratd_x_{name}_3Nm3", 0, float(T)])
        tcross.append([tag, "eal_T_instant_3Nm3_terminal", 0, float(d["T21"][-1])])

    for lag in (0, 1):
        a = d["p_kbar"][:n - lag] if lag else d["p_kbar"]
        b = d["Pstd"][lag:] if lag else d["Pstd"]
        m = min(len(a), len(b))
        if m:
            pcross.append([tag, ".p(Ry/bohr3->kbar) vs stdout P", lag,
                           float(np.abs(a[:m] - b[:m]).max())])

wcsv(os.path.join(ANA, "temperature_source_crosscheck.csv"),
     ["run", "comparison", "lag_steps", "max_abs_diff_or_value_K"], tcross)
wcsv(os.path.join(ANA, "pressure_source_crosscheck.csv"),
     ["run", "comparison", "lag_steps", "max_abs_diff_kbar"], pcross)

# ---------- statistics ----------
def classify_T(T, t_ps):
    r100 = rolling(T, 100)
    r = r100[np.isfinite(r100)]
    if len(T) < 50:
        return "INCOMPLETE"
    if not np.all(np.isfinite(T)):
        return "NUMERICAL_FAILURE"
    if T.max() > 2000 and len(r) and r[-1] > 1000:
        return "PHYSICAL_HEATING"
    if len(r) and r[-1] < 30:
        return "PHYSICAL_COOLING"
    bl = np.array(blocks(T))
    # "near 300 K" requires the long-time block means to sit within +/-25% of tempw.
    if np.all(bl > 0.75 * 300.0) and np.all(bl < 1.25 * 300.0):
        return "STABLE_NEAR_300K"
    return "BOUNDED_BUT_OFFSET"


def classify_cell(vol, conds, dmins):
    if len(vol) < 10:
        return "INCOMPLETE"
    rel = abs(vol[-1] - vol[0]) / vol[0]
    if vol.min() <= 0 or max(conds) > 5 or min(dmins) < 1.5:
        return "SINGULAR_OR_UNSTABLE_CELL"
    if rel > 0.25:
        return "DRIFTING_PR_CELL"
    return "STABLE_PR_CELL"


tstat_rows, pstat_rows, sstat_rows, cstat_rows, estat_rows, blk_rows = [], [], [], [], [], []
SUM = {}
for tag in RUNS:
    d = D[tag]
    T, t_ps, P = d["T21"], d["t_ps"], d["p_kbar"]
    st, sp = stats(T), stats(P)
    tcls = classify_T(T, t_ps)

    # cell metrics
    vol = np.array([abs(np.linalg.det(m)) for m in d["cmats"]])
    conds, svs = [], []
    for m in d["cmats"]:
        c, s = cond_and_sv(m)
        conds.append(c); svs.append(s)
    conds = np.array(conds)
    H0 = d["cmats"][0]
    strains = np.array([strain6(m, H0) for m in d["cmats"]])
    ct_ps = (d["csteps"] - d["csteps"][0]) * DT_FS / 1000.0
    dmins, dmaxs, nns = [], [], []
    idx_sample = range(0, len(d["cmats"]), max(1, len(d["cmats"]) // 200))
    for i in idx_sample:
        # .cel columns are the lattice vectors; min_max_cc wants them as rows.
        Hrows = d["cmats"][i].T
        if i < len(d["ppos"]):
            # .pos holds Cartesian bohr in the same frame -> fractional
            fr = d["ppos"][i] @ np.linalg.inv(Hrows)
        else:
            fr = frac0
        a, b, c = min_max_cc(Hrows, fr)
        dmins.append(a); dmaxs.append(b); nns.append(c)
    ccls = classify_cell(vol, conds, dmins)

    # energy conservation: edyn is the extended PR total (ekint + ut), press=0
    ed = d["edyn"]
    e_drift_total = float(ed[-1] - ed[0])
    e_drift_ps = e_drift_total / max(t_ps[-1], 1e-12)
    e_drift_ev_atom_ps = e_drift_ps * 13.605693122994 / NAT
    span = float(np.abs(ed - ed[0]).max())
    if abs(e_drift_ev_atom_ps) < 0.01:
        ecls = "GOOD_CONSERVATION"
    elif abs(e_drift_ev_atom_ps) < 0.05:
        ecls = "MODERATE_DRIFT"
    else:
        ecls = "SEVERE_DRIFT"

    SUM[tag] = {
        "label": d["label"], "version": d["version"],
        "completion": completeness[tag]["classification"],
        "job_done": d["job_done"], "steps": int(len(T)),
        "temperature_instant_3Nm3": st,
        "temperature_classification": tcls,
        "temperature_drift_K_per_ps": drift(t_ps, T),
        "temperature_initialized_K": float(T[0]),
        "first_above_600K": first_above(t_ps, T, 600),
        "first_above_1000K": first_above(t_ps, T, 1000),
        "first_above_2000K": first_above(t_ps, T, 2000),
        "fraction_above_600K": frac_above(T, 600),
        "fraction_above_1000K": frac_above(T, 1000),
        "fraction_above_2000K": frac_above(T, 2000),
        "longest_run_above_600K": longest_run(T, 600),
        "longest_run_above_1000K": longest_run(T, 1000),
        "longest_run_above_2000K": longest_run(T, 2000),
        "pressure_kbar": sp,
        "pressure_GPa": {k: (v / 10.0 if isinstance(v, float) and k != "n" else v)
                         for k, v in sp.items()},
        "pressure_drift_GPa_per_ps": (drift(t_ps, P / 10.0)),
        "stress_mean_kbar": {n_: float(np.nanmean(d["Scomp"][:, i, j]))
                             for n_, (i, j) in zip(("sxx", "syy", "szz", "syz", "sxz", "sxy"),
                                                   ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)))}
        if len(d["Scomp"]) else {},
        "cell": {"volume_initial_bohr3": float(vol[0]), "volume_final_bohr3": float(vol[-1]),
                 "volume_mean_bohr3": float(vol.mean()),
                 "volume_change_percent": float(100 * (vol[-1] - vol[0]) / vol[0]),
                 "determinant_min_bohr3": float(vol.min()),
                 "condition_number_max": float(conds.max()),
                 "singular_values_final": [float(x) for x in svs[-1]],
                 "strain6_final": [float(x) for x in strains[-1]],
                 "min_CC_bohr": float(min(dmins)), "max_CC_bohr": float(max(dmaxs)),
                 "mean_NN_CC_bohr": float(np.mean(nns))},
        "cell_classification": ccls,
        "energy": {"quantity": "edyn from .e col3 = ekint + ut (extended PR total, press=0)",
                   "initial_Ry": float(ed[0]), "final_Ry": float(ed[-1]),
                   "drift_total_Ry": e_drift_total,
                   "drift_Ry_per_ps": e_drift_ps,
                   "drift_eV_per_atom_per_ps": e_drift_ev_atom_ps,
                   "max_excursion_Ry": span,
                   "classification": ecls},
    }

    tstat_rows.append([tag, d["label"], tcls, st["n"], st["initial"], st["mean"], st["std"],
                       st["median"], st["min"], st["max"], st["final"], st["first100_mean"],
                       st["final100_mean"], st["final200_mean"], SUM[tag]["temperature_drift_K_per_ps"],
                       SUM[tag]["fraction_above_600K"], SUM[tag]["longest_run_above_600K"]])
    pstat_rows.append([tag, d["label"], sp["n"], sp["initial"], sp["mean"], sp["std"], sp["median"],
                       sp["min"], sp["max"], sp["final"], sp["first100_mean"], sp["final100_mean"],
                       sp["final200_mean"], sp["mean"] / 10.0, sp["final200_mean"] / 10.0,
                       SUM[tag]["pressure_drift_GPa_per_ps"]])
    sm = SUM[tag]["stress_mean_kbar"]
    if sm:
        norm = [sm["sxx"], sm["syy"], sm["szz"]]
        sstat_rows.append([tag, sm["sxx"], sm["syy"], sm["szz"], sm["syz"], sm["sxz"], sm["sxy"],
                           float(max(norm) - min(norm)),
                           float(np.mean([abs(sm["syz"]), abs(sm["sxz"]), abs(sm["sxy"])]))])
    cc = SUM[tag]["cell"]
    cstat_rows.append([tag, ccls, cc["volume_initial_bohr3"], cc["volume_final_bohr3"],
                       cc["volume_change_percent"], cc["determinant_min_bohr3"],
                       cc["condition_number_max"], cc["min_CC_bohr"], cc["max_CC_bohr"],
                       cc["mean_NN_CC_bohr"]] + cc["strain6_final"])
    en = SUM[tag]["energy"]
    estat_rows.append([tag, en["quantity"], en["initial_Ry"], en["final_Ry"], en["drift_total_Ry"],
                       en["drift_Ry_per_ps"], en["drift_eV_per_atom_per_ps"],
                       en["max_excursion_Ry"], en["classification"]])
    for i, (bt, bp) in enumerate(zip(blocks(T), blocks(P))):
        blk_rows.append([tag, i + 1, i * 100, min((i + 1) * 100 - 1, len(T) - 1), bt, bp])

wcsv(os.path.join(ANA, "temperature_statistics.csv"),
     ["run", "label", "classification", "n", "initial_K", "mean_K", "std_K", "median_K", "min_K",
      "max_K", "final_K", "first100_mean_K", "final100_mean_K", "final200_mean_K",
      "drift_K_per_ps", "fraction_above_600K", "longest_run_above_600K"], tstat_rows)
wcsv(os.path.join(ANA, "pressure_statistics.csv"),
     ["run", "label", "n", "initial_kbar", "mean_kbar", "std_kbar", "median_kbar", "min_kbar",
      "max_kbar", "final_kbar", "first100_mean_kbar", "final100_mean_kbar", "final200_mean_kbar",
      "mean_GPa", "final200_mean_GPa", "drift_GPa_per_ps"], pstat_rows)
wcsv(os.path.join(ANA, "stress_statistics.csv"),
     ["run", "mean_sxx_kbar", "mean_syy_kbar", "mean_szz_kbar", "mean_syz_kbar", "mean_sxz_kbar",
      "mean_sxy_kbar", "normal_anisotropy_kbar", "mean_abs_shear_kbar"], sstat_rows)
wcsv(os.path.join(ANA, "cell_statistics.csv"),
     ["run", "classification", "volume_initial_bohr3", "volume_final_bohr3", "volume_change_percent",
      "determinant_min_bohr3", "condition_number_max", "min_CC_bohr", "max_CC_bohr", "mean_NN_CC_bohr",
      "strain_xx", "strain_yy", "strain_zz", "strain_yz", "strain_xz", "strain_xy"], cstat_rows)
wcsv(os.path.join(ANA, "energy_drift_statistics.csv"),
     ["run", "quantity", "initial_Ry", "final_Ry", "drift_total_Ry", "drift_Ry_per_ps",
      "drift_eV_per_atom_per_ps", "max_excursion_Ry", "classification"], estat_rows)
wcsv(os.path.join(ANA, "block_statistics.csv"),
     ["run", "block", "nst_start", "nst_end", "temperature_mean_K", "pressure_mean_kbar"], blk_rows)

# ---------- Phase 13: comparison ----------
a, b = D["qe74"], D["qe75"]
n = min(len(a["T21"]), len(b["T21"]))
dT = np.abs(a["T21"][:n] - b["T21"][:n])
dP = np.abs(a["p_kbar"][:n] - b["p_kbar"][:n])
nzT = np.nonzero(dT > 0)[0]
nzP = np.nonzero(dP > 0)[0]
identical_init = bool(dT[0] == 0.0 and a["eka"][0] == b["eka"][0])
cmp_ = {
    "identical_initial_kinetic_energy": identical_init,
    "initial_eka_Ry": {"qe74": float(a["eka"][0]), "qe75": float(b["eka"][0])},
    "initial_T_3Nm3_K": {"qe74": float(a["T21"][0]), "qe75": float(b["T21"][0])},
    "first_differing_step_temperature": int(nzT[0]) if len(nzT) else None,
    "first_differing_step_pressure": int(nzP[0]) if len(nzP) else None,
    "identical_prefix_steps_temperature": int(nzT[0]) if len(nzT) else n,
    "max_abs_T_difference_K": float(dT.max()),
    "max_abs_P_difference_kbar": float(dP.max()),
    "mean_T_delta_K": abs(SUM["qe74"]["temperature_instant_3Nm3"]["mean"]
                          - SUM["qe75"]["temperature_instant_3Nm3"]["mean"]),
    "final200_T_delta_K": abs(SUM["qe74"]["temperature_instant_3Nm3"]["final200_mean"]
                              - SUM["qe75"]["temperature_instant_3Nm3"]["final200_mean"]),
    "mean_P_delta_kbar": abs(SUM["qe74"]["pressure_kbar"]["mean"] - SUM["qe75"]["pressure_kbar"]["mean"]),
    "T_fluctuation_scale_K": float(np.sqrt((SUM["qe74"]["temperature_instant_3Nm3"]["std"] ** 2
                                            + SUM["qe75"]["temperature_instant_3Nm3"]["std"] ** 2) / 2)),
    "P_fluctuation_scale_kbar": float(np.sqrt((SUM["qe74"]["pressure_kbar"]["std"] ** 2
                                               + SUM["qe75"]["pressure_kbar"]["std"] ** 2) / 2)),
    "same_temperature_class": SUM["qe74"]["temperature_classification"] == SUM["qe75"]["temperature_classification"],
    "same_cell_class": SUM["qe74"]["cell_classification"] == SUM["qe75"]["cell_classification"],
    "volume_change_delta_percent": abs(SUM["qe74"]["cell"]["volume_change_percent"]
                                       - SUM["qe75"]["cell"]["volume_change_percent"]),
    "tv_step0_reporting_difference": {
        "qe74_tv_nst0_K": float(a["tv_T"][0]), "qe75_tv_nst0_K": float(b["tv_T"][0]),
        "explanation": ("af7ebfc adds a tnew initialisation after vcinit in PW/src/vcsmd.f90 "
                        "(absent in pristine a11913a), so QE 7.5 reports 2*ekint/(3(N+1)k_B) "
                        "at nst=0 while QE 7.4.1 reports the uninitialised 0.0. Diagnostic only: "
                        "the dynamics are bit-identical for the first ~700 steps."),
    },
}
temp_similar = cmp_["same_temperature_class"] and cmp_["mean_T_delta_K"] <= max(2 * cmp_["T_fluctuation_scale_K"], 50)
press_similar = cmp_["mean_P_delta_kbar"] <= max(2 * cmp_["P_fluctuation_scale_kbar"], 20)
cmp_["temperature_broadly_similar"] = bool(temp_similar)
cmp_["pressure_broadly_similar"] = bool(press_similar)

ca, cb = SUM["qe74"]["temperature_classification"], SUM["qe75"]["temperature_classification"]
if completeness["qe74"]["classification"] != "COMPLETE" or completeness["qe75"]["classification"] != "COMPLETE":
    overall = "PARTIAL"
elif ca == "PHYSICAL_HEATING" and cb == "PHYSICAL_HEATING":
    overall = "BOTH_HEAT"
elif ca == "STABLE_NEAR_300K" and cb == "PHYSICAL_HEATING":
    overall = "QE74_STABLE_QE75_HEATS"
elif ca == "PHYSICAL_HEATING" and cb == "STABLE_NEAR_300K":
    overall = "QE74_HEATS_QE75_STABLE"
elif temp_similar and press_similar:
    overall = "BOTH_COMPARABLE"
elif temp_similar and not press_similar:
    overall = "BOTH_STABLE_DIFFERENT_PRESSURE"
elif not temp_similar:
    overall = "BOTH_BOUNDED_DIFFERENT_TEMPERATURE"
else:
    overall = "PARTIAL"

summary = {
    "job_id": "50160152",
    "slurm_states": {"top": "COMPLETED 0:0", "batch": "COMPLETED 0:0", "extern": "COMPLETED 0:0",
                     "50160152.0 (pw.x qe74)": "COMPLETED 0:0",
                     "50160152.1 (pw.x qe75)": "COMPLETED 0:0"},
    "remote_workdir": "/leonardo/home/userexternal/mcheng00/workspace/Diamond/"
                      "diamond8_no_velocity_qe74_vs_qe75_caimd_false_300K_20260724T092828Z",
    "overall_status": overall,
    "ensemble_label": ("ordinary Parrinello-Rahman variable-cell dynamics after native 300 K "
                       "velocity initialization; NOT NVT - no continuing thermostat is active"),
    "temperature_mode": {"keyword": "ion_temperature='initial'", "tempw_K": 300.0,
                         "continuing_thermostat_active": False},
    "input_contract": {"physical_settings_identical": ident["physical_settings_identical"],
                       "explicit_velocities_absent": no_vel["explicit_velocities_absent_in_both"],
                       "qe75_only_caimd_false": caimd_audit["qe75_only_caimd_false"],
                       "qe74_no_caimd_keyword": caimd_audit["qe74_has_no_caimd_keyword"]},
    "runs": SUM,
    "comparison": cmp_,
}
json.dump(summary, open(os.path.join(ANA, "comparison_summary_audited.json"), "w"), indent=2)
print(json.dumps({"overall": overall,
                  "qe74": [ca, SUM["qe74"]["cell_classification"], SUM["qe74"]["energy"]["classification"]],
                  "qe75": [cb, SUM["qe75"]["cell_classification"], SUM["qe75"]["energy"]["classification"]],
                  "identical_prefix": cmp_["identical_prefix_steps_temperature"],
                  "Tmean": [SUM["qe74"]["temperature_instant_3Nm3"]["mean"],
                            SUM["qe75"]["temperature_instant_3Nm3"]["mean"]],
                  "Pmean_kbar": [SUM["qe74"]["pressure_kbar"]["mean"],
                                 SUM["qe75"]["pressure_kbar"]["mean"]]}, indent=2))

# save arrays for the plotting stage
np.savez(os.path.join(EVD, "_audit_arrays.npz"),
         **{f"{tag}_{k}": np.asarray(D[tag][k], dtype=float)
            for tag in RUNS for k in ("t_ps", "T21", "p_kbar", "vcell", "edyn", "eka", "ekla")},
         qe74_vol=np.array([abs(np.linalg.det(m)) for m in D["qe74"]["cmats"]]),
         qe75_vol=np.array([abs(np.linalg.det(m)) for m in D["qe75"]["cmats"]]),
         qe74_cond=np.array([cond_and_sv(m)[0] for m in D["qe74"]["cmats"]]),
         qe75_cond=np.array([cond_and_sv(m)[0] for m in D["qe75"]["cmats"]]),
         qe74_S=D["qe74"]["Scomp"], qe75_S=D["qe75"]["Scomp"],
         qe74_cel_t=(D["qe74"]["csteps"] - D["qe74"]["csteps"][0]) * DT_FS / 1000.0,
         qe75_cel_t=(D["qe75"]["csteps"] - D["qe75"]["csteps"][0]) * DT_FS / 1000.0)
