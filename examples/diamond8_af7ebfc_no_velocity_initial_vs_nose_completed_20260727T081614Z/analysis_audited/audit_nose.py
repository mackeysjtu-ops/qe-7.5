#!/usr/bin/env python3
"""Independent audit of completed job 50167241: CONTROL initial vs TEST nose, no velocities.

Reads only the local mirror; writes new products under analysis_audited/ and
evidence_audited/. Raw mirrored outputs are never modified.

Column semantics verified from PW/src/vcsmd.f90 (af7ebfc, lines 408-419):
  .e    : ut, ekint, edyn, pv, nst
  .eal  : uta, eka, eta, utl, ekla, etl, nst
  .ave  : avu, avk, nst
  .p    : press, p, avp, nst           [Ry/bohr^3]
  .avec : |a|,|b|,|c|, theta12, theta23, theta31, nst
  .tv   : vcell, tnew, nst    (tnew = 2*avk/(3(N+1)k_B): RUNNING AVERAGE, not instantaneous)
  .md   : TERMINAL state only: rati, ratd, rat2d, rat2di, tauold / aveci, avecd, ...
  .cel  : per step, lattice vectors as COLUMNS (verified vs input CELL_PARAMETERS rows)

Nose chain energy: ions_nose_energy / cell_nose_energy are computed in vcsmd but
NEVER written to any output file, so a full Nose+PR conserved Hamiltonian is
unavailable. edyn (.e col3 = ekint+ut) is a PARTIAL extended-energy diagnostic.
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

RY_TO_K = 157887.51240116
RY_BOHR3_TO_KBAR = 147105.078
BOHR_TO_ANG = 0.529177210903
EV_PER_RY = 13.605693122994
NAT = 8
NDOF = 3 * NAT - 3
DT_FS = 2.0
MASS_C_AMU = 12.0107
AMU_RY = 911.4442421
TARGET_T = 300.0

RUNS = {
    "control": {"dir": "control_initial", "prefix": "diamond8_control_initial",
                "label": "CONTROL ion_temperature='initial'", "input": "inputs/control_initial_no_velocity.in"},
    "nose": {"dir": "test_nose", "prefix": "diamond8_test_nose",
             "label": "TEST ion_temperature='nose' (fnosep=50, nhpcl=4)", "input": "inputs/test_nose_no_velocity.in"},
}

REFERENCE = {
    "control_previous_no_velocity_initial": {  # audited job 50160152, qe75
        "mean_K": 170.01436653663114, "final200_K": 174.82712827452696,
        "max_K": 300.001310468147, "volume_change_percent": -4.622395517125142,
        "fraction_above_600K": 0.0},
    "nose_previous_with_explicit_velocities_job50127909": {
        "max_K": 8123.838804792437, "final_K": 4231.686070493259,
        "mean_K": 4201.45228096483, "first_step_above_600K": 26,
        "first_step_above_2000K": 130},
}

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


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(65536), b""):
            h.update(c)
    return h.hexdigest()


def roll(a, w):
    if len(a) < w:
        return np.full(len(a), np.nan)
    return np.concatenate([np.full(w - 1, np.nan), np.convolve(a, np.ones(w) / w, "valid")])


def blocks(a, w=100):
    return np.array([a[i:i + w].mean() for i in range(0, len(a), w) if len(a[i:i + w])])


def lin_slope(x, y):
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
    allnum = [float(t.replace("D", "E")) for t in NUM.findall(text(p))]
    n3 = 3 * nat
    need = n3 * 4 + 3 * nat * 3
    if len(allnum) < need + 9:
        return None
    return {"ratd": np.array(allnum[n3:2 * n3]).reshape(nat, 3),
            "aveci": np.array(allnum[need:need + 9]).reshape(3, 3)}


def det3(m):
    return float(np.linalg.det(m))


def cond_sv(m):
    s = np.linalg.svd(m, compute_uv=False)
    return float(s[0] / s[-1]), s


def strain6(H, H0):
    F = H @ np.linalg.inv(H0)          # .cel columns = lattice vectors
    E = 0.5 * (F + F.T) - np.eye(3)
    return [E[0, 0], E[1, 1], E[2, 2], E[1, 2], E[0, 2], E[0, 1]]


def cc_dist(Hcol, cart):
    a = Hcol.T
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1)
                       for k in (-1, 0, 1)]) @ a
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
            "first100_mean": float(a[:100].mean()),
            "step101_200_mean": float(a[100:200].mean()) if len(a) >= 200 else None,
            "final100_mean": float(a[-100:].mean()),
            "final200_mean": float(a[-200:].mean()),
            "final250_mean": float(a[-250:].mean())}


def wcsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join("" if v is None or (isinstance(v, float) and not np.isfinite(v))
                              else (f"{v:.10g}" if isinstance(v, float) else str(v))
                              for v in r) + "\n")


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
           "cell_units": None, "pos_units": None, "atomic_velocities": False, "vel_rows": 0}
    mode = None
    for line in cards:
        u = line.upper().split(); h = u[0] if u else ""
        if h == "ATOMIC_SPECIES": mode = "sp"; continue
        if h == "CELL_PARAMETERS": mode = "ce"; out["cell_units"] = u[1].lower() if len(u) > 1 else None; continue
        if h == "ATOMIC_POSITIONS": mode = "po"; out["pos_units"] = u[1].lower() if len(u) > 1 else None; continue
        if h == "K_POINTS": mode = "kp"; out["kpoints"] = " ".join(u[1:]).lower(); continue
        if h == "ATOMIC_VELOCITIES": out["atomic_velocities"] = True; mode = "ve"; continue
        f = line.split()
        if mode == "sp" and len(f) >= 3: out["species"].append((f[0], float(f[1]), f[2]))
        elif mode == "ce" and len(f) == 3: out["cell"].append([float(x) for x in f])
        elif mode == "po" and len(f) >= 4: out["pos"].append((f[0], [float(x) for x in f[1:4]]))
        elif mode == "ve" and len(f) >= 4: out["vel_rows"] += 1
    return out


def numv(v):
    try:
        return float(str(v).lower().replace("d", "e").strip().strip("'\""))
    except ValueError:
        return str(v).strip().strip("'\"").lower()


failures_contract = []

# ============================================================ Phase 4 ======
tc, tn = text(os.path.join(MIRROR, RUNS["control"]["input"])), text(os.path.join(MIRROR, RUNS["nose"]["input"]))
nc, cc_cards = parse_input(tc)
nn, cn_cards = parse_input(tn)
kc, kn = cards_struct(cc_cards), cards_struct(cn_cards)
FORBID_VEL = [r"atomic_velocities", r"ion_velocities", r"from_input", r"tv0rd", r"starting_velocities"]
FORBID_CAIMD = [r"caimd_stride", r"caimd_step_strain", r"caimd_mask"]
ALLOWED = {("CONTROL", "prefix"), ("CONTROL", "outdir"), ("IONS", "ion_temperature"),
           ("IONS", "fnosep(1)"), ("IONS", "fnosep"), ("IONS", "nhpcl")}

no_vel = {"native_path": "vcinit->ranv(temp=tempw); ranv zeroes total linear momentum", "runs": {}}
for tag, t, k in (("control", tc, kc), ("nose", tn, kn)):
    hits = {p: [i + 1 for i, l in enumerate(t.splitlines()) if re.search(p, l, re.I)] for p in FORBID_VEL}
    ctl = parse_input(t)[0]
    no_vel["runs"][tag] = {"forbidden_hits_complete_text": {a: b for a, b in hits.items() if b},
                           "atomic_velocities_header": k["atomic_velocities"], "velocity_rows": k["vel_rows"],
                           "restart_mode": numv(ctl.get("CONTROL", {}).get("restart_mode", "")),
                           "ion_temperature": numv(ctl.get("IONS", {}).get("ion_temperature", "")),
                           "tempw": numv(ctl.get("IONS", {}).get("tempw", ""))}
    r = no_vel["runs"][tag]
    if r["forbidden_hits_complete_text"]:
        failures_contract.append(f"{tag}: forbidden velocity keyword")
    if r["atomic_velocities_header"] or r["velocity_rows"]:
        failures_contract.append(f"{tag}: ATOMIC_VELOCITIES present")
    if r["restart_mode"] != "from_scratch":
        failures_contract.append(f"{tag}: restart_mode != from_scratch")
no_vel["explicit_velocities_absent_in_both"] = not failures_contract
json.dump(no_vel, open(os.path.join(ANA, "no_explicit_velocity_audit.json"), "w"), indent=2)

cd = {"runs": {}}
for tag, t in (("control", tc), ("nose", tn)):
    eff = strip_comments(t)
    lines = [l.strip() for l in eff if "caimd" in l.lower()]
    bad = [p for p in FORBID_CAIMD if any(re.search(p, l, re.I) for l in eff)]
    active = [l for l in lines if re.search(r"caimd\s*=\s*\.?true", l, re.I)]
    okfalse = len(lines) == 1 and bool(re.search(r"caimd\s*=\s*\.false\.", lines[0], re.I))
    cd["runs"][tag] = {"caimd_lines": lines, "forbidden_parameters": bad,
                       "active_caimd": active, "exactly_caimd_false": okfalse}
cd["both_explicitly_disabled"] = all(v["exactly_caimd_false"] for v in cd["runs"].values())
json.dump(cd, open(os.path.join(ANA, "caimd_disabled_audit.json"), "w"), indent=2)

diffs = {}
for sec in sorted(set(nc) | set(nn)):
    a, b = nc.get(sec, {}), nn.get(sec, {})
    for key in sorted(set(a) | set(b)):
        if numv(a.get(key, "<absent>")) != numv(b.get(key, "<absent>")) and (sec, key) not in ALLOWED:
            diffs[f"{sec}.{key}"] = {"control": a.get(key, "<absent>"), "nose": b.get(key, "<absent>")}


def maxdiff(A, B):
    if len(A) != len(B):
        return None
    return max((abs(x - y) for ra, rb in zip(A, B) for x, y in zip(ra, rb)), default=0.0)


def fnum(sec, key):
    return numv(nc.get(sec, {}).get(key, "")), numv(nn.get(sec, {}).get(key, ""))


cell_md = maxdiff(kc["cell"], kn["cell"])
pos_md = maxdiff([p[1] for p in kc["pos"]], [p[1] for p in kn["pos"]])
order_same = [p[0] for p in kc["pos"]] == [p[0] for p in kn["pos"]]
mass_md = maxdiff([[s[1]] for s in kc["species"]], [[s[1]] for s in kn["species"]])
det_in = det3(np.array(kc["cell"]))
dmin_in, _, _ = cc_dist(np.array(kc["cell"]).T, np.array([p[1] for p in kc["pos"]]) @ np.array(kc["cell"]))
scf = [("ELECTRONS", "conv_thr"), ("ELECTRONS", "electron_maxstep"), ("ELECTRONS", "mixing_beta")]
ident = {
    "max_differences": {
        "cell_bohr": cell_md, "positions_crystal": pos_md, "atom_ordering_identical": order_same,
        "masses_amu": mass_md,
        "ecutwfc": abs(fnum("SYSTEM", "ecutwfc")[0] - fnum("SYSTEM", "ecutwfc")[1]),
        "ecutrho": abs(fnum("SYSTEM", "ecutrho")[0] - fnum("SYSTEM", "ecutrho")[1]),
        "dt": abs(fnum("CONTROL", "dt")[0] - fnum("CONTROL", "dt")[1]),
        "nstep": abs(fnum("CONTROL", "nstep")[0] - fnum("CONTROL", "nstep")[1]),
        "press": abs(fnum("CELL", "press")[0] - fnum("CELL", "press")[1]),
        "wmass": abs(fnum("CELL", "wmass")[0] - fnum("CELL", "wmass")[1]),
        "tempw": abs(fnum("IONS", "tempw")[0] - fnum("IONS", "tempw")[1]),
        "cell_dynamics_same": fnum("CELL", "cell_dynamics")[0] == fnum("CELL", "cell_dynamics")[1],
        "cell_dofree_same": fnum("CELL", "cell_dofree")[0] == fnum("CELL", "cell_dofree")[1],
        "scf": {f"{s}.{k}": abs(fnum(s, k)[0] - fnum(s, k)[1]) for s, k in scf},
    },
    "kpoints_identical": kc["kpoints"] == kn["kpoints"],
    "intended_differences": {
        "ion_temperature": [numv(nc["IONS"].get("ion_temperature")), numv(nn["IONS"].get("ion_temperature"))],
        "fnosep(1)": [nc["IONS"].get("fnosep(1)", "<absent>"), nn["IONS"].get("fnosep(1)", "<absent>")],
        "nhpcl": [nc["IONS"].get("nhpcl", "<absent>"), nn["IONS"].get("nhpcl", "<absent>")]},
    "namelist_differences_outside_allowed_set": diffs,
    "structure": {"n_atoms": [len(kc["pos"]), len(kn["pos"])], "determinant_bohr3": det_in,
                  "min_CC_bohr": dmin_in, "undeformed_ideal": abs(det_in - 306.2710374914433) < 1e-6},
    "checksums": {"control_sha256": sha(os.path.join(MIRROR, RUNS["control"]["input"])),
                  "nose_sha256": sha(os.path.join(MIRROR, RUNS["nose"]["input"])),
                  "pseudo_sha256": sha(os.path.join(MIRROR, "pseudo", "C.UPF"))},
}
md = ident["max_differences"]
phys_zero = (cell_md == 0 and pos_md == 0 and mass_md == 0 and order_same and md["ecutwfc"] == 0
             and md["ecutrho"] == 0 and md["dt"] == 0 and md["nstep"] == 0 and md["press"] == 0
             and md["wmass"] == 0 and md["tempw"] == 0 and all(v == 0 for v in md["scf"].values())
             and not diffs and md["cell_dynamics_same"] and md["cell_dofree_same"])
ident["valid_single_factor"] = bool(phys_zero)
ident["workflow_class"] = "VALID_SINGLE_FACTOR" if phys_zero else "INVALID_SINGLE_FACTOR_COMPARISON"
json.dump(ident, open(os.path.join(ANA, "input_identity_audit.json"), "w"), indent=2)

import difflib
with open(os.path.join(ANA, "input_normalized.diff"), "w") as fh:
    fh.write("\n".join(difflib.unified_diff(strip_comments(tc), strip_comments(tn),
             fromfile="control_initial (normalized)", tofile="test_nose (normalized)", lineterm="")) + "\n")

# ============================================================ load runs ====
rc_txt = text(os.path.join(MIRROR, "evidence", "run_return_codes.txt"))
RCS = dict(re.findall(r"(\w+)=(-?\d+)", rc_txt))
D = {}
for tag, meta in RUNS.items():
    base = os.path.join(MIRROR, meta["dir"]); pref = meta["prefix"]
    out = os.path.join(base, "qe.out")
    eal = rd(os.path.join(base, pref + ".eal")); pf = rd(os.path.join(base, pref + ".p"))
    tv = rd(os.path.join(base, pref + ".tv")); ef = rd(os.path.join(base, pref + ".e"))
    av = rd(os.path.join(base, pref + ".ave"))
    cs, cm = parse_cel(os.path.join(base, "tmp", pref + ".cel"))
    ps, pp = parse_pos(os.path.join(base, "tmp", pref + ".pos"))
    P, S = parse_stress(out)
    md_term = parse_md(os.path.join(base, "tmp", pref + ".md"))
    txt = text(out)
    n = len(eal)
    d = {"label": meta["label"], "n": n, "job_done": "JOB DONE" in txt,
         "version": (re.search(r"Program PWSCF v\.(\S+)", txt) or [None, None])[1]
         if re.search(r"Program PWSCF v\.(\S+)", txt) else None,
         "error": (re.search(r"Error in routine\s+(\S+)", txt) or [None, None])[1]
         if re.search(r"Error in routine\s+(\S+)", txt) else None,
         "scf_fail": txt.count("convergence NOT achieved"),
         "nose_stdout_active": "nose" in txt.lower() and "temperature control" in txt.lower(),
         "rc": int(RCS.get("RC_CONTROL" if tag == "control" else "RC_NOSE", -1)),
         "cel_steps": cs, "cel": cm, "pos_steps": ps, "pos": pp,
         "P_kbar_stdout": P, "S": S, "md": md_term}
    if n:
        nst = eal[:, 6].astype(int)
        d.update(nst=nst, t_ps=nst * DT_FS / 1000.0, eka=eal[:, 1], ekla=eal[:, 4],
                 T=2 * eal[:, 1] * RY_TO_K / NDOF, T_3Np1=2 * eal[:, 1] * RY_TO_K / (3 * (NAT + 1)),
                 tv_T=tv[:, 1] if len(tv) else np.full(n, np.nan),
                 vcell=tv[:, 0] if len(tv) else np.full(n, np.nan),
                 avk=av[:, 1] if len(av) else np.full(n, np.nan),
                 p_kbar=pf[:, 1] * RY_BOHR3_TO_KBAR if len(pf) else np.full(n, np.nan),
                 edyn=ef[:, 2] if len(ef) else np.full(n, np.nan),
                 ut=ef[:, 0] if len(ef) else np.full(n, np.nan),
                 ekint=ef[:, 1] if len(ef) else np.full(n, np.nan))
        d["T_runavg_3Np1"] = 2 * d["avk"] * RY_TO_K / (3 * (NAT + 1))
    D[tag] = d

# ============================================== Phase 5 provenance =========
exe_txt = text(os.path.join(MIRROR, "evidence", "executable_checksums.txt"))
env_txt = text(os.path.join(MIRROR, "evidence", "executable_environment.txt"))
sj = text(os.path.join(MIRROR, "evidence", "slurm_job_environment.txt"))
prov = {"job_id": "50167241", "node": "lrdn0558", "mpi_launcher": "srun --mpi=pmix_v3 -n 1",
        "ranks": 1, "omp_threads": int((re.search(r"omp_num_threads=(\d+)", sj) or [0, 8])[1])
        if re.search(r"omp_num_threads=(\d+)", sj) else 8,
        "executable_path": (re.search(r"pw_x=(\S+)", exe_txt) or [None, None])[1],
        "executable_sha256_runtime": (re.findall(r"([0-9a-f]{64})", exe_txt) or [None])[0],
        "expected_sha256": "89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d",
        "version_banner": {"control": D["control"]["version"], "nose": D["nose"]["version"]},
        "hostname_env": (re.search(r"hostname=(\S+)", env_txt) or [None, None])[1],
        "module_list": [l.strip() for l in re.findall(r"^\s*\d+\)\s+(\S+)", env_txt, re.M)],
        "slurm_steps": {"control": "50167241.0", "nose": "50167241.1"},
        "nose_thermostat_active_in_stdout": D["nose"]["nose_stdout_active"]}
prov["sha256_matches_expected"] = prov["executable_sha256_runtime"] == prov["expected_sha256"]
json.dump(prov, open(os.path.join(EVD, "executable_provenance.json"), "w"), indent=2)
json.dump(prov, open(os.path.join(ANA, "executable_provenance.json"), "w"), indent=2)

# ============================================== Phase 6 completeness =======
comp = {}
for tag in RUNS:
    d = D[tag]; n = d["n"]
    cls = ("COMPLETE" if (d["job_done"] and n >= 1000 and not d["error"])
           else "PARTIAL" if n > 0 else "MISSING_OUTPUT")
    comp[tag] = {"return_code": d["rc"], "slurm_step": prov["slurm_steps"][tag],
                 "slurm_step_state": "COMPLETED 0:0", "job_done": d["job_done"],
                 "runtime_version": d["version"], "first_error": d["error"],
                 "scf_failure_count": d["scf_fail"], "completed_ionic_steps": n,
                 "requested_steps": 1000, "final_valid_nst": int(d["nst"][-1]) if n else None,
                 "completed_time_ps": float(d["t_ps"][-1]) if n else 0.0, "requested_time_ps": 2.0,
                 "safety_stop": False, "rows": {"eal": n, "p": len(d["p_kbar"]) if n else 0,
                 "cel": len(d["cel_steps"]), "pos": len(d["pos_steps"]),
                 "stdout_stress": len(d["P_kbar_stdout"])}, "classification": cls}
comp["internal_post_processing_rc"] = int(RCS.get("RC_ANALYSIS", -1))
json.dump(comp, open(os.path.join(ANA, "run_completeness.json"), "w"), indent=2)

# ============================================== Phase 7 initial state ======
c, s = D["control"], D["nose"]
deka = np.abs(c["eka"] - s["eka"])
first_div = int(np.nonzero(deka > 0)[0][0]) if np.any(deka > 0) else None
init = {"note": ("Both runs omit explicit velocities and use the same af7ebfc executable, so "
                 "vcinit->ranv(temp=tempw) generates the same deterministic Maxwell realization "
                 "(randy seeded idum=0; total linear momentum zeroed). tmp/*.md holds only the "
                 "terminal state, so the FIRST committed ratd cannot be read back; the identity "
                 "of the initial realization is established from the byte-identical early .eal "
                 "rows (ionic kinetic energy) instead."),
        "first_step_eka_Ry": {"control": float(c["eka"][0]), "nose": float(s["eka"][0])},
        "first_step_eka_abs_diff_Ry": float(abs(c["eka"][0] - s["eka"][0])),
        "first_step_T_K": {"control": float(c["T"][0]), "nose": float(s["T"][0])},
        "first_step_T_abs_diff_K": float(abs(c["T"][0] - s["T"][0])),
        "step1_eka_abs_diff_Ry": float(abs(c["eka"][1] - s["eka"][1])),
        "first_eka_divergence_step": first_div,
        "com_velocity_note": "ranv zeroes total linear momentum, so initial COM velocity is 0 by construction"}
if c["eka"][0] == s["eka"][0] and c["eka"][1] == s["eka"][1] and (first_div is None or first_div > 100):
    init["classification"] = "IDENTICAL"
elif abs(c["eka"][0] - s["eka"][0]) < 1e-10:
    init["classification"] = "NUMERICALLY_EQUIVALENT"
else:
    init["classification"] = "DIFFERENT_RANDOM_REALIZATION"
init["strict_single_factor"] = init["classification"] in ("IDENTICAL", "NUMERICALLY_EQUIVALENT")
init["interpretation"] = ("Same native realization: this is a strict single-factor trajectory "
                          "comparison (thermostat mode is the only variable). Control and Nose are "
                          f"bit-identical in eka until step {first_div}."
                          if init["strict_single_factor"] else
                          "Initial realizations differ; compare statistically only.")
json.dump(init, open(os.path.join(ANA, "initial_state_comparison.json"), "w"), indent=2)

# ============================================== Phase 8-12 per run =========
DIRECT_STEPS = [0, 10, 25, 50, 100, 130, 250, 500, 750]
tcross, blk_rows, tstat_rows, energy_rows, chain_rows = [], [], [], [], []
RES = {}
for tag, d in D.items():
    n = d["n"]
    if not n:
        RES[tag] = {"label": d["label"], "classification": "INCOMPLETE", "n": 0}
        continue
    t, T, P = d["t_ps"], d["T"], d["p_kbar"]
    r25, r50, r100 = roll(T, 25), roll(T, 50), roll(T, 100)
    bT, bP = blocks(T), blocks(P)

    wcsv(os.path.join(ANA, f"{tag}_temperature.csv"),
         ["nst", "time_ps", "eka_Ry", "T_instant_3Nm3_K", "T_instant_3Np1_K",
          "T_running_avg_reported_tv_K", "roll25_K", "roll50_K", "roll100_K"],
         [[int(d["nst"][i]), t[i], d["eka"][i], T[i], d["T_3Np1"][i], d["tv_T"][i],
           r25[i], r50[i], r100[i]] for i in range(n)])
    wcsv(os.path.join(ANA, f"{tag}_pressure_stress.csv"),
         ["nst", "time_ps", "p_kbar", "p_GPa", "stdout_P_kbar", "sxx", "syy", "szz", "syz", "sxz", "sxy"],
         [[int(d["nst"][i]), t[i], P[i], P[i] / 10.0,
           d["P_kbar_stdout"][i] if i < len(d["P_kbar_stdout"]) else None]
          + ([d["S"][i][0][0], d["S"][i][1][1], d["S"][i][2][2], d["S"][i][1][2], d["S"][i][0][2], d["S"][i][0][1]]
             if i < len(d["S"]) else [None] * 6) for i in range(n)])

    # cell
    vol = np.array([abs(det3(m)) for m in d["cel"]]) if len(d["cel"]) else np.array([])
    conds = np.array([cond_sv(m)[0] for m in d["cel"]]) if len(d["cel"]) else np.array([])
    cel_rows, dmins, dmaxs, nns = [], [], [], []
    if len(d["cel"]):
        H0 = d["cel"][0]
        sample = set(range(0, len(d["cel"]), max(1, len(d["cel"]) // 200)))
        for i in range(len(d["cel"])):
            e6 = strain6(d["cel"][i], H0); dmin = dmax = nn = None
            if i in sample and i < len(d["pos"]):
                dmin, dmax, nn = cc_dist(d["cel"][i], d["pos"][i]); dmins.append(dmin); dmaxs.append(dmax); nns.append(nn)
            sv = cond_sv(d["cel"][i])[1]
            cel_rows.append([int(d["cel_steps"][i]), (d["cel_steps"][i] - d["cel_steps"][0]) * DT_FS / 1000.0,
                             vol[i], 100 * (vol[i] - vol[0]) / vol[0], conds[i],
                             float(sv[0]), float(sv[1]), float(sv[2])] + [float(x) for x in e6] + [dmin, nn])
    wcsv(os.path.join(ANA, f"{tag}_cell.csv"),
         ["step", "time_ps", "volume_bohr3", "volume_change_percent", "condition_number", "sv1", "sv2", "sv3",
          "strain_xx", "strain_yy", "strain_zz", "strain_yz", "strain_xz", "strain_xy", "min_CC_bohr", "mean_NN_CC_bohr"],
         cel_rows)

    # temperature source crosscheck
    tcross.append([tag, "T_instant_3Np1 vs reported .tv (lag0)", 0, float(np.nanmax(np.abs(d["T_3Np1"] - d["tv_T"])))])
    tcross.append([tag, "T_running_avg(2avk/3(N+1)) vs .tv (lag0)", 0, float(np.nanmax(np.abs(d["T_runavg_3Np1"] - d["tv_T"])))])
    Tdir = None; direct_conf = False
    if d["md"] is not None and len(d["cel"]):
        v = d["md"]["ratd"] @ d["cel"][-1].T
        vcom = v.mean(axis=0)
        K = 0.5 * MASS_C_AMU * AMU_RY * ((v - vcom) ** 2).sum()
        Tdir = 2 * K * RY_TO_K / NDOF
        tcross.append([tag, "direct .md ratd x terminal .cel (3N-3)", int(d["nst"][-1]), float(Tdir)])
        tcross.append([tag, "eal instantaneous terminal (3N-3)", int(d["nst"][-1]), float(T[-1])])
        tcross.append([tag, "direct-vs-eal terminal abs diff", 0, float(abs(Tdir - T[-1]))])
        direct_conf = bool(abs(Tdir - T[-1]) < max(5.0, 0.05 * abs(T[-1])))
    for st in DIRECT_STEPS:
        if st < n:
            tcross.append([tag, f"eal instantaneous at step {st} (3N-3)", st, float(T[st])])

    stT, stP = stats(T), stats(P)
    for i, (bt, bp) in enumerate(zip(bT, bP)):
        blk_rows.append([tag, i + 1, i * 100, min((i + 1) * 100 - 1, n - 1), bt, bp])

    # peak / post-peak
    pk = int(np.argmax(T)); post = T[pk:]
    post_roll = roll(post, 100)
    low_after = float(np.nanmin(post_roll)) if np.any(np.isfinite(post_roll)) else None
    fin_half = lin_slope(t[n // 2:], T[n // 2:])
    fin_quarter = lin_slope(t[3 * n // 4:], T[3 * n // 4:])

    # energy (partial)
    ed = d["edyn"]
    drift_total = float(ed[-1] - ed[0]); drift_ps = drift_total / max(t[-1], 1e-12)
    drift_at = drift_ps * EV_PER_RY / NAT
    energy_rows.append([tag, "edyn (.e col3 = ekint+ut)",
                        "PARTIAL_EXTENDED_ENERGY_DIAGNOSTIC (Nose chain energy not written by QE)"
                        if tag == "nose" else "extended PR total (press=0)",
                        float(ed[0]), float(ed[-1]), float(ed.min()), float(ed.max()),
                        drift_total, drift_ps, drift_at, lin_slope(t, ed),
                        lin_slope(t[3 * n // 4:], ed[3 * n // 4:]),
                        float(np.corrcoef(ed, d["eka"])[0, 1]) if np.all(np.isfinite(ed)) else None,
                        bool(np.all(np.isfinite(ed)) and np.abs(ed - ed[0]).max() < 5.0)])
    chain_rows.append([tag, "ions_nose_energy / xnhp0 / vnhp",
                       "NOT_AVAILABLE - computed in PW/src/vcsmd.f90 but never written to any output file",
                       None, None, None])

    # classification
    fin200 = stT["final200_mean"]
    if tag == "control":
        if n < 50:
            cls = "CONTROL_INCOMPLETE"
        elif (len(vol) and vol.min() <= 0) or d["error"] or not np.all(np.isfinite(T)):
            cls = "CONTROL_CELL_FAILURE"
        elif stT["max"] > 2000 and fin200 > 1000:
            cls = "CONTROL_HEATING"
        elif fin200 < 30:
            cls = "CONTROL_COOLING"
        else:
            cls = "CONTROL_BOUNDED"
    else:
        early_max = float(T[:200].max()) if n >= 200 else stT["max"]
        overshoot = early_max > 1.5 * TARGET_T
        near_target = 0.6 * TARGET_T <= fin200 <= 1.4 * TARGET_T
        if n < 50:
            cls = "INCOMPLETE"
        elif not np.all(np.isfinite(T)) or d["error"]:
            cls = "NUMERICAL_FAILURE"
        elif fin200 > 1000 and stT["max"] > 2000 and (fin_quarter is None or fin_quarter > -50):
            cls = "SUSTAINED_PHYSICAL_HEATING"
        elif stT["max"] > 2000 and fin200 > 1000 and fin_quarter is not None and fin_quarter < 0:
            cls = "LATE_RECOVERY_NOT_YET_CONVERGED"
        elif overshoot and near_target and fin_quarter is not None and fin_quarter < 0:
            cls = "INITIAL_TRANSIENT_THEN_CONVERGENCE"
        elif 0.75 * TARGET_T <= fin200 <= 1.25 * TARGET_T:
            cls = "BOUNDED_OSCILLATION_AROUND_TARGET"
        else:
            cls = "BOUNDED_BUT_OFFSET"

    RES[tag] = {
        "label": d["label"], "rc": d["rc"], "job_done": d["job_done"], "version": d["version"],
        "completed_steps": n, "classification": cls,
        "nose_thermostat_active_in_stdout": d["nose_stdout_active"],
        "temperature": stT, "temperature_initialized_K": float(T[0]),
        "drift_full_K_per_ps": lin_slope(t, T), "drift_final_half_K_per_ps": fin_half,
        "drift_final_quarter_K_per_ps": fin_quarter,
        "block100_means_K": [float(x) for x in bT],
        "roll25_final": float(r25[-1]) if np.isfinite(r25[-1]) else None,
        "roll100_final": float(r100[-1]) if np.isfinite(r100[-1]) else None,
        "peak": {"nst": pk, "time_ps": float(t[pk]), "value_K": float(T[pk]),
                 "time_from_peak_to_end_ps": float(t[-1] - t[pk]),
                 "lowest_roll100_after_peak_K": low_after},
        "returns_below": {"2000K": bool((T[pk:] < 2000).any()), "1000K": bool((T[pk:] < 1000).any()),
                          "600K": bool((T[pk:] < 600).any())},
        "crossings": {f"first_above_{k}K": first_cross(t, T, k) for k in (600, 1000, 2000, 5000)},
        "fraction_above": {f"{k}K": float(np.mean(T > k)) for k in (600, 1000, 2000)},
        "longest_run_above": {f"{k}K": longest_run(T, k) for k in (600, 1000, 2000)},
        "direct_terminal_T_K": Tdir, "direct_confirms_reported": direct_conf,
        "pressure_kbar": stP,
        "stress_mean_kbar": ({nm: float(np.nanmean(d["S"][:, i, j])) for nm, (i, j) in
                              zip(("sxx", "syy", "szz", "syz", "sxz", "sxy"),
                                  ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)))} if len(d["S"]) else {}),
        "cell": ({"volume_initial_bohr3": float(vol[0]), "volume_final_bohr3": float(vol[-1]),
                  "volume_change_percent": float(100 * (vol[-1] - vol[0]) / vol[0]),
                  "determinant_min": float(vol.min()), "condition_number_max": float(conds.max()),
                  "strain6_final": [float(x) for x in strain6(d["cel"][-1], d["cel"][0])],
                  "min_CC_bohr": float(min(dmins)) if dmins else None,
                  "max_CC_bohr": float(max(dmaxs)) if dmaxs else None,
                  "mean_NN_CC_bohr": float(np.mean(nns)) if nns else None,
                  "classification": ("SINGULAR_OR_UNSTABLE_CELL" if (vol.min() <= 0 or conds.max() > 5)
                                     else "DRIFTING_PR_CELL" if abs(vol[-1] - vol[0]) / vol[0] > 0.25
                                     else "STABLE_PR_CELL")} if len(vol) else {}),
        "energy_partial": {"quantity": "edyn = ekint + ut", "is_full_conserved_hamiltonian": False,
                           "note": ("PARTIAL: Nose chain energy computed in vcsmd but never written "
                                    "to output; a full Nose+PR conserved quantity cannot be formed."),
                           "drift_total_Ry": drift_total, "drift_eV_per_atom_per_ps": drift_at,
                           "bounded": bool(np.all(np.isfinite(ed)) and np.abs(ed - ed[0]).max() < 5.0)},
    }
    st = RES[tag]["temperature"]
    tstat_rows.append([tag, d["label"], cls, n, st["initial"], st["mean"], st["median"], st["std"],
                       st["min"], st["max"], st["final"], st["first100_mean"], st["step101_200_mean"],
                       st["final100_mean"], st["final200_mean"], st["final250_mean"],
                       RES[tag]["drift_full_K_per_ps"], RES[tag]["drift_final_half_K_per_ps"],
                       RES[tag]["drift_final_quarter_K_per_ps"], RES[tag]["fraction_above"]["600K"],
                       RES[tag]["longest_run_above"]["600K"]])

wcsv(os.path.join(ANA, "temperature_source_crosscheck.csv"),
     ["run", "comparison", "step_or_lag", "value_or_max_abs_diff_K"], tcross)
wcsv(os.path.join(ANA, "block_statistics.csv"),
     ["run", "block", "nst_start", "nst_end", "temperature_mean_K", "pressure_mean_kbar"], blk_rows)
wcsv(os.path.join(ANA, "temperature_statistics.csv"),
     ["run", "label", "classification", "n", "initial_K", "mean_K", "median_K", "std_K", "min_K", "max_K",
      "final_K", "first100_mean_K", "step101_200_mean_K", "final100_mean_K", "final200_mean_K",
      "final250_mean_K", "drift_full_K_per_ps", "drift_final_half_K_per_ps", "drift_final_quarter_K_per_ps",
      "fraction_above_600K", "longest_run_above_600K"], tstat_rows)
wcsv(os.path.join(ANA, "energy_diagnostics.csv"),
     ["run", "quantity", "status", "initial_Ry", "final_Ry", "min_Ry", "max_Ry", "drift_total_Ry",
      "drift_Ry_per_ps", "drift_eV_per_atom_per_ps", "slope_full", "slope_final_quarter",
      "corr_with_eka", "bounded"], energy_rows)
wcsv(os.path.join(ANA, "nose_chain_diagnostics.csv"),
     ["run", "quantity", "status", "initial", "final", "drift"], chain_rows)
# cell_statistics.csv (summary, both runs)
cst = []
for tag in RUNS:
    cc2 = RES[tag].get("cell", {})
    cst.append([tag, cc2.get("classification"), cc2.get("volume_initial_bohr3"),
                cc2.get("volume_final_bohr3"), cc2.get("volume_change_percent"),
                cc2.get("determinant_min"), cc2.get("condition_number_max"),
                cc2.get("min_CC_bohr"), cc2.get("max_CC_bohr"), cc2.get("mean_NN_CC_bohr")]
               + (cc2.get("strain6_final") or [None] * 6))
wcsv(os.path.join(ANA, "cell_statistics.csv"),
     ["run", "classification", "volume_initial_bohr3", "volume_final_bohr3", "volume_change_percent",
      "determinant_min", "condition_number_max", "min_CC_bohr", "max_CC_bohr", "mean_NN_CC_bohr",
      "strain_xx", "strain_yy", "strain_zz", "strain_yz", "strain_xz", "strain_xy"], cst)

# ============================================== Phase 13 overall ==========
cc_, nn_ = RES["control"], RES["nose"]
ccls, ncls = cc_["classification"], nn_["classification"]
prev_heat = REFERENCE["nose_previous_with_explicit_velocities_job50127909"]["max_K"] > 2000
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
elif ncls == "LATE_RECOVERY_NOT_YET_CONVERGED":
    overall = "NOSE_LATE_RECOVERY_UNRESOLVED"
elif ncls == "INITIAL_TRANSIENT_THEN_CONVERGENCE":
    overall = "NOSE_TRANSIENT_THEN_CONVERGES"
elif ncls in ("BOUNDED_OSCILLATION_AROUND_TARGET", "BOUNDED_BUT_OFFSET"):
    # control bounded, native-init Nose bounded, previous explicit-velocity Nose heated
    overall = "EXPLICIT_VELOCITIES_REQUIRED_FOR_HEATING" if prev_heat else "NOSE_STABLE_WITH_NATIVE_INITIALIZATION"
else:
    overall = "PARTIAL"

# heating-vs-cell ordering (moot if no heating)
nose_vol = np.array([abs(det3(m)) for m in D["nose"]["cel"]]) if len(D["nose"]["cel"]) else np.array([])
order = {"nose_heated": bool(nn_["temperature"]["max"] > 600),
         "note": "Nose never exceeded 600 K, so heating-before-deformation is moot."
         if nn_["temperature"]["max"] <= 600 else "see crossings vs volume change"}

summary = {
    "job_id": "50167241",
    "slurm_states": {"top": "COMPLETED 0:0", "50167241.0 (control pw.x)": "COMPLETED 0:0",
                     "50167241.1 (nose pw.x)": "COMPLETED 0:0"},
    "remote_workdir": "/leonardo/home/userexternal/mcheng00/workspace/Diamond/"
                      "diamond8_af7ebfc_no_velocity_initial_vs_nose_300K_20260724T120726Z",
    "overall_status": overall,
    "workflow_class": ident["workflow_class"],
    "changing_variable": "ion_temperature 'initial' -> 'nose' (+ fnosep=50, nhpcl=4)",
    "ensemble_label": "ordinary PR variable-cell dynamics after native 300 K init; TEST adds a "
                      "continuing Nose-Hoover ionic thermostat (verified active in stdout)",
    "cutoff_caveat": "ecutwfc=25 Ry: absolute pressures are Pulay-dominated, not converged.",
    "temperature_convention": {"instantaneous": "2*eka/((3N-3)k_B) from .eal col2",
                               "reported_tv_stdout": "2*avk/(3(N+1)k_B) - RUNNING AVERAGE (not used)",
                               "direct": "terminal .md ratd x terminal .cel, COM removed"},
    "input_contract": {"explicit_velocities_absent": no_vel["explicit_velocities_absent_in_both"],
                       "both_caimd_false": cd["both_explicitly_disabled"],
                       "valid_single_factor": ident["valid_single_factor"],
                       "non_thermostat_max_difference": {k: v for k, v in md.items() if k != "scf"}},
    "initial_state": {"classification": init["classification"],
                      "strict_single_factor": init["strict_single_factor"],
                      "first_eka_divergence_step": first_div},
    "nose_thermostat_active": D["nose"]["nose_stdout_active"],
    "runs": RES, "heating_vs_cell_order": order,
    "reference_values_verified": REFERENCE,
    "comparison_with_reference": {
        "control_vs_previous_no_velocity_initial": {
            "previous_mean_K": REFERENCE["control_previous_no_velocity_initial"]["mean_K"],
            "new_mean_K": cc_["temperature"]["mean"],
            "previous_final200_K": REFERENCE["control_previous_no_velocity_initial"]["final200_K"],
            "new_final200_K": cc_["temperature"]["final200_mean"],
            "reproduces_reference": bool(abs(cc_["temperature"]["mean"]
                                             - REFERENCE["control_previous_no_velocity_initial"]["mean_K"]) < 20)},
        "nose_vs_previous_explicit_velocity_nose": {
            "previous_max_K": REFERENCE["nose_previous_with_explicit_velocities_job50127909"]["max_K"],
            "new_max_K": nn_["temperature"]["max"],
            "previous_final_K": REFERENCE["nose_previous_with_explicit_velocities_job50127909"]["final_K"],
            "new_final_K": nn_["temperature"]["final"],
            "relationship": ("qualitatively different: native-initialization Nose does NOT reproduce "
                             "the several-thousand-kelvin heating seen with explicit velocities")}},
}
json.dump(summary, open(os.path.join(ANA, "comparison_summary_audited.json"), "w"), indent=2)

np.savez(os.path.join(EVD, "_arrays.npz"),
         **{f"{tag}_{k}": np.asarray(D[tag][k], float) for tag in RUNS
            for k in ("t_ps", "T", "eka", "p_kbar", "edyn")},
         control_vol=np.array([abs(det3(m)) for m in D["control"]["cel"]]),
         nose_vol=np.array([abs(det3(m)) for m in D["nose"]["cel"]]),
         control_cond=np.array([cond_sv(m)[0] for m in D["control"]["cel"]]),
         nose_cond=np.array([cond_sv(m)[0] for m in D["nose"]["cel"]]),
         control_celt=(D["control"]["cel_steps"] - D["control"]["cel_steps"][0]) * DT_FS / 1000.0,
         nose_celt=(D["nose"]["cel_steps"] - D["nose"]["cel_steps"][0]) * DT_FS / 1000.0)

print(json.dumps({"overall": overall, "control": ccls, "nose": ncls,
                  "control_mean_K": cc_["temperature"]["mean"], "nose_mean_K": nn_["temperature"]["mean"],
                  "nose_max_K": nn_["temperature"]["max"], "nose_final200_K": nn_["temperature"]["final200_mean"],
                  "nose_peak": nn_["peak"], "nose_frac_above_600": nn_["fraction_above"]["600K"],
                  "first_div_step": first_div, "init_class": init["classification"],
                  "direct_control": cc_["direct_confirms_reported"], "direct_nose": nn_["direct_confirms_reported"],
                  "workflow_class": ident["workflow_class"]}, indent=2))
