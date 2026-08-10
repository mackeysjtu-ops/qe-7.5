#!/usr/bin/env python3
"""Pre-submission validation for the three-case CAIMD zero-strain 5000-step test.

Writes:
  evidence/no_explicit_velocity_audit.json
  evidence/input_identity_audit.json
  evidence/input_A_vs_B_normalized.diff
  evidence/input_A_vs_C_normalized.diff
Non-zero exit if any hard requirement fails. Standard library only.
"""
import difflib
import hashlib
import itertools
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INP = os.path.join(ROOT, "inputs")
EVID = os.path.join(ROOT, "evidence")
FILES = {
    "A": os.path.join(INP, "case_A_caimd_false_nose_5000.in"),
    "B": os.path.join(INP, "case_B_caimd_true_all_false_nose_5000.in"),
    "C": os.path.join(INP, "case_C_caimd_true_all_fixed_zero_strain_nose_5000.in"),
}
PSEUDO = os.path.join(ROOT, "pseudo", "C.UPF")
PSEUDO_EXPECTED = "9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639"

FORBID_VEL = [r"atomic_velocities", r"ion_velocities", r"from_input", r"tv0rd", r"starting_velocities"]
# CAIMD lines are the intended per-case difference
ALLOWED_KEYS = {("CONTROL", "prefix"), ("CONTROL", "outdir"),
                ("CELL", "caimd"), ("CELL", "caimd_stride")}
for i in range(1, 7):
    ALLOWED_KEYS.add(("CELL", f"caimd_step_strain({i})"))
    ALLOWED_KEYS.add(("CELL", f"caimd_mask({i})"))

failures = []


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(65536), b""):
            h.update(c)
    return h.hexdigest()


def read(p):
    return open(p, encoding="utf-8").read()


def strip_comments(t):
    return [l.split("!")[0].rstrip() for l in t.splitlines() if l.split("!")[0].strip()]


def parse(t):
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


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def min_cc(cell_rows, frac):
    best = math.inf
    shifts = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    n = len(frac)
    for a, b in itertools.combinations(range(n), 2):
        for s in shifts:
            d = [frac[b][i] + s[i] - frac[a][i] for i in range(3)]
            c = [sum(d[k] * cell_rows[k][i] for k in range(3)) for i in range(3)]
            r = math.sqrt(sum(x * x for x in c))
            if r > 1e-9:
                best = min(best, r)
    for s in shifts:
        if s != (0, 0, 0):
            c = [sum(s[k] * cell_rows[k][i] for k in range(3)) for i in range(3)]
            best = min(best, math.sqrt(sum(x * x for x in c)))
    return best


TXT = {k: read(v) for k, v in FILES.items()}
NML = {k: parse(v)[0] for k, v in TXT.items()}
CRD = {k: cards_struct(parse(v)[1]) for k, v in TXT.items()}


def main():
    os.makedirs(EVID, exist_ok=True)

    # ---- no explicit velocity ----
    nv = {"native_path": "vcinit->ranv(temp=tempw); ranv zeroes total linear momentum", "cases": {}}
    for k in FILES:
        t = TXT[k]
        raw = {p: [i + 1 for i, l in enumerate(t.splitlines()) if re.search(p, l, re.I)] for p in FORBID_VEL}
        found = {a: b for a, b in raw.items() if b}
        nv["cases"][k] = {"forbidden_hits_complete_text": found,
                          "atomic_velocities_header": CRD[k]["atomic_velocities"],
                          "velocity_rows": CRD[k]["vel_rows"],
                          "restart_mode": numv(NML[k].get("CONTROL", {}).get("restart_mode", "")),
                          "tempw": numv(NML[k].get("IONS", {}).get("tempw", "")),
                          "ion_temperature": numv(NML[k].get("IONS", {}).get("ion_temperature", ""))}
        r = nv["cases"][k]
        if found:
            failures.append(f"{k}: forbidden velocity keyword {list(found)}")
        if r["atomic_velocities_header"] or r["velocity_rows"]:
            failures.append(f"{k}: ATOMIC_VELOCITIES present")
        if r["restart_mode"] != "from_scratch":
            failures.append(f"{k}: restart_mode != from_scratch")
        if r["ion_temperature"] != "nose":
            failures.append(f"{k}: ion_temperature != nose")
        if r["tempw"] != 300.0:
            failures.append(f"{k}: tempw != 300.0")
    nv["explicit_velocities_absent_in_all"] = not failures
    nv["native_initialization_conclusion"] = ("all three cases start from native vcinit->ranv(300 K); "
                                              "no explicit velocity card or keyword present")
    json.dump(nv, open(os.path.join(EVID, "no_explicit_velocity_audit.json"), "w"), indent=2)

    # ---- CAIMD per-case settings ----
    def caimd_of(k):
        cell = NML[k].get("CELL", {})
        return {"caimd": numv(cell.get("caimd", "<absent>")),
                "caimd_stride": numv(cell.get("caimd_stride", "<absent>")),
                "mask": [numv(cell.get(f"caimd_mask({i})", "<absent>")) for i in range(1, 7)],
                "step_strain": [numv(cell.get(f"caimd_step_strain({i})", "<absent>")) for i in range(1, 7)]}
    caimd = {k: caimd_of(k) for k in FILES}
    # CASE A
    if caimd["A"]["caimd"] not in (False, "false", ".false."):
        failures.append("A: caimd is not .false.")
    # CASE B: enabled, stride 1, all masks false, all step_strain 0
    b = caimd["B"]
    if b["caimd"] not in (True, "true", ".true."):
        failures.append("B: caimd is not .true.")
    if b["caimd_stride"] != 1.0:
        failures.append("B: caimd_stride != 1")
    if any(m not in (False, "false", ".false.") for m in b["mask"]):
        failures.append("B: not all masks false")
    if any(s != 0.0 for s in b["step_strain"]):
        failures.append("B: not all step_strain zero")
    # CASE C: enabled, stride 1, all masks true, all step_strain 0
    c = caimd["C"]
    if c["caimd"] not in (True, "true", ".true."):
        failures.append("C: caimd is not .true.")
    if c["caimd_stride"] != 1.0:
        failures.append("C: caimd_stride != 1")
    if any(m not in (True, "true", ".true.") for m in c["mask"]):
        failures.append("C: not all masks true")
    if any(s != 0.0 for s in c["step_strain"]):
        failures.append("C: not all step_strain zero")
    # no nonzero strain anywhere
    for k in FILES:
        if any(isinstance(s, float) and s != 0.0 for s in caimd[k]["step_strain"]):
            failures.append(f"{k}: nonzero strain increment present")

    # ---- input identity (A as reference; ignore only allowed keys) ----
    diffs = {}
    ka = CRD["A"]
    for other in ("B", "C"):
        ko = CRD[other]
        na, no = NML["A"], NML[other]
        d = {}
        for sec in sorted(set(na) | set(no)):
            aa, bb = na.get(sec, {}), no.get(sec, {})
            for key in sorted(set(aa) | set(bb)):
                if numv(aa.get(key, "<absent>")) != numv(bb.get(key, "<absent>")) and (sec, key) not in ALLOWED_KEYS:
                    d[f"{sec}.{key}"] = {"A": aa.get(key, "<absent>"), other: bb.get(key, "<absent>")}
        if ka["cell"] != ko["cell"]:
            d["CELL_PARAMETERS"] = "differ"
        if ka["pos"] != ko["pos"]:
            d["ATOMIC_POSITIONS"] = "differ"
        if [p[0] for p in ka["pos"]] != [p[0] for p in ko["pos"]]:
            d["atom_order"] = "differ"
        if ka["species"] != ko["species"]:
            d["ATOMIC_SPECIES"] = "differ"
        if ka["kpoints"] != ko["kpoints"]:
            d["K_POINTS"] = "differ"
        if d:
            diffs[f"A_vs_{other}"] = d
            failures.append(f"A_vs_{other}: non-CAIMD differences {list(d)}")
        with open(os.path.join(EVID, f"input_A_vs_{other}_normalized.diff"), "w") as fh:
            fh.write("\n".join(difflib.unified_diff(strip_comments(TXT["A"]), strip_comments(TXT[other]),
                     fromfile="case_A (normalized)", tofile=f"case_{other} (normalized)", lineterm="")) + "\n")

    det = det3(ka["cell"])
    dmin = min_cc(ka["cell"], [p[1] for p in ka["pos"]])
    dup = any(all(abs(ka["pos"][i][1][x] - ka["pos"][j][1][x]) < 1e-8 for x in range(3))
              for i, j in itertools.combinations(range(len(ka["pos"])), 2))
    for k in FILES:
        if len(CRD[k]["pos"]) != 8:
            failures.append(f"{k}: atom count != 8")
    if det <= 0:
        failures.append("cell determinant not positive")
    if not (2.5 < dmin < 3.5):
        failures.append(f"min C-C {dmin:.4f} bohr out of range")
    if dup:
        failures.append("duplicate atoms")
    if abs(det - 306.2710374914433) > 1e-6:
        failures.append("structure is not the ideal undeformed cell")
    if sha(PSEUDO) != PSEUDO_EXPECTED:
        failures.append("pseudo checksum mismatch")

    ident = {"reference_case": "A", "caimd_settings": caimd,
             "non_caimd_differences": diffs,
             "allowed_keys": sorted(f"{s}.{k}" for s, k in ALLOWED_KEYS),
             "shared_physical": {
                 "nstep": [numv(NML[k]["CONTROL"]["nstep"]) for k in FILES],
                 "dt": [numv(NML[k]["CONTROL"]["dt"]) for k in FILES],
                 "ecutwfc": [numv(NML[k]["SYSTEM"]["ecutwfc"]) for k in FILES],
                 "ecutrho": [numv(NML[k]["SYSTEM"]["ecutrho"]) for k in FILES],
                 "fnosep(1)": [numv(NML[k]["IONS"]["fnosep(1)"]) for k in FILES],
                 "nhpcl": [numv(NML[k]["IONS"]["nhpcl"]) for k in FILES],
                 "tempw": [numv(NML[k]["IONS"]["tempw"]) for k in FILES],
                 "press": [numv(NML[k]["CELL"]["press"]) for k in FILES],
                 "wmass": [numv(NML[k]["CELL"]["wmass"]) for k in FILES],
                 "cell_dynamics": [numv(NML[k]["CELL"]["cell_dynamics"]) for k in FILES],
                 "cell_dofree": [numv(NML[k]["CELL"]["cell_dofree"]) for k in FILES]},
             "structure": {"n_atoms": [len(CRD[k]["pos"]) for k in FILES],
                           "determinant_bohr3": det, "min_CC_bohr": dmin,
                           "min_CC_angstrom": dmin * 0.529177210903, "duplicate_atoms": dup,
                           "undeformed_ideal": abs(det - 306.2710374914433) < 1e-6},
             "checksums": {k: sha(FILES[k]) for k in FILES} | {"pseudo": sha(PSEUDO),
                          "pseudo_matches": sha(PSEUDO) == PSEUDO_EXPECTED},
             "all_checks_passed": not failures}
    json.dump(ident, open(os.path.join(EVID, "input_identity_audit.json"), "w"), indent=2)

    print(json.dumps({"failures": failures, "determinant": det, "min_CC_bohr": dmin,
                      "caimd_A": caimd["A"]["caimd"], "caimd_B_masks": caimd["B"]["mask"],
                      "caimd_C_masks": caimd["C"]["mask"],
                      "all_step_strain_zero": all(all(s == 0.0 for s in caimd[k]["step_strain"])
                                                  for k in ("B", "C"))}, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
