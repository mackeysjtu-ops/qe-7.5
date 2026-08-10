#!/usr/bin/env python3
"""Pre-submission validation for the single-factor thermostat isolation.

Writes:
  evidence/no_explicit_velocity_audit.json
  evidence/caimd_disabled_audit.json
  evidence/input_identity_audit.json
  evidence/input_normalized.diff

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
CTRL = os.path.join(INP, "control_initial_no_velocity.in")
NOSE = os.path.join(INP, "test_nose_no_velocity.in")
PSEUDO = os.path.join(ROOT, "pseudo", "C.UPF")
PSEUDO_EXPECTED = "9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639"

FORBID_VEL = [r"atomic_velocities", r"ion_velocities", r"from_input", r"tv0rd",
              r"starting_velocities"]
FORBID_CAIMD = [r"caimd_stride", r"caimd_step_strain", r"caimd_mask"]
# thermostat lines are the intended changing variable
ALLOWED_KEYS = {("CONTROL", "prefix"), ("CONTROL", "outdir"),
                ("IONS", "ion_temperature"), ("IONS", "fnosep(1)"),
                ("IONS", "fnosep"), ("IONS", "nhpcl")}

failures = []


def sha256(p):
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
           "cell_units": None, "pos_units": None, "atomic_velocities": False,
           "velocity_rows": 0}
    mode = None
    for line in cards:
        u = line.upper().split(); h = u[0] if u else ""
        if h == "ATOMIC_SPECIES": mode = "sp"; continue
        if h == "CELL_PARAMETERS":
            mode = "ce"; out["cell_units"] = u[1].lower() if len(u) > 1 else None; continue
        if h == "ATOMIC_POSITIONS":
            mode = "po"; out["pos_units"] = u[1].lower() if len(u) > 1 else None; continue
        if h == "K_POINTS":
            mode = "kp"; out["kpoints"] = " ".join(u[1:]).lower(); continue
        if h == "ATOMIC_VELOCITIES":
            out["atomic_velocities"] = True; mode = "ve"; continue
        f = line.split()
        if mode == "sp" and len(f) >= 3: out["species"].append((f[0], float(f[1]), f[2]))
        elif mode == "ce" and len(f) == 3: out["cell"].append([float(x) for x in f])
        elif mode == "po" and len(f) >= 4: out["pos"].append((f[0], [float(x) for x in f[1:4]]))
        elif mode == "ve" and len(f) >= 4: out["velocity_rows"] += 1
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
        if s == (0, 0, 0):
            continue
        c = [sum(s[k] * cell_rows[k][i] for k in range(3)) for i in range(3)]
        best = min(best, math.sqrt(sum(x * x for x in c)))
    return best


def main():
    os.makedirs(EVID, exist_ok=True)
    tc, tn = read(CTRL), read(NOSE)
    nc, cc = parse(tc)
    nn, cn = parse(tn)
    kc, kn = cards_struct(cc), cards_struct(cn)

    # ---------------- no explicit velocity ----------------
    vel = {"native_initialization_path":
           "PW/src/vcsmd.f90 -> PW/src/vcsubs.f90:vcinit -> ranv(temp=tempw); "
           "ranv builds a Maxwell distribution with zero total linear momentum. "
           "input_velocities_read is FALSE because no velocity card/keyword is present.",
           "runs": {}}
    for tag, t, k in (("control", tc, kc), ("nose", tn, kn)):
        raw = {p: [i + 1 for i, l in enumerate(t.splitlines()) if re.search(p, l, re.I)]
               for p in FORBID_VEL}
        eff_lines = strip_comments(t)
        eff = {p: [i + 1 for i, l in enumerate(eff_lines) if re.search(p, l, re.I)]
               for p in FORBID_VEL}
        ctl = parse(t)[0]
        vel["runs"][tag] = {
            "forbidden_hits_complete_text": {a: b for a, b in raw.items() if b},
            "forbidden_hits_comments_stripped": {a: b for a, b in eff.items() if b},
            "atomic_velocities_header": k["atomic_velocities"],
            "velocity_rows": k["velocity_rows"],
            "restart_mode": numv(ctl.get("CONTROL", {}).get("restart_mode", "")),
            "restart_state_input": bool(
                re.search(r"startingwfc\s*=\s*'file'|startingpot\s*=\s*'file'", t, re.I)),
        }
        r = vel["runs"][tag]
        if r["forbidden_hits_complete_text"]:
            failures.append(f"{tag}: forbidden velocity keyword {list(r['forbidden_hits_complete_text'])}")
        if r["atomic_velocities_header"] or r["velocity_rows"]:
            failures.append(f"{tag}: ATOMIC_VELOCITIES card present")
        if r["restart_mode"] != "from_scratch":
            failures.append(f"{tag}: restart_mode != from_scratch")
        if r["restart_state_input"]:
            failures.append(f"{tag}: restart-state input detected")
    vel["explicit_velocities_absent_in_both"] = not failures
    json.dump(vel, open(os.path.join(EVID, "no_explicit_velocity_audit.json"), "w"), indent=2)

    # ---------------- caimd disabled ----------------
    cd = {"policy": "both inputs must contain exactly caimd=.false. and no active CAIMD keyword",
          "runs": {}}
    for tag, t in (("control", tc), ("nose", tn)):
        eff = strip_comments(t)
        lines = [l.strip() for l in eff if "caimd" in l.lower()]
        bad = [p for p in FORBID_CAIMD if any(re.search(p, l, re.I) for l in eff)]
        active = [l for l in lines if re.search(r"caimd\s*=\s*\.?true", l, re.I)]
        false_ok = len(lines) == 1 and bool(re.search(r"caimd\s*=\s*\.false\.", lines[0], re.I))
        cd["runs"][tag] = {"caimd_lines": lines, "forbidden_parameters": bad,
                           "active_caimd": active, "exactly_caimd_false": false_ok}
        if bad:
            failures.append(f"{tag}: forbidden CAIMD parameter {bad}")
        if active:
            failures.append(f"{tag}: active caimd=.true.")
        if not false_ok:
            failures.append(f"{tag}: expected exactly one 'caimd = .false.' line, got {lines}")
    cd["both_explicitly_disabled"] = all(v["exactly_caimd_false"] for v in cd["runs"].values())
    json.dump(cd, open(os.path.join(EVID, "caimd_disabled_audit.json"), "w"), indent=2)

    # ---------------- input identity ----------------
    diffs = {}
    for sec in sorted(set(nc) | set(nn)):
        a, b = nc.get(sec, {}), nn.get(sec, {})
        for key in sorted(set(a) | set(b)):
            va, vb = a.get(key, "<absent>"), b.get(key, "<absent>")
            if numv(va) != numv(vb) and (sec, key) not in ALLOWED_KEYS:
                diffs[f"{sec}.{key}"] = {"control": va, "nose": vb}
    if diffs:
        failures.append(f"non-thermostat namelist differences: {list(diffs)}")

    def maxdiff(A, B):
        if len(A) != len(B):
            return None
        m = 0.0
        for ra, rb in zip(A, B):
            for x, y in zip(ra, rb):
                m = max(m, abs(x - y))
        return m

    cell_md = maxdiff(kc["cell"], kn["cell"])
    pos_md = maxdiff([p[1] for p in kc["pos"]], [p[1] for p in kn["pos"]])
    order_same = [p[0] for p in kc["pos"]] == [p[0] for p in kn["pos"]]
    mass_md = maxdiff([[s[1]] for s in kc["species"]], [[s[1]] for s in kn["species"]])
    dmin = min_cc(kc["cell"], [p[1] for p in kc["pos"]])
    det = det3(kc["cell"])

    def fnum(sec, key):
        return numv(nc.get(sec, {}).get(key, "")), numv(nn.get(sec, {}).get(key, ""))

    scf_keys = [("ELECTRONS", "conv_thr"), ("ELECTRONS", "electron_maxstep"),
                ("ELECTRONS", "mixing_beta")]
    ident = {
        "max_differences": {
            "cell_bohr": cell_md, "positions_crystal": pos_md,
            "atom_ordering_identical": order_same, "masses_amu": mass_md,
            "ecutwfc": abs(fnum("SYSTEM", "ecutwfc")[0] - fnum("SYSTEM", "ecutwfc")[1]),
            "ecutrho": abs(fnum("SYSTEM", "ecutrho")[0] - fnum("SYSTEM", "ecutrho")[1]),
            "timestep_dt": abs(fnum("CONTROL", "dt")[0] - fnum("CONTROL", "dt")[1]),
            "nstep": abs(fnum("CONTROL", "nstep")[0] - fnum("CONTROL", "nstep")[1]),
            "press": abs(fnum("CELL", "press")[0] - fnum("CELL", "press")[1]),
            "wmass": abs(fnum("CELL", "wmass")[0] - fnum("CELL", "wmass")[1]),
            "scf_parameters": {f"{s}.{k}": abs(fnum(s, k)[0] - fnum(s, k)[1]) for s, k in scf_keys},
            "tempw": abs(fnum("IONS", "tempw")[0] - fnum("IONS", "tempw")[1]),
        },
        "kpoints": {"control": kc["kpoints"], "nose": kn["kpoints"],
                    "identical": kc["kpoints"] == kn["kpoints"]},
        "units": {"cell": [kc["cell_units"], kn["cell_units"]],
                  "positions": [kc["pos_units"], kn["pos_units"]],
                  "identical": (kc["cell_units"] == kn["cell_units"]
                                and kc["pos_units"] == kn["pos_units"])},
        "intended_differences": {
            "ion_temperature": [numv(nc["IONS"].get("ion_temperature")),
                                numv(nn["IONS"].get("ion_temperature"))],
            "fnosep(1)": [nc["IONS"].get("fnosep(1)", "<absent>"),
                          nn["IONS"].get("fnosep(1)", "<absent>")],
            "nhpcl": [nc["IONS"].get("nhpcl", "<absent>"), nn["IONS"].get("nhpcl", "<absent>")],
            "prefix": [nc["CONTROL"].get("prefix"), nn["CONTROL"].get("prefix")],
        },
        "namelist_differences_outside_allowed_set": diffs,
        "structure": {
            "n_atoms": {"control": len(kc["pos"]), "nose": len(kn["pos"])},
            "cell_determinant_bohr3": det,
            "determinant_positive": det > 0,
            "min_periodic_CC_bohr": dmin,
            "min_periodic_CC_angstrom": dmin * 0.529177210903,
            "undeformed_ideal_structure": abs(det - 306.2710374914433) < 1e-6,
        },
        "checksums": {"control_sha256": sha256(CTRL), "nose_sha256": sha256(NOSE),
                      "pseudo_sha256": sha256(PSEUDO),
                      "pseudo_expected": PSEUDO_EXPECTED,
                      "pseudo_matches": sha256(PSEUDO) == PSEUDO_EXPECTED},
    }
    md = ident["max_differences"]
    for label, v in (("cell", cell_md), ("positions", pos_md), ("masses", mass_md),
                     ("ecutwfc", md["ecutwfc"]), ("ecutrho", md["ecutrho"]),
                     ("dt", md["timestep_dt"]), ("nstep", md["nstep"]),
                     ("press", md["press"]), ("wmass", md["wmass"]),
                     ("tempw", md["tempw"])):
        if v is None or v != 0:
            failures.append(f"non-zero physical difference in {label}: {v}")
    for k, v in md["scf_parameters"].items():
        if v != 0:
            failures.append(f"non-zero SCF difference in {k}: {v}")
    if not order_same:
        failures.append("atom ordering differs")
    if len(kc["pos"]) != 8 or len(kn["pos"]) != 8:
        failures.append("atom count is not 8")
    if not ident["structure"]["undeformed_ideal_structure"]:
        failures.append("structure is not the ideal undeformed reference cell")
    if not ident["checksums"]["pseudo_matches"]:
        failures.append("pseudopotential checksum mismatch")
    ident["all_checks_passed"] = not failures
    json.dump(ident, open(os.path.join(EVID, "input_identity_audit.json"), "w"), indent=2)

    with open(os.path.join(EVID, "input_normalized.diff"), "w") as fh:
        fh.write("\n".join(difflib.unified_diff(
            strip_comments(tc), strip_comments(tn),
            fromfile="control_initial_no_velocity.in (normalized)",
            tofile="test_nose_no_velocity.in (normalized)", lineterm="")) + "\n")

    print(json.dumps({"failures": failures,
                      "max_diffs": {k: v for k, v in md.items() if k != "scf_parameters"},
                      "min_CC_bohr": dmin, "determinant": det,
                      "intended_differences": ident["intended_differences"]}, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
