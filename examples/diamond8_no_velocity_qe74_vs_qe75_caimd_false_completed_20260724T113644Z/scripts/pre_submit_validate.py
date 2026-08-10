#!/usr/bin/env python3
"""Pre-submission validation for the no-explicit-velocity QE 7.4.1 vs QE 7.5 comparison.

Writes:
  evidence/no_explicit_velocity_audit.json
  evidence/caimd_keyword_audit.json
  evidence/input_identity_check.json
  evidence/input_normalized.diff

Exit status is non-zero if any hard requirement fails.
Standard library only; no network access.
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

QE74 = os.path.join(INP, "qe74_no_velocity.in")
QE75 = os.path.join(INP, "qe75_caimd_false_no_velocity.in")
PSEUDO = os.path.join(ROOT, "pseudo", "C.UPF")

EXPECTED_PSEUDO_SHA = "9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639"

FORBIDDEN_VELOCITY_PATTERNS = [
    r"atomic_velocities",
    r"ion_velocities",
    r"from_input",
    r"tv0rd",
    r"starting_velocities",
]
FORBIDDEN_CAIMD_PATTERNS = [r"caimd_stride", r"caimd_step_strain", r"caimd_mask"]

failures = []


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def strip_comments(text):
    out = []
    for line in text.splitlines():
        s = line.split("!")[0].rstrip()
        if s.strip():
            out.append(s)
    return out


def parse_namelists(text):
    """Return {namelist: {key: value}} plus the card section lines."""
    nml = {}
    cards = []
    current = None
    for line in strip_comments(text):
        st = line.strip()
        if st.startswith("&"):
            current = st[1:].strip().upper()
            nml[current] = {}
            continue
        if st == "/":
            current = None
            continue
        if current is not None:
            for part in st.split(","):
                part = part.strip()
                if not part or "=" not in part:
                    continue
                k, v = part.split("=", 1)
                nml[current][k.strip().lower()] = v.strip().rstrip(",").strip()
        else:
            cards.append(st)
    return nml, cards


def parse_cards(cards):
    """Extract species, cell, positions, k-points from card lines."""
    out = {"species": [], "cell": [], "positions": [], "kpoints": None,
           "cell_units": None, "position_units": None,
           "has_atomic_velocities": False}
    mode = None
    for line in cards:
        up = line.strip().upper()
        head = up.split()[0] if up.split() else ""
        if head == "ATOMIC_SPECIES":
            mode = "species"; continue
        if head == "CELL_PARAMETERS":
            mode = "cell"; out["cell_units"] = up.split()[1].lower() if len(up.split()) > 1 else None; continue
        if head == "ATOMIC_POSITIONS":
            mode = "pos"; out["position_units"] = up.split()[1].lower() if len(up.split()) > 1 else None; continue
        if head == "K_POINTS":
            mode = "kpoints"; out["kpoints"] = " ".join(up.split()[1:]).lower(); continue
        if head == "ATOMIC_VELOCITIES":
            out["has_atomic_velocities"] = True; mode = "vel"; continue
        f = line.split()
        if mode == "species" and len(f) >= 3:
            out["species"].append((f[0], float(f[1]), f[2]))
        elif mode == "cell" and len(f) == 3:
            out["cell"].append([float(x) for x in f])
        elif mode == "pos" and len(f) >= 4:
            out["positions"].append((f[0], [float(x) for x in f[1:4]]))
    return out


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def min_periodic_distance(cell_rows, frac):
    """Minimum periodic C-C distance in bohr. cell_rows are lattice vectors as rows."""
    best = float("inf")
    n = len(frac)
    shifts = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    for a, b in itertools.combinations(range(n), 2):
        fa, fb = frac[a][1], frac[b][1]
        for s in shifts:
            d = [fb[i] + s[i] - fa[i] for i in range(3)]
            cart = [sum(d[k] * cell_rows[k][i] for k in range(3)) for i in range(3)]
            r = math.sqrt(sum(c * c for c in cart))
            if r > 1e-9:
                best = min(best, r)
    # also self-images
    for a in range(n):
        for s in shifts:
            if s == (0, 0, 0):
                continue
            cart = [sum(s[k] * cell_rows[k][i] for k in range(3)) for i in range(3)]
            r = math.sqrt(sum(c * c for c in cart))
            best = min(best, r)
    return best


def num(v):
    try:
        return float(str(v).lower().replace("d", "e").strip().strip("'\""))
    except ValueError:
        return None


def norm_val(v):
    f = num(v)
    if f is not None:
        return f
    return str(v).strip().strip("'\"").lower()


def main():
    os.makedirs(EVID, exist_ok=True)
    t74, t75 = read(QE74), read(QE75)
    n74, c74 = parse_namelists(t74)
    n75, c75 = parse_namelists(t75)
    k74, k75 = parse_cards(c74), parse_cards(c75)

    # ---------------- no explicit velocity audit ----------------
    vel = {"expected_native_initialization_path":
           "PW/src/vcsmd.f90 -> PW/src/vcsubs.f90:vcinit -> ranv(temp=tempw) Maxwell "
           "distribution; input_velocities_read is FALSE because no ATOMIC_VELOCITIES "
           "card and no ion_velocities keyword are present.",
           "common_initialization_mode": "ion_temperature='initial' with tempw=300.0",
           "runs": {}}
    for name, text, cards in (("qe74", t74, k74), ("qe75", t75, k75)):
        # Check the COMPLETE raw text (case-insensitive) and, separately, the
        # effective input with comments stripped.
        raw_hits = {p: [i + 1 for i, l in enumerate(text.splitlines())
                        if re.search(p, l, re.IGNORECASE)]
                    for p in FORBIDDEN_VELOCITY_PATTERNS}
        raw_found = {p: v for p, v in raw_hits.items() if v}
        eff_lines = strip_comments(text)
        eff_hits = {p: [i + 1 for i, l in enumerate(eff_lines)
                        if re.search(p, l, re.IGNORECASE)]
                    for p in FORBIDDEN_VELOCITY_PATTERNS}
        eff_found = {p: v for p, v in eff_hits.items() if v}
        found = raw_found
        vel["runs"][name] = {
            "forbidden_keyword_hits_complete_text": raw_found,
            "forbidden_keyword_hits_comments_stripped": eff_found,
            "has_atomic_velocities_card": cards["has_atomic_velocities"],
            "restart_mode": norm_val(dict(itertools.chain(*[d.items() for d in
                                    ([parse_namelists(text)[0].get("CONTROL", {})])]))
                                     .get("restart_mode", "")),
            "ion_temperature": norm_val(parse_namelists(text)[0].get("IONS", {}).get("ion_temperature", "")),
            "tempw": num(parse_namelists(text)[0].get("IONS", {}).get("tempw", "")),
        }
        if found:
            failures.append("%s: forbidden velocity keyword(s) present: %s" % (name, list(found)))
        if cards["has_atomic_velocities"]:
            failures.append("%s: ATOMIC_VELOCITIES card present" % name)
        if vel["runs"][name]["restart_mode"] != "from_scratch":
            failures.append("%s: restart_mode is not from_scratch" % name)
        if vel["runs"][name]["ion_temperature"] != "initial":
            failures.append("%s: ion_temperature is not 'initial'" % name)
        if vel["runs"][name]["tempw"] != 300.0:
            failures.append("%s: tempw is not 300.0" % name)
    vel["all_checks_passed"] = not failures
    json.dump(vel, open(os.path.join(EVID, "no_explicit_velocity_audit.json"), "w"), indent=2)

    # ---------------- CAIMD keyword audit ----------------
    caimd = {"policy": "QE 7.5 input must contain only caimd=.false.; "
                       "QE 7.4.1 input must contain no CAIMD keyword at all.",
             "runs": {}}
    for name, text in (("qe74", t74), ("qe75", t75)):
        low_lines = [l for l in strip_comments(text)]
        caimd_lines = [l.strip() for l in low_lines if "caimd" in l.lower()]
        bad = [p for p in FORBIDDEN_CAIMD_PATTERNS
               if any(re.search(p, l, re.IGNORECASE) for l in low_lines)]
        active = [l for l in caimd_lines if re.search(r"caimd\s*=\s*\.?true", l, re.IGNORECASE)]
        caimd["runs"][name] = {"caimd_lines": caimd_lines,
                               "forbidden_caimd_parameters": bad,
                               "active_caimd_lines": active}
        if bad:
            failures.append("%s: forbidden CAIMD parameter(s): %s" % (name, bad))
        if active:
            failures.append("%s: active caimd=.true. present" % name)
    if caimd["runs"]["qe74"]["caimd_lines"]:
        failures.append("qe74: contains CAIMD keyword(s), not supported by QE 7.4.1")
    q75 = caimd["runs"]["qe75"]["caimd_lines"]
    if len(q75) != 1 or not re.search(r"caimd\s*=\s*\.false\.", q75[0], re.IGNORECASE):
        failures.append("qe75: expected exactly one 'caimd = .false.' line, got %r" % q75)
    caimd["all_checks_passed"] = not failures
    json.dump(caimd, open(os.path.join(EVID, "caimd_keyword_audit.json"), "w"), indent=2)

    # ---------------- input identity ----------------
    ALLOWED_DIFF_KEYS = {("CONTROL", "prefix"), ("CONTROL", "outdir"), ("CELL", "caimd")}
    diffs = {}
    for sec in sorted(set(n74) | set(n75)):
        a, b = n74.get(sec, {}), n75.get(sec, {})
        for key in sorted(set(a) | set(b)):
            va, vb = a.get(key, "<absent>"), b.get(key, "<absent>")
            if norm_val(va) != norm_val(vb) and (sec, key) not in ALLOWED_DIFF_KEYS:
                diffs["%s.%s" % (sec, key)] = {"qe74": va, "qe75": vb}
    if diffs:
        failures.append("namelist differences outside the allowed set: %s" % list(diffs))

    def cmp_group(label, va, vb):
        same = va == vb
        if not same:
            failures.append("%s differs between inputs" % label)
        return {"identical": same, "qe74": va, "qe75": vb}

    cell = k74["cell"]
    frac = k74["positions"]
    determinant = det3(cell)
    dmin = min_periodic_distance(cell, frac)
    dup = False
    for i, j in itertools.combinations(range(len(frac)), 2):
        if all(abs(frac[i][1][k] - frac[j][1][k]) < 1e-8 for k in range(3)):
            dup = True
    if determinant <= 0:
        failures.append("cell determinant is not positive")
    if not (2.5 < dmin < 3.5):
        failures.append("minimum periodic C-C distance %.4f bohr outside expected diamond range" % dmin)
    if dup:
        failures.append("duplicate atomic positions detected")
    if len(frac) != 8:
        failures.append("atom count is not 8")

    identity = {
        "cell": cmp_group("cell", k74["cell"], k75["cell"]),
        "cell_units": cmp_group("cell_units", k74["cell_units"], k75["cell_units"]),
        "positions": cmp_group("positions", k74["positions"], k75["positions"]),
        "position_units": cmp_group("position_units", k74["position_units"], k75["position_units"]),
        "species_and_masses": cmp_group("species", k74["species"], k75["species"]),
        "kpoints": cmp_group("kpoints", k74["kpoints"], k75["kpoints"]),
        "namelist_differences_outside_allowed_set": diffs,
        "allowed_differences": sorted("%s.%s" % t for t in ALLOWED_DIFF_KEYS),
        "structure_checks": {
            "n_atoms": len(frac),
            "cell_determinant_bohr3": determinant,
            "cell_determinant_positive": determinant > 0,
            "min_periodic_CC_distance_bohr": dmin,
            "min_periodic_CC_distance_angstrom": dmin * 0.529177210903,
            "duplicate_atoms": dup,
            "gamma_point": (k74["kpoints"] or "").startswith("gamma"),
            "structure_provenance": "ideal unstrained Diamond8 conventional cell, "
                                    "a=3.567 Angstrom (DIAMOND8_SETUP.json); "
                                    "not a terminal state from a failed trajectory",
        },
        "matched_physical_fields": {
            "cutoffs": {"ecutwfc": num(n74["SYSTEM"]["ecutwfc"]), "ecutrho": num(n74["SYSTEM"]["ecutrho"])},
            "timestep_qe": num(n74["CONTROL"]["dt"]),
            "timestep_fs": 2.0,
            "nstep": int(num(n74["CONTROL"]["nstep"])),
            "temperature_mode": norm_val(n74["IONS"]["ion_temperature"]),
            "tempw": num(n74["IONS"]["tempw"]),
            "cell_dynamics": norm_val(n74["CELL"]["cell_dynamics"]),
            "cell_dofree": norm_val(n74["CELL"]["cell_dofree"]),
            "press": num(n74["CELL"]["press"]),
            "wmass": num(n74["CELL"]["wmass"]),
            "conv_thr": num(n74["ELECTRONS"]["conv_thr"]),
            "mixing_beta": num(n74["ELECTRONS"]["mixing_beta"]),
            "electron_maxstep": int(num(n74["ELECTRONS"]["electron_maxstep"])),
        },
        "checksums": {
            "qe74_input_sha256": sha256(QE74),
            "qe75_input_sha256": sha256(QE75),
            "pseudo_sha256": sha256(PSEUDO),
            "pseudo_sha256_expected": EXPECTED_PSEUDO_SHA,
            "pseudo_matches_expected": sha256(PSEUDO) == EXPECTED_PSEUDO_SHA,
        },
    }
    if sha256(PSEUDO) != EXPECTED_PSEUDO_SHA:
        failures.append("pseudopotential checksum mismatch")
    identity["all_checks_passed"] = not failures
    json.dump(identity, open(os.path.join(EVID, "input_identity_check.json"), "w"), indent=2)

    # ---------------- normalized diff ----------------
    def normalize(text):
        return [l for l in strip_comments(text)]
    d = list(difflib.unified_diff(normalize(t74), normalize(t75),
                                  fromfile="qe74_no_velocity.in (comments stripped)",
                                  tofile="qe75_caimd_false_no_velocity.in (comments stripped)",
                                  lineterm=""))
    with open(os.path.join(EVID, "input_normalized.diff"), "w") as fh:
        fh.write("\n".join(d) + "\n")

    print(json.dumps({
        "failures": failures,
        "n_atoms": len(frac),
        "determinant_bohr3": determinant,
        "min_CC_bohr": dmin,
        "min_CC_angstrom": dmin * 0.529177210903,
        "namelist_diffs_outside_allowed": diffs,
        "normalized_diff_lines": len(d),
    }, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
