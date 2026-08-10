#!/usr/bin/env python3
"""Recompute key CAIMD validation metrics from raw evidence copies.

Reads only files under ../assets/data (raw copies made from the local
validation workflows) and writes ../assets/data/recomputed_metrics.json.

Temperature convention: T = 2*eka / (ndof * k_B[Ry]) with ndof = 3*(nat-1) = 21
for the 8-atom diamond cell, matching the QE PW vc-md `.eal` ionic kinetic
energy column (columns: uta eka eta utl ekla etl nst).

No QE source file, raw workflow output, or existing report is modified.
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "assets", "data")

K_BOLTZMANN_RY = 8.617333262e-5 / 13.605693122994  # eV/K / (Ry in eV)
NDOF = 21          # 3*(nat-1), nat=8
DT_PS = 0.002      # 2 fs ionic timestep in every Diamond8 run considered here


def parse_eal(path):
    """Return list of (step, eka_Ry, temperature_K)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.replace("D", "E").split()
            if len(parts) != 7:
                continue
            eka = float(parts[1])
            step = int(parts[6])
            temp = 2.0 * eka / (NDOF * K_BOLTZMANN_RY)
            rows.append((step, eka, temp))
    return rows


def parse_caimd(path):
    """Return dict of column arrays from a 26-column .caimd stream."""
    cols = ("ionic_step time_ps events event_flag "
            "nom_xx nom_yy nom_zz nom_yz nom_xz nom_xy "
            "act_xx act_yy act_zz act_yz act_xz act_xy "
            "s_xx s_yy s_zz s_yz s_xz s_xy "
            "temperature_K volume velocity_jump_max velocity_jump_rms").split()
    data = {c: [] for c in cols}
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 26:
                continue
            for i, c in enumerate(cols):
                if c == "event_flag":
                    data[c].append(parts[i] == "T")
                elif c in ("ionic_step", "events"):
                    data[c].append(int(parts[i]))
                else:
                    data[c].append(float(parts[i]))
    return data


def parse_cel(path):
    """Return list of (step, 3x3 matrix) from a QE vc-md .cel history."""
    frames = []
    with open(path) as fh:
        lines = [ln.split() for ln in fh if ln.strip()]
    i = 0
    while i + 3 < len(lines) + 1 and i + 3 <= len(lines):
        head = lines[i]
        if len(head) != 2:
            i += 1
            continue
        step = int(head[0])
        mat = [[float(x) for x in lines[i + r]] for r in (1, 2, 3)]
        frames.append((step, mat))
        i += 4
    return frames


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def cond2(m):
    """2-norm condition number via singular values of a 3x3 matrix."""
    import numpy as np
    s = np.linalg.svd(np.array(m), compute_uv=False)
    return float(s[0] / s[-1])


def temp_stats(rows, skip_first=0):
    temps = [t for _, _, t in rows[skip_first:]]
    n = len(temps)
    mean = sum(temps) / n
    var = sum((t - mean) ** 2 for t in temps) / n
    out = {
        "rows": n,
        "initial_K": rows[0][2],
        "mean_K": mean,
        "std_K": math.sqrt(var),
        "min_K": min(temps),
        "max_K": max(temps),
        "final_K": temps[-1],
    }
    for thr in (600.0, 2000.0):
        first = next((s for s, _, t in rows if t > thr), None)
        out["first_step_above_%dK" % int(thr)] = first
        out["first_time_above_%dK_ps" % int(thr)] = (
            None if first is None else first * DT_PS)
    if n >= 500:
        out["final500_mean_K"] = sum(temps[-500:]) / 500.0
    return out


def main():
    result = {"convention": {
        "k_boltzmann_ry": K_BOLTZMANN_RY, "ndof": NDOF, "dt_ps": DT_PS,
        "temperature": "T = 2*eka/(ndof*kB_Ry), eka from .eal column 2",
    }}

    # --- Leonardo af7ebfc ordinary no-CAIMD run (job 50127909) -------------
    leo = parse_eal(os.path.join(DATA, "leonardo_af7ebfc.eal"))
    result["leonardo_af7ebfc_ordinary"] = temp_stats(leo)

    frames = parse_cel(os.path.join(DATA, "leonardo_af7ebfc.cel"))
    v0 = abs(det3(frames[0][1]))
    vols = [abs(det3(m)) for _, m in frames]
    conds = [cond2(m) for _, m in frames]
    result["leonardo_af7ebfc_cell"] = {
        "frames": len(frames),
        "volume_first_bohr3": v0,
        "volume_final_bohr3": vols[-1],
        "volume_change_percent_final": 100.0 * (vols[-1] - v0) / v0,
        "volume_change_percent_max_abs": 100.0 * max(
            abs(v - v0) / v0 for v in vols),
        "condition_number_max": max(conds),
    }

    # --- Local ordinary Control A (controls ABC workflow) ------------------
    a = parse_eal(os.path.join(DATA, "control_A.eal"))
    result["control_A_ordinary_local"] = temp_stats(a)

    # --- Four-DOF tests 1 and 2 (active CAIMD + Nose) ----------------------
    for name in ("test1", "test2"):
        rows = parse_eal(os.path.join(DATA, "%s.eal" % name))
        result["four_dof_%s" % name] = temp_stats(rows)

    # --- Large-strain continuation .caimd stream ---------------------------
    c = parse_caimd(os.path.join(DATA, "cont_0p30_diamond8_continuous.caimd"))
    xx_err = [abs(a - n) for a, n in zip(c["act_xx"], c["nom_xx"])]
    result["continuation_0p30_caimd"] = {
        "rows": len(c["ionic_step"]),
        "events_final": c["events"][-1],
        "nominal_xx_final": c["nom_xx"][-1],
        "actual_xx_final": c["act_xx"][-1],
        "xx_abs_error_max": max(xx_err),
        "velocity_jump_max": max(c["velocity_jump_max"]),
        "last_actual_strains": {
            k: c["act_" + k][-1] for k in ("xx", "yy", "zz", "yz", "xz", "xy")},
        "last_temperature_K": c["temperature_K"][-1],
    }

    out = os.path.join(DATA, "recomputed_metrics.json")
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
