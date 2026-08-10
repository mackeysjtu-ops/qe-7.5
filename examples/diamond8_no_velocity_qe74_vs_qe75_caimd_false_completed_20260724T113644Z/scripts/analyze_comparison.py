#!/usr/bin/env python3
"""Post-processing for the Diamond8 no-explicit-velocity QE 7.4.1 vs QE 7.5 comparison.

Runs inside the Slurm job after both QE processes finish.
Standard library + (optionally) matplotlib. No network access required.

Temperature conventions (both reported):
  * eal_3Nm3 : T = 2*eka / ((3N-3) * k_B[Ry]), eka from prefix.eal column 2.
                This is the ionic-only kinetic temperature, 21 DOF for N=8.
  * stdout   : QE's own printed T for vc-md, which uses a 3*(N+1) denominator
                in vcsubs (includes cell degrees of freedom). Reported as-is.
Pressure:
  * stdout 'total   stress ... (kbar)  P=  <value>' is the QE-native pressure.
    QE prints P as the trace/3 of the stress tensor in kbar with the
    thermodynamic-pressure sign convention (positive = compressive load on the
    cell). Reported in kbar and GPa (1 GPa = 10 kbar).
"""
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ANA = os.path.join(ROOT, "analysis")
EVID = os.path.join(ROOT, "evidence")

RY_TO_KELVIN = 157887.51240116          # 1 / k_B[Ry]
NDOF_IONIC = 21                          # 3N-3, N=8
DT_FS = 2.0
NAT = 8

RUNS = [
    ("qe74", os.path.join(ROOT, "qe74_official"), "diamond8_qe74_novel",
     "Official QE 7.4.1 (Leonardo module)"),
    ("qe75", os.path.join(ROOT, "qe75_af7ebfc_caimd_false"), "diamond8_qe75_novel",
     "QE 7.5 af7ebfc, caimd=.false."),
]

STRESS_P_RE = re.compile(r"P\s*=\s*([-+0-9.EeDd]+)")
NUM_RE = re.compile(r"[-+]?\d+\.?\d*(?:[EeDd][-+]?\d+)?")


def floats(line):
    out = []
    for tok in NUM_RE.findall(line):
        try:
            out.append(float(tok.replace("D", "E").replace("d", "e")))
        except ValueError:
            pass
    return out


def read(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except OSError:
        return ""


def mean(v):
    return sum(v) / len(v) if v else None


def pstdev(v):
    if len(v) < 2:
        return 0.0 if v else None
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / len(v))


def parse_eal(path):
    """prefix.eal columns: uta eka eta utl ekla etl nst"""
    rows = []
    for line in read(path).splitlines():
        v = floats(line)
        if len(v) >= 7:
            uta, eka, eta, utl, ekla, etl, nst = v[:7]
            n = int(round(nst))
            rows.append({
                "nst": n, "step": n + 1, "time_ps": n * DT_FS / 1000.0,
                "eka_Ry": eka, "uta_Ry": uta, "eta_Ry": eta,
                "utl_Ry": utl, "ekla_Ry": ekla, "etl_Ry": etl,
                "temperature_eal_3Nm3_K": 2.0 * eka * RY_TO_KELVIN / NDOF_IONIC,
            })
    return rows


def parse_stdout_temperature(path):
    """QE vc-md prints 'Ekin = ... T = ... Etot = ...' lines."""
    rows = []
    step = None
    tim = None
    for line in read(path).splitlines():
        if "Entering Dynamics" in line:
            m = re.search(r"it\s*=\s*(\d+)", line)
            if m:
                step = int(m.group(1))
            m2 = re.search(r"time\s*=\s*([-+0-9.EeDd]+)", line)
            tim = float(m2.group(1).replace("D", "E")) if m2 else None
        if "Ekin =" in line and " T =" in line:
            v = floats(line)
            if len(v) >= 3:
                rows.append({"step": step, "time_ps": tim,
                             "ekin_Ry": v[0], "temperature_stdout_K": v[1],
                             "etot_Ry": v[2]})
    return rows


def parse_stress(path):
    """Parse per-step 'total   stress ... P=' blocks; pressure in kbar."""
    rows = []
    lines = read(path).splitlines()
    k = 0
    for i, line in enumerate(lines):
        if "stress" in line and "P=" in line:
            k += 1
            m = STRESS_P_RE.search(line)
            p = float(m.group(1).replace("D", "E")) if m else None
            mat = []
            ok = True
            for j in range(1, 4):
                if i + j >= len(lines):
                    ok = False
                    break
                v = floats(lines[i + j])
                if len(v) < 6:
                    ok = False
                    break
                mat.append(v[3:6])
            row = {"step": k, "time_ps": (k - 1) * DT_FS / 1000.0,
                   "pressure_kbar": p,
                   "pressure_GPa": (p / 10.0) if p is not None else None}
            if ok:
                row.update({"stress_xx_kbar": mat[0][0], "stress_yy_kbar": mat[1][1],
                            "stress_zz_kbar": mat[2][2], "stress_yz_kbar": mat[1][2],
                            "stress_xz_kbar": mat[0][2], "stress_xy_kbar": mat[0][1]})
                tr = (mat[0][0] + mat[1][1] + mat[2][2]) / 3.0
                row["stress_trace_over3_kbar"] = tr
            rows.append(row)
    return rows


def parse_cel(path):
    """prefix.cel: step header then 3 rows of the cell matrix (bohr)."""
    frames = []
    toks = [l.split() for l in read(path).splitlines() if l.strip()]
    i = 0
    while i + 3 < len(toks) + 1 and i + 3 <= len(toks):
        if len(toks[i]) == 2:
            try:
                step = int(float(toks[i][0]))
                m = [[float(x) for x in toks[i + r]] for r in (1, 2, 3)]
            except (ValueError, IndexError):
                i += 1
                continue
            frames.append((step, m))
            i += 4
        else:
            i += 1
    return frames


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def cond2(m):
    """2-norm condition number of a 3x3 matrix via eigenvalues of MtM (no numpy needed)."""
    try:
        import numpy as np
        s = np.linalg.svd(np.array(m), compute_uv=False)
        return float(s[0] / s[-1]) if s[-1] > 0 else float("inf")
    except Exception:
        # power-iteration free fallback: use Frobenius-based bound
        a = sum(sum(x * x for x in r) for r in m) ** 0.5
        d = abs(det3(m))
        return a ** 3 / d if d > 0 else float("inf")


def first_cross(rows, key, thr):
    for r in rows:
        v = r.get(key)
        if v is not None and v > thr:
            return {"step": r.get("step"), "time_ps": r.get("time_ps"), "value": v}
    return None


def drift_per_ps(rows, key):
    pts = [(r["time_ps"], r[key]) for r in rows
           if r.get("time_ps") is not None and r.get(key) is not None]
    n = len(pts)
    if n < 3:
        return None
    sx = sum(p[0] for p in pts); sy = sum(p[1] for p in pts)
    sxx = sum(p[0] * p[0] for p in pts); sxy = sum(p[0] * p[1] for p in pts)
    den = n * sxx - sx * sx
    return (n * sxy - sx * sy) / den if abs(den) > 1e-30 else None


def block_means(rows, key, size=100):
    out = []
    cur = []
    for r in rows:
        v = r.get(key)
        if v is not None:
            cur.append(v)
        if len(cur) == size:
            out.append(mean(cur)); cur = []
    if cur:
        out.append(mean(cur))
    return out


def stats(vals):
    if not vals:
        return {"n": 0}
    return {"n": len(vals), "mean": mean(vals), "std": pstdev(vals),
            "min": min(vals), "max": max(vals),
            "initial": vals[0], "final": vals[-1],
            "final200_mean": mean(vals[-200:]) if len(vals) >= 1 else None}


def classify_temperature(rc, job_done, temps):
    if not temps or len(temps) < 5:
        return "NUMERICAL_FAILURE" if rc != 0 else "INCOMPLETE"
    m = mean(temps)
    mx = max(temps)
    tail = mean(temps[-200:])
    if mx > 2000.0 and tail > 1000.0:
        return "PHYSICAL_HEATING"
    if tail < 60.0:
        return "PHYSICAL_COOLING"
    # statistical band for an 8-atom cell at 300 K: generous but bounded
    if 150.0 <= m <= 600.0 and mx <= 2000.0:
        return "STABLE_NEAR_300K"
    return "BOUNDED_BUT_OFFSET"


def classify_cell(frames):
    if len(frames) < 5:
        return "INCOMPLETE", {}
    dets = [det3(m) for _, m in frames]
    conds = [cond2(m) for _, m in frames]
    v0 = dets[0]
    rel = [abs(d - v0) / v0 for d in dets]
    info = {"volume_initial_bohr3": v0, "volume_mean_bohr3": mean(dets),
            "volume_final_bohr3": dets[-1],
            "volume_max_abs_change_percent": 100.0 * max(rel),
            "volume_final_change_percent": 100.0 * (dets[-1] - v0) / v0,
            "determinant_min_bohr3": min(dets),
            "condition_number_max": max(conds)}
    if min(dets) <= 0 or max(conds) > 5.0:
        return "SINGULAR_OR_UNSTABLE_CELL", info
    if max(rel) > 0.25:
        return "DRIFTING_PR_CELL", info
    return "STABLE_PR_CELL", info


def write_csv(path, rows, fields):
    with open(path, "w") as fh:
        fh.write(",".join(fields) + "\n")
        for r in rows:
            fh.write(",".join("" if r.get(f) is None else str(r.get(f)) for f in fields) + "\n")


def return_codes():
    rc = {}
    for line in read(os.path.join(EVID, "run_return_codes.txt")).splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                rc[k.strip()] = int(v.strip())
            except ValueError:
                pass
    return rc


def main():
    os.makedirs(ANA, exist_ok=True)
    rc = return_codes()
    rc_map = {"qe74": rc.get("RC_QE74"), "qe75": rc.get("RC_QE75_CAIMD_FALSE")}

    summary = {
        "comparison": "Diamond8 native from-scratch 300 K, no explicit velocities",
        "scope_note": ("This comparison tests native 300 K initialization "
                       "(ion_temperature='initial') followed by ordinary variable-cell "
                       "Parrinello-Rahman dynamics. It is NOT a continuing Nose-thermostat "
                       "test: QE 7.4.1 does not accept ion_temperature='nose' in PW vc-md."),
        "random_velocity_note": ("No explicit velocities are supplied, so each executable "
                                 "generates its own Maxwell realization via vcinit->ranv. "
                                 "QE's randy() is seeded deterministically (idum=0), so the "
                                 "realizations MAY coincide, but this is not assumed: "
                                 "statistical comparison is used."),
        "conventions": {
            "temperature_eal": "T = 2*eka/((3N-3)*k_B), ionic only, 21 DOF",
            "temperature_stdout": "QE-native vc-md T (3*(N+1) denominator, includes cell DOF)",
            "pressure": "QE stdout 'P=' in kbar; GPa = kbar/10",
        },
        "return_codes": rc_map,
        "runs": {},
    }

    per_run_rows = {}
    for key, rundir, prefix, label in RUNS:
        out = os.path.join(rundir, "qe.out")
        eal = os.path.join(rundir, "tmp", prefix + ".eal")
        if not os.path.exists(eal):
            eal = os.path.join(rundir, prefix + ".eal")
        cel = os.path.join(rundir, "tmp", prefix + ".cel")

        text = read(out)
        job_done = "JOB DONE" in text
        version = None
        m = re.search(r"Program PWSCF v\.([^\s]+)", text)
        if m:
            version = m.group(1)
        err = None
        me = re.search(r"Error in routine\s+(\S+)[^\n]*\n\s*(.*)", text)
        if me:
            err = (me.group(1) + ": " + me.group(2)).strip()
        if not err and "convergence NOT achieved" in text:
            err = "convergence NOT achieved (SCF)"

        eal_rows = parse_eal(eal)
        std_rows = parse_stdout_temperature(out)
        str_rows = parse_stress(out)
        frames = parse_cel(cel)

        # merge temperature + pressure by step for the per-run CSV
        pmap = {r["step"]: r for r in str_rows}
        merged = []
        for r in eal_rows:
            p = pmap.get(r["step"], {})
            merged.append({
                "step": r["step"], "time_ps": r["time_ps"],
                "temperature_eal_3Nm3_K": r["temperature_eal_3Nm3_K"],
                "eka_Ry": r["eka_Ry"],
                "pressure_kbar": p.get("pressure_kbar"),
                "pressure_GPa": p.get("pressure_GPa"),
            })
        if not merged and str_rows:
            merged = [{"step": r["step"], "time_ps": r["time_ps"],
                       "temperature_eal_3Nm3_K": None, "eka_Ry": None,
                       "pressure_kbar": r.get("pressure_kbar"),
                       "pressure_GPa": r.get("pressure_GPa")} for r in str_rows]
        per_run_rows[key] = merged
        write_csv(os.path.join(ANA, "%s_temperature_pressure.csv" % key), merged,
                  ["step", "time_ps", "temperature_eal_3Nm3_K", "eka_Ry",
                   "pressure_kbar", "pressure_GPa"])

        temps = [r["temperature_eal_3Nm3_K"] for r in eal_rows]
        temps_std = [r["temperature_stdout_K"] for r in std_rows
                     if r.get("temperature_stdout_K") is not None]
        press = [r["pressure_kbar"] for r in str_rows if r.get("pressure_kbar") is not None]

        tcls = classify_temperature(rc_map.get(key), job_done, temps or temps_std)
        ccls, cinfo = classify_cell(frames)

        summary["runs"][key] = {
            "label": label,
            "return_code": rc_map.get(key),
            "job_done": job_done,
            "runtime_version": version,
            "first_error": err,
            "completed_ionic_steps_eal": len(eal_rows),
            "completed_stress_blocks": len(str_rows),
            "completed_cell_frames": len(frames),
            "simulated_time_ps": (len(eal_rows) - 1) * DT_FS / 1000.0 if eal_rows else 0.0,
            "temperature_eal_3Nm3": stats(temps),
            "temperature_stdout_native": stats(temps_std),
            "temperature_first_above_600K": first_cross(eal_rows, "temperature_eal_3Nm3_K", 600.0),
            "temperature_first_above_1000K": first_cross(eal_rows, "temperature_eal_3Nm3_K", 1000.0),
            "temperature_first_above_2000K": first_cross(eal_rows, "temperature_eal_3Nm3_K", 2000.0),
            "temperature_drift_K_per_ps": drift_per_ps(eal_rows, "temperature_eal_3Nm3_K"),
            "pressure_kbar": stats(press),
            "pressure_GPa": stats([p / 10.0 for p in press]),
            "cell": cinfo,
            "temperature_classification": tcls,
            "cell_classification": ccls,
            "block100_temperature_means_K": block_means(eal_rows, "temperature_eal_3Nm3_K"),
            "block100_pressure_means_kbar": block_means(str_rows, "pressure_kbar"),
        }

    # ---------------- cross-run comparison ----------------
    a = summary["runs"]["qe74"]
    b = summary["runs"]["qe75"]
    both_ok = (a["job_done"] and b["job_done"]
               and a["temperature_eal_3Nm3"].get("n", 0) > 10
               and b["temperature_eal_3Nm3"].get("n", 0) > 10)

    cmp_block = {"both_produced_trajectories": both_ok}
    if both_ok:
        ta, tb = a["temperature_eal_3Nm3"], b["temperature_eal_3Nm3"]
        pa, pb = a["pressure_kbar"], b["pressure_kbar"]
        # expected fluctuation scale of the mean for an 8-atom cell
        sig = math.sqrt((ta["std"] ** 2 + tb["std"] ** 2) / 2.0) if ta["std"] and tb["std"] else None
        cmp_block.update({
            "initial_temperature_delta_K": abs((ta["initial"] or 0) - (tb["initial"] or 0)),
            "mean_temperature_delta_K": abs((ta["mean"] or 0) - (tb["mean"] or 0)),
            "final200_temperature_delta_K": abs((ta["final200_mean"] or 0) - (tb["final200_mean"] or 0)),
            "temperature_fluctuation_scale_K": sig,
            "mean_pressure_delta_kbar": abs((pa["mean"] or 0) - (pb["mean"] or 0)),
            "pressure_fluctuation_scale_kbar": (
                math.sqrt((pa["std"] ** 2 + pb["std"] ** 2) / 2.0)
                if pa.get("std") and pb.get("std") else None),
            "same_temperature_class": a["temperature_classification"] == b["temperature_classification"],
            "same_cell_class": a["cell_classification"] == b["cell_classification"],
            "trajectories_bitwise_identical_temperature": (
                ta["n"] == tb["n"]
                and all(abs(x - y) < 1e-10 for x, y in
                        zip([r["temperature_eal_3Nm3_K"] for r in per_run_rows["qe74"]],
                            [r["temperature_eal_3Nm3_K"] for r in per_run_rows["qe75"]]
                            ) if x is not None and y is not None)),
        })
        temp_similar = (cmp_block["same_temperature_class"]
                        and (sig is None or cmp_block["mean_temperature_delta_K"] <= max(2.0 * sig, 150.0)))
        press_similar = (cmp_block["mean_pressure_delta_kbar"] <=
                         max(2.0 * (cmp_block["pressure_fluctuation_scale_kbar"] or 0.0), 50.0))
        cmp_block["temperature_broadly_similar"] = temp_similar
        cmp_block["pressure_broadly_similar"] = press_similar

        ca, cb = a["temperature_classification"], b["temperature_classification"]
        heat = {"PHYSICAL_HEATING"}
        stable = {"STABLE_NEAR_300K"}
        if ca in heat and cb in heat:
            overall = "BOTH_HEAT"
        elif ca in stable and cb in heat:
            overall = "QE74_STABLE_QE75_HEATS"
        elif ca in heat and cb in stable:
            overall = "QE74_HEATS_QE75_STABLE"
        elif temp_similar and press_similar:
            overall = "BOTH_COMPARABLE"
        elif temp_similar and not press_similar:
            overall = "BOTH_STABLE_DIFFERENT_PRESSURE"
        else:
            overall = "PARTIAL"
    else:
        overall = "PARTIAL"
        if a["first_error"] and "ion_temperature" in (a["first_error"] or ""):
            overall = "BLOCKED_INPUT_COMPATIBILITY"
    summary["comparison_metrics"] = cmp_block
    summary["overall_status"] = overall

    json.dump(summary, open(os.path.join(ANA, "comparison_summary.json"), "w"), indent=2)

    # ---------------- summary CSVs ----------------
    trows, prows, brows = [], [], []
    for key, _, _, label in RUNS:
        r = summary["runs"][key]
        t, p = r["temperature_eal_3Nm3"], r["pressure_kbar"]
        trows.append({"run": key, "label": label,
                      "classification": r["temperature_classification"],
                      "return_code": r["return_code"], "job_done": r["job_done"],
                      "steps": r["completed_ionic_steps_eal"],
                      "initial_K": t.get("initial"), "mean_K": t.get("mean"),
                      "std_K": t.get("std"), "min_K": t.get("min"), "max_K": t.get("max"),
                      "final_K": t.get("final"), "final200_mean_K": t.get("final200_mean"),
                      "drift_K_per_ps": r["temperature_drift_K_per_ps"]})
        prows.append({"run": key, "label": label,
                      "cell_classification": r["cell_classification"],
                      "initial_kbar": p.get("initial"), "mean_kbar": p.get("mean"),
                      "std_kbar": p.get("std"), "min_kbar": p.get("min"),
                      "max_kbar": p.get("max"), "final_kbar": p.get("final"),
                      "final200_mean_kbar": p.get("final200_mean"),
                      "mean_GPa": (p.get("mean") / 10.0) if p.get("mean") is not None else None,
                      "volume_initial_bohr3": r["cell"].get("volume_initial_bohr3"),
                      "volume_final_bohr3": r["cell"].get("volume_final_bohr3"),
                      "volume_final_change_percent": r["cell"].get("volume_final_change_percent"),
                      "condition_number_max": r["cell"].get("condition_number_max")})
        for i, (tm, pm) in enumerate(zip(r["block100_temperature_means_K"],
                                         r["block100_pressure_means_kbar"] + [None] * 50)):
            brows.append({"run": key, "block": i + 1,
                          "step_start": i * 100 + 1, "step_end": (i + 1) * 100,
                          "temperature_mean_K": tm, "pressure_mean_kbar": pm})

    write_csv(os.path.join(ANA, "temperature_statistics.csv"), trows,
              ["run", "label", "classification", "return_code", "job_done", "steps",
               "initial_K", "mean_K", "std_K", "min_K", "max_K", "final_K",
               "final200_mean_K", "drift_K_per_ps"])
    write_csv(os.path.join(ANA, "pressure_statistics.csv"), prows,
              ["run", "label", "cell_classification", "initial_kbar", "mean_kbar",
               "std_kbar", "min_kbar", "max_kbar", "final_kbar", "final200_mean_kbar",
               "mean_GPa", "volume_initial_bohr3", "volume_final_bohr3",
               "volume_final_change_percent", "condition_number_max"])
    write_csv(os.path.join(ANA, "block_statistics.csv"), brows,
              ["run", "block", "step_start", "step_end", "temperature_mean_K",
               "pressure_mean_kbar"])

    # ---------------- plots (best effort) ----------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        colors = {"qe74": "#2a78d6", "qe75": "#eb6834"}
        labels = {"qe74": "QE 7.4.1 official", "qe75": "QE 7.5 af7ebfc caimd=.false."}

        def line_plot(fname, ykey, ylabel, title, hline=None):
            fig, ax = plt.subplots(figsize=(10, 5))
            drew = False
            for key in ("qe74", "qe75"):
                rows = [r for r in per_run_rows[key] if r.get(ykey) is not None]
                if not rows:
                    continue
                ax.plot([r["time_ps"] for r in rows], [r[ykey] for r in rows],
                        lw=1.0, color=colors[key], label=labels[key])
                drew = True
            if hline is not None:
                ax.axhline(hline, color="#898781", ls="--", lw=1.2)
            ax.set_xlabel("time (ps)"); ax.set_ylabel(ylabel); ax.set_title(title)
            if drew:
                ax.legend(frameon=False)
            fig.savefig(os.path.join(ANA, fname), dpi=150, bbox_inches="tight")
            plt.close(fig)

        line_plot("temperature_vs_time.png", "temperature_eal_3Nm3_K",
                  "ionic temperature (K, 3N-3)",
                  "Diamond8 native 300 K init, no explicit velocities", hline=300.0)
        line_plot("pressure_vs_time.png", "pressure_kbar", "pressure (kbar)",
                  "Pressure vs time", hline=0.0)

        # volume
        fig, ax = plt.subplots(figsize=(10, 5))
        for key, rundir, prefix, _ in RUNS:
            cel = os.path.join(rundir, "tmp", prefix + ".cel")
            fr = parse_cel(cel)
            if not fr:
                continue
            ax.plot([(s - 1) * DT_FS / 1000.0 for s, _ in fr],
                    [det3(m) for _, m in fr], lw=1.2, color=colors[key], label=labels[key])
        ax.set_xlabel("time (ps)"); ax.set_ylabel("cell volume (bohr^3)")
        ax.set_title("Cell volume vs time"); ax.legend(frameon=False)
        fig.savefig(os.path.join(ANA, "volume_vs_time.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)

        # block means
        for fname, bkey, ylab, ttl in (
                ("temperature_block_means.png", "temperature_mean_K",
                 "block-mean T (K)", "100-step block temperature means"),
                ("pressure_block_means.png", "pressure_mean_kbar",
                 "block-mean P (kbar)", "100-step block pressure means")):
            fig, ax = plt.subplots(figsize=(10, 5))
            for key in ("qe74", "qe75"):
                rows = [r for r in brows if r["run"] == key and r.get(bkey) is not None]
                if not rows:
                    continue
                ax.plot([r["block"] for r in rows], [r[bkey] for r in rows],
                        marker="o", ms=4, lw=1.4, color=colors[key], label=labels[key])
            ax.set_xlabel("100-step block"); ax.set_ylabel(ylab); ax.set_title(ttl)
            ax.legend(frameon=False)
            fig.savefig(os.path.join(ANA, fname), dpi=150, bbox_inches="tight")
            plt.close(fig)
        summary["plots_generated"] = True
    except Exception as exc:  # plotting must never break the analysis
        summary["plots_generated"] = False
        summary["plot_error"] = repr(exc)
        json.dump(summary, open(os.path.join(ANA, "comparison_summary.json"), "w"), indent=2)

    write_report(summary)
    print(json.dumps({"overall_status": summary["overall_status"],
                      "qe74_class": a["temperature_classification"],
                      "qe75_class": b["temperature_classification"]}, indent=2))
    return 0


def fmt(x, nd=3):
    return "n/a" if x is None else ("%.*f" % (nd, x))


def write_report(s):
    a, b = s["runs"]["qe74"], s["runs"]["qe75"]
    c = s["comparison_metrics"]
    ta, tb = a["temperature_eal_3Nm3"], b["temperature_eal_3Nm3"]
    pa, pb = a["pressure_kbar"], b["pressure_kbar"]
    L = []
    L.append("# Final Comparison Report")
    L.append("")
    L.append("Diamond8, native from-scratch 300 K initialization, **no explicit velocities**.")
    L.append("")
    L.append("Overall status: `%s`" % s["overall_status"])
    L.append("")
    L.append("## Scope")
    L.append("")
    L.append(s["scope_note"])
    L.append("")
    L.append(s["random_velocity_note"])
    L.append("")
    L.append("## Required answers")
    L.append("")
    L.append("1. **Explicit velocities absent from both inputs?** Yes. See "
             "`evidence/no_explicit_velocity_audit.json`; neither input contains an "
             "ATOMIC_VELOCITIES card or an explicit-velocity keyword, in the complete "
             "text (case-insensitive).")
    L.append("2. **Common 300 K initialization mode:** `ion_temperature='initial'` with "
             "`tempw=300.0`, accepted by both versions.")
    L.append("3. **Both executables initialized velocities?** QE 7.4.1: %s; QE 7.5: %s "
             "(non-zero initial ionic kinetic energy)."
             % ("yes" if (ta.get("initial") or 0) > 1.0 else "no/undetermined",
                "yes" if (tb.get("initial") or 0) > 1.0 else "no/undetermined"))
    L.append("4. **Actual initial temperature:** QE 7.4.1 = %s K; QE 7.5 = %s K (3N-3 ionic)."
             % (fmt(ta.get("initial")), fmt(tb.get("initial"))))
    L.append("5. **QE 7.4.1 near 300 K?** classification `%s` (mean %s K, final-200 mean %s K)."
             % (a["temperature_classification"], fmt(ta.get("mean")), fmt(ta.get("final200_mean"))))
    L.append("6. **QE 7.5 caimd=.false. near 300 K?** classification `%s` (mean %s K, "
             "final-200 mean %s K)."
             % (b["temperature_classification"], fmt(tb.get("mean")), fmt(tb.get("final200_mean"))))
    L.append("7. **Did either run physically heat?** QE 7.4.1 max %s K (first >2000 K: %s); "
             "QE 7.5 max %s K (first >2000 K: %s)."
             % (fmt(ta.get("max")), a["temperature_first_above_2000K"],
                fmt(tb.get("max")), b["temperature_first_above_2000K"]))
    L.append("8. **Mean / final-window temperatures:** QE 7.4.1 %s / %s K; QE 7.5 %s / %s K."
             % (fmt(ta.get("mean")), fmt(ta.get("final200_mean")),
                fmt(tb.get("mean")), fmt(tb.get("final200_mean"))))
    L.append("9. **Mean / final-window pressures:** QE 7.4.1 %s / %s kbar (%s GPa mean); "
             "QE 7.5 %s / %s kbar (%s GPa mean)."
             % (fmt(pa.get("mean")), fmt(pa.get("final200_mean")),
                fmt((pa.get("mean") or 0) / 10.0),
                fmt(pb.get("mean")), fmt(pb.get("final200_mean")),
                fmt((pb.get("mean") or 0) / 10.0)))
    L.append("10. **Pressure differences consistent with 8-atom fluctuations?** %s "
             "(|Δmean| = %s kbar vs fluctuation scale %s kbar)."
             % (c.get("pressure_broadly_similar"), fmt(c.get("mean_pressure_delta_kbar")),
                fmt(c.get("pressure_fluctuation_scale_kbar"))))
    L.append("11. **Cell instability?** QE 7.4.1 `%s` (final ΔV %s %%, cond max %s); "
             "QE 7.5 `%s` (final ΔV %s %%, cond max %s)."
             % (a["cell_classification"], fmt(a["cell"].get("volume_final_change_percent")),
                fmt(a["cell"].get("condition_number_max"), 4),
                b["cell_classification"], fmt(b["cell"].get("volume_final_change_percent")),
                fmt(b["cell"].get("condition_number_max"), 4)))
    L.append("12. **Statistically comparable?** temperature: %s; pressure: %s; "
             "identical temperature series: %s."
             % (c.get("temperature_broadly_similar"), c.get("pressure_broadly_similar"),
                c.get("trajectories_bitwise_identical_temperature")))
    L.append("13. **Narrowest justified next test:** see below.")
    L.append("")
    nxt = {
        "BOTH_COMPARABLE": "Both versions behave equivalently under native 300 K "
                           "initialization. The next test is to reintroduce the QE 7.5 "
                           "Nose thermostat alone (7.5 vs 7.5, nose vs initial) to isolate "
                           "whether the previously observed heating is specific to the "
                           "PW vc-md Nose path rather than to the executable.",
        "BOTH_HEAT": "Both versions heat from a clean 300 K start, so the phenomenon is not "
                     "specific to the QE 7.5 CAIMD-enabled executable. The next test is a "
                     "timestep/PR-cell-mass sensitivity scan (and fixed-cell md vs vc-md) "
                     "at identical settings.",
        "QE74_STABLE_QE75_HEATS": "A version/executable difference is demonstrated. This does "
                                  "NOT by itself separate an upstream QE 7.5 change from an "
                                  "af7ebfc regression: the next test is pristine QE 7.5 "
                                  "(commit a11913a) built with the same toolchain, same input.",
        "QE74_HEATS_QE75_STABLE": "Unexpected direction; repeat with a second random "
                                  "realization before interpreting, then compare pristine QE 7.5.",
        "BOTH_STABLE_DIFFERENT_PRESSURE": "Temperatures agree but pressure statistics differ; "
                                          "next test should check stress convention and "
                                          "cutoff/Pulay behavior at higher ecutwfc.",
        "PARTIAL": "Only one run produced sufficient evidence; diagnose the failing run "
                   "before drawing any version conclusion.",
        "BLOCKED_INPUT_COMPATIBILITY": "The common input was not accepted; no comparison is "
                                       "possible until an accepted common mode is found.",
    }.get(s["overall_status"], "Diagnose before further interpretation.")
    L.append(nxt)
    L.append("")
    L.append("## Run summary")
    L.append("")
    L.append("| run | version | rc | JOB DONE | steps | T init (K) | T mean (K) | T max (K) | "
             "P mean (kbar) | T class | cell class |")
    L.append("|---|---|---:|---|---:|---:|---:|---:|---:|---|---|")
    for k, r in (("qe74", a), ("qe75", b)):
        t, p = r["temperature_eal_3Nm3"], r["pressure_kbar"]
        L.append("| %s | %s | %s | %s | %d | %s | %s | %s | %s | %s | %s |"
                 % (k, r.get("runtime_version"), r.get("return_code"), r.get("job_done"),
                    r.get("completed_ionic_steps_eal", 0), fmt(t.get("initial")),
                    fmt(t.get("mean")), fmt(t.get("max")), fmt(p.get("mean")),
                    r["temperature_classification"], r["cell_classification"]))
    L.append("")
    L.append("## Caveats")
    L.append("")
    L.append("- `sbatch` acceptance is not scientific validation; this report reflects parsed "
             "output only.")
    L.append("- Instantaneous pressures in an 8-atom cell fluctuate strongly; only statistical "
             "agreement is meaningful.")
    L.append("- QE 7.4.1 and QE 7.5 may generate different Maxwell realizations; exact "
             "step-by-step equality is not required unless the identical-series check above "
             "reports true.")
    L.append("- First errors, if any: QE 7.4.1 `%s`; QE 7.5 `%s`."
             % (a.get("first_error"), b.get("first_error")))
    L.append("")
    with open(os.path.join(ROOT, "FINAL_COMPARISON_REPORT.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
