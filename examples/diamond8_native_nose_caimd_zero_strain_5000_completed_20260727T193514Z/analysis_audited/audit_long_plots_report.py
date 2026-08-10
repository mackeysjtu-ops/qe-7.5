#!/usr/bin/env python3
"""Plots + audited report for job 50434898 (three 5000-step native-init Nose + CAIMD)."""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MIRROR = sys.argv[1]
ANA = os.path.join(MIRROR, "analysis_audited")
EVD = os.path.join(MIRROR, "evidence_audited")
S = json.load(open(os.path.join(ANA, "comparison_summary_audited.json")))
Z = np.load(os.path.join(EVD, "_arrays.npz"))
NAT = 8
CASES = ["A", "B", "C"]
COL = {"A": "#2a78d6", "B": "#1baf7a", "C": "#eb6834"}
LB = {"A": "CASE A caimd=.false.", "B": "CASE B no-op", "C": "CASE C active zero-strain"}
SRC = ("Source: job 50434898 mirror; T=2*eka/((3N-3)k_B) from .eal; P from .p; cell from .cel; "
       "strain/residual from .caimd. ecutwfc=25 Ry -> pressures not converged.")
plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.color": "#e1e0d9", "grid.linewidth": .8,
                     "axes.axisbelow": True, "axes.unicode_minus": False, "savefig.dpi": 150,
                     "figure.facecolor": "white", "axes.facecolor": "white"})


def roll(a, w):
    if len(a) < w:
        return np.full(len(a), np.nan)
    return np.concatenate([np.full(w - 1, np.nan), np.convolve(a, np.ones(w) / w, "valid")])


def blk(a, w):
    return np.array([a[i:i + w].mean() for i in range(0, len(a), w) if len(a[i:i + w])])


def fin(ax, title, xl, yl, hline=None, note=SRC, legend=True, loc="best"):
    if hline is not None:
        ax.axhline(hline, color="#898781", ls="--", lw=1.2)
    ax.set_title(title, fontsize=11.5); ax.set_xlabel(xl); ax.set_ylabel(yl)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if legend:
        ax.legend(frameon=False, fontsize=9, loc=loc)
    ax.figure.text(0.01, -0.02, note, fontsize=7.2, color="#52514e")


def savef(fig, name):
    fig.savefig(os.path.join(ANA, name), bbox_inches="tight"); plt.close(fig); print("wrote", name)


t = {k: Z[f"{k}_t_ps"] for k in CASES}
T = {k: Z[f"{k}_T"] for k in CASES}
eka = {k: Z[f"{k}_eka"] for k in CASES}
P = {k: Z[f"{k}_p_kbar"] for k in CASES}
E = {k: Z[f"{k}_edyn"] for k in CASES}
V = {k: Z[f"{k}_vol"] for k in CASES}
CT = {k: Z[f"{k}_celt"] for k in CASES}

# 1 instantaneous
fig, ax = plt.subplots(figsize=(11, 5))
for k in CASES:
    ax.plot(t[k], T[k], lw=.4, color=COL[k], alpha=.8, label=LB[k])
fin(ax, "Instantaneous ionic temperature, 10 ps (native init + Nose)", "time (ps)", "T (K, 3N-3)", hline=300)
savef(fig, "temperature_instantaneous_all_cases.png")

# 2 rolling 500
fig, ax = plt.subplots(figsize=(11, 5))
for k in CASES:
    ax.plot(t[k], roll(T[k], 500), lw=1.7, color=COL[k], label=LB[k])
fin(ax, "Temperature, 500-step rolling mean (heating onset ~3-4 ps in A/B)", "time (ps)", "rolling mean T (K)", hline=300)
savef(fig, "temperature_rolling500_all_cases.png")

# 3 block1000
fig, ax = plt.subplots(figsize=(11, 5))
for k in CASES:
    b = blk(T[k], 1000)
    ax.plot(np.arange(1, len(b) + 1) * 2 - 1, b, "o-", ms=7, lw=1.6, color=COL[k], label=LB[k])
fin(ax, "Temperature, 1000-step (2 ps) block means", "block centre (ps)", "block-mean T (K)", hline=300)
savef(fig, "temperature_block1000_all_cases.png")

# 4 last 5 ps
fig, ax = plt.subplots(figsize=(11, 5))
for k in CASES:
    m = t[k] >= (t[k][-1] - 5.0)
    ax.plot(t[k][m], T[k][m], lw=.5, color=COL[k], alpha=.85, label=LB[k])
fin(ax, "Temperature, last 5 ps", "time (ps)", "T (K, 3N-3)", hline=300)
savef(fig, "temperature_last5ps_all_cases.png")

# 5 block uncertainty vs 300 K
fig, ax = plt.subplots(figsize=(9.5, 5.5))
for i, k in enumerate(CASES):
    u = S["runs"][k]["block_uncertainty_1000"]
    bm, sem = u.get("block_mean"), u.get("sem")
    if bm is not None:
        ax.errorbar([i], [bm], yerr=[sem or 0], fmt="o", ms=10, color=COL[k], capsize=6, label=LB[k])
ax.axhline(300, color="#898781", ls="--", lw=1.3); ax.axhspan(270, 330, color="#0ca30c", alpha=.08)
ax.text(2.3, 300, "300 K", color="#52514e", fontsize=9)
ax.set_xticks(range(3)); ax.set_xticklabels(CASES)
fin(ax, "Final-half block mean +/- SEM vs 300 K target band", "case", "late-time mean T (K)", legend=True, loc="upper center")
savef(fig, "temperature_block_uncertainty_vs_300K.png")

# 6 definition crosscheck (CASE A eka -> T with 3 divisors)
fig, ax = plt.subplots(figsize=(11, 5))
for div, lab, c in ((21, "3N-3 = 21 (= ndega, authoritative)", "#2a78d6"),
                    (24, "3N = 24", "#eb6834"), (27, "3(N+1) = 27", "#1baf7a")):
    ax.plot(t["A"], roll(2 * eka["A"] * 157887.51240116 / div, 500), lw=1.5, color=c, label=lab)
fin(ax, "CASE A temperature under different DOF divisors (500-step rolling)", "time (ps)", "T (K)", hline=300,
    note=SRC + " ndega=get_ndof()=3N-3=21; no divisor choice turns a below/above-target eka into 300 K.")
savef(fig, "temperature_definition_crosscheck.png")

# 7 ionic KE
fig, ax = plt.subplots(figsize=(11, 5))
for k in CASES:
    ax.plot(t[k], eka[k], lw=.4, color=COL[k], alpha=.8, label=LB[k])
fin(ax, "Ionic kinetic energy", "time (ps)", "eka (Ry)")
savef(fig, "ionic_kinetic_energy_all_cases.png")

# 8 pressure
fig, ax = plt.subplots(figsize=(11, 5))
for k in CASES:
    ax.plot(t[k], P[k], lw=.4, color=COL[k], alpha=.8, label=LB[k])
fin(ax, "MD pressure (NOT converged; diagnostic only)", "time (ps)", "P (kbar)", hline=0)
savef(fig, "pressure_all_cases.png")

# 9 volume
fig, ax = plt.subplots(figsize=(11, 5))
for k in CASES:
    ax.plot(CT[k], V[k], lw=1.3, color=COL[k], label=LB[k])
fin(ax, "Cell volume: A/B contract ~9%, C held rigid by zero-strain CAIMD", "time (ps)", "volume (bohr^3)")
savef(fig, "volume_all_cases.png")

# 10 energy drift
fig, ax = plt.subplots(figsize=(11, 5))
for k in CASES:
    ax.plot(t[k], (E[k] - E[k][0]) * 13.605693122994 / NAT, lw=1.2, color=COL[k], label=LB[k])
fin(ax, "Partial extended-energy drift (edyn=ekint+ut; excludes Nose chain)", "time (ps)",
    "drift per atom (eV/atom)", hline=0, note=SRC + " Nose chain energy not written by QE -> partial.")
savef(fig, "energy_drift_all_cases.png")

# CASE C strain plots
mir_C = os.path.join(MIRROR, "case_C_caimd_zero_strain", "tmp", "diamond8_case_C_caimd_zero.caimd")
import re
NUM = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")
rows = []
for line in open(mir_C, errors="replace"):
    if line.lstrip().startswith("#"):
        continue
    parts = line.split()
    if len(parts) == 26:
        rows.append([float(x.replace("D", "E")) for x in (parts[:3] + ["0"] + parts[4:])])
arr = np.array(rows)
ct = arr[:, 1]
comps = ["xx", "yy", "zz", "yz", "xz", "xy"]
act = arr[:, 10:16]
nom = arr[:, 4:10]
vjmax = arr[:, 24]

fig, ax = plt.subplots(figsize=(11, 5))
for i, (cmp_, c) in enumerate(zip(comps, ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"))):
    ax.plot(ct, act[:, i], lw=1.0, color=c, label=f"actual {cmp_}")
fin(ax, "CASE C: all six actual strain components (target = 0)", "time (ps)", "actual strain (tensor)")
savef(fig, "case_C_all_six_strains.png")

fig, ax = plt.subplots(figsize=(11, 5))
err = np.max(np.abs(nom - act), axis=1)
pos = err > 0
if pos.any():
    ax.plot(ct[pos], err[pos], lw=1.2, color="#eb6834")
    ax.set_yscale("log")
    from matplotlib.ticker import FuncFormatter, LogLocator
    ax.yaxis.set_major_locator(LogLocator(base=10.0))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: "1e%d" % int(round(np.log10(v))) if v > 0 else ""))
else:
    ax.plot(ct, err, lw=1.2, color="#eb6834")
fin(ax, "CASE C: max |nominal - actual| strain error (all six; roundoff scale)", "time (ps)",
    "max strain error", legend=False)
savef(fig, "case_C_strain_tracking_error.png")

fig, ax = plt.subplots(figsize=(11, 5))
pos = vjmax > 0
if pos.any():
    ax.plot(ct[pos], vjmax[pos], lw=1.2, color="#eb6834")
    ax.set_yscale("log")
    from matplotlib.ticker import FuncFormatter, LogLocator
    ax.yaxis.set_major_locator(LogLocator(base=10.0))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: "1e%d" % int(round(np.log10(v))) if v > 0 else ""))
else:
    ax.plot(ct, vjmax, lw=1.2, color="#eb6834")
fin(ax, "CASE C: velocity-transform residual (velocity_jump_max)", "time (ps)",
    "velocity_jump_max (a.u.)", legend=False)
savef(fig, "case_C_velocity_transform_residual.png")

# A vs B no-op
fig, ax = plt.subplots(figsize=(11, 5))
m = min(len(T["A"]), len(T["B"]))
ax.plot(t["A"][:m], np.abs(T["A"][:m] - T["B"][:m]), lw=.6, color="#2a78d6")
ax.set_yscale("symlog", linthresh=1e-6)
fin(ax, "CASE A vs B temperature difference (bit-identical to step ~701, then chaotic)", "time (ps)",
    "|T_A - T_B| (K)", legend=False)
savef(fig, "case_A_vs_B_noop_difference.png")

# ---------------------------------------------------------------- report ----
R = S["runs"]
rr = S["temperature_reasonableness"]
ab = json.load(open(os.path.join(ANA, "case_A_vs_B_noop_audit.json")))
cC = json.load(open(os.path.join(ANA, "case_C_zero_strain_audit.json")))
nose = json.load(open(os.path.join(ANA, "nose_activation_audit.json")))
init = json.load(open(os.path.join(ANA, "initial_state_comparison.json")))


def f(x, n=1):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{n}f}"


L = []
L.append("# Audited Final Report: Diamond8 5000-step native Nose + CAIMD zero-strain")
L.append("")
L.append("Job **50434898**, QE 7.5 `af7ebfc`. Three 10 ps native-initialization Nose trajectories, "
         "no explicit velocities: CASE A `caimd=.false.`, CASE B enabled no-op, CASE C active "
         "all-component zero-strain hold.")
L.append("")
L.append(f"**Overall status: `{S['overall_status']}`**")
L.append("")
L.append("| case | T class | mean K | final-5ps K | final-2.5ps K | max K | late-half block K +/- SEM | >600K frac |")
L.append("|---|---|---:|---:|---:|---:|---:|---:|")
for k in CASES:
    r = R[k]; st = r["temperature"]; u = r["block_uncertainty_1000"]
    L.append(f"| {k} | {r['temperature_classification']} | {f(st['mean'])} | {f(r['final_5ps_mean_K'])} | "
             f"{f(r['final_2p5ps_mean_K'])} | {f(st['max'])} | {f(u.get('block_mean'))} +/- {f(u.get('sem'))} | "
             f"{r['fraction_above']['600K']:.3f} |")
L.append("")
L.append("## Headline")
L.append("")
L.append(S["headline"])
L.append("")
L.append("This overturns the earlier 2 ps picture (job 50167241) that both control and Nose stayed "
         "bounded near 170 K: the heating onset is at ~3-4 ps, beyond that window. The three 1000-step "
         "block means for CASE A are "
         f"{[round(x) for x in R['A']['block1000_means_K']]} K and for CASE C "
         f"{[round(x) for x in R['C']['block1000_means_K']]} K.")
L.append("")
L.append("## The 300 K question (three separate answers)")
L.append("")
L.append("**A. Are instantaneous fluctuations reasonable for 8 atoms?** Yes. Large step-to-step "
         "swings are expected for a 3N-3 = 21 DOF system; the instantaneous standard deviation is not "
         "itself a defect.")
L.append(f"**B. Is the late-time mean statistically compatible with 300 K?** **No, in every case.** "
         f"CASE A final-half block mean {f(R['A']['block_uncertainty_1000'].get('block_mean'))} +/- "
         f"{f(R['A']['block_uncertainty_1000'].get('sem'))} K, CASE B "
         f"{f(R['B']['block_uncertainty_1000'].get('block_mean'))} +/- "
         f"{f(R['B']['block_uncertainty_1000'].get('sem'))} K, CASE C "
         f"{f(R['C']['block_uncertainty_1000'].get('block_mean'))} +/- "
         f"{f(R['C']['block_uncertainty_1000'].get('sem'))} K. None overlaps 300 K within 2 SEM.")
L.append("**C. Is the thermostat scientifically validated for future CAIMD loading?** No. The ionic "
         "Nose is active every step with the correct target but does not maintain 300 K.")
L.append("")
L.append("**Plain-language verdict (温度是否合理):** the fluctuations are expected for 8 atoms, but "
         "the ionic subsystem is not held at the requested 300 K in any case. CASE A/B run away to "
         f"~{f(R['A']['final_5ps_mean_K'],0)} K driven by the free Parrinello-Rahman cell (which "
         "contracts ~9%); CASE C, with the cell frozen by zero-strain CAIMD, is far milder "
         f"(~{f(R['C']['final_5ps_mean_K'],0)} K) but still above 300 K and rising. The trajectories "
         "are numerically stable but thermodynamically at the wrong temperature.")
L.append("")
L.append("## Required answers")
L.append("")
L.append(f"1. **All three complete 5000 steps?** Yes - " +
         ", ".join(f"{k}={S['run_completeness'][k]['completed_ionic_steps']}({S['run_completeness'][k]['classification']})" for k in CASES) + ".")
L.append(f"2. **Explicit velocities absent?** {S['input_contract']['explicit_velocities_absent']} - no "
         "ATOMIC_VELOCITIES/ion_velocities/from_input in any input's complete text; all from_scratch.")
L.append(f"3. **Initial velocity states identical?** Native realization is the same (first-step eka = "
         f"0.019951 Ry in all three). A vs B: `{init['pairs']['A_vs_B']['classification']}` "
         f"(bit-identical to step {init['pairs']['A_vs_B']['first_eka_divergence_step']}). "
         f"A/B vs C: `{init['pairs']['A_vs_C']['classification']}` - C diverges at step "
         f"{init['pairs']['A_vs_C']['first_eka_divergence_step']} because its active projection acts "
         f"from the start.")
L.append("4. **Authoritative temperature quantity:** instantaneous ionic kinetic temperature "
         "`T = 2*eka/((3N-3)k_B)` with 3N-3 = 21, `eka` from `.eal` column 2. This matches "
         "`ndega = get_ndof() = 21`, the Nose target DOF. The `.tv`/stdout temperature is a running "
         "average and is NOT used. Terminal direct `.md` velocity reconstruction confirms `.eal` "
         f"(A within {f(abs((R['A']['direct_terminal_T_K'] or 0)-R['A']['temperature']['final']),2)} K).")
L.append(f"5. **Full-run means:** A {f(R['A']['temperature']['mean'])}, B {f(R['B']['temperature']['mean'])}, "
         f"C {f(R['C']['temperature']['mean'])} K.")
L.append(f"6. **Final-5 ps means:** A {f(R['A']['final_5ps_mean_K'])}, B {f(R['B']['final_5ps_mean_K'])}, "
         f"C {f(R['C']['final_5ps_mean_K'])} K.")
L.append(f"7. **Final-2.5 ps means:** A {f(R['A']['final_2p5ps_mean_K'])}, B {f(R['B']['final_2p5ps_mean_K'])}, "
         f"C {f(R['C']['final_2p5ps_mean_K'])} K.")
L.append(f"8. **1000-step block means (K):** A {[round(x) for x in R['A']['block1000_means_K']]}; "
         f"B {[round(x) for x in R['B']['block1000_means_K']]}; C {[round(x) for x in R['C']['block1000_means_K']]}.")
L.append(f"9. **Block-based uncertainties (final-half, SEM):** A {f(R['A']['block_uncertainty_1000'].get('sem'))}, "
         f"B {f(R['B']['block_uncertainty_1000'].get('sem'))}, C {f(R['C']['block_uncertainty_1000'].get('sem'))} K.")
L.append("10. **Late-time means compatible with 300 K?** No for all three (see answer B above).")
L.append(f"11. **Instantaneous fluctuations reasonable for 8 atoms?** Yes (std ~"
         f"{f(rr['A']['instantaneous_std_K'],0)}-{f(rr['C']['instantaneous_std_K'],0)} K); large swings "
         f"are expected and are not the problem.")
L.append("12. **Factor-of-two / DOF issue?** No. `ndega = get_ndof() = 3N-3 = 21` equals the "
         "measurement divisor, so a working thermostat would read ~300 K. No divisor choice (21/24/27) "
         "reconciles the observed eka with 300 K - the temperatures are genuinely off-target, not a "
         "definitional artifact. (This is a real change from the earlier ~170 K discussion: at 10 ps "
         "the system is HOT, not cold.)")
L.append("13. **Is Nose applied every step?** Yes. `tnosep=.true.` throughout; the `vcmove` Nose "
         "friction term acts unconditionally each step; the velocity-rescaling branch is gated "
         "`.and. .not. tnosep`, so Nose is not gated behind `ntcheck`.")
L.append("14. **Does Nose exchange measurable energy with the ions?** Indirectly yes (the trajectory "
         "is driven away from a pure microcanonical path), but the Nose chain energy is not written by "
         "QE, so a direct energy-exchange measurement is unavailable (evidence limit stated).")
L.append(f"15. **Active but ineffective, or inactive?** `{nose['classification']}` - active with the "
         "correct target but ineffective: it neither warms the cold early ions to 300 K nor prevents "
         "the later PR-cell-driven runaway.")
L.append(f"16. **Is CASE B a valid long-time no-op?** Yes - `{ab['classification']}`; bit-identical to "
         f"CASE A until step {ab.get('eka_first_divergence_step')}, then statistically equivalent "
         f"(both heat to ~5000 K); no CAIMD event fired (`any_event={ab.get('caimd_B_any_event')}`).")
L.append(f"17. **Is CASE C stable under active zero-strain projection?** Yes - "
         f"`{cC['classification']}`: cell held perfectly rigid (V=306.27 bohr^3, cond=1.0000 constant), "
         f"determinant always positive, max skew norm {cC.get('max_skew_norm'):.1e}.")
L.append(f"18. **CASE C maximum strain error:** {cC.get('max_abs_nominal_minus_actual'):.2e} (machine "
         f"precision); max velocity-transform residual {cC.get('max_velocity_jump_max'):.2e}.")
L.append(f"19. **Does active CAIMD change temperature behavior?** **Yes, materially.** A and B heat to "
         f"~{f(R['A']['final_5ps_mean_K'],0)}/{f(R['B']['final_5ps_mean_K'],0)} K; C only to "
         f"~{f(R['C']['final_5ps_mean_K'],0)} K. Freezing the cell (all six strains held at zero) "
         "removes the dominant PR-cell energy-injection channel and suppresses the heating ~8x.")
L.append("20. **Suitable for a nonzero-strain test?** Not yet. The zero-strain projection is "
         "numerically flawless, but the underlying ordinary `vc-md + Nose` baseline does not hold "
         "300 K, so a loading test would run on an uncontrolled thermal background.")
L.append("21. **Narrowest justified next step:** see recommendation below.")
L.append("22. **Sources or raw outputs modified?** No. Remote workflow read once (rsync pull); all new "
         "products under `analysis_audited/`, `evidence_audited/` and this report. No QE source touched, "
         "no job submitted/cancelled/polled, no commit/push.")
L.append("")
L.append("## Why A/B heat and C does not (mechanism)")
L.append("")
L.append("The heating is Parrinello-Rahman cell-driven. In A/B the free cell contracts ~9% "
         f"(V 306->{f(R['A']['cell']['volume_final_bohr3'],0)} bohr^3, min volume "
         f"~{f(R['A']['cell']['volume_min'],0)}) and pumps energy into the ions faster than the "
         "ionic Nose (which controls ions only; `cell_temperature='not_controlled'`) can remove it. "
         "CASE C holds all six symmetric strain components at zero, so the cell is rigid and that "
         "energy channel is closed; the residual mild heating (~520 K) is what the ionic Nose alone "
         "leaves uncontrolled.")
L.append("")
L.append("## Minimum next step (recommendation)")
L.append("")
L.append("Do **not** proceed to nonzero-strain scientific loading. The immediate diagnostic is to "
         "isolate the thermostat from the barostat: run a **fixed-cell `calculation='md'` + Nose** "
         "test (same native init, same `fnosep=50`, `nhpcl=4`, `dt=2 fs`, no PR cell, no CAIMD) to "
         "check whether the ionic Nose alone holds 300 K when there is no PR-cell energy source. "
         "In parallel, a conservative ordinary `vc-md + Nose` parameter check is warranted (the "
         "`wmass=2e6` PR cell mass and 2 fs step may make the PR cell dynamics the instability): try "
         "lower `fnosep`, shorter `dt`, and a cell thermostat / different `wmass`. Only after a "
         "stable 300 K ordinary baseline is established should the zero-strain and then small "
         "nonzero-strain CAIMD tests be repeated.")
L.append("")
L.append("## Provenance & caveats")
L.append("")
L.append(f"- Job 50434898, node lrdn1374; steps .0/.1/.2 all `COMPLETED 0:0`; A/B/C RC=0, `JOB DONE`.")
L.append(f"- Executable SHA256 {S['executable']['sha256_runtime']} (matches expected); banner v.7.5.")
L.append("- `ecutwfc=25 Ry`: absolute pressures are Pulay-dominated, not converged; used only for "
         "comparison and stability.")
L.append("- `.tv`/stdout temperature is a running average; instantaneous stats use `2*eka/((3N-3)k_B)`.")
L.append("- Nose chain conserved-energy terms are not written by QE; the extended-energy check is "
         "partial.")
L.append("")
with open(os.path.join(MIRROR, "FINAL_LONG_ZERO_STRAIN_REPORT_AUDITED.md"), "w") as fh:
    fh.write("\n".join(L) + "\n")
print("wrote FINAL_LONG_ZERO_STRAIN_REPORT_AUDITED.md")
