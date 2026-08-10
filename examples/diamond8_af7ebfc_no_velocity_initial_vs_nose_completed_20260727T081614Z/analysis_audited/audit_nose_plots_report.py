#!/usr/bin/env python3
"""Plots + audited report for job 50167241 (CONTROL initial vs TEST nose)."""
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

C = {"control": "#2a78d6", "nose": "#eb6834"}
LB = {"control": "CONTROL initial", "nose": "TEST nose (fnosep=50, nhpcl=4)"}
SRC = ("Source: job 50167241 mirror; T = 2*eka/((3N-3)k_B) from .eal; P from .p col2 (kbar); "
       "cell from .cel. ecutwfc=25 Ry -> pressures not converged.")
plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.color": "#e1e0d9",
                     "grid.linewidth": .8, "axes.axisbelow": True, "axes.unicode_minus": False,
                     "savefig.dpi": 150, "figure.facecolor": "white", "axes.facecolor": "white"})


def roll(a, w):
    if len(a) < w:
        return np.full(len(a), np.nan)
    return np.concatenate([np.full(w - 1, np.nan), np.convolve(a, np.ones(w) / w, "valid")])


def blk(a, w=100):
    return np.array([a[i:i + w].mean() for i in range(0, len(a), w) if len(a[i:i + w])])


def base(ax, title, xl, yl, note=SRC, legend=True, loc="best"):
    ax.set_title(title, fontsize=11.5); ax.set_xlabel(xl); ax.set_ylabel(yl)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if legend:
        ax.legend(frameon=False, fontsize=9, loc=loc)
    ax.figure.text(0.01, -0.02, note, fontsize=7.3, color="#52514e")


def save(fig, name):
    fig.savefig(os.path.join(ANA, name), bbox_inches="tight"); plt.close(fig)
    print("wrote", name)


t = {k: Z[f"{k}_t_ps"] for k in C}
T = {k: Z[f"{k}_T"] for k in C}
P = {k: Z[f"{k}_p_kbar"] for k in C}
E = {k: Z[f"{k}_edyn"] for k in C}
V = {k: Z[f"{k}_vol"] for k in C}
CN = {k: Z[f"{k}_cond"] for k in C}
CT = {k: Z[f"{k}_celt"] for k in C}

# 1 temperature
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], T[k], lw=.85, color=C[k], alpha=.9, label=LB[k])
ax.axhline(300, color="#898781", ls="--", lw=1.2)
ax.text(0.02, 305, "300 K target", fontsize=8.6, color="#52514e")
ax.set_ylim(40, 330)
base(ax, "Instantaneous ionic temperature: initial vs Nose (native init, no velocities)",
     "time (ps)", "T (K, 3N-3)", loc="lower right")
save(fig, "temperature_control_vs_nose.png")

# 2 rolling 25
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], roll(T[k], 25), lw=1.4, color=C[k], label=LB[k])
ax.axhline(300, color="#898781", ls="--", lw=1.2)
base(ax, "Temperature, 25-step rolling mean", "time (ps)", "rolling mean T (K)")
save(fig, "temperature_rolling25_control_vs_nose.png")

# 3 rolling 100
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], roll(T[k], 100), lw=1.7, color=C[k], label=LB[k])
ax.axhline(300, color="#898781", ls="--", lw=1.2)
base(ax, "Temperature, 100-step rolling mean", "time (ps)", "rolling mean T (K)")
save(fig, "temperature_rolling100_control_vs_nose.png")

# 4 block 100
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    b = blk(T[k]); ax.plot(np.arange(1, len(b) + 1), b, "o-", ms=5, lw=1.5, color=C[k], label=LB[k])
ax.axhline(300, color="#898781", ls="--", lw=1.2)
base(ax, "Temperature, 100-step block means", "100-step block", "block-mean T (K)")
save(fig, "temperature_block100_control_vs_nose.png")

# 5 nose temperature after peak
fig, ax = plt.subplots(figsize=(10, 5))
pk = int(np.argmax(T["nose"]))
ax.plot(t["nose"], T["nose"], lw=.7, color="#eb6834", alpha=.5, label="instantaneous")
ax.plot(t["nose"], roll(T["nose"], 100), lw=1.8, color="#b0431a", label="100-step rolling mean")
ax.axvline(t["nose"][pk], color="#898781", ls=":", lw=1.2)
ax.axhline(300, color="#898781", ls="--", lw=1.0)
ax.text(t["nose"][pk] + 0.02, 285, f"peak {T['nose'][pk]:.0f} K @ {t['nose'][pk]:.2f} ps (= initial)",
        fontsize=8.4, color="#52514e")
base(ax, "Nose temperature after its peak (no heating excursion: peak is the initial value)",
     "time (ps)", "T (K, 3N-3)")
save(fig, "nose_temperature_after_peak.png")

# 6 ionic kinetic energy
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], Z[f"{k}_eka"], lw=.9, color=C[k], label=LB[k])
base(ax, "Ionic kinetic energy (eka)", "time (ps)", "eka (Ry)")
save(fig, "ionic_kinetic_energy_control_vs_nose.png")

# 7 pressure
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], P[k], lw=.85, color=C[k], alpha=.9, label=LB[k])
ax.axhline(0, color="#898781", ls="--", lw=1.2)
base(ax, "MD pressure (target press=0; NOT a converged pressure)", "time (ps)", "P (kbar)")
save(fig, "pressure_control_vs_nose.png")

# 8 volume
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(CT[k], V[k], lw=1.4, color=C[k], label=LB[k])
base(ax, "Parrinello-Rahman cell volume", "time (ps)", "volume (bohr^3)")
save(fig, "volume_control_vs_nose.png")

# 9 condition number
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(CT[k], CN[k], lw=1.4, color=C[k], label=LB[k])
base(ax, "Cell condition number (1.0 = perfectly conditioned)", "time (ps)", "cond(H)")
save(fig, "cell_condition_control_vs_nose.png")

# 10 energy drift
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], (E[k] - E[k][0]) * 13.605693122994 / NAT, lw=1.3, color=C[k], label=LB[k])
ax.axhline(0, color="#898781", ls="--", lw=1.2)
base(ax, "Partial extended-energy drift (edyn = ekint + ut; excludes Nose chain)",
     "time (ps)", "drift per atom (eV/atom)",
     note=SRC + " Nose chain energy is not written by QE, so this is a PARTIAL diagnostic.")
save(fig, "energy_drift_control_vs_nose.png")

# 11 nose chain energy terms (not available -> honest placeholder)
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t["nose"], (E["nose"] - E["nose"][0]) * 13.605693122994 / NAT, lw=1.3, color="#eb6834",
        label="edyn drift (partial; excludes chain)")
ax.text(0.5, 0.5, "Nose chain energy terms (xnhp0, vnhp, ions_nose_energy)\n"
                  "are computed in PW/src/vcsmd.f90 but NOT written to any\n"
                  "output file, so the chain kinetic/potential terms and the\n"
                  "full Nose+PR conserved Hamiltonian cannot be reconstructed.",
        transform=ax.transAxes, ha="center", va="center", fontsize=9.5, color="#52514e",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f6f6f4", ec="#c3c2b7"))
base(ax, "Nose chain energy terms (unavailable in QE output)", "time (ps)",
     "partial edyn drift (eV/atom)")
save(fig, "nose_chain_energy_terms.png")

# 12 temperature vs volume (nose)
fig, ax = plt.subplots(figsize=(8.4, 6))
vt = np.interp(t["nose"], CT["nose"], V["nose"])
sc = ax.scatter(vt, T["nose"], c=t["nose"], cmap="viridis", s=8, alpha=.7)
ax.axhline(300, color="#898781", ls="--", lw=1.0)
cb = fig.colorbar(sc, ax=ax); cb.set_label("time (ps)")
base(ax, "Nose: temperature vs cell volume (colour = time)", "volume (bohr^3)", "T (K, 3N-3)",
     legend=False)
save(fig, "temperature_vs_volume_nose.png")

# ----------------------------------------------------------------- report ---
cc_, nn_ = S["runs"]["control"], S["runs"]["nose"]
ct, nt = cc_["temperature"], nn_["temperature"]
cp, npr = cc_["pressure_kbar"], nn_["pressure_kbar"]
REF = S["reference_values_verified"]
ii = S["initial_state"]


def f(x, n=3):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{n}f}"


L = []
L.append("# Audited Final Thermostat Isolation Report")
L.append("")
L.append("Diamond8, QE 7.5 `af7ebfc`, **no explicit velocities**, `caimd=.false.` in both runs. "
         "Single changing variable: `ion_temperature='initial'` (CONTROL) vs `'nose'` with "
         "`fnosep=50`, `nhpcl=4` (TEST). Job **50167241**.")
L.append("")
L.append(f"**Overall status: `{S['overall_status']}`**")
L.append("")
L.append(f"- CONTROL: `{cc_['classification']}` · NOSE: `{nn_['classification']}` · "
         f"workflow: `{S['workflow_class']}`")
L.append("")
L.append("## Headline result")
L.append("")
L.append("**The Nose run does not heat.** With native velocity initialization (no `ATOMIC_VELOCITIES`) "
         "the Nose-Hoover trajectory is bit-identical to the control for the first "
         f"{ii['first_eka_divergence_step']} ionic steps and settles onto the same ~170 K plateau; "
         f"its maximum instantaneous temperature is {f(nt['max'],3)} K — exactly the step-0 value — "
         "and it never crosses 600 K. This is the opposite of the earlier explicit-velocity Nose run "
         f"(job 50127909), which reached {REF['nose_previous_with_explicit_velocities_job50127909']['max_K']:.0f} K. "
         "The Nose thermostat was confirmed active in the TEST stdout "
         f"(`{S['nose_thermostat_active']}`). Therefore the earlier several-thousand-kelvin heating "
         "**requires the explicit-velocity input path**, not the Nose thermostat by itself, and not "
         "active CAIMD (disabled in every run).")
L.append("")
L.append("The plateau sits near ~170 K rather than 300 K because the run starts at the ideal lattice "
         "(a potential-energy minimum) with all energy kinetic, so equipartition halves the kinetic "
         "temperature; at these parameters the Nose chain does not pull it back to 300 K within 2 ps "
         "(it barely acts — first divergence from the control only at 1.058 ps). That is a weak-coupling "
         "/ short-time observation, **not** heating.")
L.append("")
L.append("## Required answers")
L.append("")
L.append("1. **Job ID and final Slurm states:** `50167241` — top-level, `batch`, `extern` all "
         "`COMPLETED 0:0`; step `.0` (control `pw.x`) and `.1` (nose `pw.x`) both `COMPLETED 0:0`; "
         "elapsed 00:03:44 on `lrdn0558`.")
L.append(f"2. **Did both QE calculations complete?** Yes. Both printed `JOB DONE`, returned 0, and "
         f"wrote 1000 `.eal`/`.p`/`.tv`/`.cel` rows. Classifications `{cc_['classification']}` / "
         f"`{nn_['classification']}`; SCF failures 0; runtime version `{cc_['version']}` (both).")
L.append(f"3. **Explicit velocities absent from both inputs?** Yes — "
         f"`{S['input_contract']['explicit_velocities_absent']}`: no `ATOMIC_VELOCITIES` header, "
         f"no velocity rows, no `ion_velocities`/`from_input` token in the complete text of either "
         f"input, both `restart_mode='from_scratch'`.")
L.append(f"4. **`caimd=.false.` in both?** Yes — `{S['input_contract']['both_caimd_false']}` "
         f"(exactly one `caimd = .false.` line each; no `caimd_stride`/`caimd_step_strain`/"
         f"`caimd_mask`; no `caimd=.true.`).")
L.append(f"5. **All non-thermostat settings identical?** Yes — `{S['input_contract']['valid_single_factor']}`. "
         f"Maximum difference 0 for cell, positions, ordering, masses, cutoffs, timestep, `nstep`, "
         f"`press`, `wmass`, `tempw`, `cell_dynamics`, `cell_dofree`, k-points and every SCF "
         f"parameter. Only `prefix`, `ion_temperature`, `fnosep(1)` and `nhpcl` differ.")
L.append(f"6. **Same native initial velocity realization?** Yes — initial-state class "
         f"`{ii['classification']}`. Step-0 and step-1 ionic kinetic energies are bit-identical and "
         f"the `.eal` first diverges only at step {ii['first_eka_divergence_step']} (1.058 ps). "
         f"`ranv` zeroes total linear momentum, so both start with zero COM velocity. This is a "
         f"strict single-factor trajectory comparison.")
L.append(f"7. **Did the control reproduce the previous bounded result?** Yes. New control mean "
         f"{f(ct['mean'])} K, final-200 {f(ct['final200_mean'])} K, max {f(ct['max'],3)} K, 0 records "
         f">600 K, volume change {f(cc_['cell']['volume_change_percent'])}%. Previous audited "
         f"reference (job 50160152): mean {REF['control_previous_no_velocity_initial']['mean_K']:.3f} K, "
         f"final-200 {REF['control_previous_no_velocity_initial']['final200_K']:.3f} K, volume "
         f"{REF['control_previous_no_velocity_initial']['volume_change_percent']:.3f}%. "
         f"Reproduced within small-system variation.")
L.append(f"8. **Nose peak temperature and when:** {f(nn_['peak']['value_K'],3)} K at step "
         f"{nn_['peak']['nst']} (t = {f(nn_['peak']['time_ps'])} ps) — i.e. the peak IS the initial "
         f"value; there is no heating excursion above it.")
L.append(f"9. **Nose 100-step block means (K):** "
         f"{[round(x,1) for x in nn_['block100_means_K']]}.")
L.append(f"10. **Did the temperature decline after the peak?** There is no high-temperature peak to "
         f"decline from; the rolling mean stays in a ~160-177 K band the whole run (lowest 100-step "
         f"rolling mean after step 0 = {f(nn_['peak']['lowest_roll100_after_peak_K'])} K).")
L.append(f"11. **Did it return below 2000/1000/600 K?** It never rose above any of them: "
         f"fractions above 600/1000/2000 K = "
         f"{nn_['fraction_above']['600K']}/{nn_['fraction_above']['1000K']}/"
         f"{nn_['fraction_above']['2000K']}.")
L.append(f"12. **Final-100 / final-200 / final-250 means:** {f(nt['final100_mean'])} / "
         f"{f(nt['final200_mean'])} / {f(nt['final250_mean'])} K (Nose); "
         f"{f(ct['final100_mean'])} / {f(ct['final200_mean'])} / {f(ct['final250_mean'])} K (control).")
L.append(f"13. **Final-half and final-quarter temperature slopes (Nose):** "
         f"{f(nn_['drift_final_half_K_per_ps'])} K/ps and {f(nn_['drift_final_quarter_K_per_ps'])} K/ps "
         f"(small, sign-alternating — consistent with a bounded plateau, not a trend).")
L.append("14. **Evidence of convergence toward 300 K?** No secular approach to 300 K within 2 ps; "
         "the trajectory is bounded near ~170 K. The thermostat is simply too weakly/slowly coupled "
         "at these parameters over this window to reach the target — but there is no divergence.")
L.append(f"15. **Did direct `.md` velocities confirm the reported temperature?** Yes. Terminal direct "
         f"T (from `.md` `ratd` x terminal `.cel`, COM removed): control {f(cc_['direct_terminal_T_K'])} K "
         f"vs `.eal` {f(ct['final'])} K; nose {f(nn_['direct_terminal_T_K'])} K vs {f(nt['final'])} K. "
         f"Confirmed (control={cc_['direct_confirms_reported']}, nose={nn_['direct_confirms_reported']}); "
         f"no diagnostic mismatch.")
L.append(f"16. **Did temperature growth precede cell deformation?** Moot — there was no temperature "
         f"growth. Nose never exceeded 600 K.")
L.append(f"17. **Was the cell stable?** Yes. Both `{cc_['cell']['classification']}`: volume change "
         f"control {f(cc_['cell']['volume_change_percent'])}% / nose "
         f"{f(nn_['cell']['volume_change_percent'])}%, max condition number "
         f"{f(cc_['cell']['condition_number_max'],5)} / {f(nn_['cell']['condition_number_max'],5)}, "
         f"min C-C {f(cc_['cell']['min_CC_bohr'])} / {f(nn_['cell']['min_CC_bohr'])} bohr.")
L.append(f"18. **Was the Nose/PR extended energy bounded?** The available quantity `edyn` "
         f"(`.e` col3 = ekint+ut) is bounded (control drift "
         f"{cc_['energy_partial']['drift_eV_per_atom_per_ps']:.2e} eV/atom/ps, nose "
         f"{nn_['energy_partial']['drift_eV_per_atom_per_ps']:.2e} eV/atom/ps). "
         f"**This is a PARTIAL diagnostic**: `ions_nose_energy` and the chain variables are computed "
         f"in `vcsmd` but never written to any output file, so the full Nose+PR conserved Hamiltonian "
         f"cannot be formed. No secular energy injection is visible.")
L.append("19. **Are explicit velocities required to reproduce the heating?** On this evidence, "
         "**yes** — Nose with native initialization is bounded, while the earlier run that heated "
         "differed only by supplying explicit `ATOMIC_VELOCITIES`. This is the direct isolation the "
         "task set out to make.")
L.append("20. **Is active CAIMD required?** No. CAIMD was disabled (`caimd=.false.`) in both runs "
         "here and in the earlier heating run.")
L.append("21. **Correct conclusion — transient convergence, late unresolved recovery, or sustained "
         "heating?** **None of those** — there is no heating at all. The native-initialization Nose "
         "trajectory is bounded from the start; the earlier heating is attributable to the "
         "explicit-velocity path.")
L.append("22. **Minimum next test:** see below.")
L.append("23. **Were any source files or raw remote outputs modified?** No. The remote workflow was "
         "read once (rsync pull); all new products are under `analysis_audited/`, `evidence_audited/` "
         "and this report. No QE source touched, no job submitted/cancelled/polled, no commit/push.")
L.append("")
L.append("## Comparison table")
L.append("")
L.append("| quantity | CONTROL (`initial`) | TEST (`nose`) |")
L.append("|---|---:|---:|")
for lab, a, b in (
        ("classification", cc_["classification"], nn_["classification"]),
        ("Nose active in stdout", "n/a", S["nose_thermostat_active"]),
        ("steps / JOB DONE", f"{cc_['completed_steps']} / {cc_['job_done']}",
         f"{nn_['completed_steps']} / {nn_['job_done']}"),
        ("initialized T (K)", f(cc_["temperature_initialized_K"], 4), f(nn_["temperature_initialized_K"], 4)),
        ("mean T (K)", f(ct["mean"]), f(nt["mean"])),
        ("final-200 mean T (K)", f(ct["final200_mean"]), f(nt["final200_mean"])),
        ("max T (K)", f(ct["max"], 3), f(nt["max"], 3)),
        ("min T (K)", f(ct["min"]), f(nt["min"])),
        ("fraction >600 K", cc_["fraction_above"]["600K"], nn_["fraction_above"]["600K"]),
        ("drift full (K/ps)", f(cc_["drift_full_K_per_ps"]), f(nn_["drift_full_K_per_ps"])),
        ("drift final-quarter (K/ps)", f(cc_["drift_final_quarter_K_per_ps"]),
         f(nn_["drift_final_quarter_K_per_ps"])),
        ("mean P (kbar)", f(cp["mean"], 2), f(npr["mean"], 2)),
        ("volume change (%)", f(cc_["cell"]["volume_change_percent"]), f(nn_["cell"]["volume_change_percent"])),
        ("max cond(H)", f(cc_["cell"]["condition_number_max"], 5), f(nn_["cell"]["condition_number_max"], 5)),
        ("terminal direct-vel T (K)", f(cc_["direct_terminal_T_K"]), f(nn_["direct_terminal_T_K"])),
):
    L.append(f"| {lab} | {a} | {b} |")
L.append("")
L.append("## Three-way context")
L.append("")
L.append("| run | thermostat | velocities | outcome |")
L.append("|---|---|---|---|")
L.append("| job 50160152 (control ref) | initial | native | bounded, ~170 K |")
L.append(f"| **this, control** | initial | native | bounded, {f(ct['mean'],0)} K mean |")
L.append(f"| **this, TEST** | **nose** | native | **bounded, {f(nt['mean'],0)} K mean, max {f(nt['max'],0)} K** |")
L.append(f"| job 50127909 | nose | **explicit** | **heated to "
         f"{REF['nose_previous_with_explicit_velocities_job50127909']['max_K']:.0f} K** |")
L.append("")
L.append("The only factor that flips the outcome from bounded to heated is the explicit-velocity "
         "input path.")
L.append("")
L.append("## Minimum next test (not submitted)")
L.append("")
L.append("Because native-initialization Nose is bounded while the previous explicit-velocity Nose "
         "heated, move the variable to the velocity path (single factor at a time), keeping "
         "`caimd=.false.`: (a) `ion_temperature='initial'` **with** an explicit `ATOMIC_VELOCITIES` "
         "card and **no** Nose; then (b) Nose **plus** the same explicit velocities. Audit the "
         "first-step ionic kinetic energy, the Nose chain kinetic input `ekin2nhp`, and the Nose "
         "chain variables (`xnhp0`, `vnhp`) at initialization when velocities are supplied — that "
         "initialization is the remaining structural difference against job 50127909.")
L.append("")
L.append("## Caveats")
L.append("")
L.append("- `ecutwfc=25 Ry` retained from the diagnostic deck; absolute pressures are "
         "Pulay-dominated and not converged. Pressure/volume are used only for comparison and cell "
         "stability.")
L.append("- The `.tv`/stdout temperature is a running average (`2*avk/(3(N+1)k_B)`); all "
         "instantaneous statistics use `2*eka/((3N-3)k_B)` from `.eal`.")
L.append("- `tmp/*.vel` is identically zero in the PR path and is not used; `tmp/*.md` holds only "
         "the terminal state, so direct velocity reconstruction is available at the final frame only.")
L.append("- The Nose chain conserved-energy terms are not written by QE, so the extended-energy "
         "check is a partial diagnostic (documented above).")
L.append("- 2 ps is short; conclusions cover this window. The result is a clear negative for heating, "
         "not a claim that Nose reaches 300 K given more time.")
L.append("")
with open(os.path.join(MIRROR, "FINAL_THERMOSTAT_ISOLATION_REPORT_AUDITED.md"), "w") as fh:
    fh.write("\n".join(L) + "\n")
print("wrote FINAL_THERMOSTAT_ISOLATION_REPORT_AUDITED.md")
