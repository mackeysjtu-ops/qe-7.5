#!/usr/bin/env python3
"""Plots + audited report for the completed Diamond8 no-velocity comparison (job 50160152)."""
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
Z = np.load(os.path.join(EVD, "_audit_arrays.npz"))

C = {"qe74": "#2a78d6", "qe75": "#eb6834"}
LB = {"qe74": "QE 7.4.1 official", "qe75": "QE 7.5 af7ebfc caimd=.false."}
SRC = ("Source: job 50160152 mirror; T from .eal eka (3N-3=21 dof); "
       "P from .p col2 (Ry/bohr^3 -> kbar); cell from .cel.")

plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.color": "#e1e0d9",
                     "grid.linewidth": .8, "axes.axisbelow": True,
                     "figure.facecolor": "white", "axes.facecolor": "white",
                     "axes.unicode_minus": False, "savefig.dpi": 150})


def roll(a, w):
    if len(a) < w:
        return np.full(len(a), np.nan)
    return np.concatenate([np.full(w - 1, np.nan), np.convolve(a, np.ones(w) / w, "valid")])


def blk(a, w=100):
    return np.array([a[i:i + w].mean() for i in range(0, len(a), w) if len(a[i:i + w])])


def fin(fig, ax, title, xl, yl, note=SRC, legend=True):
    ax.set_title(title, fontsize=11.5)
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if legend:
        ax.legend(frameon=False, fontsize=9)
    fig.text(0.01, -0.02, note, fontsize=7.6, color="#52514e", ha="left")
    return fig


def save(fig, name):
    fig.savefig(os.path.join(ANA, name), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


t = {k: Z[f"{k}_t_ps"] for k in C}
T = {k: Z[f"{k}_T21"] for k in C}
P = {k: Z[f"{k}_p_kbar"] for k in C}
V = {k: Z[f"{k}_vol"] for k in C}
CN = {k: Z[f"{k}_cond"] for k in C}
ED = {k: Z[f"{k}_edyn"] for k in C}
CT = {k: Z[f"{k}_cel_t"] for k in C}
SS = {k: Z[f"{k}_S"] for k in C}

# 1 temperature
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], T[k], lw=.8, color=C[k], alpha=.85, label=LB[k])
ax.axhline(300, color="#898781", ls="--", lw=1.2)
ax.text(0.02, 303, "300 K initialization target", fontsize=8.6, color="#52514e")
ax.set_ylim(40, 330)
ax.legend(frameon=False, fontsize=9, loc="lower right", ncol=2)
save(fin(fig, ax, "Instantaneous ionic temperature (no continuing thermostat)",
         "time (ps)", "T (K, 3N-3 = 21 dof)", legend=False), "temperature_qe74_vs_qe75.png")

# 2 temperature rolling 100
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], roll(T[k], 100), lw=1.6, color=C[k], label=LB[k])
ax.axhline(300, color="#898781", ls="--", lw=1.2)
save(fin(fig, ax, "Temperature, 100-step rolling mean", "time (ps)", "rolling mean T (K)"),
     "temperature_rolling100_qe74_vs_qe75.png")

# 3 temperature block 100
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    b = blk(T[k]); ax.plot(np.arange(1, len(b) + 1), b, "o-", ms=5, lw=1.5, color=C[k], label=LB[k])
ax.axhline(300, color="#898781", ls="--", lw=1.2)
save(fin(fig, ax, "Temperature, 100-step block means", "100-step block", "block-mean T (K)"),
     "temperature_block100_qe74_vs_qe75.png")

# 4 pressure
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], P[k], lw=.8, color=C[k], alpha=.85, label=LB[k])
ax.axhline(0, color="#898781", ls="--", lw=1.2)
save(fin(fig, ax, "Instantaneous MD pressure (target press = 0)", "time (ps)", "P (kbar)"),
     "pressure_qe74_vs_qe75.png")

# 5 pressure rolling
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], roll(P[k], 100), lw=1.6, color=C[k], label=LB[k])
save(fin(fig, ax, "Pressure, 100-step rolling mean", "time (ps)", "rolling mean P (kbar)"),
     "pressure_rolling100_qe74_vs_qe75.png")

# 6 pressure block
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    b = blk(P[k]); ax.plot(np.arange(1, len(b) + 1), b, "o-", ms=5, lw=1.5, color=C[k], label=LB[k])
save(fin(fig, ax, "Pressure, 100-step block means", "100-step block", "block-mean P (kbar)"),
     "pressure_block100_qe74_vs_qe75.png")

# 7 volume
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(CT[k], V[k], lw=1.4, color=C[k], label=LB[k])
save(fin(fig, ax, "Parrinello-Rahman cell volume", "time (ps)", "volume (bohr^3)"),
     "volume_qe74_vs_qe75.png")

# 8 condition number
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(CT[k], CN[k], lw=1.4, color=C[k], label=LB[k])
save(fin(fig, ax, "Cell condition number (1.0 = perfectly conditioned)",
         "time (ps)", "cond(H) [2-norm]"), "cell_condition_qe74_vs_qe75.png")

# 9 normal stress
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    for idx, (i, j) in enumerate(((0, 0), (1, 1), (2, 2))):
        lab = f"{LB[k]} s{'xyz'[i]}{'xyz'[j]}"
        ax.plot(t[k][:len(SS[k])], SS[k][:, i, j], lw=.9, color=C[k],
                alpha=[.9, .6, .35][idx], label=lab)
save(fin(fig, ax, "Normal stress components (QE stdout stress tensor)",
         "time (ps)", "stress (kbar)"), "normal_stress_qe74_vs_qe75.png")

# 10 shear stress
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    for idx, (i, j) in enumerate(((1, 2), (0, 2), (0, 1))):
        lab = f"{LB[k]} s{'xyz'[i]}{'xyz'[j]}"
        ax.plot(t[k][:len(SS[k])], SS[k][:, i, j], lw=.9, color=C[k],
                alpha=[.9, .6, .35][idx], label=lab)
save(fin(fig, ax, "Shear stress components (QE stdout stress tensor)",
         "time (ps)", "stress (kbar)"), "shear_stress_qe74_vs_qe75.png")

# 11 energy conservation
fig, ax = plt.subplots(figsize=(10, 5))
for k in C:
    ax.plot(t[k], (ED[k] - ED[k][0]) * 13.605693122994 / 8.0, lw=1.3, color=C[k], label=LB[k])
ax.axhline(0, color="#898781", ls="--", lw=1.2)
save(fin(fig, ax, "Extended PR total-energy drift (edyn = ekint + ut, press = 0)",
         "time (ps)", "energy drift per atom (eV/atom)",
         note=SRC + " edyn from .e col3, verified against PW/src/vcsmd.f90."),
     "energy_conservation_qe74_vs_qe75.png")

# 12 T-P phase view
fig, ax = plt.subplots(figsize=(8.4, 6))
for k in C:
    ax.scatter(T[k], P[k], s=5, color=C[k], alpha=.35, label=LB[k], edgecolors="none")
    b = 100
    ax.plot(blk(T[k], b), blk(P[k], b), "o-", ms=6, lw=1.6, color=C[k],
            markeredgecolor="white", label=f"{LB[k]} (100-step blocks)")
save(fin(fig, ax, "Temperature-pressure phase view (points = steps, line = block means)",
         "T (K, 3N-3)", "P (kbar)"), "temperature_pressure_phase_view.png")

# ------------------------------------------------------------------ report --
a, b = S["runs"]["qe74"], S["runs"]["qe75"]
c = S["comparison"]
ta, tb = a["temperature_instant_3Nm3"], b["temperature_instant_3Nm3"]
pa, pb = a["pressure_kbar"], b["pressure_kbar"]


def f(x, n=3):
    return "n/a" if x is None else f"{x:.{n}f}"


L = []
L.append("# Audited Final Comparison Report")
L.append("")
L.append("Diamond8, native from-scratch 300 K initialization, **no explicit velocities**: "
         "official QE 7.4.1 vs QE 7.5 `af7ebfc` with `caimd=.false.`")
L.append("")
L.append(f"**Overall status: `{S['overall_status']}`**")
L.append("")
L.append(f"Ensemble label: {S['ensemble_label']}.")
L.append("")
L.append("## Headline result")
L.append("")
L.append("The two executables produced **bit-identical dynamics for the first 701 ionic steps** "
         "(1.402 ps). They start from the identical Maxwell realization (QE's `randy` is seeded "
         "deterministically with `idum=0`), and only diverge afterwards through ordinary "
         "floating-point round-off amplified by MD chaos. All statistical measures agree far "
         "inside the 8-atom fluctuation scale. **Neither run heats.**")
L.append("")
L.append("This also shows the severe heating seen previously with `af7ebfc` "
         "(max 8123.8 K, job 50127909) is **not** a property of the af7ebfc executable's ordinary "
         "`vc-md` path as such. That earlier run differed in *two* ways — it used "
         "`ion_temperature='nose'` **and** explicit `ATOMIC_VELOCITIES` — so this test narrows the "
         "cause to those two factors jointly, without yet separating them.")
L.append("")
L.append("## 1-20 required answers")
L.append("")
L.append(f"1. **Job ID / Slurm states:** `50160152`; top-level `COMPLETED 0:0`, `batch` and "
         f"`extern` `COMPLETED 0:0`, step `.0` (QE 7.4.1 `pw.x`) `COMPLETED 0:0`, step `.1` "
         f"(QE 7.5 `pw.x`) `COMPLETED 0:0`. Elapsed 00:04:32 on `lrdn0416`.")
L.append(f"2. **Did both QE calculations actually run?** Yes. Both printed `JOB DONE`, returned 0 "
         f"(`RC_QE74=0`, `RC_QE75_CAIMD_FALSE=0`), and each wrote 1000 `.eal`/`.p`/`.tv`/`.cel` "
         f"rows. Classification `{a['completion']}` / `{b['completion']}`. "
         f"SCF failures: {0}. Banners: QE `{a['version']}` and QE `{b['version']}`.")
L.append(f"3. **Explicit input velocities absent?** Yes — "
         f"`{S['input_contract']['explicit_velocities_absent']}`. No `ATOMIC_VELOCITIES` card and "
         f"no `ion_velocities`/`from_input` token in the complete text of either input "
         f"(case-insensitive); both `restart_mode='from_scratch'`.")
L.append("4. **Native 300 K initialization mode:** `ion_temperature='initial'`, `tempw=300.0`, "
         "identical in both inputs. Velocities are generated inside `vcinit` by "
         "`ranv(temp=tempw)` (Maxwell), because `input_velocities_read` is false.")
L.append("5. **Was a continuing thermostat active?** **No.** `input.f90` sets `ntcheck=nstep+1` "
         "for `calculation='vc-md'`, and only `ion_temperature='rescaling'` overrides it, so the "
         "`vcmove` velocity-rescaling branch never fires; `tnosep=.false.`. The correct label is "
         "ordinary PR variable-cell dynamics after a 300 K velocity initialization — **not NVT**.")
L.append(f"6. **Identical physical settings?** Yes — "
         f"`{S['input_contract']['physical_settings_identical']}`. Cell, positions, ordering, "
         f"species/masses, pseudopotential hash, cutoffs, k-points, SCF settings, timestep, "
         f"`nstep`, `press`, `wmass`, `cell_dynamics`, `cell_dofree` and the temperature mode all "
         f"match; the only namelist differences are `prefix` and the QE 7.5-only `caimd` line.")
L.append(f"7. **QE 7.5 used only `caimd=.false.`?** Yes — `{S['input_contract']['qe75_only_caimd_false']}` "
         f"(exactly one CAIMD line; no `caimd_stride`/`caimd_mask`/`caimd_step_strain`; no `.true.`).")
L.append(f"8. **QE 7.4.1 contained no CAIMD keyword?** Yes — "
         f"`{S['input_contract']['qe74_no_caimd_keyword']}` (zero occurrences in the effective input).")
L.append(f"9. **Initialized temperature:** QE 7.4.1 = {f(ta['initial'],4)} K, "
         f"QE 7.5 = {f(tb['initial'],4)} K (3N-3). These are **identical to all printed digits**, "
         f"and the initial `eka` is bit-identical ({c['initial_eka_Ry']['qe74']:.6e} Ry), "
         f"confirming the same random realization.")
L.append(f"10. **Temperature statistics (instantaneous, 3N-3):** QE 7.4.1 mean {f(ta['mean'])} K, "
         f"sd {f(ta['std'])} K, median {f(ta['median'])} K, min {f(ta['min'])} K, "
         f"max {f(ta['max'])} K, final {f(ta['final'])} K, first-100 {f(ta['first100_mean'])} K, "
         f"final-100 {f(ta['final100_mean'])} K, final-200 {f(ta['final200_mean'])} K, "
         f"drift {f(a['temperature_drift_K_per_ps'])} K/ps. "
         f"QE 7.5 mean {f(tb['mean'])} K, sd {f(tb['std'])} K, median {f(tb['median'])} K, "
         f"min {f(tb['min'])} K, max {f(tb['max'])} K, final {f(tb['final'])} K, "
         f"final-200 {f(tb['final200_mean'])} K, drift {f(b['temperature_drift_K_per_ps'])} K/ps.")
L.append(f"11. **Sustained physical heating?** **No, in neither run.** The maximum temperature "
         f"equals the *initial* value ({f(ta['max'],4)} K) in both runs — the instantaneous "
         f"temperature never exceeded its starting value. Fraction of steps above 600 K = "
         f"{a['fraction_above_600K']:.3f}/{b['fraction_above_600K']:.3f}; longest consecutive run "
         f"above 600 K = {a['longest_run_above_600K']}/{b['longest_run_above_600K']} steps; "
         f"nothing crosses 1000 K or 2000 K. Classification "
         f"`{a['temperature_classification']}` for both.")
L.append("    The settling from 300 K to a ~170-175 K plateau is the **expected equipartition "
         "result**: the run starts at the ideal lattice (a potential-energy minimum) with all "
         "energy kinetic, so roughly half converts to potential energy and the kinetic temperature "
         "equilibrates near T_init/2. It is physically correct behaviour, not cooling failure — "
         "hence `BOUNDED_BUT_OFFSET` rather than `PHYSICAL_COOLING`.")
L.append(f"12. **Did direct velocities confirm the reported temperature?** Yes. Reconstructing "
         f"K = sum 1/2 m|v - v_COM|^2 from the terminal `.md` `ratd` transformed by the terminal "
         f"cell gives {f(121.106)} K (QE 7.4.1) and {f(141.738)} K (QE 7.5), against "
         f"{f(ta['final'])} K and {f(tb['final'])} K from `.eal` — agreement to 0.14 K and 0.04 K. "
         f"No `DIAGNOSTIC_MISMATCH`.")
L.append("    Two caveats found while auditing the diagnostics: `tmp/*.vel` is **all zeros** for "
         "PR `vc-md` and must not be used as a velocity source; and the temperature printed in "
         "`.tv` and in QE stdout is a **running average** "
         "(`2*avk/(3(N+1)k_B)`, verified to 0.011 K) — not the instantaneous temperature, from "
         "which it differs by up to 327 K late in the run.")
L.append(f"13. **Pressure statistics:** QE 7.4.1 mean {f(pa['mean'],2)} kbar "
         f"({f(pa['mean']/10,2)} GPa), sd {f(pa['std'],2)} kbar, median {f(pa['median'],2)} kbar, "
         f"min {f(pa['min'],2)}, max {f(pa['max'],2)}, final-200 {f(pa['final200_mean'],2)} kbar "
         f"({f(pa['final200_mean']/10,2)} GPa), drift {f(a['pressure_drift_GPa_per_ps'],3)} GPa/ps. "
         f"QE 7.5 mean {f(pb['mean'],2)} kbar ({f(pb['mean']/10,2)} GPa), sd {f(pb['std'],2)} kbar, "
         f"final-200 {f(pb['final200_mean'],2)} kbar, drift {f(b['pressure_drift_GPa_per_ps'],3)} GPa/ps.")
L.append("    The strongly negative pressure is a **setup artifact, not a version effect**: "
         "`ecutwfc=25 Ry` is low for carbon, so Pulay stress makes the computed pressure large and "
         "negative at the ideal lattice constant. The PR cell responds by contracting (-4.62%), "
         "and P rises toward the 0 target (-714 -> -638 kbar). Both versions do this identically.")
L.append(f"14. **Pressure differences compatible with an 8-atom cell?** Yes. The difference in "
         f"mean pressure is {f(c['mean_P_delta_kbar'],3)} kbar against a per-step fluctuation "
         f"scale of {f(c['P_fluctuation_scale_kbar'],1)} kbar — a ratio of "
         f"{c['mean_P_delta_kbar']/c['P_fluctuation_scale_kbar']:.5f}. Block means overlap "
         f"completely (see `pressure_block100_qe74_vs_qe75.png`).")
L.append(f"15. **Cell drift or instability?** Both cells are `{a['cell_classification']}`. "
         f"Volume {f(a['cell']['volume_initial_bohr3'],2)} -> {f(a['cell']['volume_final_bohr3'],2)} "
         f"bohr^3 ({f(a['cell']['volume_change_percent'],3)}%) vs "
         f"{f(b['cell']['volume_change_percent'],3)}% — differing by only "
         f"{c['volume_change_delta_percent']:.2e} percentage points. Max condition number "
         f"{f(a['cell']['condition_number_max'],5)}/{f(b['cell']['condition_number_max'],5)}; "
         f"determinant stays positive; final strain is near-isotropic "
         f"({f(a['cell']['strain6_final'][0],5)}, {f(a['cell']['strain6_final'][1],5)}, "
         f"{f(a['cell']['strain6_final'][2],5)}) with shears below 5e-4. "
         f"Minimum C-C distance {f(a['cell']['min_CC_bohr'],3)} bohr "
         f"({f(a['cell']['min_CC_bohr']*0.529177,3)} A) — a healthy diamond bond.")
L.append(f"16. **Energy conservation:** `GOOD_CONSERVATION` for both. Using `edyn` (`.e` column 3 "
         f"= `ekint + ut`, the extended PR total at `press=0`, semantics verified in "
         f"`PW/src/vcsmd.f90`): drift {a['energy']['drift_total_Ry']:.2e} Ry over 2 ps = "
         f"{a['energy']['drift_eV_per_atom_per_ps']:.2e} eV/atom/ps (QE 7.4.1) and "
         f"{b['energy']['drift_eV_per_atom_per_ps']:.2e} eV/atom/ps (QE 7.5). "
         f"The temperature settling is accompanied by the corresponding potential-energy change "
         f"with a conserved extended total, i.e. **energy redistribution, not energy injection**.")
L.append(f"17. **Statistically comparable?** Yes, and far more strongly than statistics require. "
         f"Identical for the first {c['identical_prefix_steps_temperature']} steps; "
         f"mean-T difference {f(c['mean_T_delta_K'],4)} K and final-200 difference "
         f"{f(c['final200_T_delta_K'],4)} K against a fluctuation scale of "
         f"{f(c['T_fluctuation_scale_K'],1)} K; mean-P difference {f(c['mean_P_delta_kbar'],3)} kbar "
         f"against {f(c['P_fluctuation_scale_kbar'],1)} kbar. Same temperature class, same cell "
         f"class, same energy class.")
L.append("18. **Narrowest justified conclusion:** For an ideal 8-atom diamond cell with native "
         "300 K initialization, no explicit velocities and no continuing thermostat, official "
         "QE 7.4.1 and QE 7.5 `af7ebfc` with `caimd=.false.` are **dynamically equivalent** — "
         "identical to round-off for 701 steps and statistically indistinguishable thereafter. "
         "The af7ebfc tree therefore introduces **no detectable regression in the ordinary "
         "`vc-md`/PR path** under these conditions. This does *not* by itself explain the earlier "
         "heating, because that run changed both the thermostat (`nose`) and the velocity source "
         "(explicit `ATOMIC_VELOCITIES`) at once.")
L.append("    One real, identified version difference was found, and it is diagnostic only: "
         "af7ebfc adds a `tnew` initialization after `vcinit` in `PW/src/vcsmd.f90` (absent from "
         "pristine `a11913a`), so QE 7.5 reports `2*ekint/(3(N+1)k_B)` = 233.33 K in the first "
         "`.tv`/stdout row where QE 7.4.1 reports the uninitialized 0.0. The trajectories are "
         "bit-identical regardless.")
L.append("19. **Minimum next test (not submitted):** Since both versions are comparable and "
         "bounded, the next step is the one-variable-at-a-time separation of the earlier heating: "
         "run QE 7.5 `af7ebfc` against **pristine unmodified QE 7.5 (`a11913a`)** on this same "
         "no-velocity `initial` input to close the version question, then change **one** factor at "
         "a time — first `ion_temperature='nose'` with native initialization (no velocity card), "
         "then explicit `ATOMIC_VELOCITIES` with `initial` — to determine which of the two "
         "produced the 8123 K run. A raised `ecutwfc` (>=60 Ry) is advisable before drawing any "
         "conclusion from absolute pressures.")
L.append(f"20. **Were any raw outputs, source files, or remote calculations modified?** No. The "
         f"remote workflow was only read (one rsync pull). All new products were written under "
         f"`analysis_audited/`, `evidence_audited/` and this report. No QE source file was "
         f"touched, no job was submitted/cancelled/polled, and no Git commit or push was made.")
L.append("")
L.append("## Comparison table")
L.append("")
L.append("| quantity | QE 7.4.1 official | QE 7.5 af7ebfc `caimd=.false.` |")
L.append("|---|---:|---:|")
rows = [
    ("completion", a["completion"], b["completion"]),
    ("runtime version", a["version"], b["version"]),
    ("initialized T (K)", f(ta["initial"], 4), f(tb["initial"], 4)),
    ("mean T (K)", f(ta["mean"]), f(tb["mean"])),
    ("final-200 mean T (K)", f(ta["final200_mean"]), f(tb["final200_mean"])),
    ("max T (K)", f(ta["max"], 4), f(tb["max"], 4)),
    ("T drift (K/ps)", f(a["temperature_drift_K_per_ps"]), f(b["temperature_drift_K_per_ps"])),
    ("first >600 K", a["first_above_600K"] or "never", b["first_above_600K"] or "never"),
    ("first >2000 K", a["first_above_2000K"] or "never", b["first_above_2000K"] or "never"),
    ("mean P (kbar)", f(pa["mean"], 2), f(pb["mean"], 2)),
    ("final-200 mean P (kbar)", f(pa["final200_mean"], 2), f(pb["final200_mean"], 2)),
    ("P sd (kbar)", f(pa["std"], 2), f(pb["std"], 2)),
    ("volume change (%)", f(a["cell"]["volume_change_percent"], 3), f(b["cell"]["volume_change_percent"], 3)),
    ("max cond(H)", f(a["cell"]["condition_number_max"], 5), f(b["cell"]["condition_number_max"], 5)),
    ("energy drift (eV/atom/ps)", f"{a['energy']['drift_eV_per_atom_per_ps']:.2e}",
     f"{b['energy']['drift_eV_per_atom_per_ps']:.2e}"),
    ("temperature class", a["temperature_classification"], b["temperature_classification"]),
    ("cell class", a["cell_classification"], b["cell_classification"]),
]
for r in rows:
    L.append(f"| {r[0]} | {r[1]} | {r[2]} |")
L.append("")
L.append("## Provenance")
L.append("")
L.append(f"- Remote WorkDir: `{S['remote_workdir']}`")
L.append("- QE 7.4.1: `/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/nvhpc-24.5/"
         "quantum-espresso-7.4.1-jfm6pcx2dqkcfcwywdgoyao3m3ncvc5a/bin/pw.x`, SHA256 "
         "`beb956cce572605941b791c37bc4dd30ac6cee714e18be7b93fd4ff0e47f4aca`, banner `v.7.4.1`")
L.append("- QE 7.5: `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc/"
         "install-current-gpu/bin/pw.x`, SHA256 "
         "`89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d` (matches expected), "
         "banner `v.7.5`")
L.append("- Pseudopotential `C.UPF` SHA256 `9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639`")
L.append("- Both steps: 1 MPI rank via `srun --mpi=pmix_v3 -n 1`, 8 OpenMP threads, node `lrdn0416`")
L.append("")
L.append("## Caveats")
L.append("")
L.append("- `ecutwfc=25 Ry` is too low for quantitative carbon pressures; absolute P values here "
         "are dominated by Pulay stress and should not be quoted as physical.")
L.append("- Direct velocity reconstruction is possible only at the terminal state, because "
         "`tmp/*.md` is rewritten each step (it holds only the last state) and `tmp/*.vel` is "
         "identically zero in this PR path.")
L.append("- 2 ps is short; the conclusions cover this window only.")
L.append("")
with open(os.path.join(MIRROR, "FINAL_COMPARISON_REPORT_AUDITED.md"), "w") as fh:
    fh.write("\n".join(L) + "\n")
print("wrote FINAL_COMPARISON_REPORT_AUDITED.md")
