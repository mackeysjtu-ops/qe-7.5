# Audited Final Comparison Report

Diamond8, native from-scratch 300 K initialization, **no explicit velocities**: official QE 7.4.1 vs QE 7.5 `af7ebfc` with `caimd=.false.`

**Overall status: `BOTH_COMPARABLE`**

Ensemble label: ordinary Parrinello-Rahman variable-cell dynamics after native 300 K velocity initialization; NOT NVT - no continuing thermostat is active.

## Headline result

The two executables produced **bit-identical dynamics for the first 701 ionic steps** (1.402 ps). They start from the identical Maxwell realization (QE's `randy` is seeded deterministically with `idum=0`), and only diverge afterwards through ordinary floating-point round-off amplified by MD chaos. All statistical measures agree far inside the 8-atom fluctuation scale. **Neither run heats.**

This also shows the severe heating seen previously with `af7ebfc` (max 8123.8 K, job 50127909) is **not** a property of the af7ebfc executable's ordinary `vc-md` path as such. That earlier run differed in *two* ways — it used `ion_temperature='nose'` **and** explicit `ATOMIC_VELOCITIES` — so this test narrows the cause to those two factors jointly, without yet separating them.

## 1-20 required answers

1. **Job ID / Slurm states:** `50160152`; top-level `COMPLETED 0:0`, `batch` and `extern` `COMPLETED 0:0`, step `.0` (QE 7.4.1 `pw.x`) `COMPLETED 0:0`, step `.1` (QE 7.5 `pw.x`) `COMPLETED 0:0`. Elapsed 00:04:32 on `lrdn0416`.
2. **Did both QE calculations actually run?** Yes. Both printed `JOB DONE`, returned 0 (`RC_QE74=0`, `RC_QE75_CAIMD_FALSE=0`), and each wrote 1000 `.eal`/`.p`/`.tv`/`.cel` rows. Classification `COMPLETE` / `COMPLETE`. SCF failures: 0. Banners: QE `7.4.1` and QE `7.5`.
3. **Explicit input velocities absent?** Yes — `True`. No `ATOMIC_VELOCITIES` card and no `ion_velocities`/`from_input` token in the complete text of either input (case-insensitive); both `restart_mode='from_scratch'`.
4. **Native 300 K initialization mode:** `ion_temperature='initial'`, `tempw=300.0`, identical in both inputs. Velocities are generated inside `vcinit` by `ranv(temp=tempw)` (Maxwell), because `input_velocities_read` is false.
5. **Was a continuing thermostat active?** **No.** `input.f90` sets `ntcheck=nstep+1` for `calculation='vc-md'`, and only `ion_temperature='rescaling'` overrides it, so the `vcmove` velocity-rescaling branch never fires; `tnosep=.false.`. The correct label is ordinary PR variable-cell dynamics after a 300 K velocity initialization — **not NVT**.
6. **Identical physical settings?** Yes — `True`. Cell, positions, ordering, species/masses, pseudopotential hash, cutoffs, k-points, SCF settings, timestep, `nstep`, `press`, `wmass`, `cell_dynamics`, `cell_dofree` and the temperature mode all match; the only namelist differences are `prefix` and the QE 7.5-only `caimd` line.
7. **QE 7.5 used only `caimd=.false.`?** Yes — `True` (exactly one CAIMD line; no `caimd_stride`/`caimd_mask`/`caimd_step_strain`; no `.true.`).
8. **QE 7.4.1 contained no CAIMD keyword?** Yes — `True` (zero occurrences in the effective input).
9. **Initialized temperature:** QE 7.4.1 = 300.0013 K, QE 7.5 = 300.0013 K (3N-3). These are **identical to all printed digits**, and the initial `eka` is bit-identical (1.995100e-02 Ry), confirming the same random realization.
10. **Temperature statistics (instantaneous, 3N-3):** QE 7.4.1 mean 170.090 K, sd 40.964 K, median 172.376 K, min 55.076 K, max 300.001 K, final 121.241 K, first-100 172.096 K, final-100 173.839 K, final-200 175.203 K, drift 1.755 K/ps. QE 7.5 mean 170.014 K, sd 40.487 K, median 172.338 K, min 55.076 K, max 300.001 K, final 141.774 K, final-200 174.827 K, drift 1.530 K/ps.
11. **Sustained physical heating?** **No, in neither run.** The maximum temperature equals the *initial* value (300.0013 K) in both runs — the instantaneous temperature never exceeded its starting value. Fraction of steps above 600 K = 0.000/0.000; longest consecutive run above 600 K = 0/0 steps; nothing crosses 1000 K or 2000 K. Classification `BOUNDED_BUT_OFFSET` for both.
    The settling from 300 K to a ~170-175 K plateau is the **expected equipartition result**: the run starts at the ideal lattice (a potential-energy minimum) with all energy kinetic, so roughly half converts to potential energy and the kinetic temperature equilibrates near T_init/2. It is physically correct behaviour, not cooling failure — hence `BOUNDED_BUT_OFFSET` rather than `PHYSICAL_COOLING`.
12. **Did direct velocities confirm the reported temperature?** Yes. Reconstructing K = sum 1/2 m|v - v_COM|^2 from the terminal `.md` `ratd` transformed by the terminal cell gives 121.106 K (QE 7.4.1) and 141.738 K (QE 7.5), against 121.241 K and 141.774 K from `.eal` — agreement to 0.14 K and 0.04 K. No `DIAGNOSTIC_MISMATCH`.
    Two caveats found while auditing the diagnostics: `tmp/*.vel` is **all zeros** for PR `vc-md` and must not be used as a velocity source; and the temperature printed in `.tv` and in QE stdout is a **running average** (`2*avk/(3(N+1)k_B)`, verified to 0.011 K) — not the instantaneous temperature, from which it differs by up to 327 K late in the run.
13. **Pressure statistics:** QE 7.4.1 mean -708.77 kbar (-70.88 GPa), sd 48.81 kbar, median -724.54 kbar, min -793.16, max -540.10, final-200 -638.30 kbar (-63.83 GPa), drift 7.008 GPa/ps. QE 7.5 mean -708.70 kbar (-70.87 GPa), sd 48.64 kbar, final-200 -637.96 kbar, drift 7.024 GPa/ps.
    The strongly negative pressure is a **setup artifact, not a version effect**: `ecutwfc=25 Ry` is low for carbon, so Pulay stress makes the computed pressure large and negative at the ideal lattice constant. The PR cell responds by contracting (-4.62%), and P rises toward the 0 target (-714 -> -638 kbar). Both versions do this identically.
14. **Pressure differences compatible with an 8-atom cell?** Yes. The difference in mean pressure is 0.067 kbar against a per-step fluctuation scale of 48.7 kbar — a ratio of 0.00138. Block means overlap completely (see `pressure_block100_qe74_vs_qe75.png`).
15. **Cell drift or instability?** Both cells are `STABLE_PR_CELL`. Volume 306.27 -> 292.11 bohr^3 (-4.622%) vs -4.622% — differing by only 8.39e-05 percentage points. Max condition number 1.00124/1.00123; determinant stays positive; final strain is near-isotropic (-0.01595, -0.01580, -0.01520) with shears below 5e-4. Minimum C-C distance 2.824 bohr (1.495 A) — a healthy diamond bond.
16. **Energy conservation:** `GOOD_CONSERVATION` for both. Using `edyn` (`.e` column 3 = `ekint + ut`, the extended PR total at `press=0`, semantics verified in `PW/src/vcsmd.f90`): drift 2.43e-03 Ry over 2 ps = 2.07e-03 eV/atom/ps (QE 7.4.1) and 1.92e-03 eV/atom/ps (QE 7.5). The temperature settling is accompanied by the corresponding potential-energy change with a conserved extended total, i.e. **energy redistribution, not energy injection**.
17. **Statistically comparable?** Yes, and far more strongly than statistics require. Identical for the first 701 steps; mean-T difference 0.0752 K and final-200 difference 0.3756 K against a fluctuation scale of 40.7 K; mean-P difference 0.067 kbar against 48.7 kbar. Same temperature class, same cell class, same energy class.
18. **Narrowest justified conclusion:** For an ideal 8-atom diamond cell with native 300 K initialization, no explicit velocities and no continuing thermostat, official QE 7.4.1 and QE 7.5 `af7ebfc` with `caimd=.false.` are **dynamically equivalent** — identical to round-off for 701 steps and statistically indistinguishable thereafter. The af7ebfc tree therefore introduces **no detectable regression in the ordinary `vc-md`/PR path** under these conditions. This does *not* by itself explain the earlier heating, because that run changed both the thermostat (`nose`) and the velocity source (explicit `ATOMIC_VELOCITIES`) at once.
    One real, identified version difference was found, and it is diagnostic only: af7ebfc adds a `tnew` initialization after `vcinit` in `PW/src/vcsmd.f90` (absent from pristine `a11913a`), so QE 7.5 reports `2*ekint/(3(N+1)k_B)` = 233.33 K in the first `.tv`/stdout row where QE 7.4.1 reports the uninitialized 0.0. The trajectories are bit-identical regardless.
19. **Minimum next test (not submitted):** Since both versions are comparable and bounded, the next step is the one-variable-at-a-time separation of the earlier heating: run QE 7.5 `af7ebfc` against **pristine unmodified QE 7.5 (`a11913a`)** on this same no-velocity `initial` input to close the version question, then change **one** factor at a time — first `ion_temperature='nose'` with native initialization (no velocity card), then explicit `ATOMIC_VELOCITIES` with `initial` — to determine which of the two produced the 8123 K run. A raised `ecutwfc` (>=60 Ry) is advisable before drawing any conclusion from absolute pressures.
20. **Were any raw outputs, source files, or remote calculations modified?** No. The remote workflow was only read (one rsync pull). All new products were written under `analysis_audited/`, `evidence_audited/` and this report. No QE source file was touched, no job was submitted/cancelled/polled, and no Git commit or push was made.

## Comparison table

| quantity | QE 7.4.1 official | QE 7.5 af7ebfc `caimd=.false.` |
|---|---:|---:|
| completion | COMPLETE | COMPLETE |
| runtime version | 7.4.1 | 7.5 |
| initialized T (K) | 300.0013 | 300.0013 |
| mean T (K) | 170.090 | 170.014 |
| final-200 mean T (K) | 175.203 | 174.827 |
| max T (K) | 300.0013 | 300.0013 |
| T drift (K/ps) | 1.755 | 1.530 |
| first >600 K | never | never |
| first >2000 K | never | never |
| mean P (kbar) | -708.77 | -708.70 |
| final-200 mean P (kbar) | -638.30 | -637.96 |
| P sd (kbar) | 48.81 | 48.64 |
| volume change (%) | -4.622 | -4.622 |
| max cond(H) | 1.00124 | 1.00123 |
| energy drift (eV/atom/ps) | 2.07e-03 | 1.92e-03 |
| temperature class | BOUNDED_BUT_OFFSET | BOUNDED_BUT_OFFSET |
| cell class | STABLE_PR_CELL | STABLE_PR_CELL |

## Provenance

- Remote WorkDir: `/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_no_velocity_qe74_vs_qe75_caimd_false_300K_20260724T092828Z`
- QE 7.4.1: `/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/nvhpc-24.5/quantum-espresso-7.4.1-jfm6pcx2dqkcfcwywdgoyao3m3ncvc5a/bin/pw.x`, SHA256 `beb956cce572605941b791c37bc4dd30ac6cee714e18be7b93fd4ff0e47f4aca`, banner `v.7.4.1`
- QE 7.5: `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x`, SHA256 `89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d` (matches expected), banner `v.7.5`
- Pseudopotential `C.UPF` SHA256 `9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639`
- Both steps: 1 MPI rank via `srun --mpi=pmix_v3 -n 1`, 8 OpenMP threads, node `lrdn0416`

## Caveats

- `ecutwfc=25 Ry` is too low for quantitative carbon pressures; absolute P values here are dominated by Pulay stress and should not be quoted as physical.
- Direct velocity reconstruction is possible only at the terminal state, because `tmp/*.md` is rewritten each step (it holds only the last state) and `tmp/*.vel` is identically zero in this PR path.
- 2 ps is short; the conclusions cover this window only.

