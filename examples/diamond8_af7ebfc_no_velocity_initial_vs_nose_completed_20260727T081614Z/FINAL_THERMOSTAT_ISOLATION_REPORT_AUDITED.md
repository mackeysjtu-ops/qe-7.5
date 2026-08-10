# Audited Final Thermostat Isolation Report

Diamond8, QE 7.5 `af7ebfc`, **no explicit velocities**, `caimd=.false.` in both runs. Single changing variable: `ion_temperature='initial'` (CONTROL) vs `'nose'` with `fnosep=50`, `nhpcl=4` (TEST). Job **50167241**.

**Overall status: `EXPLICIT_VELOCITIES_REQUIRED_FOR_HEATING`**

- CONTROL: `CONTROL_BOUNDED` · NOSE: `BOUNDED_BUT_OFFSET` · workflow: `VALID_SINGLE_FACTOR`

## Headline result

**The Nose run does not heat.** With native velocity initialization (no `ATOMIC_VELOCITIES`) the Nose-Hoover trajectory is bit-identical to the control for the first 529 ionic steps and settles onto the same ~170 K plateau; its maximum instantaneous temperature is 300.001 K — exactly the step-0 value — and it never crosses 600 K. This is the opposite of the earlier explicit-velocity Nose run (job 50127909), which reached 8124 K. The Nose thermostat was confirmed active in the TEST stdout (`True`). Therefore the earlier several-thousand-kelvin heating **requires the explicit-velocity input path**, not the Nose thermostat by itself, and not active CAIMD (disabled in every run).

The plateau sits near ~170 K rather than 300 K because the run starts at the ideal lattice (a potential-energy minimum) with all energy kinetic, so equipartition halves the kinetic temperature; at these parameters the Nose chain does not pull it back to 300 K within 2 ps (it barely acts — first divergence from the control only at 1.058 ps). That is a weak-coupling / short-time observation, **not** heating.

## Required answers

1. **Job ID and final Slurm states:** `50167241` — top-level, `batch`, `extern` all `COMPLETED 0:0`; step `.0` (control `pw.x`) and `.1` (nose `pw.x`) both `COMPLETED 0:0`; elapsed 00:03:44 on `lrdn0558`.
2. **Did both QE calculations complete?** Yes. Both printed `JOB DONE`, returned 0, and wrote 1000 `.eal`/`.p`/`.tv`/`.cel` rows. Classifications `CONTROL_BOUNDED` / `BOUNDED_BUT_OFFSET`; SCF failures 0; runtime version `7.5` (both).
3. **Explicit velocities absent from both inputs?** Yes — `True`: no `ATOMIC_VELOCITIES` header, no velocity rows, no `ion_velocities`/`from_input` token in the complete text of either input, both `restart_mode='from_scratch'`.
4. **`caimd=.false.` in both?** Yes — `True` (exactly one `caimd = .false.` line each; no `caimd_stride`/`caimd_step_strain`/`caimd_mask`; no `caimd=.true.`).
5. **All non-thermostat settings identical?** Yes — `True`. Maximum difference 0 for cell, positions, ordering, masses, cutoffs, timestep, `nstep`, `press`, `wmass`, `tempw`, `cell_dynamics`, `cell_dofree`, k-points and every SCF parameter. Only `prefix`, `ion_temperature`, `fnosep(1)` and `nhpcl` differ.
6. **Same native initial velocity realization?** Yes — initial-state class `IDENTICAL`. Step-0 and step-1 ionic kinetic energies are bit-identical and the `.eal` first diverges only at step 529 (1.058 ps). `ranv` zeroes total linear momentum, so both start with zero COM velocity. This is a strict single-factor trajectory comparison.
7. **Did the control reproduce the previous bounded result?** Yes. New control mean 169.858 K, final-200 174.046 K, max 300.001 K, 0 records >600 K, volume change -4.622%. Previous audited reference (job 50160152): mean 170.014 K, final-200 174.827 K, volume -4.622%. Reproduced within small-system variation.
8. **Nose peak temperature and when:** 300.001 K at step 0 (t = 0.000 ps) — i.e. the peak IS the initial value; there is no heating excursion above it.
9. **Nose 100-step block means (K):** [172.1, 168.2, 168.5, 172.0, 170.0, 164.2, 166.9, 168.5, 176.6, 169.9].
10. **Did the temperature decline after the peak?** There is no high-temperature peak to decline from; the rolling mean stays in a ~160-177 K band the whole run (lowest 100-step rolling mean after step 0 = 159.390 K).
11. **Did it return below 2000/1000/600 K?** It never rose above any of them: fractions above 600/1000/2000 K = 0.0/0.0/0.0.
12. **Final-100 / final-200 / final-250 means:** 169.877 / 173.237 / 172.538 K (Nose); 171.471 / 174.046 / 173.185 K (control).
13. **Final-half and final-quarter temperature slopes (Nose):** 9.671 K/ps and -10.550 K/ps (small, sign-alternating — consistent with a bounded plateau, not a trend).
14. **Evidence of convergence toward 300 K?** No secular approach to 300 K within 2 ps; the trajectory is bounded near ~170 K. The thermostat is simply too weakly/slowly coupled at these parameters over this window to reach the target — but there is no divergence.
15. **Did direct `.md` velocities confirm the reported temperature?** Yes. Terminal direct T (from `.md` `ratd` x terminal `.cel`, COM removed): control 143.427 K vs `.eal` 143.436 K; nose 144.614 K vs 144.623 K. Confirmed (control=True, nose=True); no diagnostic mismatch.
16. **Did temperature growth precede cell deformation?** Moot — there was no temperature growth. Nose never exceeded 600 K.
17. **Was the cell stable?** Yes. Both `STABLE_PR_CELL`: volume change control -4.622% / nose -4.622%, max condition number 1.00123 / 1.00123, min C-C 2.824 / 2.824 bohr.
18. **Was the Nose/PR extended energy bounded?** The available quantity `edyn` (`.e` col3 = ekint+ut) is bounded (control drift 1.88e-03 eV/atom/ps, nose 1.82e-03 eV/atom/ps). **This is a PARTIAL diagnostic**: `ions_nose_energy` and the chain variables are computed in `vcsmd` but never written to any output file, so the full Nose+PR conserved Hamiltonian cannot be formed. No secular energy injection is visible.
19. **Are explicit velocities required to reproduce the heating?** On this evidence, **yes** — Nose with native initialization is bounded, while the earlier run that heated differed only by supplying explicit `ATOMIC_VELOCITIES`. This is the direct isolation the task set out to make.
20. **Is active CAIMD required?** No. CAIMD was disabled (`caimd=.false.`) in both runs here and in the earlier heating run.
21. **Correct conclusion — transient convergence, late unresolved recovery, or sustained heating?** **None of those** — there is no heating at all. The native-initialization Nose trajectory is bounded from the start; the earlier heating is attributable to the explicit-velocity path.
22. **Minimum next test:** see below.
23. **Were any source files or raw remote outputs modified?** No. The remote workflow was read once (rsync pull); all new products are under `analysis_audited/`, `evidence_audited/` and this report. No QE source touched, no job submitted/cancelled/polled, no commit/push.

## Comparison table

| quantity | CONTROL (`initial`) | TEST (`nose`) |
|---|---:|---:|
| classification | CONTROL_BOUNDED | BOUNDED_BUT_OFFSET |
| Nose active in stdout | n/a | True |
| steps / JOB DONE | 1000 / True | 1000 / True |
| initialized T (K) | 300.0013 | 300.0013 |
| mean T (K) | 169.858 | 169.696 |
| final-200 mean T (K) | 174.046 | 173.237 |
| max T (K) | 300.001 | 300.001 |
| min T (K) | 55.076 | 55.076 |
| fraction >600 K | 0.0 | 0.0 |
| drift full (K/ps) | 1.065 | 0.587 |
| drift final-quarter (K/ps) | -3.257 | -10.550 |
| mean P (kbar) | -708.75 | -708.75 |
| volume change (%) | -4.622 | -4.622 |
| max cond(H) | 1.00123 | 1.00123 |
| terminal direct-vel T (K) | 143.427 | 144.614 |

## Three-way context

| run | thermostat | velocities | outcome |
|---|---|---|---|
| job 50160152 (control ref) | initial | native | bounded, ~170 K |
| **this, control** | initial | native | bounded, 170 K mean |
| **this, TEST** | **nose** | native | **bounded, 170 K mean, max 300 K** |
| job 50127909 | nose | **explicit** | **heated to 8124 K** |

The only factor that flips the outcome from bounded to heated is the explicit-velocity input path.

## Minimum next test (not submitted)

Because native-initialization Nose is bounded while the previous explicit-velocity Nose heated, move the variable to the velocity path (single factor at a time), keeping `caimd=.false.`: (a) `ion_temperature='initial'` **with** an explicit `ATOMIC_VELOCITIES` card and **no** Nose; then (b) Nose **plus** the same explicit velocities. Audit the first-step ionic kinetic energy, the Nose chain kinetic input `ekin2nhp`, and the Nose chain variables (`xnhp0`, `vnhp`) at initialization when velocities are supplied — that initialization is the remaining structural difference against job 50127909.

## Caveats

- `ecutwfc=25 Ry` retained from the diagnostic deck; absolute pressures are Pulay-dominated and not converged. Pressure/volume are used only for comparison and cell stability.
- The `.tv`/stdout temperature is a running average (`2*avk/(3(N+1)k_B)`); all instantaneous statistics use `2*eka/((3N-3)k_B)` from `.eal`.
- `tmp/*.vel` is identically zero in the PR path and is not used; `tmp/*.md` holds only the terminal state, so direct velocity reconstruction is available at the final frame only.
- The Nose chain conserved-energy terms are not written by QE, so the extended-energy check is a partial diagnostic (documented above).
- 2 ps is short; conclusions cover this window. The result is a clear negative for heating, not a claim that Nose reaches 300 K given more time.

