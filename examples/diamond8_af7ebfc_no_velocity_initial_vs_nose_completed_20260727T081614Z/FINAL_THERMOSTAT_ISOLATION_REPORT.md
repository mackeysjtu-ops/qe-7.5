# Final Thermostat Isolation Report

Diamond8, QE 7.5 `af7ebfc`, **no explicit velocities**, `caimd=.false.` in both runs.
Single changing variable: `ion_temperature='initial'` (CONTROL) vs `'nose'` with `fnosep=50`, `nhpcl=4` (TEST).

**Overall status: `NOSE_STABLE_WITH_NATIVE_INITIALIZATION`**

- CONTROL classification: `CONTROL_BOUNDED`
- NOSE classification: `BOUNDED_BUT_OFFSET`

## Required answers

1. **Explicit velocities absent from both inputs?** Yes — see `evidence/no_explicit_velocity_audit.json`: no `ATOMIC_VELOCITIES` header, no velocity rows, no explicit-velocity keyword anywhere in the complete text, no restart-state input, both `restart_mode='from_scratch'`.
2. **`caimd=.false.` in both?** Yes — `evidence/caimd_disabled_audit.json`: exactly one `caimd = .false.` line each, no `caimd_stride`/`caimd_step_strain`/`caimd_mask`, no `caimd=.true.`.
3. **All non-thermostat settings identical?** Yes — `evidence/input_identity_audit.json` reports maximum difference 0 for cell, positions, ordering, masses, cutoffs, timestep, `nstep`, `press`, `wmass`, `tempw` and every SCF parameter. Only `prefix`, `ion_temperature`, `fnosep(1)` and `nhpcl` differ.
4. **Same native initial velocity realization?** Not identical — first-step ionic kinetic energy 0.01995100 Ry (control) vs 0.01995100 Ry (nose), |diff| = 0.0; step-2 position max |diff| = 2.9993368899638995e-15 bohr. Initial realizations differ: the random realization is an additional uncontrolled difference; compare statistically.
5. **Did the control remain bounded?** `CONTROL_BOUNDED` — mean 169.858 K, final-200 174.046 K, max 300.001 K, min 55.076 K, drift 1.065 K/ps.
6. **Did Nose produce an early overshoot?** max 300.001 K; first >600 K at None.
7. **Did the rolling/block mean later return toward 300 K?** final-200 mean 173.237 K; 100-step rolling mean at the end 169.877 K; late-half block slope 2.101 K per 100-step block (negative = declining).
8. **Did Nose cross 600 / 1000 / 2000 K?** 600 K: None; 1000 K: None; 2000 K: None.
9. **How long did high-temperature intervals persist?** longest consecutive run above 600 K = 0 steps, above 1000 K = 0, above 2000 K = 0; fractions {'600K': 0.0, '1000K': 0.0, '2000K': 0.0}.
10. **Did direct `.md` velocities confirm the reported temperature?** control: direct 143.427 K vs `.eal` 143.436 K; nose: direct 144.614 K vs 144.623 K. Confirmed: control=True, nose=True.
11. **Did heating precede significant cell deformation?** {'first_600K_time_ps': None, 'note': 'Nose never crossed 600 K'}
12. **Was the extended-energy behavior bounded?** control edyn drift 1.880e-03 eV/atom/ps (bounded=True); nose 1.824e-03 eV/atom/ps (bounded=True). **This is a PARTIAL quantity**: `ions_nose_energy` is computed in `vcsmd` but never written to an output file, so the full Nose+PR conserved Hamiltonian cannot be formed; a change in `edyn` for the Nose run is expected by construction.
13. **Are explicit velocities required to reproduce the earlier heating?** On this evidence, **plausibly yes** — Nose alone with native initialization did not reproduce sustained heating, so the explicit-velocity path becomes the prime remaining suspect. This is an inference from two runs, not yet a proof.
14. **Is active CAIMD required?** No — CAIMD was disabled in both runs here, and in the earlier heating run as well.
15. **Narrowest justified next test:** see the recommendation below.

## Run comparison

| quantity | CONTROL (`initial`) | TEST (`nose`) |
|---|---:|---:|
| classification | CONTROL_BOUNDED | BOUNDED_BUT_OFFSET |
| return code | 0 | 0 |
| JOB DONE | True | True |
| steps | 1000 | 1000 |
| initialized T (K) | 300.0013 | 300.0013 |
| mean T (K) | 169.858 | 169.696 |
| final-200 mean T (K) | 174.046 | 173.237 |
| max T (K) | 300.001 | 300.001 |
| min T (K) | 55.076 | 55.076 |
| T drift (K/ps) | 1.065 | 0.587 |
| mean P (kbar) | -708.75 | -708.75 |
| final-200 P (kbar) | -638.19 | -638.23 |
| volume change (%) | -4.622 | -4.622 |
| max cond(H) | 1.00123 | 1.00123 |
| min C-C (bohr) | 2.824 | 2.824 |

## Comparison with previously verified results

- Previous `initial`/no-velocity control (audited): mean 170.014 K, final-200 174.827 K, volume change -4.622%. New control: mean 169.858 K, final-200 174.046 K, volume -4.622%.
- Previous Nose **with explicit velocities** (job 50127909): max 8123.839 K, final 4231.686 K, first >600 K at step 26. New Nose **without** explicit velocities: max 300.001 K, final 144.623 K.

## Safety

No external runtime killer was armed: each run is ~2 minutes of wall time and QE's own singular-cell / non-finite guards are the effective safety mechanism. The extreme conditions were evaluated post-hoc:

- CONTROL: {'nonpositive_determinant': False, 'condition_number_gt_1e8': False, 'nonfinite_values': False, 'qe_singular_cell': False, 'T_above_20000K_for_10_consecutive': False, 'any_triggered': False}
- TEST: {'nonpositive_determinant': False, 'condition_number_gt_1e8': False, 'nonfinite_values': False, 'qe_singular_cell': False, 'T_above_20000K_for_10_consecutive': False, 'any_triggered': False}

## Recommended minimum next test (not submitted)

Nose is not sufficient on its own here, so move the variable to the velocity path: (a) `ion_temperature='initial'` **with** an explicit `ATOMIC_VELOCITIES` card and no Nose; then (b) Nose **plus** explicit velocities. Inspect how the Nose chain is initialized when velocities are supplied (`vcinit` velocity branch vs `ranv`), since that is the remaining difference against job 50127909.

## Caveats

- `ecutwfc=25 Ry` is retained deliberately from the reference deck; absolute pressures are dominated by Pulay stress and are **not** converged diamond pressures. Pressure and volume are used only to compare the two runs and to judge cell stability.
- The reported `.tv`/stdout temperature is a **running average**; all instantaneous statistics above use `2*eka/((3N-3) k_B)` from `.eal`.
- `tmp/*.vel` is not used: it is identically zero in this PR path.
- `tmp/*.md` holds only the terminal state, so direct velocity reconstruction is available at the final frame only.
- `sbatch` acceptance is not scientific validation; this report reflects parsed output.

