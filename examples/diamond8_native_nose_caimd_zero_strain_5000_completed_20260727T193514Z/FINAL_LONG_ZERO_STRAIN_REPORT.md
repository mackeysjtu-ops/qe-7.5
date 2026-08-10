# Final Long Zero-Strain CAIMD Report

Diamond8, QE 7.5 `af7ebfc`, three 5000-step (10 ps) native-initialization Nose trajectories. No explicit velocities in any case.

**Overall status: `ACTIVE_CAIMD_ZERO_STRAIN_UNSTABLE`**

| case | CAIMD | T class | mean K | final-5ps K | final-2.5ps K | max K | >600K frac |
|---|---|---|---:|---:|---:|---:|---:|
| A | .false. | STABLE_ABOVE_TARGET | 2844.6 | 4906.0 | 4786.4 | 8981.7 | 0.5662 |
| B | .true. | STABLE_ABOVE_TARGET | 2833.8 | 4966.9 | 4957.5 | 8384.4 | 0.5576 |
| C | .true. | STABLE_ABOVE_TARGET | 384.1 | 523.0 | 599.7 | 1185.1 | 0.1566 |

## Required answers

1. **All three completed 5000 steps?** A=5000(COMPLETE) B=5000(COMPLETE) C=5000(COMPLETE).
2. **Explicit velocities absent in all cases?** True; see `analysis/no_explicit_velocity_audit.json` (no ATOMIC_VELOCITIES / ion_velocities / from_input in the complete text of any input).
3. **Native velocity initialization?** Yes — all three omit velocities and run `vcinit->ranv(300 K)`.
4. **Identical initial velocity realizations?** all-identical=False; pairs A_vs_B=IDENTICAL, A_vs_C=NUMERICALLY_EQUIVALENT, B_vs_C=NUMERICALLY_EQUIVALENT.
5. **CASE A approach 300 K?** `STABLE_ABOVE_TARGET` (mean 2844.6 K, final-5ps 4906.0 K).
6. **CASE B approach 300 K?** `STABLE_ABOVE_TARGET` (mean 2833.8 K, final-5ps 4966.9 K).
7. **CASE C approach 300 K?** `STABLE_ABOVE_TARGET` (mean 384.1 K, final-5ps 523.0 K).
8. **final-5-ps / final-2.5-ps mean T:** A 4906.0/4786.4; B 4966.9/4957.5; C 523.0/599.7 K.
9. **1000-step block means (K):** A [169.9, 221.3, 4062.9, 5031.1, 4737.6]; B [170.0, 221.3, 3831.8, 4999.2, 4946.5]; C [184.7, 261.9, 364.9, 490.2, 618.9].
10. **Late-time slopes consistent with a plateau?** final-quarter slope A -64.674, B -105.739, C 87.239 K/ps.
11. **Any case cross 600/1000/2000 K?** A: 600={'nst': 2045, 'time_ps': 4.09, 'value': 601.2807579414463}; B: 600={'nst': 2085, 'time_ps': 4.17, 'value': 612.5885112105578}; C: 600={'nst': 2623, 'time_ps': 5.246, 'value': 601.0852781641877}.
12. **CASE B a long-time no-op vs A?** consistent=True, mean-T delta 10.777 K, any event flag False, eka first divergence step 701.
13. **CASE C keep all six strains at zero?** nominal all zero=True; max |actual| strain from .cel 0.000.
14. **CASE C max strain-tracking error:** 2.220446049250313e-16.
15. **CASE C determinant always positive?** True (range [306.27103783218973, 306.27103783218973] bohr^3).
16. **CAIMD velocity-transform residuals bounded?** CASE C max velocity_jump_max 0.0, max skew norm 2.1464554419741031e-16.
17. **Did active CAIMD alter temperature control?** mean-T C vs A delta 2460.444 K; classifications A/B/C = STABLE_ABOVE_TARGET/STABLE_ABOVE_TARGET/STABLE_ABOVE_TARGET.
18. **Nose/PR energy bounded?** partial edyn drift (eV/atom/ps): A 2.02e-02, B 1.85e-02, C 1.27e-02; bounded A/B/C=True/True/True. PARTIAL: Nose chain energy is not written by QE.
19. **Stable-below / approaching / at target?** ACTIVE_CAIMD_ZERO_STRAIN_UNSTABLE.
20. **Minimum next step:** see below.
21. **Sources or previous raw outputs modified?** No — this job only wrote new files in its own workflow directory.

## Minimum next step (recommendation)

CASE C is unstable: diagnose the instability (cell singularity, residual growth) before any loading test.

## Caveats

- `ecutwfc=25 Ry`: absolute pressures are Pulay-dominated, not converged; pressure/volume are used only for software and stability diagnostics.
- The `.tv`/stdout temperature is a running average; all instantaneous statistics use `2*eka/((3N-3)k_B)` from `.eal`.
- `tmp/*.vel` is zero in the PR path and not used; `tmp/*.md` holds only the terminal state, so direct velocity reconstruction is at the final frame only.
- Nose chain conserved-energy terms are not written by QE; the extended-energy check is a partial diagnostic.

