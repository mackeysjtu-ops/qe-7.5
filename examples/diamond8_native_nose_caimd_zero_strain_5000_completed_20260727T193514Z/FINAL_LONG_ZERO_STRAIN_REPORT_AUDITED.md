# Audited Final Report: Diamond8 5000-step native Nose + CAIMD zero-strain

Job **50434898**, QE 7.5 `af7ebfc`. Three 10 ps native-initialization Nose trajectories, no explicit velocities: CASE A `caimd=.false.`, CASE B enabled no-op, CASE C active all-component zero-strain hold.

**Overall status: `ACTIVE_CAIMD_CHANGES_TEMPERATURE_RESPONSE`**

| case | T class | mean K | final-5ps K | final-2.5ps K | max K | late-half block K +/- SEM | >600K frac |
|---|---|---:|---:|---:|---:|---:|---:|
| A | SUSTAINED_HEATING | 2844.6 | 4906.0 | 4786.4 | 8981.7 | 4894.7 +/- 73.5 | 0.566 |
| B | SUSTAINED_HEATING | 2833.8 | 4966.9 | 4957.5 | 8384.4 | 4948.0 +/- 48.0 | 0.558 |
| C | SUSTAINED_HEATING | 384.1 | 523.0 | 599.7 | 1185.1 | 547.1 +/- 68.2 | 0.157 |

## Headline

All three cases COMPLETE 5000 steps. CASE A (control) and CASE B (no-op) are statistically identical: both sit ~170 K for ~2-3 ps then run away to ~5000 K as the free Parrinello-Rahman cell contracts ~9%. CASE C (active zero-strain) holds the cell perfectly rigid (V=306.27 bohr^3, cond=1.0 constant) and heats far less (final-5 ps ~520 K). The ionic Nose thermostat is active every step with the correct target (ndega=3N-3=21) but does NOT hold 300 K in any case; the earlier '~170 K bounded' picture was an artifact of the 2 ps window - the heating onset is at ~3-4 ps.

This overturns the earlier 2 ps picture (job 50167241) that both control and Nose stayed bounded near 170 K: the heating onset is at ~3-4 ps, beyond that window. The three 1000-step block means for CASE A are [170, 221, 4063, 5031, 4738] K and for CASE C [185, 262, 365, 490, 619] K.

## The 300 K question (three separate answers)

**A. Are instantaneous fluctuations reasonable for 8 atoms?** Yes. Large step-to-step swings are expected for a 3N-3 = 21 DOF system; the instantaneous standard deviation is not itself a defect.
**B. Is the late-time mean statistically compatible with 300 K?** **No, in every case.** CASE A final-half block mean 4894.7 +/- 73.5 K, CASE B 4948.0 +/- 48.0 K, CASE C 547.1 +/- 68.2 K. None overlaps 300 K within 2 SEM.
**C. Is the thermostat scientifically validated for future CAIMD loading?** No. The ionic Nose is active every step with the correct target but does not maintain 300 K.

**Plain-language verdict (温度是否合理):** the fluctuations are expected for 8 atoms, but the ionic subsystem is not held at the requested 300 K in any case. CASE A/B run away to ~4906 K driven by the free Parrinello-Rahman cell (which contracts ~9%); CASE C, with the cell frozen by zero-strain CAIMD, is far milder (~523 K) but still above 300 K and rising. The trajectories are numerically stable but thermodynamically at the wrong temperature.

## Required answers

1. **All three complete 5000 steps?** Yes - A=5000(COMPLETE_5000), B=5000(COMPLETE_5000), C=5000(COMPLETE_5000).
2. **Explicit velocities absent?** True - no ATOMIC_VELOCITIES/ion_velocities/from_input in any input's complete text; all from_scratch.
3. **Initial velocity states identical?** Native realization is the same (first-step eka = 0.019951 Ry in all three). A vs B: `IDENTICAL` (bit-identical to step 701). A/B vs C: `NUMERICALLY_EQUIVALENT` - C diverges at step 10 because its active projection acts from the start.
4. **Authoritative temperature quantity:** instantaneous ionic kinetic temperature `T = 2*eka/((3N-3)k_B)` with 3N-3 = 21, `eka` from `.eal` column 2. This matches `ndega = get_ndof() = 21`, the Nose target DOF. The `.tv`/stdout temperature is a running average and is NOT used. Terminal direct `.md` velocity reconstruction confirms `.eal` (A within 0.49 K).
5. **Full-run means:** A 2844.6, B 2833.8, C 384.1 K.
6. **Final-5 ps means:** A 4906.0, B 4966.9, C 523.0 K.
7. **Final-2.5 ps means:** A 4786.4, B 4957.5, C 599.7 K.
8. **1000-step block means (K):** A [170, 221, 4063, 5031, 4738]; B [170, 221, 3832, 4999, 4947]; C [185, 262, 365, 490, 619].
9. **Block-based uncertainties (final-half, SEM):** A 73.5, B 48.0, C 68.2 K.
10. **Late-time means compatible with 300 K?** No for all three (see answer B above).
11. **Instantaneous fluctuations reasonable for 8 atoms?** Yes (std ~2472-191 K); large swings are expected and are not the problem.
12. **Factor-of-two / DOF issue?** No. `ndega = get_ndof() = 3N-3 = 21` equals the measurement divisor, so a working thermostat would read ~300 K. No divisor choice (21/24/27) reconciles the observed eka with 300 K - the temperatures are genuinely off-target, not a definitional artifact. (This is a real change from the earlier ~170 K discussion: at 10 ps the system is HOT, not cold.)
13. **Is Nose applied every step?** Yes. `tnosep=.true.` throughout; the `vcmove` Nose friction term acts unconditionally each step; the velocity-rescaling branch is gated `.and. .not. tnosep`, so Nose is not gated behind `ntcheck`.
14. **Does Nose exchange measurable energy with the ions?** Indirectly yes (the trajectory is driven away from a pure microcanonical path), but the Nose chain energy is not written by QE, so a direct energy-exchange measurement is unavailable (evidence limit stated).
15. **Active but ineffective, or inactive?** `CONFIRMED_ACTIVE_BUT_INEFFECTIVE` - active with the correct target but ineffective: it neither warms the cold early ions to 300 K nor prevents the later PR-cell-driven runaway.
16. **Is CASE B a valid long-time no-op?** Yes - `LONG_TIME_STATISTICAL_NOOP`; bit-identical to CASE A until step 701, then statistically equivalent (both heat to ~5000 K); no CAIMD event fired (`any_event=False`).
17. **Is CASE C stable under active zero-strain projection?** Yes - `ACTIVE_ZERO_STRAIN_STABLE`: cell held perfectly rigid (V=306.27 bohr^3, cond=1.0000 constant), determinant always positive, max skew norm 2.1e-16.
18. **CASE C maximum strain error:** 2.22e-16 (machine precision); max velocity-transform residual 0.00e+00.
19. **Does active CAIMD change temperature behavior?** **Yes, materially.** A and B heat to ~4906/4967 K; C only to ~523 K. Freezing the cell (all six strains held at zero) removes the dominant PR-cell energy-injection channel and suppresses the heating ~8x.
20. **Suitable for a nonzero-strain test?** Not yet. The zero-strain projection is numerically flawless, but the underlying ordinary `vc-md + Nose` baseline does not hold 300 K, so a loading test would run on an uncontrolled thermal background.
21. **Narrowest justified next step:** see recommendation below.
22. **Sources or raw outputs modified?** No. Remote workflow read once (rsync pull); all new products under `analysis_audited/`, `evidence_audited/` and this report. No QE source touched, no job submitted/cancelled/polled, no commit/push.

## Why A/B heat and C does not (mechanism)

The heating is Parrinello-Rahman cell-driven. In A/B the free cell contracts ~9% (V 306->280 bohr^3, min volume ~226) and pumps energy into the ions faster than the ionic Nose (which controls ions only; `cell_temperature='not_controlled'`) can remove it. CASE C holds all six symmetric strain components at zero, so the cell is rigid and that energy channel is closed; the residual mild heating (~520 K) is what the ionic Nose alone leaves uncontrolled.

## Minimum next step (recommendation)

Do **not** proceed to nonzero-strain scientific loading. The immediate diagnostic is to isolate the thermostat from the barostat: run a **fixed-cell `calculation='md'` + Nose** test (same native init, same `fnosep=50`, `nhpcl=4`, `dt=2 fs`, no PR cell, no CAIMD) to check whether the ionic Nose alone holds 300 K when there is no PR-cell energy source. In parallel, a conservative ordinary `vc-md + Nose` parameter check is warranted (the `wmass=2e6` PR cell mass and 2 fs step may make the PR cell dynamics the instability): try lower `fnosep`, shorter `dt`, and a cell thermostat / different `wmass`. Only after a stable 300 K ordinary baseline is established should the zero-strain and then small nonzero-strain CAIMD tests be repeated.

## Provenance & caveats

- Job 50434898, node lrdn1374; steps .0/.1/.2 all `COMPLETED 0:0`; A/B/C RC=0, `JOB DONE`.
- Executable SHA256 89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d (matches expected); banner v.7.5.
- `ecutwfc=25 Ry`: absolute pressures are Pulay-dominated, not converged; used only for comparison and stability.
- `.tv`/stdout temperature is a running average; instantaneous stats use `2*eka/((3N-3)k_B)`.
- Nose chain conserved-energy terms are not written by QE; the extended-energy check is partial.

