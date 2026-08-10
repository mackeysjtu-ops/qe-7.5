# Final Comparison Report

Diamond8, native from-scratch 300 K initialization, **no explicit velocities**.

Overall status: `BOTH_COMPARABLE`

## Scope

This comparison tests native 300 K initialization (ion_temperature='initial') followed by ordinary variable-cell Parrinello-Rahman dynamics. It is NOT a continuing Nose-thermostat test: QE 7.4.1 does not accept ion_temperature='nose' in PW vc-md.

No explicit velocities are supplied, so each executable generates its own Maxwell realization via vcinit->ranv. QE's randy() is seeded deterministically (idum=0), so the realizations MAY coincide, but this is not assumed: statistical comparison is used.

## Required answers

1. **Explicit velocities absent from both inputs?** Yes. See `evidence/no_explicit_velocity_audit.json`; neither input contains an ATOMIC_VELOCITIES card or an explicit-velocity keyword, in the complete text (case-insensitive).
2. **Common 300 K initialization mode:** `ion_temperature='initial'` with `tempw=300.0`, accepted by both versions.
3. **Both executables initialized velocities?** QE 7.4.1: yes; QE 7.5: yes (non-zero initial ionic kinetic energy).
4. **Actual initial temperature:** QE 7.4.1 = 300.001 K; QE 7.5 = 300.001 K (3N-3 ionic).
5. **QE 7.4.1 near 300 K?** classification `STABLE_NEAR_300K` (mean 170.090 K, final-200 mean 175.203 K).
6. **QE 7.5 caimd=.false. near 300 K?** classification `STABLE_NEAR_300K` (mean 170.014 K, final-200 mean 174.827 K).
7. **Did either run physically heat?** QE 7.4.1 max 300.001 K (first >2000 K: None); QE 7.5 max 300.001 K (first >2000 K: None).
8. **Mean / final-window temperatures:** QE 7.4.1 170.090 / 175.203 K; QE 7.5 170.014 / 174.827 K.
9. **Mean / final-window pressures:** QE 7.4.1 -712.448 / -642.175 kbar (-71.245 GPa mean); QE 7.5 -712.379 / -641.833 kbar (-71.238 GPa mean).
10. **Pressure differences consistent with 8-atom fluctuations?** True (|Δmean| = 0.068 kbar vs fluctuation scale 48.359 kbar).
11. **Cell instability?** QE 7.4.1 `STABLE_PR_CELL` (final ΔV -4.622 %, cond max 1.0012); QE 7.5 `STABLE_PR_CELL` (final ΔV -4.622 %, cond max 1.0012).
12. **Statistically comparable?** temperature: True; pressure: True; identical temperature series: False.
13. **Narrowest justified next test:** see below.

Both versions behave equivalently under native 300 K initialization. The next test is to reintroduce the QE 7.5 Nose thermostat alone (7.5 vs 7.5, nose vs initial) to isolate whether the previously observed heating is specific to the PW vc-md Nose path rather than to the executable.

## Run summary

| run | version | rc | JOB DONE | steps | T init (K) | T mean (K) | T max (K) | P mean (kbar) | T class | cell class |
|---|---|---:|---|---:|---:|---:|---:|---:|---|---|
| qe74 | 7.4.1 | 0 | True | 1000 | 300.001 | 170.090 | 300.001 | -712.448 | STABLE_NEAR_300K | STABLE_PR_CELL |
| qe75 | 7.5 | 0 | True | 1000 | 300.001 | 170.014 | 300.001 | -712.379 | STABLE_NEAR_300K | STABLE_PR_CELL |

## Caveats

- `sbatch` acceptance is not scientific validation; this report reflects parsed output only.
- Instantaneous pressures in an 8-atom cell fluctuate strongly; only statistical agreement is meaningful.
- QE 7.4.1 and QE 7.5 may generate different Maxwell realizations; exact step-by-step equality is not required unless the identical-series check above reports true.
- First errors, if any: QE 7.4.1 `None`; QE 7.5 `None`.

