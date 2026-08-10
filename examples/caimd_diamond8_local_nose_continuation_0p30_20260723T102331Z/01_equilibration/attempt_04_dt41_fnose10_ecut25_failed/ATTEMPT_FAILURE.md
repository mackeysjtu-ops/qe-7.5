# Attempt 04 Failure: gentler Nose, improved but still hot

This fourth equilibration attempt used:

- `dt=41.341381006`, validated as `0.002 ps` per ionic step.
- `ecutwfc/ecutrho=25/120 Ry`.
- `fnosep=10 THz`.

The trajectory stayed geometrically finite and did not show the earlier cell
runaway. It is still invalid because the ionic temperature drifted above the
accepted local validation range:

- Temperature classification: `FAIL`.
- Temperature mean: `420.9892 K`.
- Temperature max: `901.0 K`.
- Final 200-step mean temperature: `678.6355 K`.
- Final volume: `261.6983068184196 bohr^3`.
- Observed time increment: `0.002 ps`.

Diagnosis: `fnosep=10 THz` is much more stable than 100 THz, but too weak for
the continuing PR stress-release heating in this local Diamond8 setup. The next
attempt increases only the ionic Nose frequency to `20 THz`.
