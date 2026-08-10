# Attempt 03 Failure: correct 2 fs timestep, unstable ionic temperature

This third equilibration attempt corrected the QE timestep to
`dt=41.341381006`, and the validator confirmed an observed stdout increment of
`0.002 ps` per ionic step. The Parrinello-Rahman geometry remained finite and
bounded through 500 steps.

The run is still invalid because the ionic Nose temperature control failed:

- Temperature classification: `FAIL`.
- Parsed stdout temperature max: `9843.3 K`.
- Parsed stdout temperature mean: `8896.8916 K`.
- Final 200-step mean temperature: `9331.734 K`.
- Final volume: `311.2205229265596 bohr^3`.
- Final minimum periodic C-C distance reported by the local validator:
  `2.0018905642932054 bohr`.
- Observed time increment: `0.002 ps`.

Diagnosis: this is no longer a PR cell runaway. QE printed only 9 time steps per
Nose oscillation for `fnosep=100 THz` at the corrected 2 fs step. The follow-up
attempt therefore uses a gentler `fnosep=10 THz` and the QE diamond test-suite
carbon cutoff level `ecutwfc/ecutrho=25/120 Ry` instead of the very cheap
`12/48 Ry` screen value.
