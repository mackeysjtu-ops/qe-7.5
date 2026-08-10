# Attempt 01 Failure: `wmass=2000`

Generated: 2026-07-23T10:23:31Z workflow, diagnosed during local-only rerun.

This first equilibration attempt used the Diamond-8 Nose+Parrinello-Rahman
input generated with `wmass=2000.0`, inherited from the short-step local
CAIMD+Nose harness. Quantum ESPRESSO itself completed the requested 500 ionic
steps with return code 0 and `JOB DONE`, but the trajectory is not physically or
numerically acceptable.

Key validation evidence from `VALIDATION.json`:

- `trajectory_blocks`: 500 cell, 500 position, 500 velocity blocks.
- Temperature source: stdout `T=...K`.
- Temperature classification: `FAIL`.
- Temperature max: `8990.9 K`.
- Temperature mean: `3602.85625 K`.
- Final volume: `824403536.2501643 bohr^3`.
- Final minimum periodic C-C distance reported by the local validator:
  `196.38870831656803 bohr`.

Diagnosis: this is a Parrinello-Rahman cell runaway under the required 2 fs
local validation timestep, not a QE source-code defect. The follow-up attempt
therefore preserves the same scientific input model but scales `wmass` from
the local harness `dt=5.0 a.u.` value by `(82.682762012/5.0)^2`, preserving the
approximate `dt^2/wmass` update scale for the larger mandated timestep.
