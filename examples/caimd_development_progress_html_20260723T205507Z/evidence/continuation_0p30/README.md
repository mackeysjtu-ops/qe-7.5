# Local Diamond-8 CAIMD + Nose Continuation To 0.30

This directory contains a completely local macOS validation of the QE 7.5
CAIMD implementation at commit `af7ebfc`.

The loading ensemble is thermostatted constrained-strain Parrinello-Rahman
dynamics. It is not unconstrained full NpT because laboratory `xx` strain is
externally driven by CAIMD.

The loading rate is `0.05 ps^-1`, which is 50 times faster than
`0.001 ps^-1`. It is used only for local functional validation. The target
global engineering compression of about 0.30 is a software stress test, not a
quantitative high-strain stress curve.

The QE input `dt` is `41.341381006000`, validated against QE stdout as a 0.002 ps
per-step increment. The Parrinello-Rahman cell mass is `2000000.00000000`. The
starting estimate from the local harness scaling would be
`136728.78267866`, then this workflow uses the heavier value to keep PR
cell kinetic modes from draining the ionic temperature during the short local
equilibration.

The ionic Nose setup is `fnosep=50.0 THz` with `nhpcl=4`.
The carbon cutoff is `ecutwfc/ecutrho=25.0/120.0 Ry`,
matching the bundled QE Diamond+C.UPF test-suite level.

The Nose thermostat target remains `300.0 K`. Initial velocities
are warm-started at `450.0 K` under the 3N-3
ionic convention because archived zero-displacement Diamond8 starts at 300 K
equilibrated to a cold ionic kinetic temperature near 180 K.

The ordinary Nose+PR equilibration is `500` steps. The archived
attempt directories document the 500-step screens and the stable 1500-step
single-chain extension that motivated the four-link Nose-Hoover chain.

Run commands are recorded in `evidence/commands.log`. Completed stages write
`*_COMPLETE.json` markers and can be skipped by rerunning the workflow.
