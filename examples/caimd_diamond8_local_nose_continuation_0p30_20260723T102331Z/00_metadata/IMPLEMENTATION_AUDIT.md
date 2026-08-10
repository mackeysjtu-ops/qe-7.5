# Implementation Audit

Generated: 2026-07-23T11:19:07Z

Scope: local-only audit for the Diamond-8 CAIMD + Nose continuation validation.

## Confirmed Source Behavior

- Active CAIMD allows `ion_temperature='nose'` and `ion_temperature='not_controlled'`; unsupported ionic thermostats remain rejected by the active-CAIMD allowlist documented in `CAIMD_IMPLEMENTATION.md`.
- PW `calculation='vc-md'` with `cell_dynamics='pr'` maps to `calc='cd'` in `PW/src/input.f90`, so variable-cell MD enters `PW/src/vcsmd.f90` and `PW/src/vcsubs.f90`.
- `PW/src/input.f90:1144-1148` sets `tnosep=.true.` for `ion_temperature='nose'`; `PW/src/input.f90:2017-2020` initializes ionic Nose variables with `get_ndof()`.
- `PW/src/dynamics_module.f90:2517-2524` gives the no-fixed-coordinate ionic degree count as `3*nat - 3`.
- `PW/src/vcsmd.f90:334-354` applies the existing PW `vc-md` ionic Nose update around `vcmove`.
- `PW/src/vcsubs.f90:708-715` applies the Nose friction term to atomic velocity histories, not to CAIMD-imposed strain.
- Active CAIMD diagnostic temperature is computed from ionic kinetic energy, not atomic-plus-cell kinetic energy: `PW/src/caimd_module.f90:201-213`.
- `PW/src/vcsubs.f90:1143-1151` refreshes projected ionic kinetic quantities after final CAIMD projection; `vc_refresh_projected_nose_ekin` refreshes Nose chain kinetic inputs from the committed projected `ratd`.
- `PW/src/caimd_module.f90:321-349` transforms `ratd`, `rat2d`, and `rat2di` by `inverse(H_projected) * H_candidate`, preserving physical Cartesian ionic velocity across CAIMD cell projection.
- `.caimd` rows are premove diagnostics. The `temperature_K` column reports the last committed ionic kinetic temperature when available; the first row may fall back to the target/input temperature.
- Native CAIMD restart remains unsupported. `caimd_configure` rejects restart mode, and `vcsmd` rejects existing MD history for active CAIMD at start of run.
- Segmented continuation in this workflow is therefore controlled external continuation with `restart_mode='from_scratch'`, predecessor terminal cell, predecessor fractional positions, and predecessor physical Cartesian velocities. It resets PR history, Nose extended variables, CAIMD event count, and the local CAIMD reference cell.

## Loading Ensemble Name

The loading runs are thermostatted constrained-strain Parrinello-Rahman dynamics. They are not unconstrained full NpT because laboratory `xx` strain is externally driven by active CAIMD.

## High-Strain Limitation

The target global engineering compression of about 0.30 is far outside the nominal small-strain regime of CAIMD V1. This local run is a software stress test, thermostat-under-loading test, and continuation-workflow test. It must not be used as a quantitative high-strain stress curve.
