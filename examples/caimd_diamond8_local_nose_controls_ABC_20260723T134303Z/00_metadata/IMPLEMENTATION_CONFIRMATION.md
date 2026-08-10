# Implementation Confirmation

Generated: 2026-07-23T13:53:34Z

This file records a read-only source inspection of the current QE tree. No
implementation file was modified.

## Ordinary `vc-md` Nose Call Path

- `PW/src/input.f90` maps `calculation='vc-md'` plus `cell_dynamics='pr'` to
  `calc='cd'` and enables MD, forces, and variable-cell motion
  (PW/src/input.f90:621-637).
- `PW/src/run_pwscf.f90` computes forces and stress, writes the CAIMD pre-move
  diagnostic if enabled, then calls `move_ions`
  (PW/src/run_pwscf.f90:232-257).
- `PW/src/move_ions.f90` dispatches variable-cell shape MD to `vcsmd`
  (PW/src/move_ions.f90:386-393).
- `PW/src/vcsmd.f90` calls `ions_nosevel` before `vcmove` and `ions_noseupd`
  after `vcmove` when `tnosep` is active
  (PW/src/vcsmd.f90:334-351).
- `PW/src/vcsubs.f90` applies the ionic Nose friction term to `rat2d` before
  the PR/Beeman update
  (PW/src/vcsubs.f90:696-712).

## CAIMD Allowlist And All-False No-Op

- Active CAIMD is defined as `caimd=.true.` with at least one mask component
  true. Only active CAIMD is restricted to `ion_temperature='not_controlled'`
  or `'nose'`, and it rejects cell thermostats
  (PW/src/input.f90:1739-1756).
- `caimd_controls_cell()` returns `caimd_enabled_ .AND. ANY(caimd_mask_)`
  (PW/src/caimd_module.f90:113-115). Therefore the
  all-false Control B is expected to write CAIMD diagnostics because
  `caimd_enabled_` is true, but it must not enter the projection/tangent path.

## Projection Location And State Representation

- Active final projection is called near the end of `vcmove` after the candidate
  PR cell update and before the final metric/energy refresh
  (PW/src/vcsubs.f90:1120-1150).
- `ratd` is the current fractional/lattice velocity history, `rat2d` is the
  current fractional/lattice acceleration, and `rat2di` is the previous
  acceleration history. During initial velocity input, Cartesian
  `ATOMIC_VELOCITIES` are temporarily held in `rat2di`, then converted to
  `ratd` through `sigma/vcell`
  (PW/src/vcsubs.f90:185-265).
- `caimd_transform_atomic_history` transforms `ratd`, `rat2d`, and `rat2di`
  from the PR candidate basis into the projected basis while preserving
  physical Cartesian velocities to roundoff
  (PW/src/caimd_module.f90:321-349).

## Energy And Nose Synchronization

- After active projection, `vc_refresh_projected_energies` recomputes ionic
  kinetic energy `eka`, cell kinetic energy `ekla`, pressure, and ionic
  temperature from the projected `ratd` and metric `g`
  (PW/src/vcsubs.f90:1159-1241).
- When Nose is active, `vc_refresh_projected_nose_ekin` refreshes `ekin2nhp`
  from the same projected `ratd`/`g` before `ions_noseupd`
  (PW/src/vcsubs.f90:1245-1278).

## `.caimd`, `.eal`, And DOF

- `.caimd` is written by `caimd_write_premove` from `run_pwscf` before
  `move_ions`, so each row describes the state before the following ionic move.
  Its temperature is `caimd_temperature_` if available, otherwise the input
  temperature fallback at the first row
  (PW/src/run_pwscf.f90:241-257;
  PW/src/caimd_module.f90:215-267).
- `caimd_note_ionic_temperature` defines `.caimd temperature_K` as
  `2*eka/(max(3*(nat-1),1)*k_B)` with `eka` from `vcsmd`
  (PW/src/caimd_module.f90:201-213).
- `.eal` is written in `vcsmd` as
  `uta, eka, eta, utl, ekla, etl, nst` with format `6D12.5,I6`
  (PW/src/vcsmd.f90:397-411;
  PW/src/vcsmd.f90:498-506).
- For the unconstrained 8-atom system the ionic DOF convention is
  `3*nat - 3 = 21`, matching `dynamics_module::get_ndof`
  (PW/src/dynamics_module.f90:2503-2527).

## Velocity History Caveat

The PR `.md` file is the authoritative terminal history for `ratd`; the `.vel`
trajectory blocks are also inspected in this workflow, but they are reported as
diagnostic only unless they agree with `.eal`/terminal `.md` kinetic energy.
