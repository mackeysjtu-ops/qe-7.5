# Project Guidance

This repository is an initialized Quantum ESPRESSO 7.5 source tree. The untouched
source baseline is preserved by commit `a11913ac6cd37f39cf31056e8e4823ec7a5de49c`
with subject `Original QE 7.5 source baseline`.

For constrained AIMD work, keep QE source changes narrow and follow the
implemented design notes in `CAIMD_IMPLEMENTATION.md`. The first implementation
report is `CAIMD_V1_IMPLEMENTATION_REPORT.md`.

## Validated Source Map

- `Modules/input_parameters.f90`: declares PW namelist variables for
  `&CONTROL`, `&IONS`, and `&CELL`.
- `Modules/read_namelists.f90`: sets defaults, validates allowed namelist
  values, broadcasts input variables, and applies calculation-dependent defaults.
- `PW/src/input.f90`: maps `calculation='vc-md'` and `cell_dynamics='pr'` to
  the PR variable-cell MD switch `calc='cd'`, forces stress evaluation for
  moving-cell runs, converts input pressure to internal units, initializes
  `cell_dofree`, and sets the default cell mass.
- `PW/src/pwscf.f90` and `PW/src/run_pwscf.f90`: main `pw.x` execution path,
  electronic loop, force and stress calls, XML step output, and call to
  `move_ions`.
- `PW/src/move_ions.f90`: dispatches MD algorithms; variable-cell shape MD is
  routed to `vcsmd`.
- `PW/src/vcsmd.f90`: bridge between PW globals and the variable-cell MD
  routines. It converts Cartesian/fractional quantities, scales stress by
  volume, calls `vcinit` on the first step and `vcmove` later, writes the
  updated cell, volume, reciprocal lattice, and positions back to PW globals,
  and reads/writes the MD history file.
- `PW/src/vcsubs.f90`: PR/Wentzcovitch Beeman integrator. This is where stress,
  kinetic virial, external pressure, cell mass, `iforceh`, cell velocities, and
  cell accelerations enter the cell equation of motion.
- `PW/src/forces.f90`, `PW/src/stress.f90`, and `PW/src/pwcom.f90`: compute and
  store atomic forces and stress in `force_mod`.
- `Modules/cell_base.f90` and `Modules/deviatoric.f90`: existing cell variables,
  `iforceh`, `cell_dofree`, pressure, fixed-volume/fixed-area handling, and
  direct matrix-element cell constraints.

## CAIMD Boundaries

Implemented `&CELL` input names are:

- `caimd`: logical user switch.
- `caimd_stride`: integer stride between imposed deformation events.
- `caimd_step_strain(6)`: nominal per-MD-step strain increment.
- `caimd_mask(6)`: logical mask for controlled strain components.

The six-component order is `(xx, yy, zz, yz, xz, xy)`.

The first implementation applies only to `calculation='vc-md'` with
`cell_dynamics='pr'`, i.e. `calc='cd'`. Atomic forces, SCF/electronic
minimization, stress calculation, thermostats, and ordinary atomic motion remain
behaviorally unchanged.

Input parsing stays in the existing input modules. Reusable CAIMD strain,
rotation, output, and history projection logic lives in
`PW/src/caimd_module.f90`. The validated insertion boundary is the PR Beeman
path in `PW/src/vcsubs.f90`, with `PW/src/vcsmd.f90` and `PW/src/input.f90`
used for state wiring and restart/history checks. Projection of only the
returned `avec` is insufficient because `vcmove` uses cell velocity history in
the metric coupling that updates atomic acceleration, velocity, and position.

At imposed events, preserve fractional positions and preserve peculiar
Cartesian atomic velocities by transforming `ratd` with
`inverse(H_projected) * H_candidate`; do not zero atomic thermal motion.
Transform atomic acceleration histories consistently. When CAIMD is active with
any true mask, suppress the antisymmetric cell increment and project `avecd`,
`avec2d`, and `avec2di`; when all mask values are false, CAIMD is an exact
no-op.

The first implementation rejects restart, stale MD history, non-`all`
`cell_dofree`, fixed-volume/fixed-area modes, non-PR cell dynamics, and
non-`vc-md` calculations.

## Build Notes

Implementation files are `Modules/input_parameters.f90`,
`Modules/read_namelists.f90`, `PW/src/input.f90`,
`PW/src/caimd_module.f90`, `PW/src/run_pwscf.f90`, `PW/src/vcsubs.f90`,
`PW/src/vcsmd.f90`, `PW/src/Makefile`, and `PW/CMakeLists.txt`.

The CMake serial `qe_pw_exe` build passed out of tree. A traditional
`./configure` in a temporary source copy succeeded, but `make pw` failed before
reaching `PW/src` in a pre-existing macOS `cpp`/`LAXlib` preprocessing path
(`cc: error: no input files`).

Do not reorganize the QE directory structure for this feature.

## CAIMD Warnings

- `cell_dofree`/`iforceh` constrains direct 3x3 cell-matrix elements, not the
  six physical symmetric strain components requested for CAIMD. The first
  implementation must reject non-`all` `cell_dofree`.
- The shear convention must be explicit. The current plan treats off-diagonal
  `caimd_step_strain` values as tensor shear components; engineering shear would
  differ by a factor of two.
- Stress/strain output must pair stress with the cell on which that stress was
  calculated. In the current call path, XML step output happens before
  `move_ions`, while `.cel` trajectory output happens after the cell move.
- Resetting masked PR cell velocities or accelerations is a driven-dynamics
  choice, not ordinary energy-conserving NpT dynamics. Document it in any future
  input help and tests.
- QE's baseline PR path does not enforce a triangular or explicitly symmetric
  cell representation. CAIMD must remove internal antisymmetric increments when
  active, while preserving the all-false-mask regression case.
