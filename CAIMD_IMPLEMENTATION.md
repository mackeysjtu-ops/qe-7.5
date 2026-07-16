# CAIMD V1 Implementation Notes

This document is the stable software and algorithm specification for the
Quantum ESPRESSO 7.5 constrained ab initio molecular dynamics (CAIMD) V1
implementation in this repository.

CAIMD V1 is a driven small-strain variable-cell AIMD protocol implemented only
for PW `vc-md` with Parrinello-Rahman cell dynamics. It constrains selected
nominal strain components to a prescribed schedule while leaving ordinary QE
electronic minimization, atomic force computation, stress computation, and
unmasked PR cell response unchanged.

## User Interface

The CAIMD input names live in the PW/CP shared `&CELL` namelist:

```fortran
caimd                  ! logical, default .false.
caimd_stride           ! integer, default 1
caimd_step_strain(6)   ! real, default 0.0
caimd_mask(6)          ! logical, default .false.
```

The six-component order is:

```text
(xx, yy, zz, yz, xz, xy)
```

Off-diagonal values are tensor shear components. Engineering shear is twice the
reported tensor shear.

An imposed event increment is:

```text
event_increment(i) = caimd_stride * caimd_step_strain(i)
```

Each event component must satisfy:

```text
abs(event_increment(i)) <= 1.0d-2
```

This is a per-event hard limit. CAIMD V1 also emits a one-time warning on the
I/O rank when the maximum absolute cumulative nominal masked strain first
exceeds `1.0d-2`. The cumulative warning is diagnostic only; it does not alter
the trajectory.

## Supported Boundary

CAIMD V1 is supported only when all of the following hold:

- PW executable, not CP;
- `calculation='vc-md'`;
- `cell_dynamics='pr'`;
- `cell_dofree='all'`;
- `restart_mode='from_scratch'`;
- no fixed-volume or fixed-area cell modes;
- no existing/stale MD history file;
- no fixed atomic coordinate component when `caimd=.true.` and any
  `caimd_mask` component is true;
- `ion_temperature='not_controlled'` when `caimd=.true.` and any mask component
  is true;
- `cell_temperature='not_controlled'` when `caimd=.true.` and any mask
  component is true.

When `caimd=.false.` or all `caimd_mask` values are false, CAIMD must be an
exact no-op. In particular, inactive/all-false CAIMD must not reject ordinary
fixed atomic coordinates or thermostats beyond baseline QE behavior.

CP does not implement CAIMD. Because CP can parse the shared `&CELL` variables,
`caimd=.true.` in CP input must fail explicitly instead of being silently
ignored.

## Ignored `&CELL` Policy

Non-variable-cell PW runs normally ignore `&CELL`. CAIMD is an exception only to
prevent silent loss of a user request:

- no `&CELL`: ordinary non-variable-cell input is unchanged;
- an immediate valid ignored `&CELL` with `caimd=.false.`: accepted and ignored;
- an immediate ignored `&CELL` with `caimd=.true.`: explicit error requiring
  PW `calculation='vc-md'`;
- an immediate malformed ignored `&CELL`: normal namelist error;
- a `&CELL` placed after card data: not scanned by the ignored-`&CELL` probe,
  so card parsing owns the outcome.

The probe must inspect only the next significant record after the prior
namelist. It must not search through `ATOMIC_SPECIES`, `CELL_PARAMETERS`,
`ATOMIC_POSITIONS`, `K_POINTS`, or later card data.

## Physical Convention

QE stores the cell matrix as columns of lattice vectors. CAIMD uses the total
deformation gradient

```text
F = H * inverse(H_ref)
```

and the small symmetric strain

```text
E = sym(F - I)
```

where `H_ref` is the CAIMD reference cell recorded at the start of the run.
At each imposed event, masked components of the target nominal strain are
advanced by `caimd_stride * caimd_step_strain`. Unmasked components come from
the ordinary PR candidate cell.

The projected cell is reconstructed as:

```text
H_projected = (I + E_projected) * H_ref
```

This is a first-version small-strain protocol, not a finite-strain constitutive
framework.

## PR Coupling And History Projection

The validated insertion point is inside the PR Beeman path in
`PW/src/vcsubs.f90`, with state bridging in `PW/src/vcsmd.f90`.

Projection of only the returned cell matrix is insufficient. The PR path also
uses cell velocity and acceleration histories in the metric coupling that
updates atomic acceleration, velocity, and position. Therefore CAIMD V1:

- projects tangent cell histories for masked symmetric strain components;
- suppresses antisymmetric internal cell increments while active;
- projects `avecd`, `avec2d`, and `avec2di`;
- transforms atomic velocity and acceleration histories consistently.

At a projection event, CAIMD preserves fractional atomic positions and
preserves peculiar Cartesian atomic velocities. The fractional velocity-like
history arrays are transformed as:

```text
ratd_projected  = inverse(H_projected) * H_candidate * ratd_candidate
rat2d_projected = inverse(H_projected) * H_candidate * rat2d_candidate
rat2di_projected = inverse(H_projected) * H_candidate * rat2di_candidate
```

The same transformation principle applies to the stored acceleration histories.
Atomic thermal motion is not zeroed.

## Stress And Atomic Forces

CAIMD V1 does not alter force or stress calculation.

The ordinary PW path is:

- `PW/src/pwscf.f90`: input and top-level run dispatch;
- `PW/src/run_pwscf.f90`: electronic step, force and stress calls, step output,
  and `move_ions`;
- `PW/src/forces.f90`: atomic force driver and fixed-coordinate masking;
- `PW/src/stress.f90`: stress driver and stress tensor sum;
- `PW/src/pwcom.f90`: `force_mod` storage for force and stress;
- `PW/src/move_ions.f90`: MD dispatch;
- `PW/src/vcsmd.f90` and `PW/src/vcsubs.f90`: variable-cell MD state transfer
  and PR integration.

The PR cell equation still sees the ordinary stress tensor, kinetic virial,
external pressure, cell mass, and cell constraints. CAIMD constrains selected
resulting cell strain components after the ordinary PR candidate update rather
than modifying the stress tensor or atomic forces.

## Pressure

Input `press` is read in GPa from `&CELL` and converted in `PW/src/input.f90`
to QE internal pressure units:

```fortran
press_ = press / ry_kbar
```

The pressure term enters the variable-cell PR equations through the existing
cell dynamics machinery in `PW/src/vcsubs.f90`. CAIMD V1 leaves this pressure
path unchanged.

## Diagnostics

CAIMD writes a premove diagnostic file:

```text
tmp_dir/prefix.caimd
```

The file has 26 columns:

```text
ionic_step
time_ps
imposed_events_done
following_move_is_event
nominal_xx nominal_yy nominal_zz nominal_yz nominal_xz nominal_xy
actual_xx actual_yy actual_zz actual_yz actual_xz actual_xy
stress_xx_GPa stress_yy_GPa stress_zz_GPa stress_yz_GPa stress_xz_GPa stress_xy_GPa
temperature_K
volume_bohr3
velocity_jump_max
velocity_jump_rms
```

Important semantics:

- `stress_*_GPa` is the stress computed for the pre-move cell of that ionic
  step.
- `actual_*` is the small strain of that same pre-move cell relative to the
  CAIMD reference cell.
- `following_move_is_event` tells whether the move after the logged stress is
  scheduled to apply an imposed strain event.
- `temperature_K` is a CAIMD diagnostic derived from the ionic kinetic-energy
  quantity passed through the existing `vcsmd` path when available; otherwise it
  falls back to the current input/target temperature value available to the
  premove writer. It is not a new thermostat variable and does not include a
  redesigned CAIMD temperature definition.
- `velocity_jump_max` and `velocity_jump_rms` are residual diagnostics of the
  enforced velocity-history transformation. At a correctly transformed event
  they should be at roundoff scale; they are not physical velocity jumps to be
  interpreted as a thermostat or impulse.

Because QE XML step output occurs before `move_ions` and `.cel` trajectory
output occurs after the move, stress-strain analysis should pair `.caimd`
stress/strain rows with the premove cell, not blindly with the following `.cel`
row.

## Restart Policy

Restart is intentionally unsupported in V1. A valid restart would require
persisting and validating extra CAIMD protocol state:

- event count and stride phase;
- cumulative imposed strain;
- reference cell;
- history projection convention;
- output/log append state.

The implementation rejects restart mode, stale MD history before cleanup, and
unexpected MD history at the start of `vcsmd`.

## Numerical Checks

CAIMD V1 requires:

- positive `caimd_stride`;
- finite `caimd_step_strain` values;
- finite event increments;
- per-event increment bound `abs(stride * step_strain) <= 1.0d-2`;
- positive projected-cell determinant above the CAIMD determinant floor;
- projected-cell condition number below the CAIMD condition bound;
- explicit error paths for singular or ill-conditioned 3x3 inversions.

## Implementation Files

Changed source files for CAIMD V1 are:

- `Modules/input_parameters.f90`
- `Modules/read_namelists.f90`
- `PW/src/input.f90`
- `PW/src/caimd_module.f90`
- `PW/src/run_pwscf.f90`
- `PW/src/vcsmd.f90`
- `PW/src/vcsubs.f90`
- `PW/src/Makefile`
- `PW/CMakeLists.txt`

Validation support is in:

- `test-suite/caimd/caimd_validation.py`

## Validation Expectations

Focused validation should include:

- ordinary PW no-CAIMD behavior;
- `caimd=.false.` and all-false mask no-op regressions;
- active zero-increment constraints;
- stride scheduling;
- tensor shear convention;
- non-orthogonal-cell left-multiplication behavior;
- antisymmetric/rotation suppression diagnostics;
- `.caimd` 26-column schema and finite diagnostic columns;
- explicit rejections for unsupported calculation, non-PR cell dynamics,
  non-`all` `cell_dofree`, restart/stale history, invalid stride, oversized
  event increment, fixed atomic coordinates, thermostats/rescaling, and CP
  `caimd=.true.`;
- parser behavior for ignored non-variable-cell `&CELL` without card scanning.

Leonardo-specific build and validation instructions are intentionally kept in
`LEONARDO_BUILD.md`, not in this algorithm document.
