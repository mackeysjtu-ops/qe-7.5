# CAIMD Thermostat Integration Audit

Status: implementation audit complete; source changes follow this document.

Repository state at audit start:

- Repository root:
  `/Users/mengqicheng/Downloads/Learning/qe-7.5`
- Branch: `main`
- HEAD: `13ba2ebbb43ff296305ae42183699fda995079da`
  (`Validate 300 K vc-md initialization`)
- Initial `git status --short --branch`: `## main`

## Source Path

The active CAIMD boundary remains PW `calculation='vc-md'` with
`cell_dynamics='pr'`. In `PW/src/input.f90`, that input is mapped to
`calc='cd'`, `lmd=.true.`, and `lmovecell=.true.`. The ionic move is dispatched
as:

```text
PW/src/run_pwscf.f90
  force/stress/XML step output/CAIMD premove output
  -> PW/src/move_ions.f90
     -> PW/src/vcsmd.f90
        -> PW/src/vcsubs.f90:vcinit, first step
        -> PW/src/vcsubs.f90:vcmove, later steps
```

`PW/src/move_ions.f90` calls `dynamics_module: langevin_md` only for
`calc(1:1)=='l'`, the fixed-cell Langevin branch. PW `vc-md`/PR does not use
that branch. It always reaches `vcsmd` for variable-cell shape dynamics.

## Thermostat Modes In PW `vc-md`

The input parser accepts many `ion_temperature` strings, but they are not all
implemented in the PW `vc-md`/PR integrator.

- `not_controlled`: supported existing CAIMD behavior. No ionic thermostat is
  applied.
- `nose`: implemented in the actual `vc-md` call graph. `input.f90` sets
  `tnosep=.true.`, `pos_iosys` initializes the ionic Nose-Hoover variables, and
  `vcsmd` calls `ions_nosevel` before `vcmove` and `ions_noseupd` after
  `vcmove`. Inside `vcmove`, the Nose friction term is applied to ionic
  accelerations in lattice coordinates as a function of the current
  lattice-coordinate velocity `ratd`.
- `rescaling`: partly implemented in `vcmove`, but not as an ionic-only
  thermostat. The rescaling test uses the averaged total kinetic quantity
  `avk`, and for dynamics cases with `calc(2:2)=='d'` it rescales both
  `ratd` and `avecd`. This couples to PR cell velocity and violates the CAIMD
  requirement that ionic temperature control must not thermostat imposed or
  unmasked cell degrees of freedom.
- `berendsen`, `svr`, and `andersen`: implemented for the fixed-cell
  `dynamics_module` path through `apply_thermostat`, but PW `vc-md` does not
  call that routine. In `vc-md`, these keywords set parser-level variables and
  the starting temperature, but do not apply the corresponding thermostat during
  `vcsmd/vcmove`.
- `initial`: initializes velocities when no explicit velocity card is supplied;
  it is not a sustained thermostat.
- `rescale-v`, `rescale-T`, and `reduce-T`: handled by
  `dynamics_module.apply_thermostat` for fixed-cell MD, but not by the
  `vcsmd/vcmove` variable-cell path.
- Langevin: unavailable for PW `vc-md`/PR in this source. `vc-md` ignores
  non-`beeman` `ion_dynamics`, and the Langevin implementation is reached only
  by fixed-cell `calc='ld'` or `calc='ls'`.

## Selected Allowlist

The active-CAIMD allowlist should be deliberately narrow:

- `not_controlled`
- `nose`

`nose` is the minimum scientifically valid thermostat already present in the
real PW `vc-md`/PR call graph. It is deterministic, canonical extended-system
temperature control and acts on ionic velocities only through the existing
Nose-Hoover friction term. It does not require inventing new input keywords or
duplicating fixed-cell thermostat logic.

All other `ion_temperature` modes must continue to reject early for active
CAIMD with a precise message.

`cell_temperature` remains restricted to `not_controlled`. This task is about
ionic temperature control only.

## Velocity And Basis Semantics

`vcsmd` converts PW Cartesian atomic positions and forces to lattice
coordinates before calling `vcinit` or `vcmove`. In the variable-cell routines:

- `avec` is the physical cell matrix in bohr.
- `rat` is fractional/lattice atomic position.
- `ratd` is the lattice-coordinate atomic velocity-like history.
- The physical Cartesian ionic velocity is `avec * ratd`.
- `rat2d` and `rat2di` are lattice-coordinate acceleration histories.
- `avecd`, `avec2d`, and `avec2di` are PR cell velocity/acceleration histories.

The existing CAIMD event projection changes the PR candidate cell to the
projected constrained cell and transforms `ratd`, `rat2d`, and `rat2di` by
`inverse(H_projected) * H_candidate`. That preserves the physical Cartesian
ionic velocity across the basis change. Nose friction is applied before this
projection, and the transformed histories must be the ones committed to the
MD history file.

For active CAIMD with Nose, the kinetic quantities used for the Nose update
must be refreshed after the final projected cell and transformed `ratd` are
known. Otherwise the thermostat update may use kinetic energy from a different
cell/velocity representation than the committed history.

## Temperature Convention

The scientifically relevant thermostat temperature is ionic kinetic
temperature, not PR cell kinetic temperature. For active CAIMD diagnostics, the
26-column `.caimd` schema should remain unchanged, but the `temperature_K`
field should be documented as the ionic kinetic temperature recorded from the
last completed move when available.

The first premove row is written before `vcinit` has built the `vcsmd` velocity
history, so it may fall back to the input target temperature. Subsequent rows
lag the thermostat update by one ionic move because `.caimd` output is written
before `move_ions`.

## Continuation State

Full CAIMD restart remains unsupported. Exact Nose-Hoover continuation would
also require persisting and restoring the Nose variables:

- `xnhp0`, `xnhpm`, `xnhpp`
- `vnhp`
- chain dimensions and atom-to-chain mapping
- accumulated random-number state is not relevant to deterministic Nose, but
  would be relevant for stochastic thermostats if ever supported.

For segmented validation without formal restart, the terminal cell, atomic
positions, and physical Cartesian velocities can be extracted into a new
`from_scratch` segment, but that is not an exact Nose extended-state restart.

## Implementation Decision

Implement active CAIMD compatibility for `ion_temperature='nose'` only. The
source patch should:

- replace the blanket active-CAIMD thermostat rejection with the allowlist above;
- keep cell thermostats rejected;
- refresh projected ionic kinetic energy and Nose chain kinetic inputs after
  CAIMD final projection;
- record CAIMD diagnostic temperature from ionic kinetic energy;
- keep disabled and all-false-mask CAIMD as no-op regressions;
- leave restart unsupported.
