# Project Guidance

This repository is a Git-initialized Quantum ESPRESSO 7.5 source tree. The
untouched source baseline is commit
`a11913ac6cd37f39cf31056e8e4823ec7a5de49c`
(`Original QE 7.5 source baseline`).

## Canonical CAIMD Documents

- `CAIMD_IMPLEMENTATION.md`: stable CAIMD V1 algorithm, source map, input
  semantics, runtime boundaries, and validation expectations.
- `LEONARDO_BUILD.md`: repeatable Leonardo GPU/MPI build, deployment,
  incremental rebuild, validation harness, evidence, and recovery manual.
- `CAIMD_V1_IMPLEMENTATION_REPORT.md`: historical local implementation report.
- `CAIMD_LEONARDO_GPU_VALIDATION_REPORT.md`: historical Leonardo validation
  evidence and later hardening checkpoints.

Read `CAIMD_IMPLEMENTATION.md` before changing CAIMD source behavior. Read
`LEONARDO_BUILD.md` before any Leonardo transfer, build, Slurm allocation, or
remote validation.

## Source Boundaries

CAIMD V1 is a PW-only constrained AIMD mode for:

- `calculation='vc-md'`;
- `cell_dynamics='pr'`;
- `cell_dofree='all'`;
- `restart_mode='from_scratch'`;
- no fixed-volume/fixed-area cell modes;
- no fixed atomic coordinates when any CAIMD mask component is active;
- no ionic or cell thermostat/rescaling when any CAIMD mask component is active.

The implementation must remain an exact no-op when `caimd=.false.` or when all
`caimd_mask` values are false. CP does not implement CAIMD; `caimd=.true.` in
CP input must fail explicitly.

## Stable Source Map

- `Modules/input_parameters.f90`: declares PW/CP namelist variables, including
  CAIMD `&CELL` names.
- `Modules/read_namelists.f90`: defaults, validation, broadcast, ignored
  `&CELL` handling, and CP CAIMD rejection.
- `PW/src/input.f90`: maps PW `vc-md`/PR into `calc='cd'`, initializes CAIMD
  runtime state, applies unsupported-boundary rejections, pressure conversion,
  and cell mass setup.
- `PW/src/run_pwscf.f90`: electronic step, force/stress calls, XML step output,
  CAIMD premove logging, and `move_ions`.
- `PW/src/move_ions.f90`: dispatches variable-cell MD to `vcsmd`.
- `PW/src/vcsmd.f90`: bridges PW globals to variable-cell MD routines, rejects
  unsupported CAIMD history, updates cell/positions, and records temperature
  diagnostics.
- `PW/src/vcsubs.f90`: Parrinello-Rahman/Wentzcovitch Beeman integrator and
  CAIMD projection insertion point.
- `PW/src/caimd_module.f90`: CAIMD state, strain projection, history
  transformation, diagnostics, checks, and `.caimd` output.
- `PW/src/forces.f90`, `PW/src/stress.f90`, `PW/src/pwcom.f90`: ordinary QE
  atomic force and stress computation/storage; CAIMD V1 must not alter these.

## Development Rules

- Keep CAIMD changes narrow and preserve ordinary QE behavior outside the
  documented boundary.
- Do not redesign CAIMD physics while hardening parser, validation, or
  documentation behavior.
- Preserve the six-component strain order `(xx, yy, zz, yz, xz, xy)` and the
  tensor shear convention.
- Preserve fractional positions and peculiar Cartesian atomic velocities across
  imposed cell events.
- Keep restart unsupported until CAIMD protocol state is persisted and tested.
- Do not reorganize the QE directory structure for this feature.
- Use the persistent harness in `test-suite/caimd/caimd_validation.py` for
  focused CAIMD validation. It must not download pseudopotentials.

## Git And Validation

Treat pre-existing uncommitted changes as user-owned work. Do not reset, clean,
or overwrite unrelated files. Validate locally before Leonardo whenever
practical, and keep local source changes tracked before deploying them remotely.

Do not push to any remote unless explicitly requested.
