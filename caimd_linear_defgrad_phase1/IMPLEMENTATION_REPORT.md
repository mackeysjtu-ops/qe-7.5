# CAIMD Linear Deformation-Gradient Phase 1 Implementation Report

Status: `CAIMD_LINEAR_DEFGRAD_PHASE1_IMPLEMENTED`

## Repository Identity

- Repository root: `/Users/mengqicheng/Downloads/Learning/qe-7.5`
- Active branch: `main`
- HEAD before change: `af7ebfc757c32e0b4f050ed9f2208e01020b2ee3`
- Source identity: local QE 7.5 CAIMD tree corresponding to the
  `qe-7.5-caimd-thermostat-af7ebfc` executable lineage.
- Initial Git state: tracked tree clean, with pre-existing untracked
  `examples/` result/workflow directories.
- No commit, tag, push, Slurm submission, scientific AIMD validation, or
  production executable overwrite was performed.

## Existing CAIMD Conventions Recorded

- Existing source map inspected: `Modules/input_parameters.f90`,
  `Modules/read_namelists.f90`, `PW/src/input.f90`,
  `PW/src/caimd_module.f90`, `PW/src/vcsubs.f90`, `PW/src/vcsmd.f90`,
  `PW/src/run_pwscf.f90`, and `test-suite/caimd/caimd_validation.py`.
- Existing input variables were `caimd`, `caimd_stride`,
  `caimd_step_strain(6)`, and `caimd_mask(6)` in `&CELL`.
- Existing event convention is unchanged: events occur at ionic steps
  `1, 1+stride, 1+2*stride, ...`.
- Existing six-component legacy order is unchanged:
  `(xx, yy, zz, yz, xz, xy)`.
- QE stores the cell matrix with lattice vectors as columns.  A physical
  affine relation `r_new = A * r_old` is represented as
  `H_new = A * H_old`.
- Existing CAIMD preserves fractional positions and transforms fractional
  velocity-like histories so peculiar Cartesian ionic velocities are preserved
  across imposed cell events.

## Files Changed

- `Modules/input_parameters.f90`
- `Modules/read_namelists.f90`
- `PW/src/input.f90`
- `PW/src/caimd_module.f90`
- `PW/Doc/INPUT_PW.def`
- `CAIMD_IMPLEMENTATION.md`
- `test-suite/caimd/caimd_linear_defgrad_math.py`
- `examples/caimd_linear_defgrad_phase1/README.md`
- `examples/caimd_linear_defgrad_phase1/local_z_engineering_compression.in.fragment`
- `examples/caimd_linear_defgrad_phase1/local_x_on_local_z_simple_shear.in.fragment`

Evidence files are under `caimd_linear_defgrad_phase1/evidence/`.

## Input API Added

```fortran
caimd_loading_mode = 'legacy_strain'   ! default
caimd_loading_mode = 'linear_defgrad'
caimd_step_defgrad(3,3) = 0.0d0        ! default
```

`caimd_step_defgrad` supports indexed namelist syntax such as
`caimd_step_defgrad(1,3)=2.0d-6` and `caimd_step_defgrad(3,3)=-2.0d-6`.

`legacy_strain` remains the default.  Existing inputs that omit
`caimd_loading_mode` and use `caimd_step_strain/caimd_mask` keep the previous
behavior.

## Mathematical Definition Implemented

Let `m` be the one-based CAIMD event index, `s=caimd_stride`,
`D=caimd_step_defgrad`, and `I` be the identity matrix.

```text
F_load(m) = I + (m * s) D
F_load(0) = I
A_event(m) = F_load(m) * inverse(F_load(m-1))
```

At a `linear_defgrad` event, the current PR candidate cell is updated by:

```text
H_projected = A_event(m) * H_PR_candidate
```

The implementation intentionally does not use `F_load(m) * H_reference` as a
replacement cell, so ordinary PR cell evolution accumulated between events is
preserved.

The Green-Lagrange diagnostic helper computes:

```text
E = 0.5 * (transpose(F_load) * F_load - I)
```

This diagnostic is not used as the driver in `linear_defgrad` mode.

## Position, Velocity, And PR-Cell Behavior

- Fractional atomic positions are preserved when the cell is affinely mapped.
- Fractional velocity, acceleration, and previous-acceleration histories are
  transformed through the existing CAIMD path to preserve peculiar Cartesian
  ionic velocities.
- In `linear_defgrad` mode, cell velocity and acceleration histories are
  left-multiplied by the same `A_event` as the cell at the event.
- The legacy symmetric tangent projection is skipped in `linear_defgrad` mode.
  Therefore zero entries of `caimd_step_defgrad` do not freeze the corresponding
  PR cell degrees of freedom.
- `caimd_mask(6)` remains a legacy symmetric-strain control.  In
  `linear_defgrad` mode it does not define loading and does not freeze zero
  deformation-gradient entries.

## Runtime Validation

The implementation rejects:

- unknown `caimd_loading_mode`;
- non-positive `caimd_stride`;
- non-finite `caimd_step_strain`;
- non-finite `caimd_step_defgrad`;
- non-finite or oversized legacy event strain increments;
- nonzero `caimd_step_defgrad` in `legacy_strain` mode;
- nonzero `caimd_step_strain` in `linear_defgrad` mode;
- unsupported CAIMD calculation scope;
- restart mode other than `from_scratch`;
- non-`all` `cell_dofree`;
- fixed-volume or fixed-area modes;
- fixed atomic coordinates for active CAIMD controls;
- unsupported ionic or cell thermostats for active CAIMD controls;
- singular, nearly singular, non-finite, or non-positive `F_load` and
  `A_event` determinants in `linear_defgrad` mode.

The determinant floor remains the existing CAIMD double-precision threshold:
`1.0d-18`.

## Output And Provenance

Initialization prints the selected loading mode, stride,
`caimd_step_defgrad` matrix, event convention, and restart support status.

The legacy `prefix.caimd` data rows remain 26 columns for parser compatibility.
Header lines are additive and include the selected loading mode.

For `linear_defgrad`, event matrices are written to:

```text
tmp_dir/prefix.caimd_defgrad
```

Rows record ionic step, event index, `F_load` before and after the event,
`A_event`, determinants, Green-Lagrange diagnostic components, and velocity
residual diagnostics.

## Deterministic Test Results

Command log paths:

- `caimd_linear_defgrad_phase1/evidence/test_py_compile.log`
- `caimd_linear_defgrad_phase1/evidence/test_linear_defgrad_math.log`
- `caimd_linear_defgrad_phase1/evidence/test_local_smoke.log`
- `caimd_linear_defgrad_phase1/evidence/test_caimd_validation_list_cases.log`
- `caimd_linear_defgrad_phase1/evidence/test_git_diff_check.log`

Results:

- Python syntax checks passed.
- `test-suite/caimd/caimd_linear_defgrad_math.py` passed.
- Local short `pw.x` smoke tests for D13 parsing, D33 parsing, zero-D identity
  events, short compression, short shear, and invalid determinant rejection
  passed.
- Existing CAIMD validation harness case listing passed.
- `git diff --check` passed.

Key deterministic values verified:

- `F_load(0)=I`.
- Explicit zero-D `linear_defgrad` inputs are active CAIMD controls and write
  identity `F_load`/`A_event` diagnostics.
- Simple shear with `stride=10`, `D(1,3)=2.0e-6`, `m=10000` gives
  `F_load(1,3)=0.20` and `det(F_load)=1.0`.
- Uniaxial compression with `stride=10`, `D(3,3)=-2.0e-6`, `m=7500` gives
  `F_load(3,3)=0.85` and `det(F_load)=0.85`.
- Simple-shear Green-Lagrange diagnostic gives `E_xz=E_zx=0.10` and
  `E_zz=0.02`.
- Compression Green-Lagrange diagnostic gives `E_zz=-0.13875`.
- Incremental compression map satisfies
  `A_event(3,3)=[1-m*stride*dc]/[1-(m-1)*stride*dc]`.
- Repeated simple-shear incremental maps recover the cumulative linear shear
  within floating-point tolerance.
- Positive and negative shear and normal loading are supported when the
  determinant remains positive.
- Singular/inverted compression, simultaneous legacy and defgrad loading
  definitions, and invalid mode are rejected.
- Legacy default mode remains `legacy_strain`.
- Applying `A_event` to a perturbed current cell is not equivalent to resetting
  to `F_load * H_reference`.
- Stride-10 event indices remain `1, 11, 21, ...`.

## Build Results

Command log paths:

- `caimd_linear_defgrad_phase1/evidence/build_configure.log`
- `caimd_linear_defgrad_phase1/evidence/build_qe_pw_exe.log`

Build tree:

```text
/private/tmp/qe_linear_defgrad_build
```

Commands:

```text
cmake -S . -B /private/tmp/qe_linear_defgrad_build -DQE_ENABLE_MPI=OFF -DQE_ENABLE_OPENMP=OFF -DQE_ENABLE_TEST=OFF
cmake --build /private/tmp/qe_linear_defgrad_build --target qe_pw_exe -j 4
```

Result: local serial `qe_pw_exe` build succeeded and produced:

```text
/private/tmp/qe_linear_defgrad_build/bin/pw.x
```

The build emitted existing QE/macOS compiler/linker warnings, but no
linear-defgrad CAIMD compile or link failure.

## Evidence Artifacts

- Complete material patch:
  `caimd_linear_defgrad_phase1/evidence/linear_defgrad_phase1.patch`
- Pre-change hashes:
  `caimd_linear_defgrad_phase1/evidence/pre_hashes.sha256`
- Post-change hashes:
  `caimd_linear_defgrad_phase1/evidence/post_hashes.sha256`
- Test logs:
  `caimd_linear_defgrad_phase1/evidence/test_*.log`
- Local smoke summary:
  `caimd_linear_defgrad_phase1/evidence/local_smoke_summary.json`
- Build logs:
  `caimd_linear_defgrad_phase1/evidence/build_*.log`

## Backward Compatibility

`legacy_strain` is the default.  The legacy small symmetric-strain path, mask
semantics, six-component order, event timing, tangent projection, `.caimd`
26-column data rows, and old `caimd_step_strain` validation are preserved for
legacy inputs.

The new `linear_defgrad` path is active when `caimd=.true.` and
`caimd_loading_mode='linear_defgrad'`, including the zero-matrix
`caimd_step_defgrad` case.  A zero-D input therefore behaves as an explicit
identity-loading CAIMD control and writes identity deformation-gradient
diagnostics.

## Known Limitations

- Restart remains unsupported.
- No scientific AIMD validation was run in this phase.
- No Miller-index parsing or structure rotation was implemented.
- No schedule-file parsing or arbitrary nonlinear loading paths were added.
- No thermostat, stress, force, timestep, or PR equation redesign was done.
- The generated PW input `.txt/.xml/.html` documentation was not regenerated
  from `PW/Doc/INPUT_PW.def` in this phase.

## Explicit Phase 2 Requirements

The following are not completed in this task:

- 64-atom 3000 K short AIMD validation;
- `[111]`-oriented structure generation and verification;
- `(111)[11-2]` shear-direction verification using atomic coordinates;
- partial-cell-response scientific validation;
- comparison with the reference stress-strain curve;
- long-run restart and segmentation;
- 75000-step compression;
- 100000-step shear;
- production readiness.
