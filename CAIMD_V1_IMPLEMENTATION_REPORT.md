# CAIMD V1 Implementation Report

Status: SUCCESS

This report records the first working constrained AIMD implementation in this
Quantum ESPRESSO 7.5 source tree. The implementation commit that contains this
file is created after the report is staged; its hash is reported in the final
Codex response.

## Repository And Commit Baseline

- Repository root: `/Users/mengqicheng/Downloads/Learning/qe-7.5`
- Branch: `main`
- Source baseline commit: `a11913ac6cd37f39cf31056e8e4823ec7a5de49c`
  (`Original QE 7.5 source baseline`)
- Documentation commit: `f7256b18f21913f506144e256fc7e8986fbb46e4`
  (`Document CAIMD implementation design`)
- Initial tracked/staged diffs before work: none.
- Initial untracked files before work: `AGENTS.md`,
  `CAIMD_IMPLEMENTATION.md`

## Outcome

Implemented CAIMD V1 for PW `calculation='vc-md'` with
`cell_dynamics='pr'`, `cell_dofree='all'`, and `restart_mode='from_scratch'`.

Implemented input variables in `&CELL`:

- `caimd`: logical, default `.FALSE.`
- `caimd_stride`: positive integer, default `1`
- `caimd_step_strain(6)`: small tensor-strain increment per MD step
- `caimd_mask(6)`: logical mask for controlled strain components

Component order is:

```text
(xx, yy, zz, yz, xz, xy)
```

Off-diagonal components are tensor shear. Engineering shear is twice the tensor
component.

## Implemented Source Changes

- `Modules/input_parameters.f90`
  - Declares CAIMD input variables and adds them to `&CELL`.
- `Modules/read_namelists.f90`
  - Sets defaults, broadcasts CAIMD inputs, validates positive stride, finite
    strain values, and `max(abs(caimd_stride*caimd_step_strain)) <= 1.0d-2`.
- `PW/src/input.f90`
  - Configures CAIMD after QE resolves `calculation`, `cell_dynamics`,
    `cell_dofree`, pressure, fixed-volume/fixed-area state, and restart mode.
  - Rejects stale `tmp_dir/prefix.md` history before from-scratch cleanup can
    delete it.
- `PW/src/caimd_module.f90`
  - New PW-local module with CAIMD runtime state, six-component strain mapping,
    event schedule, determinant/condition checks, reference-cell strain
    projection, tangent-space history projection, velocity-preserving atomic
    history transform, and pre-move output.
- `PW/src/vcsubs.f90`
  - Projects PR cell velocity and acceleration histories inside `vcinit` and
    `vcmove`.
  - Applies final strain-and-hold projection before the committed cell state is
    returned and written to MD history.
- `PW/src/vcsmd.f90`
  - Rejects unsupported start-of-run MD history and records CAIMD temperature
    diagnostics.
- `PW/src/run_pwscf.f90`
  - Writes pre-move CAIMD stress-strain rows after `stress(sigma)` and before
    `move_ions`.
- `PW/src/Makefile`
  - Adds `caimd_module.o`.
- `PW/CMakeLists.txt`
  - Adds `src/caimd_module.f90`.
- `AGENTS.md` and `CAIMD_IMPLEMENTATION.md`
  - Updated to reflect the delivered implementation and validation status.

## Algorithmic Decisions

- Active CAIMD means `caimd=.TRUE.` and at least one `caimd_mask` component is
  true.
- If all mask values are false, CAIMD writes diagnostics but skips all
  projection; the PR trajectory remains an exact no-op regression.
- The first event is applied on the first ionic move satisfying
  `MOD(istep-1, caimd_stride) == 0`.
- Event increment is `caimd_stride * caimd_step_strain`, but only masked
  components are accumulated as imposed nominal strain.
- Masked components are controlled against the initial CAIMD reference cell:
  `Hnew = (I + Enew) * Href`. This avoids finite-product drift between nominal
  and actual small strain.
- Unmasked symmetric strain components are taken from the PR candidate cell.
- Active CAIMD suppresses the antisymmetric cell increment and projects
  `avecd`, `avec2d`, and `avec2di` in laboratory tangent space.
- Fractional atomic positions are preserved across the projected cell change.
- Peculiar Cartesian atomic velocity `u = H * ratd` is preserved by transforming
  `ratd`, `rat2d`, and `rat2di` with
  `inverse(H_projected) * H_candidate`.

## Output

CAIMD V1 writes a text stream at:

```text
tmp_dir/prefix.caimd
```

The stream is written on the I/O rank only, before the ionic/cell move, so
stress and strain are aligned with the cell on which the stress was computed.
No XML schema changes were made.

Columns:

```text
ionic_step time_ps imposed_events_done following_move_is_event
nominal_xx nominal_yy nominal_zz nominal_yz nominal_xz nominal_xy
actual_xx actual_yy actual_zz actual_yz actual_xz actual_xy
stress_xx_GPa stress_yy_GPa stress_zz_GPa stress_yz_GPa stress_xz_GPa stress_xy_GPa
temperature_K volume_bohr3 velocity_jump_max velocity_jump_rms
```

## Rejections And Limits

CAIMD V1 rejects:

- `restart_mode='restart'` or equivalent restart state;
- stale MD history `tmp_dir/prefix.md`;
- non-`all` `cell_dofree`;
- fixed-volume or fixed-area modes;
- non-PR cell dynamics;
- non-`vc-md` calculations;
- non-positive `caimd_stride`;
- non-finite `caimd_step_strain`;
- event components larger than `1.0d-2`;
- projected cells with non-positive/tiny determinant or excessive condition
  number.

Known V1 limitations:

- No restart support and no persisted CAIMD protocol state.
- No finite-strain measure or rotating-load-frame protocol.
- No XML/QEXSD schema integration.
- No user-facing generated input documentation update in `PW/Doc`.
- Validation is local serial/small-cell smoke and regression testing, not a
  production scientific validation campaign.

## Build Evidence

Modified CMake build:

```text
cmake -S . -B /private/tmp/qe_caimd_v1_build -DQE_ENABLE_MPI=OFF -DQE_ENABLE_OPENMP=OFF -DQE_ENABLE_TEST=OFF
cmake --build /private/tmp/qe_caimd_v1_build --target qe_pw_exe -j 4
```

Result: succeeded. The executable was
`/private/tmp/qe_caimd_v1_build/bin/pw.x`.

Baseline CMake build from commit
`a11913ac6cd37f39cf31056e8e4823ec7a5de49c`:

```text
git worktree add --detach /private/tmp/qe_caimd_v1_baseline_src a11913ac6cd37f39cf31056e8e4823ec7a5de49c
cmake -S /private/tmp/qe_caimd_v1_baseline_src -B /private/tmp/qe_caimd_v1_baseline_build -DQE_ENABLE_MPI=OFF -DQE_ENABLE_OPENMP=OFF -DQE_ENABLE_TEST=OFF
cmake --build /private/tmp/qe_caimd_v1_baseline_build --target qe_pw_exe -j 4
```

Result: succeeded. The executable was
`/private/tmp/qe_caimd_v1_baseline_build/bin/pw.x`.

Traditional Makefile/configure attempt in a temporary source copy:

```text
rsync -a --exclude=.git --exclude='*.o' --exclude='*.mod' --exclude='*.a' --exclude='*.x' /Users/mengqicheng/Downloads/Learning/qe-7.5/ /private/tmp/qe_caimd_v1_make_src/
cd /private/tmp/qe_caimd_v1_make_src
./configure
make pw -j4
```

Result: `./configure` succeeded. `make pw` failed before reaching `PW/src`, in
the pre-existing `LAXlib` preprocessing path on macOS:

```text
cpp ... laxlib.h -o laxlib.fh
cc: error: no input files
make[1]: *** [laxlib.fh] Error 1
make: *** [libla] Error 1
```

This failure occurred before the CAIMD source or `PW/src/Makefile` object list
was compiled.

## Runtime Validation Evidence

All runtime tests used tiny local C-cell inputs with `C.UPF`, `ecutwfc=16 Ry`,
`ecutrho=64 Ry`, Gaussian smearing, `cell_dynamics='pr'`, `press=100 kbar`,
`wmass=20`, `dt=5`, and zero initial ionic velocity
(`ion_temperature='not_controlled'`) for deterministic comparisons.

Validation command: an inline Python harness ran baseline and modified `pw.x`
executables in isolated directories under `/private/tmp/qe_caimd_v1_tests`.

Passed cases: 16 total.

Expected successful cases:

- baseline ordinary `vc-md`;
- modified ordinary `vc-md`;
- explicit `caimd=.false.`;
- `caimd=.true.` with all mask values false;
- fixed `xx` imposed strain;
- fully fixed zero-strain cell with atomic motion;
- uniaxial stride-2 events;
- tensor `xy` shear;
- non-orthogonal cell with fixed `xx` strain.

Expected rejection cases:

- stale existing MD history;
- restart mode;
- `cell_dofree='xy'`;
- `cell_dynamics='w'`;
- `caimd_stride=0`;
- excessive event increment.

Key metrics from the successful validation run:

```text
disabled_vs_baseline_max_cell_bohr: 0.0
explicit_false_vs_default_max_cell_bohr: 0.0
all_false_vs_default_max_cell_bohr: 0.0
fixed_xx_max_actual_minus_nominal: 8.814563662307151e-17
fully_fixed_max_cell_change_bohr: 0.0
stride2_flags: [True, False, True, False, True, False]
stride2_final_events: 3
shear_xy_final_tensor: 0.0004
shear_xy_final_engineering: 0.0008
nonorth_max_abs_xx_strain: 0.0
max_velocity_jump_reported: 1.21217487951527e-19
cases: 16
expected_rejections: existing_history, reject_restart, reject_dofree_xy,
  reject_wentzcovitch, reject_stride0, reject_event_large
```

Static check:

```text
git diff --check
```

Result: passed with no output.

## Problems Found And Fixed During Implementation

- Initial MD-history rejection in `vcsmd` rejected the second step of the same
  from-scratch run. Fixed by only rejecting an existing MD history at
  start-of-run (`istep == 0`) in `vcsmd`.
- QE from-scratch cleanup deletes `tmp_dir/prefix.md` before `vcsmd` opens it.
  Fixed by adding a pre-cleanup stale-history check in `tmpdir_iosys`.
- The first projection rule compounded incremental strain as repeated
  `(I + delta_E)` products, causing small finite-product drift between nominal
  and actual strain. Fixed by targeting cumulative masked strain against the
  initial CAIMD reference cell.
- Deterministic disabled-mode regression required zero initial velocity; random
  velocity initialization produced harmless 1e-7 Bohr trajectory differences
  between independent runs.

## Locality And Cleanup

No network access, Leonardo access, Slurm commands, remote pushes, or dependency
installs were used.

Temporary paths used for evidence capture:

- `/private/tmp/qe_caimd_v1_build`
- `/private/tmp/qe_caimd_v1_baseline_src`
- `/private/tmp/qe_caimd_v1_baseline_build`
- `/private/tmp/qe_caimd_v1_tests`
- `/private/tmp/qe_caimd_v1_make_src`
- `/private/tmp/qe_caimd_v1_make_build`

These temporary paths are removed before the final implementation commit.

## Remaining Risks

- Makefile/configure could not be fully validated on this macOS environment
  because of the unrelated `LAXlib` preprocessing failure.
- CAIMD V1 is intentionally a constrained driven-dynamics mode, not ordinary
  energy-conserving NpT dynamics.
- Scientific production use still needs larger physical tests, timestep/strain
  sensitivity checks, and comparison against the intended stress-strain
  protocol.
