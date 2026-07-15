# CAIMD V1 Implementation Report

Status: SUCCESS

This report records the first working constrained AIMD implementation in this
Quantum ESPRESSO 7.5 source tree. The implementation commit that contains this
file is created after the report is staged; its hash is reported in the final
Codex response.

## macOS Validation Closure Checkpoint

- Status: SUCCESS
- Current validation phase: Complete -- macOS validation closure passed after
  fixing a non-variable-cell CAIMD rejection gap.
- Tests completed:
  - Verified repository root `/Users/mengqicheng/Downloads/Learning/qe-7.5`.
  - Verified branch `main`.
  - Verified `HEAD=df2047eb4387a028d49b4ae419a0836c080c3f83`.
  - Verified commit count `3`.
  - Verified clean working tree and clean index before edits.
  - Read `AGENTS.md`, `CAIMD_IMPLEMENTATION.md`, and this report.
  - Verified `ATOMIC_VELOCITIES { a.u }` syntax and units from
    `PW/Doc/INPUT_PW.txt` and `Modules/read_cards.f90`.
  - Configured and built current serial CMake `qe_pw_exe`.
  - Configured and built baseline serial CMake `qe_pw_exe`.
  - Ran deterministic nonzero-velocity baseline/default/false/all-false
    regressions.
  - Ran fully fixed zero-strain cell with nonzero thermal atomic motion.
  - Ran stride-2 event continuity at full and half event increments.
  - Ran rotation-suppression checks for orthogonal, non-orthogonal, and
    anisotropic cells.
  - Validated the `.caimd` stress-strain stream, restart/reinitialization
    behavior, invalid-input rejections, and CAIMD determinant/inverse/condition
    helper paths.
- Exact latest command and result:
  - The final inline Python validation harness completed successfully with
    `"failures": []`.
- Current modified files:
  - `Modules/read_namelists.f90`
  - `CAIMD_IMPLEMENTATION.md`
  - `CAIMD_V1_IMPLEMENTATION_REPORT.md`
- Defects found:
  - Non-variable-cell inputs can currently hide `caimd=.true.` inside an ignored
    `&CELL` namelist. A tested `calculation='scf'` input with
    `&CELL caimd=.true.` ran to `JOB DONE` instead of failing with the CAIMD
    `vc-md/pr` compatibility error.
- Fixes applied:
  - Added a guarded ignored-`&CELL` probe in `Modules/read_namelists.f90` so
    non-variable-cell inputs containing `caimd=.true.` are rejected explicitly,
    while ordinary `scf` inputs with no `&CELL` and ignored `&CELL
    caimd=.false.` inputs still run.
- Builds completed:
  - Current CMake build succeeded:
    `cmake -S . -B /private/tmp/qe_caimd_mac_current_build -DQE_ENABLE_MPI=OFF -DQE_ENABLE_OPENMP=OFF -DQE_ENABLE_TEST=OFF`
    then
    `cmake --build /private/tmp/qe_caimd_mac_current_build --target qe_pw_exe -j 4`.
  - Baseline CMake build succeeded from
    `a11913ac6cd37f39cf31056e8e4823ec7a5de49c`:
    `cmake -S /private/tmp/qe_caimd_mac_baseline_src -B /private/tmp/qe_caimd_mac_baseline_build -DQE_ENABLE_MPI=OFF -DQE_ENABLE_OPENMP=OFF -DQE_ENABLE_TEST=OFF`
    then
    `cmake --build /private/tmp/qe_caimd_mac_baseline_build --target qe_pw_exe -j 4`.
  - Incremental current rebuild after the rejection fix succeeded:
    `cmake --build /private/tmp/qe_caimd_mac_current_build --target qe_pw_exe -j 4`.
- Known blockers:
  - None. `sysctl -n machdep.cpu.brand_string hw.model hw.optional.arm64` was
    blocked by the sandbox, so hardware metadata will be collected with other
    local commands where available.
- Exact next action:
  - Inspect the final diff, run static checks, clean temporary helper files, and
    create commit 4 `Validate CAIMD V1 on macOS`.

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

## macOS Validation Closure Evidence

Closure validation used local Apple Silicon serial CMake builds only:

- Current executable:
  `/private/tmp/qe_caimd_mac_current_build/bin/pw.x`.
- Baseline executable from `a11913ac6cd37f39cf31056e8e4823ec7a5de49c`:
  `/private/tmp/qe_caimd_mac_baseline_build/bin/pw.x`.
- Temporary test root:
  `/private/tmp/qe_caimd_mac_tests`.

The final inline Python harness completed with `"failures": []`. It verified:

- ordinary baseline/default/explicit-false/all-false regressions with
  deterministic nonzero `ATOMIC_VELOCITIES`;
- fully fixed zero-strain cell: `max_cell_change=0.0`, nonzero atomic motion
  (`max_position_displacement_scaled=0.019814352828091283`), and nonzero
  velocity metric (`min_cart_speed_metric=0.07332121111929343`);
- stride-2 event schedule `[T,F,T,F,T,F,T,F]` with event counts
  `[0,1,1,2,2,3,3,4]`;
- full and half event increments with velocity-jump diagnostics at roundoff
  (`max_velocity_jump <= 1.5332934166833742e-19`);
- orthogonal, non-orthogonal, and anisotropic rotation suppression with
  `max_skew_fro <= 1.4334732947067358e-17`;
- `.caimd` stream shape and units: 8 rows, 26 columns including the logical
  event flag, stress agreement with stdout within `4.955300604549961e-4 GPa`,
  and premove volume agreement within `3.519616029734607e-7 bohr^3`;
- `.caimd` reinitialization after deleting only stale MD history: first run
  wrote 3 rows, second run wrote 2 rows rather than appending;
- invalid input rejections for non-`vc-md` CAIMD, non-PR cell dynamics,
  non-`all`/shape/2Dshape `cell_dofree`, zero/negative stride, oversized event
  increment, NaN/Infinity strain values, restart mode, and stale MD history;
- temporary helper-driver coverage for CAIMD singular inverse and
  ill-conditioned-cell error paths.

During closure, one additional bug was found and fixed: `calculation='scf'`
with `&CELL caimd=.true.` was previously ignored and ran as a normal `scf`.
`Modules/read_namelists.f90` now probes ignored non-variable-cell `&CELL`
namelists and rejects that CAIMD request explicitly while preserving ordinary
`scf` behavior for no-`&CELL` and `caimd=.false.` cases.

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
- macOS closure found that ignored non-variable-cell `&CELL caimd=.true.`
  inputs were not rejected. Fixed in `Modules/read_namelists.f90` by probing
  ignored `&CELL` namelists, rejecting active CAIMD there, and restoring the
  input stream for ordinary non-variable-cell inputs.

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
- `/private/tmp/qe_caimd_mac_current_build`
- `/private/tmp/qe_caimd_mac_baseline_src`
- `/private/tmp/qe_caimd_mac_baseline_build`
- `/private/tmp/qe_caimd_mac_tests`
- `/private/tmp/qe_caimd_mac_probe_vel`
- `/private/tmp/qe_caimd_helper_driver.f90`
- `/private/tmp/qe_caimd_helper_driver`
- `/private/tmp/qe_caimd_condition_driver.f90`
- `/private/tmp/qe_caimd_condition_driver`

The macOS closure temporary paths were removed before commit 4.

## Remaining Risks

- Makefile/configure could not be fully validated on this macOS environment
  because of the unrelated `LAXlib` preprocessing failure.
- CAIMD V1 is intentionally a constrained driven-dynamics mode, not ordinary
  energy-conserving NpT dynamics.
- Scientific production use still needs larger physical tests, timestep/strain
  sensitivity checks, and comparison against the intended stress-strain
  protocol.
