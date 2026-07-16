# CAIMD Leonardo GPU Validation Report

Status: SUCCESS

This file is the persistent local checkpoint and final report for compiling and
validating the QE 7.5 CAIMD V1 implementation as an MPI-enabled NVIDIA GPU build
on the Leonardo Booster partition.

## Final Summary

- Status: SUCCESS.
- Resumed after interruption from the completed remote deployment and CMake
  build state; no source-code changes were made in this validation task.
- Local repository: `/Users/mengqicheng/Downloads/Learning/qe-7.5`, branch
  `main`, starting HEAD `b53f59316802cb8145a72bbe3a1904161bc55826`.
- Remote deployment root:
  `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919`.
- Current source hash:
  `b53f59316802cb8145a72bbe3a1904161bc55826`; baseline source hash:
  `a11913ac6cd37f39cf31056e8e4823ec7a5de49c`.
- Leonardo allocations used: `49490363` on `lrdn0241`, `49522699` on
  `lrdn0064`, and `49525164` on `lrdn2026`. Allocation `49525164` was released
  after validation and `sacct` reported `CANCELLED+`, elapsed `01:00:52`, end
  time `2026-07-16T11:12:03`.
- Build/runtime stack: `profile/chem-phys`, NVHPC `24.5`, HPCX MPI `2.19`,
  CMake `3.27.9`, NVIDIA A100-SXM-64GB, driver `535.274.02`, driver CUDA
  `12.2`, FFTW/OpenBLAS/ScaLAPACK/LibXC from the NVHPC/HPCX module stack.
- Numerical-library workaround: static `libopenblas.a` for BLAS/LAPACK,
  deployment-local stripped `libscalapack.so`, stripped FFTW shared libraries,
  local FFTW SONAME symlinks, local stripped `libopenblas.so.0`, and
  build/install RPATH plus loader `LD_LIBRARY_PATH` pointing at
  `$ROOT/local-libs/lib`.
- Executables:
  - current installed GPU/MPI `pw.x`:
    `$ROOT/install-current-gpu/bin/pw.x`, SHA256
    `4b316dad82d2cb18155fa62ed4957aef041355955b4a7a5687f15fe7fbcd0b35`;
  - current build `pw.x`:
    `$ROOT/build-current-gpu/bin/pw.x`, same SHA256;
  - untouched baseline GPU/MPI `pw.x`:
    `$ROOT/build-baseline-gpu/bin/pw.x`, SHA256
    `7804c655a1ea6c157e09f402ef4064940edea08558317d1b529623fe8825d67f`;
  - isolated traditional `configure && make pw` observation `pw.x`:
    `$ROOT/validation/traditional/current-traditional-20260716T090743Z/src/PW/src/pw.x`,
    SHA256 `86dbe69314ec2a176f09d832d331bf6a6acb31b573c75032b69680646fd2a2be`.
- Linkage inspection: final current and baseline CMake executables resolved
  HPCX MPI, NVHPC CUDA libraries, deployment-local ScaLAPACK, and
  deployment-local FFTW with no missing libraries.
- Valid startup: current one-rank, current two-rank, and baseline one-rank SCF
  inputs all completed with `JOB DONE`, zero return code, and no MPI abort.
- A100 activity:
  - current `gpu_scf8` one-rank run: 7 samples, 1190 MiB peak memory, 2% peak
    GPU utilization, 18.30 W power delta;
  - active CAIMD `stride2_xx` two-rank run: 20 samples, 2334 MiB peak memory,
    98% peak GPU utilization, 37.54 W power delta.
- MPI launcher: one-rank runs use Slurm `srun --mpi=pmix_v3`; two-rank QE runs
  use one Slurm task wrapping `mpiexec --oversubscribe -n 2`, which avoided the
  direct two-task QE reachability failure seen during harness development.
- Persistent validation matrix: final harness checksum
  `5461da3d9070b8fafa28e4f85643e844b46c23d94ea22c8244f3f51a1ff5bef2`;
  34 cases, zero failures. Final summary:
  `$ROOT/validation/analysis/summary.json`, SHA256
  `03401cc13b390d1a907d66221bd9e48de9e02b46abe081db906f22d1db85b53c`.
- Pseudopotential provenance: `src-current/pseudo/C.UPF` copied to
  `validation/pseudopotentials/C.UPF`, element C, SHA256
  `9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639`.
- Baseline/no-op regression: structured post-processing in
  `$ROOT/validation/analysis/structured_regression.json` found exact equality
  for baseline vs current default, current default vs `caimd=.false.`, and
  current default vs all-false mask across parsed total energies, total force
  norms, SCF correction norms, total stress, printed cells, positions,
  volumes, `.ave/.avec/.e/.eal/.p/.tv`, and `tmp/*.md`.
- Current one-rank vs two-rank default `vc-md`: small MPI floating-point
  differences only; max_abs total energy `4.785e-05 Ry`, stress `0.81 kbar`,
  cell `7.47e-07 bohr`, positions `8.45e-08 crystal`, MD history
  `2.40e-04`.
- Ordinary current-QE smoke tests passed: SCF without `&CELL`, SCF with ignored
  `&CELL caimd=.false.`, fixed-cell MD, relax, vc-relax, and ordinary vc-md.
- Parser observations were recorded without fixes: ignored `&CELL` after cards,
  `&CELL caimd=.true.` in SCF, and malformed after-card `&CELL` cases abort
  through QE input/error paths; deferred parser issues remain deferred.
- Active CAIMD tests passed for fixed zero normal strain, fully fixed zero
  strain with nonzero atomic motion, stride-2 uniaxial strain-and-hold, tensor
  `xy` shear, non-orthogonal cell, anisotropic cell, and a two-rank GPU active
  case. `.caimd` files had 26 columns, no duplicate rows, finite temperature,
  finite velocity residuals, and MPI single-writer behavior.
- Rotation metrics from active cases show antisymmetric increments at numerical
  noise: fully fixed zero strain had max `E` and `W` both `0`; stride-2
  two-rank GPU had max `||W||_F=3.80e-08` and max skew-to-strain ratio
  `1.82e-04`; shear `xy` had max `||W||_F=6.76e-09` and ratio `3.42e-05`.
- Representative MPI rejection tests passed: non-`vc-md`, non-PR cell
  dynamics, non-`all` `cell_dofree`, restart, stale MD history, zero stride,
  and negative stride all terminated with nonzero status and no misleading
  successful `.caimd` output.
- Traditional build observation: isolated
  `$ROOT/validation/traditional/current-traditional-20260716T090743Z/src`
  `./configure MPIF90=mpif90 CC=mpicc F77=mpifort FC=mpifort --enable-parallel`
  followed by `make pw -j8` completed with `traditional_build_rc=0`. Log SHA256
  is `f2528fd4d8340e860c52f5492f5bc424f0f4c2e73505f1f81d706785fddcf44c`.
- No push was performed. Deferred CAIMD architectural review issues remain
  open and were not fixed during this validation task.

## Live Checkpoint

- Status: SUCCESS
- Current phase: all requested Leonardo GPU/MPI validation, traditional build
  observation, documentation, and local commit preparation completed.
- Local repository path: `/Users/mengqicheng/Downloads/Learning/qe-7.5`
- Local HEAD: `b53f59316802cb8145a72bbe3a1904161bc55826`
- Code-bearing CAIMD source commit: `7aa61919fe1b23ec4e80dbbccae550cde5c81509`
- Deferred-issue documentation commit:
  `b53f59316802cb8145a72bbe3a1904161bc55826`
- Remote deployment root:
  `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919`
- Current remote source commit: `b53f59316802cb8145a72bbe3a1904161bc55826`
- Baseline source commit: `a11913ac6cd37f39cf31056e8e4823ec7a5de49c`
- Current Slurm allocation or last allocation ID: `49525164` (released after
  validation; `sacct` state `CANCELLED+`);
  previous allocation `49522699` failed after invalid empty-stdin `pw.x` smoke
  attempt.
- Node hostname: `lrdn2026.leonardo.local`
- Module stack: `profile/chem-phys`, `nvhpc/24.5`, `hpcx-mpi/2.19`,
  `fftw/3.3.10--hpcx-mpi--2.19--nvhpc--24.5`,
  `libxc/7.0.0--nvhpc--24.5`, `openblas/0.3.26--nvhpc--24.5`,
  `netlib-scalapack/2.2.0--hpcx-mpi--2.19--nvhpc--24.5`,
  `cmake/3.27.9`
- Current build status: current GPU/MPI `pw.x` built successfully and was
  relinked successfully with local ScaLAPACK and FFTW soname dependencies.
- Baseline build status: baseline GPU/MPI `pw.x` built successfully and
  inspected.
- Latest resume note, 2026-07-16:
  - Local preflight still matches the requested resume state: branch `main`,
    HEAD `b53f59316802cb8145a72bbe3a1904161bc55826`, commit count `5`, no
    tracked or staged changes, and only this report untracked.
  - Allocation `49525164` was running on `lrdn2026` during resumed validation
    and was released after all remote validation completed.
  - The first persistent validation harness run created
    `validation/analysis/summary.json` and
    `evidence/tests/caimd_leonardo_validation_full.log`.
  - Verified from that run: current and baseline one-rank valid SCF startup
    both completed with `JOB DONE`; current one-rank A100 activity was observed
    on `gpu_scf8` with 7 monitor samples, 1188 MiB peak GPU memory, 9% peak GPU
    utilization, and 13.49 W power delta; one-rank baseline/default,
    explicit-false, and all-false `vc-md` total energies matched exactly in
    the available output-energy comparison.
  - The first failed cases were
    `current_startup_scf_2r`, `current_vcmd_default_2r`,
    `ordinary_md_fixed`, `active_fixed_xx_zero_1r`,
    `active_stride2_xx_1r`, `active_shear_xy_1r`,
    `active_stride2_xx_2r_gpu`, `active_reinit_first`, and
    `active_reinit_second`.
  - Inspection showed these were validation-harness/input issues rather than a
    QE source defect: direct two-rank QE under `srun` hit an HPCX reachability
    problem, while a one-task Slurm step wrapping
    `mpiexec --oversubscribe -n 2` worked; fixed-cell `calculation='md'`
    cannot use `ion_dynamics='beeman'`; active small-cell CAIMD cases require
    `nosym=.true., noinv=.true.` to avoid symmetry checks during imposed shear
    or anisotropic strain; the nonorthogonal case needed a looser/longer SCF
    setup before dynamics.
  - A patched validation harness exists locally at
    `/private/tmp/caimd_leonardo_validation.py` and is ready to upload. It uses
    `mpiexec --oversubscribe -n 2` inside a one-task Slurm step for QE
    two-rank cases, uses `verlet` for fixed-cell MD, adds `nosym/noinv` for
    active CAIMD inputs, and relaxes SCF thresholds for the small validation
    cells.
  - The first patched rerun exited before the test matrix because the harness
    stdout fallback parser treated unrelated three-number output lines as cell
    matrices and hit a singular matrix in the rotation metric. The harness has
    been corrected locally to parse only explicit `CELL_PARAMETERS` blocks from
    stdout and to skip singular matrix pairs in rotation analysis. This is a
    validation parser fix, not a QE source change.
  - The second patched rerun completed the 34-case matrix with only
    `ordinary_md_fixed` failing. The output showed fixed-cell MD entered
    dynamics and moved atoms, then aborted in `checkallsym` because the
    deterministic input velocities did not preserve the original symmetry
    operations. The fixed-cell MD smoke input has been corrected locally to use
    `nosym=.true., noinv=.true.`; this is a validation-input correction, not a
    QE source change.
  - The third patched rerun completed successfully with zero failures across
    34 cases. The final machine summary is
    `validation/analysis/summary.json` with SHA256
    `03401cc13b390d1a907d66221bd9e48de9e02b46abe081db906f22d1db85b53c`.
  - Structured no-op regression post-processing was saved as
    `validation/analysis/structured_regression.json`. It compares total
    energies, total force norms, SCF correction norms, total stress, printed
    cell matrices, printed atomic positions, volumes, QE MD history files
    `.ave/.avec/.e/.eal/.p/.tv`, and `tmp/*.md`. Baseline vs current default,
    current default vs `caimd=.false.`, and current default vs all-false mask
    are exactly identical for all parsed quantities. Current one-rank vs
    two-rank differs only at small MPI floating-point levels, e.g. total-energy
    max_abs `4.785e-05 Ry`, stress max_abs `0.81 kbar`, cell max_abs
    `7.47e-07 bohr`, and MD history max_abs `2.40e-04`.
- Completed tests:
  - Local preflight verified the expected starting state before documentation
    edits: branch `main`, HEAD
    `7aa61919fe1b23ec4e80dbbccae550cde5c81509`, commit count `4`, clean working
    tree, clean index, and no required Git remote.
  - Read `AGENTS.md`, `CAIMD_IMPLEMENTATION.md`, and
    `CAIMD_V1_IMPLEMENTATION_REPORT.md`.
  - Recorded deferred post-cluster review issues in
    `CAIMD_V1_IMPLEMENTATION_REPORT.md`.
  - Created documentation-only commit
    `b53f59316802cb8145a72bbe3a1904161bc55826`
    (`Record deferred CAIMD review issues`).
  - Created remote deployment directory skeleton and remote phase marker
    `phase00_local_preflight.done`.
  - Transferred Git bundle
    `qe-7.5-caimd-v1-main-b53f593.bundle` to the deployment root and verified
    SHA256 `94e3c394ab649dbb7afa813969294f630d8cc3ba8cd849476ed319da63ed9b0a`.
  - Cloned remote `src-current` at
    `b53f59316802cb8145a72bbe3a1904161bc55826`.
  - Cloned remote `src-baseline` at
    `a11913ac6cd37f39cf31056e8e4823ec7a5de49c`.
  - Verified both remote source trees have clean Git status and no `.DS_Store`
    or local compiled build artifacts.
  - Created remote phase marker `phase01_deployment.done`.
  - Recorded reference QE 7.4.1 module environment and module show output.
  - Verified selected build stack: NVHPC 24.5 compilers, HPCX MPI 2.19 wrappers
    using `nvfortran`/`nvc`, CMake 3.27.9, NVHPC bundled CUDA compiler 12.4,
    NVHPC-built FFTW/OpenBLAS/ScaLAPACK/LibXC.
  - Created remote phase marker `phase02_environment.done`.
  - Obtained interactive Booster GPU allocation `49490363` on
    `lrdn0241.leonardo.local`.
  - Recorded GPU visibility: `NVIDIA A100-SXM-64GB`,
    `CUDA_VISIBLE_DEVICES=0`, NVIDIA driver `535.274.02`, driver CUDA `12.2`.
  - Inspected QE 7.5 CMake GPU/MPI options and NVHPC GPU architecture handling.
  - Configured current GPU/MPI build successfully with NVHPC CUDA flags
    `-cuda;-gpu=cc80`, MPI enabled, OpenACC/OpenMP enabled, ScaLAPACK enabled,
    LibXC enabled, FFTW3 vendor selected, static OpenBLAS for BLAS/LAPACK, and
    `QE_GPU_ARCHS=sm_80`.
  - Obtained replacement interactive Booster GPU allocation `49522699` on
    `lrdn0064.leonardo.local` after allocation `49490363` ended.
  - Recorded resumed allocation evidence in
    `evidence/environment/phase02_resume_allocation_49522699.log`: one
    `NVIDIA A100-SXM-64GB`, `CUDA_VISIBLE_DEVICES=0`, NVIDIA driver
    `535.274.02`, driver CUDA `12.2`, and the selected NVHPC/HPCX module stack.
  - Backed up remote build scripts with suffix
    `.pre_stripped_scalapack.20260716T075306Z`, then updated
    `load_build_environment.sh`, `configure_current_gpu.sh`, and
    `configure_baseline_gpu.sh` to use
    `$ROOT/local-libs/lib/libscalapack.so`, record script checksums, preserve
    prior configure command/cache evidence, and set build/install RPATH.
  - Reconfigured the current GPU/MPI build successfully in allocation
    `49522699`; CMake reported `Found SCALAPACK` at the deployment-local
    stripped `libscalapack.so` plus static OpenBLAS.
  - Confirmed from `PW/CMakeFiles/qe_pw_exe.dir/link.txt` that the ScaLAPACK
    path had moved to the local stripped copy and that the next failing shared
    module libraries were FFTW (`libfftw3.so` and `libfftw3_omp.so`).
  - Recorded FFTW module metadata in
    `evidence/environment/fftw_module_show_49522699.txt`; the module path is
    `/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/nvhpc-24.5/fftw-3.3.10-ncl5ohbh5nie4oxbsddysah3osrwm3zf`.
  - Backed up the ScaLAPACK-workaround scripts with suffix
    `.pre_stripped_fftw.20260716T080145Z`, then updated the build scripts to
    create stripped local copies of `libfftw3.so` and `libfftw3_omp.so`, pass
    them through `FFTW3_DOUBLE` and `FFTW3_DOUBLE_OPENMP`, and record script
    checksums in `evidence/builds/script_update_stripped_fftw.sha256`.
  - Reconfigured the current GPU/MPI build successfully with CMake reporting
    `Found FFTW3` at the deployment-local stripped `libfftw3.so` and
    `libfftw3_omp.so`.
  - Built current GPU/MPI `pw.x` successfully; `phase03_current_build.done`
    exists, and build/install executable SHA256 is
    `00e0b11f1f3e133eb9d39ceeff89cd65f6c0738c9e52f3eab2968d608394ef21`.
  - Inspected current executable: it is an x86-64 ELF linked with MPI and
    NVHPC CUDA runtime libraries (`libcufft`, `libcurand`, `libcusolver`,
    `libcublas`, and `libcudart` appear in `ldd`/`readelf` evidence).
  - The current executable started far enough to print the QE 7.5 banner,
    current Git commit `b53f59316802cb8145a72bbe3a1904161bc55826`, MPI/OpenMP
    process layout, and memory report before the empty-input namelist error.
  - Backed up configure scripts with suffix
    `.pre_local_soname_deps.20260716T080822Z`, then updated them to create a
    stripped local `libopenblas.so.0`, FFTW soname symlinks
    `libfftw3.so.3`/`libfftw3_omp.so.3`, and `-Wl,-rpath-link,$ROOT/local-libs/lib`;
    checksum evidence is in
    `evidence/builds/script_update_local_soname_deps.sha256`.
  - Obtained replacement allocation `49525164` on
    `lrdn2026.leonardo.local` and recorded node/GPU/module evidence in
    `evidence/environment/phase04_allocation_49525164.log`.
  - Reconfigured and relinked current GPU/MPI `pw.x` successfully with
    `-Wl,-rpath-link,$ROOT/local-libs/lib`; the final link no longer emitted
    the earlier OpenBLAS/FFTW/ScaLAPACK debug-section linker warnings.
  - Inspected the relinked current executable in
    `evidence/builds/current_pw_x_soname_inspection.log`; SHA256 is
    `4b316dad82d2cb18155fa62ed4957aef041355955b4a7a5687f15fe7fbcd0b35`,
    ScaLAPACK and FFTW sonames resolve through `local-libs`, and CUDA/MPI
    dependencies are visible in `ldd`/`readelf`.
  - Configured and built untouched baseline GPU/MPI `pw.x` successfully from
    `a11913ac6cd37f39cf31056e8e4823ec7a5de49c`; phase marker
    `evidence/state/phase04_baseline_build.done` exists.
  - Performed final executable inspection in
    `evidence/builds/final_executable_inspection.log`: current build and
    install checksums match
    `4b316dad82d2cb18155fa62ed4957aef041355955b4a7a5687f15fe7fbcd0b35`;
    baseline checksum is
    `7804c655a1ea6c157e09f402ef4064940edea08558317d1b529623fe8825d67f`;
    both executables resolve HPCX MPI, NVHPC CUDA libraries, local
    ScaLAPACK, and local FFTW with no missing libraries.
  - Installed persistent remote validation harness
    `scripts/caimd_leonardo_validation.py`; syntax check passed and checksum
    evidence is in `evidence/tests/caimd_leonardo_validation_script.sha256`.
- Failed tests:
  - The attempted current startup smoke command
    `timeout 30s "$ROOT/install-current-gpu/bin/pw.x" -h < /dev/null` was
    invalid for this build: `pw.x` entered normal stdin parsing, failed because
    no `&CONTROL` namelist was present, invoked `MPI_ABORT`, and ended Slurm
    allocation `49522699` with state `FAILED` and exit code `0:9`.
- Failures encountered:
  - Phase 1 deployment script initially failed before cloning because
    `git bundle verify` was run outside a Git repository on Leonardo and
    returned `error: need a repository to verify a bundle`.
  - Phase 2 reference-environment script initially hung at
    `pw.x --version`; the remote `pw.x` process and wrapper script were killed,
    and the probe was changed to `timeout 20s pw.x --version < /dev/null ||
    true`.
  - First current GPU CMake configure failed because `FindBLAS/FindLAPACK`
    did not auto-detect the loaded NVHPC OpenBLAS module. OpenBLAS and
    ScaLAPACK library paths were inspected and the configure scripts were
    updated to pass `BLAS_LIBRARIES`, `LAPACK_LIBRARIES`, and
    `SCALAPACK_LIBRARIES` explicitly.
  - Second current GPU CMake configure found the OpenBLAS shared library but
    failed at the LAPACK `zhpev` link check because `/usr/bin/ld` could not read
    compressed debug sections in `libopenblas.so`. A direct link test confirmed
    shared OpenBLAS failed while static `libopenblas.a` linked successfully.
    Configure scripts were updated to use static OpenBLAS for BLAS/LAPACK.
  - Current GPU build then compiled through `qe_pw` and failed only at final
    `pw.x` link because `/usr/bin/ld` could not read compressed debug sections
    in module `libscalapack.so`.
  - A deployment-root copy of `libscalapack.so` stripped with
    `objcopy --strip-debug` linked successfully in a small MPI/ScaLAPACK test
    when paired with static OpenBLAS.
  - After switching ScaLAPACK to the stripped local copy, the next final link
    failed on module `libfftw3.so` with the same compressed debug section
    `ld` error.
  - The successful current final link still emitted non-fatal linker warnings
    because the local stripped `libscalapack.so` has a dynamic dependency on
    `libopenblas.so.0`; this will be cleaned up by adding a stripped local
    `libopenblas.so.0` and relinking.
  - The invalid empty-stdin current startup smoke triggered `MPI_ABORT` and
    ended allocation `49522699`; subsequent validation must use a valid input
    file with `-inp`, never an empty stdin help probe.
- Source changes made:
  - Documentation-only deferred-issue record in
    `CAIMD_V1_IMPLEMENTATION_REPORT.md`.
- Exact latest successful command:
  - `ssh -o BatchMode=yes -o ConnectTimeout=15 leonardo 'python3 - <<PY ...'`
    inspected the first persistent validation summary and confirmed the
    failures above plus the valid one-rank startup/GPU/no-op evidence.
- Exact failed command, if any:
  - In allocation job `49490363`,
    `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/scripts/build_current_gpu.sh`
  - In allocation job `49522699`,
    `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/scripts/build_current_gpu.sh`
    failed at final link on module `libfftw3.so`.
  - In allocation job `49522699`,
    `timeout 30s "$ROOT/install-current-gpu/bin/pw.x" -h < /dev/null`
    was an invalid smoke probe and ended the allocation via `MPI_ABORT`.
- Exact next action:
  - No remote validation action remains. Complete local documentation review,
    run local Git checks, create the local validation commit, and do not push.
- Relevant log paths:
  - Local report: `CAIMD_LEONARDO_GPU_VALIDATION_REPORT.md`
  - Remote checkpoint, after Phase 1 begins:
    `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/evidence/REMOTE_CHECKPOINT.md`
  - Phase 1 deployment log:
    `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/evidence/state/phase01_deployment.log`
  - Phase 2 reference environment log:
    `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/evidence/environment/phase02_reference_environment.log`
  - Phase 2 build environment log:
    `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/evidence/environment/phase02_build_environment.log`
  - Phase 3 allocation log:
    `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/evidence/environment/phase03_allocation.log`
  - Current configure log:
    `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/evidence/builds/current_configure.log`

## Local Environment

- Checkpoint initialized: `2026-07-15T18:19:25Z`
- Local OS: macOS 26.5.2 build 25F84
- Kernel: `Darwin chengmengqideMacBook-Air-2.local 25.5.0 Darwin Kernel Version 25.5.0: Tue Jun  9 22:26:46 PDT 2026; root:xnu-12377.121.10~1/RELEASE_ARM64_T8103 arm64`
- Git remote policy: no push; no remote Git repository is required.

## Planned Remote Layout

```text
/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/
├── src-current/
├── src-baseline/
├── build-current-gpu/
├── build-baseline-gpu/
├── install-current-gpu/
├── validation/
│   ├── inputs/
│   ├── runs/
│   ├── analysis/
│   └── pseudopotentials/
├── scripts/
└── evidence/
    ├── REMOTE_CHECKPOINT.md
    ├── state/
    ├── environment/
    ├── builds/
    └── tests/
```

## Planned Build Matrix

- Current source: local `main` after the documentation-only checkpoint commit.
- Baseline source: untouched QE 7.5 commit
  `a11913ac6cd37f39cf31056e8e4823ec7a5de49c`.
- Target executable: `pw.x`.
- Required environment: coherent Leonardo NVHPC + HPCX MPI + CUDA stack, with
  FFTW, BLAS/LAPACK, ScaLAPACK, and CMake selected after inspecting the
  installed QE 7.4.1 module.
- Build mode: out-of-tree CMake MPI/GPU builds; traditional configure/make is a
  secondary evidence item after CMake builds succeed.

## Deferred Issues

Deferred source-scan concerns are recorded in
`CAIMD_V1_IMPLEMENTATION_REPORT.md` under
`Deferred Post-Cluster Review Issues`. They are intentionally not fixed during
this initial Leonardo build and smoke-validation task. Cluster observations will
be used to scope a later dedicated task, and successful cluster validation does
not close those issues.
