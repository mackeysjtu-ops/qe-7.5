# Leonardo CAIMD Build And Validation Manual

This document is the repeatable Leonardo GPU/MPI workflow for the QE 7.5 CAIMD
V1 source tree. It is a procedure manual, not a historical validation report.
Historical evidence belongs in `CAIMD_LEONARDO_GPU_VALIDATION_REPORT.md`.

## Canonical Remote Root

The validated Leonardo deployment root is:

```bash
ROOT=/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919
```

Expected layout:

```text
$ROOT/
├── src-current/
├── src-baseline/
├── build-current-gpu/
├── build-baseline-gpu/
├── install-current-gpu/
├── local-libs/
├── scripts/
├── validation/
│   ├── inputs/
│   ├── runs/
│   ├── analysis/
│   └── pseudopotentials/
└── evidence/
    ├── REMOTE_CHECKPOINT.md
    ├── state/
    ├── environment/
    ├── builds/
    └── tests/
```

Do not overwrite completed validation evidence, baseline source/build trees,
or local stripped libraries unless a recovery step explicitly requires it and
the old file has been archived.

## Environment

All build and run commands should source the loader:

```bash
source "$ROOT/scripts/load_build_environment.sh"
```

The validated module family is:

```text
profile/chem-phys
nvhpc/24.5
hpcx-mpi/2.19
fftw/3.3.10--hpcx-mpi--2.19--nvhpc--24.5
libxc/7.0.0--nvhpc--24.5
openblas/0.3.26--nvhpc--24.5
netlib-scalapack/2.2.0--hpcx-mpi--2.19--nvhpc--24.5
cmake/3.27.9
```

The loader also sets:

```bash
export CAIMD_REMOTE_ROOT="$ROOT"
export LD_LIBRARY_PATH="$ROOT/local-libs/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OMPI_MCA_coll_hcoll_enable=0
```

Run the Python validation harness with the Leonardo Python 3.11 environment:

```bash
source "$HOME/venvs/mace_py311/bin/activate"
python --version
```

The system/default `python3` can be too old for the tracked harness.

## Checkpoints

Before any operation expected to run longer than about five minutes, update:

- local checkpoint section in `CAIMD_LEONARDO_GPU_VALIDATION_REPORT.md`;
- remote checkpoint `$ROOT/evidence/REMOTE_CHECKPOINT.md`.

Each checkpoint should record:

- local HEAD and working-tree state;
- remote source HEAD;
- phase;
- modified files;
- latest successful command;
- latest failed command, if any;
- Slurm allocation and node, if applicable;
- evidence/log paths;
- next action.

## Local Source Rule

All source, documentation, and harness changes must be tracked locally before
remote deployment. Do not edit the remote source tree directly.

Preferred local pre-deployment checks:

```bash
git status --short
git diff --check
python3 -m py_compile test-suite/caimd/caimd_validation.py
python3 test-suite/caimd/caimd_validation.py --help
python3 test-suite/caimd/caimd_validation.py --list-cases
```

Build or smoke-test locally first when practical.

## Deploy Current Source With Git Bundle

Use a Git bundle or patch generated from the local repository. Do not deploy
source with broad `rsync` for CAIMD validation.

Local:

```bash
git status --short
git rev-parse HEAD
git bundle create /private/tmp/qe-caimd-current.bundle HEAD ^a11913ac6cd37f39cf31056e8e4823ec7a5de49c
scp -o BatchMode=yes -o ConnectTimeout=15 /private/tmp/qe-caimd-current.bundle leonardo:/tmp/qe-caimd-current.bundle
```

Remote:

```bash
ROOT=/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919
cd "$ROOT/src-current"
git status --short
git fetch /tmp/qe-caimd-current.bundle HEAD:refs/heads/caimd-deployed
git checkout caimd-deployed
git rev-parse HEAD
git status --short
```

The remote `src-current` tree must match the exact local commit under test.
Do not modify `src-baseline` during current-source deployment.

## Acquire And Record A GPU Allocation

Use an A100-capable allocation consistent with the current project/account
policy. After allocation starts, record:

```bash
hostname
date -u
module list
nvidia-smi
git -C "$ROOT/src-current" rev-parse HEAD
```

Save this evidence under `$ROOT/evidence/environment/`.

## Configure Current GPU Build

Use the existing configure helper:

```bash
"$ROOT/scripts/configure_current_gpu.sh"
```

The helper is expected to:

- source the loader;
- create `$ROOT/build-current-gpu` and `$ROOT/install-current-gpu`;
- preserve and reuse `$ROOT/local-libs/lib`;
- copy and strip module ScaLAPACK/OpenBLAS/FFTW shared libraries when needed;
- configure CMake with NVHPC MPI wrappers, CUDA/OpenACC, LibXC, FFTW3,
  ScaLAPACK, static OpenBLAS for BLAS/LAPACK, and RPATH/RPATH-link to
  `$ROOT/local-libs/lib`;
- write command, cache, and checksum evidence under `$ROOT/evidence/builds/`.

If the configure command changes, archive the old command and CMake cache under
`$ROOT/evidence/builds/` before replacing them.

## Incremental Build Current `pw.x`

Inside the allocation:

```bash
source "$ROOT/scripts/load_build_environment.sh"
cmake --build "$ROOT/build-current-gpu" --target qe_pw_exe -j8
```

Before installing or replacing the runnable binary, archive the previous one:

```bash
ts=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$ROOT/evidence/builds/archived-binaries"
if [ -x "$ROOT/install-current-gpu/bin/pw.x" ]; then
  cp -p "$ROOT/install-current-gpu/bin/pw.x" \
    "$ROOT/evidence/builds/archived-binaries/pw.x.pre-current-$ts"
  sha256sum "$ROOT/evidence/builds/archived-binaries/pw.x.pre-current-$ts" \
    > "$ROOT/evidence/builds/archived-binaries/pw.x.pre-current-$ts.sha256"
fi
cmake --install "$ROOT/build-current-gpu" --prefix "$ROOT/install-current-gpu"
sha256sum "$ROOT/install-current-gpu/bin/pw.x" \
  > "$ROOT/evidence/builds/current_install_pw_x.sha256"
```

Inspect dynamic dependencies:

```bash
ldd "$ROOT/install-current-gpu/bin/pw.x"
readelf -d "$ROOT/install-current-gpu/bin/pw.x"
```

Save dependency evidence under `$ROOT/evidence/builds/`.

## Optional CP Build

If CAIMD boundary changes touch CP-visible parsing, build the CP executable when
practical:

```bash
source "$ROOT/scripts/load_build_environment.sh"
cmake --build "$ROOT/build-current-gpu" --target qe_cpv_exe -j8
```

Record whether `cp.x` was built, skipped, or unavailable. CAIMD remains PW-only;
CP validation is limited to ordinary CP smoke where feasible and explicit
`caimd=.true.` rejection.

## Startup Smoke

Do not use empty stdin or `pw.x -h` as the GPU smoke. Some MPI/GPU builds enter
normal input parsing and can abort the allocation on empty stdin.

Use a valid input generated by the persistent harness, or a minimal valid SCF
input with `-inp`.

## Persistent Validation Harness

The tracked harness is:

```bash
test-suite/caimd/caimd_validation.py
```

Local help/list checks:

```bash
python3 test-suite/caimd/caimd_validation.py --help
python3 test-suite/caimd/caimd_validation.py --list-cases
```

Leonardo run from the deployed source:

```bash
cd "$ROOT/src-current"
source "$ROOT/scripts/load_build_environment.sh"
source "$HOME/venvs/mace_py311/bin/activate"
python test-suite/caimd/caimd_validation.py \
  --mode leonardo \
  --allocation "$SLURM_JOB_ID" \
  --current-exe "$ROOT/install-current-gpu/bin/pw.x" \
  --baseline-exe "$ROOT/build-baseline-gpu/bin/pw.x" \
  --cp-exe "$ROOT/install-current-gpu/bin/cp.x" \
  --pseudo "$ROOT/src-current/pseudo/C.UPF" \
  --work-root "$ROOT/validation"
```

The harness:

- generates deterministic PW and optional CP inputs;
- copies a local pseudopotential and refuses missing pseudopotentials;
- creates an isolated timestamped validation root;
- supports `--cases`, `--groups`, and `--list-cases`;
- handles expected rejections;
- writes `analysis/summary.json`;
- records executable checksums, commands, ranks, logs, paths, and diagnostics;
- does not download data.

Important groups:

```text
A_startup_gpu
B_noop
C_ordinary
D_parser
E_cp
F_active
G_policies
H_caimd_output
I_rejections
J_rotation
```

## Evidence To Record

For each final validation pass, record:

- local commit hash deployed;
- remote `src-current` hash;
- current `pw.x` checksum;
- baseline `pw.x` checksum;
- optional `cp.x` checksum;
- module list and GPU/node evidence;
- configure/build logs;
- harness command;
- harness `summary.json` path and checksum;
- pass/fail case list;
- skipped cases and why.

## Recovery Notes

If configure fails:

- inspect the CMake log and cache first;
- do not delete the build tree unless a clean configure is explicitly needed;
- archive the old configure command/cache before changing flags.

If link fails on module libraries:

- verify `local-libs` stripped copies and SONAME symlinks;
- verify `-Wl,-rpath-link,$ROOT/local-libs/lib`;
- keep original module-library checksums in evidence.

If a Slurm allocation fails:

- record `sacct` state and exit code;
- do not call a failed submission a failed scientific validation unless QE ran
  a valid input and produced interpretable output.

If validation fails:

- inspect the exact case output in the isolated run root;
- separate build/runtime failures from expected scientific or parser
  rejections;
- update both checkpoints before retrying a long build or run.

## Completion Policy

A Leonardo revalidation is complete only when:

- local source state is committed and identified;
- remote source matches that commit;
- current executable was rebuilt or the reason for not rebuilding is recorded;
- dependency inspection has no missing libraries;
- the persistent harness finished with no unexpected failures, or all remaining
  failures are explicitly documented as unresolved;
- local and remote evidence paths are recorded.

Never push from this workflow unless the user explicitly requests a push.
