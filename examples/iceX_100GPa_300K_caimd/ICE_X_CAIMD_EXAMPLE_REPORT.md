# Ice X CAIMD Example Report

Status: PARTIAL

The example package, bounded Leonardo convergence tests, short pilot jobs,
analysis scripts, and final manual-submit deformation templates are prepared.
The remaining limitation is deployment/runtime validation: historical Leonardo
pilots were collected with an earlier CAIMD executable whose PW
`vc-md`/Parrinello-Rahman initialization path ignored explicit
`ATOMIC_VELOCITIES` for variable-cell MD. This source tree now contains a local
source-level fix and baseline/current validation for explicit PR velocity
initialization; the rebuilt Leonardo executable and short Ice X runtime
validation are still pending.

## Local Preflight

- Repository root: `/Users/mengqicheng/Downloads/Learning/qe-7.5`
- Branch: `main`
- Requested starting HEAD:
  `0282a6e6f4e7d330f68ef4a2732777d8343fd716`
- Starting status: no tracked or staged diffs; only untracked local source file
  `iceX.data`.
- Git remotes: none.
- Local source file: `/Users/mengqicheng/Downloads/Learning/qe-7.5/iceX.data`
- Source SHA256:
  `103ccf1fe76e3b3a737914a4eeddb4016d6fb8b308129516541ddc5de08aa6c6`
- QE/CAIMD source modifications: none.
- Push: none.

## Structure Conversion

- Converter: `tools/convert_icex_data.py`
- Original copy: `structure/iceX_original.data`
- QE structure card: `structure/iceX_qe_structure.inc`
- Structure validation: `structure/structure_validation.json` and
  `structure/structure_validation.md`

Validated facts:

- 144 atoms total.
- 48 O and 96 H atoms.
- Type mapping: LAMMPS type 1=O, type 2=H from Masses-section comments and
  masses.
- Cell volume: `451.706300063683 Angstrom^3`.
- LAMMPS tilt factors: `xy=0`, `xz=0`, `yz=0`.
- Maximum roundtrip error: `8.881784197001252e-16 Angstrom`.
- Minimum distances: H-H `1.880901 Angstrom`, H-O `1.151811 Angstrom`,
  O-O `2.303624 Angstrom`.
- No rotation, symmetrization, relaxation, or atom-order changes.

## Initial Velocities

- Generator: `scripts/generate_300K_velocities.py`
- Velocity card: `structure/iceX_velocities_300K.inc`
- Metadata: `structure/iceX_velocities_300K.json`
- Seed: `20260718`
- QE-style degrees of freedom: `3N-3 = 429`
- Final temperature by the generator's QE kinetic-energy formula:
  `299.9999999999999 K`
- Maximum residual COM velocity after rescale:
  `1.466377223952049e-19`

Historical limitation: Leonardo active CAIMD pilots reported max `.caimd`
temperatures of only `0.0035-0.0093 K`. Source inspection and untouched QE 7.5
baseline tests show that old `vcsmd -> vcinit` did not receive the parsed
velocity card and instead initialized through `ranv(temp=0, ...)` when active
CAIMD required `ion_temperature='not_controlled'`.

## 300 K vc-md velocity initialization checkpoint

- Status: `LOCAL FIX VALIDATED; LEONARDO VALIDATION PENDING`
- Starting HEAD:
  `20a9f2f55f9164cdcc5777f43f02891baaa8bfee`
- Current HEAD:
  `20a9f2f55f9164cdcc5777f43f02891baaa8bfee`
- Current phase: local clean current/baseline builds and baseline/current
  comparisons completed; Leonardo deployment and short Ice X validation are
  next.
- Root-cause classification: `D`, an untouched QE 7.5 PW
  `vc-md`/Parrinello-Rahman initialization omission. The velocity card is
  parsed and fixed-cell MD consumes it, but `vcsmd -> vcinit` does not receive
  `tv0rd` or `dynamics_module%vel` and instead calls `ranv(temp, ...)`.
- Current executable path:
  `/private/tmp/qe_velocity_current_build/bin/pw.x`
- Current executable SHA256:
  `aa988325721d7438abe38fe84d27d648957c4fc17711e185689ca6f7d7ebdb20`
- Baseline executable path:
  `/private/tmp/qe_velocity_baseline_build/bin/pw.x`
- Baseline executable SHA256:
  `c5df6bee268f8b43549e9d16ffae264d45474668aece91739fea1f0565506952`
- Tests completed:
  local tiny fixed-cell MD explicit velocities, ordinary `vc-md` explicit
  velocities, `caimd=.false.`, all-false CAIMD, active zero-strain CAIMD,
  active one-event CAIMD, no-card `vc-md`, explicit-zero-card `vc-md`, and an
  independent Ice X velocity-card reconstruction. After the source change, clean
  current-only K_velocity and full harness runs passed with failures `[]`; clean
  current plus untouched QE 7.5 baseline K_velocity and full harness runs also
  passed with failures `[]`.
- Evidence paths:
  `/private/tmp/qe_vcmd_velocity_diag/current-pre/analysis/summary.json`,
  `/private/tmp/qe_velocity_nstep1`, and
  `/private/tmp/iceX_velocities_300K_regen.{inc,json}`,
  `/private/tmp/qe_velocity_harness_post/k-velocity-3/analysis/summary.json`,
  `/private/tmp/qe_velocity_harness_post/full-current/analysis/summary.json`,
  `/private/tmp/qe_velocity_clean_validation/k-velocity-with-baseline/analysis/summary.json`,
  and `/private/tmp/qe_velocity_clean_validation/full-with-baseline/analysis/summary.json`.
- Key evidence: untouched QE 7.5 baseline `vc-md` with a nonzero explicit
  velocity card saves zero first-step PR `ratd` and zero temperature, while
  baseline fixed-cell MD consumes the same card with input/file velocity mismatch
  `2.99781900692242e-17`. The repaired current `vc-md` explicit case initializes
  PR `ratd_max_abs=2.9673590504451037e-05` with input-to-initialized mismatch
  `2.710505431213761e-20`; no-card and explicit-zero-card `vc-md` controls
  remain zero. The Ice X card independently reconstructs to
  `299.999999999994 K` with `3N-3` degrees of freedom and
  `295.862068965511 K` under the PR `3*(nat+1)` temperature denominator used by
  `.tv`/`.caimd`.
- Files modified: `PW/src/vcsmd.f90`, `PW/src/vcsubs.f90`,
  `test-suite/caimd/caimd_validation.py`, and this report.
- Latest successful command:
  `python3 test-suite/caimd/caimd_validation.py --current-exe /private/tmp/qe_velocity_current_build/bin/pw.x --baseline-exe /private/tmp/qe_velocity_baseline_build/bin/pw.x --pseudo pseudo/C.UPF --work-root /private/tmp/qe_velocity_clean_validation --run-id full-with-baseline --timeout 180`
- Latest failed command: an initial parallel `shasum` raced ahead of
  `/private/tmp/iceX_velocities_300K_regen.{inc,json}` creation; rerun
  succeeded and both regenerated files matched the tracked files.
- Active or latest Slurm job: none for this checkpoint.
- Exact next action: commit the local source/test/doc checkpoint, deploy the
  exact source state to Leonardo, rebuild GPU/MPI `pw.x`, update Slurm checksum
  guards, and run short remote validation before any final deformation job.

## Remote Deployment And Provenance

- Remote host: `login02.leonardo.local`
- Remote user: `mcheng00`
- Remote root:
  `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919`
- Remote example:
  `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/example_iceX_100GPa_300K_caimd`
- `pw.x`:
  `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/install-current-gpu/bin/pw.x`
- `pw.x` SHA256:
  `c9588df640013f098a33e33e1eec6e1f0820ed53f7488efffb0a83ec7ed25b48`
- `cp.x` SHA256:
  `2723af23df8306b5c4a52e380bd6a618606077198592fa1b6d03f8e287ba5882`
- `ldd` inspection found no reported `not found` entries in the inspected
  output.
- Pseudopotential symlinks:
  - `pseudo/H.upf -> /leonardo/home/userexternal/mcheng00/psd/H.upf`
  - `pseudo/O.upf -> /leonardo/home/userexternal/mcheng00/psd/O.upf`

## Completed Validation

- Local static validation: PASS.
- Remote static validation: PASS before short jobs; final archive validation is
  performed after the final local commit.
- Bounded cutoff/k-point/reference jobs: completed.
- Bounded zero-strain and dt/wmass pilots: completed.
- Bounded active tension/compression/shear pilots: completed.
- `.caimd` analysis: 4 active pilot files parsed, 20 rows each, no parser flags.
- Final deformation jobs: not submitted.

## Key Numerical Findings

- Reference SCF pressure: `99.623 GPa`, target deviation `-0.377 GPa`.
- Reference max atomic force: `0.003139 eV/A`.
- Optional `vc-relax`: prepared but not run or used.
- Cutoff: 65 Ry is the highest tested pilot setting; cutoff convergence remains
  unresolved because 55->65 Ry still changes energy by `4.711 meV/atom`.
- K-points: 2x2x2 differs from 3x3x3 by `0.040 meV/atom` and `0.005 GPa`.
- Selected pilot settings: `ecutwfc=65 Ry`, `ecutrho=520 Ry`, `2x2x2`,
  `dt=5.0`, `wmass=65.7122994644 amu`, `caimd_stride=5`,
  `caimd_step_strain=+/-5.0d-5`.

## Leonardo Jobs

| Job ID | Purpose | State | Node | Elapsed |
| --- | --- | --- | --- | --- |
| 49700757 | first cutoff attempt, exposed matrix stdin bug | COMPLETED 0:0 | lrdn3448 | 00:00:36 |
| 49700800 | reference SCF | COMPLETED 0:0 | lrdn0026 | 00:00:58 |
| 49700996 | remaining convergence matrix | COMPLETED 0:0 | lrdn0452 | 00:04:31 |
| 49701674 | zero-strain plus dt/wmass matrix | COMPLETED 0:0 | lrdn0026 | 00:26:40 |
| 49701713 | active strain pilot matrix | COMPLETED 0:0 | lrdn0062 | 00:19:24 |

Peak GPU monitor values:

- reference SCF: 52%, 6084 MiB, 244.83 W.
- convergence matrix: 75%, 6098 MiB, 263.12 W.
- parameter matrix: 82%, 10390 MiB, 310.64 W.
- strain matrix: 77%, 10390 MiB, 320.08 W.

## Final Manual Templates

Prepared but not submitted:

- `deformation/tension`
- `deformation/compression`
- `deformation/shear`

Each contains `qe.in`, `run.slurm`, `README.md`, `expected_files.txt`,
`analyze.sh`, and `CASE_METADATA.json`. The final Slurm scripts use the custom
CAIMD QE 7.5 executable and validate the executable and pseudopotential symlinks
before launch.

Manual commands:

```bash
ROOT=/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919

cd "$ROOT/example_iceX_100GPa_300K_caimd/deformation/tension"
sbatch run.slurm

cd "$ROOT/example_iceX_100GPa_300K_caimd/deformation/compression"
sbatch run.slurm

cd "$ROOT/example_iceX_100GPa_300K_caimd/deformation/shear"
sbatch run.slurm
```

## Evidence Paths

- `convergence/analysis/convergence_summary.json`
- `convergence/analysis/convergence_summary.md`
- `convergence/analysis/stress_strain.csv`
- `convergence/analysis/stress_strain_summary.json`
- `convergence/analysis/caimd_stress_strain_summary.md`
- `evidence/matrix-49701674.tsv`
- `evidence/matrix-49701713.tsv`
- `evidence/LEONARDO_JOB_SUMMARY.md`

## Next Action

Before submitting final deformation jobs, deploy and validate a rebuilt CAIMD
`pw.x` that includes the PW `vc-md`/PR explicit-velocity initialization fix, and
update each Slurm checksum guard to that executable.
