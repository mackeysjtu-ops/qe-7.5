# Final Validation Report

Generated: 2026-07-23T13:57:16Z

Status: `BASELINE_NOSE_FAIL`

## Outcome

Control A already exhibits physical instability in ordinary PW vc-md + Nose; this does not isolate a CAIMD-specific bug.

## Repository And Local Inputs

- Repository root: `/Users/mengqicheng/Downloads/Learning/qe-7.5`
- Branch/HEAD: `main` / `af7ebfc`
- Executable: `/private/tmp/qe_diamond8_af7ebfc_build/bin/pw.x`
- Executable SHA256: `7f7f5f815f53983c3b41dfab1570d96e8f8816d56cdcc574a3a2e33fcf492303`
- Pseudopotential: `/Users/mengqicheng/Downloads/Learning/qe-7.5/pseudo/C.UPF`
- Pseudopotential SHA256: `9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639`
- QE source modified: no
- Git commit created: no
- Leonardo/SSH/Slurm/network/git push: not used

## Prior Evidence Rechecked

- Four-DOF final report exists: False
- 0.30 continuation final report exists: False
- Common terminal velocity temperature from raw `.md`: 559.089239666311 K
- Common `.eal` final-500 mean: 296.671459825860 K
- Previous Test1 max CAIMD temperature: 8503.57805468553 K
- Previous Test2 max CAIMD temperature: 9657.714403976866 K

## Common State

- Source step: 2000
- Original velocity temperature: 559.089239666311 K
- Removed COM velocity au: [-9.970255415540918e-17, 3.7882150961585613e-16, -9.652109840552203e-17]
- Velocity scaling factor: 7.3252096767418784e-01
- Final scaled temperature: 300.000000000000 K
- Ionic DOF: 21
- Minimum periodic C-C distance: 2.683309253044 bohr

## Parameters

- `calculation='vc-md'`, `cell_dynamics='pr'`, `cell_dofree='all'`, `press=0.0`
- `ion_temperature='nose'`, `tempw=300`, `fnosep(1)=50`, `nhpcl=4`, `wmass=2.0e6`
- timestep: 2.0 fs; requested steps: 5000; requested time: 10.0 ps
- cutoffs: ecutwfc=25.0 Ry, ecutrho=120.0 Ry; Gamma point

## Control Results

- Control A: FAIL; steps=664; JOB DONE=False
  EAL mean/final500/max K: 3843.933566439214 / 4770.521264699599 / 8047.601691661173
  Terminal direct `.md` T K: 5064.727517264306
  CAIMD T mean K: None
- Control B: not run or no validation available.
- Control C: not run or no validation available.

## No-Op And Projection Diagnostics

- A-vs-B maximum numeric difference: None (None)
- A-vs-B no-op classified OK: False
- Control C maximum cell-strain residual: None
- Control C velocity jump max: None
- Earliest diverging observable: `{'control': 'A', 'observable': 'eal ionic temperature', 'step': 26}`

## Temperature Alignment

- `.caimd` rows are pre-move diagnostics written before `move_ions`.
- `.eal` rows are written in `vcsmd` after the PR/Beeman update with `nst=istep-1`.
- `analysis/temperature_crosscheck.csv` records zero/one-step offset checks where `.caimd` exists.

## Required Questions

1. Ordinary PW `vc-md` + Nose stable? not stable/classification FAIL
2. All-false CAIMD valid no-op under Nose? not determined because Control A failed first
3. Active all-true zero-strain CAIMD stable? not run
4. High temperature physically present in velocities? See terminal `.md` direct temperatures above and `temperature_crosscheck.csv`.
5. First divergent control: {'control': 'A', 'observable': 'eal ionic temperature', 'step': 26}
6. Earliest differing physical variable: {'control': 'A', 'observable': 'eal ionic temperature', 'step': 26}
7. Likely issue class: BASELINE_NOSE_FAIL

## Local Checks

- Python syntax check: 0
- Shell syntax check: 0
- `git diff --check`: 0
- Tracked source diff empty: True
- Targeted local CAIMD+Nose harness return code: 0
- Full local CAIMD harness: not_run

## Output Files

- `analysis/temperature_statistics.csv`
- `analysis/temperature_crosscheck.csv`
- `analysis/cell_statistics.csv`
- `analysis/nose_extended_energy_statistics.csv`
- `analysis/control_A_vs_B_differences.csv`
- `analysis/control_C_projection_statistics.csv`
- `analysis/DIAGNOSTIC_CLASSIFICATION.json`
