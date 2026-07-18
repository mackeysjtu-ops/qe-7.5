# Tensor xy shear

Target state:

- Pressure: 100 GPa, written as `press = 1000.0` because QE PW input uses kbar.
- Nominal initial temperature: deterministic 300 K velocity card from `../../structure/iceX_velocities_300K.inc`.
- Active thermostats: none. Active CAIMD uses `ion_temperature='not_controlled'` and `cell_temperature='not_controlled'`.
- Current CAIMD V1 caveat: QE 7.5 `vc-md`/PR with `not_controlled` initializes the variable-cell MD temperature from the namelist path, not from the explicit `ATOMIC_VELOCITIES` card used here. The velocity card is retained for reproducibility and future source work, but the current Leonardo pilot evidence starts near 0 K.

Structure and provenance:

- Source: `../../structure/iceX_original.data`.
- Source SHA256: `103ccf1fe76e3b3a737914a4eeddb4016d6fb8b308129516541ddc5de08aa6c6`.
- Converted structure: `../../structure/iceX_qe_structure.inc`, no rotation, relaxation, symmetrization, or atom-order changes.
- Pseudopotentials are expected as remote symlinks at `../../pseudo/H.upf` and `../../pseudo/O.upf`.

Parameters:

- `ecutwfc = 65.0 Ry`
- `ecutrho = 520.0 Ry`
- `K_POINTS automatic 2 2 2 0 0 0`
- `dt = 5.0000` Rydberg atomic time units
- `wmass = 65.7122994644 amu`
- `caimd_mask(xx, yy, zz, yz, xz, xy) = (F, F, F, F, F, T)`
- `caimd_step_strain(xx, yy, zz, yz, xz, xy) = (0.0e+00, 0.0e+00, 0.0e+00, 0.0e+00, 0.0e+00, 5.0e-05)`
- `caimd_stride = 5`
- Per-event masked nominal increment: `0.00025`
- `nstep = 100`
- Expected events: `20`
- Final nominal strain: +0.005 nominal xy tensor strain, equivalent to +0.010 engineering shear
- Loading direction/convention: positive xy tensor shear; engineering gamma_xy = 2 * tensor_E_xy

Expected outputs:

- `slurm-<jobid>.out`
- `slurm-<jobid>.err`
- `nvidia-smi-<jobid>.log`
- `qe.out`
- `tmp/icex_100gpa_300k_shear.caimd`
- `icex_100gpa_300k_shear.caimd`
- `../../convergence/analysis/stress_strain.csv`
- `../../convergence/analysis/stress_strain_summary.json`
- common analysis outputs under `../../convergence/analysis/` after `./analyze.sh`

This is a bounded pilot template, not a production-converged mechanical protocol.
Review `../../PARAMETER_SELECTION.md` and `../../ICE_X_CAIMD_EXAMPLE_REPORT.md`
before submitting.
