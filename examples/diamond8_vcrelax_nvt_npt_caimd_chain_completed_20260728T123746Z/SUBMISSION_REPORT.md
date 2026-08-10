# Diamond8 vc-relax / NVT / NPT / CAIMD chain — submission report

Submitted 2026-07-28T12:08:47Z. Remote root: /leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_vcrelax_nvt_npt_caimd_chain_20260728T120614Z

## Jobs (strict afterok chain)
| # | Job ID | Name | Script | Dependency |
|---|---|---|---|---|
| 1 | 50558068 | d8_vcrelax | scripts/01_vc_relax.slurm | (none) |
| 2 | 50558069 | d8_fixed_nvt | scripts/02_fixed_nvt.slurm | afterok:50558068 |
| 3 | 50558070 | d8_free_npt | scripts/03_free_npt.slurm | afterok:50558069 |
| 4 | 50558071 | d8_caimd_zn | scripts/04_caimd_zero_nonzero.slurm | afterok:50558070 |

## Contract
- 8 carbon atoms only (no H/O); ecutwfc 80 Ry, ecutrho 320 Ry; K_POINTS automatic 2 2 2 0 0 0; NC LDA C.UPF.
- Executable SHA256 89d4fda8...10e7d; pseudo SHA256 9c06aa64...c978639 (both verified).
- JOB1 vc-relax (ibrav=1, cell_dofree=volume, press=0). JOB2 fixed-cell md + ionic Nose 300 K.
  JOB3 vc-md caimd=.false. + ionic Nose 300 K + cell Nose (cell_temperature=nose, temph=300, fnoseh=1) + cell_dofree=volume, press=0.
  JOB4 (ibrav=0, cell_dofree=all): Stage 4A active CAIMD all-6-masks zero-strain (gated); Stage 4B xx-only increment (step_strain(1)=5.0d-6, stride=10).
- No ATOMIC_VELOCITIES / ion_velocities anywhere; native velocity initialization.

## Future inspection
sacct -j 50558068,50558069,50558070,50558071 --format=JobID,JobName%22,State,ExitCode,Elapsed,NodeList%16,Reason
