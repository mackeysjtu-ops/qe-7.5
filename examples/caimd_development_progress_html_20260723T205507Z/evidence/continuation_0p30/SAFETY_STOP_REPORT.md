# Safety Stop Report

Generated: 2026-07-23T11:21:26Z local macOS run.

## Status

The workflow stopped during the active continuous CAIMD loading stage. Segmented
continuation runs and final analysis plots/CSVs were not executed because the
continuous trajectory became numerically and physically invalid.

## Last Valid Prerequisite

The ordinary Nose+Parrinello-Rahman equilibration passed with:

- `nstep=500`.
- `dt=41.341381006`, observed as `0.002 ps` per step.
- `ecutwfc/ecutrho=25/120 Ry`.
- `ion_temperature='nose'`, `tempw=300 K`, `fnosep=50 THz`, `nhpcl=4`.
- Initial velocity warm start: 450 K under the 3N-3 ionic convention.
- `cell_dynamics='pr'`, `cell_dofree='all'`, `press=0`, `wmass=2.0e6 AMU`.
- Validation source for equilibration temperature: `.eal` ionic kinetic energy.
- Equilibration validation: `valid=true`, temperature classification `PASS`,
  ionic 3N-3 mean `251.672 K`, final 200-step mean `252.092 K`.
- Extracted terminal velocity temperature: `317.176 K` under the 3N-3
  convention.

## Failed Stage

Continuous active CAIMD run:

- Directory: `02_continuous_0p00_to_0p30`.
- Input: `qe.in`.
- Output: `qe.out`.
- QE return code: `1`.
- `JOB DONE`: absent.
- Completed trajectory blocks: 264 cell/position/velocity blocks.
- CAIMD diagnostic rows: 265.
- Last CAIMD event count: 53.
- Last nominal `xx` strain: `-0.0265`.
- Last actual `xx` strain: `-0.0265` within numerical tolerance.
- CAIMD `xx_abs_error_max`: about `2.26e-16`.
- CAIMD velocity residual max: about `1.64e-14`.

The constrained `xx` projection itself stayed accurate, but the unconstrained
cell degrees of freedom became unstable:

- Last reported actual strains from `.caimd` were approximately
  `xx=-0.0265`, `yy=0.5581`, `zz=0.8005`, `yz=-0.4742`,
  `xz=1.0797`, `xy=-0.4621`.
- Last CAIMD temperature: about `7.00e10 K`.
- Last pressure from validation: `11776.88 kbar`.
- Minimum periodic C-C distance reported by the local trajectory validator:
  `0.2563 bohr`.
- QE stopped with:
  `Error in routine caimd_inverse3 (1): caimd projected cell is singular`.

## Interpretation

This is a genuine safety stop for the requested local validation conditions.
The active CAIMD projection maintained the imposed lab `xx` strain accurately,
but the remaining unconstrained PR cell degrees of freedom evolved into a
singular projected cell under the local high-rate loading and large residual
stress. Continuing to segmented validation or producing stress/strain plots
from this trajectory would be misleading.

## Not Executed

- Three 1000-step segmented continuation runs to global `-0.10`, `-0.20`, and
  `-0.30`.
- Boundary validation between segments.
- Continuous-vs-segmented final comparison.
- Final plots and CSV analysis.

No remote, Slurm, Leonardo, network, or push operations were used.
