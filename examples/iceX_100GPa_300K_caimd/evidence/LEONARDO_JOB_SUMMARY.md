# Leonardo Job Summary

Remote example:

`/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/example_iceX_100GPa_300K_caimd`

## Jobs

| Job ID | Purpose | State | Exit | Node | Elapsed | Evidence |
| --- | --- | --- | --- | --- | --- | --- |
| 49700757 | first convergence attempt; exposed matrix stdin consumption bug | COMPLETED | 0:0 | lrdn3448 | 00:00:36 | `evidence/matrix-49700757.tsv` on remote |
| 49700800 | reference SCF | COMPLETED | 0:0 | lrdn0026 | 00:00:58 | `reference/scf/qe.out` on remote |
| 49700996 | remaining cutoff/k-point convergence matrix | COMPLETED | 0:0 | lrdn0452 | 00:04:31 | `convergence/analysis/convergence_summary.json` |
| 49701674 | zero-strain plus dt/wmass pilot matrix | COMPLETED | 0:0 | lrdn0026 | 00:26:40 | `evidence/matrix-49701674.tsv` |
| 49701713 | active tension/compression/shear pilot matrix | COMPLETED | 0:0 | lrdn0062 | 00:19:24 | `evidence/matrix-49701713.tsv` |

## GPU Monitor Peaks

| Log | Samples | Peak GPU % | Peak memory MiB | Peak power W |
| --- | ---: | ---: | ---: | ---: |
| `evidence/nvidia-smi-matrix-49700757.log` | 1 | 0 | 0 | 59.64 |
| `reference/scf/nvidia-smi-49700800.log` | 2 | 52 | 6084 | 244.83 |
| `evidence/nvidia-smi-matrix-49700996.log` | 9 | 75 | 6098 | 263.12 |
| `evidence/nvidia-smi-matrix-49701674.log` | 54 | 82 | 10390 | 310.64 |
| `evidence/nvidia-smi-matrix-49701713.log` | 39 | 77 | 10390 | 320.08 |

Final tension, compression, and shear jobs were not submitted.
