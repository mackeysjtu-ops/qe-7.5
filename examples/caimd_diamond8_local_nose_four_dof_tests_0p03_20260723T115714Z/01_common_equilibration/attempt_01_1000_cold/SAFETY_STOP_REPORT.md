# Safety Stop Report

Generated: 2026-07-23T12:09:08Z

Status: FAIL

Run directory: `/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_four_dof_tests_0p03_20260723T115714Z/01_common_equilibration`

Prefix: `diamond8_common_equil`

Expected steps: `1000`

Completed trajectory blocks: `{'cel': 1000, 'pos': 1000, 'vel': 1000}`

Return code: `0`

JOB DONE: `True`

First failing observable:

```json
null
```

QE error lines:

```json
[]
```

Temperature classification: `FAIL`

Cell determinant stats:

```json
{
  "count": 1000,
  "final_500_mean": 297.917964015947,
  "finite": true,
  "first_500_mean": 305.04602993622575,
  "max": 306.2710188038857,
  "mean": 301.4819969760864,
  "min": 292.2090962869458,
  "std": 4.227798679810431
}
```

Cell condition stats:

```json
{
  "count": 1000,
  "final_500_mean": 1.0006794772141323,
  "finite": true,
  "first_500_mean": 1.0000957851092593,
  "max": 1.0011001247574336,
  "mean": 1.0003876311616957,
  "min": 1.0000000015367205,
  "std": 0.000342784335582633
}
```

The partial outputs are preserved in this directory. This test must not be
treated as successful unless it is rerun and validated from a clean input state.
