# Safety Stop Report

Generated: 2026-07-23T12:34:21Z

Status: FAIL

Run directory: `/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_four_dof_tests_0p03_20260723T115714Z/02_test1_zero_load_full_pr`

Prefix: `diamond8_test1_zero_load_full_pr`

Expected steps: `5000`

Completed trajectory blocks: `{'cel': 3576, 'pos': 3576, 'vel': 3576}`

Return code: `None`

JOB DONE: `False`

First failing observable:

```json
{
  "observable": "temperature_K",
  "reason": "sustained above 600 K",
  "step_index": 26,
  "value": 606.8279664883704
}
```

QE error lines:

```json
[]
```

Temperature classification: `FAIL`

Cell determinant stats:

```json
{
  "count": 3576,
  "final_500_mean": 291.3128822205318,
  "finite": true,
  "first_500_mean": 257.8805625613096,
  "max": 298.3777184010486,
  "mean": 269.42298923750747,
  "min": 257.81512173868316,
  "std": 11.749476265696863
}
```

Cell condition stats:

```json
{
  "count": 3576,
  "final_500_mean": 1.1354275883588507,
  "finite": true,
  "first_500_mean": 1.0030654239357482,
  "max": 1.158303384975943,
  "mean": 1.0617766572254488,
  "min": 1.0025927740818137,
  "std": 0.04703654555555345
}
```

The partial outputs are preserved in this directory. This test must not be
treated as successful unless it is rerun and validated from a clean input state.
