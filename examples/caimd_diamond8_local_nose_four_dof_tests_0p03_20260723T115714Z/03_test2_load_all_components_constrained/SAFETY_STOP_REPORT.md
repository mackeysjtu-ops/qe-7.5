# Safety Stop Report

Generated: 2026-07-23T12:36:03Z

Status: FAIL

Run directory: `/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_four_dof_tests_0p03_20260723T115714Z/03_test2_load_all_components_constrained`

Prefix: `diamond8_test2_all_components_constrained`

Expected steps: `5000`

Completed trajectory blocks: `{'cel': 779, 'pos': 779, 'vel': 779}`

Return code: `None`

JOB DONE: `False`

First failing observable:

```json
{
  "observable": "temperature_K",
  "reason": "sustained above 600 K",
  "step_index": 26,
  "value": 603.2009759377513
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
  "count": 779,
  "final_500_mean": 257.0101515768595,
  "finite": true,
  "first_500_mean": 257.44176949817535,
  "max": 257.82078881694946,
  "mean": 257.22595159932246,
  "min": 256.6295854657888,
  "std": 0.34787682826051464
}
```

Cell condition stats:

```json
{
  "count": 779,
  "final_500_mean": 1.0047143644382828,
  "finite": true,
  "first_500_mean": 1.0036350697026046,
  "max": 1.0059271412368667,
  "mean": 1.0042062765134734,
  "min": 1.003137206809314,
  "std": 0.0008747801018310189
}
```

The partial outputs are preserved in this directory. This test must not be
treated as successful unless it is rerun and validated from a clean input state.
