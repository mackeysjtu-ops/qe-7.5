# Final Validation Report

Status: BLOCKED

Workflow directory: `/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_four_dof_tests_0p03_20260723T115714Z`

Branch/HEAD: `main` / `af7ebfc Support Nose thermostat with active CAIMD`

Executable: `/private/tmp/qe_diamond8_af7ebfc_build/bin/pw.x`

Executable SHA256: `7f7f5f815f53983c3b41dfab1570d96e8f8816d56cdcc574a3a2e33fcf492303`

Pseudopotential: `/Users/mengqicheng/Downloads/Learning/qe-7.5/pseudo/C.UPF`

Pseudopotential SHA256: `9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639`

No Leonardo, Slurm, SSH, network access, or Git push was used by this workflow.

No QE source code was modified by this workflow.

## Common Setup

- timestep: `2.0 fs`
- test duration: `10.0 ps`
- strain rate for Tests 2-4: `0.003 ps^-1`
- target final nominal xx strain: `-0.03`
- Nose parameters: `fnosep=50.0 THz`, `nhpcl=4`
- `wmass=2000000.0`
- common-state source step: `2000`
- common-state velocity temperature: `559.0892396663107` K

## Common Equilibration

```json
{
  "block500_max_mean": 296.67145982585964,
  "block500_min_mean": 177.6277339630961,
  "count": 2000,
  "drift_K_per_ps": 38.22305303269425,
  "final_500_mean": 296.67145982585964,
  "finite": true,
  "first_500_mean": 177.6277339630961,
  "fraction_above_600K": 0.0,
  "fraction_below_150K": 0.163,
  "longest_above_2000K_steps": 0,
  "longest_above_600K_steps": 0,
  "longest_below_150K_steps": 8,
  "max": 559.1473474923312,
  "mean": 215.59478242367328,
  "min": 34.07813993873525,
  "rolling100_max": 385.10283002355357,
  "rolling100_min": 157.05326486030933,
  "rolling25_max": 410.11657269075823,
  "rolling25_min": 152.94205448023183,
  "std": 74.42361401504472
}
```

## Test Comparison

| test | completed | event | final nominal xx | final actual xx | T mean | T std | final500 T | det min | cond max | velocity jump | class |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| test1 | 3576 | 358 | 0 | 0 | 4611.46 | 1178.8 | 4448.57 | 257.815 | 1.1583 | 2.62009e-18 | FAIL |
| test2 | 779 | 78 | -0.00468 | -0.00468 | 4954.88 | 1532.72 | 5245.41 | 256.63 | 1.00593 | 1.73472e-18 | FAIL |
| test3 | None | None |  |  |  |  |  |  |  |  | None |
| test4 | None | None |  |  |  |  |  |  |  |  | None |

## Old Failure Comparison

The previous 0.05 ps^-1 local stress test reached only about -0.0265 nominal xx before QE stopped with `caimd projected cell is singular`.

```json
{
  "caimd_path": "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_continuation_0p30_20260723T102331Z/02_continuous_0p00_to_0p30/tmp/diamond8_continuous.caimd",
  "first_abs_actual_component_ge_0.05": [
    {
      "component": "zz",
      "ionic_step": 185,
      "nominal_xx": -0.018500000000000013,
      "threshold": 0.05,
      "time_ps": 0.3680000682817448,
      "value": 0.0518628406245667
    },
    {
      "component": "xz",
      "ionic_step": 198,
      "nominal_xx": -0.020000000000000014,
      "threshold": 0.05,
      "time_ps": 0.3940000731059985,
      "value": -0.051037993383567175
    },
    {
      "component": "yy",
      "ionic_step": 219,
      "nominal_xx": -0.022000000000000016,
      "threshold": 0.05,
      "time_ps": 0.43600008089902376,
      "value": 0.05301615416793548
    },
    {
      "component": "yz",
      "ionic_step": 229,
      "nominal_xx": -0.023000000000000017,
      "threshold": 0.05,
      "time_ps": 0.45600008460998814,
      "value": -0.05233899919464141
    },
    {
      "component": "xy",
      "ionic_step": 237,
      "nominal_xx": -0.024000000000000018,
      "threshold": 0.05,
      "time_ps": 0.4720000875787596,
      "value": 0.051279070828451054
    }
  ],
  "first_abs_actual_component_ge_0.1": [
    {
      "component": "zz",
      "ionic_step": 191,
      "nominal_xx": -0.019000000000000013,
      "threshold": 0.1,
      "time_ps": 0.38000007050832346,
      "value": 0.10245880108284555
    },
    {
      "component": "xz",
      "ionic_step": 210,
      "nominal_xx": -0.021000000000000015,
      "threshold": 0.1,
      "time_ps": 0.4180000775591558,
      "value": -0.10166018577219099
    },
    {
      "component": "yy",
      "ionic_step": 232,
      "nominal_xx": -0.023500000000000017,
      "threshold": 0.1,
      "time_ps": 0.46200008572327744,
      "value": 0.10408518108487064
    },
    {
      "component": "yz",
      "ionic_step": 244,
      "nominal_xx": -0.02450000000000002,
      "threshold": 0.1,
      "time_ps": 0.4860000901764347,
      "value": -0.10356209883174816
    },
    {
      "component": "xy",
      "ionic_step": 264,
      "nominal_xx": -0.02650000000000002,
      "threshold": 0.1,
      "time_ps": 0.5260000975983635,
      "value": -0.19549750539721616
    }
  ],
  "first_abs_actual_component_ge_0.25": [
    {
      "component": "zz",
      "ionic_step": 204,
      "nominal_xx": -0.020500000000000015,
      "threshold": 0.25,
      "time_ps": 0.40600007533257715,
      "value": 0.2527554069216642
    },
    {
      "component": "yy",
      "ionic_step": 255,
      "nominal_xx": -0.02550000000000002,
      "threshold": 0.25,
      "time_ps": 0.5080000942584956,
      "value": 0.25388020435115544
    },
    {
      "component": "yz",
      "ionic_step": 258,
      "nominal_xx": -0.02600000000000002,
      "threshold": 0.25,
      "time_ps": 0.5140000953717848,
      "value": -0.2590158898491714
    },
    {
      "component": "xz",
      "ionic_step": 260,
      "nominal_xx": -0.02600000000000002,
      "threshold": 0.25,
      "time_ps": 0.5180000961139777,
      "value": 0.27631015095899514
    },
    {
      "component": "xy",
      "ionic_step": 265,
      "nominal_xx": -0.02650000000000002,
      "threshold": 0.25,
      "time_ps": 0.5280000979694599,
      "value": -0.46212053710344453
    }
  ],
  "first_abs_actual_component_ge_0.5": [
    {
      "component": "zz",
      "ionic_step": 219,
      "nominal_xx": -0.022000000000000016,
      "threshold": 0.5,
      "time_ps": 0.43600008089902376,
      "value": 0.5040686142741149
    },
    {
      "component": "xz",
      "ionic_step": 262,
      "nominal_xx": -0.02650000000000002,
      "threshold": 0.5,
      "time_ps": 0.5220000968561705,
      "value": 0.5966940686762194
    },
    {
      "component": "yz",
      "ionic_step": 264,
      "nominal_xx": -0.02650000000000002,
      "threshold": 0.5,
      "time_ps": 0.5260000975983635,
      "value": -0.5143531519603401
    },
    {
      "component": "yy",
      "ionic_step": 265,
      "nominal_xx": -0.02650000000000002,
      "threshold": 0.5,
      "time_ps": 0.5280000979694599,
      "value": 0.5581433215010438
    }
  ],
  "first_temperature_ge_2000K": {
    "actual_components": {
      "xx": -0.000500000000000056,
      "xy": 7.33487167305583e-08,
      "xz": 6.685559823204113e-08,
      "yy": 1.000460805311576e-07,
      "yz": -9.872664778175537e-08,
      "zz": 4.427092630265861e-07
    },
    "ionic_step": 3,
    "nominal_xx": -0.0005,
    "temperature_K": 90977.49352344994,
    "time_ps": 0.004000000742192878
  },
  "first_temperature_ge_600K": {
    "actual_components": {
      "xx": -0.000500000000000056,
      "xy": 7.33487167305583e-08,
      "xz": 6.685559823204113e-08,
      "yy": 1.000460805311576e-07,
      "yz": -9.872664778175537e-08,
      "zz": 4.427092630265861e-07
    },
    "ionic_step": 3,
    "nominal_xx": -0.0005,
    "temperature_K": 90977.49352344994,
    "time_ps": 0.004000000742192878
  },
  "largest_stepwise_unmasked_growth": [
    {
      "at_ionic_step": 265,
      "component": "zz",
      "max_step_to_step_delta": 0.45231870294152166
    },
    {
      "at_ionic_step": 264,
      "component": "xz",
      "max_step_to_step_delta": 0.41570581634677173
    },
    {
      "at_ionic_step": 265,
      "component": "xy",
      "max_step_to_step_delta": 0.26662303170622836
    },
    {
      "at_ionic_step": 265,
      "component": "yy",
      "max_step_to_step_delta": 0.09512512846089716
    },
    {
      "at_ionic_step": 264,
      "component": "yz",
      "max_step_to_step_delta": 0.05609443338061082
    }
  ],
  "last_row": {
    "actual_components": {
      "xx": -0.026499999999999968,
      "xy": -0.46212053710344453,
      "xz": 1.0796891802659618,
      "yy": 0.5581433215010438,
      "yz": -0.4741702477636549,
      "zz": 0.8004691398698094
    },
    "event_count": 53,
    "ionic_step": 265,
    "nominal_components": {
      "xx": -0.02650000000000002,
      "xy": 0.0,
      "xz": 0.0,
      "yy": 0.0,
      "yz": 0.0,
      "zz": 0.0
    },
    "temperature_K": 70017846364.13759,
    "time_ps": 0.5280000979694599,
    "velocity_jump_max": 1.640126728885195e-14,
    "volume_bohr3": 237.37055141113925
  },
  "previous_root": "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_continuation_0p30_20260723T102331Z",
  "rate_ps_inv": 0.05,
  "row_count": 265
}
```

## Interpretation

At least one later test failed; its `SAFETY_STOP_REPORT.md` identifies the first unstable observable. Earlier valid tests remain independent evidence because every test starts from the same common state.

This run distinguishes thermostat behavior from free PR-cell behavior by comparing fully constrained loading, transverse-normal-only relaxation, and full unmasked PR response at the same 0.003 ps^-1 rate.
