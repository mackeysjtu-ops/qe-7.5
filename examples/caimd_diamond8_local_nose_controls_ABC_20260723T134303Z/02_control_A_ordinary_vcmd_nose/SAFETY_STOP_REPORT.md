# Safety Stop Report

Generated: 2026-07-23T13:57:15Z

Status: FAIL

Control: A

Run directory: `/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose`

Prefix: `diamond8_control_A_ordinary_vcmd_nose`

Completed trajectory blocks: 664

JOB DONE: False

Primary temperature source: .eal ionic kinetic energy

EAL temperature stats:

```json
{
  "count": 664,
  "finite": true,
  "mean": 3843.933566439214,
  "std": 2000.0237223448755,
  "min": 184.5629835450013,
  "max": 8047.601691661173,
  "first_500_mean": 3558.3748340024267,
  "final_500_mean": 4770.521264699599,
  "final": 5064.730659738254,
  "drift_per_ps": 3351.0491708631257,
  "rolling25_min": 242.28214700740045,
  "rolling25_max": 5412.335807042567,
  "rolling100_min": 242.28214700740045,
  "rolling100_max": 4940.54386093579,
  "block500_min_mean": 3558.3748340024267,
  "block500_max_mean": 4714.539458014783,
  "fraction_below_150K": 0.0,
  "fraction_above_600K": 0.8192771084337349,
  "longest_below_150K_steps": 0,
  "longest_above_600K_steps": 541,
  "longest_above_2000K_steps": 189,
  "first_step_above_600K": 26,
  "first_step_above_2000K": 130
}
```

Safety context:

```json
{
  "control": "A",
  "elapsed_s": 55.119303941726685,
  "preserved_outputs": [
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/qe.out",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/qe.err",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/diamond8_control_A_ordinary_vcmd_nose.e",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/diamond8_control_A_ordinary_vcmd_nose.eal",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/diamond8_control_A_ordinary_vcmd_nose.ave",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/diamond8_control_A_ordinary_vcmd_nose.avec",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/diamond8_control_A_ordinary_vcmd_nose.p",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/diamond8_control_A_ordinary_vcmd_nose.tv",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/tmp/diamond8_control_A_ordinary_vcmd_nose.md",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/tmp/diamond8_control_A_ordinary_vcmd_nose.cel",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/tmp/diamond8_control_A_ordinary_vcmd_nose.pos",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/tmp/diamond8_control_A_ordinary_vcmd_nose.vel",
    "/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z/02_control_A_ordinary_vcmd_nose/tmp/diamond8_control_A_ordinary_vcmd_nose.for"
  ],
  "returncode": -15,
  "safety": {
    "consecutive_records": 20,
    "control": "A",
    "detected_utc": "2026-07-23T13:54:45Z",
    "first_threshold_step": 130,
    "first_threshold_temperature_K": 2533.4179132851405,
    "latest_step": 664,
    "latest_temperature_K": 5064.730659738254,
    "prefix": "diamond8_control_A_ordinary_vcmd_nose",
    "reason": "temperature above 2000 K for sustained records"
  },
  "status": "SAFETY_STOP"
}
```

QE error lines:

```json
[]
```

The partial outputs are preserved in this directory and are not treated as a
successful 5000-step control.
