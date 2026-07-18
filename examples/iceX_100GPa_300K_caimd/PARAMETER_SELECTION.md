# Parameter Selection

This file records the parameter rationale and bounded Leonardo evidence used to
prepare the Ice X CAIMD example. The settings are pilot/example settings, not
publication-quality production parameters.

## Structure

- Source file: local `iceX.data`.
- Source SHA256:
  `103ccf1fe76e3b3a737914a4eeddb4016d6fb8b308129516541ddc5de08aa6c6`.
- Atom count: 144 total, O=48, H=96.
- Type mapping: LAMMPS type 1=O and type 2=H from Masses-section comments and
  masses.
- Cell: orthorhombic, rows in Angstrom:
  - `9.2145100000 0.0000000000 0.0000000000`
  - `0.0000000000 6.5156430000 0.0000000000`
  - `0.0000000000 0.0000000000 7.5236160000`
- Volume: `451.706300063683 Angstrom^3`.
- Density under Angstrom interpretation: `3.178909003080 g/cm^3`.
- The conversion applies no rotation, symmetrization, relaxation, or atom-order
  changes.

## Pseudopotentials

| Element | Remote path | SHA256 | Type | XC | Recommended cutoff |
| --- | --- | --- | --- | --- | --- |
| H | `/leonardo/home/userexternal/mcheng00/psd/H.upf` | `c7182caf52c4cad9fe912262ea40feaf4c3b757204189a14cb365869e3e0888e` | PAW, `is_ultrasoft="T"` | PBE | not present |
| O | `/leonardo/home/userexternal/mcheng00/psd/O.upf` | `fb96915f25537a73000b2db2fde27610f287b923f62c2ff52a674f651f394203` | PAW, `is_ultrasoft="T"` | PBE | not present |

Both files are UPF v2.0.1 ATOMPAW scalar-relativistic PBE datasets. No
machine-readable recommended cutoff was found, so the cutoff choice is based on
bounded testing. The example uses `ecutrho/ecutwfc = 8` because these UPFs carry
PAW/USPP-like flags.

## Pressure

Target pressure is `100 GPa`. QE PW input documentation states that `press` is
in kbar, and `PW/src/input.f90` converts with `press_ = press / ry_kbar`.
Therefore the inputs use:

```text
press = 1000.0
```

because `100 GPa = 1000 kbar`.

## Cutoff Matrix

Cutoff tests used gamma-only k-points and `ecutrho = 8 * ecutwfc`.

| ecutwfc Ry | ecutrho Ry | Energy Ry | Pressure GPa | Total force Ry/bohr | Max force eV/A | Delta to 65 Ry meV/atom |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 35 | 280 | -2109.44021125 | 95.258 | 0.008658 | 0.1844 | 23.063 |
| 45 | 360 | -2109.57776344 | 100.584 | 0.009089 | 0.1857 | 10.067 |
| 55 | 440 | -2109.63444610 | 99.811 | 0.009149 | 0.1856 | 4.711 |
| 65 | 520 | -2109.68430794 | 100.322 | 0.009073 | 0.1853 | 0.000 |

Status: unresolved for strict production convergence. The bounded matrix did not
bracket the default 1 meV/atom criterion by 65 Ry. The final examples retain
65/520 Ry as the highest tested pilot setting, not as a production-converged
cutoff.

## K-Point Matrix

K-point tests used `ecutwfc=65 Ry`, `ecutrho=520 Ry`.

| Grid | Energy Ry | Pressure GPa | Total force Ry/bohr | Max force eV/A | Delta to 3x3x3 meV/atom |
| --- | ---: | ---: | ---: | ---: | ---: |
| gamma | -2109.68430794 | 100.322 | 0.009073 | 0.1853 | 10.380 |
| 2x2x2 | -2109.79374555 | 99.623 | 0.000277 | 0.0031 | 0.040 |
| 3x3x3 | -2109.79416654 | 99.618 | 0.000031 | 0.0015 | 0.000 |

Selected bounded-pilot grid: `2 2 2 0 0 0`. It is close to 3x3x3 in total
energy and pressure for this supercell at the tested cutoff. A conservative
production candidate would rerun sensitivity checks with 3x3x3 and an expanded
cutoff range.

## Reference State

Reference SCF uses the original converted structure, `ecutwfc=65 Ry`,
`ecutrho=520 Ry`, and 2x2x2 k-points.

| Quantity | Value |
| --- | ---: |
| Total energy | -2109.79374555 Ry |
| Pressure | 99.623 GPa |
| Deviation from target | -0.377 GPa |
| Total force | 0.000277 Ry/bohr |
| Max atomic force | 0.003139 eV/A |
| Wall time | 48.80 s |
| Peak GPU activity | 52%, 6084 MiB, 244.83 W |

The pressure is within 5 GPa of 100 GPa and max force is below 0.05 eV/A, so the
optional 0 K `vc-relax` template was prepared but not run and not used.

## Time Step And Cell Mass

PW `dt` is in Rydberg atomic time units. QE prints femtoseconds as
`dt_fs = dt * 2.0e3 * au_ps`, so `dt=5.0` is approximately `0.241889 fs`.

For PR `vc-md`, QE derives the default `wmass` as:

```text
wmass_default = 0.75 * total_atomic_mass / pi^2
```

For this 144-atom cell, `wmass_default = 65.7122994644 amu`.

The parameter screen completed six short cases. These historical Leonardo pilot
runs used the earlier deployed CAIMD executable whose PW `vc-md`/PR
initialization path did not honor the explicit 300 K velocity card. They are
therefore near-0 K software stability screens, not 300 K MD sensitivity tests.

| Case | dt | wmass amu | Status | Final pressure GPa | Max force eV/A | T range K | Wall time |
| --- | ---: | ---: | --- | ---: | ---: | --- | --- |
| zero strain | 5.0 | 65.7122994644 | ok | 99.677 | 0.0068 | 0.0-0.0 | 4m22s |
| dt5 default | 5.0 | 65.7122994644 | ok | 99.677 | 0.0068 | 0.0-0.0 | 4m22s |
| dt7.5 default | 7.5 | 65.7122994644 | ok | 99.731 | 0.0103 | 0.0-0.0 | 4m23s |
| dt10 default | 10.0 | 65.7122994644 | ok | 99.801 | 0.0139 | 0.0-0.0 | 4m25s |
| dt5 half wmass | 5.0 | 32.8561497322 | ok | 99.720 | 0.0084 | 0.0-0.0 | 4m22s |
| dt5 double wmass | 5.0 | 131.4245989289 | ok | 99.654 | 0.0057 | 0.0-0.0 | 4m22s |

Selected bounded-pilot values:

- `dt = 5.0`
- `wmass = 65.7122994644 amu`

These are conservative short-run software settings. The explicit-velocity
initialization issue now has local source-level validation, Leonardo GPU/MPI
harness validation, and a short Ice X nstep=1 smoke showing `295.86 K`
initialization. The historical pilot matrix itself remains near-0 K evidence and
should not be reinterpreted as 300 K MD sensitivity data.

## CAIMD Strain Schedule

CAIMD component order is:

```text
(xx, yy, zz, yz, xz, xy)
```

The pilot strain cases used `nstep=20`, `caimd_stride=5`, and 4 imposed events.
No `.caimd` parser flags were raised.

| Case | Step strain | Mask | Final nominal strain | Final actual strain | Final stress component GPa | Max T K |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| tension xx small | +1.0e-5 | xx | +0.0002 xx | +0.0002 xx | stress_xx 99.437 | 0.00381 |
| tension xx selected | +5.0e-5 | xx | +0.0010 xx | +0.0010 xx | stress_xx 98.522 | 0.00737 |
| compression xx selected | -5.0e-5 | xx | -0.0010 xx | -0.0010 xx | stress_xx 100.821 | 0.00354 |
| shear xy selected | +5.0e-5 | xy | +0.0010 tensor xy | +0.0010 tensor xy | stress_xy -0.395 | 0.00932 |

Selected final-template strain settings:

- `caimd_stride = 5`
- tension: `caimd_step_strain(1) = +5.0d-5`
- compression: `caimd_step_strain(1) = -5.0d-5`
- shear: `caimd_step_strain(6) = +5.0d-5`
- event increment: `caimd_stride * caimd_step_strain = +/-2.5d-4`
- final `nstep = 100`, expected 20 events
- final tensor strain magnitude: `0.005`
- final engineering shear for xy: `gamma_xy = 2 * tensor_E_xy = 0.010`

These are pilot-scale deformation templates, not production strain-rate
parameters.

## Leonardo Jobs And GPU Activity

| Job ID | Purpose | State | Node | Elapsed | Peak GPU activity |
| --- | --- | --- | --- | --- | --- |
| 49700757 | first cutoff attempt, exposed matrix stdin bug | COMPLETED 0:0 | lrdn3448 | 00:00:36 | 0%, 0 MiB, 59.64 W |
| 49700800 | reference SCF | COMPLETED 0:0 | lrdn0026 | 00:00:58 | 52%, 6084 MiB, 244.83 W |
| 49700996 | remaining convergence matrix | COMPLETED 0:0 | lrdn0452 | 00:04:31 | 75%, 6098 MiB, 263.12 W |
| 49701674 | zero-strain plus dt/wmass matrix | COMPLETED 0:0 | lrdn0026 | 00:26:40 | 82%, 10390 MiB, 310.64 W |
| 49701713 | active strain pilot matrix | COMPLETED 0:0 | lrdn0062 | 00:19:24 | 77%, 10390 MiB, 320.08 W |

Final tension, compression, and shear jobs were not submitted.

## Current Parameter Status

| Parameter | Selected value | Status |
| --- | --- | --- |
| `press` | `1000.0` kbar | verified unit conversion |
| `ecutwfc` | `65 Ry` | highest tested pilot setting; production convergence unresolved |
| `ecutrho` | `520 Ry` | paired with 65 Ry at ratio 8 |
| k-grid | `2 2 2 0 0 0` | bounded-pilot setting |
| structure | original converted Ice X | validated; no relaxation used |
| initial velocity card | deterministic 300 K | generated/validated; supported input path and Leonardo Ice X nstep=1 initialization validated at 295.86 K |
| `dt` | `5.0` | conservative pilot setting; final 300 K deformation runs not submitted |
| `wmass` | `65.7122994644 amu` | QE default formula; final 300 K deformation runs not submitted |
| `caimd_stride` | `5` | pilot setting |
| strain increment | `+/-5.0d-5` tensor component per MD step | pilot setting |
