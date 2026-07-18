# Ice X 100 GPa / nominal 300 K CAIMD example

This directory is a reproducible Quantum ESPRESSO 7.5 CAIMD V1 example package
for the local Ice X LAMMPS structure `iceX.data`.

The package contains a validated structure conversion, pseudopotential
provenance, bounded Leonardo convergence and pilot evidence, analysis scripts,
and three final manually submittable deformation templates. It does not modify
QE or CAIMD source code.

## Important Temperature Caveat

A deterministic 300 K `ATOMIC_VELOCITIES { a.u }` card is generated and included
in all MD inputs. The velocity metadata verifies 300 K after center-of-mass
removal and rescaling.

However, the current QE 7.5 `vc-md`/Parrinello-Rahman path used by active CAIMD
does not actually initialize the variable-cell MD trajectory from that velocity
card when `ion_temperature='not_controlled'`. Source inspection shows:

- active CAIMD rejects ionic thermostats/rescaling and allows only
  `ion_temperature='not_controlled'`;
- that mode sets the internal `temperature` variable to `0.0`;
- `vcsmd` calls `vcinit(..., temperature, ...)`;
- `vcinit` calls `ranv(temp, ...)` to initialize variable-cell MD velocities.

The completed Leonardo pilot outputs therefore start and remain near 0 K, not
300 K. The final templates keep the 300 K velocity card for reproducibility and
for future source work, but with the current source they are best understood as
nominal-300 K input decks whose active `vc-md` trajectory initializes near 0 K.

## Layout

- `structure/`: original LAMMPS data copy, QE structure card, 300 K velocity
  card, and validation evidence.
- `pseudo/`: remote pseudopotential provenance and checksums. UPF binaries are
  not stored in this Git example.
- `convergence/`: bounded cutoff and k-point SCF inputs plus analysis outputs.
- `reference/`: static SCF, optional bounded `vc-relax`, and zero-strain MD
  controls.
- `pilot/`: bounded dt/wmass and strain-increment pilot inputs.
- `deformation/`: final manual-submit tension, compression, and shear templates.
- `scripts/`: conversion, generation, validation, analysis, and Slurm helpers.
- `evidence/`: local validation records and compact Leonardo job evidence.

## Remote Pseudopotential Setup

On Leonardo, after extracting this example under:

```bash
ROOT=/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919
cd "$ROOT/example_iceX_100GPa_300K_caimd"
ln -sfn /leonardo/home/userexternal/mcheng00/psd/H.upf pseudo/H.upf
ln -sfn /leonardo/home/userexternal/mcheng00/psd/O.upf pseudo/O.upf
```

The Slurm scripts verify the expected custom CAIMD `pw.x` checksum and require
these pseudo symlinks before launching QE.

## Completed Bounded Validation

The Leonardo short validation completed:

- cutoff matrix: 35, 45, 55, 65 Ry at gamma;
- k-point matrix: gamma, 2x2x2, 3x3x3 at 65 Ry;
- reference SCF at 100 GPa input pressure;
- zero-strain and five dt/wmass `vc-md` pilot cases;
- short active tension, compression, and tensor xy shear CAIMD pilots.

See `PARAMETER_SELECTION.md`, `ICE_X_CAIMD_EXAMPLE_REPORT.md`, and
`convergence/analysis/`.

## Final Deformation Templates

These final jobs were intentionally not submitted by this preparation workflow.
Submit them manually only after reviewing the evidence and the temperature
caveat:

```bash
ROOT=/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919

cd "$ROOT/example_iceX_100GPa_300K_caimd/deformation/tension"
sbatch run.slurm

cd "$ROOT/example_iceX_100GPa_300K_caimd/deformation/compression"
sbatch run.slurm

cd "$ROOT/example_iceX_100GPa_300K_caimd/deformation/shear"
sbatch run.slurm
```

Each final directory contains `qe.in`, `run.slurm`, `README.md`,
`expected_files.txt`, `analyze.sh`, and `CASE_METADATA.json`.

## Local Regeneration

From the QE source root:

```bash
python3 examples/iceX_100GPa_300K_caimd/tools/convert_icex_data.py
python3 examples/iceX_100GPa_300K_caimd/scripts/generate_300K_velocities.py
python3 examples/iceX_100GPa_300K_caimd/scripts/prepare_example_inputs.py
python3 examples/iceX_100GPa_300K_caimd/scripts/write_case_slurm.py
python3 examples/iceX_100GPa_300K_caimd/scripts/validate_example.py
```
