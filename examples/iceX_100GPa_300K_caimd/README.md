# Ice X 100 GPa / 300 K CAIMD example

This directory is a reproducible Quantum ESPRESSO 7.5 CAIMD V1 example for an
Ice X LAMMPS structure at 100 GPa with a deterministic nominal 300 K initial
velocity set.

The example is intentionally a bounded validation and manual-submission package.
It does not claim an isothermal stress-strain curve, a publication-quality strain
rate, or fully converged material response. Active CAIMD inputs use
`ion_temperature='not_controlled'` and `cell_temperature='not_controlled'`; the
300 K value is only the initial ionic kinetic temperature.

## Layout

- `structure/`: original LAMMPS data copy, QE structure card, 300 K velocity
  card, and structure/velocity validation evidence.
- `pseudo/`: remote pseudopotential provenance and checksums. UPF binaries are
  not stored in this Git example.
- `convergence/`: bounded cutoff and k-point SCF inputs plus analysis outputs.
- `reference/`: static SCF, optional bounded `vc-relax`, and zero-strain MD
  controls.
- `pilot/`: bounded dt/wmass and strain-increment pilot inputs.
- `deformation/`: final manual-submit tension, compression, and shear templates.
- `scripts/`: conversion, generation, validation, analysis, and Slurm helpers.
- `evidence/`: local and remote validation records.

## Remote pseudopotential setup

On Leonardo, after extracting this example under:

```bash
ROOT=/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919
cd "$ROOT/example_iceX_100GPa_300K_caimd"
ln -sfn /leonardo/home/userexternal/mcheng00/psd/H.upf pseudo/H.upf
ln -sfn /leonardo/home/userexternal/mcheng00/psd/O.upf pseudo/O.upf
```

The Slurm scripts verify the expected `pw.x` checksum and require these pseudo
symlinks before launching QE.

## Bounded validation jobs

The convergence matrix is capped at four cutoffs and three k-point grids:

```bash
cd "$ROOT/example_iceX_100GPa_300K_caimd"
sbatch scripts/run_convergence.slurm
```

The pilot matrix is capped at one zero-strain control, five dt/wmass cases, and
four short active strain cases:

```bash
cd "$ROOT/example_iceX_100GPa_300K_caimd"
sbatch --export=CASES_FILE=pilot/pilot_cases.txt scripts/run_convergence.slurm
```

## Final deformation templates

These final jobs are intentionally not submitted by this preparation workflow.
Submit them manually only after reviewing the bounded validation evidence:

```bash
cd "$ROOT/example_iceX_100GPa_300K_caimd/deformation/tension"
sbatch run.slurm

cd "$ROOT/example_iceX_100GPa_300K_caimd/deformation/compression"
sbatch run.slurm

cd "$ROOT/example_iceX_100GPa_300K_caimd/deformation/shear"
sbatch run.slurm
```

## Local regeneration

From the QE source root:

```bash
python3 examples/iceX_100GPa_300K_caimd/tools/convert_icex_data.py
python3 examples/iceX_100GPa_300K_caimd/scripts/generate_300K_velocities.py
python3 examples/iceX_100GPa_300K_caimd/scripts/prepare_example_inputs.py
python3 examples/iceX_100GPa_300K_caimd/scripts/write_case_slurm.py
python3 examples/iceX_100GPa_300K_caimd/scripts/validate_example.py
```
