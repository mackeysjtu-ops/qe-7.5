# Ice X CAIMD Example Report

Status: IN PROGRESS

This report is the local framework checkpoint for the Ice X 100 GPa / nominal
300 K CAIMD example. It will be updated after deployment and bounded Leonardo
validation jobs.

## Local Preflight

- Repository root: `/Users/mengqicheng/Downloads/Learning/qe-7.5`
- Branch: `main`
- Starting HEAD: `0282a6e6f4e7d330f68ef4a2732777d8343fd716`
- Starting status: no tracked or staged diffs; only untracked local source file
  `iceX.data`.
- Git remotes: none.
- Local source file: `/Users/mengqicheng/Downloads/Learning/qe-7.5/iceX.data`
- Source SHA256: `103ccf1fe76e3b3a737914a4eeddb4016d6fb8b308129516541ddc5de08aa6c6`

## Structure Conversion

- Converter: `tools/convert_icex_data.py`
- Original copy: `structure/iceX_original.data`
- QE structure card: `structure/iceX_qe_structure.inc`
- Structure validation: `structure/structure_validation.json` and
  `structure/structure_validation.md`

Validated facts:

- 144 atoms total.
- 48 O and 96 H atoms.
- Type mapping: LAMMPS type 1=O, type 2=H from Masses-section comments and
  masses.
- Cell volume: `451.706300063683 Angstrom^3`.
- LAMMPS tilt factors: `xy=0`, `xz=0`, `yz=0`.
- Maximum roundtrip error: `8.881784197001252e-16 Angstrom`.
- No rotation, symmetrization, relaxation, or atom-order changes.

## Initial Velocities

- Generator: `scripts/generate_300K_velocities.py`
- Velocity card: `structure/iceX_velocities_300K.inc`
- Metadata: `structure/iceX_velocities_300K.json`
- Seed: `20260718`
- QE-style degrees of freedom: `3N-3 = 429`
- Final temperature by QE kinetic-energy formula: `299.9999999999999 K`
- Maximum residual COM velocity after rescale:
  `1.466377223952049e-19 Bohr/Ry atomic-time units`

## Remote Preflight Already Verified

- Remote host: `login02.leonardo.local`
- Remote user: `mcheng00`
- Remote root:
  `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919`
- `pw.x`:
  `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-v1-7aa61919/install-current-gpu/bin/pw.x`
- `pw.x` SHA256:
  `c9588df640013f098a33e33e1eec6e1f0820ed53f7488efffb0a83ec7ed25b48`
- `cp.x` SHA256:
  `2723af23df8306b5c4a52e380bd6a618606077198592fa1b6d03f8e287ba5882`
- `ldd` inspection found deployment-local FFTW/ScaLAPACK/OpenBLAS/MPI/CUDA/NVHPC
  libraries and no reported `not found` entries in the inspected output.
- Pseudopotentials:
  - H checksum
    `c7182caf52c4cad9fe912262ea40feaf4c3b757204189a14cb365869e3e0888e`
  - O checksum
    `fb96915f25537a73000b2db2fde27610f287b923f62c2ff52a674f651f394203`

## Prepared Inputs

Prepared bounded matrices:

- 4 cutoff SCF cases.
- 3 k-point SCF cases.
- 1 reference SCF case.
- 1 optional bounded `vc-relax` template.
- 1 zero-strain MD control.
- 5 dt/wmass short pilots.
- 4 active CAIMD strain pilots.
- 3 final 24 h manual-submit deformation templates.

Active CAIMD inputs use:

- `calculation='vc-md'`
- `restart_mode='from_scratch'`
- `cell_dynamics='pr'`
- `cell_dofree='all'`
- `ion_temperature='not_controlled'`
- `cell_temperature='not_controlled'`
- explicit deterministic `ATOMIC_VELOCITIES { a.u }`

## Pending

- Local static validation.
- Local commit of the framework checkpoint.
- Exact example archive creation and remote extraction.
- Remote pseudo symlink creation.
- Remote checkpoint file.
- Bounded Leonardo convergence and pilot jobs.
- Analysis updates and final manual-submit package report.

Final tension/compression/shear jobs have not been submitted.
