# Ice X structure conversion validation

- Source SHA256: `103ccf1fe76e3b3a737914a4eeddb4016d6fb8b308129516541ddc5de08aa6c6`
- Structure include SHA256: `f3474180cd30cd95f32d13c9348fe1316f430028dbb398ffa68c0b69b2f92ebd`
- Atom count: 144 total, O=48, H=96
- Stoichiometry: 48 H2O units; H:O=96:48
- Cell volume: 451.706300063683 Angstrom^3
- Estimated density: 3.178909003080 g/cm^3
- Cell rows in Angstrom:
  - `9.2145100000 0.0000000000 0.0000000000`
  - `0.0000000000 6.5156430000 0.0000000000`
  - `0.0000000000 0.0000000000 7.5236160000`
- LAMMPS tilt factors xy/xz/yz: {'xy': 0.0, 'xz': 0.0, 'yz': 0.0}
- Maximum cartesian roundtrip error: 8.882e-16 Angstrom
- Minimum fractional separation norm: 0.124999918607
- Minimum H-H: 1.880901168145 Angstrom between LAMMPS atom ids [75, 117]
- Minimum H-O: 1.151810612721 Angstrom between LAMMPS atom ids [3, 49]
- Minimum O-O: 2.303624158081 Angstrom between LAMMPS atom ids [27, 41]

## Interpretation

The conversion preserves the LAMMPS cell vectors, atom order, and atom positions. No rotation, symmetrization, relaxation, coordinate rounding beyond text formatting, or species remapping beyond the Masses-section H/O identification was applied.
Coordinates are interpreted as Angstrom because this gives a physically plausible high-pressure ice density and O-H/O-O separations; interpreting the same numbers as Bohr would make the density and bond lengths unphysical for the stated Ice X cell.
