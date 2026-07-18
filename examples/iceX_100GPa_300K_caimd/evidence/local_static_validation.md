# Local static validation

Overall: PASS

- PASS `required_files`: 23 required files present
- PASS `structure`: 144 atoms, 48 O, 96 H, source/structure hashes and roundtrip validated
- PASS `velocities`: deterministic 300 K velocities, ndof=429, COM removed
- PASS `qe_inputs`: all QE inputs pass static CAIMD boundary and pressure checks
- PASS `slurm`: manual deformation Slurm scripts have required Leonardo resources and checks
- PASS `pseudo`: no local UPF binaries present; remote symlinks documented
