# Local static validation

Overall: PASS

- PASS `required_files`: 33 required files present
- PASS `structure`: 144 atoms, 48 O, 96 H, source/structure hashes and roundtrip validated
- PASS `velocities`: deterministic 300 K velocities, ndof=429, COM removed
- PASS `qe_inputs`: all QE inputs pass static CAIMD boundary and pressure checks
- PASS `velocity_initialization_evidence`: runtime evidence confirms explicit PR velocity initialization
- PASS `slurm`: manual deformation Slurm scripts have required Leonardo resources and checks
- PASS `pseudo`: no regular UPF binaries present; symlinks allowed/documented: []
