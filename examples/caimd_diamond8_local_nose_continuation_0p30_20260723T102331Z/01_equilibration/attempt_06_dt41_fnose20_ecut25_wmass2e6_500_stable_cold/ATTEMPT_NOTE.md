# Attempt 06: stable 500-step screen, extended for equilibration

This 500-step screen used:

- `dt=41.341381006`, validated as `0.002 ps` per ionic step.
- `ecutwfc/ecutrho=25/120 Ry`.
- `fnosep=20 THz`.
- `wmass=2.0e6 AMU`.

The trajectory passed the non-temperature stability checks: SCF completed,
`JOB DONE` was present, 500 cell/position/velocity blocks were written, the
cell determinant remained positive and sane, and the volume changed only from
about `306.27 bohr^3` to `302.62 bohr^3`.

The ionic-only `.eal` temperature was still below the local validation target:
mean `177.81 K`, final 200-step mean `181.95 K`. This is not a safety stop;
it is an under-equilibrated but stable screen. Per the task instruction, the
next run extends the ordinary Nose+PR equilibration from 500 to 1500 steps.
