# Attempt 08: four-link Nose-Hoover chain, still cold

This 500-step screen used:

- `dt=41.341381006`, validated as `0.002 ps` per ionic step.
- `ecutwfc/ecutrho=25/120 Ry`.
- `fnosep=50 THz`.
- `nhpcl=4`.
- `wmass=2.0e6 AMU`.
- Initial velocities scaled to 300 K under the 3N-3 ionic convention.

The run remained numerically stable but the ionic-only `.eal` temperature still
settled near the same cold kinetic range as the single-chain attempts:

- `.eal` ionic 3N-3 mean temperature: `176.6614 K`.
- `.eal` ionic 3N-3 final 200-step mean: `179.0862 K`.
- Final volume: `302.6230930560058 bohr^3`.

Diagnosis: changing the Nose chain/frequency alone does not solve the
cold-start effect from a perfect zero-displacement Diamond8 structure. The next
screen keeps the Nose target `tempw=300 K` but uses a 450 K initial velocity
warm start so that kinetic and potential energy sharing can settle near the
target ionic temperature during the short local equilibration.
