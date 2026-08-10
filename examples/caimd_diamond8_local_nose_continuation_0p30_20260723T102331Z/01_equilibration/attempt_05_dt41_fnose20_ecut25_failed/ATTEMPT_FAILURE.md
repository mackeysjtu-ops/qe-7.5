# Attempt 05 Failure: ionic-only temperature is too cold

This fifth equilibration attempt used:

- `dt=41.341381006`, validated as `0.002 ps` per ionic step.
- `ecutwfc/ecutrho=25/120 Ry`.
- `fnosep=20 THz`.
- `wmass=136728.78267866 AMU`.

With the corrected validator interpretation, the relevant non-CAIMD
temperature source is the ionic kinetic energy `eka` from `diamond8_equil.eal`,
not the stdout `T=` value that averages total kinetic contributions in the PR
path. The run remains invalid because the ionic-only temperature drifted cold:

- `.eal` ionic 3N-3 temperature mean: about `176.84 K`.
- `.eal` ionic 3N-3 final temperature: about `146.78 K`.
- `.eal` ionic 3N-3 final 200-step mean: about `187.63 K`.
- Final volume: `261.86449730305395 bohr^3`.
- Observed time increment: `0.002 ps`.

Diagnosis: increasing `fnosep` from 10 to 20 THz did not materially alter the
energy transfer from ions into PR cell kinetic modes. The next attempt
therefore increases the PR cell mass to slow cell motion during the short local
equilibration while keeping `cell_dynamics='pr'` and `cell_dofree='all'`.
