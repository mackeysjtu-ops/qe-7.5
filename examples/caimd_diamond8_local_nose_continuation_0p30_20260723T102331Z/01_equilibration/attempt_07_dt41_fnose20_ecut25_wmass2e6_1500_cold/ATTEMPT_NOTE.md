# Attempt 07: stable 1500-step extension, still cold

This 1500-step extension used the same stable geometry parameters as the
500-step screen:

- `dt=41.341381006`, validated as `0.002 ps` per ionic step.
- `ecutwfc/ecutrho=25/120 Ry`.
- `fnosep=20 THz`.
- `nhpcl=1`.
- `wmass=2.0e6 AMU`.

The run remained geometrically and numerically stable through 1500 steps, but
the ionic-only `.eal` temperature stayed below target:

- `.eal` ionic 3N-3 mean temperature: `186.6262 K`.
- `.eal` ionic 3N-3 final 200-step mean: `189.2243 K`.
- Final volume: `276.5579092492622 bohr^3`.
- Observed time increment: `0.002 ps`.

Diagnosis: simple extension is not enough for this single-chain Nose setup to
thermalize the initially undisplaced hard Diamond8 solid at the requested local
2 fs timestep. The next screen uses a four-thermostat Nose-Hoover chain and a
stronger but bounded `fnosep=50 THz`.
