# Implementation Confirmation

Generated: 2026-07-23T12:10:34Z

Scope: local-only confirmation for active CAIMD plus ionic Nose thermostat in
the current Quantum ESPRESSO source tree.

## Confirmed Behavior

- Active CAIMD supports `ion_temperature='nose'` and `ion_temperature='not_controlled'`.
- Unsupported active-CAIMD ionic thermostat modes remain rejected in
  `PW/src/input.f90` by the active CAIMD allowlist.
- PW `calculation='vc-md'` with `cell_dynamics='pr'` maps to `calc='cd'`, so
  variable-cell MD enters `PW/src/vcsmd.f90` and `PW/src/vcsubs.f90`.
- Ionic Nose setup uses QE's `get_ndof()` convention; for no fixed coordinates
  this gives `3*nat - 3`.
- `vcsmd` applies `ions_nosevel` and `ions_noseupd` around `vcmove`.
- The Nose friction term acts on ionic velocity histories (`ratd`/related
  histories), not on the imposed CAIMD strain schedule itself.
- CAIMD projection transforms `ratd`, `rat2d`, and `rat2di` so the committed
  physical Cartesian ionic velocity is preserved across a projected cell event.
- After final CAIMD projection, `vc_refresh_projected_nose_ekin` synchronizes
  projected ionic kinetic energy with the Nose chain kinetic inputs.
- `.caimd` column `temperature_K` is based on ionic kinetic energy. Rows are
  premove diagnostics; the first row may fall back to the input target
  temperature before a committed ionic kinetic temperature exists.
- CAIMD component order is `(xx, yy, zz, yz, xz, xy)`.
- `caimd_mask(i)=.true.` means component `i` is imposed or held by CAIMD.
- Unmasked symmetric components retain the ordinary Parrinello-Rahman candidate
  response before final projection.

No source redesign or source modification was performed during this audit.
