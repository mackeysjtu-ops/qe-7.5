# Attempt 02 Failure: `dt=82.682762012`

This second equilibration attempt used the scaled Parrinello-Rahman cell mass
`wmass=546915.13071464` but retained `dt=82.682762012`.

Quantum ESPRESSO completed 500 ionic steps with return code 0 and `JOB DONE`,
but the run is invalid for two independent reasons:

- The stdout trajectory time advances by `0.004 ps` per ionic step. Therefore
  `dt=82.682762012` corresponds to 4 fs in this PW `vc-md` path, not the 2 fs
  timestep required by the task.
- The trajectory still undergoes a clear runaway: the local validation reported
  final volume `10181335.385172607 bohr^3`, maximum parsed stdout temperature
  `9888.0 K`, and final minimum periodic C-C distance `36.40807066573568 bohr`.

Diagnosis: this attempt cannot be used as the requested 2 fs local
thermostatted PR equilibration. The follow-up workflow sets
`dt=41.341381006`, which should correspond to a 0.002 ps/step stdout increment.
