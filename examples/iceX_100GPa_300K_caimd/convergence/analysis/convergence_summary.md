# Ice X bounded convergence summary

- Completed cases: 18
- Pending validation cases: []
- Final not-submitted cases: ['deformation/compression', 'deformation/shear', 'deformation/tension']
- Optional not-run cases: ['reference/optional_vc_relax']

| Case | Status | Ecut Ry | K-grid | Energy Ry | Pressure GPa | Force Ry/bohr |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| convergence/ecut/ecut_35Ry | ok | 35.0 | [1, 1, 1] | -2109.44021125 | 95.25800000000001 | 0.008658 |
| convergence/ecut/ecut_45Ry | ok | 45.0 | [1, 1, 1] | -2109.57776344 | 100.584 | 0.009089 |
| convergence/ecut/ecut_55Ry | ok | 55.0 | [1, 1, 1] | -2109.6344461 | 99.811 | 0.009149 |
| convergence/ecut/ecut_65Ry | ok | 65.0 | [1, 1, 1] | -2109.68430794 | 100.322 | 0.009073 |
| convergence/kpoints/k_2x2x2 | ok | 65.0 | [2, 2, 2] | -2109.79374555 | 99.623 | 0.000277 |
| convergence/kpoints/k_3x3x3 | ok | 65.0 | [3, 3, 3] | -2109.79416654 | 99.618 | 3.1e-05 |
| convergence/kpoints/k_gamma | ok | 65.0 | [1, 1, 1] | -2109.68430794 | 100.322 | 0.009073 |
| deformation/compression | not_submitted | 65.0 | [2, 2, 2] |  |  |  |
| deformation/shear | not_submitted | 65.0 | [2, 2, 2] |  |  |  |
| deformation/tension | not_submitted | 65.0 | [2, 2, 2] |  |  |  |
| pilot/compression/xx_stepm5e-5_stride5 | ok | 65.0 | [2, 2, 2] | -2109.77188913 | 100.19 | 0.004908 |
| pilot/parameter_screen/dt10_wmass_default | ok | 65.0 | [2, 2, 2] | -2109.78708206 | 99.801 | 0.000311 |
| pilot/parameter_screen/dt5_wmass_default | ok | 65.0 | [2, 2, 2] | -2109.79196836 | 99.67699999999999 | 0.000101 |
| pilot/parameter_screen/dt5_wmass_double | ok | 65.0 | [2, 2, 2] | -2109.79284746 | 99.654 | 0.000142 |
| pilot/parameter_screen/dt5_wmass_half | ok | 65.0 | [2, 2, 2] | -2109.79026555 | 99.72 | 6.3e-05 |
| pilot/parameter_screen/dt7p5_wmass_default | ok | 65.0 | [2, 2, 2] | -2109.78985226 | 99.731 | 0.000117 |
| pilot/shear/xy_step5e-5_stride5 | ok | 65.0 | [2, 2, 2] | -2109.79158027 | 99.684 | 0.004842 |
| pilot/tension/xx_step1e-5_stride5 | ok | 65.0 | [2, 2, 2] | -2109.79639249 | 99.562 | 0.001037 |
| pilot/tension/xx_step5e-5_stride5 | ok | 65.0 | [2, 2, 2] | -2109.81254132 | 99.147 | 0.004962 |
| reference/optional_vc_relax | optional_not_run | 65.0 | [2, 2, 2] |  |  |  |
| reference/scf | ok | 65.0 | [2, 2, 2] | -2109.79374555 | 99.623 | 0.000277 |
| reference/zero_strain | ok | 65.0 | [2, 2, 2] | -2109.79196836 | 99.67699999999999 | 0.000101 |

