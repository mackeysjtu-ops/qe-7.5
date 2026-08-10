# Audited Final Comparison Report

## Overall Status

Overall status: `PARTIAL`.

The official QE 7.4.1 process did execute, but it failed during input parsing with `unknown ion_temperature nose` before any MD trajectory was produced. The af7ebfc QE 7.5 process executed the same physical, non-CAIMD input and completed all 1000 requested steps; its `.eal` kinetic energy and terminal `.md` velocities show physical heating.

## Required Answers

1. Exact Slurm Job ID and terminal state: `50127909`; top-level state `COMPLETED` exit `0:0`. Step `50127909.0` failed, step `50127909.1` completed.
2. Did both QE processes execute? Yes. QE74 launched and failed in `iosys`; QE75 launched and completed.
3. Did both inputs contain no CAIMD keywords? Yes.
4. Were the physical inputs identical? Yes. Only `prefix`/`outdir` differ.
5. Was the initial direct-velocity temperature 300 K in both inputs? QE74 `300.000000000000 K`; QE75 `300.000000000000 K`.
6. Did official QE 7.4.1 physically heat? Not determined: no MD was accepted or run because the official module rejected `ion_temperature='nose'`.
7. Did af7ebfc physically heat with CAIMD absent? Yes. Classification `PHYSICAL_HEATING`; max `.eal`/21-DOF temperature `8123.839 K`, final `4231.686 K`.
8. Did `.eal`, QE stdout, and direct velocities agree? `.eal` and terminal direct `.md` velocities agree within `0.057 K`. QE stdout/`.tv` show the same qualitative heating but a different temperature observable/convention, ending near `3328.6 K` rather than the `.eal`/direct `4231.7 K`.
9. Which run first crossed 600 K? Only af7ebfc has a trajectory; it crossed at step `26` time `0.05` ps.
10. Which run first crossed 2000 K? Only af7ebfc has a trajectory; it crossed at step `130` time `0.258` ps.
11. Did cell instability precede or follow heating? Heating begins before significant cell deformation or cell singularity. Cell volume changed only `0.302%` by the final frame; max condition number `1.017`.
12. Is CAIMD activation required to reproduce the af7ebfc failure? No for the af7ebfc executable: the input contains no CAIMD keywords. Whether official QE 7.4.1 shows the same heating is unresolved because it did not accept the Nose input.
13. Narrowest justified conclusion: this completed job is a `PARTIAL` comparison. It proves ordinary non-CAIMD af7ebfc heats from the common state, but it does not test ordinary Nose dynamics in official QE 7.4.1 because that executable rejects the requested input.
14. Minimum next test: pristine unmodified QE 7.5 using the exact same input, executable/module path recorded, before varying Nose/PR parameters.
15. Were any raw outputs, QE source files, or existing remote calculations modified? No. Only local `analysis/`, `evidence/`, and this audited report in the local mirror were written.

## Key Evidence

- Remote workflow: `/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_ordinary_nose_qe74_vs_qe75_af7ebfc_20260723T175625Z`
- Local mirror: `/Users/mengqicheng/Downloads/Learning/qe-7.5/examples/diamond8_ordinary_nose_qe74_vs_af7ebfc_completed_20260723T185307Z`
- QE74 executable: `/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/nvhpc-24.5/quantum-espresso-7.4.1-jfm6pcx2dqkcfcwywdgoyao3m3ncvc5a/bin/pw.x` SHA256 `beb956cce572605941b791c37bc4dd30ac6cee714e18be7b93fd4ff0e47f4aca` runtime version `7.4.1`
- af7ebfc executable: `/leonardo/home/userexternal/mcheng00/software/qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x` SHA256 `89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d` runtime version `7.5`
- Input SHA256: QE74 `fb36b256662bc48b4486ab0f71e71e6c66324e91de590229c6cb1c8cf2ef35d6`; QE75 `33b584fc5ecccae84138d4bdce6de85872a218a8693deb2544bc80ed293ed51d`; pseudo `9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639`

## Caveats

- The printed `tmp/*.vel` file is not used as the vc-md physical velocity history because it remains at the initial PW global velocity values; terminal physical velocities were reconstructed from `.md` `ratd` and the final `.cel` matrix.
- The official QE 7.4.1 result is a parser/runtime compatibility failure, not a thermostat stability result.
