# Diamond8 vc-relax — QE 7.5 (custom, CAIMD off) vs QE 7.4.1 (official) — audited

**Overall status: `BOTH_VCRELAX_VALID_AND_CONSISTENT`**
**Downstream decision: `DOWNSTREAM_READY_USING_QE75_STRUCTURE`**

Read-only result audit. No job was submitted, rerun, cancelled, resubmitted, or rebuilt; no source,
executable, pseudopotential, or raw output was modified; no build/commit/push. Exactly one grouped
`sacct` query was issued. Scope confined to `/leonardo/home/userexternal/mcheng00/workspace/Diamond/`;
the **SinglePhase** H/O / water / fixed-O / λ-integration / 3000 K workflows were **excluded entirely**
(not inspected or cited).

- Remote WorkDir: `/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_vcrelax_qe75_vs_qe74_ibrav_20260728T205518Z`
- Local mirror: `examples/diamond8_vcrelax_qe75_vs_qe74_completed_20260728T234619Z/`

## 1. Slurm (single grouped `sacct`)

| Job | ID | Name | State | Exit | Elapsed | Node | pw.x step |
|---|---|---|---|---|---|---|---|
| JOB75 | 50605866 | d8_vcr75_ibrav | COMPLETED | 0:0 | 00:00:37 | lrdn2889 | COMPLETED 0:0 (3 s) |
| JOB74 | 50605867 | d8_vcr74_ibrav | COMPLETED | 0:0 | 00:00:36 | lrdn2931 | COMPLETED 0:0 (3 s) |

Independent jobs, no dependency. Neither RUNNING nor PENDING. Both scripts exited 0 (internal gate
passed: `GATE_qe75=0`, `GATE_qe74=0`).

## 2. Provenance (verified, all match expected)

| | QE 7.5 (JOB75) | QE 7.4.1 (JOB74) |
|---|---|---|
| banner | v.7.5 | v.7.4.1 |
| executable | `…/qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x` | `/leonardo/prod/spack/06/…/quantum-espresso-7.4.1-jfm6pcx2…/bin/pw.x` |
| exe SHA256 | `89d4fda8…a10e7d` ✓ | `beb956cc…f4aca` ✓ |
| module | (custom build env) | `quantum-espresso/7.4.1--hpcx-mpi--2.19--nvhpc--24.5-openblas` |
| input SHA256 | `5e1e6e94…3660` | `11925e5a…8571` |

Pseudopotential `C.UPF` SHA256 `9c06aa64…978639` ✓ (matches). Local repo `main @ af7ebfc757…` ✓.

## 3. Input contract (both cases, verified independently)

8 C atoms, no H/O; `calculation='vc-relax'`; `ibrav=1`; `cell_dofree='ibrav'`; **`cell_dofree='volume'`
absent**; **no CAIMD**; **no explicit velocities**; `ecutwfc=80` / `ecutrho=320` Ry; `K_POINTS automatic
2 2 2 0 0 0`; `press=0.0` kbar; `conv_thr=1.0d-8`; same initial cell (`celldm(1)=6.7406530895` bohr),
same atomic positions/ordering, same C.UPF. **Normalized physical-input diff is EMPTY** (ignoring only
comments, `prefix`, `outdir`) — see `analysis_audited/normalized_input_diff.txt`.

## 4. Per-version results (independent parse of raw `qe.out`)

Both versions ran SCF and converged BFGS explicitly: `bfgs converged in 4 scf cycles and 2 bfgs steps`
+ `End of BFGS Geometry Optimization` (not a max-step stop). Pressure trajectory (kbar):
`-154.33 → 221.52 → -4.50 → 0.25 → 0.25`. Volume (bohr³): `306.271 → 282.798 → 295.861 → 295.559 → 295.559`.

| quantity | QE 7.5 | QE 7.4.1 |
|---|---|---|
| classification | **QE75_VCRELAX_VALID** | **QE74_VCRELAX_VALID** |
| final pressure | 0.25 kbar (0.025 GPa) | 0.25 kbar (0.025 GPa) |
| final volume | 295.5588 bohr³ (43.7973 Å³) | 295.5588 bohr³ (43.7973 Å³) |
| final lattice const a | 6.661131 bohr (3.524919 Å) | 6.661131 bohr (3.524919 Å) |
| final max atomic force | 3.98e-7 Ry/bohr (1.02e-5 eV/Å) | 3.98e-7 Ry/bohr (1.02e-5 eV/Å) |
| final total energy | −91.50602495 Ry | −91.50602495 Ry |
| energy per atom | −11.438253 Ry (−155.62536 eV) | −11.438253 Ry (−155.62536 eV) |
| BFGS steps | 2 | 2 |
| cell determinant | 295.5588 bohr³ (>0) | 295.5588 bohr³ (>0) |
| cell condition number | 1.0 | 1.0 |
| min C–C distance | 1.526335 Å | 1.526335 Å |
| cubic / shear | cubic, no shear | cubic, no shear |
| SCF/NaN/Inf/fatal | none | none |

Each satisfies every `VALID` criterion (COMPLETED 0:0, RC 0, JOB DONE, explicit BFGS convergence, 8 C,
positive finite determinant, max force ≤ 1e-4 Ry/bohr, |P| ≤ 5 kbar, cubic, no non-finite values).

## 5. QE 7.5 vs QE 7.4.1 comparison

| metric | value | tolerance | pass |
|---|---|---|---|
| final pressure difference | 0.00 kbar | ≤5 kbar each | ✓ |
| final volume rel. diff | 0.0 % | ≤0.5 % | ✓ |
| final lattice-const rel. diff | 0.0 % | ≤0.2 % | ✓ |
| normalized Cartesian RMS diff | 4.3e-10 Å | ≤1.0e-3 Å | ✓ |
| max single-atom diff | 6.1e-10 Å | ≤5.0e-3 Å | ✓ |
| energy/atom diff | ~0 eV | diagnostic | — |
| relaxation-step diff | 0 | diagnostic | — |
| both cells cubic | yes | required | ✓ |

The two runs are **bit-identical to all printed digits** except a ~1e-10 fractional round-off on the
four tetrahedral atoms. **Consistency verdict: `CONSISTENT`.**

## 6. Gate audit

Both Slurm scripts exited 0 and both results are genuinely `VALID`, so exit 0 correctly reflects a
converged, in-tolerance vc-relax. **Verdict: `GATES_CORRECT`.**

## 7. Required answers

1. **Both Slurm jobs completed?** Yes — both COMPLETED, exit 0:0.
2. **Both QE calculations ran SCF?** Yes — both printed SCF convergence and `JOB DONE`; pw.x steps COMPLETED 0:0.
3. **Did `cell_dofree='ibrav'` avoid the previous failure?** Yes — no `Isotropic expansion…ibrav=1` error in either run; the earlier `cell_dofree='volume'` crash is absent.
4. **QE 7.5 explicitly converged?** Yes — `bfgs converged in 4 scf cycles and 2 bfgs steps` + `End of BFGS Geometry Optimization`.
5. **QE 7.4.1 explicitly converged?** Yes — identical convergence markers.
6. **Final P / V / a / force / energy:** both versions — P 0.25 kbar (0.025 GPa); V 295.5588 bohr³ (43.7973 Å³); a 6.661131 bohr (3.524919 Å); max force 3.98e-7 Ry/bohr (1.02e-5 eV/Å); E −91.50602495 Ry (−155.62536 eV/atom).
7. **Both final cells cubic?** Yes — diagonal cells, off-diagonal 0, condition number 1.0.
8. **How different are the final structures?** Position RMS 4.3e-10 Å, max atom 6.1e-10 Å, ΔV 0 %, Δa 0 % — negligible round-off.
9. **Consistent under the tolerances?** Yes — `CONSISTENT` (all structural tolerances satisfied).
10. **Gates behaved correctly?** Yes — `GATES_CORRECT`.
11. **Downstream chain ready to be prepared?** Yes — but only prepared, not submitted in this task.
12. **Which structure downstream?** The **QE 7.5** final relaxed structure (`qe75/analysis/final_structure.json`; audited copy `analysis_audited/qe75_final_structure.json`).
13. **Single minimum next action:** In a later, submission-authorized session, prepare the downstream NVT/NPT/CAIMD chain from the QE 7.5 relaxed structure (a = 6.661131 bohr) — do not submit it here.

## 8. Safety confirmation

No job submitted/rerun/cancelled/resubmitted/rebuilt; no dependency altered; no source/executable/
pseudopotential/raw output modified; no build/commit/push. One grouped `sacct` only. New products written
only under `analysis_audited/`, `evidence_audited/`, and this report. SinglePhase untouched. The
downstream NVT/NPT/CAIMD chain was **not** submitted.
