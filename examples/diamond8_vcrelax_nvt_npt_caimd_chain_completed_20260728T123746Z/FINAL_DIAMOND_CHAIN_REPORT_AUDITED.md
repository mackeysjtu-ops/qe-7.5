# Diamond8 vc-relax / NVT / NPT / CAIMD dependency chain — completed-run audit

**Overall status: `VCRELAX_FAILED`**  (earliest failed stage: JOB 1 vc-relax; JOBS 2–4 never ran)

Read-only result audit. No job was submitted, rerun, cancelled, resubmitted, rebuilt, committed, or
pushed. Scope was confined strictly to
`/leonardo/home/userexternal/mcheng00/workspace/Diamond/`. The unrelated **SinglePhase** H/O
lambda-integration workflow (1152 atoms, 3000 K, fixed-O, jobs 50468684/85/86/91, water, λ-integration,
1152-atom / 3000 K trees) was **excluded entirely** — not inspected, cited, or used.

- Remote WorkDir: `/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_vcrelax_nvt_npt_caimd_chain_20260728T120614Z`
- Local mirror: `examples/diamond8_vcrelax_nvt_npt_caimd_chain_completed_20260728T123746Z/`
- Submission timestamp (from `SUBMISSION_REPORT.md`): 2026-07-28T12:08:47Z

## 1. Slurm chain (single grouped `sacct` query)

| Stage | job name | Job ID | dependency | Slurm state | exit | elapsed |
|---|---|---|---|---|---|---|
| JOB 1 vc-relax | d8_vcrelax | 50558068 | (none) | **FAILED** | 1:0 | 00:00:06 |
| JOB 2 fixed NVT | d8_fixed_nvt | 50558069 | afterok:50558068 | CANCELLED | 0:0 | — |
| JOB 3 free NPT | d8_free_npt | 50558070 | afterok:50558069 | CANCELLED | 0:0 | — |
| JOB 4 CAIMD 4A/4B | d8_caimd_zn | 50558071 | afterok:50558070 | CANCELLED (Reason=Dependency) | 0:0 | — |

JOB 1 exited non-zero after 6 s; its `pw.x` step failed (exit 1:0) after ~2 s. JOBS 2–4 were CANCELLED
because their `afterok` dependency was never satisfied. No job is still active.

## 2. Root cause of JOB 1 failure

`01_vc_relax/CRASH`:

```
Error in routine cell_dofree (1):
Isotropic expansion is only allowed for ibrav=1; i.e. for simple cubic
```

QE aborted **at input parsing**, before any SCF or BFGS step. Cause (verified by read-only source
inspection at the audited commit `af7ebfc`):

- `01_vc_relax.in` requests `cell_dofree='volume'` with `ibrav=1`.
- In `Modules/cell_base.f90`, `init_dofree` guards the `'volume'` case with `IF (ibrav /= 1) CALL errore(...)`,
  reading the **module-level** `ibrav`.
- In `PW/src/input.f90`, `init_dofree(cell_dofree)` is called (line ~1306) **before** the module `ibrav`
  is populated by `cell_base_init` (invoked later inside `structure_iosys`). So at the guard's call site
  the module `ibrav` is not yet 1, and QE aborts even though the `&SYSTEM` namelist correctly sets `ibrav=1`.

**Conclusion:** the `ibrav=1 + cell_dofree='volume'` design used for JOB 1 does not work in this pw.x
build. This is a workflow-input defect, not a hardware/convergence failure. Classification:
`VCRELAX_NUMERICAL_FAILURE` (parse-time abort; no pressure/volume/force data produced).

## 3. Provenance (verified read-only, all match expected)

- QE version banner: **v.7.5**.
- Executable `.../qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x`
  SHA256 `89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d` — **matches**.
- Carbon pseudopotential SHA256 `9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639`
  — **matches** (norm-conserving LDA C).
- JOB 1 input: exactly **8 C atoms** (no H/O), `ecutwfc=80` Ry, `ecutrho=320` Ry,
  `K_POINTS automatic 2 2 2 0 0 0`, `calculation='vc-relax'`, `press=0.0`.
- Local repo `/Users/mengqicheng/Downloads/Learning/qe-7.5` at `af7ebfc757c32e0b4f050ed9f2208e01020b2ee3`.

## 4. Per-stage status

| Stage | Job ID | ran | classification |
|---|---|---|---|
| JOB 1 vc-relax | 50558068 | yes (crashed at parse) | `VCRELAX_NUMERICAL_FAILURE` |
| JOB 2 fixed NVT | 50558069 | no | NOT_RUN (NVT_INCOMPLETE) |
| JOB 3 free NPT | 50558070 | no | NOT_RUN (NPT_INCOMPLETE) |
| JOB 4A CAIMD zero | 50558071 | no | NOT_RUN (CAIMD_ZERO_INCOMPLETE) |
| JOB 4B CAIMD nonzero | 50558071 | no | NOT_RUN (moot — JOB 4 never started) |

Only the JOB 1 input (`inputs/generated/01_vc_relax.in`) was ever generated; the runtime-generated
inputs for JOBS 2–4 do not exist because those jobs were cancelled before running. There are **no
explicit atomic velocities** anywhere (JOB 1 is vc-relax with no velocity card; the MD inputs that would
carry velocities were never produced). No structure handoff (01→02→03→04, or zero→nonzero) ever occurred.
No plots were generated — there is no trajectory/history data (QE aborted before the first SCF step; the
JOB 1 pressure/volume/force history CSVs are empty by construction).

## 5. Gate and dependency audit

- JOB 1 returned non-zero (QE code 1) and the Slurm job exited FAILED, so its `afterok` was correctly
  never satisfied.
- The dependency graph is the correct **`afterok`** chain (JOB2←JOB1, JOB3←JOB2, JOB4←JOB3), not
  `afterany`. Upstream failure cascaded exactly as intended: JOBS 2–4 were CANCELLED, none ran.
- **Verdict: gates and dependencies behaved correctly.** There was no dependency/gate defect; the only
  defect is inside the JOB 1 input.

## 6. Required answers

1. vc-relax reached 0 kbar: **No** — crashed at input parsing; no pressure computed.
2. Fixed-cell NVT maintained 300 K: **N/A** — JOB 2 never ran.
3. Fixed-cell NVT pressure: **N/A** — JOB 2 never ran.
4. Free-cell NPT maintained 300 K: **N/A** — JOB 3 never ran.
5. Free-cell NPT maintained 0 kbar: **N/A** — JOB 3 never ran.
6. Ionic Nose active: **N/A** — no MD stage ran.
7. Cell Nose active (NPT): **N/A** — JOB 3 never ran.
8. CAIMD disabled in JOB 3: **N/A** — JOB 3 never ran.
9. CAIMD zero-strain stable: **N/A** — JOB 4 never ran.
10. Small nonzero strain applied correctly: **N/A** — JOB 4 never ran.
11. Earliest failed stage: **vc-relax (JOB 1)**.
12. Ready for a larger/longer test: **No** — the vc-relax baseline never succeeded.
13. QE source modified: **No.**
14. Remote WorkDir: `.../Diamond/diamond8_vcrelax_nvt_npt_caimd_chain_20260728T120614Z`.
15. Four Job IDs: 50558068 / 50558069 / 50558070 / 50558071.
16. Slurm states/exit codes: FAILED(1:0) / CANCELLED / CANCELLED / CANCELLED.
17. Dependency-chain verdict: **CORRECT** (`afterok`, upstream failure cascaded).
18. Executable SHA: `89d4fda8…` — matches; pseudo SHA `9c06aa64…` — matches.
19. System identity: **8 C atoms**, ecutwfc 80 / ecutrho 320 Ry, k-points 2×2×2 — matches expected.
20. JOB 1 crash cause: `cell_dofree='volume'` hits the `ibrav/=1` guard in `init_dofree` at parse time
    (module `ibrav` not yet set), despite input `ibrav=1`.
21. JOBS 2–4: **not run** (correctly CANCELLED by `afterok`).
22. Explicit velocities: **none** present (no MD input generated).
23. Minimum single next action: fix JOB 1 so cell relaxation does not use the broken
    `cell_dofree='volume'` path in this build (use a `cell_dofree` accepted by this pw.x, or a build in
    which the `ibrav` guard sees `ibrav=1`), then resubmit the chain. **Recommendation only — not applied.**

## 7. Safety confirmation

No job was submitted, rerun, cancelled, resubmitted, or rebuilt; no dependency was altered. No QE
source, executable, pseudopotential, or remote raw output was modified. No build, commit, or push
occurred. Exactly one grouped `sacct` query was issued. Only new audit files under `analysis_audited/`,
`evidence_audited/`, and this report were written; mirrored raw outputs were not edited. The SinglePhase
workflow was not touched.
