# Diamond8 vc-relax / NVT / NPT / CAIMD dependency chain — audit

**Overall status: `DIAMOND_CHAIN_NOT_SUBMITTED`**

Read-only recovery audit, 2026-07-28. Scope was confined strictly to
`/leonardo/home/userexternal/mcheng00/workspace/Diamond/`. The unrelated SinglePhase H/O
lambda-integration workflow (1152 atoms, 3000 K, fixed-O, jobs 50468684/85/86/91) was **excluded
entirely** — not inspected, cited, or used.

## 1. Finding

The four-job Diamond CAIMD dependency chain was **prepared but never submitted**. The prior session
ended (usage limit, then redirected to an unrelated task) **before** the workflow was transferred to
Leonardo and before the four `sbatch` commands were issued.

### Evidence of non-submission (all read-only)

- **No matching workflow root.** `ls /leonardo/.../Diamond/diamond8_vcrelax_nvt_npt_caimd_chain_*`
  → no matching directory. The Diamond workspace contains several *other* diamond workflows
  (initial-vs-nose, native-nose-caimd-zero-strain-5000, qe74-vs-qe75, diamond111 studies) but **not**
  the vc-relax/NVT/NPT/CAIMD chain.
- **No submission records.** No `DEPENDENCY_CHAIN.json`, no `01_vc_relax.slurm` /
  `04_caimd_zero_nonzero.slurm`, and no `submitted_jobs`/`WORKFLOW_METADATA.json` for this chain
  exist anywhere under `Diamond/` (`find -maxdepth 3`).
- **No Slurm accounting.** One grouped `sacct --starttime 2026-07-27 --name=d8_vcrelax,d8_fixed_nvt,d8_free_npt,d8_caimd_zn`
  returned **zero rows**. The four jobs were never submitted, so no Job IDs exist. (The only recent
  Diamond job, 50434898, is the earlier `diamond8_native_nose_caimd_zero_strain_5000` workflow — a
  different test, not part of this chain.)
- **No Job IDs were fabricated.** Because nothing was submitted, `JOB1..JOB4` IDs are `null`.

### Preparation that does exist (local only)

The prior session had fully **built and locally validated** the chain in a scratchpad staging
directory. Those artifacts are copied here under `prepared_not_submitted/` as a record:

- `scripts/01_vc_relax.slurm`, `02_fixed_nvt.slurm`, `03_free_npt.slurm`, `04_caimd_zero_nonzero.slurm`
- `scripts/analyze_stage.py`, `extract_final_structure.py`, `build_overall_report.py`
- `inputs/generated/01_vc_relax.in` (ibrav=1 cubic, 8 C atoms, ecutwfc 80 / ecutrho 320 Ry,
  `K_POINTS automatic 2 2 2 0 0 0`, `vc-relax`, `cell_dofree='volume'`, `press=0.0`)
- `pseudo/C.UPF`

These were never `rsync`-transferred to Leonardo and never submitted, so **no calculation ran** for
this chain.

## 2. Provenance (verified read-only, matches expected)

- Local source repo: `/Users/mengqicheng/Downloads/Learning/qe-7.5`, branch `main`,
  HEAD `af7ebfc757c32e0b4f050ed9f2208e01020b2ee3` — matches expected.
- Executable `/leonardo/.../qe-7.5-caimd-thermostat-af7ebfc/install-current-gpu/bin/pw.x` exists,
  SHA256 `89d4fda827d2cc22a3b4c1670e90dc40ce629af42096c785b55c70c355a10e7d` — matches expected.
- Carbon pseudopotential `C.UPF` SHA256 `9c06aa642d22981491b32a18ec129086fe8366164984930389e4fc071c978639`
  — matches expected (norm-conserving LDA C).

## 3. Per-stage status

| Stage | intended job name | intended dependency | Job ID | Slurm state | status |
|---|---|---|---|---|---|
| JOB 1 vc-relax | d8_vcrelax | (none) | — | — | NOT_SUBMITTED |
| JOB 2 fixed NVT | d8_fixed_nvt | afterok:JOB1 | — | — | NOT_SUBMITTED |
| JOB 3 free NPT | d8_free_npt | afterok:JOB2 | — | — | NOT_SUBMITTED |
| JOB 4 CAIMD 4A/4B | d8_caimd_zn | afterok:JOB3 | — | — | NOT_SUBMITTED |

No stage ran; phases 3–10 (Slurm accounting, runtime contract, per-stage temperature/pressure/CAIMD
audits, gate/dependency audit) are **not applicable** because there is no calculation output to audit.
No plots were generated (no data). The intended dependency graph in the prepared scripts was the
correct `afterok` chain (JOB2←JOB1, JOB3←JOB2, JOB4←JOB3), not `afterany`, but it was never exercised.

## 4. Overall verdict and recommended next step

**`DIAMOND_CHAIN_NOT_SUBMITTED`.** The chain exists only as validated local preparation; nothing was
deployed or run on Leonardo, so the vc-relax/NVT/NPT/CAIMD baseline remains unvalidated.

**Minimum recommended next step (single):** in a fresh, submission-authorized session, transfer the
prepared workflow (`prepared_not_submitted/`) to a new
`/leonardo/.../Diamond/diamond8_vcrelax_nvt_npt_caimd_chain_<UTC>` root and submit the four jobs as the
intended `afterok` chain (`sbatch --parsable` with `--dependency=afterok:` between stages). Then, in a
later session, re-run this audit against the produced outputs.

## 5. Safety confirmation

No job was submitted, rerun, cancelled, or resubmitted; no dependencies were altered. No QE source,
executable, or remote raw output was modified; no build, commit, or push occurred. Only local audit
files under this mirror were written. The SinglePhase workflow was not touched.
