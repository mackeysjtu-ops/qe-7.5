# Diamond8 vc-relax → NVT → NPT → CAIMD chain — audit report

**OVERALL STATUS: `BLOCKED_AMBIGUOUS_WORKFLOW`** (Phase 1 stop condition)

Audit date: 2026-07-28T06:46:56Z. Read-only. No job submitted, no rerun, no
modification of any source, executable, or raw result.

## 1. Why the audit stopped at Phase 1

The four-job dependent workflow described in the task **was not found on
Leonardo**. The blocking condition is *absence*, not ambiguity between several
candidates: zero candidate roots exist.

Verified read-only on `login01.leonardo.local` as `mcheng00`
(remote time 2026-07-28T06:44:06Z):

| Evidence required by Phase 1 | Found |
|---|---|
| `.../Diamond/diamond8_vcrelax_nvt_npt_caimd_chain_*` | **0 matches** |
| `scripts/01_vc_relax.slurm` | **0 matches** (home + scratch, maxdepth 6) |
| `scripts/02_fixed_nvt.slurm` | **0 matches** |
| `scripts/03_free_npt.slurm` | **0 matches** |
| `scripts/04_caimd_zero_nonzero.slurm` | **0 matches** |
| `DEPENDENCY_CHAIN.json` | **0 matches** |
| four stage directories | **0 matches** |
| `WORKFLOW_METADATA.json` | 3 matches, all unrelated studies |

Also searched: `/leonardo_work/ICT26_CMSP_0`, all `*chain*` directories
(1 match, an unrelated CARTACE/LJ7 run), all `*vcrelax*` / `*vc_relax*` /
`*npt_caimd*` directories (2 matches, both inside the QE source tree), every
directory modified since 2026-07-26, and the local
`qe-7.5/examples/` tree (no `*chain*` or `*vcrelax*` workflow).

Full search record and rejection rationale:
`evidence_audited/workflow_selection.json`.

## 2. What does exist from the relevant period

The newest Diamond workflow is
`/leonardo/home/userexternal/mcheng00/workspace/Diamond/diamond8_native_nose_caimd_zero_strain_5000_20260727T151734Z`
(created 2026-07-27T15:26:03Z). Its own `WORKFLOW_METADATA.json` describes it as:

> three 5000-step (10 ps) native-init Nose trajectories: CASE A `caimd=.false.`
> (control), CASE B `caimd=.true.` all masks `.false.` (no-op), CASE C
> `caimd=.true.` all masks `.true.` `step_strain=0` (zero-strain hold)

That is a **three-case single-stage study**, not a four-job dependent chain. It
has no vc-relax stage, no fixed-cell NVT stage, no free-cell NPT stage, and no
dependency metadata; its `scripts/` holds `generate_inputs.py`,
`pre_submit_validate.py`, `run_long_zero_strain.slurm`,
`analyze_long_zero_strain.py`. It was already mirrored locally as
`examples/diamond8_native_nose_caimd_zero_strain_5000_completed_20260727T193514Z`.

Phase 1 explicitly forbids selecting a directory solely because it is the
newest, and this one does not match the described design on any structural
criterion. It was therefore **rejected, not audited as a substitute**.

Interpretation (stated as inference, not fact): the completed calculations the
user has seen on Leonardo are most plausibly this Jul 27 three-case study. There
is no filesystem evidence that the vc-relax → NVT → NPT → CAIMD chain was ever
created, let alone submitted.

## 3. Slurm accounting

**Not queried.** The task authorizes exactly one grouped `sacct` query and only
*after* all four Job IDs have been recovered. Zero Job IDs were recovered, so the
precondition was never satisfied and no scheduler query was issued.

## 4. Answers to the primary questions

Questions 1–10 (vc-relax pressure, NVT temperature, NVT pressure, NPT
temperature/pressure, ionic and cell Nose activity, volume boundedness, CAIMD
zero strain, xx strain, first failed stage, larger-cell readiness) are **all
unanswerable**: they require raw output of a workflow that does not exist.
No value is reported, because any number would be fabricated.

Phases 2–14 were not executed for the same reason. The only physical statement
that can be made without data is the one that is definitional, not empirical:

* **D. Is the pressure observed in fixed-cell NVT meaningful as a controlled
  ensemble quantity?** No. With a fixed cell the barostat degrees of freedom are
  absent, so pressure is a diagnostic observable only — it is not controlled.
  This follows from the ensemble definition and holds independently of any run.

## 5. Local repository state (read-only)

```
root:   /Users/mengqicheng/Downloads/Learning/qe-7.5
branch: main
HEAD:   af7ebfc757c32e0b4f050ed9f2208e01020b2ee3   (matches expected)
```

Not modified by this audit. Executable and pseudopotential SHA-256 verification
(Phase 4) was not performed, because it is only meaningful against the runtime
evidence of the selected workflow.

## 6. Safety contract

| Action | Count |
|---|---|
| jobs submitted | 0 |
| jobs cancelled / resubmitted | 0 |
| QE reruns | 0 |
| dependencies changed | 0 |
| QE source modified / rebuilt | no / no |
| remote raw outputs modified | no |
| commits / pushes | 0 / 0 |
| waiting, sleeping, watching, polling | none |

Remote access was limited to read-only inspection (`ls`, `find`, `cat`,
`hostname`, `date`).

## 7. Minimum next step

Exactly one, and it is a *location* step, not a scientific one: **establish
whether the four-job chain was ever created and submitted**, before any audit is
attempted again.

If you authorize one scheduler query, the decisive check is:

```bash
sacct -S 2026-07-26 -u mcheng00 \
  --format=JobID,JobName%36,State,ExitCode,Elapsed,Start,End,WorkDir%80
```

* If it shows four dependent jobs with a `WorkDir` that no longer exists, the
  chain ran and its directory was moved or removed — recover from that path.
* If it shows no such jobs, the previous session stopped **before** submission,
  the chain must be prepared and submitted afresh, and the completed
  calculations the user observed are the Jul 27 three-case study.

Until that is settled, no vc-relax, NVT, NPT, or CAIMD conclusion can be drawn.
