# QE 7.5 Local Modification Audit

Audit date: 2026-08-10 (Europe/Rome)

## 1. Repository

| Item | Verified value |
| --- | --- |
| Local path | `/Users/mengqicheng/Downloads/Learning/qe-7.5` |
| Repository root | `/Users/mengqicheng/Downloads/Learning/qe-7.5` |
| QE version | Quantum ESPRESSO 7.5 (`CMakeLists.txt` and `include/qe_version.h`) |
| Local baseline | `a11913ac6cd37f39cf31056e8e4823ec7a5de49c` (`Original QE 7.5 source baseline`) |
| Official QE tag | `qe-7.5`, official GitLab commit `770a0b2d12928a67048e2f3da8d10d057e52179e` |
| Current branch at audit start | `main` |
| HEAD at audit start | `af7ebfc757c32e0b4f050ed9f2208e01020b2ee3` (`Support Nose thermostat with active CAIMD`) |
| Initial remotes | None |
| Initial tags | None |
| Initial remote branches | None |

The official references used were:

- <https://gitlab.com/QEF/q-e/-/tags/qe-7.5>
- <https://www.quantum-espresso.org/download/software/qe-7.5-ReleasePack.tar.gz>

The downloaded official ReleasePack had SHA-256
`f046f924ea16830197c2dc0ef46e8ec6d015a1d3f8105e8d6430d2272cd33554`.
A file-by-file comparison proved that every common file in the local baseline
matches that ReleasePack.  The only baseline difference is that the local tree
omits the ReleasePack's `KCW/examples/.gitignore` and the `example01`,
`example01_sym`, `example02`, `example03`, `example04`, `example05`,
`example05.1`, `example06`, and `test_io` KCW example directories.  Therefore
`a11913a...` is the correct project-local comparison base, while
`770a0b2...` is the upstream tag commit; they are not the same Git tree.

Initial working-tree state was six modified tracked files, no staged files,
and 17 untracked status entries containing source/test material and generated
results.  No reset, clean, rebase, branch deletion, tag deletion, history
rewrite, or force operation was performed.

## 2. Git Branches

| Branch | Latest commit at audit start | Unique commits | Status | Notes |
| --- | --- | ---: | --- | --- |
| `main` | `af7ebfc` | 21 after `a11913a` | Active development | Only local branch; linear history; not previously published |

There were no local branches other than `main`, no remote-tracking branches,
and no Git tags.  Consequently there are no unmerged branch-only commits to
recover and no additional development branch to publish.  Two
`refs/codex/turn-diffs/checkpoints/...` references point to a tree snapshot,
not commits or branches; they are local recovery metadata and are not treated
as branches for GitHub synchronization.  A stale prunable temporary worktree
registration at `/private/tmp/qe_velocity_baseline_src` was observed but not
removed because it contains no unique commit and cleanup was not required to
preserve the project.

The 21 project commits after the local QE baseline form one coherent CAIMD
development line: design, first strain-and-hold implementation, macOS and
Leonardo validation, parser and boundary hardening, the persistent harness,
the Ice X framework, explicit `vc-md` velocity support, and active-CAIMD Nose
thermostat support.

## 3. Modified Source Code

The complete QE core-source difference from `a11913a...`, including the
working-tree linear-deformation-gradient extension, is limited to the ten
files below.  No custom change was found in `UtilXlib`, `FFTXlib`, `LAXlib`,
`KS_Solvers`, `upflib`, `XClib`, `forces.f90`, `stress.f90`, or `pwcom.f90`.

| File | Module / routine | Main modification | Purpose / class |
| --- | --- | --- | --- |
| `Modules/input_parameters.f90` | shared PW/CP input parameters and `&CELL` namelist | Declares `caimd`, legacy strain fields, `caimd_loading_mode`, and `caimd_step_defgrad(3,3)` | New feature; MD/constraint/I/O |
| `Modules/read_namelists.f90` | CELL defaults, reads, broadcasts, validation, ignored-CELL probe | Defaults and broadcasts CAIMD fields; validates modes and finite values; rejects CP CAIMD; safely probes ignored `&CELL` blocks | New feature and parser bug hardening |
| `PW/CMakeLists.txt` | PW source list | Adds `src/caimd_module.f90` | Build integration |
| `PW/src/Makefile` | PW object list | Adds `caimd_module.o` | Build integration |
| `PW/src/caimd_module.f90` | new `caimd_module` | Owns CAIMD state, validation, strain projection, finite-deformation event maps, history/velocity transforms, diagnostics, and `.caimd`/`.caimd_defgrad` output | New numerical method; MD constraint; I/O |
| `PW/src/input.f90` | `control_iosys`, `iosys_end`, stale-history checks | Configures CAIMD for PW `vc-md`/PR, rejects unsupported restart/cell/atom/thermostat boundaries, and treats `linear_defgrad` as an active controller even without a legacy mask | New feature; validation/bug fix |
| `PW/src/run_pwscf.f90` | PW ionic-step driver | Writes pre-move CAIMD stress/strain diagnostics after electronic force/stress evaluation | I/O/diagnostics |
| `PW/src/vcsmd.f90` | `vcsmd` bridge | Preserves supported explicit input velocities, rejects unsupported history, and records final ionic-temperature diagnostics | MD bug fix and validation |
| `PW/src/vcsubs.f90` | `vcinit`, `vcmove`, projected-energy helpers | Inserts projection into the PR Beeman path; transforms atomic/cell histories; refreshes projected energy and Nose kinetic inputs | Numerical method; MD constraint |
| `PW/Doc/INPUT_PW.def` | PW input reference | Documents all CAIMD inputs, mode semantics, matrix indexing, units, event timing, and loading-axis convention | Documentation/I/O |

### CAIMD V1 (`legacy_strain`)

The committed implementation is a PW-only constrained AIMD mode for
`calculation='vc-md'`, `cell_dynamics='pr'`, `cell_dofree='all'`, and
`restart_mode='from_scratch'`.  It projects selected components in the stable
six-component order `(xx, yy, zz, yz, xz, xy)`, uses tensor shear, preserves
fractional atomic positions and peculiar Cartesian velocities across imposed
events, and leaves ordinary QE force and stress evaluation unchanged.  It is
designed as a runtime no-op when CAIMD is disabled or when the legacy mask is
all false.

The thermostat hardening in `af7ebfc` permits the existing PW `vc-md` Nose
path with active CAIMD and refreshes energy variables after the final
projection.  Other ionic thermostats and all cell thermostats remain rejected
for active CAIMD.

### Linear deformation-gradient phase 1 (`linear_defgrad`)

The initial uncommitted source work adds a second loading definition:

```text
F_load(m)  = I + (m * caimd_stride) * caimd_step_defgrad
A_event(m) = F_load(m) * inverse(F_load(m-1))
H_new      = A_event(m) * H_PR_candidate
```

The incremental map is left-multiplied onto the current PR candidate cell so
ordinary PR evolution between events is preserved.  Cell velocity and
acceleration histories receive the same affine map at events.  Zero entries in
`caimd_step_defgrad` do not freeze unrelated PR degrees of freedom.  The
diagnostic Green-Lagrange tensor is written to `prefix.caimd_defgrad`; it does
not drive the cell.

Input examples for local-z compression and local-x/local-z shear, a
deterministic math test, a six-case local `pw.x` smoke test, and Leonardo
campaign-preparation tools accompany this implementation.  The code is a
phase-1 implementation: it has current compile and small-runtime evidence but
no long finite-deformation scientific validation.

### Custom-method search and code-quality findings

Searches across the baseline-to-current diff and all selected untracked source
did not find a custom LJ7, standalone collective-variable/CV, `mu`, bias, or
custom force implementation.  Matches for generic LJ/force/constraint terms
outside the ten files above are ordinary QE 7.5 source.  The only occurrence of
LJ7 in the research archives is an unrelated audit note, not code in this
repository.

No temporary debug print was found in the modified source.  The
`caimd_print_configuration` and CAIMD output writes are intentional,
documented runtime diagnostics.  No alternate live source implementation was
found.  Saved `.patch`/`.diff` files are recovery/evidence snapshots of the
same implementation, not competing compiled versions.  Two unused imports in
`PW/src/input.f90` introduced by the phase-1 work were removed during this
audit; this is a nonfunctional cleanup confirmed by an incremental rebuild and
repeat regression run.

## 4. Tests

Status values below are deliberately limited to `PASS`, `FAIL`, `UNKNOWN`,
and `NOT RUN`.  A `JOB DONE` marker alone was never treated as scientific
success.

| Test | Input / evidence | Code tested | Result | Status |
| --- | --- | --- | --- | --- |
| Current clean local CMake build | Current tree; `/private/tmp/qe_caimd_audit_build.yV9kyQ` | All current PW/CAIMD Fortran | Serial `qe_pw_exe` compiled and linked; `pw.x` SHA-256 `99bca404807f5c46bebf4dd9c9b941149324c2b6753482c24dd35f574e66377c` | PASS |
| Current deterministic finite-deformation math | `test-suite/caimd/caimd_linear_defgrad_math.py` | Event indexing, matrix multiplication orientation, compression/shear formulae, source formula | `CAIMD linear_defgrad deterministic math tests passed` | PASS |
| Current `linear_defgrad` runtime smoke | `caimd_linear_defgrad_phase1/local_smoke_linear_defgrad.py`; generated two-atom inputs | New parser, event maps, output schema, determinant rejection | D13, D33, identity, 12-step compression, 12-step shear passed; invalid determinant rejected as expected | PASS |
| Current CAIMD regression harness | `test-suite/caimd/caimd_validation.py`; current rebuilt `pw.x` | Parser, no-op, active projection, Nose, histories, output, rejection boundaries | 42 cases executed, 6 skipped, zero failures | PASS |
| Current baseline and CP subcases | Same harness | Baseline comparison and CP rejection | Four baseline cases and two CP cases skipped because separate baseline `pw.x`/`cp.x` executables were not built in this audit | NOT RUN |
| Historical local CAIMD V1 closure | `CAIMD_V1_IMPLEMENTATION_REPORT.md` | Initial CAIMD V1 | CMake builds and 16 small runtime cases passed; traditional make failed before CAIMD compilation in legacy LAXlib preprocessing | PASS |
| Historical Leonardo CAIMD hardening | `CAIMD_LEONARDO_GPU_VALIDATION_REPORT.md`, stored summaries | GPU/MPI CAIMD V1, policies, velocity initialization | Final 34-case, 29-case, and 40-case harness checkpoints had zero failures; Slurm jobs `49543531`, `49723755`, `49728066`, and `49725763` completed. Earlier attempt `49542425` failed and was superseded | PASS |
| Ice X 100 GPa / 300 K framework | `examples/iceX_100GPa_300K_caimd` | Structure, convergence, pilot CAIMD, velocity initialization | Static/convergence/pilot and velocity checks have completed evidence; final tension/compression/shear production jobs were not submitted | PASS |
| Ice X final deformation production | Prepared deformation inputs and Slurm scripts | Long scientific tension/compression/shear runs | Explicitly not submitted | NOT RUN |
| Diamond-8 continuation to 0.30 | `examples/caimd_diamond8_local_nose_continuation_0p30_20260723T102331Z` | Long active CAIMD + Nose continuation | Equilibration completed, continuous loading stopped after 264 blocks with a singular projected cell; segmented continuation was not run | FAIL |
| Diamond-8 controls A/B/C | `examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z` | Ordinary `vc-md` Nose baseline and controls | Control A failed after 664 steps; controls B/C were not run. A targeted small CAIMD+Nose harness did pass | FAIL |
| Diamond-8 four-DOF campaign | `examples/caimd_diamond8_local_nose_four_dof_tests_0p03_20260723T115714Z` | Zero-load and constrained-loading temperature behavior | Common equilibration passed; tests 1/2 failed high-temperature safety criteria; tests 3/4 were not run | FAIL |
| Short no-velocity initial-vs-Nose isolation | `examples/diamond8_af7ebfc_no_velocity_initial_vs_nose_completed_20260727T081614Z` | Short thermostat isolation | Both 1000-step runs completed and were bounded in the short window; later long tests supersede any heating conclusion | PASS |
| Native Nose / active zero-strain 5000-step test | `examples/diamond8_native_nose_caimd_zero_strain_5000_completed_20260727T193514Z` | Long thermostat response with CAIMD false/no-op/active zero strain | All three runs reached `JOB DONE`, but temperatures were scientifically inconsistent with the 300 K target (roughly 4900 K in free PR and 523 K active zero strain) | FAIL |
| QE 7.4 vs QE 7.5 CAIMD-false, no velocity | `examples/diamond8_no_velocity_qe74_vs_qe75_caimd_false_completed_20260724T113644Z` | Ordinary/no-op version comparison | Both 1000-step runs completed and were bit-identical through the first 701 steps; low-cutoff limitation recorded | PASS |
| Ordinary Nose QE 7.4 vs `af7ebfc` | `examples/diamond8_ordinary_nose_qe74_vs_af7ebfc_completed_20260723T185307Z` | Cross-version thermostat comparison | QE 7.4 rejected the `nose` input; QE 7.5 completed but heated above 8000 K; comparison objective unresolved | FAIL |
| Older ordinary-Nose placeholder | `examples/diamond8_ordinary_nose_qe74_vs_qe75_af7ebfc_20260723T175625Z` | Prepared comparison scripts | Report remains `PENDING`; no valid completed result in this directory | UNKNOWN |
| Ambiguous chain audit | `examples/diamond8_vcrelax_nvt_npt_caimd_chain_AUDIT_BLOCKED_20260728T064656Z` | Intended vc-relax/NVT/NPT/CAIMD chain | Original workflow identity and scheduler evidence were insufficient | UNKNOWN |
| Prepared vc-relax/NVT/NPT/CAIMD chain | `examples/diamond8_vcrelax_nvt_npt_caimd_chain_completed_20260728T070746Z` | Four-stage workflow | Scripts and inputs prepared; no submission | NOT RUN |
| Submitted chain with failed vc-relax | `examples/diamond8_vcrelax_nvt_npt_caimd_chain_completed_20260728T123746Z` | Four-stage workflow | Job `50558068` failed during vc-relax input/initialization; dependent jobs `50558069`-`50558071` never ran | FAIL |
| QE 7.5 vs QE 7.4 vc-relax | `examples/diamond8_vcrelax_qe75_vs_qe74_completed_20260728T234619Z` | Ordinary vc-relax baseline | Jobs `50605866` and `50605867` completed, converged, and agreed | PASS |

No Leonardo connection, transfer, build, allocation, submission, monitoring,
restart, cancellation, or result modification was performed during this audit.
No expensive production calculation was started.

## 5. Analysis / Python Tools

All relevant tools were retained.  The main groups are:

- `test-suite/caimd/caimd_validation.py`: persistent generated-input software
  harness for local and Leonardo validation.
- `test-suite/caimd/caimd_linear_defgrad_math.py`: deterministic source/math
  checks for the new loading mode.
- `caimd_linear_defgrad_phase1/local_smoke_linear_defgrad.py`,
  `leonardo_campaign_tools.py`, and
  `setup_leonardo_linear_defgrad_campaign.py`: runtime smoke and future
  campaign preparation.
- `examples/iceX_100GPa_300K_caimd/scripts` and `tools`: input preparation,
  velocity generation, convergence/strain analysis, validation, and Slurm
  script generation.
- `examples/caimd_diamond8_*_20260723*/scripts`: continuation, terminal-state
  extraction, segmented validation, controls A/B/C, and four-DOF workflows.
- `examples/diamond8_*_2026072*/scripts` and `analysis_audited`: input
  generation, pre-submit gates, output parsers, audit recomputation, plotting,
  vc-relax extraction, and report generation.
- `examples/caimd_development_progress_html_20260723T205507Z/scripts`:
  reproducible metric recomputation and figure generation for the self-contained
  HTML progress report.

All 41 Python files under the relevant CAIMD/test archives passed
`python3 -m py_compile` in this audit.  Shell and Slurm scripts were preserved
as research workflow records; they were not submitted.

## 6. Generated / Temporary Files

Before ignore rules, the untracked research/result material contained roughly
738.5 MiB of files.  The new ignore policy leaves about 29.3 MiB (943 files
before adding this report) as the candidate source/evidence set and excludes
about 709.2 MiB (1460 files) of reproducible runtime products.

Kept for Git:

- all modified QE source and input documentation;
- all CAIMD test inputs, fragments, Python/shell/Slurm scripts, reports, JSON,
  CSV, NPZ, HTML/SVG/PNG analyses, hashes, patches, and compact harness output;
- expected-rejection `CRASH` files and small named harness `.out` files;
- compact scientific comparisons needed to understand both passes and fails.

Left local and ignored, without deleting or moving anything:

- Python bytecode/cache files;
- root build products and machine-specific `make.inc`;
- timestamped campaign `tmp`/`.save` trees, charge density, wavefunctions, and
  XML restart state;
- large `qe.out`, `qe.err`, Slurm stdout/stderr, duplicate `*.UPF` files, and
  raw `.ave`, `.avec`, `.caimd`, `.cel`, `.e`, `.eal`, `.for`, `.mdtrj`, `.p`,
  `.pos`, `.tv`, and `.vel` streams in the bulky timestamped campaigns.

The excluded files remain at their original local paths.  No scientific result,
trajectory, restart, input, source file, or unique analysis was deleted.

## 7. Repository Cleanup

A root `.gitignore` was added with narrow rules for operating-system/editor
artifacts, Python caches, local QE build products, and the bulky timestamped
campaign runtime files listed above.  The rules deliberately keep compact
analysis and failure evidence.  The largest newly selected file is about
1.3 MiB; no selected file approaches GitHub's 100 MiB object limit.

The complete selected source and untracked candidate set was scanned for
private-key headers, GitHub/PAT token patterns, AWS access keys, passwords,
API keys, access tokens, and secret keys.  No match was found.  No credentials,
SSH private keys, or token files are selected for Git.

## 8. GitHub Sync

| Item | Verified state after the content push |
| --- | --- |
| Repository | `mackeysjtu-ops/qe-7.5` |
| URL | <https://github.com/mackeysjtu-ops/qe-7.5> |
| Visibility | **PRIVATE**, directly verified on the repository page before and after push |
| Remote | `origin` -> `git@github.com:mackeysjtu-ops/qe-7.5.git` |
| SSH identity | Verified as `mackeysjtu-ops` |
| Pushed branch | `main`; it tracks `origin/main` |
| Synchronized content commit | `d9dad0437b0a739fe06a1b721d3ecd262151cfcb` (`Preserve CAIMD linear deformation and validation archive`) |
| Other branches | None; no other local branch exists |
| Tags | None |
| Force operations | None permitted or used |

The existing `gh` credential for `mackeysjtu-ops` was invalid.  A GitHub CLI
device reauthorization was attempted, but the page requested an additional
Workflow permission and kept its authorization action disabled; the request
was cancelled rather than bypassed.  Repository ownership and PRIVATE
visibility were therefore verified in the existing signed-in GitHub session,
and authenticated SSH is used for Git transport.

The ordinary, non-force push created `main` on GitHub and set the local branch
to track `origin/main`.  The refreshed repository page showed `Private`, one
branch, zero tags, and content commit `d9dad0437b0a739fe06a1b721d3ecd262151cfcb`.
This report-state update is a subsequent metadata-only commit; its exact
`main`/`origin/main` SHA-1 is recorded in the final task response because a Git
commit cannot embed its own SHA-1 without changing that SHA-1.

## 9. Remaining Issues

1. `linear_defgrad` has compile, deterministic-math, and short-runtime evidence,
   but no long finite-deformation scientific validation, no 64-atom/3000 K
   validation, and no Leonardo GPU/MPI validation.
2. CAIMD restart remains intentionally unsupported because protocol state is
   not persisted.
3. The long native-Nose and active-zero-strain studies expose unresolved
   temperature-control behavior.  A successful short harness is software
   closure, not proof of 300 K scientific control.
4. The 0.30 continuation reached a singular-cell safety stop; controls A/B/C
   and the four-DOF campaign remain incomplete or failed.
5. Final Ice X tension/compression/shear production calculations were never
   submitted.
6. Six current harness cases were not run because separate baseline `pw.x` and
   `cp.x` executables were not built.  Historical local/Leonardo reports cover
   these boundaries, but they were not repeated against this exact dirty-tree
   executable.
7. The `gh` CLI credential remains invalid; GitHub web verification plus SSH
   transport are used for this synchronization.
8. Bulky raw runtime results remain local and intentionally are not mirrored to
   GitHub.  Their compact analyses, scripts, inputs, and reports are preserved
   in Git.

## Validation Commands Executed In This Audit

Key local commands, in addition to Git inventory/diff/log/branch scans, were:

```bash
diff -qr <export-of-a11913a> <official-qe-7.5-ReleasePack>
git diff --check
python3 -m py_compile <41 relevant Python files>
python3 test-suite/caimd/caimd_linear_defgrad_math.py --source-root .
cmake -S . -B /private/tmp/qe_caimd_audit_build.yV9kyQ \
  -DQE_ENABLE_MPI=OFF -DQE_ENABLE_OPENMP=OFF -DQE_ENABLE_TEST=OFF
cmake --build /private/tmp/qe_caimd_audit_build.yV9kyQ \
  --target qe_pw_exe -j 4
python3 caimd_linear_defgrad_phase1/local_smoke_linear_defgrad.py \
  /private/tmp/qe_caimd_audit_build.yV9kyQ/bin/pw.x
python3 test-suite/caimd/caimd_validation.py --mode local \
  --current-exe /private/tmp/qe_caimd_audit_build.yV9kyQ/bin/pw.x \
  --pseudo pseudo/C.UPF \
  --work-root /private/tmp/qe_caimd_audit_validation \
  --run-id audit-final-20260810 --timeout 180
ssh -o BatchMode=yes -o ConnectTimeout=15 -T git@github.com
```

The local harness summary is at
`/private/tmp/qe_caimd_audit_validation/audit-final-20260810/analysis/summary.json`.
It records 42 executed cases, 6 skips, and an empty failure list.  The current
smoke summary is preserved at
`caimd_linear_defgrad_phase1/evidence/local_smoke_summary.json`.
