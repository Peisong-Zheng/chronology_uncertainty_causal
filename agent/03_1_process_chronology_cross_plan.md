# Simulation realization by Z age-scheme cross: compact diagnostic

Authorized on 10 November 2026. The root script is `03_1_check_process_chronology_cross.py`. This branch uses the completed 03 experiment as a read-only source. Its question is whether the observed X-Y decision changes follow the simulated value history, the assigned Z age scheme, or a particular combination. The outcome is the X-Y relation; recovering the entire three-variable network is outside this branch's scope.

## What is crossed

A **simulation realization** is one generated X,Y,Z history under the fixed linear equations. An **age scheme** is one saved uncertain observed chronology for Z together with its ten ordered chronology members. The members were drawn around that observed chronology. Reassignment uses the complete saved age arrays, without generating new processes, observed centers or chronology members.

The two systems remain **Common driver only** (Z to X and Z to Y) and **Common driver + direct X to Y**. Memories remain Z/X/Y = 500/200/200 years, common-driver gains two, direct gain zero/one, and every direct equation delay 20 years. Sampling positions, innovations, measurement noise, analysis grid, estimator settings and chronology model are unchanged. The age schemes retain C-type 40-year proposal SD and the original correlation length; realized error is not bounded by that SD.

Realizations 1, 2 and 3 use the existing seeds 20261009, 20261301 and 20261302. Age schemes 1, 2 and 3 are the Z chronology packages originally associated with those same seeds. All source realizations have the same sampling positions and true sample ages, making the pointwise reassignment well defined.

For realization r and scheme a, keep the X/Y analysis inputs from realization r's irregular accurate-age reference. Assign the ages from scheme a to the sampled Z values from realization r, then interpolate Z onto the unchanged 29-1 kyr grid of 1401 points. Member m always denotes the same saved chronology curve within a scheme. The fixed observation mask and source true-age arrays are checked, not assumed. Whole curves retain their original ordering and coverage; no sorting, clipping or extrapolation is allowed.

| Simulation realization | Z age scheme 1 | Z age scheme 2 | Z age scheme 3 |
|---|---|---|---|
| Realization 1 | Reuse completed 03 | New conditional analysis | New conditional analysis |
| Realization 2 | New conditional analysis | Reuse completed 03 | New conditional analysis |
| Realization 3 | New conditional analysis | New conditional analysis | Reuse completed 03 |

This matrix is run for both systems. Only Z has uncertain ages in the cross. The previously completed X/Y-only and all-three conditions remain in 03; they are not expanded here. Both original-uniform and irregular accurate-age references are retained, with the latter used for comparisons in figures and paired changes.

## Scale and reuse

- Crossed uncertain inputs: 2 systems x 3 realizations x 3 Z schemes x 10 members = **180**.
- Accurate-age references: 2 systems x 3 realizations x 2 references = **12**.
- Total branch cases: **192**. The 60 diagonal cases and 12 references reuse completed original fits in full.
- New off-diagonal cases: **120**, each requiring three conditional calls (PCMCI+, LKIF and GC), giving **360 new method calls**.
- X/Y inputs and estimator seeds are unchanged by the Z assignment. All four pairwise methods, including CCM, reuse each realization/system's irregular accurate-age evidence. They are not recalculated to create duplicate evidence.
- Every case retains seven complete native archives and 26 directional rows, giving **1344 archives and 4992 rows**. Of the archives, 984 are exact copies and 360 are newly computed. Source failures and unresolved evidence are preserved. A 26-row edge table is the final case commit.

The saved estimator seed stays fixed within each realization, matching 03. Thus a row pattern describes a realization under its original estimator seed; it is not a pure variance-component estimate for dynamics alone. Every age-scheme comparison within that row holds the estimator seed fixed.

## Scientific outputs

The primary outcome is exact X-Y recovery: neither X/Y arrow for Common driver only, and only X to Y for Common driver + direct X to Y, with no unresolved evidence. All requested cases remain in the denominator. Unresolved, failed and pending outcomes remain distinguishable from definite wrong-arrow states. Directional retention, original native scores, within-cell support with usable/requested counts, and score changes relative to the same realization's irregular accurate reference are saved.

Native three-variable fits necessarily retain Z-related fields, but no full-network recovery metric or Z-edge performance sweep is produced. Existing 03 findings are not replaced or relabeled. The source case ID, realization seed, age-scheme seed, member index and reuse origin are explicit in every case.

## Figures

1. **00_cross_design:** distinguish simulated value histories from age assignments; show three accurate-age Z histories over the preselected 17-13 kyr window, the reuse/new 3x3 matrix, and the actual age-center errors with the range of the ten saved members. Ages increase to the left. The member range is descriptive, not a confidence interval.
2. **01_xy_recovery_matrix:** two system rows by three conditional-method columns. Each panel is a 3x3 matrix: simulation realizations are rows, Z age schemes are columns. Cells show correct/10 and unresolved counts, with unavailable counts if needed. The original diagonal is outlined. Row labels give that realization/method's accurate-age correct/1 result. All panels use the same 0-100% color scale. Pairwise evidence is invariant to the Z assignment and remains in tables rather than being repeated as additional curves.

Use descriptive system names, unique panel labels and direct count annotations. Captions explain both factors and the requested denominator. A low recovery fraction can arise from different state changes; those remain available in the saved X-Y state table. No diagram suggests nine independent histories.

## Validation and execution

Read the source generated series with round-trip float parsing, and use the same ordinary CSV parsing as 03 at estimator input. Independently reconstruct the diagonal Z inputs and require equality with the saved 03 inputs before copying their fits. Check the full arrays actually read by estimators when reusing pairwise evidence. Check source settings, packages, numerical platform, computation digest, source files and estimator seeds; do not edit old metadata to allow a mismatch.

The preflight in `agent/verify_process_chronology_cross.py` checks all 192 prepared identities, complete age-array reuse, unchanged X/Y, common grids, paired systems, diagonal evidence, copied native archives, failure retention and resume/input/source guards. One off-diagonal case compares all three conditional native outputs in direct and Windows-spawn calls. Serialization and numerical arrays are checked, not only significant arrows.

Use 32 outer workers with one numerical thread each, the existing native GPDC null caches, and no nested estimator pools. No new GPU method or worker benchmark is required. Run only one instance for the branch folder. Output locations are `data/processed/03_1_check_process_chronology_cross/` and `figures/03_1_check_process_chronology_cross/`. Runtime production code does not depend on `tmp/`.

After completion, audit all committed cases and native evidence, exact source copies, summary reconstruction and input guards. Review the two actual figures and captions. Confirm hashes, sizes and modification times of all pre-existing scientific files against the pre-run snapshot. Preserve complete unfavorable and unresolved outcomes rather than selectively repeating them.

## Interpretation limits

A pattern following a row suggests realization-associated sensitivity; a pattern following a column suggests age-scheme-associated sensitivity; a particular cell can indicate their combination. These are descriptive contrasts among three existing realizations and three existing age schemes. They do not constitute nine independent climate histories, a calibrated population error rate, or a formal decomposition of process and chronology variance. The cross is an exploratory diagnostic motivated by the initial 03 result; it does not provide independent confirmation of a universal method ranking.

## Completed execution

The 10 November run completed all 192 cases in 799.423 seconds using 32 workers. It performed the planned 360 new conditional calls and retained 984 exact source-native copies. All 4992 directional rows have successful fit status; unresolved scientific decisions remain explicit. Minimum sampled available physical memory was 11.26 GiB.

The preflight passed at `data/benchmarks/process_chronology_cross_checks/20261110T123801313973Z/report.json`. The complete artifact audit and PNG/PDF figure review passed in `data/logs/process_chronology_cross_20261110T124145323892Z/`. Every matrix cell and accurate-age reference was checked against the saved directional evidence. The [results report](03_1_process_chronology_cross_results.md) explains the six matrices, accurate-age exceptions and the limits of the descriptive crossing; [captions](figure_review/07_process_chronology_cross.md) accompany the two figures.

The final preservation audit confirmed unchanged hashes, sizes and modification times for all 26,464 pre-existing scientific files. No existing 02 or 03 result or figure was overwritten, and no scientific fit was repeated during auditing.
