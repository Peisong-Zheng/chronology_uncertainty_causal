# Figure captions: realization-by-age-scheme crossing

Reviewed on 10 November 2026 against the completed 192-case run of `03_1_check_process_chronology_cross.py`. Both final PNGs and rendered one-page PDFs were inspected. The figures assess the complete X–Y decision only; they do not assess recovery of the full three-variable network.

## Figure 00 — Crossed design and actual age assignments

**Files:** `00_cross_design.png` and `00_cross_design.pdf`

**Caption.** Crossed design separating variation among simulated histories from sensitivity to the assigned Z chronology. (a) A realization supplies one X, Y and Z history generated under fixed linear dynamics and intrinsic memories; an age scheme supplies one uncertain observed Z chronology and its ten ordered chronology members. (b) The three accurate-age Z inputs are shown over the fixed 17–13 kyr BP window; Z is identical between the two prescribed systems for a given realization. (c) Each realization receives every age scheme, with outlined diagonal cells reusing the original 03 results and the remaining cells requiring new conditional analyses; each cell contains ten members and the same crossing is used for both systems. (d) Lines show observed chronology errors relative to the true sample ages, and shaded bands show the pointwise minimum and maximum errors among the ten saved members, whose curves were drawn around the observed chronology; the bands are descriptive ranges, not confidence intervals, and positive errors assign samples older ages. Only Z timing changes within a realization, so X/Y inputs and pairwise results, including CCM, remain unchanged; the crossing comprises three simulated histories and three reused age schemes, not nine independent histories.

## Figure 01 — Exact X–Y recovery for every realization and age scheme

**Files:** `01_xy_recovery_matrix.png` and `01_xy_recovery_matrix.pdf`

**Caption.** Exact recovery of the prescribed X–Y relation after conditioning on Z with uncertain ages. Panels (a–c) show Common driver only, for which the correct decision retains neither X/Y direction; panels (d–f) show Common driver + direct X → Y, for which only X → Y is correct, with PCMCI+, LKIF and GC in successive columns. Within each panel, rows identify simulation realizations and columns identify complete Z age schemes; colour gives the percentage of all ten requested members with the correct complete X–Y state, the first cell label gives the exact correct/requested count, and U gives the number with unresolved evidence. Row labels marked Acc give the corresponding single irregular accurate-age decision for that realization and method, with unresolved, failed or pending reference outcomes explicitly identified where present; outlined cells reuse the original 03 pairing. A decision with unresolved evidence is never counted as correct, and failed or pending cases remain in the requested denominator, with F/P counts shown where present; saved state tables distinguish definite missed or extra arrows. Row and column patterns are descriptive comparisons among three existing histories and three existing age schemes, with the estimator seed fixed within each realization, and do not estimate independent variance components or a population error rate.

## Reading notes for the result summary

- Read across a row to compare complete Z age schemes for the same simulated history and the same estimator seed.
- Read down a column to compare the same complete age scheme across the three simulated histories, each retaining its original estimator seed.
- Compare each row with its own Acc decision before describing deterioration or improvement. Acc 0/1 is one incorrect or unresolved reference decision, not an estimated zero population success probability.
- A low correct count need not mean that the retained arrows are consistently wrong: U identifies unresolved evidence, while the saved complete states distinguish specific arrow errors.
- The original diagonal alone paired each history with only one age scheme. The off-diagonal cells reveal whether that pattern persists under reassignment; they are additional comparisons of existing histories, not new independent simulations.

## Completed figure review

The final exports are in `figures/03_1_check_process_chronology_cross/`. Figure 00 has four labelled panels; Figure 01 has six labelled scientific panels and a separate common colour bar. Text, cell annotations, reference labels, diagonal outlines and legends are readable in both PNG and PDF exports, without clipped text or material overlap.

All 54 cross-cell counts and 18 accurate-reference counts were checked against the saved state summaries and independently reconstructed from 1,116 corresponding X/Y directional rows. Every cross-cell denominator is ten and each accurate-reference denominator is one. No failed or pending outcomes occur in these displayed cases; unresolved counts remain visible and are excluded from exact recovery. The two exceptional accurate-reference decisions remain displayed: PCMCI+ is unresolved for Realization 2 in Common driver + direct X → Y, and LKIF is incorrect for Realization 1 in that system.

PDF review used the existing Windows PDF renderer without installing packages or changing the scientific environment. The audit, final export hashes and all checked cell/reference counts are recorded in `data/logs/process_chronology_cross_20261110T124145323892Z/figure_review.json`. The two branch figures were re-exported once from unchanged code and saved inputs during the canvas audit; no fits, scientific tables, root plotting source or original 03 figures were changed.
