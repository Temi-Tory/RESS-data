# Corrections applied to the manuscript draft (2026-09-06)

Base file: `validation/probability/newress.zip` → `main.tex` (the current, applied revision —
1568 lines, compiles to 43pp / 978,773 bytes unmodified). Working copy edited in the scratchpad,
recompiled after every change (`pdflatex`, no errors). **The corrected file is not yet copied
back into `newress.zip` / the working repo — it sits in this `manuscript/` folder pending your
review.** After edits below: 44pp / 979,920 bytes.

## Applied

1. **Dropped the "27–28" drone maxcond figure (2 occurrences: §5.4.3, §5.5).** This number was
   never a real measurement — the identification run demanded ~25.5GB on a 16GB machine and never
   completed (`validation/probability/notes/CORPUS_CAMPAIGN_HANDBACK.md` §G;
   `validation/fresh_20260816/MASTER_FINDINGS.md` corrections ledger #1). Replaced both instances
   with the honest measured-memory statement (identification itself exceeded the memory budget).
2. **Added the dPrPm lower/upper bound columns to Table `tab:grid_accuracy`**, resolving the
   in-manuscript author placeholder `[Author note: carry the dPrPm bound columns over from the
   original submission's table.]`. Values read directly from Tong & Tien (2019), Table 2 (page 7
   of the source PDF, `.../Examples of other ppl written works/tong-tien-2019-....pdf`) — not
   retyped from memory. Added one sentence noting the sole nonzero dPrPm gap (node 16, 0.00048)
   coincides with the source study's own hardest-to-bound node, and that IPA carries no gap there.
3. **Grid runtime table (`tab:grid_runtime`): 2.378ms/4.94MB → 0.928ms/1.72MB**, per
   `MASTER_FINDINGS.md`'s "REFRESHED" line (#2 in its confirmation table). Added a caption caveat
   that these are machine-specific figures, not for cross-machine comparison.
4. **Interval-vs-BDD one-shot timing (§5.3.1 prose + `tab:interval_timing`): "2.9 to 95" →
   "3.9 to 99.9"**, and every per-family cell in the table updated to the fresh values in
   `validation/probability/data/interval_bdd_vs_ipa_timing.csv` (this is the current, Aug-30-dated
   authoritative copy of that CSV — confirmed by direct inspection, not assumed).
5. **`tab:structured`: grid unique-diamond count 41→39, KarlNetwork 147→145** — confirmed cosmetic
   drift per `MASTER_FINDINGS.md` #4; grid's 39 also independently confirmed against
   `validation/probability/data/complexity_validation.csv` directly (Karl's 145 taken from the
   ledger, not independently re-derived this session — flag if you want it re-measured).
6. **`tab:drone`: vtol-dense-decentralized `|C|_max` 17→16.** Confirmed fresh 2026-08-30
   (`validation/probability/drone_ksweep_remeasure_results.csv`, per `INDEX.md` §3 item 2) — a
   real, expected change from the `is_det` generalisation fix, not noise. Updated the two
   downstream prose mentions of "15–17" (§5.4.3, §5.5) to "16" (both dense configurations now
   measure identically).
7. **§4.3 cost-formula wording**: added "is therefore a sound a-priori **upper bound**... It is an
   upper bound rather than the realised cost because memoisation... reduces the work actually
   performed... observed to be a small fraction of $W$ in practice" — per `MASTER_FINDINGS.md`
   correction #2 (the unqualified "computable in advance"/"definite cost" framing overclaims;
   realised work is measured at 1–5% of the upper bound in practice).
8. **§5 preamble**: added the shared BDD/IPA measurement-environment sentence (same process, same
   machine, dynamic reordering unless stated) and the degeneracy-preserving synthetic-widening
   methods sentence — per `MASTER_FINDINGS.md` correction #9 and this session's
   `reviewer_response/JUDGMENT_CALLS.md` #2.

9. **Data and Software Availability section**: split the single package-only link into a code
   statement (package, General registry + GitHub, version pinned) and a data statement (new
   Zenodo deposit, DOI placeholder, MIT licence, README-mapped) — reflecting your decision to
   mint a separate paper-scoped deposit rather than pointing reviewers only at the package repo.
10. **RETRACTED then RE-APPLIED with clean data, same session (2026-09-06, later still).**
    A proper re-measurement was done: two scripts (`task2_pbox_steps_v2_one.jl`,
    `pbox_cost_vs_diamonds_v2_one.jl`) rewritten so each Julia process performs exactly ONE timed
    propagation on the real target, with JIT warm-up done on a tiny throwaway 5-node network
    instead of the target itself. Both validated on a single run before being trusted for a full
    batch. **Final clean steps-scaling curve**: $1.5424$s / $8.6771$s / $51.0163$s / $339.0249$s at
    steps 25/50/100/200 (`data/task2_pbox_steps_v2/`), exponent $\approx2.6$ overall, with the
    *local* exponent increasing across the range (2.49, 2.56, 2.73 between successive points) —
    reinstated into §5.3.2 and §5.5. As a bonus, the companion cost-vs-diamond-count experiment
    (`data/pbox_cost_vs_diamonds_v2/`, all 14 networks) found a genuinely clean result: propagation
    time at fixed steps=50 tracks `measured_ops` (the realised, memoisation-adjusted work count
    already used for the exact/interval cost story), not the raw `sum_2^C` bound — ratio
    time/`measured_ops` stayed within roughly 0.14–0.26 s/op across an 89× range in `measured_ops`
    (`bridge_5`, the single highest-`sum_2^C` network in the batch at 4072, finished faster than
    two lower-`sum_2^C`-but-higher-`measured_ops` networks — the cleanest possible demonstration
    that `sum_2^C` is not the driver). This directly strengthens `manuscript/PBOX_COST_MECHANISM_DRAFT.md`
    (updated) — recommend folding this into the §4.4 addition as a headline point, not a footnote.
    Originally applied (then retracted, now reconfirmed): p-box steps-scaling
    curve (§5.3.2 prose + §5.5 limitations paragraph) "grows quadratically... 2.7s/8.3s/110s at
    steps 50/200/800" → "grows super-linearly (empirical exponent ≈2.8)... 1.4s/7.0s/45.8s/300.8s
    at steps 25/50/100/200", backed by `pre-write final/data/task2_pbox_steps/`. **This has been
    reverted** (both occurrences now read "grows super-linearly with the discretisation level
    [NUMBER PENDING RE-VERIFICATION]", with an inline `%% PENDING RE-VERIFICATION` comment).
    Reason: `task2_pbox_steps_confirm.jl` runs a warm-up plus eight timed propagations, of
    increasing size, **all in one Julia process** — exactly the pattern already documented
    elsewhere in this project as producing polluted timings (`MASTER_FINDINGS.md`: "float corner
    legs ran after the interval leg's multi-GB cache in the same process: 2,291s vs 81s
    fresh-process reference... **process-hygiene rule, now twice learned: one timing measurement
    per fresh process**"). A live follow-up experiment (p-box cost vs. diamond count across other
    networks) surfaced the same pattern directly — a warm-up and its immediately-following timed
    call showing timings inconsistent with a clean JIT-then-fast story — which is what prompted
    re-examining Task 2's own script and finding it has the identical structural flaw. **The old
    "quadratic, 2.7/8.3/110s" claim was already known wrong (a different bug); this session's
    replacement number cannot currently be trusted either, for a different reason, and neither
    should go in the manuscript until a proper one-measurement-per-fresh-process re-run is done.**
    The qualitative claim ("p-box cost grows faster than linearly with discretisation level, and
    is a controllable, not-free, dial") is very likely still correct — only the precise exponent
    and absolute times are in doubt. **Note for a follow-up pass**: the same stale "quadratic"
    claim also still appears (not manuscript-facing) in
    `validation/certified_bound_threshold_sweep_steps200.jl`, `validation/rc_pbox_cvx.jl`, and
    `validation/probability/notes/PBOX_DILEMMA_SUMMARY.md`.

## Deliberately NOT applied here — flagged as open decisions, not silently done

- **Adversarial (fanin-k/mesh-w) table/figure** — confirmed absent from the manuscript entirely
  (no "adversarial"/"fanin"/"mesh" string anywhere in `main.tex`), despite fully validated fresh
  data existing. This is a genuine scope decision (add a new subsection vs. leave it out), not a
  correction — see `../README.md`.
- **§5.4.1/§5.4.3 additional disclosure/rigour-narrative sentences** recommended in
  `reviewer_response/JUDGMENT_CALLS.md` #1 and #2 (the "proxies, here is exactly how" framing and
  the "we caught our own warm-up mistake" sentence) — these are voice/rhetorical choices for you
  to accept or reject, not inserted unilaterally.
- **§4.3/§4.1 credal-network citations** (`JUDGMENT_CALLS.md` #3) — needs two new bibliography
  entries first; not inserted pending your sign-off on the recommendation.
- Editorial/proofreading pass (reviewer comment 8, first table) — not attempted; flagged as still
  owed in the tracker.

## Verification performed
- Compiled the untouched original (`pdflatex`, from a fresh unzip) — 43pp, byte-identical to the
  zip's own `main.pdf` (978,773 bytes) — confirms the baseline is exactly what `newress.zip` ships.
- Compiled the edited copy after every change batch — zero `pdflatex` errors throughout; final
  version 44pp, 979,920 bytes (net +1 page from the expanded grid-accuracy table and the added
  prose, nothing else moved).
- Every numeric change above is cited to a specific artifact file (named in each bullet), per the
  project's "fresh-or-it-didn't-happen" rule — none are hand-typed from memory or from a prior
  session's prose alone.
