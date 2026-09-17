# Three flagged judgment calls — recommendations for sign-off

Prior sessions drafted full responses to these three points but flagged their own doubt about
them inline (quoted in the audit). Per the user's instruction, here is a recommendation for each,
clearly marked as a recommendation, not settled text. Accept, edit, or reject each independently.

---

## 1. Drone case study — "just another generated network that fits a story"?

**The risk**: reviewer #5 (first table) already criticised the original six-Pareto case study for
being disconnected from real reliability insight; a reviewer could look at the *rebuilt* three
proxy networks and suspect they were reverse-engineered to produce a flattering result.

**Recommendation: lean into full disclosure rather than minimising the construction.** The
rebuilt case study is, on the evidence gathered, the more defensible of the two options because
every input is either (a) taken directly from the source study's own numbers/rules, or (b) one
single, explicitly flagged extension (the weather-derating interval), not several invented
parameters as in the original six-Pareto version. The honest move is to say this plainly and
early in §5.4.1, not to let a reviewer discover the proxy nature and wonder what else was
shaped to fit:

> "The source study's optimised network layouts are not published (its own data-availability
> statement: 'available upon reasonable request'). We therefore do not reproduce them; we
> construct three configurations as explicit, labelled *proxies* for the qualitative character
> of three described trade-off points, built entirely from the source study's own stated rules
> and figures, with exactly one flagged extension. We report this construction in full (Table X)
> so that the proxy nature of the topology — as distinct from the reliability-model inputs, which
> are traceable — is not left implicit."

This is close to the wording already in `RESS_edit_proposals.md` EDIT 15 §5.4.2 — the
recommendation is to make sure this disclosure survives into the final draft verbatim and is not
softened, and to say the same thing again, briefly, in the reviewer response letter itself
(comment 5 and R2 comment 6 responses) rather than only in the manuscript. A reviewer who reads
"we built proxies and here is exactly how, because the real ones aren't public" is reassured, not
suspicious; a reviewer who has to infer it is not.

**One more thing worth doing** (small, not yet in the plan's task list — optional): add one
sentence noting that the *qualitative* redundancy finding (K=16 is where sifted/naive BDD both
fail while IPA still completes) does not depend on the proxy construction being exactly right —
it is a structural property of "how many alternate routes are provisioned," which would hold for
the real optimised topology too if it shared a similar order of redundancy. This decouples the
tractability finding from the proxy-fidelity question a reviewer might otherwise conflate.

---

## 2. BDD-vs-IPA comparison — is the environment fair?

**The risk**: reviewer #2 already caught one uncontrolled comparison (the 134× dPrPm figure); a
second comparison that turns out to be unfair (different language, different hardware, a cold vs.
warm timing asymmetry) would be far more damaging the second time.

**Recommendation: state the environment once, explicitly, as its own short paragraph, and reuse
it everywhere.** The ingredients are already scattered across the validation notes but not
consolidated into one manuscript-facing statement. Recommended text for the §5 preamble (this
already exists in `RESS_edit_proposals.md` EDIT 11 in part — extend it to name the BDD side too):

> "All decision-diagram comparisons in this section use CUDD, called from the same Julia process
> and the same machine as the proposed method, with dynamic variable reordering (sifting) enabled
> unless stated otherwise; a fixed-order (unsifted) build is additionally reported wherever it
> changes the conclusion. All runtimes, for both methods, were measured after a discarded warm-up
> run of the same computation, so that JIT/program-initialisation cost is excluded from every
> figure reported."

This one paragraph pre-empts the fairness question for every subsequent number (grid, corpus,
drone) rather than requiring a per-table caveat. It is directly supported by the validation
record: the campaign notes document a real instance where an *unwarmed* comparison flipped the
conclusion (BDD looked faster than IPA until warm-up was added; corrected, IPA was 14.38× faster)
— this is exactly the kind of mistake the paragraph above forecloses, and the manuscript should
say plainly that this check was done, not just that the number is now warm. Suggested addition to
§5.4.3: "An initial unwarmed measurement of this comparison reversed the conclusion; the
discrepancy was traced to JIT/compilation cost and eliminated by the warm-up protocol stated
above — we note this explicitly because it is exactly the class of methodological error the
reviewers rightly flagged in the original dPrPm comparison, and we wanted the same scrutiny
applied to our own new comparison."  Naming your own near-miss is stronger than silently fixing
it — it demonstrates the rigour reviewer #2 asked for, rather than merely asserting it.

Also recommend keeping the scope sentence already drafted ("no claim is made... nor that the
crossover point has been located precisely") — it is the right level of honesty and should not be
strengthened into a general "BDD cannot handle drone-scale networks" claim.

---

## 3. "Broader applicability to Bayesian networks" (R3 comment 8) — safe to claim?

**The risk**: claiming the method "transfers" to Bayesian-network inference invites a PGM-literate
reviewer to ask why, if so, IPA isn't just cutset conditioning with extra steps — undermining the
novelty argument the manuscript makes elsewhere.

**Recommendation: the existing draft response is already appropriately narrow — strengthen it
with the credal-network complexity literature instead of broadening it.** The already-drafted
answer ("the machinery transfers to source-to-node reachability queries on other directed acyclic
probabilistic models... the capability that does not transfer back is imprecise propagation")
is the right scope: transfer the *structural* relationship (cutset-conditioning specialisation),
not a claim that IPA "does Bayesian inference." The literature check done today
(`literature/CITATION_STATUS.md` item 1) gives this a much sharper, citable form. Recommended
replacement paragraph for §4.3/response to R3.8:

> "The relation to exact inference in probabilistic graphical models is structural, not
> incidental: the method is a specialisation of cutset conditioning (Pearl 1988), in the same
> width-governed complexity class as junction-tree inference (Lauritzen & Spiegelhalter 1988), so
> the underlying machinery — local separating sets, memoised conditional maps, independent-branch
> factorisation — transfers to source-to-node reachability queries on any directed acyclic
> probabilistic model, Bayesian networks included. What does not transfer automatically is the
> paper's second contribution, exact imprecise propagation: for classical (precise) Bayesian
> networks, bounded treewidth guarantees polynomial exact inference, but this guarantee is known
> not to survive the move to imprecise (credal) parameters — bounded-treewidth credal-network
> inference remains NP-hard in general (Mauá and Cozman 2020), and the only known tractable
> relaxation is an approximation scheme, not an exact one (Fagiuoli and Zaffalon 1998, for the
> polytree/binary-variable special case). The present method's width-governed *exact* imprecise
> result is therefore specific to the arithmetic structure of the reachability problem and the
> explicit diamond-conditioning construction, not a free consequence of bounded width that would
> automatically extend to general credal-network inference."

This is a stronger, not weaker, claim than the current draft, because it is backed by a real
complexity-theoretic contrast rather than an assertion — and it directly answers "why hasn't
someone already done this for Bayesian networks" with a citable reason (credal-net inference
doesn't get the same free pass that precise inference does). Needs Fagiuoli & Zaffalon (1998) and
Mauá & Cozman (2020) added to the bibliography (see `literature/CITATION_STATUS.md`).

---

## Net effect if all three recommendations are accepted
- §5.4.1 gains a short, explicit "these are proxies, here is exactly how" disclosure (mostly
  already drafted, just needs to survive un-softened).
- §5 preamble gains one paragraph naming the shared BDD/IPA measurement environment; §5.4.3 gains
  one sentence disclosing the caught warm-up bug as evidence of rigour, not hiding it.
- §4.3 (and the R3.8 tracker entry) gains two new citations and a sharper, more defensible
  broader-applicability paragraph.
- Two new bibliography entries: Fagiuoli & Zaffalon (1998), Mauá & Cozman (2020).
