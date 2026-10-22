# Vignette drafts for `vignettes/null-distribution.Rmd`

Drafts only. Nothing here has been applied to the vignette. Written for
Carly to edit and land in a separate commit.

---

## Draft 1 (Item 3): the "Empirical p-value" paragraph

**Currently in the vignette** (the paragraph after the formula):

> The `+1` in both numerator and denominator is not smoothing. It is
> the correct exact expression for a Monte Carlo test where the null
> samples are drawn with replacement, as derived by Phipson and Smyth
> (2010). Using the more intuitive $k/n$ estimator systematically
> understates the true p-value by roughly $1/n$ and can produce
> impossible p-values of exactly zero. The corrected form has a
> minimum achievable value of $1/(n+1)$, which is the smallest p-value
> that a Monte Carlo test with $n$ samples can honestly report.

**Proposed replacement:**

> We compute the p-value with the $+1$ in both numerator and
> denominator, and we use it as a conservative convention rather than as
> an exact result. The reasoning is that the observed statistic is itself
> one plausible member of the reference set: if the anchor is nothing
> special, $F_{\text{obs}}$ is just another draw alongside the $n$ random
> anchors. Counting it keeps the p-value strictly above zero, so we never
> report an impossible $p = 0$, and it never gives a smaller value than
> the plain $k/n$ estimator. The random anchors are sampled without
> replacement from the finite set of integer genomic positions, not drawn
> with replacement, so this is a different resampling scheme from the one
> Phipson and Smyth (2010) analyze. We cite their paper for the case
> against reporting zero p-values from Monte Carlo samples, and for the
> $(k+1)/(n+1)$ form that comes out of it, not as a proof that the form is
> exact here. The smallest value the test can report is $1/(n+1)$.

### Editorial notes for Carly

- I reversed the framing from "this is the correct exact expression" to
  "this is a conservative convention". The original sentence is the part
  that contradicts the code (`sample.int(..., replace = FALSE)` in
  `null_distribution_F`), and it also overclaims: Phipson and Smyth's
  derivation is for p-values from randomly drawn permutations of a
  permutation test, not for anchor positions.
- "Never gives a smaller value than $k/n$" is true for any $k \le n$ and
  I checked it: $(k+1)/(n+1) \ge k/n$ exactly when $k \le n$. I kept the
  claim to that safe form and did not say "conservative with respect to
  the true p-value", because I cannot show that. Whether the observed
  anchor really is exchangeable with random anchors under the null is
  an assumption, and the vignette's "Assumptions of the random-anchor
  null" section is where that assumption belongs. I did not add it to
  the draft to keep the length close to the original; you may want a
  one-line pointer.
- I dropped the original's "systematically understates the true p-value by
  roughly $1/n$" sentence. That is accurate for permutation tests with
  replacement but I could not justify it for our scheme.
- Related, outside this paragraph: the "Construction" section says
  `null_distribution_F()` "samples $N$ positions uniformly". The code
  samples without replacement, so "uniformly without replacement" would be
  more accurate. It also says failed fits "are dropped"; the count is now
  returned (`n_null_failed`), which the section could mention.
- Another effect worth a sentence somewhere: with 1 kb bins, anchors in the
  same bin give almost the same F, and nearby anchors give correlated F, so
  the effective null size is smaller than $n$. I put this in the roxygen for
  `null_distribution_F`; I left it out of the draft paragraph.
- I am not certain the phrase "alongside the $n$ random anchors" is the way
  you would want to explain the +1. Please change it to your own wording.

---

## Draft 2 (Item 4): the "cannot occur by construction" cell

**A labeling note first.** The session prompt calls this the "upper-right"
cell. In the vignette's table, the cell that says "Cannot occur by
construction" is **large LMM p / small empirical p**, which is the
*lower-left* cell, and the text below the table already calls it "the
lower-left cell". The upper-right cell (small LMM p / large empirical p)
is the one the text describes as "most nuanced", and it is correct as
written. I wrote the draft for the cell whose content matches the
prompt's description, the one that makes the false claim. If you meant a
different cell, tell me.

**Currently in the table (lower-left cell):**

> Cannot occur by construction: if the observed F is unremarkable, it
> cannot be extreme relative to a null derived from data with the same
> structure.

**Proposed replacement for the cell:**

> Possible, and informative. The effect at this anchor is not individually
> significant, but F at random anchors is systematically lower still, so
> this anchor stands out from the rest of the genome. Interpret as
> "unusual for this genome, not strong on its own".

**The paragraph below the table currently says:**

> The lower-left cell is precluded because the empirical p-value uses
> $F_{\text{obs}}$ as the reference; if that reference is not
> extreme, no null can make it look extreme.

That paragraph is also wrong, and I think it has to change with the cell or the
vignette will contradict itself. **Proposed replacement:**

> The lower-left cell is not precluded. The two p-values compare
> $F_{\text{obs}}$ to different references: the LMM p-value to an $F$
> distribution with Kenward-Roger degrees of freedom, the empirical
> p-value to the $F$ values observed at other anchors in this genome. If
> the group effect is weak everywhere, $F_{\text{obs}}$ can be modest
> (say, an LMM p of 0.2) and still exceed nearly all of the null, because
> the null values are smaller. The combination says the anchor is unusual
> relative to the rest of the genome even though the effect at it is not
> individually significant. We would treat it as a prompt to look
> further, not as evidence of an effect, since the group difference
> itself is not supported at this anchor.

### Editorial notes for Carly

- The logic of the original claim treats the null F values as if they
  could never be smaller than the observed F when the observed F is
  unremarkable. They are built from fits at other anchors on the same
  data, so nothing ties them to the nominal $F$ reference. When group
  effects are weak everywhere, null F values cluster near small numbers
  and a modest $F_{\text{obs}}$ can sit in their upper tail.
- I chose "informative" and "prompt to look further" deliberately. I want
  to avoid implying the result supports a biological claim about the
  anchor; the LMM p-value says the group difference is not itself supported.
  You may prefer a more cautious or a more enthusiastic reading.
- I am unsure how often this cell arises in practice. My intuition is that
  it is uncommon when the null is roughly $F$-distributed and more likely
  when the null is concentrated at small values, but I have no evidence;
  the simulation study would give it. I avoided saying "rare" in the draft for
  that reason.
- The table cell is deliberately short to fit the column. The long form is
  in the paragraph.
- The following paragraph in the vignette says the upper-right case
  "suggests one of two possibilities". I left that text alone; it reads
  correctly.

---

## Applied in Session 3 (Item 4 only)

Authorized by Carly after Session 2. Both edits are now in the working tree
of `vignettes/null-distribution.Rmd` (uncommitted; left for your editorial
pass). The text is exactly Draft 2 above, with no changes.

- **Lower-left table cell** ("Large LMM p" row, "Small empirical p"
  column): replaced "Cannot occur by construction ..." with "Possible, and
  informative. ...". This confirms the cell identification in the labeling
  note: the lower-left cell was the right one.
- **Paragraph after the table**: "The lower-left cell is precluded ..."
  replaced with "The lower-left cell is not precluded. ..." (full text as
  in Draft 2).

Not applied: Draft 1 (the "Empirical p-value" paragraph, Item 3). The
vignette still says the `+1` form is "the correct exact expression" for
samples "drawn with replacement", which now contradicts the roxygen for
`null_distribution_F` and `empirical_pvalue`. The Construction-section
wording notes (without replacement, `n_null_failed`) are also still open.
The "Assumptions" section also still says failed fits "are dropped and the
effective n is reported" which is roughly right but could point to
`n_null_failed`.

---

## Applied in Session 3, second pass (Item 3)

Draft 1 is now in `vignettes/null-distribution.Rmd`, replacing the "The `+1`
... is not smoothing" paragraph under "Empirical p-value". Text is verbatim
from Draft 1; it still needs your editorial pass (see the notes under
Draft 1). Uncommitted.

Still not touched: the "Construction" section's "samples $N$ positions
uniformly" (code is without replacement) and "are dropped" (count is now
returned as `n_null_failed`), and the matching line in "Assumptions".
