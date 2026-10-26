# Methods framing memo

**Status.** Input document for the methods paper's methods section.
Sections 1 through 6 and Section 8 are appropriate for excerpting into
the lab-internal methods document that supersedes the earlier
Kruskal-Wallis specification. Section 7 is paper-specific.

**Last updated.** 2026-10-24, after Tier 2 simulation study completion.

**Scope note.** This memo covers methods only. Introduction, results,
discussion, abstract, and figure design are deferred to subsequent
paper drafting sessions. Mathematical notation uses standard inline
LaTeX for pandoc rendering.

## 1. The method and its scope

CovDecayAnalyzer analyzes how sequencing coverage of a reference genome
changes with distance from a user-supplied anchor position, and tests
whether that pattern differs between experimental groups. The method is
designed for whole-genome shotgun sequencing data from dsDNA organisms
with a localized replication origin, where coverage near the origin is
enhanced during active replication relative to coverage at distant
positions. The anchor represents the position where coverage enhancement
is expected to be maximal, typically a known or putative replication
origin.

For a given sample $i$ and a set of binned, median-normalized coverage
values $y_{ij}$ at bin $j$, the method tests whether the slope of
$\log y_{ij}$ against distance from the anchor differs systematically
between the experimental groups to which samples belong. The null
hypothesis is that all groups share the same decay slope; the
alternative is that at least one group's slope differs. A secondary
test addresses whether the decay pattern at the supplied anchor is
specific to that position versus a general property of the reference
genome.

The method applies to linear reference coordinates. Users with circular
physical genomes represented in linear references can run the method
directly; the method does not treat position $L$ as adjacent to
position $1$. This limitation is noted here and discussed as a bounded
case in the paper's discussion.

## 2. Preprocessing

The method accepts per-base coverage counts from `breseq`-generated
output tables (`Exported.coverage.tab`) or any equivalent per-base
coverage schema (position and depth columns with implicit integer
position ordering).

The genome is partitioned into non-overlapping bins of fixed size $B$.
The default bin size is $B = 5000$ bp; a smaller bin size of
$B = 1000$ bp is supported and produces calibrated results in the
simulation study (Section 7). For each sample $i$ and bin $j$ spanning
reference positions $[(j-1)B + 1, jB]$, the per-base coverage values
within the bin are averaged to produce a mean coverage $c_{ij}$.

To account for differences in sequencing depth between samples, mean
coverage per bin is median-normalized within each sample:

$$y_{ij} = \frac{c_{ij}}{\mathrm{median}_j(c_{ij})}$$

The log-transformed normalized coverage is computed with a small
pseudocount to accommodate zero-coverage bins:

$$\ell_{ij} = \log(y_{ij} + \epsilon), \quad \epsilon = 10^{-5}$$

The distance from each bin's midpoint to the anchor position $a$ is
computed as:

$$d_j = |m_j - a|, \quad m_j = (j - 0.5)B$$

Optionally, the analysis can be restricted to a window of width $W$
around the anchor, excluding bins with $d_j > W$. The default is no
window restriction.

## 3. Primary test: linear mixed-effects model

The primary test is a linear mixed-effects model fit to all per-bin
log-normalized coverage values across all samples and groups
simultaneously. Let $g(i)$ denote the group to which sample $i$
belongs, with $G$ total groups and $n_g$ samples in group $g$. The
model is:

$$\ell_{ij} = \beta_{0,g(i)} + \beta_{1,g(i)} d_j + u_{0,i} + u_{1,i} d_j + e_{ij}$$

where $\beta_{0,g}$ and $\beta_{1,g}$ are group-level intercept and
decay slope parameters, $u_{0,i}$ and $u_{1,i}$ are sample-level random
intercept and random slope with joint distribution:

$$\begin{pmatrix} u_{0,i} \\ u_{1,i} \end{pmatrix} \sim \mathcal{N}\left(\mathbf{0}, \begin{pmatrix} \tau_0^2 & \rho \tau_0 \tau_1 \\ \rho \tau_0 \tau_1 & \tau_1^2 \end{pmatrix}\right)$$

and $e_{ij} \sim \mathcal{N}(0, \sigma^2)$ are independent residuals.
The random intercept $u_{0,i}$ captures sample-to-sample variation in
overall coverage level after median normalization, and the random
slope $u_{1,i}$ captures sample-to-sample variation in decay rate
within a group.

The null hypothesis is
$H_0: \beta_{1,1} = \beta_{1,2} = \cdots = \beta_{1,G}$,
that is, all groups share the same decay slope. The test statistic is
the $F$-statistic from a Wald test on the group-by-distance interaction
term, with Kenward-Roger small-sample degree-of-freedom correction.
The $p$-value is from the $F$ distribution with Kenward-Roger-
approximated numerator and denominator degrees of freedom.

Implementation uses `lmerTest::lmer` for the fit with REML estimation,
and `pbkrtest::KRmodcomp` for the Kenward-Roger correction. The
exported function is `test_decay_slopes`.

### 3.1 On the random intercept

Median normalization constrains each sample's intercept to the
sample-level mean, so the random intercept $u_{0,i}$ absorbs
sample-level variation that median normalization leaves. The estimated
correlation $\rho$ between random intercept and random slope is
typically close to $-1$ (median $-0.993$ across simulated cells; see
Section 7), and the covariance matrix is often reported as singular by
`lme4::isSingular`. Despite this near-degeneracy, a cautionary
parameterization that drops the random intercept entirely (random
slope only; Section 5.4) exhibits severe Type I error inflation to
$0.507$ (95% CI $0.493$ to $0.521$ at $R = 5000$, reference cell), with
additional bin-size dependence ($+0.037$ increase from 5 kb to 1 kb
binning at $R = 1000$). The random intercept is therefore retained in
the primary model; the Kenward-Roger correction correctly accommodates
the near-degenerate covariance estimate.

## 4. Secondary test: random-anchor null distribution

The primary test asks whether groups differ in decay slope at the
user-supplied anchor. A separate question is whether the decay pattern
at the anchor is specific to that position, or whether a comparable
pattern would emerge at randomly chosen positions on the reference.
The random-anchor null test addresses this question and is a
contribution unique to this package in the comparator set.

Given the observed $F$-statistic $F_{\text{obs}}$ from the primary test
at anchor $a$, draw $K$ random anchor positions $a_k$ uniformly without
replacement from the allowed set of integer positions on the reference.
If an exclusion radius $r$ is specified, random anchors are drawn from
positions outside $[a - r, a + r]$. The default $r$ is $0$ when no
window is specified, and $W$ when a window of width $W$ is specified.
For each random anchor $a_k$, the distances $d_j^{(k)} = |m_j - a_k|$
are recomputed and the primary LMM is refit, producing an
$F$-statistic $F_k$.

The empirical $p$-value is calculated as:

$$p_{\text{emp}} = \frac{1 + \#\{k : F_k \geq F_{\text{obs}}\}}{K + 1}$$

The $+1$ offset in numerator and denominator provides a conservative
floor and avoids the impossible value $p = 0$ when all $F_k$ are below
$F_{\text{obs}}$. This convention is widely used for Monte Carlo
$p$-values and is appropriate for the without-replacement sampling
employed here; it is more conservative than exact calibration would
require.

Fit failures at individual random anchors (singular models that fail
to converge) are counted and reported as an attribute on the returned
object, but are excluded from the empirical $p$-value calculation. In
the simulation study (Section 7), fit failure rates were below $1\%$
at $R = 5000$ across all tested cells.

Implementation functions: `null_distribution_F` for the Monte Carlo
null, and `test_decay_significance` for the combined test.

## 5. Alternative tests provided by the package

The package provides four alternative tests for direct comparison, with
explicit recommendations for the conditions under which each is
preferred (Section 6).

### 5.1 Pooled-variance two-stage ANOVA on per-sample slopes

For each sample $i$, fit an ordinary least squares regression
$\ell_{ij} = \alpha_i + s_i d_j + e_{ij}$ to obtain a per-sample slope
estimate $s_i$. The $n = \sum_g n_g$ slope estimates are then tested
with a pooled-variance one-way ANOVA (for $G \geq 3$) or a
pooled-variance $t$-test (for $G = 2$) under the null that all
$\mathbb{E}(s_i)$ are equal across groups.

The pooled two-stage test operates on slope estimates, not on raw
coverage values. OLS slope estimates on 20 to 50 binned observations
per sample are approximately normal under the generative conditions
tested in Section 7, by the usual central limit argument. The pooled
two-stage test is calibrated at $n$ as small as 2 per group in the
parameter space tested, with Type I in the same range as the primary
LMM across all cells.

Implementation: `base::lm` for per-sample fits and `base::t.test` or
`base::aov` for the second stage. Function:
`test_decay_slopes_two_stage(variant = "pooled")`.

### 5.2 Welch's two-stage test

Identical to pooled two-stage at the first stage. At the second stage,
uses Welch's variance-corrected $t$-test for $G = 2$ or Welch's ANOVA
for $G \geq 3$. Appropriate in principle when group-level slope
variance differs between groups.

In simulations, Welch's two-stage loses power at $G = 2$ with small
$n$ relative to pooled (power $0.17$ versus $0.39$ at
$n = 3, \Delta_{\text{std}} = 0.5$; see Section 7). The power loss is
large relative to the variance-robustness gain for the usage pattern
this package supports. Welch's two-stage is retained in the package as
a methodological alternative but is not recommended for the common
case.

Function: `test_decay_slopes_two_stage(variant = "welch")`.

### 5.3 Kruskal-Wallis plus Dunn's post-hoc

Non-parametric alternative. Per-sample slopes are computed as in
Section 5.1 and then analyzed with the Kruskal-Wallis omnibus test on
ranks, with Dunn's post-hoc for pairwise comparisons under Holm
multiplicity adjustment.

Three variants are supported.

**Chi-square $p$-value.** The KW statistic is referred to a chi-squared
distribution with $G - 1$ degrees of freedom. The chi-square
approximation is miscalibrated at small $n$ because the KW statistic
has a discrete sampling distribution whose exact tail probabilities
deviate from the continuous approximation. Simulation under the null
at $G = 2, n = 3$ shows the chi-square-approximated test rejects at
nominal $\alpha = 0.05$ with probability $0.0996$, approximately
double the intended rate.

**Permutation $p$-value.** The KW statistic is calibrated by exact
enumeration (when the number of distinct group-label assignments is
below $2 \times 10^5$) or by random permutation otherwise. Correct in
principle, but has a power floor at small $n$: at $G = 2, n = 3$, the
minimum achievable $p$-value exceeds $\alpha = 0.05$ because of the
small number of distinct permutations, so the test cannot reject
regardless of effect size (simulation confirms Type I of $0$ at this
cell).

**Permutation Dunn's.** Pairwise rank-sum differences are calibrated
by permutation, with Holm multiplicity adjustment. Simulations show
permutation Dunn's agrees with normal-approximation Dunn's within
Monte Carlo noise at $G = 3, n = 5, \Delta_{\text{std}} = 1$ (rejection
rates $0.962$ versus $0.958$, respectively), indicating that Dunn's
chi-square approximation is adequate at the post-hoc level even when
the omnibus KW chi-square is miscalibrated.

The KW+Dunn method has been used historically in the laboratory
domain this package serves. It is retained in the package for
reproducibility of prior analyses, with explicit documentation of its
small-$n$ limitations.

Functions: `test_decay_slopes_kw(variant = "chisq")`,
`test_decay_slopes_kw(variant = "perm")`, with separate functions for
the two Dunn's variants.

### 5.4 Slope-only cautionary model

The primary LMM with the random intercept $u_{0,i}$ removed, retaining
only the random slope $u_{1,i}$. Documented in Section 3.1 as a
cautionary comparator demonstrating that the random intercept absorbs
real correlation structure. Not recommended for data analysis; included
for pedagogical purposes.

Function: `test_decay_slopes(random_structure = "slope_only")`.

## 6. Recommendations by study design

The simulation study (Section 7) characterizes the primary test and
all alternatives across a parameter space spanning plausible scenarios
for dsDNA viruses with localized replication origins. The following
recommendations are derived from that characterization.

**Balanced designs with $n \geq 3$ per group, no known
heteroscedasticity, dropout rate below $2\%$.** Use the primary LMM
(`test_decay_slopes`). The pooled-variance two-stage ANOVA
(`test_decay_slopes_two_stage(variant = "pooled")`) is approximately
equivalent in simulations: paired rejection differences within $0.006$
across the 63-cell figure grid at $R = 1000$.

**Unbalanced designs.** Either the primary LMM or the pooled two-stage.
Simulations across 30 unbalanced cells show paired rejection
differences within $0.002$.

**Known or suspected heteroscedastic residual variance between groups
at $n = 3$.** Use the pooled two-stage. The primary LMM is calibrated
but conservative under these conditions, with up to $0.041$ paired
rejection difference below pooled (at $G = 2, n = 3,
\Delta_{\text{std}} = 0.5$ with residual SD ratio 2 and
$\sigma_{\text{atom}} = 0.5$). At $n \geq 5$ or at homoscedastic
residuals, the two methods agree within $0.015$.

**High dropout rates ($\geq 2\%$ at 5 kb run lengths).** Preprocess to
address dropouts before applying either test. Both the primary LMM and
the pooled two-stage are below nominal in simulations (LMM Type I
$0.005$ to $0.011$; pooled Type I $0.020$ to $0.043$). The LMM power
loss relative to pooled is up to $0.112$ at effect cells. Neither
method is calibrated under these conditions, and no comparator in the
package is.

**Smallest designs ($G = 2, n = 2$).** Primary LMM, which achieves
$0.381$ power at $\Delta_{\text{std}} = 1$ with correct $\alpha = 0.05$
Type I. Welch's two-stage loses power substantially at small $n$
(Section 5.2); the KW permutation variant cannot reject at $n \leq 3$
because of its power floor (Section 5.3); the KW chi-square variant is
anticonservative at these sizes.

**Anchor near reference ends ($0.0 L$ or $1.0 L$, single-arm
geometry).** Mild Type I elevation at $G = 3, n = 3$ ($0.063$ to
$0.070$), affecting all tested methods similarly. Interior anchors
($0.1 L$ to $0.9 L$) show no elevation in simulations.

**Assessing anchor specificity.** Use `test_decay_significance` with
the random-anchor null distribution (Section 4). No comparator provides
an equivalent test.

**Historical reproducibility.** Users reproducing analyses from the
earlier Kruskal-Wallis plus Dunn literature should use
`test_decay_slopes_kw` with the specific variant that was originally
applied, documenting its small-$n$ limitations where applicable.

## 7. Simulation study

Characterization of the primary test and alternatives was performed
across a parameter space defined over:

- Samples per group: $n \in \{2, 3, 4, 5, 6, 8, 10\}$
- Number of groups: $G \in \{2, 3, 4\}$
- Standardized effect size:
  $\Delta_{\text{std}} = \Delta / s_{\text{ref}} \in \{0, 0.5, 1.0, 1.5\}$
- Coefficient of variation: $\mathrm{cv} = \tau / s_{\text{ref}}$ in a
  range from $0.1$ to $0.5$
- Residual noise amplitude: $\sigma_{\text{atom}}$ in a range from
  $0.05$ to $0.5$
- Genome length: $L \in \{50, 100, 238, 500\}$ kb
- Reference decay rate: $s_{\text{ref}} \in
  \{3 \times 10^{-6}, 6 \times 10^{-6}, 1.2 \times 10^{-5}\}$ per bp
- Anchor position: $a/L \in \{0.1, 0.25, 0.5, 0.75, 0.9\}$, with
  sensitivity cells at the extremes $0.0$ and $1.0$
- Bin size: $B \in \{1, 5\}$ kb (main grid at 5 kb; 1 kb as paired
  sensitivity at the reference cell)
- Group-effect pattern: "one group differs from a shared reference"

Sensitivity analyses perturbed the baseline model with autocorrelated
residuals, zero-coverage dropouts at varying rates and run lengths,
heteroscedastic residual variance between groups, heteroscedastic
between-sample slope variance, and anchor positions at the reference
extremes.

The study ran approximately $800{,}000$ simulated datasets analyzed
by seven arm-variants (primary LMM, slope-only cautionary, two
two-stage parametric variants, three KW+Dunn variants). Infrastructure
and seeds are fixed; the complete pipeline is reproducible via
`targets::tar_make()` from the `simulations/` directory.

Primary findings referenced throughout Sections 3, 4, 5, and 6:

- Primary LMM Type I: $0.043$ to $0.054$ across all 21 null cells of
  the main grid at $R = 5000$ ($135{,}000$ datasets).
- Primary LMM power at $G = 2, n = 5, \Delta_{\text{std}} = 1$:
  $0.985$ (95% CI $0.98$ to $0.99$).
- Smallest $n$ for $80\%$ power ($G = 2$, reference): $n = 10$ at
  $\Delta_{\text{std}} = 0.5$, $n = 4$ at $1.0$, $n = 3$ at $1.5$.
- Approximate collapse of power onto $\Delta / \tau_{\text{eff}}$,
  where $\tau_{\text{eff}}$ combines between-sample slope variance and
  within-sample residual noise. Deviations across $s_{\text{ref}}$
  levels match the noncentral-$F$ prediction with this explicit
  $\tau_{\text{eff}}$ definition. The sample size calculator provided
  by the package indexes on $\Delta / \tau_{\text{eff}}$.

Full simulation results are archived at `simulations/results/tier2/`
with per-replicate metadata enabling independent reanalysis.

## 8. Reproducibility

Simulation code, generator, analysis arms, and infrastructure live at
`simulations/` in the project repository, with dependencies pinned via
`renv.lock`. The full pipeline is reconstructable via
`targets::tar_make()` from the `simulations/` directory.

Summary tables used in Sections 3 through 7 are committed in
`simulations/results/tier2/`. Per-replicate records are either
committed (Tier 1 exploratory and A2-3 focused cells) or archived
externally with a committed checksum manifest (Tier 2 production;
archive location documented in the repository README).

Session logs documenting each simulation run are at `session-logs/`.
The session log for the most recent Tier 2 run is
`session-logs/session-2026-10-23-a2-4.md`, with reruns at
`session-logs/session-2026-10-24-s3-rerun.md` and
`session-logs/session-2026-10-24-s3-rerun-baseline.md`.

Package source: `https://github.com/zebragirl314/CovDecayAnalyzer`.
Package tests pass R CMD check with zero errors, zero warnings, and
zero notes at the time of this writing. Overall line coverage is above
$91\%$.
