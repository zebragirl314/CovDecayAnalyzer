# CovDecayAnalyzer

Coverage decay analysis around viral genome anchor positions, developed
for Human Cytomegalovirus (HCMV) TB40/E analyses of read coverage around
the lytic origin of replication (OriLyt) and at genome start and end
positions. The methods generalize to any linear reference genome with a
defined anchor coordinate.

## What the package does

Given per-base sequencing coverage input, the package:

1. Bins and median-normalizes coverage per sample.
2. Fits log-linear decay slopes with distance from a user-supplied
   anchor position (OriLyt, genome start, genome end, or arbitrary
   coordinate).
3. Compares slopes across treatment groups using a linear mixed-effects
   model with a Kenward-Roger small-sample F test for the
   distance-by-group interaction, and pairwise slope contrasts via
   emmeans.
4. Assesses whether the observed group difference is stronger than
   would be expected at a random genomic position by building a null
   distribution of interaction F-statistics from random anchors and
   reporting an empirical p-value.
5. Optionally summarizes viral reads-per-million (RPM) across
   treatment groups.

## Installation

Not yet on Bioconductor. From GitHub:

```r
# install.packages("remotes")
remotes::install_github("zebragirl314/CovDecayAnalyzer")
```

## Input data

The package reads tab-delimited per-base coverage files. Required
columns are `position`, `unique_top_cov`, `unique_bot_cov`,
`redundant_top_cov`, and `redundant_bot_cov`. Any additional columns
are permitted and ignored.

The package does not assume any particular on-disk layout. You tell
it where each sample's file lives via a `coverage_file` column in
your sample metadata table:

```r
meta <- data.frame(
    sample        = c("WT_1", "WT_2", "FV_1", "FV_2"),
    group         = c("WT",   "WT",   "FV",   "FV"),
    coverage_file = c(
        "/path/to/WT_1_coverage.tab",
        "/path/to/WT_2_coverage.tab",
        "/path/to/FV_1_coverage.tab",
        "/path/to/FV_2_coverage.tab"
    )
)
```

If your files come from breseq's mutation identification stage in the
standard layout (`<sample>/08_mutation_identification/Exported.coverage.tab`),
you can build the `coverage_file` column in one line with
`breseq_coverage_paths()`:

```r
meta$coverage_file <- breseq_coverage_paths(
    sample_names = meta$sample,
    base_dir     = "/path/to/breseq/output"
)
```

For any other tool or layout, populate `coverage_file` however you
like (`list.files()`, `file.path()`, a lookup table, etc.).

## Trying the pipeline on shipped example data

The package ships four simulated samples in `inst/extdata/` with a
deliberately strong group difference in decay rate, so a first-time
user can verify the pipeline works end-to-end before pointing it at
real data:

```r
library(CovDecayAnalyzer)

extdata_dir <- system.file("extdata", package = "CovDecayAnalyzer")
meta <- read.csv(
    file.path(extdata_dir, "sample_metadata.csv"),
    stringsAsFactors = FALSE
)
meta$coverage_file <- breseq_coverage_paths(meta$sample, extdata_dir)

combined <- load_all_samples(
    sample_metadata = meta,
    bin_size        = 100,
    genome_length   = 2000
)

## Visual: two clearly separated decay curves
plot_coverage_decay(combined, anchor_position = 1000, window = 1000)

## Formal test: strongly significant distance x group interaction
res <- test_decay_slopes(combined, anchor_position = 1000)
res
```

## Minimal example on real data

```r
library(CovDecayAnalyzer)

## 1. Build metadata with a coverage_file column
meta <- data.frame(
    sample = c("s1", "s2", "s3", "s4"),
    group  = c("A", "A", "B", "B")
)
meta$coverage_file <- breseq_coverage_paths(
    sample_names = meta$sample,
    base_dir     = "/path/to/breseq/output"
)

## 2. Read and normalize
combined <- load_all_samples(
    sample_metadata = meta,
    bin_size        = 1000,
    genome_length   = 238080
)

## 3. Diagnostic per-sample slopes
slopes <- fit_decay_slopes(
    combined,
    anchor_position = 93000,   # OriLyt coordinate; verify vs your reference
    window          = 50000
)
plot_decay_slopes(slopes)

## 4. Formal group comparison (LMM + Kenward-Roger F)
res <- test_decay_slopes(
    combined,
    anchor_position = 93000,
    window          = 50000
)
res              # prints F, p, per-group slopes, pairwise contrasts
res$pairwise     # Tukey-adjusted pairwise slope contrasts

## 5. Anchor-specificity test
sig <- test_decay_significance(
    combined,
    anchor_position = 93000,
    genome_length   = 238080,
    n_random        = 500,
    window          = 50000,
    exclude_radius  = 50000,
    seed            = 2026,
    BPPARAM         = BiocParallel::MulticoreParam(workers = 8)
)
sig                          # observed F, empirical p
plot_null_distribution(sig)
```

For a full walkthrough with the shipped simulated data, see
`vignette("intro", package = "CovDecayAnalyzer")`. For the statistical
methods behind `test_decay_slopes()` and `test_decay_significance()`,
see `vignette("null-distribution", package = "CovDecayAnalyzer")`.

## Method summary

* Coverage normalization: median-scaling per sample followed by
  `log(x + epsilon)` transform.
* Group comparison: linear mixed-effects model
  `log_norm_cov ~ distance * group + (1 + distance | sample)` fit by
  REML, with the `distance:group` interaction tested by the
  Kenward-Roger F statistic (Kenward & Roger 1997) using `lmerTest`
  (Kuznetsova et al. 2017) on top of `lme4` (Bates et al. 2015).
* Post-hoc contrasts: Tukey-adjusted pairwise slope contrasts via
  `emmeans::emtrends()` (Lenth 2024).
* Anchor specificity: empirical p-value comparing the observed F to a
  null of F-statistics from N random anchor positions, using the
  `(k+1)/(n+1)` continuity correction of Phipson & Smyth (2010).

## References

* Bates D, Machler M, Bolker B, Walker S (2015). Fitting linear
  mixed-effects models using lme4. *Journal of Statistical Software*
  67(1):1-48. doi:10.18637/jss.v067.i01
* Deatherage DE, Barrick JE (2014). Identification of mutations in
  laboratory-evolved microbes from next-generation sequencing data
  using breseq. *Methods in Molecular Biology* 1151:165-188.
  doi:10.1007/978-1-4939-0554-6_12
* Kenward MG, Roger JH (1997). Small sample inference for fixed
  effects from restricted maximum likelihood. *Biometrics*
  53(3):983-997. doi:10.2307/2533558
* Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package:
  Tests in linear mixed effects models. *Journal of Statistical
  Software* 82(13):1-26. doi:10.18637/jss.v082.i13
* Lenth RV (2024). emmeans: Estimated Marginal Means, aka Least-Squares
  Means. R package. https://CRAN.R-project.org/package=emmeans
* Phipson B, Smyth GK (2010). Permutation p-values should never be
  zero: calculating exact p-values when permutations are randomly
  drawn. *Statistical Applications in Genetics and Molecular Biology*
  9(1), Article 39. doi:10.2202/1544-6115.1585

## Authors

Olivia Daigle, Noelle Kosarek, Carly Bobak (maintainer).

## License

Artistic-2.0

## Development

To work on the package source:

```r
install.packages(c("devtools", "roxygen2"))

devtools::document()   # regenerate NAMESPACE and man/*.Rd from roxygen tags
devtools::test()       # run tests/testthat/*
devtools::check()      # full R CMD check
devtools::install()    # install from the local checkout
```

To regenerate the simulated example data in `inst/extdata/`:

```r
setwd("/path/to/CovDecayAnalyzer")
source("inst/scripts/simulate_extdata.R")
```
