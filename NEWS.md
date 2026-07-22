# CovDecayAnalyzer 0.1.0

Initial development release. Full end-to-end workflow now available in
package functions with roxygen-generated documentation.

## Coverage input and normalization

The loading API is layout-agnostic: the sample metadata table carries a
`coverage_file` column with a direct path per sample, so the package
does not assume any particular on-disk directory structure. A separate
`breseq_coverage_paths()` helper is provided for the common case of
breseq output in the standard layout.

* `read_coverage_tab()` reads a tab-delimited per-base coverage file
  and computes total per-base coverage.
* `bin_coverage()` aggregates per-base coverage into fixed-width bins.
* `normalize_coverage()` median-normalizes binned coverage and applies
  the `log(x + epsilon)` transform used by the decay model.
* `load_sample_coverage()` composes read + bin + normalize for one
  file at a user-supplied path.
* `load_all_samples()` iterates the above over a metadata table with
  a `coverage_file` column, with error tolerance and optional
  parallelism via BiocParallel.
* `breseq_coverage_paths()` builds the `coverage_file` column from a
  sample-name vector and a breseq output root directory. The only
  breseq-specific function in the package.

## Distance calculation and per-sample slopes

* `add_distance_from_anchor()` computes absolute distance from any
  anchor coordinate; supports the OriLyt, genome start, and genome end
  analyses in a single function.
* `fit_decay_slopes()` fits per-sample log-linear decay regressions
  as a diagnostic view.

## Formal group comparison

* `test_decay_slopes()` fits a linear mixed-effects model with a
  distance-by-group interaction, tested by the Kenward-Roger F
  approximation for small-sample inference. Returns per-group slope
  estimates and Tukey-adjusted pairwise contrasts via emmeans.

## Anchor specificity

* `null_distribution_F()` refits the mixed model at N random anchor
  positions to build a null distribution of interaction F-statistics.
* `empirical_pvalue()` computes the Monte Carlo `(k+1)/(n+1)`
  continuity-corrected p-value of Phipson & Smyth (2010).
* `test_decay_significance()` end-to-end wrapper combining the anchor
  fit with a random-anchor null distribution.

## Visualization

* `plot_coverage_decay()` LOESS-smoothed coverage vs distance by group.
* `plot_decay_slopes()` boxplot of per-sample decay slopes by group.
* `plot_null_distribution()` histogram of random-anchor null with
  observed F-statistic marked.

## Viral RPM

* `compute_viral_rpm()` vectorized RPM calculation.
* `plot_rpm_by_group()` mean-plus-error-bar barplot by treatment group.

## S3 classes

* `decay_test` (from `test_decay_slopes`) with `print()` method.
* `decay_significance` (from `test_decay_significance`) with `print()`
  method.

## Included data

* `inst/extdata/{A1,A2,B1,B2}/08_mutation_identification/Exported.coverage.tab`:
  four simulated coverage files matching breseq's schema, deliberately
  designed with a strong group difference in decay rate so example
  runs produce clearly interpretable results.
* `inst/extdata/sample_metadata.csv`: matching metadata table.
* `inst/scripts/simulate_extdata.R`: R recipe for regenerating the
  simulated dataset.

## Vignettes

* `intro`: end-to-end walkthrough of the pipeline on the shipped
  simulated data.
* `null-distribution`: methods reference covering the linear
  mixed-effects test and the random-anchor null distribution, with
  the assumptions, small-sample correction, and interpretation of
  observed-vs-null p-value combinations.
