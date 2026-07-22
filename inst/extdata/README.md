# inst/extdata

Simulated example data shipped with the package. All files here are
generated synthetically. No real experimental data is included, so
these files can be redistributed and used in public examples freely.

## Files

```
inst/extdata/
├── A1/08_mutation_identification/Exported.coverage.tab
├── A2/08_mutation_identification/Exported.coverage.tab
├── B1/08_mutation_identification/Exported.coverage.tab
├── B2/08_mutation_identification/Exported.coverage.tab
└── sample_metadata.csv
```

The directory layout matches what breseq produces, purely so users can
see how `breseq_coverage_paths()` works against a realistic tree. The
package itself does not depend on this layout: `load_all_samples()`
reads whatever paths you put in the metadata's `coverage_file` column.

## Sample design

Four samples across two treatment groups:

| Sample | Group  | Decay rate | Description                          |
|--------|--------|------------|--------------------------------------|
| A1     | WT-GCV | 5e-4 / bp  | Slow coverage decay from anchor      |
| A2     | WT-GCV | 5e-4 / bp  | Slow coverage decay from anchor      |
| B1     | WT+GCV | 3e-3 / bp  | Fast coverage decay from anchor      |
| B2     | WT+GCV | 3e-3 / bp  | Fast coverage decay from anchor      |

The simulated genome is 2000 bp long with the anchor at position 1000.
Coverage follows `total_cov = 500 * exp(-decay * |position - 1000|)`
plus Poisson noise, split ~50/50 between `unique_top_cov` and
`unique_bot_cov`. Redundant coverage columns are zero.

## Expected results

Because the group difference is deliberately strong, running the full
pipeline on these samples produces a highly significant group by
distance interaction, an empirical p-value at the floor
`1 / (n_random + 1)`, and clearly separated loess curves in
`plot_coverage_decay()`. This makes the shipped example useful as a
sanity check that the pipeline is installed and working correctly.

## Regeneration

The recipe lives in `inst/scripts/simulate_extdata.R`. To regenerate:

```r
setwd("/path/to/CovDecayAnalyzer")
source("inst/scripts/simulate_extdata.R")
```

Because R and Python use different RNG streams, regenerated files
will be structurally similar but not bit-identical to the ones
currently shipped. Analytical behavior is unchanged.

## Accessing the files from R

```r
library(CovDecayAnalyzer)

extdata_dir <- system.file("extdata", package = "CovDecayAnalyzer")

## Read one file directly
cov <- read_coverage_tab(
    file.path(extdata_dir, "A1", "08_mutation_identification",
              "Exported.coverage.tab")
)

## Or load all four samples through the full pipeline. Because the
## layout matches breseq's convention, the path helper works.
meta <- read.csv(file.path(extdata_dir, "sample_metadata.csv"))
meta$coverage_file <- breseq_coverage_paths(meta$sample, extdata_dir)
combined <- load_all_samples(
    sample_metadata = meta,
    bin_size        = 100,
    genome_length   = 2000
)
```
