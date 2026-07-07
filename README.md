# CovDecayAnalyzer

A packaged Shiny app for analyzing sequencing coverage decay relative to an OriLyt (origin of lytic replication) position across sample groups.

## Capabilities

- **Customizable Setup** — configure input directory, genome length, bin size, OriLyt position, and sample metadata
- **Coverage Plots** — normalized coverage vs. distance from OriLyt, LOESS-smoothed, with optional region filtering
- **Decay Slopes** — per-sample log-linear decay slopes, boxplots by group, and statistical tests (ANOVA, Kruskal-Wallis, pairwise permutation)
- **Sample Explorer** — per-sample coverage profile with linear fit overlay

## Requirements

``` r
install.packages(c(
  "shiny", "bslib", "data.table", "ggplot2",
  "broom", "rcompanion"
))
```

## Input Data

For each sample, the app expects a coverage file at:

```         
<base_dir>/<sample>/08_mutation_identification/Exported.coverage.tab
```

with columns:

- `position`
- `unique_top_cov`, `unique_bot_cov`, `redundant_top_cov`, `redundant_bot_cov`

## Sample Metadata Format

Pasted as CSV in the Setup tab:

``` csv
sample,group
10_WTyesGCV_4_S24,WT+GCV
...
```

## Running the App

``` r
CovDecayAnalyzer::covApp()
```

Then set the base directory, load metadata, and click **Load & Process Data**.

## Processing Pipeline

1.  Total coverage = sum of unique/redundant top/bottom coverage
2.  Coverage binned by configurable bin size
3.  Normalized to per-sample median coverage
4.  Distance from OriLyt computed per bin
5.  Log-transformed (`log(norm_cov + 1e-5)`) for linear decay modeling
6.  Per-sample decay slope fit via `lm(log_norm_cov ~ dist_from_ori)`

## Downloadable Outputs

- Coverage plot (TIFF) and binned coverage table (CSV)
- Decay slope boxplot (TIFF) and slope table (CSV)
