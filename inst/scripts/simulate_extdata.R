## Simulate example coverage.tab files shipped in inst/extdata/.
##
## Running this script from the package source root regenerates the
## four example samples: A1, A2 (group WT-GCV, slow coverage decay) and
## B1, B2 (group WT+GCV, fast coverage decay). The files are written
## into inst/extdata/<sample>/08_mutation_identification/ so that the
## example directory layout matches breseq's real output structure.
##
## Design goals for the simulated data:
##
## 1. Realistic schema: all eight columns produced by breseq
##    (unique_top_cov, unique_bot_cov, redundant_top_cov,
##    redundant_bot_cov, raw_redundant_top_cov, raw_redundant_bot_cov,
##    e_value, position).
##
## 2. Small enough to bundle: 2000 bp per sample x 4 samples produces
##    ~270 KB of extdata total.
##
## 3. Clear, visible group effect. WT+GCV samples decay six times
##    faster than WT-GCV samples, so plot_coverage_decay() shows an
##    obvious separation and test_decay_slopes() returns a highly
##    significant interaction.
##
## To regenerate the shipped files, run from the package source root:
##
##   setwd("/path/to/CovDecayAnalyzer")
##   source("inst/scripts/simulate_extdata.R")
##
## The current shipped files were produced with an equivalent Python
## script using the same design parameters. Re-running this R script
## will produce structurally similar (but not bit-identical) files
## because R and Python use different RNG streams; the analytical
## behavior is unchanged.

set.seed(2026)

genome_length <- 2000L
anchor        <- 1000L
base_coverage <- 500L

configs <- list(
    list(sample = "A1", group = "WT-GCV", decay = 5e-4, mult = 1.00),
    list(sample = "A2", group = "WT-GCV", decay = 5e-4, mult = 0.95),
    list(sample = "B1", group = "WT+GCV", decay = 3e-3, mult = 1.05),
    list(sample = "B2", group = "WT+GCV", decay = 3e-3, mult = 0.90)
)

outdir_root <- "inst/extdata"

positions <- seq_len(genome_length)
distance  <- abs(positions - anchor)

for (cfg in configs) {
    subdir <- file.path(outdir_root, cfg$sample, "08_mutation_identification")
    dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

    expected_total <- base_coverage * cfg$mult * exp(-cfg$decay * distance)
    total <- rpois(genome_length, lambda = pmax(expected_total, 1e-6))

    ## ~50/50 split across the two unique strands
    unique_top <- rbinom(genome_length, size = total, prob = 0.5)
    unique_bot <- total - unique_top

    ## Redundant coverage is typically zero in these regions
    redundant_top     <- integer(genome_length)
    redundant_bot     <- integer(genome_length)
    raw_redundant_top <- integer(genome_length)
    raw_redundant_bot <- integer(genome_length)

    ## Plausible per-position e-values
    e_values <- round(runif(genome_length, min = 20000, max = 40000), 1)

    out <- data.frame(
        unique_top_cov        = unique_top,
        unique_bot_cov        = unique_bot,
        redundant_top_cov     = redundant_top,
        redundant_bot_cov     = redundant_bot,
        raw_redundant_top_cov = raw_redundant_top,
        raw_redundant_bot_cov = raw_redundant_bot,
        e_value               = e_values,
        position              = positions
    )

    outpath <- file.path(subdir, "Exported.coverage.tab")
    write.table(
        out, file = outpath,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    message("wrote ", outpath)
}

## Companion metadata table
meta <- data.frame(
    sample = vapply(configs, `[[`, character(1L), "sample"),
    group  = vapply(configs, `[[`, character(1L), "group")
)
write.csv(meta, file.path(outdir_root, "sample_metadata.csv"),
    row.names = FALSE, quote = FALSE)
message("wrote ", file.path(outdir_root, "sample_metadata.csv"))
