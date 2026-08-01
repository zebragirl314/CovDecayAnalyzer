#' CovDecayAnalyzer: Coverage Decay Analysis Around Viral Genome Anchor Positions
#'
#' Tools for quantifying sequencing coverage decay around anchor positions
#' in viral genomes. The package was developed for Human Cytomegalovirus
#' (HCMV) TB40/E analyses of read coverage around the lytic origin of
#' replication (OriLyt) and at genome start and end positions, but the
#' methods are applicable to any linear reference genome with a defined
#' anchor coordinate.
#'
#' @section Motivation:
#'
#' In lytic HCMV infection, DNA replication initiates at OriLyt and
#' proceeds bidirectionally. When a mixed population of viruses is
#' sequenced during replication, coverage is higher near the origin
#' (because more copies of the origin-proximal sequence exist in the
#' population) and decays with distance from the origin. The rate of this
#' decay is informative about replication dynamics, and comparing decay
#' rates between treatment groups (for example, drug-treated versus
#' untreated, or wild-type versus mutant virus) can reveal
#' treatment-induced changes in replication behavior.
#'
#' @section Typical workflow:
#'
#' 1. Read and normalize coverage. [load_all_samples()] iterates over a
#'    metadata table of samples, reads each sample's breseq per-base
#'    coverage file, bins coverage into fixed-width windows, and
#'    median-normalizes each sample.
#'
#' 2. Fit and inspect per-sample decay slopes. [fit_decay_slopes()]
#'    fits a simple log-linear regression per sample as a diagnostic
#'    view and as input to per-group boxplots.
#'
#' 3. Perform the formal group comparison. [test_decay_slopes()] fits a
#'    linear mixed-effects model to all bins across all samples
#'    simultaneously and tests the distance-by-group interaction with
#'    a Kenward-Roger small-sample F test.
#'
#' 4. Assess anchor specificity. [test_decay_significance()] repeats the
#'    mixed-model fit at many random anchor positions across the genome
#'    to build a null distribution of interaction F-statistics, then
#'    reports an empirical p-value for the anchor of scientific interest.
#'
#' @section Input assumptions:
#'
#' The package assumes per-base coverage input in the format produced by
#' breseq during its mutation identification stage, located at
#' `<sample>/08_mutation_identification/Exported.coverage.tab`. Required
#' columns are `position`, `unique_top_cov`, `unique_bot_cov`,
#' `redundant_top_cov`, and `redundant_bot_cov`. Additional columns
#' produced by newer breseq versions (`raw_redundant_top_cov`,
#' `raw_redundant_bot_cov`, `e_value`) are permitted but ignored.
#'
#' @references
#'
#' Deatherage DE, Barrick JE (2014). Identification of mutations in
#' laboratory-evolved microbes from next-generation sequencing data using
#' breseq. *Methods in Molecular Biology* 1151:165-188.
#' \doi{10.1007/978-1-4939-0554-6_12}
#'
#' Bates D, Machler M, Bolker B, Walker S (2015). Fitting linear
#' mixed-effects models using lme4. *Journal of Statistical Software*
#' 67(1):1-48. \doi{10.18637/jss.v067.i01}
#'
#' Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package:
#' Tests in linear mixed effects models. *Journal of Statistical Software*
#' 82(13):1-26. \doi{10.18637/jss.v082.i13}
#'
#' Kenward MG, Roger JH (1997). Small sample inference for fixed effects
#' from restricted maximum likelihood. *Biometrics* 53(3):983-997.
#' \doi{10.2307/2533558}
#'
#' Phipson B, Smyth GK (2010). Permutation p-values should never be zero:
#' calculating exact p-values when permutations are randomly drawn.
#' *Statistical Applications in Genetics and Molecular Biology* 9(1),
#' Article 39. \doi{10.2202/1544-6115.1585}
#'
#' @keywords internal
#' @importFrom utils globalVariables
"_PACKAGE"


## Silence R CMD check notes about data.table non-standard evaluation
## and ggplot2 tidy-eval pronouns. Every column name referenced inside
## `dt[, ...]` expressions must be listed; `.data` is included because
## ggplot2 uses it as a pronoun for column lookup inside aes().
globalVariables(c(
    ".", ".data", "bin", "bin_end", "bin_mid", "bin_start",
    "distance_from_anchor", "distance_scaled", "estimate",
    "group", "log_norm_cov", "mean_cov", "norm_cov",
    "position", "redundant_bot_cov", "redundant_top_cov",
    "sample", "slope", "se_slope", "term", "total_cov",
    "unique_bot_cov", "unique_top_cov",
    "F_stat", "p_value", "null_F", "viral_rpm"
))
