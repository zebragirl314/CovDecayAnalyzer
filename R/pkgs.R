#' @import shiny
#' @import bslib
#' @import ggplot2
#' @import data.table
#' @importFrom broom tidy
#' @importFrom rcompanion pairwisePermutationTest
#' @importFrom utils read.csv globalVariables
#' @importFrom stats aov kruskal.test lm median sd shapiro.test

utils::globalVariables(c("bin", "bin_end", "bin_mid", "bin_start", "dist_from_ori", "estimate", "group",
                       "log_norm_cov", "mean_cov", "norm_cov", "position", "redundant_bot_cov",
                       "redundant_top_cov", "slope", "term", "total_cov", "unique_bot_cov", "unique_top_cov", "."))

