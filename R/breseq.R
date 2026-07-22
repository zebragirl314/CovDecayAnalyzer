## =========================================================================
## breseq-specific helpers.
##
## The rest of the package is layout-agnostic. This file contains the only
## code that assumes anything about breseq's on-disk output structure, so
## the coupling to breseq is contained and easy to audit or replace.
## =========================================================================


#' Build coverage file paths for a breseq output tree
#'
#' Convenience helper that constructs full paths to per-sample
#' `Exported.coverage.tab` files under a directory of breseq output.
#' The result is intended to be assigned to the `coverage_file` column
#' of a metadata table before passing to [load_all_samples()].
#'
#' @details
#'
#' In a standard breseq run, each sample's mutation identification
#' outputs are written to
#' `<base_dir>/<sample>/08_mutation_identification/`, with the
#' per-base coverage table at `Exported.coverage.tab` within that
#' directory. This function assembles those paths.
#'
#' If your breseq output uses non-default subdirectory or filename
#' conventions (older or forked breseq versions, custom pipelines),
#' override `subdir` or `filename` accordingly.
#'
#' If you do not have breseq output at all, you do not need this
#' function: fill the `coverage_file` column of your metadata table
#' by any other means (for example, `list.files()`, `file.path()`, or
#' a `dplyr::mutate()` call) and pass the metadata to
#' [load_all_samples()].
#'
#' Paths are constructed as strings and returned without checking that
#' the files actually exist. Missing files are reported by
#' [load_all_samples()] when it tries to read them.
#'
#' @param sample_names Character vector. Sample identifiers, matching
#'   the directory names immediately under `base_dir`.
#' @param base_dir Character scalar. Root directory containing
#'   per-sample subdirectories.
#' @param subdir Character scalar. Subdirectory within each sample
#'   folder containing the coverage file. Default
#'   `"08_mutation_identification"`.
#' @param filename Character scalar. Name of the coverage file within
#'   `subdir`. Default `"Exported.coverage.tab"`.
#'
#' @return A character vector of paths, one per element of
#'   `sample_names`, in the same order.
#'
#' @seealso [load_all_samples()], [read_coverage_tab()].
#'
#' @examples
#' ## Using the shipped simulated example data.
#' extdata_dir <- system.file("extdata", package = "CovDecayAnalyzer")
#' breseq_coverage_paths(
#'     sample_names = c("A1", "A2", "B1", "B2"),
#'     base_dir     = extdata_dir
#' )
#'
#' ## Typical usage: build the `coverage_file` column of a metadata
#' ## table before calling load_all_samples().
#' meta <- data.frame(
#'     sample = c("A1", "A2", "B1", "B2"),
#'     group  = c("WT-GCV", "WT-GCV", "WT+GCV", "WT+GCV"),
#'     stringsAsFactors = FALSE
#' )
#' meta$coverage_file <- breseq_coverage_paths(meta$sample, extdata_dir)
#' head(meta)
#'
#' @export
breseq_coverage_paths <- function(sample_names,
                                base_dir,
                                subdir = "08_mutation_identification",
                                filename = "Exported.coverage.tab") {
    if (!is.character(sample_names) || length(sample_names) == 0L ||
        anyNA(sample_names)) {
        stop("`sample_names` must be a non-empty character vector ",
            "with no NA values.")
    }
    .validate_string(base_dir, "base_dir")
    .validate_string(subdir, "subdir")
    .validate_string(filename, "filename")

    file.path(base_dir, sample_names, subdir, filename)
}
