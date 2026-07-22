## =========================================================================
## Coverage I/O: read per-base coverage, bin, and median-normalize.
##
## The package does not assume any particular on-disk layout for the input
## coverage files. Callers pass the full path to each sample's file
## explicitly, either directly (for one sample) or via a `coverage_file`
## column in a metadata table (for many samples).
##
## Users whose files sit in the layout produced by breseq's mutation
## identification stage can build the `coverage_file` column with the
## helper `breseq_coverage_paths()` in one line.
## =========================================================================

#' Read a per-base coverage file
#'
#' Reads a tab-delimited per-base coverage file and returns a
#' [data.table::data.table] with total per-base coverage computed as
#' the sum of the four strand- and uniqueness-specific coverage columns.
#'
#' @details
#'
#' The file must be tab-delimited with a header row and contain at
#' minimum the columns `position`, `unique_top_cov`, `unique_bot_cov`,
#' `redundant_top_cov`, and `redundant_bot_cov`. Any additional columns
#' are read but ignored.
#'
#' This schema matches the `Exported.coverage.tab` file produced by
#' breseq's mutation identification stage. Files produced by other
#' tools can also be read directly if they use the same column names,
#' or via a shim that renames columns before writing to disk.
#'
#' Total coverage is computed as
#' `unique_top_cov + unique_bot_cov + redundant_top_cov + redundant_bot_cov`.
#'
#' @param path Character scalar. Absolute or relative path to a
#'   coverage file. The file must exist and be readable.
#'
#' @return A [data.table::data.table] with one row per genomic position.
#'   All columns from the input file are preserved. One additional
#'   column is added:
#'   \describe{
#'     \item{`total_cov`}{Integer. Sum of `unique_top_cov`,
#'       `unique_bot_cov`, `redundant_top_cov`, and `redundant_bot_cov`.}
#'   }
#'
#' @references
#' Deatherage DE, Barrick JE (2014). Identification of mutations in
#' laboratory-evolved microbes from next-generation sequencing data
#' using breseq. *Methods in Molecular Biology* 1151:165-188.
#' \doi{10.1007/978-1-4939-0554-6_12}
#'
#' @seealso
#' [bin_coverage()] for the next step in the workflow;
#' [load_sample_coverage()] for a convenience wrapper that runs read,
#' bin, and normalize together;
#' [breseq_coverage_paths()] for building file paths when working with
#' breseq output.
#'
#' @examples
#' ## Read one of the four simulated samples shipped with the package.
#' cov_file <- system.file(
#'     "extdata", "A1", "08_mutation_identification",
#'     "Exported.coverage.tab",
#'     package = "CovDecayAnalyzer"
#' )
#' cov <- read_coverage_tab(cov_file)
#' head(cov)
#' summary(cov$total_cov)
#'
#' @importFrom data.table fread ':='
#' @export
read_coverage_tab <- function(path) {
    if (!is.character(path) || length(path) != 1L || is.na(path)) {
        stop("`path` must be a single non-NA character string.")
    }
    if (!file.exists(path)) {
        stop(sprintf("Coverage file not found: %s", path))
    }

    dt <- data.table::fread(path)

    required_cols <- c(
        "position",
        "unique_top_cov", "unique_bot_cov",
        "redundant_top_cov", "redundant_bot_cov"
    )
    missing_cols <- setdiff(required_cols, names(dt))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "Coverage file %s is missing required column(s): %s",
            path, paste(missing_cols, collapse = ", ")
        ))
    }

    dt[, total_cov := unique_top_cov + unique_bot_cov +
        redundant_top_cov + redundant_bot_cov]
    dt[]
}


#' Bin per-base coverage into fixed-width windows
#'
#' Aggregates per-base coverage into equally-sized non-overlapping bins
#' by taking the mean of `total_cov` within each bin.
#'
#' @details
#'
#' Binning reduces the working table size roughly `bin_size`-fold and
#' smooths local noise while preserving large-scale trends. Bins are
#' constructed by integer division: base at position `p` is assigned to
#' bin `((p - 1) %/% bin_size) + 1`. The final bin is right-truncated
#' at `genome_length`.
#'
#' Smaller bins retain more spatial detail but produce noisier per-bin
#' means. Larger bins are smoother but may average across distinct
#' local features. The default of 1000 bp works well for HCMV-sized
#' genomes (~238 kb) and gives approximately 238 bins per sample.
#'
#' @param cov A [data.table::data.table] as returned by
#'   [read_coverage_tab()]. Must contain columns `position` (positive
#'   integers) and `total_cov` (non-negative numeric).
#' @param bin_size Positive integer. Bin width in base pairs. Default
#'   1000.
#' @param genome_length Positive integer. Reference genome length in
#'   base pairs. Used to clip the final bin end.
#'
#' @return A [data.table::data.table] with one row per bin:
#'   \describe{
#'     \item{`bin`}{Integer. Bin index, 1-based.}
#'     \item{`mean_cov`}{Numeric. Mean total coverage in the bin.}
#'     \item{`bin_start`, `bin_end`}{Integer. Inclusive 1-based
#'       coordinates.}
#'     \item{`bin_mid`}{Numeric. Bin midpoint used as the bin's
#'       representative coordinate for distance calculations.}
#'   }
#'
#' @seealso [normalize_coverage()], [load_sample_coverage()].
#'
#' @examples
#' set.seed(1)
#' cov <- data.table::data.table(
#'     position  = 1:5000,
#'     total_cov = rpois(5000, lambda = 30)
#' )
#' bin_coverage(cov, bin_size = 1000, genome_length = 5000)
#'
#' @importFrom data.table data.table ':='
#' @export
bin_coverage <- function(cov, bin_size = 1000L, genome_length) {
    .validate_positive_integer(bin_size, "bin_size")
    .validate_positive_integer(genome_length, "genome_length")
    if (!inherits(cov, "data.table")) {
        stop("`cov` must be a data.table (e.g. from read_coverage_tab()).")
    }
    if (!all(c("position", "total_cov") %in% names(cov))) {
        stop("`cov` must contain columns `position` and `total_cov`.")
    }

    ## Non-mutating copy so callers can reuse cov downstream.
    dt <- data.table::copy(cov)
    dt[, bin := ((position - 1L) %/% bin_size) + 1L]

    binned <- dt[, list(mean_cov = mean(total_cov)), by = bin]
    binned[, bin_start := (bin - 1L) * bin_size + 1L]
    binned[, bin_end := pmin(bin_start + bin_size - 1L, genome_length)]
    binned[, bin_mid := (bin_start + bin_end) / 2]
    binned[]
}


#' Median-normalize binned coverage
#'
#' Divides each bin's mean coverage by the sample-level median coverage
#' across bins, yielding a scale-invariant `norm_cov` value comparable
#' across samples with different sequencing depths. Also produces
#' `log_norm_cov = log(norm_cov + epsilon)` for downstream log-linear
#' decay modeling.
#'
#' @details
#'
#' Median rather than mean normalization is used because sequencing
#' coverage is overdispersed and can be strongly influenced by
#' high-coverage regions (such as the origin during replication) that
#' would bias the mean. After normalization, `norm_cov` has median
#' exactly 1 within each sample.
#'
#' The `epsilon` addition inside `log()` prevents `log(0)` for
#' zero-coverage bins. The default `1e-5` is well below any meaningful
#' coverage value at typical viral sequencing depths.
#'
#' If the sample's median coverage is zero (essentially no reads
#' recovered), `norm_cov` and `log_norm_cov` are set to `NA_real_` and
#' a warning is issued.
#'
#' @param binned A [data.table::data.table] as returned by
#'   [bin_coverage()]. Must contain a `mean_cov` column.
#' @param epsilon Small positive numeric added inside `log()`.
#'   Default `1e-5`.
#'
#' @return The input table with two columns added:
#'   \describe{
#'     \item{`norm_cov`}{Numeric. `mean_cov` divided by per-sample median.}
#'     \item{`log_norm_cov`}{Numeric. `log(norm_cov + epsilon)`.}
#'   }
#'
#' @seealso [bin_coverage()], [load_sample_coverage()].
#'
#' @examples
#' set.seed(1)
#' binned <- data.table::data.table(
#'     bin       = 1:20,
#'     mean_cov  = rpois(20, lambda = 30),
#'     bin_start = seq(1, by = 1000, length.out = 20),
#'     bin_end   = seq(1000, by = 1000, length.out = 20),
#'     bin_mid   = seq(500.5, by = 1000, length.out = 20)
#' )
#' normalize_coverage(binned)
#'
#' @importFrom data.table ':='
#' @importFrom stats median
#' @export
normalize_coverage <- function(binned, epsilon = 1e-5) {
    if (!inherits(binned, "data.table")) {
        stop("`binned` must be a data.table (e.g. from bin_coverage()).")
    }
    if (!"mean_cov" %in% names(binned)) {
        stop("`binned` must contain a `mean_cov` column.")
    }
    if (!is.numeric(epsilon) || length(epsilon) != 1L ||
        !is.finite(epsilon) || epsilon <= 0) {
        stop("`epsilon` must be a single positive finite number.")
    }

    med <- stats::median(binned$mean_cov, na.rm = TRUE)
    if (!is.finite(med) || med == 0) {
        binned[, norm_cov := NA_real_]
        binned[, log_norm_cov := NA_real_]
        warning("Sample has zero or non-finite median coverage; ",
            "`norm_cov` and `log_norm_cov` set to NA.")
        return(binned[])
    }

    binned[, norm_cov := mean_cov / med]
    binned[, log_norm_cov := log(norm_cov + epsilon)]
    binned[]
}


#' Load, bin, and normalize coverage for a single sample
#'
#' Convenience wrapper that composes [read_coverage_tab()],
#' [bin_coverage()], and [normalize_coverage()] for one coverage file
#' and stamps the sample name onto the result. Returns a tibble
#' suitable for row-binding across many samples.
#'
#' @details
#'
#' This function is layout-agnostic: it reads the file at the path you
#' give it. To process many samples, use [load_all_samples()] which
#' takes a metadata table with a `coverage_file` column and iterates
#' this function under the hood with error tolerance and optional
#' parallelism.
#'
#' @param path Character scalar. Path to the coverage file for this
#'   sample.
#' @param sample_name Character scalar. Value to write into the
#'   returned `sample` column. Used to identify the sample in
#'   downstream analyses.
#' @param bin_size Positive integer. Bin width in base pairs. See
#'   [bin_coverage()]. Default 1000.
#' @param genome_length Positive integer. Reference genome length in
#'   base pairs.
#' @param epsilon Positive numeric. Log-transform offset. See
#'   [normalize_coverage()]. Default `1e-5`.
#'
#' @return A tibble with one row per bin and the columns from
#'   [normalize_coverage()] plus:
#'   \describe{
#'     \item{`sample`}{Character. The `sample_name` argument.}
#'   }
#'
#' @seealso [load_all_samples()] for multi-sample loading.
#'
#' @examples
#' cov_file <- system.file(
#'     "extdata", "A1", "08_mutation_identification",
#'     "Exported.coverage.tab",
#'     package = "CovDecayAnalyzer"
#' )
#' out <- load_sample_coverage(
#'     path          = cov_file,
#'     sample_name   = "A1",
#'     bin_size      = 100,
#'     genome_length = 2000
#' )
#' head(out)
#'
#' @importFrom data.table ':='
#' @importFrom tibble as_tibble
#' @export
load_sample_coverage <- function(path,
                                sample_name,
                                bin_size = 1000L,
                                genome_length,
                                epsilon = 1e-5) {
    .validate_string(path, "path")
    .validate_string(sample_name, "sample_name")

    cov <- read_coverage_tab(path)
    binned <- bin_coverage(cov, bin_size = bin_size,
        genome_length = genome_length)
    binned <- normalize_coverage(binned, epsilon = epsilon)
    binned[, sample := sample_name]
    tibble::as_tibble(binned)
}


#' Load, bin, and normalize coverage for a set of samples
#'
#' Given a sample metadata table with a `coverage_file` column giving
#' the path to each sample's coverage file, loads every sample's file,
#' bins and median-normalizes each sample independently, and returns
#' one combined tibble with sample and group metadata attached. This is
#' the standard entry point to the package.
#'
#' @details
#'
#' The metadata table is the single source of truth for what to load
#' and where each file lives. This decouples the package from any
#' particular on-disk layout: files can live anywhere, be named
#' anything, and come from any tool that writes the required columns
#' (`position`, `unique_top_cov`, `unique_bot_cov`, `redundant_top_cov`,
#' `redundant_bot_cov`).
#'
#' Users with breseq output in the standard layout can build the
#' `coverage_file` column with [breseq_coverage_paths()] in one line
#' (see examples). Users with any other layout can fill the column
#' directly with `file.path()` or `list.files()`.
#'
#' Samples whose file is missing, unreadable, or does not contain the
#' required columns are skipped with a warning rather than aborting
#' the whole load.
#'
#' Parallel processing is delegated to [BiocParallel::bplapply()]. The
#' default [BiocParallel::SerialParam()] processes samples
#' sequentially. To use multiple cores, pass a
#' [BiocParallel::MulticoreParam] or [BiocParallel::SnowParam] instance.
#'
#' @param sample_metadata A data.frame or tibble with, at minimum,
#'   three columns:
#'   \describe{
#'     \item{`sample`}{Character. Sample identifier (must be unique).}
#'     \item{`group`}{Character. Treatment or condition label used in
#'       downstream group comparisons.}
#'     \item{`coverage_file`}{Character. Full path to the sample's
#'       coverage file.}
#'   }
#'   Additional columns are preserved in the output and can be used
#'   for downstream filtering or plotting.
#' @param bin_size Positive integer. Bin width in base pairs. Default
#'   1000.
#' @param genome_length Positive integer. Reference genome length in
#'   base pairs.
#' @param epsilon Positive numeric. Log-transform offset. Default
#'   `1e-5`.
#' @param BPPARAM A `BiocParallel::BiocParallelParam` instance
#'   controlling parallelism. Default [BiocParallel::SerialParam()].
#'
#' @return A tibble with one row per (sample, bin), columns from
#'   [normalize_coverage()] plus:
#'   \describe{
#'     \item{`sample`}{Character. Sample name.}
#'     \item{`group`}{Character. Group label.}
#'     \item{...}{Any additional columns from `sample_metadata` (except
#'       `coverage_file`, which is dropped from the output).}
#'   }
#'
#' @seealso
#' [breseq_coverage_paths()] for building `coverage_file` from a breseq
#' output tree;
#' [add_distance_from_anchor()] and [fit_decay_slopes()] for the next
#' steps of the workflow.
#'
#' @examples
#' ## Example 1: breseq output in standard layout, using the shipped
#' ## example data.
#' extdata_dir <- system.file("extdata", package = "CovDecayAnalyzer")
#' meta <- read.csv(
#'     file.path(extdata_dir, "sample_metadata.csv"),
#'     stringsAsFactors = FALSE
#' )
#' meta$coverage_file <- breseq_coverage_paths(
#'     sample_names = meta$sample,
#'     base_dir     = extdata_dir
#' )
#' combined <- load_all_samples(
#'     sample_metadata = meta,
#'     bin_size        = 100,
#'     genome_length   = 2000
#' )
#' head(combined)
#'
#' ## Example 2: arbitrary layout, paths filled in directly.
#' custom_meta <- data.frame(
#'     sample        = meta$sample,
#'     group         = meta$group,
#'     coverage_file = meta$coverage_file,
#'     stringsAsFactors = FALSE
#' )
#' combined2 <- load_all_samples(
#'     sample_metadata = custom_meta,
#'     bin_size        = 100,
#'     genome_length   = 2000
#' )
#' identical(combined, combined2)
#'
#' @importFrom data.table rbindlist as.data.table
#' @importFrom tibble as_tibble
#' @importFrom BiocParallel bplapply SerialParam
#' @export
load_all_samples <- function(sample_metadata,
                            bin_size = 1000L,
                            genome_length,
                            epsilon = 1e-5,
                            BPPARAM = BiocParallel::SerialParam()) {
    if (!is.data.frame(sample_metadata)) {
        stop("`sample_metadata` must be a data.frame or tibble.")
    }
    required_meta <- c("sample", "group", "coverage_file")
    missing_meta <- setdiff(required_meta, names(sample_metadata))
    if (length(missing_meta) > 0L) {
        stop(sprintf(
            "`sample_metadata` is missing required column(s): %s. ",
            paste(missing_meta, collapse = ", ")),
            "For breseq output, build `coverage_file` with ",
            "`breseq_coverage_paths()`.")
    }

    meta_dt <- data.table::as.data.table(sample_metadata)
    if (anyDuplicated(meta_dt$sample) > 0L) {
        stop("`sample_metadata$sample` contains duplicated values.")
    }

    ## Internal per-sample calls: work directly with data.table to
    ## avoid the per-sample tibble round-trip that load_sample_coverage
    ## would trigger.
    loaded <- BiocParallel::bplapply(
        seq_len(nrow(meta_dt)),
        function(i) {
            sname <- meta_dt$sample[i]
            fpath <- meta_dt$coverage_file[i]
            tryCatch(
                {
                    cov <- read_coverage_tab(fpath)
                    binned <- bin_coverage(cov, bin_size = bin_size,
                        genome_length = genome_length)
                    binned <- normalize_coverage(binned, epsilon = epsilon)
                    binned[, sample := sname]
                    binned
                },
                error = function(e) {
                    warning(sprintf(
                        "Sample '%s' skipped: %s",
                        sname, conditionMessage(e)
                    ))
                    NULL
                }
            )
        },
        BPPARAM = BPPARAM
    )

    loaded <- loaded[!vapply(loaded, is.null, logical(1L))]
    if (length(loaded) == 0L) {
        stop("No samples produced usable coverage data. ",
            "Check the paths in `sample_metadata$coverage_file`.")
    }

    combined <- data.table::rbindlist(loaded, use.names = TRUE, fill = TRUE)

    ## Drop the coverage_file column from the metadata before merging;
    ## it is not useful downstream and would clutter the returned table.
    meta_join <- meta_dt[, setdiff(names(meta_dt), "coverage_file"),
        with = FALSE]
    combined <- merge(combined, meta_join, by = "sample",
        all.x = TRUE, sort = FALSE)
    tibble::as_tibble(combined)
}
