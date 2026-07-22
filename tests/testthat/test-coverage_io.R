## Tests for coverage_io.R.
## Synthetic coverage tables are built on the fly rather than reading
## real data, so tests are portable.

.make_fake_cov <- function(n_positions = 500L, seed = 1L) {
    set.seed(seed)
    data.table::data.table(
        position          = seq_len(n_positions),
        unique_top_cov    = rpois(n_positions, lambda = 20),
        unique_bot_cov    = rpois(n_positions, lambda = 20),
        redundant_top_cov = rpois(n_positions, lambda = 5),
        redundant_bot_cov = rpois(n_positions, lambda = 5)
    )
}

.write_fake_cov <- function(path, n_positions = 500L, seed = 1L) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(.make_fake_cov(n_positions, seed), path, sep = "\t")
    path
}


test_that("read_coverage_tab reads a well-formed file and adds total_cov", {
    tmp <- tempfile(fileext = ".tab")
    on.exit(unlink(tmp), add = TRUE)

    data.table::fwrite(.make_fake_cov(100L), tmp, sep = "\t")
    dt <- read_coverage_tab(tmp)

    expect_s3_class(dt, "data.table")
    expect_true("total_cov" %in% names(dt))
    expect_equal(
        dt$total_cov,
        dt$unique_top_cov + dt$unique_bot_cov +
            dt$redundant_top_cov + dt$redundant_bot_cov
    )
})


test_that("read_coverage_tab errors on missing columns", {
    tmp <- tempfile(fileext = ".tab")
    on.exit(unlink(tmp), add = TRUE)

    bad <- data.table::data.table(position = 1:10, unique_top_cov = 1:10)
    data.table::fwrite(bad, tmp, sep = "\t")

    expect_error(read_coverage_tab(tmp), "missing required column")
})


test_that("bin_coverage produces correct bin coordinates", {
    cov <- .make_fake_cov(2500L)
    cov[, total_cov := unique_top_cov + unique_bot_cov +
        redundant_top_cov + redundant_bot_cov]

    binned <- bin_coverage(cov, bin_size = 1000L, genome_length = 2500L)

    expect_s3_class(binned, "data.table")
    expect_equal(nrow(binned), 3L)
    expect_equal(binned$bin_start, c(1L, 1001L, 2001L))
    expect_equal(binned$bin_end,   c(1000L, 2000L, 2500L))
    expect_equal(binned$bin_mid,   c(500.5, 1500.5, 2250.5))
})


test_that("normalize_coverage divides by the per-sample median", {
    cov <- .make_fake_cov(2000L)
    cov[, total_cov := unique_top_cov + unique_bot_cov +
        redundant_top_cov + redundant_bot_cov]
    binned <- bin_coverage(cov, bin_size = 500L, genome_length = 2000L)

    normed <- normalize_coverage(binned)

    expect_true(all(c("norm_cov", "log_norm_cov") %in% names(normed)))
    expect_equal(
        stats::median(normed$norm_cov, na.rm = TRUE),
        1,
        tolerance = 1e-8
    )
})


test_that("load_sample_coverage reads a file from a direct path", {
    tmp <- tempfile(fileext = ".tab")
    on.exit(unlink(tmp), add = TRUE)
    .write_fake_cov(tmp, n_positions = 5000L, seed = 1L)

    out <- load_sample_coverage(
        path          = tmp,
        sample_name   = "test_sample",
        bin_size      = 1000L,
        genome_length = 5000L
    )

    expect_s3_class(out, "tbl_df")
    expect_true(all(c("sample", "norm_cov", "log_norm_cov") %in%
        names(out)))
    expect_true(all(out$sample == "test_sample"))
})


test_that("load_all_samples reads files from the coverage_file column", {
    tdir <- tempfile()
    dir.create(tdir)
    on.exit(unlink(tdir, recursive = TRUE), add = TRUE)

    f1 <- .write_fake_cov(file.path(tdir, "s1.tab"),
        n_positions = 5000L, seed = 1L)
    f2 <- .write_fake_cov(file.path(tdir, "s2.tab"),
        n_positions = 5000L, seed = 2L)

    meta <- data.frame(
        sample        = c("s1", "s2"),
        group         = c("g1", "g2"),
        coverage_file = c(f1, f2),
        stringsAsFactors = FALSE
    )

    combined <- load_all_samples(
        sample_metadata = meta,
        bin_size        = 1000L,
        genome_length   = 5000L
    )

    expect_s3_class(combined, "tbl_df")
    expect_setequal(unique(combined$sample), c("s1", "s2"))
    expect_true(all(c("group", "norm_cov", "log_norm_cov") %in%
        names(combined)))
    ## coverage_file should NOT propagate into the output
    expect_false("coverage_file" %in% names(combined))
})


test_that("load_all_samples skips missing files with a warning", {
    tdir <- tempfile()
    dir.create(tdir)
    on.exit(unlink(tdir, recursive = TRUE), add = TRUE)

    f1 <- .write_fake_cov(file.path(tdir, "s1.tab"),
        n_positions = 5000L, seed = 1L)

    meta <- data.frame(
        sample        = c("s1", "s2"),
        group         = c("g1", "g2"),
        coverage_file = c(f1, file.path(tdir, "does_not_exist.tab")),
        stringsAsFactors = FALSE
    )

    withCallingHandlers(
        combined <- load_all_samples(
            sample_metadata = meta,
            bin_size        = 1000L,
            genome_length   = 5000L
        ),
        warning = function(w) invokeRestart("muffleWarning")
    )

    expect_setequal(unique(combined$sample), "s1")
})


test_that("load_all_samples errors clearly when coverage_file is missing", {
    meta <- data.frame(sample = "s1", group = "g1")
    expect_error(
        load_all_samples(meta, bin_size = 1000L, genome_length = 5000L),
        "coverage_file"
    )
})
