## =========================================================================
## Sobol low-discrepancy sequence (up to 7 dimensions), for the Tier 2 power
## surface design (A2-4). No Sobol package is installed in the renv library
## and installing one is a dependency decision, so this is a small base-R
## implementation: Gray-code construction with the Joe-Kuo (new-joe-kuo-6)
## direction numbers for dimensions 1 to 7, 30-bit integers. Deterministic,
## unscrambled; the first point (the origin) is dropped.
##
## Validated in tests/testthat/test-a24.R: every dimension is dyadically
## stratified, and the 2-D projections of the first 400 points are spread
## evenly over a 4 x 4 grid.
## =========================================================================

## s (degree), a (polynomial coefficient bits), initial m_1..m_s, per
## dimension 2..7; dimension 1 is the van der Corput sequence.
.sobol_dirs <- list(
    list(s = 1L, a = 0L, m = 1L),
    list(s = 2L, a = 1L, m = c(1L, 3L)),
    list(s = 3L, a = 1L, m = c(1L, 3L, 1L)),
    list(s = 3L, a = 2L, m = c(1L, 1L, 1L)),
    list(s = 4L, a = 1L, m = c(1L, 1L, 3L, 3L)),
    list(s = 4L, a = 4L, m = c(1L, 3L, 5L, 13L))
)
.SOBOL_BITS <- 30L

.sobol_v <- function(dir) {
    B <- .SOBOL_BITS
    if (is.null(dir)) return(2^(B - seq_len(B)))  # van der Corput
    s <- dir$s
    m <- numeric(B)
    m[seq_len(s)] <- dir$m
    for (k in seq.int(s + 1L, B)) {
        mk <- bitwXor(as.integer(2^s * m[k - s]), as.integer(m[k - s]))
        for (i in seq_len(s - 1L)) {
            if (bitwAnd(bitwShiftR(dir$a, s - 1L - i), 1L) == 1L) {
                mk <- bitwXor(mk, as.integer(2^i * m[k - i]))
            }
        }
        m[k] <- mk
    }
    m * 2^(B - seq_len(B))
}

#' First `n` points (origin excluded) of the `d`-dimensional Sobol sequence.
#' @return Numeric matrix `n` x `d` in (0, 1).
sobol_points <- function(n, d) {
    if (d < 1L || d > 7L) stop("`d` must be between 1 and 7.")
    B <- .SOBOL_BITS
    V <- lapply(seq_len(d), function(j) {
        .sobol_v(if (j == 1L) NULL else .sobol_dirs[[j - 1L]])
    })
    x <- integer(d)
    out <- matrix(NA_real_, n, d)
    for (i in seq_len(n)) {
        ## Rightmost zero bit of (i - 1) (1-indexed).
        c <- 1L
        v <- i - 1L
        while (bitwAnd(v, 1L) == 1L) {
            v <- bitwShiftR(v, 1L)
            c <- c + 1L
        }
        for (j in seq_len(d)) {
            x[j] <- bitwXor(x[j], as.integer(V[[j]][c]))
        }
        out[i, ] <- x / 2^B
    }
    out
}
