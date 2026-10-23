## =========================================================================
## Scenario definitions, scenario IDs, seed derivation, and tier grids.
##
## A scenario is the set of *generative* parameters. Bin size is an
## *analysis* choice, so it is not part of the scenario ID: the same
## simulated atoms are analyzed at 5 kb and 1 kb (paired, PLAN.md 2.1).
## A grid row is a scenario plus a bin size.
## =========================================================================

ATOM_SIZE <- 1000L

## Reference values (PLAN.md section 4.1). These are centers of plausible
## ranges, not estimates for any organism.
scenario_defaults <- function() {
    list(
        L            = 238080,
        s_ref        = 6e-6,
        anchor_frac  = 0.5,
        n            = 5L,        # samples per group; scalar or length G
        G            = 2L,
        Delta        = 0,
        tau          = 2e-6,
        sigma_atom   = 0.15,
        depth_median = 200,
        depth_sd_log = 0.4,
        atom_size    = ATOM_SIZE,
        sens         = list()
    )
}

## Fields that define the data-generating process (and so the scenario ID).
.scenario_id_fields <- c("L", "s_ref", "anchor_frac", "n", "G", "Delta",
    "tau", "sigma_atom", "depth_median", "depth_sd_log", "atom_size")

make_scenario <- function(..., bin_size = 5000L) {
    sc <- utils::modifyList(scenario_defaults(), list(...))
    sc$G <- as.integer(sc$G)
    sc$n <- as.integer(sc$n)
    if (!length(sc$n) %in% c(1L, sc$G)) {
        stop("`n` must have length 1 or length G.")
    }
    if (any(sc$n < 2L)) {
        stop("Every group needs at least 2 samples.")
    }
    if (sc$G < 2L) {
        stop("`G` must be at least 2.")
    }
    if (sc$anchor_frac < 0 || sc$anchor_frac > 1) {
        stop("`anchor_frac` must be in [0, 1].")
    }
    sc$bin_size <- as.integer(bin_size)
    if (sc$bin_size %% sc$atom_size != 0L) {
        stop("`bin_size` must be a multiple of the atom size.")
    }
    sc$anchor <- sc$anchor_frac * sc$L
    sc$scenario_id <- scenario_id(sc)
    class(sc) <- c("cda_scenario", "list")
    sc
}

scenario_id <- function(sc) {
    parts <- vapply(.scenario_id_fields, function(f) {
        paste0(f, "=", paste(format(sc[[f]], digits = 15, trim = TRUE),
            collapse = ","))
    }, character(1L))
    sens <- if (length(sc$sens)) paste(deparse(sc$sens), collapse = "") else ""
    key <- paste(c(parts, paste0("sens=", sens)), collapse = ";")
    paste0("sc", substr(digest::digest(key, algo = "md5", serialize = FALSE),
        1L, 10L))
}

## Seed from a hash of (scenario_id, replicate_id): any replicate is
## regenerable alone, and adding replicates never disturbs existing ones.
derive_seed <- function(scenario_id, replicate_id) {
    h <- digest::digest(paste(scenario_id, replicate_id, sep = "|"),
        algo = "xxhash32", serialize = FALSE)
    max(1L, strtoi(substr(h, 2L, 8L), base = 16L))
}

## Run `expr` under a fixed seed and restore the caller's RNG state, so a
## replicate never perturbs (or depends on) the global stream.
with_seed <- function(seed, expr) {
    has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old <- if (has_seed) get(".Random.seed", envir = globalenv()) else NULL
    on.exit({
        if (is.null(old)) {
            if (exists(".Random.seed", envir = globalenv(),
                inherits = FALSE)) {
                rm(".Random.seed", envir = globalenv())
            }
        } else {
            assign(".Random.seed", old, envir = globalenv())
        }
    })
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
        sample.kind = "Rejection")
    force(expr)
}

## Q9 group-effect pattern: one group differs from the shared reference.
## G = 2: (0, Delta); G = 3: (0, 0, Delta); G = 4: (0, 0, 0, Delta).
group_deltas <- function(G, Delta) {
    c(rep(0, G - 1L), Delta)
}

n_per_group <- function(sc) {
    if (length(sc$n) == 1L) rep(sc$n, sc$G) else sc$n
}

## Effective between-sample slope SD (PLAN.md 2.5), per bp:
## sqrt(tau^2 + sigma_bin^2 / Sxx), Sxx over the bin distances of a sample.
tau_eff <- function(sc, bin_size = sc$bin_size) {
    geo <- bin_geometry(sc$L, bin_size)
    d <- abs(geo$bin_mid - sc$anchor)
    sxx <- sum((d - mean(d))^2)
    sigma_bin <- sc$sigma_atom / sqrt(bin_size / sc$atom_size)
    sqrt(sc$tau^2 + sigma_bin^2 / sxx)
}

## Scenario (as list) from one grid row.
scenario_from_row <- function(row) {
    row <- as.list(row)
    args <- row[intersect(names(row), c(.scenario_id_fields, "bin_size"))]
    args$atom_size <- NULL
    do.call(make_scenario, args)
}

scenario_to_row <- function(sc) {
    tibble::tibble(
        scenario_id  = sc$scenario_id,
        L            = sc$L,
        s_ref        = sc$s_ref,
        anchor_frac  = sc$anchor_frac,
        n            = sc$n[1L],
        G            = sc$G,
        Delta        = sc$Delta,
        tau          = sc$tau,
        sigma_atom   = sc$sigma_atom,
        depth_median = sc$depth_median,
        depth_sd_log = sc$depth_sd_log,
        bin_size     = sc$bin_size,
        tau_eff      = tau_eff(sc),
        delta_std    = sc$Delta / tau_eff(sc)
    )
}

## One row per scenario x bin size from a list of make_scenario() calls.
rows_to_grid <- function(scs) {
    do.call(rbind, lapply(scs, scenario_to_row))
}

## One row per (scenario, bin size, batch) with the replicate range.
expand_batches <- function(grid, R, batch_size = 100L) {
    starts <- seq.int(1L, R, by = batch_size)
    out <- lapply(seq_len(nrow(grid)), function(i) {
        tibble::tibble(
            grid[rep(i, length(starts)), ],
            batch_id  = seq_along(starts),
            rep_start = starts,
            rep_end   = pmin(starts + batch_size - 1L, R)
        )
    })
    do.call(rbind, out)
}

## ---------------------------------------------------------------------
## Tier grids. Only the pilot is defined in A2-1; later tiers are built in
## A2-2 onward (PLAN.md 4.2).
## ---------------------------------------------------------------------

## Tier 0 pilot cells. `purposes` tags which pilot measurements a cell
## feeds; a cell can feed several.
pilot_grid <- function() {
    cells <- list()
    add <- function(purpose, ...) {
        sc <- make_scenario(...)
        cells[[length(cells) + 1L]] <<- list(sc = sc, purpose = purpose)
    }

    ## Timing: smallest design, largest design at 5 kb, and the largest
    ## realistic design at 1 kb.
    add("timing_small", L = 50000, G = 2L, n = 2L, bin_size = 5000L)
    add("timing_large", L = 500000, G = 4L, n = 10L, bin_size = 5000L)
    add("timing_large_1kb", L = 500000, G = 4L, n = 10L, bin_size = 1000L)

    ## Singular and failure rates under the baseline null, reference
    ## geometry, n x G.
    for (G in 2:4) {
        for (n in c(2L, 3L, 5L, 10L)) {
            add("null", G = G, n = n)
        }
    }

    ## M3 preliminary: Type I at three anchor positions.
    for (a in c(0.1, 0.5, 0.9)) {
        add("m3", G = 2L, n = 5L, anchor_frac = a)
    }

    ## Collapse sanity (PLAN.md 2.5), n = 5, G = 2. Pair A varies s_ref,
    ## which does not enter tau_eff in this generator (ratio 1). Pair B
    ## varies tau, which changes tau_eff by a known factor, and rescales
    ## Delta to hold Delta / tau_eff matched.
    ref <- make_scenario(G = 2L, n = 5L)
    deltas_a <- c(3e-6, 4.5e-6)
    for (d in deltas_a) {
        add("collapse_ref", G = 2L, n = 5L, Delta = d)
        add("collapse_sref", G = 2L, n = 5L, Delta = d, s_ref = 1.2e-5)
        tb <- make_scenario(G = 2L, n = 5L, tau = 1e-6)
        d_b <- d * tau_eff(tb) / tau_eff(ref)
        add("collapse_tau", G = 2L, n = 5L, Delta = d_b, tau = 1e-6)
    }

    grid <- rows_to_grid(lapply(cells, `[[`, "sc"))
    grid$purposes <- vapply(cells, `[[`, character(1L), "purpose")
    keep <- !duplicated(grid[c("scenario_id", "bin_size")])
    dup_tags <- tapply(grid$purposes, paste(grid$scenario_id, grid$bin_size),
        function(x) paste(unique(x), collapse = ";"))
    grid <- grid[keep, ]
    grid$purposes <- unname(dup_tags[paste(grid$scenario_id, grid$bin_size)])
    grid
}

build_grid <- function(tier) {
    switch(tier,
        pilot = pilot_grid(),
        stop("Tier '", tier, "' is not defined in A2-1; see PLAN.md 4.2.")
    )
}
