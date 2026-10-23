## KW on per-sample slopes with chi-square p-value (the traditional route),
## plus Dunn-Holm ungated and gated on KW p < 0.05. Returns 3 rows:
## arm = "kw_chisq", "kw_chisq_dunn_ungated", "kw_chisq_dunn_gated".
analyze_kw_chisq <- function(data,
                             scenario_id = attr(data, "scenario_id"),
                             replicate_id = attr(data, "replicate_id"),
                             bin_size = attr(data, "bin_size")) {
    ensure_package()
    kw_fun <- function(slope, group) {
        k <- stats::kruskal.test(slope, group)
        list(p = k$p.value, H = unname(k$statistic))
    }
    .kw_rows("kw_chisq", data, kw_fun, scenario_id, replicate_id, bin_size)
}
