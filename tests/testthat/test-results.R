# Tests for results() and filterResults() — Phase 4

# ─── Helper ──────────────────────────────────────────────────────────────────

.make_tested_object <- function() {
  set.seed(123L)
  n_sites <- 15L
  methyl_mat <- matrix(
    c(rep(0.9, n_sites), rep(0.9, n_sites), rep(0.2, n_sites)),
    nrow = n_sites, ncol = 3L,
    dimnames = list(NULL, c("ctrl_1", "ctrl_2", "treat_1"))
  )
  # Make last 5 sites non-differential
  methyl_mat[(n_sites - 4L):n_sites, 3L] <- 0.9
  cov_mat <- matrix(10L,
    nrow = n_sites, ncol = 3L,
    dimnames = dimnames(methyl_mat)
  )
  site_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr_sim", n_sites),
    ranges = IRanges::IRanges(start = seq_len(n_sites) * 50L, width = 1L),
    strand = rep("+", n_sites),
    mod_type = factor(rep("6mA", n_sites), levels = c("4mC", "5mC", "6mA")),
    motif = rep("GATC", n_sites)
  )
  GenomeInfoDb::seqinfo(site_gr) <- GenomeInfoDb::Seqinfo(
    seqnames = "chr_sim",
    seqlengths = 100000L,
    isCircular = FALSE
  )
  cd <- S4Vectors::DataFrame(
    sample_name = c("ctrl_1", "ctrl_2", "treat_1"),
    condition   = c("control", "control", "treatment"),
    replicate   = 1:3,
    row.names   = c("ctrl_1", "ctrl_2", "treat_1")
  )
  rse <- SummarizedExperiment::SummarizedExperiment(
    assays     = list(methylation = methyl_mat, coverage = cov_mat),
    rowRanges  = site_gr,
    colData    = cd
  )
  obj <- new("commaData", rse)
  diffMethyl(obj, formula = ~condition)
}

# ─── results() ────────────────────────────────────────────────────────────────

test_that("results: returns a data.frame", {
  dm <- .make_tested_object()
  res <- results(dm)
  expect_s3_class(res, "data.frame")
})

test_that("results: one row per site in object", {
  dm <- .make_tested_object()
  res <- results(dm)
  expect_equal(nrow(res), nrow(dm))
})

test_that("results: contains required columns", {
  dm <- .make_tested_object()
  res <- results(dm)
  required <- c(
    "chrom", "position", "strand", "mod_type",
    "dm_pvalue", "dm_padj", "dm_delta_beta"
  )
  expect_true(all(required %in% colnames(res)))
})

test_that("results: contains per-condition mean_beta columns", {
  dm <- .make_tested_object()
  res <- results(dm)
  expect_true("dm_mean_beta_control" %in% colnames(res))
  expect_true("dm_mean_beta_treatment" %in% colnames(res))
})

test_that("results: error when diffMethyl not yet run", {
  data(comma_example_data)
  expect_error(results(comma_example_data), "run diffMethyl")
})

test_that("results: mod_type filter works", {
  data(comma_example_data)
  dm <- diffMethyl(comma_example_data, formula = ~condition)
  res <- results(dm, mod_type = "6mA")
  expect_true(all(res$mod_type == "6mA"))
})

test_that("results: mod_type filter reduces rows", {
  data(comma_example_data)
  dm <- diffMethyl(comma_example_data, formula = ~condition)
  res_all <- results(dm)
  res_6ma <- results(dm, mod_type = "6mA")
  expect_lt(nrow(res_6ma), nrow(res_all))
})

test_that("results: error on invalid mod_type", {
  dm <- .make_tested_object()
  expect_error(results(dm, mod_type = "9mX"), "not found")
})

test_that("results: combined site filters keep original row alignment", {
  skip_if_not_installed("limma")
  obj <- .make_two_modtype_fixture(
    n_6ma = 8L,
    n_5mc = 6L,
    sample_names = c("ctrl_1", "ctrl_2", "treat_1", "treat_2"),
    conditions = c("control", "control", "treatment", "treatment")
  )
  dm <- diffMethyl(obj, formula = ~condition, method = "quasi_f")
  all_res <- results(dm)
  info <- as.data.frame(siteInfo(dm))
  expected <- which(
    info$mod_type == "5mC" &
      info$motif == "CCWGG" &
      info$mod_context == "5mC_CCWGG"
  )

  res <- results(
    dm,
    mod_type = "5mC",
    motif = "CCWGG",
    mod_context = "5mC_CCWGG"
  )

  expect_equal(rownames(res), as.character(expected))
  expect_equal(res$site_key, all_res$site_key[expected])
  expect_equal(res$dm_pvalue, all_res$dm_pvalue[expected])
})

test_that("results: validates filters against the current site subset", {
  skip_if_not_installed("limma")
  obj <- .make_two_modtype_fixture(
    n_6ma = 8L,
    n_5mc = 6L,
    sample_names = c("ctrl_1", "ctrl_2", "treat_1", "treat_2"),
    conditions = c("control", "control", "treatment", "treatment")
  )
  dm <- diffMethyl(obj, formula = ~condition, method = "quasi_f")

  expect_error(
    results(dm, mod_type = "6mA", motif = "CCWGG"),
    "'motif' value\\(s\\) not found"
  )
})

test_that("results: dm_pvalue in [0, 1] for non-NA sites", {
  dm <- .make_tested_object()
  res <- results(dm)
  pv <- res$dm_pvalue[!is.na(res$dm_pvalue)]
  expect_true(all(pv >= 0 & pv <= 1))
})

test_that("results: dm_delta_beta in [-1, 1] for non-NA sites", {
  dm <- .make_tested_object()
  res <- results(dm)
  db <- res$dm_delta_beta[!is.na(res$dm_delta_beta)]
  expect_true(all(db >= -1 & db <= 1))
})

test_that("results: known fixture preserves row order and result values", {
  dm <- .make_tested_object()
  res <- results(dm)

  expect_equal(res$position, seq(50L, 750L, by = 50L))
  expect_equal(as.character(res$chrom), rep("chr_sim", 15L))
  expect_equal(as.character(res$strand), rep("+", 15L))
  expect_equal(as.character(res$mod_type), rep("6mA", 15L))
  expect_equal(res$motif, rep("GATC", 15L))

  expect_equal(res$dm_delta_beta, c(rep(-0.7, 10L), rep(0, 5L)))
  expect_equal(res$dm_mean_beta_control, rep(0.9, 15L))
  expect_equal(
    res$dm_mean_beta_treatment,
    c(rep(0.2, 10L), rep(0.9, 5L))
  )
  expect_true(all(res$dm_padj[seq_len(10L)] < 0.001))
  expect_equal(res$dm_padj[11:15], rep(1, 5L))
})

# ─── filterResults() ──────────────────────────────────────────────────────────

test_that("filterResults: returns a data.frame", {
  dm <- .make_tested_object()
  sig <- filterResults(dm, padj = 1, delta_beta = 0)
  expect_s3_class(sig, "data.frame")
})

test_that("filterResults: subset of results()", {
  dm <- .make_tested_object()
  res <- results(dm)
  sig <- filterResults(dm, padj = 1, delta_beta = 0)
  expect_lte(nrow(sig), nrow(res))
})

test_that("filterResults: all returned sites meet padj threshold", {
  dm <- .make_tested_object()
  sig <- filterResults(dm, padj = 0.5)

  expect_equal(nrow(sig), 10L)
  expect_equal(sig$position, seq(50L, 500L, by = 50L))
  expect_true(all(sig$dm_padj <= 0.5))
})

test_that("filterResults: all returned sites meet delta_beta threshold", {
  dm <- .make_tested_object()
  sig <- filterResults(dm, padj = 1, delta_beta = 0.1)

  expect_equal(nrow(sig), 10L)
  expect_equal(sig$position, seq(50L, 500L, by = 50L))
  expect_equal(sig$dm_delta_beta, rep(-0.7, 10L))
  expect_true(all(abs(sig$dm_delta_beta) >= 0.1))
})

test_that("filterResults: tight thresholds return 0 rows", {
  dm <- .make_tested_object()
  sig <- filterResults(dm, padj = 0, delta_beta = 1)
  expect_equal(nrow(sig), 0L)
})

test_that("filterResults: very loose thresholds return all non-NA sites", {
  dm <- .make_tested_object()
  res <- results(dm)
  n_testable <- sum(!is.na(res$dm_padj) & !is.na(res$dm_delta_beta))
  sig <- filterResults(dm, padj = 1, delta_beta = 0)
  expect_equal(nrow(sig), n_testable)
})

test_that("filterResults: error when diffMethyl not run", {
  data(comma_example_data)
  expect_error(filterResults(comma_example_data), "run diffMethyl")
})

test_that("filterResults: no NA in dm_padj of returned rows", {
  dm <- .make_tested_object()
  sig <- filterResults(dm, padj = 1, delta_beta = 0)
  expect_true(!any(is.na(sig$dm_padj)))
})

test_that("filterResults: no NA in dm_delta_beta of returned rows", {
  dm <- .make_tested_object()
  sig <- filterResults(dm, padj = 1, delta_beta = 0)
  expect_true(!any(is.na(sig$dm_delta_beta)))
})

test_that(
  paste(
    "filterResults: delta_beta=0 threshold retains more sites than",
    "delta_beta=0.5"
  ),
  {
    dm <- .make_tested_object()
    sig_d0 <- filterResults(dm, padj = 1, delta_beta = 0)
    sig_d05 <- filterResults(dm, padj = 1, delta_beta = 0.5)
    expect_gte(nrow(sig_d0), nrow(sig_d05))
  }
)

test_that("filterResults: both thresholds applied simultaneously (AND logic)", {
  dm <- .make_tested_object()
  padj_thresh <- 0.5
  db_thresh <- 0.1
  sig <- filterResults(dm, padj = padj_thresh, delta_beta = db_thresh)

  expect_equal(nrow(sig), 10L)
  expect_equal(sig$position, seq(50L, 500L, by = 50L))
  expect_true(all(sig$dm_padj <= padj_thresh))
  expect_true(all(abs(sig$dm_delta_beta) >= db_thresh))
})

test_that("filterResults: padj=0 returns an empty data frame", {
  dm <- .make_tested_object()
  sig <- filterResults(dm, padj = 0, delta_beta = 0)
  expect_equal(nrow(sig), 0L)
  expect_s3_class(sig, "data.frame")
})

test_that("filterResults: NA thresholds error rather than returning NA rows", {
  dm <- .make_tested_object()
  expect_error(filterResults(dm, padj = NA_real_), "padj")
  expect_error(filterResults(dm, delta_beta = NA_real_), "delta_beta")
})

test_that("filterResults: non-finite or negative thresholds error", {
  dm <- .make_tested_object()
  expect_error(filterResults(dm, padj = Inf), "padj")
  expect_error(filterResults(dm, padj = -0.1), "padj")
  expect_error(filterResults(dm, delta_beta = Inf), "delta_beta")
  expect_error(filterResults(dm, delta_beta = -0.1), "delta_beta")
})

test_that(
  "results: data.frame and GRanges paths align for selected result layer",
  {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 12L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(
      obj,
      formula = ~condition,
      method = "quasi_f",
      result_name = "quasi_f.loose",
      min_coverage = 0L
    )
    dm <- diffMethyl(
      dm,
      formula = ~condition,
      method = "quasi_f",
      result_name = "quasi_f.empty",
      min_coverage = 1000L
    )

    df <- results(dm, name = "quasi_f.loose")
    gr <- results(dm, name = "quasi_f.loose", as = "GRanges")
    result_cols <- c("dm_pvalue", "dm_padj", "dm_delta_beta")
    gr_result_data <- as.data.frame(GenomicRanges::mcols(gr)[, result_cols])
    df_result_data <- df[, result_cols]
    rownames(gr_result_data) <- NULL
    rownames(df_result_data) <- NULL

    expect_equal(as.character(GenomeInfoDb::seqnames(gr)), df$chrom)
    expect_equal(GenomicRanges::start(gr), df$position)
    expect_equal(as.character(GenomicRanges::strand(gr)), df$strand)
    expect_equal(gr_result_data, df_result_data)

    sig <- filterResults(
      dm,
      padj = 1,
      delta_beta = 0,
      name = "quasi_f.loose"
    )
    keep <- !is.na(df$dm_padj) &
      !is.na(df$dm_delta_beta) &
      df$dm_padj <= 1 &
      abs(df$dm_delta_beta) >= 0
    expected <- df[
      keep, ,
      drop = FALSE
    ]

    expect_gt(nrow(sig), 0L)
    expect_equal(sig$position, expected$position)
    sig_result_data <- sig[, result_cols]
    expected_result_data <- expected[, result_cols]
    rownames(sig_result_data) <- NULL
    rownames(expected_result_data) <- NULL
    expect_equal(sig_result_data, expected_result_data)
  }
)
