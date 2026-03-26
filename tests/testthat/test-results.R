# Tests for results() and filterResults() — Phase 4

# ─── Helper ──────────────────────────────────────────────────────────────────

.make_tested_object <- function() {
    set.seed(123L)
    n_sites  <- 15L
    site_keys <- paste0("chr_sim:", seq_len(n_sites) * 50L, ":+:6mA:GATC")
    methyl_mat <- matrix(
        c(rep(0.9, n_sites), rep(0.9, n_sites), rep(0.2, n_sites)),
        nrow = n_sites, ncol = 3L,
        dimnames = list(site_keys, c("ctrl_1", "ctrl_2", "treat_1"))
    )
    # Make last 5 sites non-differential
    methyl_mat[(n_sites - 4L):n_sites, 3L] <- 0.9
    cov_mat <- matrix(10L, nrow = n_sites, ncol = 3L,
                      dimnames = dimnames(methyl_mat))
    rd <- S4Vectors::DataFrame(
        chrom    = rep("chr_sim", n_sites),
        position = seq_len(n_sites) * 50L,
        strand   = rep("+", n_sites),
        mod_type = rep("6mA", n_sites),
        motif    = rep("GATC", n_sites),
        row.names = site_keys
    )
    cd <- S4Vectors::DataFrame(
        sample_name = c("ctrl_1", "ctrl_2", "treat_1"),
        condition   = c("control", "control", "treatment"),
        replicate   = 1:3,
        row.names   = c("ctrl_1", "ctrl_2", "treat_1")
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = methyl_mat, coverage = cov_mat),
        rowData = rd,
        colData = cd
    )
    obj <- new("commaData", se,
               genomeInfo = c(chr_sim = 100000L),
               annotation = GenomicRanges::GRanges(),
               motifSites = GenomicRanges::GRanges())
    diffMethyl(obj, formula = ~ condition)
}

# ─── results() ────────────────────────────────────────────────────────────────

test_that("results: returns a data.frame", {
    dm  <- .make_tested_object()
    res <- results(dm)
    expect_s3_class(res, "data.frame")
})

test_that("results: one row per site in object", {
    dm  <- .make_tested_object()
    res <- results(dm)
    expect_equal(nrow(res), nrow(dm))
})

test_that("results: contains required columns", {
    dm   <- .make_tested_object()
    res  <- results(dm)
    required <- c("chrom", "position", "strand", "mod_type",
                  "dm_pvalue", "dm_padj", "dm_delta_beta")
    expect_true(all(required %in% colnames(res)))
})

test_that("results: contains per-condition mean_beta columns", {
    dm  <- .make_tested_object()
    res <- results(dm)
    expect_true("dm_mean_beta_control"   %in% colnames(res))
    expect_true("dm_mean_beta_treatment" %in% colnames(res))
})

test_that("results: error when diffMethyl not yet run", {
    data(comma_example_data)
    expect_error(results(comma_example_data), "run diffMethyl")
})

test_that("results: mod_type filter works", {
    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition)
    res <- results(dm, mod_type = "6mA")
    expect_true(all(res$mod_type == "6mA"))
})

test_that("results: mod_type filter reduces rows", {
    data(comma_example_data)
    dm      <- diffMethyl(comma_example_data, formula = ~ condition)
    res_all <- results(dm)
    res_6mA <- results(dm, mod_type = "6mA")
    expect_lt(nrow(res_6mA), nrow(res_all))
})

test_that("results: error on invalid mod_type", {
    dm <- .make_tested_object()
    expect_error(results(dm, mod_type = "9mX"), "not found")
})

test_that("results: dm_pvalue in [0, 1] for non-NA sites", {
    dm  <- .make_tested_object()
    res <- results(dm)
    pv  <- res$dm_pvalue[!is.na(res$dm_pvalue)]
    expect_true(all(pv >= 0 & pv <= 1))
})

test_that("results: dm_delta_beta in [-1, 1] for non-NA sites", {
    dm  <- .make_tested_object()
    res <- results(dm)
    db  <- res$dm_delta_beta[!is.na(res$dm_delta_beta)]
    expect_true(all(db >= -1 & db <= 1))
})

# ─── filterResults() ──────────────────────────────────────────────────────────

test_that("filterResults: returns a data.frame", {
    dm  <- .make_tested_object()
    sig <- filterResults(dm, padj = 1, delta_beta = 0)
    expect_s3_class(sig, "data.frame")
})

test_that("filterResults: subset of results()", {
    dm  <- .make_tested_object()
    res <- results(dm)
    sig <- filterResults(dm, padj = 1, delta_beta = 0)
    expect_lte(nrow(sig), nrow(res))
})

test_that("filterResults: all returned sites meet padj threshold", {
    dm  <- .make_tested_object()
    sig <- filterResults(dm, padj = 0.5)
    if (nrow(sig) > 0) {
        expect_true(all(sig$dm_padj <= 0.5))
    }
})

test_that("filterResults: all returned sites meet delta_beta threshold", {
    dm  <- .make_tested_object()
    sig <- filterResults(dm, padj = 1, delta_beta = 0.1)
    if (nrow(sig) > 0) {
        expect_true(all(abs(sig$dm_delta_beta) >= 0.1))
    }
})

test_that("filterResults: tight thresholds return 0 rows", {
    dm  <- .make_tested_object()
    sig <- filterResults(dm, padj = 0, delta_beta = 1)
    expect_equal(nrow(sig), 0L)
})

test_that("filterResults: very loose thresholds return all non-NA sites", {
    dm  <- .make_tested_object()
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
    dm  <- .make_tested_object()
    sig <- filterResults(dm, padj = 1, delta_beta = 0)
    expect_true(!any(is.na(sig$dm_padj)))
})

test_that("filterResults: no NA in dm_delta_beta of returned rows", {
    dm  <- .make_tested_object()
    sig <- filterResults(dm, padj = 1, delta_beta = 0)
    expect_true(!any(is.na(sig$dm_delta_beta)))
})

test_that("filterResults: delta_beta=0 threshold retains more sites than delta_beta=0.5", {
    dm      <- .make_tested_object()
    sig_d0  <- filterResults(dm, padj = 1, delta_beta = 0)
    sig_d05 <- filterResults(dm, padj = 1, delta_beta = 0.5)
    expect_gte(nrow(sig_d0), nrow(sig_d05))
})

test_that("filterResults: both thresholds applied simultaneously (AND logic)", {
    dm <- .make_tested_object()
    padj_thresh <- 0.5
    db_thresh   <- 0.1
    sig <- filterResults(dm, padj = padj_thresh, delta_beta = db_thresh)
    if (nrow(sig) > 0) {
        expect_true(all(sig$dm_padj       <= padj_thresh))
        expect_true(all(abs(sig$dm_delta_beta) >= db_thresh))
    }
})

test_that("filterResults: padj=0 returns an empty data frame", {
    dm  <- .make_tested_object()
    sig <- filterResults(dm, padj = 0, delta_beta = 0)
    expect_equal(nrow(sig), 0L)
    expect_s3_class(sig, "data.frame")
})
