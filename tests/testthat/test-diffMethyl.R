# Tests for diffMethyl() — Phase 4 differential methylation analysis

# ─── Helpers ─────────────────────────────────────────────────────────────────

.make_dm_data <- function(n_sites = 20L, n_ctrl = 2L, n_treat = 1L) {
    set.seed(99L)
    n_samp  <- n_ctrl + n_treat
    site_keys <- paste0("chr_sim:", seq_len(n_sites) * 100L, ":+:6mA")

    # First half of sites: differentially methylated (ctrl ~0.9, treat ~0.2)
    # Second half: not (both ~0.5)
    n_diff    <- n_sites %/% 2L
    beta_ctrl <- c(rep(0.9, n_diff), rep(0.5, n_sites - n_diff)) +
                 matrix(rnorm(n_ctrl * n_sites, 0, 0.05), nrow = n_sites)
    beta_treat <- c(rep(0.2, n_diff), rep(0.5, n_sites - n_diff)) +
                  rnorm(n_treat * n_sites, 0, 0.05)

    methyl_mat <- cbind(
        matrix(pmax(0, pmin(1, beta_ctrl)), nrow = n_sites, ncol = n_ctrl),
        matrix(pmax(0, pmin(1, beta_treat)), nrow = n_sites, ncol = n_treat)
    )
    colnames(methyl_mat) <- c(
        paste0("ctrl_",  seq_len(n_ctrl)),
        paste0("treat_", seq_len(n_treat))
    )
    rownames(methyl_mat) <- site_keys
    cov_mat <- matrix(30L, nrow = n_sites, ncol = n_samp,
                      dimnames = list(site_keys, colnames(methyl_mat)))

    rd <- S4Vectors::DataFrame(
        chrom    = rep("chr_sim", n_sites),
        position = seq_len(n_sites) * 100L,
        strand   = rep("+", n_sites),
        mod_type = rep("6mA", n_sites),
        is_diff  = c(rep(TRUE, n_diff), rep(FALSE, n_sites - n_diff)),
        row.names = site_keys
    )
    cd <- S4Vectors::DataFrame(
        sample_name = colnames(methyl_mat),
        condition   = c(rep("control", n_ctrl), rep("treatment", n_treat)),
        replicate   = seq_len(n_samp),
        row.names   = colnames(methyl_mat)
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = methyl_mat, coverage = cov_mat),
        rowData = rd,
        colData = cd
    )
    new("commaData", se,
        genomeInfo = c(chr_sim = 100000L),
        annotation = GenomicRanges::GRanges(),
        motifSites = GenomicRanges::GRanges())
}

# ─── Basic functionality ──────────────────────────────────────────────────────

test_that("diffMethyl: returns a commaData object", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition)
    expect_s4_class(dm, "commaData")
})

test_that("diffMethyl: dimension unchanged after running", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition)
    expect_equal(dim(dm), dim(obj))
})

test_that("diffMethyl: rowData gains dm_ result columns", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition)
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    expect_true("dm_pvalue"     %in% colnames(rd))
    expect_true("dm_padj"       %in% colnames(rd))
    expect_true("dm_delta_beta" %in% colnames(rd))
})

test_that("diffMethyl: per-condition mean_beta columns added", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition)
    rd  <- colnames(SummarizedExperiment::rowData(dm))
    expect_true("dm_mean_beta_control"   %in% rd)
    expect_true("dm_mean_beta_treatment" %in% rd)
})

test_that("diffMethyl: metadata records result_cols", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition)
    md  <- S4Vectors::metadata(dm)
    expect_true(!is.null(md$diffMethyl_result_cols))
    expect_true("dm_pvalue" %in% md$diffMethyl_result_cols)
    expect_true("dm_padj"   %in% md$diffMethyl_result_cols)
})

test_that("diffMethyl: metadata params recorded correctly", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "beta_binomial",
                      p_adjust_method = "BH")
    params <- S4Vectors::metadata(dm)$diffMethyl_params
    expect_equal(params$method, "beta_binomial")
    expect_equal(params$p_adjust_method, "BH")
})

test_that("diffMethyl: existing rowData columns preserved", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition)
    rd  <- SummarizedExperiment::rowData(dm)
    expect_true("is_diff" %in% colnames(rd))
    expect_true("mod_type" %in% colnames(rd))
})

# ─── Statistical correctness ──────────────────────────────────────────────────

test_that("diffMethyl: dm_pvalue in [0, 1] for non-NA sites", {
    obj    <- .make_dm_data()
    dm     <- diffMethyl(obj, formula = ~ condition)
    pvals  <- SummarizedExperiment::rowData(dm)$dm_pvalue
    nonNA  <- pvals[!is.na(pvals)]
    expect_true(all(nonNA >= 0 & nonNA <= 1))
})

test_that("diffMethyl: dm_padj >= dm_pvalue for all non-NA sites", {
    obj    <- .make_dm_data()
    dm     <- diffMethyl(obj, formula = ~ condition)
    rd     <- as.data.frame(SummarizedExperiment::rowData(dm))
    ok     <- !is.na(rd$dm_pvalue) & !is.na(rd$dm_padj)
    expect_true(all(rd$dm_padj[ok] >= rd$dm_pvalue[ok] - 1e-10))
})

test_that("diffMethyl: true positive enrichment in ground truth diff sites", {
    # Truly differential sites should have lower padj than non-diff sites
    obj <- .make_dm_data(n_sites = 30L, n_ctrl = 2L, n_treat = 1L)
    dm  <- diffMethyl(obj, formula = ~ condition)
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    median_padj_diff    <- median(rd$dm_padj[rd$is_diff  & !is.na(rd$dm_padj)])
    median_padj_nondiff <- median(rd$dm_padj[!rd$is_diff & !is.na(rd$dm_padj)])
    expect_true(median_padj_diff < median_padj_nondiff)
})

test_that("diffMethyl: delta_beta sign matches direction (treat - ctrl)", {
    obj <- .make_dm_data(n_sites = 10L)
    dm  <- diffMethyl(obj, formula = ~ condition)
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    # First sites are diff (control high, treatment low) → delta_beta < 0
    delta <- rd$dm_delta_beta[1:5]
    expect_true(all(delta[!is.na(delta)] < 0))
})

# ─── mod_type filtering ───────────────────────────────────────────────────────

test_that("diffMethyl: mod_type = '6mA' tests only 6mA sites", {
    data(comma_example_data)
    dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))

    # 6mA sites should have non-NA p-values
    has_6mA_pval  <- !is.na(rd$dm_pvalue[rd$mod_type == "6mA"])
    has_5mC_pval  <- !is.na(rd$dm_pvalue[rd$mod_type == "5mC"])
    expect_true(any(has_6mA_pval))
    expect_true(!any(has_5mC_pval))
})

test_that("diffMethyl: all mod types tested when mod_type = NULL", {
    data(comma_example_data)
    dm <- diffMethyl(comma_example_data, formula = ~ condition)
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))
    expect_true(any(!is.na(rd$dm_pvalue[rd$mod_type == "6mA"])))
    expect_true(any(!is.na(rd$dm_pvalue[rd$mod_type == "5mC"])))
})

# ─── min_coverage filtering ───────────────────────────────────────────────────

test_that("diffMethyl: min_coverage = 1000 → all NA p-values", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, min_coverage = 1000L)
    pvals <- SummarizedExperiment::rowData(dm)$dm_pvalue
    expect_true(all(is.na(pvals)))
})

# ─── p_adjust_method ─────────────────────────────────────────────────────────

test_that("diffMethyl: p_adjust_method = 'bonferroni' is accepted", {
    obj <- .make_dm_data()
    expect_no_error(
        diffMethyl(obj, formula = ~ condition, p_adjust_method = "bonferroni")
    )
})

test_that("diffMethyl: p_adjust_method = 'none' gives padj equal to pvalue", {
    obj  <- .make_dm_data()
    dm   <- diffMethyl(obj, formula = ~ condition, p_adjust_method = "none")
    rd   <- as.data.frame(SummarizedExperiment::rowData(dm))
    ok   <- !is.na(rd$dm_pvalue) & !is.na(rd$dm_padj)
    expect_equal(rd$dm_pvalue[ok], rd$dm_padj[ok], tolerance = 1e-10)
})

# ─── Error handling ───────────────────────────────────────────────────────────

test_that("diffMethyl: error on non-commaData input", {
    expect_error(diffMethyl(data.frame(x = 1)), "'object' must be a commaData")
})

test_that("diffMethyl: error when formula is not a formula", {
    obj <- .make_dm_data()
    expect_error(diffMethyl(obj, formula = "~ condition"), "'formula' must be a formula")
})

test_that("diffMethyl: error when formula variable not in colData", {
    obj <- .make_dm_data()
    expect_error(
        diffMethyl(obj, formula = ~ nonexistent_col),
        "not found in sample"
    )
})

test_that("diffMethyl: error when mod_type not present in object", {
    obj <- .make_dm_data()
    expect_error(diffMethyl(obj, mod_type = "4mC"), "not found in object")
})

test_that("diffMethyl: method argument must be valid", {
    obj <- .make_dm_data()
    expect_error(diffMethyl(obj, method = "bogus"), "'arg' should be one of")
})

test_that("diffMethyl: works with comma_example_data", {
    data(comma_example_data)
    expect_no_error(
        diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    )
})
