# Tests for diffMethyl() — Phase 4 differential methylation analysis

# ─── Helpers ─────────────────────────────────────────────────────────────────

.make_dm_data <- function(n_sites = 20L, n_ctrl = 2L, n_treat = 1L) {
    .make_diff_methyl_fixture(
        n_sites = n_sites,
        n_ctrl = n_ctrl,
        n_treat = n_treat
    )
}

# ─── Basic functionality ──────────────────────────────────────────────────────

test_that("diffMethyl: returns a commaData object", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    expect_s4_class(dm, "commaData")
})

test_that("diffMethyl: dimension unchanged after running", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    expect_equal(dim(dm), dim(obj))
})

test_that("diffMethyl: rowData gains dm_ result columns", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    expect_true("dm_pvalue"     %in% colnames(rd))
    expect_true("dm_padj"       %in% colnames(rd))
    expect_true("dm_delta_beta" %in% colnames(rd))
})

test_that("diffMethyl: per-condition mean_beta columns added", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd  <- colnames(SummarizedExperiment::rowData(dm))
    expect_true("dm_mean_beta_control"   %in% rd)
    expect_true("dm_mean_beta_treatment" %in% rd)
})

test_that("diffMethyl: metadata records result_cols", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    md  <- S4Vectors::metadata(dm)
    expect_true(!is.null(md$diffMethyl_result_cols))
    expect_true("dm_pvalue" %in% md$diffMethyl_result_cols)
    expect_true("dm_padj"   %in% md$diffMethyl_result_cols)
})

test_that("diffMethyl: metadata params recorded correctly", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f",
                      p_adjust_method = "BH")
    params <- S4Vectors::metadata(dm)$diffMethyl_params
    expect_equal(params$method, "quasi_f")
    expect_equal(params$p_adjust_method, "BH")
})

test_that("diffMethyl: existing rowData columns preserved", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd  <- SummarizedExperiment::rowData(dm)
    expect_true("is_diff" %in% colnames(rd))
    expect_true("mod_type" %in% colnames(rd))
})

# ─── Statistical correctness ──────────────────────────────────────────────────

test_that("diffMethyl: dm_pvalue in [0, 1] for non-NA sites", {
    obj    <- .make_dm_data()
    dm     <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    pvals  <- SummarizedExperiment::rowData(dm)$dm_pvalue
    nonNA  <- pvals[!is.na(pvals)]
    expect_true(all(nonNA >= 0 & nonNA <= 1))
})

test_that("diffMethyl: dm_padj >= dm_pvalue for all non-NA sites", {
    obj    <- .make_dm_data()
    dm     <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd     <- as.data.frame(SummarizedExperiment::rowData(dm))
    ok     <- !is.na(rd$dm_pvalue) & !is.na(rd$dm_padj)
    expect_true(all(rd$dm_padj[ok] >= rd$dm_pvalue[ok] - 1e-10))
})

test_that("diffMethyl: true positive enrichment in ground truth diff sites", {
    # Truly differential sites should have lower padj than non-diff sites
    obj <- .make_dm_data(n_sites = 30L, n_ctrl = 2L, n_treat = 1L)
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    median_padj_diff    <- median(rd$dm_padj[rd$is_diff  & !is.na(rd$dm_padj)])
    median_padj_nondiff <- median(rd$dm_padj[!rd$is_diff & !is.na(rd$dm_padj)])
    expect_true(median_padj_diff < median_padj_nondiff)
})

test_that("diffMethyl: delta_beta sign matches direction (treat - ctrl)", {
    obj <- .make_dm_data(n_sites = 10L)
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    # First sites are diff (control high, treatment low) → delta_beta < 0
    delta <- rd$dm_delta_beta[1:5]
    expect_true(all(delta[!is.na(delta)] < 0))
})

test_that("diffMethyl: limma uses observed count assays when present", {
    skip_if_not_installed("limma")
    n_sites <- 12L
    sample_names <- c("ctrl_1", "ctrl_2", "treat_1", "treat_2")
    beta <- matrix(
        0.5,
        nrow = n_sites,
        ncol = length(sample_names),
        dimnames = list(NULL, sample_names)
    )
    coverage <- matrix(
        20L,
        nrow = n_sites,
        ncol = length(sample_names),
        dimnames = dimnames(beta)
    )
    mod_counts <- matrix(
        rep(c(2L, 2L, 18L, 18L), each = n_sites),
        nrow = n_sites,
        dimnames = dimnames(beta)
    )
    canonical_counts <- coverage - mod_counts
    sample_info <- data.frame(
        sample_name = sample_names,
        condition = c("control", "control", "treatment", "treatment"),
        replicate = c(1L, 2L, 1L, 2L),
        stringsAsFactors = FALSE
    )
    obj <- .make_commaData_fixture(
        beta = beta,
        coverage = coverage,
        sample_info = sample_info,
        positions = seq_len(n_sites) * 100L,
        mod_counts = mod_counts,
        canonical_counts = canonical_counts
    )
    legacy <- obj
    SummarizedExperiment::assays(legacy) <-
        SummarizedExperiment::assays(legacy)[c("methylation", "coverage")]

    dm_observed <- suppressWarnings(
        diffMethyl(obj, formula = ~ condition, method = "limma")
    )
    dm_legacy <- suppressWarnings(
        diffMethyl(legacy, formula = ~ condition, method = "limma")
    )
    p_observed <- SummarizedExperiment::rowData(dm_observed)$dm_pvalue
    p_legacy <- SummarizedExperiment::rowData(dm_legacy)$dm_pvalue

    expect_true(any(!is.na(p_observed)))
    expect_true(median(p_observed, na.rm = TRUE) < 0.05)
    expect_true(all(is.na(p_legacy) | p_legacy > 0.5))
})

# ─── mod_type filtering ───────────────────────────────────────────────────────

test_that("diffMethyl: mod_type = '6mA' tests only 6mA sites", {
    data(comma_example_data)
    dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA",
                     method = "quasi_f")
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))

    # 6mA sites should have non-NA p-values
    has_6mA_pval  <- !is.na(rd$dm_pvalue[rd$mod_type == "6mA"])
    has_5mC_pval  <- !is.na(rd$dm_pvalue[rd$mod_type == "5mC"])
    expect_true(any(has_6mA_pval))
    expect_true(!any(has_5mC_pval))
})

test_that("diffMethyl: all mod types tested when mod_type = NULL", {
    data(comma_example_data)
    dm <- diffMethyl(comma_example_data, formula = ~ condition, method = "quasi_f")
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))
    expect_true(any(!is.na(rd$dm_pvalue[rd$mod_type == "6mA"])))
    expect_true(any(!is.na(rd$dm_pvalue[rd$mod_type == "5mC"])))
})

# ─── min_coverage filtering ───────────────────────────────────────────────────

test_that("diffMethyl: min_coverage = 1000 → all NA p-values", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, min_coverage = 1000L,
                      method = "quasi_f")
    pvals <- SummarizedExperiment::rowData(dm)$dm_pvalue
    expect_true(all(is.na(pvals)))
})

# ─── p_adjust_method ─────────────────────────────────────────────────────────

test_that("diffMethyl: p_adjust_method = 'bonferroni' produces padj >= pvalue", {
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f",
                      p_adjust_method = "bonferroni")
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    ok  <- !is.na(rd$dm_pvalue) & !is.na(rd$dm_padj)
    # Guard against vacuous pass: must have at least one non-NA row
    expect_true(any(ok))
    # Bonferroni correction should always produce padj >= pvalue
    expect_true(all(rd$dm_padj[ok] >= rd$dm_pvalue[ok] - 1e-10))
    # Bonferroni should be more conservative than BH (larger padj)
    dm_bh <- diffMethyl(obj, formula = ~ condition, method = "quasi_f",
                        p_adjust_method = "BH")
    rd_bh <- as.data.frame(SummarizedExperiment::rowData(dm_bh))
    ok_both <- ok & !is.na(rd_bh$dm_padj)
    # Guard against vacuous pass
    expect_true(any(ok_both))
    expect_true(all(rd$dm_padj[ok_both] >= rd_bh$dm_padj[ok_both] - 1e-10))
})

test_that("diffMethyl: p_adjust_method = 'none' gives padj equal to pvalue", {
    obj  <- .make_dm_data()
    dm   <- diffMethyl(obj, formula = ~ condition, method = "quasi_f",
                       p_adjust_method = "none")
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

test_that("diffMethyl: rejects multi-factor formulas before v1", {
    obj <- .make_dm_data(n_ctrl = 2L, n_treat = 2L)
    SummarizedExperiment::colData(obj)$batch <- c("b1", "b2", "b1", "b2")
    expect_error(
        diffMethyl(obj, formula = ~ condition + batch, method = "quasi_f"),
        "exactly one two-level RHS variable"
    )
})

test_that("diffMethyl: rejects interaction formulas before v1", {
    obj <- .make_dm_data(n_ctrl = 2L, n_treat = 2L)
    SummarizedExperiment::colData(obj)$batch <- c("b1", "b2", "b1", "b2")
    expect_error(
        diffMethyl(obj, formula = ~ condition:batch, method = "quasi_f"),
        "exactly one two-level RHS variable"
    )
})

test_that("diffMethyl: rejects transformed RHS terms before v1", {
    obj <- .make_dm_data(n_ctrl = 2L, n_treat = 2L)
    SummarizedExperiment::colData(obj)$batch <- c(1, 2, 1, 2)
    expect_error(
        diffMethyl(obj, formula = ~ factor(batch), method = "quasi_f"),
        "untransformed two-level RHS variable"
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

test_that("diffMethyl: returns commaData with result columns for comma_example_data", {
    data(comma_example_data)
    dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA",
                     method = "quasi_f")
    expect_s4_class(dm, "commaData")
    rd <- colnames(SummarizedExperiment::rowData(dm))
    expect_true(all(c("dm_pvalue", "dm_padj", "dm_delta_beta") %in% rd))
    # All 393 6mA sites should have non-NA results
    rd_df <- as.data.frame(SummarizedExperiment::rowData(dm))
    n_6ma <- sum(rd_df$mod_type == "6mA")
    n_6ma_result <- sum(rd_df$mod_type == "6mA" & !is.na(rd_df$dm_pvalue))
    expect_equal(n_6ma_result, n_6ma)
    # 5mC sites should have NA results (not tested)
    n_5mc_na <- sum(rd_df$mod_type == "5mC" & is.na(rd_df$dm_pvalue))
    expect_equal(n_5mc_na, sum(rd_df$mod_type == "5mC"))
})

# ─── methylKit method ─────────────────────────────────────────────────────────

test_that("diffMethyl: method='methylkit' errors with informative message if methylKit absent", {
    skip_if(requireNamespace("methylKit", quietly = TRUE),
            "methylKit is installed; skipping absent-package test")
    obj <- .make_dm_data()
    expect_error(
        diffMethyl(obj, formula = ~ condition, method = "methylkit"),
        "methylKit"
    )
})

test_that("diffMethyl: method='methylkit' message uses actual level names", {
    skip_if_not(requireNamespace("methylKit", quietly = TRUE),
                "methylKit not installed")
    obj <- .make_dm_data()
    msgs <- character(0)
    withCallingHandlers(
        diffMethyl(obj, formula = ~ condition, method = "methylkit"),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    )
    # The replacement message should name the actual levels
    expect_true(any(grepl("treatment", msgs, fixed = TRUE)))
    expect_true(any(grepl("control", msgs, fixed = TRUE)))
    # The generic methylKit "group: 0" / "group: 1" message must NOT appear
    expect_false(any(grepl("group: 0", msgs, fixed = TRUE)))
    expect_false(any(grepl("group: 1", msgs, fixed = TRUE)))
})

# ─── limma method ─────────────────────────────────────────────────────────────

test_that("diffMethyl: method='limma' returns commaData with correct columns", {
    skip_if_not_installed("limma")
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "limma")
    expect_s4_class(dm, "commaData")
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))
    expect_true(all(c("dm_pvalue", "dm_padj", "dm_delta_beta",
                      "dm_mean_beta_control", "dm_mean_beta_treatment") %in%
                        colnames(rd)))
    expect_equal(nrow(rd), nrow(SummarizedExperiment::rowData(obj)))
})

test_that("diffMethyl: method='limma' produces valid p-values in [0, 1]", {
    skip_if_not_installed("limma")
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "limma")
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    pvals <- rd$dm_pvalue[!is.na(rd$dm_pvalue)]
    expect_true(length(pvals) > 0)
    expect_true(all(pvals >= 0 & pvals <= 1))
    expect_true(all(rd$dm_padj[!is.na(rd$dm_padj)] >= rd$dm_pvalue[!is.na(rd$dm_pvalue)]))
})

test_that("diffMethyl: limma and quasi_f delta_beta values are highly correlated", {
    skip_if_not_installed("limma")
    obj  <- .make_dm_data(n_sites = 40L)
    dm_l <- diffMethyl(obj, formula = ~ condition, method = "limma")
    dm_q <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    db_l <- SummarizedExperiment::rowData(dm_l)$dm_delta_beta
    db_q <- SummarizedExperiment::rowData(dm_q)$dm_delta_beta
    ok   <- !is.na(db_l) & !is.na(db_q)
    expect_true(sum(ok) > 0)
    expect_gt(cor(db_l[ok], db_q[ok]), 0.95)
})

test_that("diffMethyl: method='limma' records alpha in metadata", {
    skip_if_not_installed("limma")
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "limma", alpha = 1.0)
    expect_equal(S4Vectors::metadata(dm)$diffMethyl_params$alpha, 1.0)
})

test_that("diffMethyl: method='limma' errors with informative message if limma absent", {
    skip_if(requireNamespace("limma", quietly = TRUE),
            "limma is installed; skipping absent-package test")
    obj <- .make_dm_data()
    expect_error(
        diffMethyl(obj, formula = ~ condition, method = "limma"),
        "limma"
    )
})

test_that("diffMethyl: non-positive alpha errors informatively", {
    skip_if_not_installed("limma")
    obj <- .make_dm_data()
    expect_error(
        diffMethyl(obj, formula = ~ condition, method = "limma", alpha = 0),
        "alpha"
    )
    expect_error(
        diffMethyl(obj, formula = ~ condition, method = "limma", alpha = -1),
        "alpha"
    )
})

# ─── quasi_f method ───────────────────────────────────────────────────────────

test_that("diffMethyl: quasi_f recovers majority of ground-truth diff sites", {
    skip_if_not_installed("limma")
    data(comma_example_data)
    dm_q  <- diffMethyl(comma_example_data, formula = ~ condition,
                        method = "quasi_f", mod_type = "6mA")
    rd    <- as.data.frame(SummarizedExperiment::rowData(dm_q))
    rd6   <- rd[rd$mod_type == "6mA", ]
    n_true_diff <- sum(rd6$is_diff, na.rm = TRUE)
    n_detected  <- sum(rd6$is_diff & !is.na(rd6$dm_pvalue) & rd6$dm_pvalue < 0.2,
                       na.rm = TRUE)
    expect_gte(n_detected, floor(n_true_diff * 0.5))
})

test_that("diffMethyl: method='quasi_f' records method in metadata", {
    skip_if_not_installed("limma")
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    expect_equal(S4Vectors::metadata(dm)$diffMethyl_params$method, "quasi_f")
})

test_that("diffMethyl: method='quasi_f' errors informatively if limma absent", {
    skip_if(requireNamespace("limma", quietly = TRUE),
            "limma is installed; skipping absent-package test")
    obj <- .make_dm_data()
    expect_error(
        diffMethyl(obj, formula = ~ condition, method = "quasi_f"),
        "limma"
    )
})

# ─── .applyMultipleTesting() direct tests ────────────────────────────────────

test_that("applyMultipleTesting: BH correction returns values in [0, 1]", {
    pvals <- c(0.01, 0.05, 0.1, 0.5, 0.9)
    padj  <- commaKit:::.applyMultipleTesting(pvals, method = "BH")
    expect_true(all(padj >= 0 & padj <= 1))
})

test_that("applyMultipleTesting: method='none' returns original p-values unchanged", {
    pvals <- c(0.01, 0.05, 0.1, 0.5, 0.9)
    padj  <- commaKit:::.applyMultipleTesting(pvals, method = "none")
    expect_equal(padj, pvals)
})

test_that("applyMultipleTesting: NA values pass through as NA", {
    pvals <- c(0.01, NA_real_, 0.1)
    padj  <- commaKit:::.applyMultipleTesting(pvals, method = "BH")
    expect_true(is.na(padj[2]))
})

test_that("applyMultipleTesting: output length equals input length", {
    pvals <- c(0.001, 0.01, 0.05, 0.1)
    padj  <- commaKit:::.applyMultipleTesting(pvals, method = "BH")
    expect_equal(length(padj), length(pvals))
})

test_that("applyMultipleTesting: bonferroni method accepted without error", {
    pvals <- c(0.01, 0.05, 0.1)
    expect_no_error(commaKit:::.applyMultipleTesting(pvals, method = "bonferroni"))
})

# ─── Edge cases ───────────────────────────────────────────────────────────────

test_that("diffMethyl: site with single condition after NA removal gets NA p-value", {
    # Set all treatment sample methylation to NA → only 'control' group present
    # for every site → GLM cannot be fitted → all p-values NA
    obj    <- .make_dm_data(n_sites = 5L, n_ctrl = 2L, n_treat = 1L)
    methyl <- SummarizedExperiment::assay(obj, "methylation")
    methyl[, "treat_1"] <- NA_real_
    SummarizedExperiment::assay(obj, "methylation") <- methyl
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    expect_true(all(is.na(rd$dm_pvalue)))
})

test_that("diffMethyl: perfect separation (ctrl=0, treat=1) does not crash", {
    # Perfect separation may cause GLM non-convergence; result should be
    # NA or a valid p-value — never an error or NaN outside [0,1].
    obj    <- .make_dm_data(n_sites = 5L)
    methyl <- SummarizedExperiment::assay(obj, "methylation")
    methyl[, "ctrl_1"]  <- 0.0
    methyl[, "ctrl_2"]  <- 0.0
    methyl[, "treat_1"] <- 1.0
    SummarizedExperiment::assay(obj, "methylation") <- methyl
    # suppress the expected GLM non-convergence warnings
    dm <- expect_no_error(suppressWarnings(
        diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    ))
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))
    pv <- rd$dm_pvalue
    expect_true(all(is.na(pv) | (pv >= 0 & pv <= 1)))
})

test_that("diffMethyl: site with zero coverage in all samples gets NA p-value", {
    obj    <- .make_dm_data(n_sites = 3L)
    # Zero out coverage for the first site across all samples
    cov    <- SummarizedExperiment::assay(obj, "coverage")
    cov[1L, ] <- 0L
    SummarizedExperiment::assay(obj, "coverage") <- cov
    SummarizedExperiment::assay(obj, "mod_counts")[1L, ] <- 0L
    SummarizedExperiment::assay(obj, "canonical_counts")[1L, ] <- 0L
    methyl <- SummarizedExperiment::assay(obj, "methylation")
    methyl[1L, ] <- NA_real_
    SummarizedExperiment::assay(obj, "methylation") <- methyl
    dm <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))
    expect_true(is.na(rd$dm_pvalue[1]))
})

test_that("diffMethyl methylkit: all-zero-coverage site does not crash calculateDiffMeth", {
    # Regression test: glm.fit inside methylKit's logReg crashes with
    # 'object of type closure is not subsettable' when a site has
    # zero coverage across ALL samples. Verify the wrapper filters
    # such sites before calling calculateDiffMeth.
    skip_if_not(requireNamespace("methylKit", quietly = TRUE), "methylKit not installed")
    obj  <- .make_dm_data(n_sites = 3L, n_ctrl = 3L, n_treat = 3L)
    cov  <- SummarizedExperiment::assay(obj, "coverage")
    cov[2L, ] <- 0L
    SummarizedExperiment::assay(obj, "coverage") <- cov
    SummarizedExperiment::assay(obj, "mod_counts")[2L, ] <- 0L
    SummarizedExperiment::assay(obj, "canonical_counts")[2L, ] <- 0L
    mth  <- SummarizedExperiment::assay(obj, "methylation")
    mth[2L, ] <- NA_real_
    SummarizedExperiment::assay(obj, "methylation") <- mth
    expect_no_error(
        dm <- suppressWarnings(
            diffMethyl(obj, formula = ~ condition, method = "methylkit")
        )
    )
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))
    # The all-zero site (row 2) is untestable and retains p = NA.
    # The critical check is no crash.
    expect_true(is.na(rd$dm_pvalue[2L]))
    # Other sites should have p-values (not NA or 1 from skip_idx)
    expect_false(is.na(rd$dm_pvalue[1L]))
    expect_false(is.na(rd$dm_pvalue[3L]))
})

# ─── Ground-truth recovery on comma_example_data ─────────────────────────────

test_that("diffMethyl: dm_delta_beta is strongly negative for is_diff 6mA sites in comma_example_data", {
    # The 30 is_diff 6mA sites have control ~0.90 and treatment ~0.25 (set.seed(42)),
    # so dm_delta_beta (treatment - control) should be substantially negative.
    # Non-diff sites have both conditions at ~0.90, so delta_beta should be near 0.
    data(comma_example_data)
    dm   <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA",
                       method = "quasi_f")
    rd   <- as.data.frame(SummarizedExperiment::rowData(dm))
    rd6  <- rd[rd$mod_type == "6mA", ]

    delta_diff    <- rd6$dm_delta_beta[rd6$is_diff    & !is.na(rd6$dm_delta_beta)]
    delta_nondiff <- rd6$dm_delta_beta[!rd6$is_diff   & !is.na(rd6$dm_delta_beta)]

    # Simulated differential sites: control ~0.90, treatment ~0.25 → delta < -0.3
    expect_lt(median(delta_diff), -0.3)
    # Non-differential sites: both conditions ~0.90 → |delta| near 0
    expect_lt(abs(median(delta_nondiff)), 0.1)
})

test_that("diffMethyl: majority of is_diff 6mA sites recovered at pvalue < 0.2 in comma_example_data", {
    # With a strong simulated signal (~0.65 delta_beta for 30 sites), diffMethyl()
    # should detect at least half of the ground-truth differentially methylated sites.
    # Note: comma_example_data has 6 samples (3 control + 3 treatment), but this
    # small simulated fixture still makes FDR-corrected padj < 0.05 too strict
    # for a recovery smoke test. Use a lenient raw p-value threshold (0.2) to
    # verify the model correctly ranks differential sites.
    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA",
                      method = "quasi_f")
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    rd6 <- rd[rd$mod_type == "6mA", ]

    n_true_diff <- sum(rd6$is_diff, na.rm = TRUE)      # 30 ground-truth sites
    n_detected  <- sum(rd6$is_diff & !is.na(rd6$dm_pvalue) & rd6$dm_pvalue < 0.2,
                       na.rm = TRUE)

    # Conservative: at least 50% recall (15/30) using raw p-value
    expect_gte(n_detected, floor(n_true_diff * 0.5))
})

test_that("diffMethyl: significant hits are enriched for is_diff sites in comma_example_data", {
    # Among sites called significant (padj < 0.05, |delta_beta| > 0.2), the majority
    # should be ground-truth is_diff = TRUE — i.e., precision >= 50%.
    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA",
                      method = "quasi_f")
    sig <- filterResults(dm, padj = 0.05, delta_beta = 0.2, mod_type = "6mA")

    if (nrow(sig) > 0L) {
        precision <- mean(sig$is_diff, na.rm = TRUE)
        expect_gte(precision, 0.5)
    } else {
        skip("No significant sites found; cannot evaluate precision.")
    }
})

# ─────────────────────────────────────────────────────────────────────────────
# mod_context parameter and per-context looping
# ─────────────────────────────────────────────────────────────────────────────

test_that("diffMethyl: mod_context parameter filters to matching contexts", {
    data(comma_example_data)
    dm_6mA <- diffMethyl(comma_example_data, formula = ~ condition,
                          mod_context = "6mA_GATC", method = "quasi_f")
    si <- siteInfo(dm_6mA)
    # Only 6mA_GATC sites in result
    expect_true(all(si$mod_context[!is.na(si$dm_pvalue)] == "6mA_GATC"))
})

test_that("diffMethyl: loops separately over each mod_context", {
    data(comma_example_data)
    # comma_example_data has 2 contexts: 6mA_GATC and 5mC_CCWGG
    dm <- diffMethyl(comma_example_data, formula = ~ condition, method = "quasi_f")
    si <- siteInfo(dm)
    # Both contexts should have results (non-NA pvalues)
    has_6mA <- any(!is.na(si$dm_pvalue[si$mod_context == "6mA_GATC"]))
    has_5mC <- any(!is.na(si$dm_pvalue[si$mod_context == "5mC_CCWGG"]))
    expect_true(has_6mA)
    expect_true(has_5mC)
})

test_that("diffMethyl: mod_context stored in metadata params", {
    data(comma_example_data)
    dm <- diffMethyl(comma_example_data, formula = ~ condition,
                     mod_context = "6mA_GATC", method = "quasi_f")
    params <- S4Vectors::metadata(dm)$diffMethyl_params
    expect_equal(params$mod_context, "6mA_GATC")
})

# ─── reference argument and factor-level tests ───────────────────────────────

# Helper: build a commaData where WT has lower methylation than HNS
.make_ref_test_data <- function(as_factor = FALSE) {
    set.seed(7L)
    n_sites <- 20L
    # WT (n=2) ~ 0.2 methylation, HNS (n=2) ~ 0.8 methylation
    # delta_beta (HNS - WT) should be positive
    methyl_mat <- cbind(
        matrix(pmax(0, pmin(1, rnorm(n_sites * 2L, 0.2, 0.03))),
               nrow = n_sites, ncol = 2L),
        matrix(pmax(0, pmin(1, rnorm(n_sites * 2L, 0.8, 0.03))),
               nrow = n_sites, ncol = 2L)
    )
    colnames(methyl_mat) <- c("wt_1", "wt_2", "hns_1", "hns_2")
    cov_mat <- matrix(30L, nrow = n_sites, ncol = 4L,
                      dimnames = list(NULL, colnames(methyl_mat)))
    cond_vals <- c("WT", "WT", "HNS", "HNS")
    if (as_factor) {
        cond_vals <- factor(cond_vals, levels = c("WT", "HNS"))
    }
    sample_info <- data.frame(
        sample_name = colnames(methyl_mat),
        condition   = cond_vals,
        replicate   = 1:4,
        stringsAsFactors = FALSE
    )
    .make_commaData_fixture(
        methyl_mat,
        cov_mat,
        sample_info,
        positions = seq_len(n_sites) * 100L
    )
}

test_that("diffMethyl: factor condition uses first factor level as reference", {
    # WT is first factor level, HNS is second; WT has lower methylation.
    # delta_beta = HNS - WT should be positive.
    obj <- .make_ref_test_data(as_factor = TRUE)
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    db  <- rd$dm_delta_beta[!is.na(rd$dm_delta_beta)]
    expect_true(mean(db) > 0,
                info = "Expected delta_beta > 0 (HNS - WT) with factor levels")
    # Verify reference recorded in metadata
    params <- S4Vectors::metadata(dm)$diffMethyl_params
    expect_equal(params$reference, "WT")
})

test_that("diffMethyl: reference argument overrides alphabetical default", {
    # Without reference, alphabetical ordering picks HNS as ref (H < W).
    # With reference = "WT", delta_beta should be positive (HNS - WT).
    obj <- .make_ref_test_data(as_factor = FALSE)
    dm  <- diffMethyl(obj, formula = ~ condition, reference = "WT",
                      method = "quasi_f")
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    db  <- rd$dm_delta_beta[!is.na(rd$dm_delta_beta)]
    expect_true(mean(db) > 0,
                info = "Expected delta_beta > 0 (HNS - WT) with reference = 'WT'")
    params <- S4Vectors::metadata(dm)$diffMethyl_params
    expect_equal(params$reference, "WT")
})

test_that("diffMethyl: invalid reference value produces informative error", {
    obj <- .make_ref_test_data(as_factor = FALSE)
    expect_error(
        diffMethyl(obj, formula = ~ condition, reference = "CTRL"),
        regexp = "'reference' value 'CTRL' not found"
    )
})

# ─── Two-level contrast contract (#135/#137) ─────────────────────────────────

.make_three_level_dm_data <- function() {
    obj <- .make_dm_data(n_sites = 12L, n_ctrl = 2L, n_treat = 2L)
    methyl_mat <- methylation(obj)
    colnames(methyl_mat) <- c("A_1", "A_2", "B_1", "C_1")
    # Make the group directions intentionally different so B-A and C-A would
    # disagree if a backend silently chose the wrong contrast.
    methyl_mat[, "A_1"] <- 0.5
    methyl_mat[, "A_2"] <- 0.5
    methyl_mat[, "B_1"] <- 0.2
    methyl_mat[, "C_1"] <- 0.9
    SummarizedExperiment::assay(obj, "methylation", withDimnames = FALSE) <- methyl_mat
    cov_mat <- SummarizedExperiment::assay(obj, "coverage", withDimnames = FALSE)
    colnames(cov_mat) <- colnames(methyl_mat)
    SummarizedExperiment::assay(obj, "coverage", withDimnames = FALSE) <- cov_mat
    SummarizedExperiment::colData(obj)$sample_name <- colnames(methyl_mat)
    SummarizedExperiment::colData(obj)$condition <- c("A", "A", "B", "C")
    rownames(SummarizedExperiment::colData(obj)) <- colnames(methyl_mat)
    obj
}

test_that("diffMethyl: errors clearly for primary variables with more than 2 levels", {
    obj <- .make_three_level_dm_data()
    expect_error(
        diffMethyl(obj, formula = ~ condition, method = "quasi_f"),
        "currently supports exactly 2 levels"
    )
})

test_that("diffMethyl: two-level factor reference defines treatment-reference direction", {
    obj <- .make_dm_data(n_sites = 10L, n_ctrl = 2L, n_treat = 2L)
    SummarizedExperiment::colData(obj)$condition <- factor(
        SummarizedExperiment::colData(obj)$condition,
        levels = c("treatment", "control")
    )
    dm <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))
    # Factor reference is treatment, so delta is control - treatment; first sites
    # have control high and treatment low.
    expect_true(all(rd$dm_delta_beta[1:5] > 0, na.rm = TRUE))
    expect_equal(S4Vectors::metadata(dm)$diffMethyl_params$reference, "treatment")
    expect_equal(S4Vectors::metadata(dm)$diffMethyl_params$treatment, "control")
})

test_that("diffMethyl: invalid p_adjust_method errors before p.adjust", {
    obj <- .make_dm_data()
    expect_error(
        diffMethyl(obj, method = "quasi_f", p_adjust_method = "bogus"),
        "p_adjust_method"
    )
})
