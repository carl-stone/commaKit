# Tests for diffMethyl() — Phase 4 differential methylation analysis

# ─── Helpers ─────────────────────────────────────────────────────────────────

.make_dm_data <- function(n_sites = 20L, n_ctrl = 2L, n_treat = 1L) {
    set.seed(99L)
    n_samp  <- n_ctrl + n_treat
    site_keys <- paste0("chr_sim:", seq_len(n_sites) * 100L, ":+:6mA:GATC")

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
        motif    = rep("GATC", n_sites),
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

test_that("diffMethyl: limma and beta_binomial delta_beta values are highly correlated", {
    skip_if_not_installed("limma")
    obj  <- .make_dm_data(n_sites = 40L)
    dm_l <- diffMethyl(obj, formula = ~ condition, method = "limma")
    dm_b <- diffMethyl(obj, formula = ~ condition, method = "beta_binomial")
    db_l <- SummarizedExperiment::rowData(dm_l)$dm_delta_beta
    db_b <- SummarizedExperiment::rowData(dm_b)$dm_delta_beta
    ok   <- !is.na(db_l) & !is.na(db_b)
    expect_true(sum(ok) > 0)
    expect_gt(cor(db_l[ok], db_b[ok]), 0.95)
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

test_that("diffMethyl: method='quasi_f' returns commaData with correct columns", {
    skip_if_not_installed("limma")
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    expect_s4_class(dm, "commaData")
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))
    expect_true(all(c("dm_pvalue", "dm_padj", "dm_delta_beta",
                      "dm_mean_beta_control", "dm_mean_beta_treatment") %in%
                        colnames(rd)))
    expect_equal(nrow(rd), nrow(SummarizedExperiment::rowData(obj)))
})

test_that("diffMethyl: method='quasi_f' produces valid p-values in [0, 1]", {
    skip_if_not_installed("limma")
    obj <- .make_dm_data()
    dm  <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    rd  <- as.data.frame(SummarizedExperiment::rowData(dm))
    pvals <- rd$dm_pvalue[!is.na(rd$dm_pvalue)]
    expect_true(length(pvals) > 0)
    expect_true(all(pvals >= 0 & pvals <= 1))
    expect_true(all(rd$dm_padj[!is.na(rd$dm_padj)] >=
                        rd$dm_pvalue[!is.na(rd$dm_pvalue)]))
})

test_that("diffMethyl: quasi_f and beta_binomial delta_beta are highly correlated", {
    skip_if_not_installed("limma")
    obj  <- .make_dm_data(n_sites = 40L)
    dm_q <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    dm_b <- diffMethyl(obj, formula = ~ condition, method = "beta_binomial")
    db_q <- SummarizedExperiment::rowData(dm_q)$dm_delta_beta
    db_b <- SummarizedExperiment::rowData(dm_b)$dm_delta_beta
    ok   <- !is.na(db_q) & !is.na(db_b)
    expect_true(sum(ok) > 0)
    expect_gt(cor(db_q[ok], db_b[ok]), 0.99)
})

test_that("diffMethyl: quasi_f has more power than beta_binomial on ground-truth data", {
    skip_if_not_installed("limma")
    data(comma_example_data)
    dm_q  <- diffMethyl(comma_example_data, formula = ~ condition,
                        method = "quasi_f", mod_type = "6mA")
    dm_b  <- diffMethyl(comma_example_data, formula = ~ condition,
                        method = "beta_binomial", mod_type = "6mA")
    is_diff <- as.data.frame(SummarizedExperiment::rowData(comma_example_data))$is_diff
    pq <- SummarizedExperiment::rowData(dm_q)$dm_pvalue[is_diff]
    pb <- SummarizedExperiment::rowData(dm_b)$dm_pvalue[is_diff]
    # quasi_f should produce lower mean p-value for true positives
    expect_lt(mean(pq, na.rm = TRUE), mean(pb, na.rm = TRUE))
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
    padj  <- comma:::.applyMultipleTesting(pvals, method = "BH")
    expect_true(all(padj >= 0 & padj <= 1))
})

test_that("applyMultipleTesting: method='none' returns original p-values unchanged", {
    pvals <- c(0.01, 0.05, 0.1, 0.5, 0.9)
    padj  <- comma:::.applyMultipleTesting(pvals, method = "none")
    expect_equal(padj, pvals)
})

test_that("applyMultipleTesting: NA values pass through as NA", {
    pvals <- c(0.01, NA_real_, 0.1)
    padj  <- comma:::.applyMultipleTesting(pvals, method = "BH")
    expect_true(is.na(padj[2]))
})

test_that("applyMultipleTesting: output length equals input length", {
    pvals <- c(0.001, 0.01, 0.05, 0.1)
    padj  <- comma:::.applyMultipleTesting(pvals, method = "BH")
    expect_equal(length(padj), length(pvals))
})

test_that("applyMultipleTesting: bonferroni method accepted without error", {
    pvals <- c(0.01, 0.05, 0.1)
    expect_no_error(comma:::.applyMultipleTesting(pvals, method = "bonferroni"))
})

# ─── .betaBinomialTest() edge cases (via diffMethyl) ─────────────────────────

test_that("diffMethyl: site with single condition after NA removal gets NA p-value", {
    # Set all treatment sample methylation to NA → only 'control' group present
    # for every site → GLM cannot be fitted → all p-values NA
    obj    <- .make_dm_data(n_sites = 5L, n_ctrl = 2L, n_treat = 1L)
    methyl <- SummarizedExperiment::assay(obj, "methylation")
    methyl[, "treat_1"] <- NA_real_
    SummarizedExperiment::assay(obj, "methylation") <- methyl
    dm  <- diffMethyl(obj, formula = ~ condition)
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
    dm <- expect_no_error(suppressWarnings(diffMethyl(obj, formula = ~ condition)))
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
    methyl <- SummarizedExperiment::assay(obj, "methylation")
    methyl[1L, ] <- NA_real_
    SummarizedExperiment::assay(obj, "methylation") <- methyl
    dm <- diffMethyl(obj, formula = ~ condition)
    rd <- as.data.frame(SummarizedExperiment::rowData(dm))
    expect_true(is.na(rd$dm_pvalue[1]))
})

# ─── Ground-truth recovery on comma_example_data ─────────────────────────────

test_that("diffMethyl: dm_delta_beta is strongly negative for is_diff 6mA sites in comma_example_data", {
    # The 30 is_diff 6mA sites have control ~0.90 and treatment ~0.25 (set.seed(42)),
    # so dm_delta_beta (treatment - control) should be substantially negative.
    # Non-diff sites have both conditions at ~0.90, so delta_beta should be near 0.
    data(comma_example_data)
    dm   <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
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
    # Note: comma_example_data has only 3 samples (2 control + 1 treatment), giving
    # 1 residual df for the quasibinomial GLM. FDR-corrected padj < 0.05 is not
    # achievable with this few replicates. We use the raw p-value at a lenient
    # threshold (0.2) to verify the model correctly ranks differential sites.
    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
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
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    sig <- filterResults(dm, padj = 0.05, delta_beta = 0.2, mod_type = "6mA")

    if (nrow(sig) > 0L) {
        precision <- mean(sig$is_diff, na.rm = TRUE)
        expect_gte(precision, 0.5)
    } else {
        skip("No significant sites found; cannot evaluate precision.")
    }
})
