# Tests for mValues()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_mval_data <- function(n_sites = 10L, n_samples = 3L,
                             seed = 99L, cov_val = 20L) {
    set.seed(seed)
    positions <- seq(1000L, by = 1000L, length.out = n_sites)
    site_keys <- paste0("chr_sim:", positions, ":+:6mA")
    betas <- matrix(
        runif(n_sites * n_samples, 0.05, 0.95),
        nrow = n_sites, ncol = n_samples,
        dimnames = list(site_keys, paste0("s", seq_len(n_samples)))
    )
    cov_mat <- matrix(cov_val, nrow = n_sites, ncol = n_samples,
                      dimnames = dimnames(betas))
    rd <- S4Vectors::DataFrame(
        chrom    = rep("chr_sim", n_sites),
        position = positions,
        strand   = rep("+", n_sites),
        mod_type = rep("6mA", n_sites),
        row.names = site_keys
    )
    cd <- S4Vectors::DataFrame(
        sample_name = paste0("s", seq_len(n_samples)),
        condition   = rep("ctrl", n_samples),
        replicate   = seq_len(n_samples),
        row.names   = paste0("s", seq_len(n_samples))
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = betas, coverage = cov_mat),
        rowData = rd,
        colData = cd
    )
    new("commaData", se,
        genomeInfo = c(chr_sim = 100000L),
        annotation = GenomicRanges::GRanges(),
        motifSites = GenomicRanges::GRanges())
}

# ─── Return type and dimensions ───────────────────────────────────────────────

test_that("mValues: returns a numeric matrix", {
    obj <- .make_mval_data()
    m <- mValues(obj)
    expect_true(is.matrix(m))
    expect_true(is.numeric(m))
})

test_that("mValues: dimensions match methylation matrix", {
    obj <- .make_mval_data()
    m <- mValues(obj)
    expect_equal(dim(m), dim(methylation(obj)))
})

test_that("mValues: dimnames match methylation matrix", {
    obj <- .make_mval_data()
    m <- mValues(obj)
    expect_equal(dimnames(m), dimnames(methylation(obj)))
})

# ─── Numeric correctness ──────────────────────────────────────────────────────

test_that("mValues: formula is correct for known values", {
    # beta = 0.8, coverage = 10, alpha = 0.5
    # m_reads = round(0.8 * 10) = 8, u_reads = 2
    # M = log2((8 + 0.5) / (2 + 0.5)) = log2(8.5 / 2.5)
    obj <- .make_mval_data(n_sites = 1L, n_samples = 1L, cov_val = 10L)
    # Override the beta to exactly 0.8
    SummarizedExperiment::assay(obj, "methylation")[1, 1] <- 0.8
    m <- mValues(obj, alpha = 0.5)
    expected <- log2((8 + 0.5) / (2 + 0.5))
    expect_equal(m[1, 1], expected)
})

test_that("mValues: M-value at beta=0.5 is near zero", {
    obj <- .make_mval_data(n_sites = 1L, n_samples = 1L, cov_val = 100L)
    SummarizedExperiment::assay(obj, "methylation")[1, 1] <- 0.5
    m <- mValues(obj, alpha = 0.5)
    expect_true(abs(m[1, 1]) < 0.1)
})

test_that("mValues: high beta gives positive M-value", {
    obj <- .make_mval_data(n_sites = 1L, n_samples = 1L, cov_val = 20L)
    SummarizedExperiment::assay(obj, "methylation")[1, 1] <- 0.9
    m <- mValues(obj, alpha = 0.5)
    expect_true(m[1, 1] > 0)
})

test_that("mValues: low beta gives negative M-value", {
    obj <- .make_mval_data(n_sites = 1L, n_samples = 1L, cov_val = 20L)
    SummarizedExperiment::assay(obj, "methylation")[1, 1] <- 0.1
    m <- mValues(obj, alpha = 0.5)
    expect_true(m[1, 1] < 0)
})

test_that("mValues: all finite for typical beta and coverage values", {
    obj <- .make_mval_data()
    m <- mValues(obj)
    expect_true(all(is.finite(m[!is.na(m)])))
})

# ─── NA propagation ───────────────────────────────────────────────────────────

test_that("mValues: NA beta propagates to NA M-value", {
    obj <- .make_mval_data()
    SummarizedExperiment::assay(obj, "methylation")[1, 1] <- NA_real_
    m <- mValues(obj)
    expect_true(is.na(m[1, 1]))
})

test_that("mValues: zero coverage produces NA M-value", {
    obj <- .make_mval_data()
    SummarizedExperiment::assay(obj, "coverage")[2, 2] <- 0L
    m <- mValues(obj)
    expect_true(is.na(m[2, 2]))
    # Other sites in same row/column are not affected
    expect_false(is.na(m[2, 1]))
})

# ─── alpha parameter ──────────────────────────────────────────────────────────

test_that("mValues: alpha = 0 errors", {
    obj <- .make_mval_data()
    expect_error(mValues(obj, alpha = 0), "positive")
})

test_that("mValues: negative alpha errors", {
    obj <- .make_mval_data()
    expect_error(mValues(obj, alpha = -1), "positive")
})

test_that("mValues: non-numeric alpha errors", {
    obj <- .make_mval_data()
    expect_error(mValues(obj, alpha = "0.5"), "positive")
})

test_that("mValues: alpha = 1 runs without error", {
    obj <- .make_mval_data()
    expect_no_error(mValues(obj, alpha = 1))
})

test_that("mValues: smaller alpha produces more extreme M-values", {
    obj <- .make_mval_data(n_sites = 1L, n_samples = 1L, cov_val = 20L)
    SummarizedExperiment::assay(obj, "methylation")[1, 1] <- 0.9
    m_large_alpha <- mValues(obj, alpha = 2)
    m_small_alpha <- mValues(obj, alpha = 0.1)
    expect_true(m_small_alpha[1, 1] > m_large_alpha[1, 1])
})

# ─── mod_type filter ──────────────────────────────────────────────────────────

test_that("mValues: mod_type filter reduces rows to matching sites", {
    data(comma_example_data)
    m_all  <- mValues(comma_example_data)
    m_6mA  <- mValues(comma_example_data, mod_type = "6mA")
    n_6mA  <- sum(SummarizedExperiment::rowData(comma_example_data)$mod_type == "6mA")
    expect_equal(nrow(m_6mA), n_6mA)
    expect_true(nrow(m_6mA) < nrow(m_all))
})

test_that("mValues: invalid mod_type errors with informative message", {
    obj <- .make_mval_data()
    expect_error(mValues(obj, mod_type = "4mC"), "not found")
})

# ─── Input validation ─────────────────────────────────────────────────────────

test_that("mValues: error on non-commaData input", {
    expect_error(mValues(data.frame(x = 1)), "commaData")
    expect_error(mValues(matrix(0.5, 3, 3)), "commaData")
})

# ─── Integration with comma_example_data ──────────────────────────────────────

test_that("mValues: works with comma_example_data", {
    data(comma_example_data)
    m <- mValues(comma_example_data)
    expect_equal(dim(m), dim(methylation(comma_example_data)))
    expect_true(is.numeric(m))
})
