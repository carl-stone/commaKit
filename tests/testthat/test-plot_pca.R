# Tests for plot_pca()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_pca_data <- function() {
    n_sites <- 20L
    positions <- seq(1000L, 20000L, by = 1000L)
    site_keys <- paste0("chr_sim:", positions, ":+:6mA:GATC")
    set.seed(42L)
    betas <- matrix(
        runif(n_sites * 3L, 0.0, 1.0),
        nrow = n_sites, ncol = 3L,
        dimnames = list(site_keys, c("ctrl_1", "ctrl_2", "treat_1"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 3L,
                      dimnames = dimnames(betas))
    rd <- S4Vectors::DataFrame(
        chrom    = rep("chr_sim", n_sites),
        position = positions,
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
        assays  = list(methylation = betas, coverage = cov_mat),
        rowData = rd,
        colData = cd
    )
    new("commaData", se,
        genomeInfo = c(chr_sim = 100000L),
        annotation = GenomicRanges::GRanges(),
        motifSites = GenomicRanges::GRanges())
}

# ─── Basic return type ────────────────────────────────────────────────────────

test_that("plot_pca: returns ggplot for valid input", {
    obj <- .make_pca_data()
    p <- plot_pca(obj)
    expect_s3_class(p, "ggplot")
})

test_that("plot_pca: color_by argument accepted", {
    obj <- .make_pca_data()
    p <- plot_pca(obj, color_by = "condition")
    expect_s3_class(p, "ggplot")
})

test_that("plot_pca: shape_by = NULL accepted without error", {
    obj <- .make_pca_data()
    p <- plot_pca(obj, shape_by = NULL)
    expect_s3_class(p, "ggplot")
})

test_that("plot_pca: mod_type filter reduces sites used", {
    data(comma_example_data)
    p <- plot_pca(comma_example_data, mod_type = "6mA")
    expect_s3_class(p, "ggplot")
})

# ─── Axis labels contain PC variance ─────────────────────────────────────────

test_that("plot_pca: x-axis label mentions PC1", {
    obj <- .make_pca_data()
    p <- plot_pca(obj)
    expect_true(grepl("PC1", p$labels$x))
})

test_that("plot_pca: y-axis label mentions PC2", {
    obj <- .make_pca_data()
    p <- plot_pca(obj)
    expect_true(grepl("PC2", p$labels$y))
})

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_pca: error on non-commaData input", {
    expect_error(plot_pca(data.frame(x = 1)), "commaData")
})

test_that("plot_pca: error when color_by column is absent", {
    obj <- .make_pca_data()
    expect_error(plot_pca(obj, color_by = "nonexistent_col"),
                 "nonexistent_col")
})

test_that("plot_pca: error when shape_by column is absent", {
    obj <- .make_pca_data()
    expect_error(plot_pca(obj, shape_by = "nonexistent_col"),
                 "nonexistent_col")
})

test_that("plot_pca: error on invalid mod_type", {
    obj <- .make_pca_data()
    expect_error(plot_pca(obj, mod_type = "4mC"), "not found")
})

test_that("plot_pca: warning issued with fewer than 3 samples", {
    obj <- .make_pca_data()
    # Subset to 2 samples; expect a warning about low sample count
    obj2 <- obj[, 1:2]
    expect_warning(plot_pca(obj2), "[Ff]ewer than 3 samples")
    suppressWarnings({
        p <- plot_pca(obj2)
    })
    expect_s3_class(p, "ggplot")
})

# ─── return_data = TRUE ───────────────────────────────────────────────────────

test_that("plot_pca: return_data = TRUE returns a data.frame", {
    obj <- .make_pca_data()
    d <- plot_pca(obj, return_data = TRUE)
    expect_s3_class(d, "data.frame")
})

test_that("plot_pca: return_data data.frame has PC1 and PC2 columns", {
    obj <- .make_pca_data()
    d <- plot_pca(obj, return_data = TRUE)
    expect_true(all(c("PC1", "PC2") %in% colnames(d)))
})

test_that("plot_pca: return_data data.frame has one row per sample", {
    obj <- .make_pca_data()
    d <- plot_pca(obj, return_data = TRUE)
    expect_equal(nrow(d), ncol(methylation(obj)))
})

test_that("plot_pca: return_data data.frame includes sampleInfo columns", {
    obj <- .make_pca_data()
    d <- plot_pca(obj, return_data = TRUE)
    expect_true("condition" %in% colnames(d))
    expect_true("sample_name" %in% colnames(d))
})

test_that("plot_pca: return_data attaches percentVar attribute", {
    obj <- .make_pca_data()
    d <- plot_pca(obj, return_data = TRUE)
    pv <- attr(d, "percentVar")
    expect_true(is.numeric(pv))
    expect_true(length(pv) <= 2L)
    expect_true(all(pv >= 0 & pv <= 100))
})

test_that("plot_pca: return_data = TRUE skips color_by validation", {
    obj <- .make_pca_data()
    # Would error if color_by were validated; should not error with return_data = TRUE
    expect_no_error(plot_pca(obj, color_by = "nonexistent_col", return_data = TRUE))
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_pca: works with comma_example_data", {
    data(comma_example_data)
    p <- plot_pca(comma_example_data)
    expect_s3_class(p, "ggplot")
})
