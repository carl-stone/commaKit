# Tests for plot_pca()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_pca_data <- function() {
    n_sites <- 20L
    positions <- seq(1000L, 20000L, by = 1000L)
    set.seed(42L)
    betas <- matrix(
        runif(n_sites * 3L, 0.0, 1.0),
        nrow = n_sites, ncol = 3L,
        dimnames = list(NULL, c("ctrl_1", "ctrl_2", "treat_1"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 3L,
                      dimnames = dimnames(betas))
    sample_info <- data.frame(
        sample_name = c("ctrl_1", "ctrl_2", "treat_1"),
        condition   = c("control", "control", "treatment"),
        replicate   = 1:3,
        stringsAsFactors = FALSE
    )
    .make_commaData_fixture(betas, cov_mat, sample_info, positions)
}

# ─── Basic return type ────────────────────────────────────────────────────────

test_that("plot_pca: returns ggplot for valid input", {
    obj <- .make_pca_data()
    p <- plot_pca(obj)
    expect_s3_class(p, "ggplot")
    # p$data has exactly one row per sample and required columns
    expect_equal(nrow(p$data), ncol(obj))
    expect_true(all(c("PC1", "PC2", "sample_name") %in% colnames(p$data)))
})

test_that("plot_pca: color_by argument accepted", {
    obj <- .make_pca_data()
    p <- plot_pca(obj, color_by = "condition")
    expect_s3_class(p, "ggplot")
    # Verify colour mapping references the color_by column
    expect_true("colour" %in% names(p$mapping))
    expect_equal(p$labels$colour, "condition")
})

test_that("plot_pca: shape_by = NULL accepted without error", {
    obj <- .make_pca_data()
    p <- plot_pca(obj, shape_by = NULL)
    expect_s3_class(p, "ggplot")
    # Verify no shape aesthetic mapped
    expect_false("shape" %in% names(p$mapping))
})

test_that("plot_pca: mod_type filter reduces sites used", {
    data(comma_example_data)
    p_unfiltered <- plot_pca(comma_example_data)
    p_filtered    <- plot_pca(comma_example_data, mod_type = "6mA")
    expect_s3_class(p_filtered, "ggplot")
    # Same samples, different sites used → different PC1 values
    d_unfiltered <- p_unfiltered$data[order(p_unfiltered$data$sample_name), ]
    d_filtered   <- p_filtered$data[order(p_filtered$data$sample_name), ]
    expect_false(identical(d_filtered$PC1, d_unfiltered$PC1))
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

test_that("plot_pca: default color falls back to sample_name without condition", {
    obj <- .make_pca_data()
    SummarizedExperiment::colData(obj)$condition <- NULL

    p <- plot_pca(obj)
    d <- plot_pca(obj, return_data = TRUE)

    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$colour, "sample_name")
    expect_true("sample_name" %in% colnames(d))
    expect_false("condition" %in% colnames(d))
})

test_that("plot_pca: explicit condition color still validates missing column", {
    obj <- .make_pca_data()
    SummarizedExperiment::colData(obj)$condition <- NULL

    expect_error(
        plot_pca(obj, color_by = "condition"),
        "'color_by' column 'condition' not found"
    )
})

test_that("plot_pca: return_data attaches percentVar attribute", {
    obj <- .make_pca_data()
    d <- plot_pca(obj, return_data = TRUE)
    pv <- attr(d, "percentVar")
    expect_true(is.numeric(pv))
    expect_true(length(pv) <= 2L)
    expect_true(all(pv >= 0 & pv <= 100))
})

test_that("plot_pca: return_data includes requested color and shape columns", {
    obj <- .make_pca_data()
    d <- plot_pca(obj,
                  color_by = "condition",
                  shape_by = "replicate",
                  return_data = TRUE)
    si <- as.data.frame(sampleInfo(obj))
    si <- si[match(d$sample_name, si$sample_name), , drop = FALSE]

    expect_true(all(c("condition", "replicate") %in% colnames(d)))
    expect_equal(d$condition, si$condition)
    expect_equal(d$replicate, si$replicate)
})

test_that("plot_pca: return_data = TRUE validates missing color_by column", {
    obj <- .make_pca_data()
    expect_error(
        plot_pca(obj, color_by = "nonexistent_col", return_data = TRUE),
        "'color_by' column 'nonexistent_col' not found"
    )
})

test_that("plot_pca: return_data = TRUE validates missing shape_by column", {
    obj <- .make_pca_data()
    expect_error(
        plot_pca(obj, shape_by = "nonexistent_col", return_data = TRUE),
        "'shape_by' column 'nonexistent_col' not found"
    )
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_pca: works with comma_example_data", {
    data(comma_example_data)
    p <- plot_pca(comma_example_data)
    expect_s3_class(p, "ggplot")
    # 6 samples → 6 rows in plot data
    expect_equal(nrow(p$data), 6L)
})
