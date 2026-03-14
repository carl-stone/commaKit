# Tests for plot_heatmap()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_heatmap_fixtures <- function() {
    n_sites   <- 15L
    positions <- seq(1000L, 15000L, by = 1000L)
    site_keys <- paste0("chr_sim:", positions, ":+:6mA")
    set.seed(10L)
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
    obj <- new("commaData", se,
               genomeInfo = c(chr_sim = 100000L),
               annotation = GenomicRanges::GRanges(),
               motifSites = GenomicRanges::GRanges())

    ## Produce results table: use synthetic dm_padj / dm_delta_beta
    set.seed(11L)
    res <- data.frame(
        chrom         = rep("chr_sim", n_sites),
        position      = positions,
        strand        = rep("+", n_sites),
        mod_type      = rep("6mA", n_sites),
        dm_pvalue     = runif(n_sites, 0, 0.1),
        dm_padj       = runif(n_sites, 0, 0.05),
        dm_delta_beta = runif(n_sites, -0.5, 0.5),
        row.names     = site_keys,
        stringsAsFactors = FALSE
    )
    list(obj = obj, res = res)
}

# ─── Basic return type ────────────────────────────────────────────────────────

test_that("plot_heatmap: returns ggplot for valid input", {
    fix <- .make_heatmap_fixtures()
    p <- plot_heatmap(fix$res, fix$obj)
    expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap: n_sites controls number of displayed sites", {
    fix <- .make_heatmap_fixtures()
    p5 <- plot_heatmap(fix$res, fix$obj, n_sites = 5L)
    # The y scale should have at most 5 levels
    bd <- ggplot2::ggplot_build(p5)
    expect_lte(length(unique(bd$data[[1L]]$y)), 5L)
})

test_that("plot_heatmap: n_sites larger than available sites clamps silently", {
    fix <- .make_heatmap_fixtures()
    # n_sites = 1000 > 15 sites available
    p <- plot_heatmap(fix$res, fix$obj, n_sites = 1000L)
    expect_s3_class(p, "ggplot")
})

# ─── NA handling ─────────────────────────────────────────────────────────────

test_that("plot_heatmap: NA beta values handled without error", {
    fix <- .make_heatmap_fixtures()
    methyl_mat <- methylation(fix$obj)
    methyl_mat[1L, "ctrl_1"] <- NA
    SummarizedExperiment::assay(fix$obj, "methylation") <- methyl_mat
    p <- plot_heatmap(fix$res, fix$obj)
    expect_s3_class(p, "ggplot")
})

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_heatmap: error on non-data.frame results", {
    fix <- .make_heatmap_fixtures()
    expect_error(plot_heatmap(list(), fix$obj))
})

test_that("plot_heatmap: error on non-commaData object", {
    fix <- .make_heatmap_fixtures()
    expect_error(plot_heatmap(fix$res, data.frame()), "commaData")
})

test_that("plot_heatmap: error when required results columns are absent", {
    fix <- .make_heatmap_fixtures()
    bad_res <- fix$res
    bad_res$dm_padj <- NULL
    expect_error(plot_heatmap(bad_res, fix$obj), "dm_padj")
})

test_that("plot_heatmap: error when no rows have non-NA padj", {
    fix <- .make_heatmap_fixtures()
    fix$res$dm_padj <- NA_real_
    expect_error(plot_heatmap(fix$res, fix$obj))
})

test_that("plot_heatmap: error on invalid n_sites", {
    fix <- .make_heatmap_fixtures()
    expect_error(plot_heatmap(fix$res, fix$obj, n_sites = 0L), "n_sites")
    expect_error(plot_heatmap(fix$res, fix$obj, n_sites = -1L), "n_sites")
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_heatmap: works with comma_example_data results", {
    data(comma_example_data)
    cd_dm <- diffMethyl(comma_example_data, ~ condition, mod_type = "6mA")
    res <- results(cd_dm)
    p <- plot_heatmap(res, cd_dm, n_sites = 20L)
    expect_s3_class(p, "ggplot")
})
