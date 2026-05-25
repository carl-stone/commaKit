# Tests for plot_methylation_distribution()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_dist_data <- function() {
    n_sites   <- 10L
    positions <- seq(1000L, 10000L, by = 1000L)
    set.seed(1L)
    betas <- matrix(
        runif(n_sites * 3L, 0.1, 0.9),
        nrow = n_sites, ncol = 3L,
        dimnames = list(NULL, c("ctrl_1", "ctrl_2", "treat_1"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 3L,
                      dimnames = dimnames(betas))
    site_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr_sim", n_sites),
        ranges   = IRanges::IRanges(start = positions, width = 1L),
        strand   = rep("+", n_sites),
        mod_type    = factor(rep("6mA", n_sites), levels = c("4mC", "5mC", "6mA")),
        motif       = rep("GATC", n_sites)
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
        assays     = list(methylation = betas, coverage = cov_mat),
        rowRanges  = site_gr,
        colData    = cd
    )
    new("commaData", rse)
}

## Object with two modification types
.make_dist_data_two_mods <- function() {
    n_6ma <- 8L; n_5mc <- 4L
    n_sites <- n_6ma + n_5mc
    positions <- seq(1000L, n_sites * 1000L, by = 1000L)
    mod_types  <- factor(c(rep("6mA", n_6ma), rep("5mC", n_5mc)),
                         levels = c("4mC", "5mC", "6mA"))
    motif_vals <- c(rep("GATC", n_6ma), rep("CCWGG", n_5mc))
    set.seed(2L)
    betas <- matrix(
        runif(n_sites * 2L, 0.1, 0.9),
        nrow = n_sites, ncol = 2L,
        dimnames = list(NULL, c("samp1", "samp2"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 2L,
                      dimnames = dimnames(betas))
    site_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr_sim", n_sites),
        ranges   = IRanges::IRanges(start = positions, width = 1L),
        strand   = rep("+", n_sites),
        mod_type    = mod_types,
        motif       = motif_vals
    )
    GenomeInfoDb::seqinfo(site_gr) <- GenomeInfoDb::Seqinfo(
        seqnames = "chr_sim",
        seqlengths = 100000L,
        isCircular = FALSE
    )
    cd <- S4Vectors::DataFrame(
        sample_name = c("samp1", "samp2"),
        condition   = c("ctrl", "treat"),
        replicate   = 1:2,
        row.names   = c("samp1", "samp2")
    )
    rse <- SummarizedExperiment::SummarizedExperiment(
        assays     = list(methylation = betas, coverage = cov_mat),
        rowRanges  = site_gr,
        colData    = cd
    )
    new("commaData", rse)
}

# ─── Basic return type ────────────────────────────────────────────────────────

test_that("plot_methylation_distribution: returns ggplot for valid input", {
    obj <- .make_dist_data()
    p <- plot_methylation_distribution(obj)
    expect_s3_class(p, "ggplot")
    d <- p$data
    # Required columns present
    expect_true("beta" %in% colnames(d))
    expect_true("sample_name" %in% colnames(d))
    # Row count: n_sites * n_samples (no NAs in fixture)
    expect_equal(nrow(d), 10L * 3L)
    # Beta values match the methylation assay
    methyl_mat <- methylation(obj)
    expect_equal(sort(d$beta), sort(as.vector(methyl_mat)))
})

test_that("plot_methylation_distribution: mod_type filter returns ggplot", {
    obj <- .make_dist_data_two_mods()
    p_unfiltered <- plot_methylation_distribution(obj)
    p <- plot_methylation_distribution(obj, mod_type = "6mA")
    expect_s3_class(p, "ggplot")
    # Filtered data has fewer rows than unfiltered (filters to one mod_type)
    expect_lt(nrow(p$data), nrow(p_unfiltered$data))
})

# ─── Faceting ─────────────────────────────────────────────────────────────────

test_that("plot_methylation_distribution: multi-mod object produces facets", {
    obj <- .make_dist_data_two_mods()
    p <- plot_methylation_distribution(obj)
    expect_s3_class(p, "ggplot")
    # facet_wrap wraps in a FacetWrap layer, not FacetNull
    expect_false(inherits(p$facet, "FacetNull"))
})

test_that("plot_methylation_distribution: single-mod object has no facets", {
    obj <- .make_dist_data()
    p <- plot_methylation_distribution(obj)
    expect_true(inherits(p$facet, "FacetNull"))
})

# ─── NA handling ─────────────────────────────────────────────────────────────

test_that("plot_methylation_distribution: NAs in beta values are silently excluded", {
    obj <- .make_dist_data()
    # Inject NAs into the methylation matrix
    methyl_mat <- methylation(obj)
    methyl_mat[1:3, "ctrl_1"] <- NA
    SummarizedExperiment::assay(obj, "methylation") <- methyl_mat
    p <- plot_methylation_distribution(obj)
    expect_s3_class(p, "ggplot")
    # Verify no NA betas remain in p$data
    expect_true(all(!is.na(p$data$beta)))
})

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_methylation_distribution: error on non-commaData input", {
    expect_error(plot_methylation_distribution(data.frame(x = 1)),
                 "commaData")
})

test_that("plot_methylation_distribution: error on invalid mod_type", {
    obj <- .make_dist_data()
    expect_error(plot_methylation_distribution(obj, mod_type = "4mC"),
                 "not found")
})

test_that("plot_methylation_distribution: error on invalid per_sample", {
    obj <- .make_dist_data()
    expect_error(plot_methylation_distribution(obj, per_sample = "yes"),
                 "per_sample")
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_methylation_distribution: works with comma_example_data", {
    data(comma_example_data)
    p <- plot_methylation_distribution(comma_example_data)
    expect_s3_class(p, "ggplot")
    expect_gt(nrow(p$data), 0L)
})
