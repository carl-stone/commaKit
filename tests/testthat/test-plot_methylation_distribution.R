# Tests for plot_methylation_distribution()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_dist_data <- function() {
    n_sites   <- 10L
    positions <- seq(1000L, 10000L, by = 1000L)
    site_keys <- paste0("chr_sim:", positions, ":+:6mA:GATC")
    set.seed(1L)
    betas <- matrix(
        runif(n_sites * 3L, 0.1, 0.9),
        nrow = n_sites, ncol = 3L,
        dimnames = list(site_keys, c("ctrl_1", "ctrl_2", "treat_1"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 3L,
                      dimnames = dimnames(betas))
    rd <- S4Vectors::DataFrame(
        chrom       = rep("chr_sim", n_sites),
        position    = positions,
        strand      = rep("+", n_sites),
        mod_type    = rep("6mA", n_sites),
        motif       = rep("GATC", n_sites),
        mod_context = rep("6mA_GATC", n_sites),
        row.names   = site_keys
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

## Object with two modification types
.make_dist_data_two_mods <- function() {
    n_6ma <- 8L; n_5mc <- 4L
    n_sites <- n_6ma + n_5mc
    positions <- seq(1000L, n_sites * 1000L, by = 1000L)
    mod_types  <- c(rep("6mA", n_6ma), rep("5mC", n_5mc))
    motif_vals <- c(rep("GATC", n_6ma), rep("CCWGG", n_5mc))
    site_keys  <- paste0("chr_sim:", positions, ":+:", mod_types, ":", motif_vals)
    set.seed(2L)
    betas <- matrix(
        runif(n_sites * 2L, 0.1, 0.9),
        nrow = n_sites, ncol = 2L,
        dimnames = list(site_keys, c("samp1", "samp2"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 2L,
                      dimnames = dimnames(betas))
    rd <- S4Vectors::DataFrame(
        chrom       = rep("chr_sim", n_sites),
        position    = positions,
        strand      = rep("+", n_sites),
        mod_type    = mod_types,
        motif       = motif_vals,
        mod_context = paste(mod_types, motif_vals, sep = "_"),
        row.names   = site_keys
    )
    cd <- S4Vectors::DataFrame(
        sample_name = c("samp1", "samp2"),
        condition   = c("ctrl", "treat"),
        replicate   = 1:2,
        row.names   = c("samp1", "samp2")
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

test_that("plot_methylation_distribution: returns ggplot for valid input", {
    obj <- .make_dist_data()
    p <- plot_methylation_distribution(obj)
    expect_s3_class(p, "ggplot")
})

test_that("plot_methylation_distribution: per_sample = FALSE returns ggplot", {
    obj <- .make_dist_data()
    p <- plot_methylation_distribution(obj, per_sample = FALSE)
    expect_s3_class(p, "ggplot")
})

test_that("plot_methylation_distribution: mod_type filter returns ggplot", {
    obj <- .make_dist_data_two_mods()
    p <- plot_methylation_distribution(obj, mod_type = "6mA")
    expect_s3_class(p, "ggplot")
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
})

test_that("plot_methylation_distribution: example data, mod_type = '6mA'", {
    data(comma_example_data)
    p <- plot_methylation_distribution(comma_example_data, mod_type = "6mA")
    expect_s3_class(p, "ggplot")
})
