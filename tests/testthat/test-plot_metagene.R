# Tests for plot_metagene()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_metagene_data <- function() {
    # 10 sites distributed within annotated genes so metagene has overlap
    n_sites   <- 10L
    positions <- c(1500L, 2000L, 2500L, 3000L, 3500L,
                   6500L, 7000L, 7500L, 8000L, 8500L)
    set.seed(5L)
    betas <- matrix(
        runif(n_sites * 2L, 0.1, 0.9),
        nrow = n_sites, ncol = 2L,
        dimnames = list(NULL, c("ctrl_1", "treat_1"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 2L,
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
        sample_name = c("ctrl_1", "treat_1"),
        condition   = c("control", "treatment"),
        replicate   = 1:2,
        row.names   = c("ctrl_1", "treat_1")
    )
    ann_gr <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = c(1000L, 6000L),
                                     end   = c(4000L, 9000L)),
        strand   = "+"
    )
    ann_gr$feature_type <- c("gene", "gene")
    ann_gr$name         <- c("geneA", "geneB")

    rse <- SummarizedExperiment::SummarizedExperiment(
        assays     = list(methylation = betas, coverage = cov_mat),
        rowRanges  = site_gr,
        colData    = cd
    )
    obj <- new("commaData", rse)
    S4Vectors::metadata(obj)$annotation <- ann_gr
    S4Vectors::metadata(obj)$motifSites <- GenomicRanges::GRanges()
    obj
}

# ─── Basic return type ────────────────────────────────────────────────────────

test_that("plot_metagene: returns ggplot with line data for valid feature type", {
    obj <- .make_metagene_data()
    p <- plot_metagene(obj, feature = "gene")
    expect_s3_class(p, "ggplot")
    d <- p$data
    # bin_center values are in [0, 1]
    expect_true(all(d$bin_center >= 0))
    expect_true(all(d$bin_center <= 1))
    # sample_name matches fixture
    expect_setequal(as.character(d$sample_name), c("ctrl_1", "treat_1"))
    # mean_beta values are finite
    expect_true(all(is.finite(d$mean_beta)))
})

test_that("plot_metagene: mod_type filter produces identical data when all sites match", {
    obj <- .make_metagene_data()
    p <- plot_metagene(obj, feature = "gene")
    p_filt <- plot_metagene(obj, feature = "gene", mod_type = "6mA")
    expect_s3_class(p_filt, "ggplot")
    # All sites are 6mA, so filtering by 6mA should produce identical data
    d_all  <- p$data
    d_filt <- p_filt$data
    expect_identical(d_all$bin_center, d_filt$bin_center)
    expect_identical(as.character(d_all$sample_name),
                     as.character(d_filt$sample_name))
    expect_identical(d_all$mean_beta, d_filt$mean_beta)
})

test_that("plot_metagene: n_bins parameter controls bin count", {
    obj <- .make_metagene_data()
    p_20 <- plot_metagene(obj, feature = "gene", n_bins = 20L)
    p_50 <- plot_metagene(obj, feature = "gene", n_bins = 50L)
    # Compute expected bin centers the same way the function does
    breaks_20 <- seq(0, 1, length.out = 21L)
    centers_20 <- (breaks_20[-length(breaks_20)] + breaks_20[-1L]) / 2
    breaks_50 <- seq(0, 1, length.out = 51L)
    centers_50 <- (breaks_50[-length(breaks_50)] + breaks_50[-1L]) / 2
    # Every observed bin_center must be present in the expected grid
    obs_20 <- round(unique(p_20$data$bin_center), 10)
    obs_50 <- round(unique(p_50$data$bin_center), 10)
    expect_true(all(obs_20 %in% round(centers_20, 10)))
    expect_true(all(obs_50 %in% round(centers_50, 10)))
})

# ─── x-axis range ─────────────────────────────────────────────────────────────

test_that("plot_metagene: x-axis spans approximately [0, 1]", {
    obj <- .make_metagene_data()
    p <- plot_metagene(obj, feature = "gene")
    bd <- ggplot2::ggplot_build(p)
    x_vals <- bd$data[[1L]]$x
    expect_true(min(x_vals) >= 0 - 1e-6)
    expect_true(max(x_vals) <= 1 + 1e-6)
})

test_that("plot_metagene: extracts beta values from annotated row order", {
    obj <- .make_metagene_data()
    SummarizedExperiment::assay(obj, "methylation")[, "ctrl_1"] <-
        seq(0.1, 1.0, length.out = nrow(obj))
    SummarizedExperiment::assay(obj, "methylation")[, "treat_1"] <-
        seq(0.2, 0.9, length.out = nrow(obj))

    expected <- plot_metagene(obj, feature = "gene", n_bins = 10L)$data
    real_annotateSites <- annotateSites

    local_mocked_bindings(
        annotateSites = function(object, features = NULL, keep = "all", ...) {
            annotated <- real_annotateSites(
                object,
                features = features,
                keep = keep,
                ...
            )
            annotated[rev(seq_len(nrow(annotated))), ]
        },
        .package = "commaKit"
    )

    observed <- plot_metagene(obj, feature = "gene", n_bins = 10L)$data
    expect_equal(observed, expected)
})

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_metagene: error on non-commaData input", {
    expect_error(plot_metagene(data.frame()), "commaData")
})

test_that("plot_metagene: error when annotation is empty", {
    obj <- .make_metagene_data()
    S4Vectors::metadata(obj)$annotation <- GenomicRanges::GRanges()
    expect_error(plot_metagene(obj, feature = "gene"),
                 "annotation")
})

test_that("plot_metagene: error when feature type not found in annotation", {
    obj <- .make_metagene_data()
    expect_error(plot_metagene(obj, feature = "promoter"),
                 "promoter")
})

test_that("plot_metagene: error on invalid mod_type", {
    obj <- .make_metagene_data()
    expect_error(plot_metagene(obj, feature = "gene", mod_type = "4mC"),
                 "not found")
})

test_that("plot_metagene: error when n_bins < 2", {
    obj <- .make_metagene_data()
    expect_error(plot_metagene(obj, feature = "gene", n_bins = 1L),
                 "n_bins")
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_metagene: works with comma_example_data", {
    data(comma_example_data)
    p <- plot_metagene(comma_example_data, feature = "gene")
    expect_s3_class(p, "ggplot")
    expect_gt(nrow(p$data), 0L)
})
