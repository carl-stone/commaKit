# Tests for plot_metagene()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_metagene_data <- function() {
    # 10 sites distributed within annotated genes so metagene has overlap
    n_sites   <- 10L
    positions <- c(1500L, 2000L, 2500L, 3000L, 3500L,
                   6500L, 7000L, 7500L, 8000L, 8500L)
    site_keys <- paste0("chr_sim:", positions, ":+:6mA")
    set.seed(5L)
    betas <- matrix(
        runif(n_sites * 2L, 0.1, 0.9),
        nrow = n_sites, ncol = 2L,
        dimnames = list(site_keys, c("ctrl_1", "treat_1"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 2L,
                      dimnames = dimnames(betas))
    rd <- S4Vectors::DataFrame(
        chrom    = rep("chr_sim", n_sites),
        position = positions,
        strand   = rep("+", n_sites),
        mod_type = rep("6mA", n_sites),
        row.names = site_keys
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

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = betas, coverage = cov_mat),
        rowData = rd,
        colData = cd
    )
    new("commaData", se,
        genomeInfo = c(chr_sim = 100000L),
        annotation = ann_gr,
        motifSites = GenomicRanges::GRanges())
}

# ─── Basic return type ────────────────────────────────────────────────────────

test_that("plot_metagene: returns ggplot for valid feature type", {
    obj <- .make_metagene_data()
    p <- plot_metagene(obj, feature = "gene")
    expect_s3_class(p, "ggplot")
})

test_that("plot_metagene: mod_type filter accepted", {
    obj <- .make_metagene_data()
    p <- plot_metagene(obj, feature = "gene", mod_type = "6mA")
    expect_s3_class(p, "ggplot")
})

test_that("plot_metagene: n_bins parameter accepted", {
    obj <- .make_metagene_data()
    p <- plot_metagene(obj, feature = "gene", n_bins = 20L)
    expect_s3_class(p, "ggplot")
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

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_metagene: error on non-commaData input", {
    expect_error(plot_metagene(data.frame()), "commaData")
})

test_that("plot_metagene: error when annotation is empty", {
    obj <- .make_metagene_data()
    obj@annotation <- GenomicRanges::GRanges()
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
})
