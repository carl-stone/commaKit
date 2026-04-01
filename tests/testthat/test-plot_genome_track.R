# Tests for plot_genome_track()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_track_data <- function() {
    n_sites   <- 10L
    positions <- seq(1000L, 10000L, by = 1000L)
    site_keys <- paste0("chr_sim:", positions, ":+:6mA:GATC")
    set.seed(3L)
    betas <- matrix(
        runif(n_sites * 2L, 0.1, 0.9),
        nrow = n_sites, ncol = 2L,
        dimnames = list(site_keys, c("ctrl_1", "treat_1"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 2L,
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
        sample_name = c("ctrl_1", "treat_1"),
        condition   = c("control", "treatment"),
        replicate   = 1:2,
        row.names   = c("ctrl_1", "treat_1")
    )
    ann_gr <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = c(2000L, 6000L),
                                     end   = c(4000L, 8000L)),
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

test_that("plot_genome_track: returns ggplot for valid chromosome", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim")
    # May return a patchwork or ggplot depending on patchwork availability
    expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("plot_genome_track: annotation = FALSE suppresses annotation track", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim", annotation = FALSE)
    expect_s3_class(p, "ggplot")
})

test_that("plot_genome_track: mod_type filter returns ggplot without error", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim",
                           mod_type = "6mA", annotation = FALSE)
    expect_s3_class(p, "ggplot")
})

# ─── Positional filtering ─────────────────────────────────────────────────────

test_that("plot_genome_track: start/end filtering reduces displayed sites", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim",
                           start = 1000L, end = 5000L, annotation = FALSE)
    expect_s3_class(p, "ggplot")
    bd <- ggplot2::ggplot_build(p)
    # All x values in layer 1 (points) should be within [1000, 5000]
    x_vals <- bd$data[[1L]]$x
    expect_true(all(x_vals >= 1000L & x_vals <= 5000L))
})

test_that("plot_genome_track: start only (no end) accepted", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim",
                           start = 3000L, annotation = FALSE)
    expect_s3_class(p, "ggplot")
})

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_genome_track: error on non-commaData input", {
    expect_error(plot_genome_track(data.frame(), chromosome = "chr_sim"),
                 "commaData")
})

test_that("plot_genome_track: error when chromosome not in genome", {
    obj <- .make_track_data()
    expect_error(plot_genome_track(obj, chromosome = "chrX"),
                 "chrX")
})

test_that("plot_genome_track: error when no sites on chromosome", {
    obj <- .make_track_data()
    # Add chr2 to genomeInfo but not to data
    new_gi <- c(chr_sim = 100000L, chr2 = 50000L)
    obj@genomeInfo <- new_gi
    expect_error(plot_genome_track(obj, chromosome = "chr2"),
                 "No methylation sites found")
})

test_that("plot_genome_track: error when start > end", {
    obj <- .make_track_data()
    expect_error(
        plot_genome_track(obj, chromosome = "chr_sim",
                          start = 5000L, end = 1000L),
        "'end' must be"
    )
})

test_that("plot_genome_track: error on invalid mod_type", {
    obj <- .make_track_data()
    expect_error(plot_genome_track(obj, chromosome = "chr_sim",
                                   mod_type = "4mC"),
                 "not found")
})

test_that("plot_genome_track: error when invalid annotation argument", {
    obj <- .make_track_data()
    expect_error(
        plot_genome_track(obj, chromosome = "chr_sim", annotation = 42),
        "'annotation' must be"
    )
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_genome_track: works with comma_example_data", {
    data(comma_example_data)
    p <- plot_genome_track(comma_example_data, chromosome = "chr_sim",
                           annotation = FALSE)
    expect_s3_class(p, "ggplot")
})
