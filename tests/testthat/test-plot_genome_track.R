# Tests for plot_genome_track()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_track_data <- function() {
    n_sites   <- 10L
    positions <- seq(1000L, 10000L, by = 1000L)
    set.seed(3L)
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
        ranges   = IRanges::IRanges(start = c(2000L, 6000L),
                                     end   = c(4000L, 8000L)),
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

# ─── Data mapping ─────────────────────────────────────────────────────────────

test_that("plot_genome_track: p$data maps exact positions and beta values", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim", annotation = FALSE)
    expect_s3_class(p, "ggplot")
    # p$data should have position, beta, sample_name, mod_type columns
    expect_true("position" %in% colnames(p$data))
    expect_true("beta" %in% colnames(p$data))
    expect_true("sample_name" %in% colnames(p$data))
    # 10 sites * 2 samples = 20 rows
    expect_equal(nrow(p$data), 20L)
    # Positions should match the fixture exactly
    expect_equal(sort(unique(p$data$position)), seq(1000L, 10000L, by = 1000L))
    # Beta values should match the methylation assay
    methyl_mat <- methylation(obj)
    expect_equal(sort(p$data$beta), sort(as.vector(methyl_mat)))
})

test_that("plot_genome_track: annotation = FALSE returns ggplot without GeomRect layer", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim", annotation = FALSE)
    expect_s3_class(p, "ggplot")
    layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomPoint" %in% layer_classes)
    expect_false("GeomRect" %in% layer_classes)
})

test_that("plot_genome_track: mod_type filter to 6mA produces identical data (all sites are 6mA)", {
    obj <- .make_track_data()
    p_all <- plot_genome_track(obj, chromosome = "chr_sim", annotation = FALSE)
    p_filt <- plot_genome_track(obj, chromosome = "chr_sim",
                                mod_type = "6mA", annotation = FALSE)
    expect_s3_class(p_filt, "ggplot")
    # All sites are 6mA, so data should be identical
    expect_equal(nrow(p_filt$data), nrow(p_all$data))
    expect_equal(p_filt$data$position, p_all$data$position)
})

# ─── Positional filtering ─────────────────────────────────────────────────────

test_that("plot_genome_track: start/end filtering retains only positions in range", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim",
                           start = 1000L, end = 5000L, annotation = FALSE)
    expect_s3_class(p, "ggplot")
    # Positions in p$data should be exactly those in [1000, 5000]
    positions_in_range <- c(1000L, 2000L, 3000L, 4000L, 5000L)
    expect_equal(sort(unique(p$data$position)), positions_in_range)
    # Row count: 5 positions * 2 samples = 10
    expect_equal(nrow(p$data), 10L)
})

test_that("plot_genome_track: start only (no end) filters positions >= start", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim",
                           start = 3000L, annotation = FALSE)
    expect_s3_class(p, "ggplot")
    # Positions should be >= 3000
    expect_true(all(p$data$position >= 3000L))
    # Exact positions: 3000, 4000, 5000, ..., 10000
    expected_pos <- seq(3000L, 10000L, by = 1000L)
    expect_equal(sort(unique(p$data$position)), expected_pos)
    # 8 positions * 2 samples = 16 rows
    expect_equal(nrow(p$data), 16L)
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
    new_si <- GenomeInfoDb::Seqinfo(
        seqnames = c("chr_sim", "chr2"),
        seqlengths = c(100000L, 50000L),
        isCircular = c(FALSE, FALSE)
    )
    GenomeInfoDb::seqinfo(obj) <- new_si
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

test_that("plot_genome_track: comma_example_data has correct row count", {
    data(comma_example_data)
    p <- plot_genome_track(comma_example_data, chromosome = "chr_sim",
                           annotation = FALSE)
    expect_s3_class(p, "ggplot")
    # 588 sites * 6 samples = 3528 rows
    expect_equal(nrow(p$data), 588L * 6L)
})
