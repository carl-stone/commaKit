## Tests for plot_tss_profile()

## ── Inline test fixtures ────────────────────────────────────────────────────

## 12 sites on chr_sim distributed around 3 gene TSS positions.
## Gene 1 (+): TSS at position 1000 → sites at 700, 900, 1000, 1100, 1300
## Gene 2 (+): TSS at position 5000 → sites at 4600, 4800, 5050, 5200
## Gene 3 (-): TSS at position 9000 → sites at 8600, 8900, 9100
## (For - strand, TSS = end of feature: gene3 spans 8000–9000, so TSS = 9000)

.make_tss_data <- function() {
    positions <- c(700L, 900L, 1000L, 1100L, 1300L,
                   4600L, 4800L, 5050L, 5200L,
                   8600L, 8900L, 9100L)
    n_sites   <- length(positions)
    set.seed(1L)
    beta_ctrl  <- round(runif(n_sites, 0.5, 1.0), 3)
    beta_treat <- round(runif(n_sites, 0.0, 0.5), 3)

    methyl_mat <- matrix(
        c(beta_ctrl, beta_treat),
        nrow = n_sites, ncol = 2,
        dimnames = list(NULL, c("ctrl_1", "treat_1"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 2,
                      dimnames = list(NULL, c("ctrl_1", "treat_1")))

    site_gr <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = positions, width = 1L),
        strand   = "+",
        mod_type    = factor("6mA", levels = c("4mC", "5mC", "6mA")),
        motif       = "GATC"
    )
    GenomeInfoDb::seqinfo(site_gr) <- GenomeInfoDb::Seqinfo(
        seqnames = "chr_sim",
        seqlengths = 100000L,
        isCircular = FALSE
    )
    cd <- S4Vectors::DataFrame(
        sample_name = c("ctrl_1", "treat_1"),
        condition   = c("control", "treatment"),
        replicate   = c(1L, 1L)
    )
    rownames(cd) <- c("ctrl_1", "treat_1")

    rse <- SummarizedExperiment::SummarizedExperiment(
        assays     = list(methylation = methyl_mat, coverage = cov_mat),
        rowRanges  = site_gr,
        colData    = cd
    )

    ## 3 genes (two + strand, one - strand)
    annot_gr <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(
            start = c(800L,  4500L, 8000L),
            end   = c(2000L, 6000L, 9000L)
        ),
        strand = c("+", "+", "-")
    )
    annot_gr$feature_type <- c("gene", "gene", "gene")
    annot_gr$name         <- c("geneA", "geneB", "geneC")

    obj <- new("commaData", rse)
    S4Vectors::metadata(obj)$annotation <- annot_gr
    S4Vectors::metadata(obj)$motifSites <- GenomicRanges::GRanges()
    obj
}

## Same as above but with regulatory features added to annotation
.make_tss_data_with_regulatory <- function() {
    obj <- .make_tss_data()
    annot_gr <- annotation(obj)

    ## Two sigma_binding features overlapping some sites
    reg_gr <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = c(880L, 4750L), end = c(920L, 4810L)),
        strand   = "+"
    )
    reg_gr$feature_type <- c("sigma_binding", "sigma_binding")
    reg_gr$name         <- c("sigma1", "sigma2")

    combined_gr <- c(annot_gr, reg_gr)
    obj2 <- new("commaData", SummarizedExperiment::SummarizedExperiment(
            assays     = list(methylation = methylation(obj), coverage = coverage(obj)),
            rowRanges  = rowRanges(obj),
            colData    = sampleInfo(obj)
        ))
    S4Vectors::metadata(obj2) <- S4Vectors::metadata(obj)
    S4Vectors::metadata(obj2)$annotation <- combined_gr
    obj2
}

## ── Basic return type ────────────────────────────────────────────────────────

test_that("returns ggplot for valid input", {
    obj <- .make_tss_data()
    p   <- plot_tss_profile(obj, feature_type = "gene", window = 500L)
    expect_s3_class(p, "ggplot")
    # Verify p$data has rel_pos within [-window, window] and beta values present
    expect_true("rel_pos" %in% colnames(p$data))
    expect_true(all(p$data$rel_pos >= -500L & p$data$rel_pos <= 500L))
    expect_true("beta" %in% colnames(p$data))
    expect_true(any(!is.na(p$data$beta)))
})

test_that("returns ggplot when mod_type is specified", {
    obj <- .make_tss_data()
    p_unfiltered <- plot_tss_profile(obj, feature_type = "gene")
    p_filtered   <- plot_tss_profile(obj, feature_type = "gene", mod_type = "6mA")
    expect_s3_class(p_filtered, "ggplot")
    # Filtering by mod_type produces same or fewer data rows
    expect_true(nrow(p_filtered$data) <= nrow(p_unfiltered$data))
})

## ── x-axis and TSS marker ───────────────────────────────────────────────────

test_that("x-axis data values are within [-window, +window]", {
    obj    <- .make_tss_data()
    win    <- 500L
    p      <- plot_tss_profile(obj, feature_type = "gene", window = win)
    ld     <- ggplot2::ggplot_build(p)$data[[1]]
    x_vals <- ld$x[!is.na(ld$x)]
    expect_true(length(x_vals) > 0L)
    expect_true(all(x_vals >= -win & x_vals <= win))
})

test_that("geom_vline at x = 0 is present", {
    obj   <- .make_tss_data()
    p     <- plot_tss_profile(obj, feature_type = "gene", window = 500L)
    ## Find geom_vline layer(s)
    layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true(any(layer_classes == "GeomVline"))
    ## Confirm the intercept is 0
    vline_data <- p$layers[[which(layer_classes == "GeomVline")[1]]]$data
    expect_equal(vline_data$xintercept, 0L)
})

## ── color_by modes ───────────────────────────────────────────────────────────

test_that("color_by = 'mod_type' returns ggplot", {
    obj <- .make_tss_data()
    p   <- plot_tss_profile(obj, feature_type = "gene", color_by = "mod_type")
    expect_s3_class(p, "ggplot")
    # Verify colour mapping references mod_type
    expect_true("colour" %in% names(p$mapping))
    expect_equal(p$labels$colour, "Modification type")
})

test_that("color_by = 'regulatory_element' with valid types returns ggplot", {
    obj <- .make_tss_data_with_regulatory()
    p   <- plot_tss_profile(obj, feature_type = "gene",
                             color_by = "regulatory_element",
                             regulatory_feature_types = "sigma_binding")
    expect_s3_class(p, "ggplot")
    # Verify regulatory_element column present in plot data
    expect_true("regulatory_element" %in% colnames(p$data))
    # Verify colour mapping is present
    expect_true("colour" %in% names(p$mapping))
})

test_that("color_by = 'regulatory_element' with absent types falls back to sample", {
    obj <- .make_tss_data()   # no regulatory features
    expect_message(
        p <- plot_tss_profile(obj, feature_type = "gene",
                              color_by = "regulatory_element",
                              regulatory_feature_types = "not_a_real_type"),
        "Falling back to color_by = 'sample'"
    )
    expect_s3_class(p, "ggplot")
})

## ── facet_by modes ───────────────────────────────────────────────────────────

test_that("facet_by = 'sample' produces FacetWrap layer", {
    obj <- .make_tss_data()
    p   <- plot_tss_profile(obj, feature_type = "gene", facet_by = "sample")
    expect_s3_class(p, "ggplot")
    expect_s3_class(p$facet, "FacetWrap")
})

test_that("facet_by = 'mod_type' produces FacetWrap layer", {
    obj <- .make_tss_data()
    p   <- plot_tss_profile(obj, feature_type = "gene", facet_by = "mod_type")
    expect_s3_class(p, "ggplot")
    expect_s3_class(p$facet, "FacetWrap")
})

## ── smooth overlay ───────────────────────────────────────────────────────────

test_that("show_smooth = TRUE adds a geom_line layer", {
    obj          <- .make_tss_data()
    p            <- suppressWarnings(
        plot_tss_profile(obj, feature_type = "gene",
                         window = 500L, show_smooth = TRUE,
                         smooth_span = 0.5)
    )
    expect_s3_class(p, "ggplot")
    layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true(any(layer_classes == "GeomLine"))
})

test_that("show_smooth = TRUE warns on numerically unstable LOESS fit", {
    obj <- .make_tss_data()
    expect_warning(
        plot_tss_profile(obj, feature_type = "gene",
                         window = 500L, show_smooth = TRUE,
                         smooth_span = 0.5),
        "LOESS.*numerical instability"
    )
})

test_that("show_smooth = TRUE with < 10 pts per group warns but still plots", {
    ## Make object with very few sites so each sample group has < 10 points
    positions  <- c(950L, 1050L)
    beta_mat   <- matrix(c(0.8, 0.9, 0.3, 0.2), nrow = 2,
                         dimnames = list(NULL, c("ctrl_1", "treat_1")))
    cov_mat    <- matrix(20L, nrow = 2, ncol = 2,
                         dimnames = list(NULL, c("ctrl_1", "treat_1")))
    site_gr <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = positions, width = 1L),
        strand   = "+",
        mod_type    = factor("6mA", levels = c("4mC", "5mC", "6mA")),
        motif       = "GATC"
    )
    GenomeInfoDb::seqinfo(site_gr) <- GenomeInfoDb::Seqinfo(
        seqnames = "chr_sim",
        seqlengths = 100000L,
        isCircular = FALSE
    )
    cd <- S4Vectors::DataFrame(sample_name = c("ctrl_1", "treat_1"),
                                condition   = c("control", "treatment"),
                                replicate   = c(1L, 1L))
    rownames(cd) <- c("ctrl_1", "treat_1")
    rse <- SummarizedExperiment::SummarizedExperiment(
        assays     = list(methylation = beta_mat, coverage = cov_mat),
        rowRanges  = site_gr, colData = cd
    )
    annot_gr <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = 800L, end = 2000L),
        strand   = "+"
    )
    annot_gr$feature_type <- "gene"
    annot_gr$name         <- "geneA"
    tiny_obj <- new("commaData", rse)
    S4Vectors::metadata(tiny_obj)$annotation <- annot_gr
    S4Vectors::metadata(tiny_obj)$motifSites <- GenomicRanges::GRanges()

    expect_warning(
        p <- plot_tss_profile(tiny_obj, feature_type = "gene",
                              window = 500L, show_smooth = TRUE),
        "Fewer than 10 data points"
    )
    expect_s3_class(p, "ggplot")
})

test_that("color_by = 'none' returns ggplot with no colour aesthetic", {
    obj <- .make_tss_data()
    p   <- plot_tss_profile(obj, feature_type = "gene", color_by = "none")
    expect_s3_class(p, "ggplot")
    ## The global aes should not contain a colour mapping
    global_aes_names <- names(p$mapping)
    expect_false("colour" %in% global_aes_names)
})

test_that("color_by = 'none' + facet_by = 'mod_type' + show_smooth returns ggplot", {
    obj <- .make_tss_data()
    p   <- suppressWarnings(
        plot_tss_profile(obj, feature_type = "gene",
                         color_by = "none", facet_by = "mod_type",
                         show_smooth = TRUE, smooth_span = 0.5)
    )
    expect_s3_class(p, "ggplot")
    # Verify no colour aesthetic
    expect_false("colour" %in% names(p$mapping))
    # Verify faceted
    expect_s3_class(p$facet, "FacetWrap")
    # Verify smooth layer present
    layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true(any(layer_classes == "GeomLine"))
})

test_that("show_smooth = TRUE with facet_by = 'sample' produces separate smooth per panel", {
    obj <- .make_tss_data()
    p   <- suppressWarnings(
        plot_tss_profile(obj, feature_type = "gene",
                         facet_by = "sample", show_smooth = TRUE,
                         smooth_span = 0.5)
    )
    expect_s3_class(p, "ggplot")
    ## Find the geom_line layer (smooth)
    layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    line_idx <- which(layer_classes == "GeomLine")[1]
    expect_true(!is.na(line_idx))
    ## The smooth data must include sample_name so ggplot can route per-panel
    smooth_data <- p$layers[[line_idx]]$data
    expect_true("sample_name" %in% names(smooth_data))
    ## One row set per sample
    expect_gte(length(unique(smooth_data$sample_name)), 1L)
})

## ── Error conditions ─────────────────────────────────────────────────────────

test_that("error on non-commaData input", {
    expect_error(plot_tss_profile(list()), "commaData")
})

test_that("error on empty annotation", {
    obj <- .make_tss_data()
    ## Replace annotation with empty GRanges
    obj2 <- new("commaData",
                SummarizedExperiment::SummarizedExperiment(
                    assays     = list(methylation = methylation(obj),
                                   coverage    = coverage(obj)),
                    rowRanges  = rowRanges(obj),
                    colData    = sampleInfo(obj)
                ))
    S4Vectors::metadata(obj2) <- S4Vectors::metadata(obj)
    S4Vectors::metadata(obj2)$annotation <- GenomicRanges::GRanges()
    expect_error(plot_tss_profile(obj2), "annotation\\(object\\) is empty")
})

test_that("error when feature_type not in annotation", {
    obj <- .make_tss_data()
    expect_error(
        plot_tss_profile(obj, feature_type = "not_a_feature"),
        "not_a_feature"
    )
})

test_that("error on invalid mod_type", {
    obj <- .make_tss_data()
    expect_error(
        plot_tss_profile(obj, mod_type = "bogus"),
        "not found"
    )
})

test_that("error when color_by = 'regulatory_element' without regulatory_feature_types", {
    obj <- .make_tss_data()
    expect_error(
        plot_tss_profile(obj, color_by = "regulatory_element"),
        "regulatory_feature_types"
    )
})

test_that("error when no sites fall within window", {
    obj <- .make_tss_data()
    ## Use a tiny window that catches nothing (sites are >= 200 bp from TSS)
    expect_error(
        plot_tss_profile(obj, feature_type = "gene", window = 1L),
        "No methylation sites found within"
    )
})

test_that("error on window < 1", {
    obj <- .make_tss_data()
    expect_error(plot_tss_profile(obj, window = 0L), "'window' must be a positive integer")
})

test_that("error on alpha outside (0, 1]", {
    obj <- .make_tss_data()
    expect_error(plot_tss_profile(obj, alpha = 0),   "'alpha'")
    expect_error(plot_tss_profile(obj, alpha = 1.5), "'alpha'")
})

## ── Integration with comma_example_data ─────────────────────────────────────

test_that("works with comma_example_data, feature_type = 'gene'", {
    data(comma_example_data)
    p <- plot_tss_profile(comma_example_data, feature_type = "gene",
                          window = 500L)
    expect_s3_class(p, "ggplot")
    # Verify plot data has expected structure
    expect_true(all(c("rel_pos", "beta", "sample_name") %in% colnames(p$data)))
    expect_true(nrow(p$data) > 0L)
})
