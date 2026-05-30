# Tests for plot_heatmap()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_heatmap_fixtures <- function() {
    n_sites   <- 15L
    positions <- seq(1000L, 15000L, by = 1000L)
    set.seed(10L)
    betas <- matrix(
        runif(n_sites * 3L, 0.0, 1.0),
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
    obj <- new("commaData", rse)

    set.seed(11L)
    res <- data.frame(
        chrom         = rep("chr_sim", n_sites),
        position      = positions,
        strand        = rep("+", n_sites),
        mod_type      = rep("6mA", n_sites),
        motif         = rep("GATC", n_sites),
        dm_pvalue     = runif(n_sites, 0, 0.1),
        dm_padj       = runif(n_sites, 0, 0.05),
        dm_delta_beta = runif(n_sites, -0.5, 0.5),
        stringsAsFactors = FALSE
    )
    list(obj = obj, res = res)
}

# ─── Data mapping ─────────────────────────────────────────────────────────────

test_that("plot_heatmap: p$data maps site_key, sample_name, and beta exactly", {
    fix <- .make_heatmap_fixtures()
    p <- plot_heatmap(fix$res, fix$obj)
    expect_s3_class(p, "ggplot")
    # p$data should have site_key, sample_name, beta columns
    expect_true("site_key" %in% colnames(p$data))
    expect_true("sample_name" %in% colnames(p$data))
    expect_true("beta" %in% colnames(p$data))
    # 15 sites * 3 samples = 45 rows (all padj non-NA)
    n_nonNA <- sum(!is.na(fix$res$dm_padj))
    expect_equal(nrow(p$data), n_nonNA * ncol(methylation(fix$obj)))
    # Sample names should match the object
    expect_equal(sort(unique(p$data$sample_name)), c("ctrl_1", "ctrl_2", "treat_1"))
})

test_that("plot_heatmap: n_sites limits p$data to top N sites by padj", {
    fix <- .make_heatmap_fixtures()
    p5 <- plot_heatmap(fix$res, fix$obj, n_sites = 5L)
    # 5 sites * 3 samples = 15 rows
    expect_equal(nrow(p5$data), 5L * 3L)
    # site_key should have exactly 5 unique values
    expect_equal(length(unique(p5$data$site_key)), 5L)
})

test_that("plot_heatmap: n_sites larger than available sites clamps to all available", {
    fix <- .make_heatmap_fixtures()
    p <- plot_heatmap(fix$res, fix$obj, n_sites = 1000L)
    expect_s3_class(p, "ggplot")
    # All 15 sites should be shown
    n_nonNA <- sum(!is.na(fix$res$dm_padj))
    expect_equal(length(unique(p$data$site_key)), n_nonNA)
    expect_equal(nrow(p$data), n_nonNA * 3L)
})

test_that("plot_heatmap: subset and sorted results use preserved rowname indices", {
    fix <- .make_heatmap_fixtures()

    subset_res <- fix$res[c(10L, 2L, 7L, 5L), , drop = FALSE]
    subset_res$dm_padj <- c(0.03, 0.01, 0.02, 0.04)
    subset_res$dm_delta_beta <- c(0.4, -0.2, 0.1, -0.5)

    p <- plot_heatmap(
        subset_res,
        fix$obj,
        n_sites = 3L,
        annotation_cols = character(0)
    )

    selected <- subset_res[order(subset_res$dm_padj), , drop = FALSE]
    selected <- selected[seq_len(3L), , drop = FALSE]
    selected <- selected[order(selected$dm_delta_beta), , drop = FALSE]

    expected_mat <- methylation(fix$obj)[as.integer(rownames(selected)), ,
                                         drop = FALSE]
    expected_site_keys <- paste(
        selected$chrom,
        selected$position,
        selected$strand,
        selected$mod_type,
        selected$motif,
        sep = ":"
    )

    expect_equal(
        as.character(p$data$site_key),
        rep(expected_site_keys, times = ncol(expected_mat))
    )
    expect_equal(
        as.character(p$data$sample_name),
        rep(colnames(expected_mat), each = nrow(expected_mat))
    )
    expect_equal(p$data$beta, as.vector(expected_mat))
})

# ─── NA handling ─────────────────────────────────────────────────────────────

test_that("plot_heatmap: NA beta values are preserved in p$data", {
    fix <- .make_heatmap_fixtures()
    methyl_mat <- methylation(fix$obj)
    methyl_mat[1L, "ctrl_1"] <- NA
    SummarizedExperiment::assay(fix$obj, "methylation") <- methyl_mat
    p <- plot_heatmap(fix$res, fix$obj)
    expect_s3_class(p, "ggplot")
    # NA beta should be present in p$data for that site-sample combo
    na_rows <- is.na(p$data$beta)
    expect_true(any(na_rows))
    # The NA should be for the site we injected it into
    na_site <- p$data$site_key[na_rows]
    na_samp <- p$data$sample_name[na_rows]
    expect_equal(na_samp, "ctrl_1")
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

test_that("plot_heatmap: error when selected rownames are blank", {
    fix <- .make_heatmap_fixtures()
    bad_res <- fix$res
    bad_res$dm_padj <- seq(0.001, 0.015, length.out = nrow(bad_res))
    attr(bad_res, "row.names") <- c("", as.character(seq(2L, nrow(bad_res))))

    expect_error(
        plot_heatmap(
            bad_res,
            fix$obj,
            n_sites = 1L,
            annotation_cols = character(0)
        ),
        "row names must be non-missing integer row indices"
    )
})

test_that("plot_heatmap: error when selected rownames are missing", {
    fix <- .make_heatmap_fixtures()
    bad_res <- fix$res
    bad_res$dm_padj <- seq(0.001, 0.015, length.out = nrow(bad_res))
    attr(bad_res, "row.names") <- c(NA_character_,
                                    as.character(seq(2L, nrow(bad_res))))

    expect_error(
        plot_heatmap(
            bad_res,
            fix$obj,
            n_sites = 1L,
            annotation_cols = character(0)
        ),
        "row names must be non-missing integer row indices"
    )
})

test_that("plot_heatmap: error when selected rownames are non-integer", {
    fix <- .make_heatmap_fixtures()
    bad_res <- fix$res
    bad_res$dm_padj <- seq(0.001, 0.015, length.out = nrow(bad_res))
    rownames(bad_res)[1L] <- "chr_sim:1000:+:6mA:GATC"

    expect_error(
        plot_heatmap(
            bad_res,
            fix$obj,
            n_sites = 1L,
            annotation_cols = character(0)
        ),
        "Non-integer row name"
    )
})

test_that("plot_heatmap: error when selected rownames are out of range", {
    fix <- .make_heatmap_fixtures()
    bad_res <- fix$res
    bad_res$dm_padj <- seq(0.001, 0.015, length.out = nrow(bad_res))
    rownames(bad_res)[1L] <- as.character(nrow(fix$obj) + 1L)

    expect_error(
        plot_heatmap(
            bad_res,
            fix$obj,
            n_sites = 1L,
            annotation_cols = character(0)
        ),
        "Out-of-range row name"
    )
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_heatmap: comma_example_data has correct row count", {
    data(comma_example_data)
    cd_dm <- diffMethyl(comma_example_data, ~ condition, mod_type = "6mA")
    res <- results(cd_dm)
    p <- plot_heatmap(res, cd_dm, n_sites = 20L)
    expect_s3_class(p, "ggplot")
    # 20 sites * 6 samples = 120 rows
    expect_equal(nrow(p$data), 20L * 6L)
})
