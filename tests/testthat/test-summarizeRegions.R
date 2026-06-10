.make_region_summary_fixture <- function() {
    beta <- matrix(
        c(0.2, 0.4, 0.6, 0.8,
          0.1, 0.5, 0.7, 0.9),
        nrow = 4L,
        ncol = 2L,
        dimnames = list(NULL, c("s1", "s2"))
    )
    coverage <- matrix(
        c(10L, 10L, 10L, 10L,
          20L, 20L, 20L, 20L),
        nrow = 4L,
        ncol = 2L,
        dimnames = dimnames(beta)
    )
    mod_counts <- matrix(
        c(2L, 4L, NA, 8L,
          2L, 10L, 14L, 18L),
        nrow = 4L,
        ncol = 2L,
        dimnames = dimnames(beta)
    )
    canonical_counts <- coverage - mod_counts
    sample_info <- data.frame(
        sample_name = c("s1", "s2"),
        condition = c("control", "treatment"),
        replicate = 1:2,
        stringsAsFactors = FALSE
    )
    .make_commaData_fixture(
        beta = beta,
        coverage = coverage,
        sample_info = sample_info,
        positions = c(100L, 200L, 300L, 1000L),
        mod_type = c("6mA", "6mA", "5mC", "6mA"),
        motif = c("GATC", "GATC", "CCWGG", "GATC"),
        mod_counts = mod_counts,
        canonical_counts = canonical_counts
    )
}

.make_regions <- function() {
    GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges = IRanges::IRanges(
            start = c(50L, 900L, 2000L),
            end = c(350L, 1100L, 2100L)
        ),
        region_label = c("first", "second", "empty")
    )
}

test_that("summarizeRegions() aggregates count evidence by region and sample", {
    obj <- .make_region_summary_fixture()
    regions <- .make_regions()

    out <- summarizeRegions(obj, regions)

    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), length(regions) * ncol(obj))
    first_a <- out[out$region_id == "region_1" & out$sample_name == "s1", ]
    expect_equal(first_a$n_sites, 2L)
    expect_equal(first_a$total_mod_counts, 6)
    expect_equal(first_a$total_canonical_counts, 14)
    expect_equal(first_a$total_valid_coverage, 20)
    expect_equal(first_a$region_methylation, 6 / 20)

    first_b <- out[out$region_id == "region_1" & out$sample_name == "s2", ]
    expect_equal(first_b$n_sites, 3L)
    expect_equal(first_b$total_mod_counts, 26)
    expect_equal(first_b$total_valid_coverage, 60)
    expect_equal(first_b$region_methylation, 26 / 60)
    expect_equal(first_b$region_region_label, "first")
})

test_that("summarizeRegions() applies min_sites and keeps empty regions", {
    obj <- .make_region_summary_fixture()
    regions <- .make_regions()

    out <- summarizeRegions(obj, regions, min_sites = 3L)

    first_a <- out[out$region_id == "region_1" & out$sample_name == "s1", ]
    expect_equal(first_a$n_sites, 2L)
    expect_true(is.na(first_a$region_methylation))

    empty_a <- out[out$region_id == "region_3" & out$sample_name == "s1", ]
    expect_equal(empty_a$n_sites, 0L)
    expect_equal(empty_a$total_mod_counts, 0)
    expect_equal(empty_a$total_valid_coverage, 0)
    expect_true(is.na(empty_a$region_methylation))
})

test_that("summarizeRegions() respects modification site filters", {
    obj <- .make_region_summary_fixture()
    regions <- .make_regions()

    out <- summarizeRegions(obj, regions, mod_type = "6mA")

    first_a <- out[out$region_id == "region_1" & out$sample_name == "s1", ]
    expect_equal(first_a$n_sites, 2L)
    expect_equal(first_a$total_mod_counts, 6)
    expect_equal(first_a$total_valid_coverage, 20)

    out_context <- summarizeRegions(obj, regions, mod_context = "5mC_CCWGG")
    first_b <- out_context[out_context$region_id == "region_1" &
                               out_context$sample_name == "s2", ]
    expect_equal(first_b$n_sites, 1L)
    expect_equal(first_b$total_mod_counts, 14)
    expect_equal(first_b$total_valid_coverage, 20)
})

test_that("summarizeRegions() requires count evidence assays", {
    obj <- .make_region_summary_fixture()
    SummarizedExperiment::assays(obj) <-
        SummarizedExperiment::assays(obj)[c("methylation", "coverage",
                                            "canonical_counts")]

    expect_error(
        summarizeRegions(obj, .make_regions()),
        "requires count evidence assays"
    )
})

test_that("summarizeRegions() validates inputs", {
    obj <- .make_region_summary_fixture()
    expect_error(summarizeRegions(obj, data.frame()), "GRanges")
    expect_error(summarizeRegions(obj, .make_regions(), min_sites = -1),
                 "min_sites")
})


test_that("summarizeRegions() validates motif and mod_context filters", {
    obj <- .make_region_summary_fixture()
    regions <- .make_regions()

    expect_error(summarizeRegions(obj, regions, motif = "BAD"),
                 "Requested motif")
    expect_error(summarizeRegions(obj, regions, mod_context = "6mA_BAD"),
                 "Requested mod_context")
})

test_that("summarizeRegions() returns a typed 0-row data frame for empty regions", {
    obj <- .make_region_summary_fixture()
    empty_regions <- GenomicRanges::GRanges(
        seqnames = character(0),
        ranges = IRanges::IRanges(start = integer(0), end = integer(0))
    )

    out <- summarizeRegions(obj, empty_regions)

    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 0L)
    expect_true(all(c("region_id", "sample_name", "n_sites",
                      "total_mod_counts", "total_valid_coverage",
                      "region_methylation", "total_canonical_counts") %in%
                    colnames(out)))
})

test_that("summarizeRegions() rejects assays with no usable count evidence", {
    obj <- .make_region_summary_fixture()
    SummarizedExperiment::assay(obj, "mod_counts")[] <- NA_integer_

    expect_error(
        summarizeRegions(obj, .make_regions()),
        "requires count evidence"
    )
})

test_that("summarizeRegions() excludes zero-coverage sites from usable counts", {
    obj <- .make_region_summary_fixture()
    SummarizedExperiment::assay(obj, "coverage")[1, 1] <- 0L
    SummarizedExperiment::assay(obj, "mod_counts")[1, 1] <- 0L
    SummarizedExperiment::assay(obj, "canonical_counts")[1, 1] <- 0L

    out <- summarizeRegions(obj, .make_regions())
    first_s1 <- out[out$region_id == "region_1" & out$sample_name == "s1", ]

    expect_equal(first_s1$n_sites, 1L)
    expect_equal(first_s1$total_mod_counts, 4)
    expect_equal(first_s1$total_valid_coverage, 10)
    expect_equal(first_s1$region_methylation, 4 / 10)
})

test_that("summarizeRegions() ignores strand for coordinate region summaries", {
    obj <- .make_region_summary_fixture()
    regions <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges = IRanges::IRanges(start = 50L, end = 350L),
        strand = "-"
    )

    out <- summarizeRegions(obj, regions, mod_type = "6mA")
    first_s1 <- out[out$region_id == "region_1" & out$sample_name == "s1", ]

    expect_equal(first_s1$n_sites, 2L)
    expect_equal(first_s1$total_mod_counts, 6)
})

test_that("summarizeRegions() reports optional component totals conservatively", {
    obj <- .make_region_summary_fixture()
    SummarizedExperiment::assay(obj, "canonical_counts")[2, 1] <- NA_integer_
    SummarizedExperiment::assay(obj, "other_mod_counts") <- matrix(
        c(1L, NA, 3L, 4L,
          2L, 2L, 2L, 2L),
        nrow = nrow(obj),
        ncol = ncol(obj),
        dimnames = dimnames(SummarizedExperiment::assay(obj, "coverage"))
    )

    out <- summarizeRegions(obj, .make_regions())
    first_s1 <- out[out$region_id == "region_1" & out$sample_name == "s1", ]
    first_s2 <- out[out$region_id == "region_1" & out$sample_name == "s2", ]

    expect_true(is.na(first_s1$total_canonical_counts))
    expect_true(is.na(first_s1$total_other_mod_counts))
    expect_equal(first_s2$total_other_mod_counts, 6)
})

test_that("summarizeRegions() rejects non-integer min_sites", {
    obj <- .make_region_summary_fixture()
    expect_error(summarizeRegions(obj, .make_regions(), min_sites = 1.5),
                 "min_sites")
})
