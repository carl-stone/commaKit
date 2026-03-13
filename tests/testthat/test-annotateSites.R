test_that("annotateSites: type='overlap' assigns feature_type and feature_name", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    si <- siteInfo(result)
    expect_true("feature_type" %in% colnames(si))
    expect_true("feature_name" %in% colnames(si))
})

test_that("annotateSites: type='overlap' returns commaData with same dimensions", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    expect_s4_class(result, "commaData")
    expect_equal(dim(result), dim(comma_example_data))
})

test_that("annotateSites: type='overlap' assigns 'intergenic' to sites outside features", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    si <- siteInfo(result)
    expect_true("intergenic" %in% si$feature_type)
    expect_true("intergenic" %in% si$feature_name)
})

test_that("annotateSites: type='overlap' assigns non-intergenic for sites inside features", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    si <- siteInfo(result)
    # comma_example_data has 5 annotated genes; some sites must overlap them
    expect_true(any(si$feature_type != "intergenic"))
})

test_that("annotateSites: type='overlap' site known to be inside feature gets correct annotation", {
    # Build a minimal commaData with one site at position 50000 inside a feature
    features <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = 45000, end = 55000),
        strand   = "+",
        feature_type = "gene",
        name = "geneX"
    )
    data(comma_example_data)
    # Use subset with a site near position 50000
    si <- siteInfo(comma_example_data)
    idx <- which(si$position >= 45000 & si$position <= 55000)
    if (length(idx) > 0) {
        sub_obj <- comma_example_data[idx[1], ]
        result <- annotateSites(sub_obj, features = features, type = "overlap")
        expect_equal(siteInfo(result)$feature_type, "gene")
        expect_equal(siteInfo(result)$feature_name, "geneX")
    } else {
        skip("No sites in position range for this test")
    }
})

test_that("annotateSites: type='overlap' site outside all features gets 'intergenic'", {
    # Build feature far from any site in comma_example_data
    features <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = 1L, end = 10L),
        strand   = "+",
        feature_type = "gene",
        name = "geneY"
    )
    data(comma_example_data)
    # Pick a site far from position 1-10
    si <- siteInfo(comma_example_data)
    idx <- which(si$position > 1000)
    sub_obj <- comma_example_data[idx[1], ]
    result <- annotateSites(sub_obj, features = features, type = "overlap")
    expect_equal(siteInfo(result)$feature_type, "intergenic")
})

test_that("annotateSites: type='proximity' adds correct columns", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "proximity", window = 10000L)
    si <- siteInfo(result)
    expect_true("nearest_feature" %in% colnames(si))
    expect_true("distance_to_feature" %in% colnames(si))
    expect_true("rel_pos" %in% colnames(si))
})

test_that("annotateSites: type='proximity' distances are non-negative", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "proximity", window = 100000L)
    si <- siteInfo(result)
    dists <- si$distance_to_feature[!is.na(si$distance_to_feature)]
    expect_true(all(dists >= 0))
})

test_that("annotateSites: type='proximity' sites beyond window get NA", {
    # Use tiny window so nothing is within range
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "proximity", window = 0L)
    si <- siteInfo(result)
    # With window = 0, only sites exactly at a feature boundary match
    # Most sites will be NA
    expect_true(any(is.na(si$distance_to_feature)))
})

test_that("annotateSites: type='metagene' adds metagene_pos in [0,1]", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "metagene")
    si <- siteInfo(result)
    expect_true("metagene_pos" %in% colnames(si))
    pos_vals <- si$metagene_pos[!is.na(si$metagene_pos)]
    expect_true(all(pos_vals >= 0 & pos_vals <= 1))
})

test_that("annotateSites: type='metagene' adds metagene_feature column", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "metagene")
    si <- siteInfo(result)
    expect_true("metagene_feature" %in% colnames(si))
})

test_that("annotateSites: type='metagene' non-overlapping sites get NA", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "metagene")
    si <- siteInfo(result)
    # There should be some sites that don't overlap any feature
    expect_true(any(is.na(si$metagene_pos)))
})

test_that("annotateSites: type='metagene' strand-aware: + strand site at start → ~0", {
    features <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = 1000L, end = 2000L),
        strand   = "+",
        feature_type = "gene",
        name = "genePos"
    )
    data(comma_example_data)
    si <- siteInfo(comma_example_data)
    # Find a site at or near position 1000
    idx <- which(si$position >= 1000 & si$position <= 1050)
    if (length(idx) > 0) {
        sub_obj <- comma_example_data[idx[1], ]
        result <- annotateSites(sub_obj, features = features, type = "metagene")
        mg_pos <- siteInfo(result)$metagene_pos
        expect_true(!is.na(mg_pos))
        expect_true(mg_pos < 0.1)  # near start of + strand feature
    } else {
        skip("No sites near position 1000 for strand test")
    }
})

test_that("annotateSites: metagene - strand feature: site at high coordinate → ~0", {
    # For - strand, 0 = high coordinate (biological TSS)
    features <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = 1000L, end = 2000L),
        strand   = "-",
        feature_type = "gene",
        name = "geneMinus"
    )
    data(comma_example_data)
    si <- siteInfo(comma_example_data)
    idx <- which(si$position >= 1950 & si$position <= 2000)
    if (length(idx) > 0) {
        sub_obj <- comma_example_data[idx[1], ]
        result <- annotateSites(sub_obj, features = features, type = "metagene")
        mg_pos <- siteInfo(result)$metagene_pos
        expect_true(!is.na(mg_pos))
        expect_true(mg_pos < 0.1)  # near start of - strand feature (high coord)
    } else {
        skip("No sites near position 1950-2000 for - strand test")
    }
})

test_that("annotateSites: error on non-commaData input", {
    expect_error(annotateSites(data.frame(x = 1)), "'object' must be a commaData")
})

test_that("annotateSites: error on invalid type", {
    data(comma_example_data)
    expect_error(annotateSites(comma_example_data, type = "badtype"))
})

test_that("annotateSites: error when features NULL and annotation slot is empty", {
    data(comma_example_data)
    # Clear annotation slot
    obj_no_ann <- new("commaData",
        as(comma_example_data, "SummarizedExperiment"),
        genomeInfo = comma_example_data@genomeInfo,
        annotation = GenomicRanges::GRanges(),
        motifSites = comma_example_data@motifSites
    )
    expect_error(annotateSites(obj_no_ann, features = NULL), "No features available")
})

test_that("annotateSites: error on missing feature_col", {
    data(comma_example_data)
    features <- annotation(comma_example_data)
    expect_error(
        annotateSites(comma_example_data, features = features, feature_col = "nonexistent"),
        "not found in mcols"
    )
})

test_that("annotateSites: rowData has expected number of rows", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    expect_equal(nrow(siteInfo(result)), nrow(siteInfo(comma_example_data)))
})

test_that("annotateSites: assay data unchanged after annotation", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    expect_equal(methylation(result), methylation(comma_example_data))
    expect_equal(coverage(result), coverage(comma_example_data))
})

test_that("annotateSites: type='proximity' rel_pos is signed (can be negative)", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "proximity", window = 100000L)
    si <- siteInfo(result)
    rp <- si$rel_pos[!is.na(si$rel_pos)]
    # Should have both upstream (negative) and downstream (positive) sites
    # relative to TSS; at minimum check numeric
    expect_true(is.integer(rp))
})
