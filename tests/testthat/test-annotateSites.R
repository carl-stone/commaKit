# ── overlap mode ──────────────────────────────────────────────────────────────

test_that("annotateSites: type='overlap' adds feature_types and feature_names columns", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    si <- siteInfo(result)
    expect_true("feature_types" %in% colnames(si))
    expect_true("feature_names" %in% colnames(si))
})

test_that("annotateSites: type='overlap' feature_types and feature_names are CharacterList", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    si <- siteInfo(result)
    expect_true(is(si$feature_types, "CharacterList"))
    expect_true(is(si$feature_names, "CharacterList"))
})

test_that("annotateSites: type='overlap' CharacterList has length equal to nrow(siteInfo)", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    si <- siteInfo(result)
    expect_equal(length(si$feature_types), nrow(si))
    expect_equal(length(si$feature_names), nrow(si))
})

test_that("annotateSites: type='overlap' returns commaData with same dimensions", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    expect_s4_class(result, "commaData")
    expect_equal(dim(result), dim(comma_example_data))
})

test_that("annotateSites: type='overlap' intergenic sites have length-0 list elements", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    si <- siteInfo(result)
    # comma_example_data has 5 genes covering a portion of chr_sim; some sites
    # will not overlap any feature
    expect_true(any(lengths(si$feature_types) == 0))
    expect_true(any(lengths(si$feature_names) == 0))
})

test_that("annotateSites: type='overlap' genic sites have length > 0 list elements", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "overlap")
    si <- siteInfo(result)
    expect_true(any(lengths(si$feature_types) > 0))
})

test_that("annotateSites: type='overlap' site inside feature gets correct type in CharacterList", {
    # Build a feature covering position 45000-55000 on chr_sim
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 45000L, end = 55000L),
        strand       = "+",
        feature_type = "gene",
        name         = "geneX"
    )
    data(comma_example_data)
    si <- siteInfo(comma_example_data)
    idx <- which(si$position >= 45000L & si$position <= 55000L)
    if (length(idx) == 0L) skip("No sites in 45000-55000 range")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, type = "overlap")
    res_si  <- siteInfo(result)
    expect_true("gene" %in% res_si$feature_types[[1L]])
    expect_true("geneX" %in% res_si$feature_names[[1L]])
})

test_that("annotateSites: type='overlap' site outside all features gets length-0 element", {
    # Feature covers only positions 1-10; pick a site far from it
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 1L, end = 10L),
        strand       = "+",
        feature_type = "gene",
        name         = "geneY"
    )
    data(comma_example_data)
    si  <- siteInfo(comma_example_data)
    idx <- which(si$position > 1000L)
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, type = "overlap")
    expect_equal(length(siteInfo(result)$feature_types[[1L]]), 0L)
    expect_equal(length(siteInfo(result)$feature_names[[1L]]), 0L)
})

test_that("annotateSites: type='overlap' site overlapping 2 features gets both in CharacterList", {
    # Two features both covering position 50000
    features <- GenomicRanges::GRanges(
        seqnames     = c("chr_sim", "chr_sim"),
        ranges       = IRanges::IRanges(start = c(45000L, 48000L),
                                        end   = c(55000L, 52000L)),
        strand       = c("+", "+"),
        feature_type = c("gene", "TF_binding_site"),
        name         = c("geneA", "tfbsB")
    )
    data(comma_example_data)
    si  <- siteInfo(comma_example_data)
    idx <- which(si$position >= 48000L & si$position <= 52000L)
    if (length(idx) == 0L) skip("No sites in 48000-52000 range for multi-overlap test")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, type = "overlap")
    res_si  <- siteInfo(result)
    expect_equal(length(res_si$feature_types[[1L]]), 2L)
    expect_true("gene" %in% res_si$feature_types[[1L]])
    expect_true("TF_binding_site" %in% res_si$feature_types[[1L]])
})

# ── proximity mode ────────────────────────────────────────────────────────────

test_that("annotateSites: type='proximity' adds correct list-typed columns", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "proximity", window = 10000L)
    si <- siteInfo(result)
    expect_true("nearby_features" %in% colnames(si))
    expect_true("distances_to_features" %in% colnames(si))
    expect_true("rel_positions" %in% colnames(si))
    expect_true(is(si$nearby_features,       "CharacterList"))
    expect_true(is(si$distances_to_features, "IntegerList"))
    expect_true(is(si$rel_positions,         "IntegerList"))
})

test_that("annotateSites: type='proximity' distances are non-negative", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "proximity", window = 100000L)
    si     <- siteInfo(result)
    all_dists <- unlist(si$distances_to_features)
    expect_true(all(all_dists >= 0L))
})

test_that("annotateSites: type='proximity' all distances are within window", {
    window <- 5000L
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "proximity", window = window)
    si     <- siteInfo(result)
    all_dists <- unlist(si$distances_to_features)
    expect_true(all(all_dists <= window))
})

test_that("annotateSites: type='proximity' sites beyond window get length-0 elements", {
    # window = 0 means only sites exactly inside a feature boundary match
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "proximity", window = 0L)
    si     <- siteInfo(result)
    expect_true(any(lengths(si$nearby_features) == 0L))
})

test_that("annotateSites: type='proximity' rel_positions contain both positive and negative values", {
    data(comma_example_data)
    result  <- annotateSites(comma_example_data, type = "proximity", window = 100000L)
    si      <- siteInfo(result)
    all_rp  <- unlist(si$rel_positions)
    expect_true(any(all_rp > 0L, na.rm = TRUE))
    expect_true(any(all_rp < 0L, na.rm = TRUE))
})

# ── metagene mode ─────────────────────────────────────────────────────────────

test_that("annotateSites: type='metagene' adds metagene_features and metagene_positions columns", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "metagene")
    si <- siteInfo(result)
    expect_true("metagene_features"  %in% colnames(si))
    expect_true("metagene_positions" %in% colnames(si))
})

test_that("annotateSites: type='metagene' columns are list-typed", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "metagene")
    si <- siteInfo(result)
    expect_true(is(si$metagene_features,  "CharacterList"))
    expect_true(is(si$metagene_positions, "NumericList") ||
                    is(si$metagene_positions, "List"))
})

test_that("annotateSites: type='metagene' all metagene_positions values are in [0, 1]", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "metagene")
    si     <- siteInfo(result)
    all_pos <- unlist(si$metagene_positions)
    expect_true(all(all_pos >= 0 & all_pos <= 1))
})

test_that("annotateSites: type='metagene' non-overlapping sites get length-0 elements", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, type = "metagene")
    si <- siteInfo(result)
    expect_true(any(lengths(si$metagene_features) == 0L))
    expect_true(any(lengths(si$metagene_positions) == 0L))
})

test_that("annotateSites: type='metagene' strand-aware: + strand site at feature start → ~0", {
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 1000L, end = 2000L),
        strand       = "+",
        feature_type = "gene",
        name         = "genePos"
    )
    data(comma_example_data)
    si  <- siteInfo(comma_example_data)
    idx <- which(si$position >= 1000L & si$position <= 1050L)
    if (length(idx) == 0L) skip("No sites near position 1000 for strand test")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, type = "metagene")
    res_si  <- siteInfo(result)
    pos_vals <- unlist(res_si$metagene_positions)
    expect_true(length(pos_vals) > 0L)
    expect_true(all(pos_vals < 0.1))  # near start of + strand feature
})

test_that("annotateSites: type='metagene' - strand site at high coordinate → ~0", {
    # For - strand, 0 = high coordinate (biological TSS)
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 1000L, end = 2000L),
        strand       = "-",
        feature_type = "gene",
        name         = "geneMinus"
    )
    data(comma_example_data)
    si  <- siteInfo(comma_example_data)
    idx <- which(si$position >= 1950L & si$position <= 2000L)
    if (length(idx) == 0L) skip("No sites near position 1950-2000 for - strand test")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, type = "metagene")
    res_si  <- siteInfo(result)
    pos_vals <- unlist(res_si$metagene_positions)
    expect_true(length(pos_vals) > 0L)
    expect_true(all(pos_vals < 0.1))  # near start of - strand feature (high coord)
})

test_that("annotateSites: type='metagene' site overlapping 2 features gets 2 positions", {
    features <- GenomicRanges::GRanges(
        seqnames     = c("chr_sim", "chr_sim"),
        ranges       = IRanges::IRanges(start = c(45000L, 48000L),
                                        end   = c(55000L, 52000L)),
        strand       = c("+", "+"),
        feature_type = c("gene", "regulatory"),
        name         = c("geneA", "regB")
    )
    data(comma_example_data)
    si  <- siteInfo(comma_example_data)
    idx <- which(si$position >= 48000L & si$position <= 52000L)
    if (length(idx) == 0L) skip("No sites in 48000-52000 range for multi-metagene test")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, type = "metagene")
    res_si  <- siteInfo(result)
    expect_equal(length(res_si$metagene_positions[[1L]]), 2L)
    expect_equal(length(res_si$metagene_features[[1L]]),  2L)
})

# ── error handling ────────────────────────────────────────────────────────────

test_that("annotateSites: error on non-commaData input", {
    expect_error(annotateSites(data.frame(x = 1)), "'object' must be a commaData")
})

test_that("annotateSites: error on invalid type", {
    data(comma_example_data)
    expect_error(annotateSites(comma_example_data, type = "badtype"))
})

test_that("annotateSites: error when features NULL and annotation slot is empty", {
    data(comma_example_data)
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

# ── data integrity ────────────────────────────────────────────────────────────

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
