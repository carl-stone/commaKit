library(GenomicRanges)
library(IRanges)

# ── keep = "all" (default unified output) ─────────────────────────────────────

test_that("annotateSites default keep='all' adds all four list columns", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data)
    si <- siteInfo(result)
    expect_true("feature_types" %in% colnames(si))
    expect_true("feature_names" %in% colnames(si))
    expect_true("rel_position"  %in% colnames(si))
    expect_true("frac_position" %in% colnames(si))
})

test_that("annotateSites keep='all' columns are correct types", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data)
    si <- siteInfo(result)
    expect_true(is(si$feature_types, "CharacterList"))
    expect_true(is(si$feature_names, "CharacterList"))
    expect_true(is(si$rel_position,  "IntegerList"))
    expect_true(is(si$frac_position, "List"))
})

test_that("annotateSites keep='all' all columns have length == nrow(siteInfo)", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data)
    si <- siteInfo(result)
    n <- nrow(si)
    expect_equal(length(si$feature_types), n)
    expect_equal(length(si$feature_names), n)
    expect_equal(length(si$rel_position),  n)
    expect_equal(length(si$frac_position), n)
})

test_that("annotateSites keep='all' frac_position is NA for outside sites", {
    data(comma_example_data)
    result  <- annotateSites(comma_example_data)
    si      <- siteInfo(result)
    rp_list <- as.list(si$rel_position)
    fp_list <- as.list(si$frac_position)
    # For each site, indices where rel_position != 0 should have NA frac_position
    for (i in seq_along(rp_list)) {
        rp <- rp_list[[i]]
        fp <- fp_list[[i]]
        if (length(rp) > 0L) {
            outside <- rp != 0L
            if (any(outside)) {
                expect_true(all(is.na(fp[outside])))
            }
        }
    }
})

test_that("annotateSites keep='all' frac_position is in [0,1] for inside sites", {
    data(comma_example_data)
    result  <- annotateSites(comma_example_data)
    si      <- siteInfo(result)
    rp_list <- as.list(si$rel_position)
    fp_list <- as.list(si$frac_position)
    # Indices with rel_position == 0 should have frac in [0,1]
    inside_fracs <- mapply(function(rp, fp) {
        if (length(rp) == 0L) return(numeric(0))
        fp[rp == 0L]
    }, rp_list, fp_list, SIMPLIFY = FALSE)
    all_inside <- unlist(inside_fracs)
    all_inside <- all_inside[!is.na(all_inside)]
    if (length(all_inside) > 0L) {
        expect_true(all(all_inside >= 0 & all_inside <= 1))
    }
})

test_that("annotateSites keep='all' rel_position is 0 for inside sites", {
    # Use a feature covering the whole 100kb chrom
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 1L, end = 100000L),
        strand       = "+",
        feature_type = "gene",
        name         = "bigGene"
    )
    data(comma_example_data)
    result <- annotateSites(comma_example_data, features = features)
    si     <- siteInfo(result)
    # All sites should be inside this feature
    all_rp <- unlist(as.list(si$rel_position))
    expect_true(all(all_rp == 0L))
})

test_that("annotateSites keep='all' rel_position is negative upstream (+ strand)", {
    # Feature at 50000-60000 on + strand; sites < 50000 should be upstream (negative)
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 50000L, end = 60000L),
        strand       = "+",
        feature_type = "gene",
        name         = "genePos"
    )
    data(comma_example_data)
    si  <- siteInfo(comma_example_data)
    # Pick a site just upstream of the feature (within window=50)
    idx <- which(si$position < 50000L & si$position >= 49960L)
    if (length(idx) == 0L) skip("No sites in 49960-49999 for upstream test")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features)
    res_si  <- siteInfo(result)
    rp <- unlist(as.list(res_si$rel_position))
    expect_true(any(rp < 0L))
})

test_that("annotateSites keep='all' rel_position is negative upstream (- strand, high coord)", {
    # Feature at 50000-60000 on - strand; sites > 60000 are upstream on - strand (negative)
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 50000L, end = 60000L),
        strand       = "-",
        feature_type = "gene",
        name         = "geneMinus"
    )
    data(comma_example_data)
    si  <- siteInfo(comma_example_data)
    idx <- which(si$position > 60000L & si$position <= 60040L)
    if (length(idx) == 0L) skip("No sites in 60001-60040 for - strand upstream test")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features)
    res_si  <- siteInfo(result)
    rp <- unlist(as.list(res_si$rel_position))
    expect_true(any(rp < 0L))
})

test_that("annotateSites keep='all' feature outside window gets length-0 lists", {
    # Feature at 1-10; site at position > 60; window=50 → no overlap
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 1L, end = 10L),
        strand       = "+",
        feature_type = "gene",
        name         = "geneY"
    )
    data(comma_example_data)
    si  <- siteInfo(comma_example_data)
    idx <- which(si$position > 60L)
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, window = 50L)
    res_si  <- siteInfo(result)
    expect_equal(length(res_si$feature_types[[1L]]), 0L)
    expect_equal(length(res_si$rel_position[[1L]]),  0L)
    expect_equal(length(res_si$frac_position[[1L]]), 0L)
})

test_that("annotateSites metadata_cols adds parallel CharacterList column", {
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 1L, end = 100000L),
        strand       = "+",
        feature_type = "gene",
        name         = "bigGene",
        extra_col    = "sigmaX"
    )
    data(comma_example_data)
    result <- annotateSites(comma_example_data, features = features,
                             metadata_cols = "extra_col")
    si <- siteInfo(result)
    expect_true("extra_col_values" %in% colnames(si))
    expect_true(is(si$extra_col_values, "CharacterList"))
})

# ── keep = "overlap" ──────────────────────────────────────────────────────────

test_that("annotateSites keep='overlap' returns only feature_types and feature_names", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, keep = "overlap")
    si <- siteInfo(result)
    expect_true("feature_types"  %in% colnames(si))
    expect_true("feature_names"  %in% colnames(si))
    expect_false("rel_position"  %in% colnames(si))
    expect_false("frac_position" %in% colnames(si))
})

test_that("annotateSites keep='overlap' intergenic sites have length-0 elements", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, keep = "overlap")
    si <- siteInfo(result)
    expect_true(any(lengths(si$feature_types) == 0))
})

test_that("annotateSites keep='overlap' genic sites have non-empty elements", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, keep = "overlap")
    si <- siteInfo(result)
    expect_true(any(lengths(si$feature_types) > 0))
})

test_that("annotateSites keep='overlap' site inside feature gets correct type and name", {
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 45000L, end = 55000L),
        strand       = "+",
        feature_type = "gene",
        name         = "geneX"
    )
    data(comma_example_data)
    si  <- siteInfo(comma_example_data)
    idx <- which(si$position >= 45000L & si$position <= 55000L)
    if (length(idx) == 0L) skip("No sites in 45000-55000 for overlap test")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, keep = "overlap")
    res_si  <- siteInfo(result)
    expect_true("gene"  %in% res_si$feature_types[[1L]])
    expect_true("geneX" %in% res_si$feature_names[[1L]])
})

test_that("annotateSites keep='overlap' site overlapping 2 features gets both", {
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
    if (length(idx) == 0L) skip("No sites in 48000-52000 for multi-overlap test")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, keep = "overlap")
    res_si  <- siteInfo(result)
    expect_equal(length(res_si$feature_types[[1L]]), 2L)
})

test_that("annotateSites keep='overlap' assay data unchanged", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, keep = "overlap")
    expect_equal(methylation(result), methylation(comma_example_data))
    expect_equal(coverage(result),    coverage(comma_example_data))
})

# ── keep = "proximity" ────────────────────────────────────────────────────────

test_that("annotateSites keep='proximity' has feature_types, feature_names, rel_position but not frac_position", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, keep = "proximity", window = 10000L)
    si <- siteInfo(result)
    expect_true("feature_types"  %in% colnames(si))
    expect_true("feature_names"  %in% colnames(si))
    expect_true("rel_position"   %in% colnames(si))
    expect_false("frac_position" %in% colnames(si))
})

test_that("annotateSites keep='proximity' rel_position is IntegerList", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, keep = "proximity", window = 10000L)
    si <- siteInfo(result)
    expect_true(is(si$rel_position, "IntegerList"))
})

test_that("annotateSites keep='proximity' sites inside feature have rel_position == 0", {
    features <- GenomicRanges::GRanges(
        seqnames     = "chr_sim",
        ranges       = IRanges::IRanges(start = 45000L, end = 55000L),
        strand       = "+",
        feature_type = "gene",
        name         = "geneX"
    )
    data(comma_example_data)
    si  <- siteInfo(comma_example_data)
    idx <- which(si$position >= 45000L & si$position <= 55000L)
    if (length(idx) == 0L) skip("No sites inside 45000-55000")
    sub_obj <- comma_example_data[idx, ]
    result  <- annotateSites(sub_obj, features = features,
                              keep = "proximity", window = 500L)
    res_si  <- siteInfo(result)
    all_rp  <- unlist(as.list(res_si$rel_position))
    expect_true(all(all_rp == 0L))
})

test_that("annotateSites keep='proximity' rel_position contains positive and negative values", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, keep = "proximity", window = 100000L)
    si     <- siteInfo(result)
    all_rp <- unlist(as.list(si$rel_position))
    expect_true(any(all_rp > 0L, na.rm = TRUE))
    expect_true(any(all_rp < 0L, na.rm = TRUE))
})

test_that("annotateSites keep='proximity' abs(rel_position) within window", {
    window <- 5000L
    data(comma_example_data)
    result <- annotateSites(comma_example_data, keep = "proximity", window = window)
    si     <- siteInfo(result)
    all_rp <- unlist(as.list(si$rel_position))
    expect_true(all(abs(all_rp) <= window))
})

# ── keep = "metagene" ─────────────────────────────────────────────────────────

test_that("annotateSites keep='metagene' has feature_types, feature_names, frac_position but not rel_position", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, keep = "metagene")
    si <- siteInfo(result)
    expect_true("feature_types"  %in% colnames(si))
    expect_true("feature_names"  %in% colnames(si))
    expect_true("frac_position"  %in% colnames(si))
    expect_false("rel_position"  %in% colnames(si))
})

test_that("annotateSites keep='metagene' all frac_position values are in [0,1]", {
    data(comma_example_data)
    result    <- annotateSites(comma_example_data, keep = "metagene")
    si        <- siteInfo(result)
    all_fracs <- unlist(as.list(si$frac_position))
    all_fracs <- all_fracs[!is.na(all_fracs)]
    if (length(all_fracs) > 0L) {
        expect_true(all(all_fracs >= 0 & all_fracs <= 1))
    }
})

test_that("annotateSites keep='metagene' non-overlapping sites get length-0 elements", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data, keep = "metagene")
    si <- siteInfo(result)
    expect_true(any(lengths(si$feature_types) == 0L))
    expect_true(any(lengths(si$frac_position) == 0L))
})

test_that("annotateSites keep='metagene' + strand site near feature start → frac near 0", {
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
    result  <- annotateSites(sub_obj, features = features, keep = "metagene")
    res_si  <- siteInfo(result)
    frac_vals <- unlist(as.list(res_si$frac_position))
    frac_vals <- frac_vals[!is.na(frac_vals)]
    expect_true(length(frac_vals) > 0L)
    expect_true(all(frac_vals < 0.1))  # near start of + strand feature
})

test_that("annotateSites keep='metagene' - strand site near high coord → frac near 0", {
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
    if (length(idx) == 0L) skip("No sites near 1950-2000 for - strand test")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, keep = "metagene")
    res_si  <- siteInfo(result)
    frac_vals <- unlist(as.list(res_si$frac_position))
    frac_vals <- frac_vals[!is.na(frac_vals)]
    expect_true(length(frac_vals) > 0L)
    expect_true(all(frac_vals < 0.1))  # near high coord = TSS on - strand
})

test_that("annotateSites keep='metagene' site overlapping 2 features gets 2 frac_position entries", {
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
    if (length(idx) == 0L) skip("No sites in 48000-52000 range")
    sub_obj <- comma_example_data[idx[1L], ]
    result  <- annotateSites(sub_obj, features = features, keep = "metagene")
    res_si  <- siteInfo(result)
    expect_equal(length(res_si$frac_position[[1L]]), 2L)
    expect_equal(length(res_si$feature_names[[1L]]),  2L)
})

# ── error handling ────────────────────────────────────────────────────────────

test_that("annotateSites errors on non-commaData input", {
    expect_error(annotateSites(data.frame(x = 1)), "'object' must be a commaData")
})

test_that("annotateSites errors on invalid keep value", {
    data(comma_example_data)
    expect_error(annotateSites(comma_example_data, keep = "badkeep"))
})

test_that("annotateSites errors when features NULL and annotation slot is empty", {
    data(comma_example_data)
    obj_no_ann <- new("commaData",
        as(comma_example_data, "SummarizedExperiment"),
        genomeInfo = comma_example_data@genomeInfo,
        annotation = GenomicRanges::GRanges(),
        motifSites = comma_example_data@motifSites
    )
    expect_error(annotateSites(obj_no_ann, features = NULL), "No features available")
})

test_that("annotateSites errors on missing feature_col", {
    data(comma_example_data)
    features <- annotation(comma_example_data)
    expect_error(
        annotateSites(comma_example_data, features = features, feature_col = "nonexistent"),
        "not found in mcols"
    )
})

test_that("annotateSites errors on metadata_cols not in mcols(features)", {
    data(comma_example_data)
    features <- annotation(comma_example_data)
    expect_error(
        annotateSites(comma_example_data, features = features,
                      metadata_cols = "nonexistent_col"),
        "metadata_cols not found"
    )
})

# ── data integrity ────────────────────────────────────────────────────────────

test_that("annotateSites rowData has same number of rows after annotation", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data)
    expect_equal(nrow(siteInfo(result)), nrow(siteInfo(comma_example_data)))
})

test_that("annotateSites assay data unchanged after annotation", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data)
    expect_equal(methylation(result), methylation(comma_example_data))
    expect_equal(coverage(result),    coverage(comma_example_data))
})

test_that("annotateSites returns commaData with same dimensions", {
    data(comma_example_data)
    result <- annotateSites(comma_example_data)
    expect_s4_class(result, "commaData")
    expect_equal(dim(result), dim(comma_example_data))
})
