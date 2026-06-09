## Tests for commaData accessor functions and subsetting methods.
## Uses comma_example_data where available, falling back to inline construction.

library(testthat)
library(SummarizedExperiment)
library(S4Vectors)
library(GenomicRanges)

# ─────────────────────────────────────────────────────────────────────────────
# Helper: build a two-mod-type commaData for testing
# ─────────────────────────────────────────────────────────────────────────────

.make_two_modtype <- function() {
    obj <- .make_two_modtype_fixture(
        n_6ma = 10L,
        n_5mc = 5L,
        sample_names = c("ctrl_1", "ctrl_2", "treat_1"),
        conditions = c("control", "control", "treatment"),
        replicate = c(1L, 2L, 1L),
        seed = 11L
    )

    # Store annotation in metadata
    ann <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = 1L, end = 500L)
    )
    GenomicRanges::mcols(ann)$feature_type <- "gene"
    GenomicRanges::mcols(ann)$name         <- "geneA"
    S4Vectors::metadata(obj)$annotation <- ann
    S4Vectors::metadata(obj)$motifSites <- GenomicRanges::GRanges()

    obj
}

# ─────────────────────────────────────────────────────────────────────────────
# methylation()
# ─────────────────────────────────────────────────────────────────────────────

test_that("methylation() returns a numeric matrix", {
    obj <- .make_two_modtype()
    m <- methylation(obj)
    expect_true(is.matrix(m))
    expect_true(is.numeric(m))
})

test_that("methylation() has correct dimensions", {
    obj <- .make_two_modtype()
    m <- methylation(obj)
    expect_equal(nrow(m), nrow(obj))
    expect_equal(ncol(m), ncol(obj))
})

test_that("methylation() values are in [0, 1] range (ignoring NA)", {
    obj <- .make_two_modtype()
    m <- methylation(obj)
    valid <- m[!is.na(m)]
    expect_true(all(valid >= 0 & valid <= 1))
})

# ─────────────────────────────────────────────────────────────────────────────
# siteCoverage()
# ─────────────────────────────────────────────────────────────────────────────

test_that("siteCoverage() returns a numeric matrix", {
    obj <- .make_two_modtype()
    cov <- siteCoverage(obj)
    expect_true(is.matrix(cov))
    expect_true(is.numeric(cov))
    # Coverage values should be non-negative integers
    valid <- cov[!is.na(cov)]
    expect_true(length(valid) > 0)
    expect_true(all(valid >= 0))
})

test_that("siteCoverage() has correct dimensions", {
    obj <- .make_two_modtype()
    expect_equal(dim(siteCoverage(obj)), dim(methylation(obj)))
})

# ─────────────────────────────────────────────────────────────────────────────
# modCounts(), canonicalCounts(), assayProvenance()
# ─────────────────────────────────────────────────────────────────────────────

test_that("modCounts() and canonicalCounts() return raw count assays", {
    obj <- .make_two_modtype()
    mod <- modCounts(obj)
    canonical <- canonicalCounts(obj)

    expect_true(is.matrix(mod))
    expect_true(is.matrix(canonical))
    expect_equal(dim(mod), dim(methylation(obj)))
    expect_equal(dim(canonical), dim(methylation(obj)))
    expect_true(all(mod + canonical <= siteCoverage(obj), na.rm = TRUE))
})

test_that("raw count accessors error clearly for legacy objects", {
    obj <- .make_two_modtype()
    SummarizedExperiment::assays(obj) <-
        SummarizedExperiment::assays(obj)[c("methylation", "coverage")]

    expect_error(modCounts(obj), "mod_counts")
    expect_error(canonicalCounts(obj), "canonical_counts")
    expect_no_error(validObject(obj))
})

test_that("assayProvenance() returns recorded layer metadata or an empty list", {
    obj <- .make_two_modtype()
    provenance <- assayProvenance(obj)

    expect_true(is.list(provenance))
    expect_true(all(c("methylation", "coverage", "mod_counts",
                      "canonical_counts") %in% names(provenance)))
    expect_equal(provenance$mod_counts$type, "reconstructed_counts")

    S4Vectors::metadata(obj)$assay_provenance <- NULL
    expect_identical(assayProvenance(obj), list())
})

# ─────────────────────────────────────────────────────────────────────────────
# sampleInfo()
# ─────────────────────────────────────────────────────────────────────────────

test_that("sampleInfo() returns a data.frame", {
    obj <- .make_two_modtype()
    si <- sampleInfo(obj)
    expect_true(is.data.frame(si))
})

test_that("sampleInfo() has one row per sample", {
    obj <- .make_two_modtype()
    expect_equal(nrow(sampleInfo(obj)), ncol(obj))
})

test_that("sampleInfo() contains required columns", {
    obj <- .make_two_modtype()
    si <- sampleInfo(obj)
    expect_true(all(c("sample_name", "condition", "replicate") %in% colnames(si)))
})

# ─────────────────────────────────────────────────────────────────────────────
# siteInfo()
# ─────────────────────────────────────────────────────────────────────────────

test_that("siteInfo() returns a DataFrame", {
    obj <- .make_two_modtype()
    expect_true(is(siteInfo(obj), "DataFrame"))
})

test_that("siteInfo() has one row per site", {
    obj <- .make_two_modtype()
    expect_equal(nrow(siteInfo(obj)), nrow(obj))
})

test_that("siteInfo() contains required columns", {
    obj <- .make_two_modtype()
    si <- siteInfo(obj)
    expect_true(all(c("chrom", "position", "strand", "mod_type") %in% colnames(si)))
})

test_that("siteInfo() site_key always includes chromosome for single-chromosome objects", {
    obj <- .make_two_modtype()
    si <- siteInfo(obj)

    expect_equal(
        si$site_key[1],
        paste(si$chrom[1], si$position[1], si$strand[1],
              as.character(si$mod_type[1]), si$motif[1], sep = ":")
    )
    expect_equal(length(strsplit(si$site_key[1], ":", fixed = TRUE)[[1]]), 5L)
})

test_that("siteInfo() site_key uses same field contract for multi-chromosome objects", {
    obj <- .make_two_modtype()
    rr <- rowRanges(obj)
    GenomeInfoDb::seqlevels(rr) <- c("chr_sim", "chr_alt")
    GenomeInfoDb::seqnames(rr)[seq_len(3L)] <- "chr_alt"
    GenomeInfoDb::seqlengths(rr) <- c(chr_sim = 100000L, chr_alt = 50000L)
    rowRanges(obj) <- rr

    si <- siteInfo(obj)

    expect_true(any(si$chrom == "chr_alt"))
    expect_true(all(vapply(
        strsplit(si$site_key, ":", fixed = TRUE),
        length,
        integer(1L)
    ) == 5L))
    expect_equal(
        si$site_key,
        paste(si$chrom, si$position, si$strand,
              as.character(si$mod_type), si$motif, sep = ":")
    )
})

# ─────────────────────────────────────────────────────────────────────────────
# modTypes()
# ─────────────────────────────────────────────────────────────────────────────

test_that("modTypes() returns a character vector", {
    obj <- .make_two_modtype()
    expect_type(modTypes(obj), "character")
})

test_that("modTypes() returns both 6mA and 5mC for two-mod-type object", {
    obj <- .make_two_modtype()
    mt <- modTypes(obj)
    expect_true("6mA" %in% mt)
    expect_true("5mC" %in% mt)
})

test_that("modTypes() returns only unique values", {
    obj <- .make_two_modtype()
    mt <- modTypes(obj)
    expect_equal(length(mt), length(unique(mt)))
})

# ─────────────────────────────────────────────────────────────────────────────
# genome()
# ─────────────────────────────────────────────────────────────────────────────

test_that("genome() returns a named integer vector", {
    obj <- .make_two_modtype()
    g <- genome(obj)
    expect_true(is.integer(g))
    expect_false(is.null(names(g)))
})

test_that("genome() returns correct chromosome sizes", {
    obj <- .make_two_modtype()
    expect_equal(genome(obj), c(chr_sim = 100000L))
})

test_that("genomeSizes() returns chromosome sizes with genome() compatibility", {
    obj <- .make_two_modtype()
    expect_equal(genomeSizes(obj), c(chr_sim = 100000L))
    expect_equal(genome(obj), genomeSizes(obj))
})

test_that("genome() returns NULL when no Seqinfo", {
    obj <- .make_two_modtype()
    rr <- rowRanges(obj)
    GenomeInfoDb::seqlengths(rr) <- NA_integer_
    rowRanges(obj) <- rr
    expect_null(genome(obj))
})

# ─────────────────────────────────────────────────────────────────────────────
# annotation()
# ─────────────────────────────────────────────────────────────────────────────

test_that("annotation() returns a GRanges", {
    obj <- .make_two_modtype()
    expect_true(is(annotation(obj), "GRanges"))
})

test_that("annotation() returns correct number of features", {
    obj <- .make_two_modtype()
    expect_equal(length(annotation(obj)), 1L)
})

# ─────────────────────────────────────────────────────────────────────────────
# motifSites()
# ─────────────────────────────────────────────────────────────────────────────

test_that("motifSites() returns a GRanges", {
    obj <- .make_two_modtype()
    expect_true(is(motifSites(obj), "GRanges"))
})

test_that("motifSites() is empty when no motif was specified", {
    obj <- .make_two_modtype()
    expect_equal(length(motifSites(obj)), 0L)
})

# ─────────────────────────────────────────────────────────────────────────────
# [ subsetting
# ─────────────────────────────────────────────────────────────────────────────

test_that("[ subsetting by site index returns correct number of sites", {
    obj <- .make_two_modtype()
    sub <- obj[1:5, ]
    expect_equal(nrow(sub), 5L)
    expect_equal(ncol(sub), ncol(obj))
})

test_that("[ subsetting preserves commaData class", {
    obj <- .make_two_modtype()
    sub <- obj[1:3, ]
    expect_true(is(sub, "commaData"))
})

test_that("[ subsetting keeps custom slots intact", {
    obj <- .make_two_modtype()
    sub <- obj[1:3, ]
    expect_equal(genome(sub), genome(obj))
    expect_equal(length(annotation(sub)), length(annotation(obj)))
})

test_that("[ subsetting by sample index returns correct number of samples", {
    obj <- .make_two_modtype()
    sub <- obj[, 1:2]
    expect_equal(ncol(sub), 2L)
    expect_equal(nrow(sub), nrow(obj))
})

test_that("[ subsetting by logical vector works", {
    obj <- .make_two_modtype()
    keep <- rowData(obj)$mod_type == "6mA"
    sub  <- obj[keep, ]
    expect_equal(nrow(sub), sum(keep))
})

# ─────────────────────────────────────────────────────────────────────────────
# filterSites() method
# ─────────────────────────────────────────────────────────────────────────────

test_that("filterSites() by mod_type returns only that mod type", {
    obj   <- .make_two_modtype()
    sub   <- filterSites(obj, mod_type = "6mA")
    types <- as.character(unique(rowData(sub)$mod_type))
    expect_equal(types, "6mA")
})

test_that("filterSites() by mod_type reduces site count", {
    obj  <- .make_two_modtype()
    sub  <- filterSites(obj, mod_type = "6mA")
    expect_lt(nrow(sub), nrow(obj))
})

test_that("filterSites() by condition returns only those samples", {
    obj  <- .make_two_modtype()
    sub  <- filterSites(obj, condition = "control")
    cond <- unique(colData(sub)$condition)
    expect_equal(cond, "control")
    expect_equal(ncol(sub), 2L)
})

test_that("filterSites() by chrom filters sites by chromosome", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, chrom = "chr_sim")
    expect_equal(nrow(sub), nrow(obj))  # all sites are on chr_sim

    sub_none <- filterSites(obj, chrom = "nonexistent_chr")
    expect_equal(nrow(sub_none), 0L)
})

test_that("filterSites() returns commaData", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, mod_type = "6mA")
    expect_true(is(sub, "commaData"))
})

# ─────────────────────────────────────────────────────────────────────────────
# comma_example_data integration tests (skipped if data not yet generated)
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# motifs() accessor
# ─────────────────────────────────────────────────────────────────────────────

test_that("motifs() returns a character vector", {
    obj <- .make_two_modtype()
    expect_type(motifs(obj), "character")
})

test_that("motifs() returns sorted unique non-NA values", {
    obj <- .make_two_modtype()
    m <- motifs(obj)
    expect_equal(m, c("CCWGG", "GATC"))
})

test_that("motifs() excludes NA values", {
    obj <- .make_two_modtype()
    rd  <- SummarizedExperiment::rowData(obj)
    rd$motif[1L] <- NA_character_
    SummarizedExperiment::rowData(obj) <- rd
    m <- motifs(obj)
    expect_false(any(is.na(m)))
})

test_that("motifs() returns empty character vector when all motifs are NA", {
    obj <- .make_two_modtype()
    rd  <- SummarizedExperiment::rowData(obj)
    rd$motif <- NA_character_
    SummarizedExperiment::rowData(obj) <- rd
    expect_equal(length(motifs(obj)), 0L)
})

# ─────────────────────────────────────────────────────────────────────────────
# filterSites() by motif
# ─────────────────────────────────────────────────────────────────────────────

test_that("filterSites() by motif returns only sites with that motif", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, motif = "GATC")
    expect_true(all(SummarizedExperiment::rowData(sub)$motif == "GATC"))
})

test_that("filterSites() by motif reduces site count", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, motif = "GATC")
    expect_lt(nrow(sub), nrow(obj))
})

test_that("filterSites() by motif='GATC' returns only 6mA sites", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, motif = "GATC")
    expect_true(all(SummarizedExperiment::rowData(sub)$mod_type == "6mA"))
})

test_that("filterSites() with both motif and mod_type composes correctly", {
    obj  <- .make_two_modtype()
    sub  <- filterSites(obj, mod_type = "6mA", motif = "GATC")
    rd   <- SummarizedExperiment::rowData(sub)
    expect_true(all(rd$mod_type == "6mA"))
    expect_true(all(rd$motif   == "GATC"))
    expect_equal(nrow(sub), 10L)  # all 10 6mA sites have motif GATC
})

test_that("filterSites() by motif returns commaData", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, motif = "GATC")
    expect_true(is(sub, "commaData"))
})

test_that("filterSites() by non-existent motif returns zero-row object", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, motif = "TTAA")
    expect_equal(nrow(sub), 0L)
})

# ─────────────────────────────────────────────────────────────────────────────
# comma_example_data integration tests (skipped if data not yet generated)
# ─────────────────────────────────────────────────────────────────────────────

test_that("comma_example_data loads and accessors work correctly", {
    skip_if_not(exists("comma_example_data") ||
                tryCatch({ data(comma_example_data); TRUE }, error = function(e) FALSE),
                message = "comma_example_data not yet generated — run data-raw/create_example_data.R")

    data(comma_example_data)
    expected_samples <- c("ctrl_1", "ctrl_2", "ctrl_3",
                          "treat_1", "treat_2", "treat_3")

    expect_true(is(comma_example_data, "commaData"))
    expect_equal(sort(modTypes(comma_example_data)), c("5mC", "6mA"))
    expect_equal(ncol(comma_example_data), 6L)
    expect_equal(nrow(comma_example_data), 588L)
    expect_equal(genome(comma_example_data), c(chr_sim = 100000L))
    expect_equal(colnames(comma_example_data), expected_samples)

    samples <- sampleInfo(comma_example_data)
    expect_equal(as.character(samples$sample_name), expected_samples)
    expect_equal(as.character(samples$condition),
                 c(rep("control", 3L), rep("treatment", 3L)))

    sites <- siteInfo(comma_example_data)
    mod_type <- as.character(sites$mod_type)
    motif <- as.character(sites$motif)
    expect_equal(sum(mod_type == "6mA" & motif == "GATC"), 393L)
    expect_equal(sum(mod_type == "5mC" & motif == "CCWGG"), 195L)
    expect_equal(range(siteCoverage(comma_example_data), na.rm = TRUE),
                 c(10, 150))
})

# ─────────────────────────────────────────────────────────────────────────────
# modContexts() accessor
# ─────────────────────────────────────────────────────────────────────────────

test_that("modContexts: returns sorted unique mod_context values", {
    obj <- .make_two_modtype()
    mc <- modContexts(obj)
    expect_equal(mc, c("5mC_CCWGG", "6mA_GATC"))
})

test_that("modContexts: works even without mod_context in mcols", {
    obj <- .make_two_modtype()
    # mod_context should not be in mcols anymore
    expect_false("mod_context" %in% colnames(GenomicRanges::mcols(rowRanges(obj))))
    # But modContexts() should still work
    mc <- modContexts(obj)
    expect_equal(mc, c("5mC_CCWGG", "6mA_GATC"))
})

test_that("modContexts: returns character vector", {
    obj <- .make_two_modtype()
    expect_type(modContexts(obj), "character")
})

test_that("modContexts: single context returns length-1 vector", {
    data(comma_example_data)
    sub <- filterSites(comma_example_data, mod_type = "6mA")
    expect_equal(modContexts(sub), "6mA_GATC")
})

test_that("modContexts: returns 'mod_type' only for NA-motif rows", {
    obj <- .make_two_modtype()
    rd  <- rowData(obj)
    rd$motif      <- NA_character_
    rowData(obj) <- rd
    mc <- modContexts(obj)
    expect_true(all(mc %in% c("6mA", "5mC")))
    expect_false(any(grepl("NA", mc)))
})

# ─────────────────────────────────────────────────────────────────────────────
# filterSites(object, mod_context = ...) filtering
# ─────────────────────────────────────────────────────────────────────────────

test_that("filterSites: mod_context filters to matching rows", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, mod_context = "6mA_GATC")
    # Use siteInfo() to check mod_context (computed on demand)
    si <- siteInfo(sub)
    expect_true(all(si$mod_context == "6mA_GATC"))
    expect_equal(nrow(sub), 10L)
})

test_that("filterSites: mod_context = '5mC_CCWGG' retains only 5mC rows", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, mod_context = "5mC_CCWGG")
    si <- siteInfo(sub)
    expect_true(all(si$mod_context == "5mC_CCWGG"))
    expect_equal(nrow(sub), 5L)
})

test_that("filterSites: mod_context with non-existent context returns 0 rows", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, mod_context = "4mC_GATC")
    expect_equal(nrow(sub), 0L)
})

test_that("filterSites: mod_context and mod_type can be combined", {
    obj <- .make_two_modtype()
    sub <- filterSites(obj, mod_context = "6mA_GATC", mod_type = "6mA")
    expect_equal(nrow(sub), 10L)
    expect_true(all(rowData(sub)$mod_type == "6mA"))
})

# ─────────────────────────────────────────────────────────────────────────────
# caller() and minCoverage()
# ─────────────────────────────────────────────────────────────────────────────

test_that("caller() returns NA for objects without stored caller", {
    obj <- .make_two_modtype()
    expect_equal(caller(obj), NA_character_)
})

test_that("minCoverage() returns NA for objects without stored min_coverage", {
    obj <- .make_two_modtype()
    expect_equal(minCoverage(obj), NA_integer_)
})

test_that("caller() returns stored caller from metadata", {
    obj <- .make_two_modtype()
    S4Vectors::metadata(obj)$caller <- "modkit"
    expect_equal(caller(obj), "modkit")
})

test_that("minCoverage() returns stored min_coverage from metadata", {
    obj <- .make_two_modtype()
    S4Vectors::metadata(obj)$min_coverage <- 10L
    expect_equal(minCoverage(obj), 10L)
})

test_that("example data has caller and min_coverage", {
    data(comma_example_data)
    expect_equal(caller(comma_example_data), "modkit")
    expect_equal(minCoverage(comma_example_data), 5L)
})

test_that("show() displays caller and min_coverage", {
    data(comma_example_data)
    out <- capture.output(show(comma_example_data))
    expect_true(any(grepl("caller:", out, fixed = TRUE)))
    expect_true(any(grepl("min_coverage:", out, fixed = TRUE)))
})

test_that("coverage() compatibility wrapper warns and returns siteCoverage", {
    obj <- .make_two_modtype()
    expect_warning(
        cov <- coverage(obj),
        regexp = "deprecated"
    )
    expect_equal(cov, siteCoverage(obj))
})

test_that("coverage() compatibility wrapper rejects ignored IRanges arguments", {
    obj <- .make_two_modtype()
    expect_error(coverage(obj, shift = 1L), regexp = "siteCoverage")
    expect_error(coverage(obj, width = 100L), regexp = "siteCoverage")
    expect_error(coverage(obj, weight = 2L), regexp = "siteCoverage")
})

test_that("subset.commaData compatibility method warns and returns filterSites", {
    obj <- .make_two_modtype()
    expect_warning(
        sub <- subset(obj, mod_type = "6mA"),
        regexp = "deprecated"
    )
    expect_equal(sub, filterSites(obj, mod_type = "6mA"))
})
