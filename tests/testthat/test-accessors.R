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
    n_6ma <- 10L; n_5mc <- 5L; n_samp <- 3L
    n_total <- n_6ma + n_5mc
    samp_names <- c("ctrl_1", "ctrl_2", "treat_1")

    keys_6ma <- paste0("chr_sim:", seq_len(n_6ma) * 100L, ":+:6mA:GATC")
    keys_5mc <- paste0("chr_sim:", seq_len(n_5mc) * 200L, ":-:5mC:CCWGG")
    all_keys <- c(keys_6ma, keys_5mc)

    methyl <- matrix(runif(n_total * n_samp, 0.1, 0.95),
                     nrow = n_total, dimnames = list(all_keys, samp_names))
    cov    <- matrix(as.integer(runif(n_total * n_samp, 10, 50)),
                     nrow = n_total, dimnames = list(all_keys, samp_names))

    rd <- S4Vectors::DataFrame(
        chrom       = rep("chr_sim", n_total),
        position    = c(seq_len(n_6ma) * 100L, seq_len(n_5mc) * 200L),
        strand      = c(rep("+", n_6ma), rep("-", n_5mc)),
        mod_type    = c(rep("6mA", n_6ma), rep("5mC", n_5mc)),
        motif       = c(rep("GATC", n_6ma), rep("CCWGG", n_5mc)),
        mod_context = c(rep("6mA_GATC", n_6ma), rep("5mC_CCWGG", n_5mc)),
        row.names   = all_keys
    )
    cd <- S4Vectors::DataFrame(
        sample_name = samp_names,
        condition   = c("control", "control", "treatment"),
        replicate   = c(1L, 2L, 1L),
        row.names   = samp_names
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = methyl, coverage = cov),
        rowData = rd, colData = cd
    )
    ann <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = 1L, end = 500L)
    )
    GenomicRanges::mcols(ann)$feature_type <- "gene"
    GenomicRanges::mcols(ann)$name         <- "geneA"

    new("commaData", se,
        genomeInfo = c(chr_sim = 100000L),
        annotation = ann,
        motifSites = GenomicRanges::GRanges())
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
# coverage()
# ─────────────────────────────────────────────────────────────────────────────

test_that("coverage() returns a matrix", {
    obj <- .make_two_modtype()
    cov <- coverage(obj)
    expect_true(is.matrix(cov))
})

test_that("coverage() has correct dimensions", {
    obj <- .make_two_modtype()
    expect_equal(dim(coverage(obj)), dim(methylation(obj)))
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
# subset() method
# ─────────────────────────────────────────────────────────────────────────────

test_that("subset() by mod_type returns only that mod type", {
    obj   <- .make_two_modtype()
    sub   <- subset(obj, mod_type = "6mA")
    types <- unique(rowData(sub)$mod_type)
    expect_equal(types, "6mA")
})

test_that("subset() by mod_type reduces site count", {
    obj  <- .make_two_modtype()
    sub  <- subset(obj, mod_type = "6mA")
    expect_lt(nrow(sub), nrow(obj))
})

test_that("subset() by condition returns only those samples", {
    obj  <- .make_two_modtype()
    sub  <- subset(obj, condition = "control")
    cond <- unique(colData(sub)$condition)
    expect_equal(cond, "control")
    expect_equal(ncol(sub), 2L)
})

test_that("subset() by chrom filters sites by chromosome", {
    obj <- .make_two_modtype()
    sub <- subset(obj, chrom = "chr_sim")
    expect_equal(nrow(sub), nrow(obj))  # all sites are on chr_sim

    sub_none <- subset(obj, chrom = "nonexistent_chr")
    expect_equal(nrow(sub_none), 0L)
})

test_that("subset() returns commaData", {
    obj <- .make_two_modtype()
    sub <- subset(obj, mod_type = "6mA")
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
# subset() by motif
# ─────────────────────────────────────────────────────────────────────────────

test_that("subset() by motif returns only sites with that motif", {
    obj <- .make_two_modtype()
    sub <- subset(obj, motif = "GATC")
    expect_true(all(SummarizedExperiment::rowData(sub)$motif == "GATC"))
})

test_that("subset() by motif reduces site count", {
    obj <- .make_two_modtype()
    sub <- subset(obj, motif = "GATC")
    expect_lt(nrow(sub), nrow(obj))
})

test_that("subset() by motif='GATC' returns only 6mA sites", {
    obj <- .make_two_modtype()
    sub <- subset(obj, motif = "GATC")
    expect_true(all(SummarizedExperiment::rowData(sub)$mod_type == "6mA"))
})

test_that("subset() with both motif and mod_type composes correctly", {
    obj  <- .make_two_modtype()
    sub  <- subset(obj, mod_type = "6mA", motif = "GATC")
    rd   <- SummarizedExperiment::rowData(sub)
    expect_true(all(rd$mod_type == "6mA"))
    expect_true(all(rd$motif   == "GATC"))
    expect_equal(nrow(sub), 10L)  # all 10 6mA sites have motif GATC
})

test_that("subset() by motif returns commaData", {
    obj <- .make_two_modtype()
    sub <- subset(obj, motif = "GATC")
    expect_true(is(sub, "commaData"))
})

test_that("subset() by non-existent motif returns zero-row object", {
    obj <- .make_two_modtype()
    sub <- subset(obj, motif = "TTAA")
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
    expect_true(is(comma_example_data, "commaData"))
    expect_equal(sort(modTypes(comma_example_data)), c("5mC", "6mA"))
    expect_equal(ncol(comma_example_data), 6L)
    expect_equal(nrow(comma_example_data), 588L)
    expect_equal(genome(comma_example_data), c(chr_sim = 100000L))
})

# ─────────────────────────────────────────────────────────────────────────────
# modContexts() accessor
# ─────────────────────────────────────────────────────────────────────────────

test_that("modContexts: returns sorted unique mod_context values", {
    obj <- .make_two_modtype()
    mc <- modContexts(obj)
    expect_equal(mc, c("5mC_CCWGG", "6mA_GATC"))
})

test_that("modContexts: returns character vector", {
    obj <- .make_two_modtype()
    expect_type(modContexts(obj), "character")
})

test_that("modContexts: single context returns length-1 vector", {
    data(comma_example_data)
    sub <- subset(comma_example_data, mod_type = "6mA")
    expect_equal(modContexts(sub), "6mA_GATC")
})

test_that("modContexts: returns 'mod_type' only for NA-motif rows", {
    obj <- .make_two_modtype()
    rd  <- rowData(obj)
    rd$motif      <- NA_character_
    rd$mod_context <- rd$mod_type   # fallback: mod_type without motif suffix
    rowData(obj) <- rd
    mc <- modContexts(obj)
    expect_true(all(mc %in% c("6mA", "5mC")))
    expect_false(any(grepl("NA", mc)))
})

# ─────────────────────────────────────────────────────────────────────────────
# subset(object, mod_context = ...) filtering
# ─────────────────────────────────────────────────────────────────────────────

test_that("subset: mod_context filters to matching rows", {
    obj <- .make_two_modtype()
    sub <- subset(obj, mod_context = "6mA_GATC")
    expect_true(all(rowData(sub)$mod_context == "6mA_GATC"))
    expect_equal(nrow(sub), 10L)
})

test_that("subset: mod_context = '5mC_CCWGG' retains only 5mC rows", {
    obj <- .make_two_modtype()
    sub <- subset(obj, mod_context = "5mC_CCWGG")
    expect_true(all(rowData(sub)$mod_context == "5mC_CCWGG"))
    expect_equal(nrow(sub), 5L)
})

test_that("subset: mod_context with non-existent context returns 0 rows", {
    obj <- .make_two_modtype()
    sub <- subset(obj, mod_context = "4mC_GATC")
    expect_equal(nrow(sub), 0L)
})

test_that("subset: mod_context and mod_type can be combined", {
    obj <- .make_two_modtype()
    sub <- subset(obj, mod_context = "6mA_GATC", mod_type = "6mA")
    expect_equal(nrow(sub), 10L)
    expect_true(all(rowData(sub)$mod_type == "6mA"))
})
