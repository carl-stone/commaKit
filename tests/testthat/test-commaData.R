## Tests for the commaData S4 class, constructor, and show() method
##
## These tests use:
##  - system.file("extdata", "example_modkit.bed", package = "commaKit") for file-based tests
##  - Direct object construction for class/validity tests (avoids file I/O overhead)

library(testthat)
library(SummarizedExperiment)
library(S4Vectors)
library(GenomicRanges)

# ── Helper: build a minimal valid commaData without file I/O ──────────────────

.make_minimal_commaData <- function(n_sites = 5L, n_samples = 2L) {
    methyl    <- matrix(runif(n_sites * n_samples, 0.1, 0.9),
                        nrow = n_sites, ncol = n_samples,
                        dimnames = list(NULL, paste0("s", seq_len(n_samples))))
    cov       <- matrix(as.integer(runif(n_sites * n_samples, 10, 50)),
                        nrow = n_sites, ncol = n_samples,
                        dimnames = list(NULL, paste0("s", seq_len(n_samples))))
    sample_info <- data.frame(
        sample_name = paste0("s", seq_len(n_samples)),
        condition   = rep(c("control", "treatment"), length.out = n_samples),
        replicate   = seq_len(n_samples),
        stringsAsFactors = FALSE
    )
    .make_commaData_fixture(
        methyl,
        cov,
        sample_info,
        positions = seq_len(n_sites) * 100L
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# Class instantiation and validity
# ─────────────────────────────────────────────────────────────────────────────

test_that("minimal commaData object is valid", {
    obj <- .make_minimal_commaData()
    expect_true(is(obj, "commaData"))
    expect_true(is(obj, "RangedSummarizedExperiment"))
    expect_no_error(validObject(obj))
})

test_that("assay matrices have no rownames", {
    obj <- .make_minimal_commaData()
    expect_null(rownames(methylation(obj)))
    expect_null(rownames(siteCoverage(obj)))
})

test_that("rowRanges has no names", {
    obj <- .make_minimal_commaData()
    expect_null(names(rowRanges(obj)))
})

test_that("siteInfo() has no rownames and includes site_key", {
    obj <- .make_minimal_commaData()
    si <- siteInfo(obj)
    expect_null(rownames(si))
    expect_true("site_key" %in% colnames(si))
    expect_equal(si$site_key[1], "chr_sim:100:+:6mA:GATC")
})

test_that("validity rejects missing rowData columns", {
    obj <- .make_minimal_commaData()
    # Remove required column
    rowData(obj)$mod_type <- NULL
    expect_error(validObject(obj), regexp = "mod_type")
})

test_that("validity rejects unrecognized mod_type values", {
    obj <- .make_minimal_commaData()
    # Assigning a value not in factor levels converts the column to character
    mcols(rowRanges(obj))$mod_type <- factor(rep("7mX", nrow(obj)), levels = c("7mX"))
    expect_error(validObject(obj), regexp = "7mX")
})

test_that("validity rejects non-factor mod_type column", {
    obj <- .make_minimal_commaData()
    # Replace factor with character vector (simulates old-style object)
    mcols(rowRanges(obj))$mod_type <- rep("6mA", nrow(obj))
    expect_error(validObject(obj), regexp = "must be a factor")
})

test_that("validity rejects rowRanges with width != 1", {
    obj <- .make_minimal_commaData()
    rr <- rowRanges(obj)
    # Widen the first range to 10bp
    IRanges::ranges(rr)[1, ] <- IRanges::IRanges(start = start(rr)[1], width = 10L)
    rowRanges(obj) <- rr
    expect_error(validObject(obj), regexp = "1-bp ranges")
})

test_that("validity rejects missing colData columns", {
    obj <- .make_minimal_commaData()
    colData(obj)$condition <- NULL
    expect_error(validObject(obj), regexp = "condition")
})


# ─────────────────────────────────────────────────────────────────────────────
# show() method
# ─────────────────────────────────────────────────────────────────────────────

test_that("show() produces output without error", {
    obj <- .make_minimal_commaData()
    expect_output(show(obj), regexp = "commaData")
    expect_output(show(obj), regexp = "sites")
    expect_output(show(obj), regexp = "samples")
})

test_that("show() lists modification types", {
    obj <- .make_minimal_commaData()
    expect_output(show(obj), regexp = "6mA")
})

test_that("show() prints motifs line with GATC", {
    obj <- .make_minimal_commaData()
    expect_output(show(obj), regexp = "motifs")
    expect_output(show(obj), regexp = "GATC")
})

test_that("show() prints 'not available' when all motifs are NA", {
    obj <- .make_minimal_commaData()
    rd  <- rowData(obj)
    rd$motif <- NA_character_
    rowData(obj) <- rd
    expect_output(show(obj), regexp = "not available")
})

test_that("validity passes when motif column contains NA values", {
    obj <- .make_minimal_commaData()
    rd  <- rowData(obj)
    rd$motif[1L] <- NA_character_
    rowData(obj) <- rd
    expect_no_error(validObject(obj))
})

test_that("validity passes when all motif values are NA", {
    obj <- .make_minimal_commaData()
    rd  <- rowData(obj)
    rd$motif <- NA_character_
    rowData(obj) <- rd
    expect_no_error(validObject(obj))
})

# ─────────────────────────────────────────────────────────────────────────────
# Constructor: argument validation
# ─────────────────────────────────────────────────────────────────────────────

test_that("commaData() errors on colData missing required columns", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    expect_error(
        commaData(
            files   = c(s1 = bed_file),
            colData = data.frame(sample_name = "s1", replicate = 1L),
            genome  = c(chr_sim = 100000L)
        ),
        regexp = "condition"
    )
})

test_that("commaData() errors on non-named files vector", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    expect_error(
        commaData(
            files   = bed_file,    # no names
            colData = data.frame(sample_name = "s1", condition = "ctrl", replicate = 1L),
            genome  = c(chr_sim = 100000L)
        ),
        regexp = "named character vector"
    )
})

test_that("commaData() errors on mismatched sample names", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    expect_error(
        commaData(
            files   = c(wrong_name = bed_file),
            colData = data.frame(sample_name = "s1", condition = "ctrl", replicate = 1L),
            genome  = c(chr_sim = 100000L)
        ),
        regexp = "wrong_name"
    )
})

test_that("commaData() errors on invalid caller", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    expect_error(
        commaData(
            files   = c(s1 = bed_file),
            colData = data.frame(sample_name = "s1", condition = "ctrl", replicate = 1L),
            genome  = c(chr_sim = 100000L),
            caller  = "unknown_caller"
        )
    )
})

# ─────────────────────────────────────────────────────────────────────────────
# Constructor: successful construction from extdata example file
# ─────────────────────────────────────────────────────────────────────────────

test_that("commaData() constructs valid object from example modkit BED", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    genome_info <- c(chr_sim = 100000L)
    cd <- commaData(
        files   = c(s1 = bed_file),
        colData = data.frame(
            sample_name = "s1",
            condition   = "control",
            replicate   = 1L,
            stringsAsFactors = FALSE
        ),
        genome  = genome_info,
        caller  = "modkit"
    )

    expect_true(is(cd, "commaData"))
    expect_no_error(validObject(cd))
    expect_true(nrow(cd) > 0)
    expect_equal(ncol(cd), 1L)
    expect_true("methylation" %in% assayNames(cd))
    expect_true("coverage" %in% assayNames(cd))
    expect_true("mod_counts" %in% assayNames(cd))
    expect_true("canonical_counts" %in% assayNames(cd))
    expect_true(all(!is.na(modCounts(cd))))
    expect_true(all(!is.na(canonicalCounts(cd))))
    expect_true(all(modCounts(cd) + canonicalCounts(cd) <= siteCoverage(cd)))
    expect_true(all(GenomeInfoDb::isCircular(GenomeInfoDb::seqinfo(cd))))
    expect_equal(
        GenomeInfoDb::seqinfo(SummarizedExperiment::rowRanges(cd)),
        commaKit:::.makeSeqinfo(genome_info)
    )
})

test_that("commaData() applies min_coverage filter correctly", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    # All sites in example file have coverage >= 18, so min_coverage = 50
    # should set all beta values to NA
    cd_high <- commaData(
        files        = c(s1 = bed_file),
        colData      = data.frame(sample_name = "s1", condition = "ctrl",
                                  replicate = 1L, stringsAsFactors = FALSE),
        genome       = c(chr_sim = 100000L),
        min_coverage = 50L
    )
    # Most sites should have NA beta (coverage < 50 in example data)
    m <- methylation(cd_high)
    expect_true(sum(is.na(m)) > 0)
})

test_that("commaData() accepts a tibble as colData without warning", {
    skip_if_not_installed("tibble")
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    tbl_cd <- tibble::tibble(
        sample_name = "s1",
        condition   = "control",
        replicate   = 1L
    )
    expect_no_warning(
        commaData(
            files   = c(s1 = bed_file),
            colData = tbl_cd,
            genome  = c(chr_sim = 100000L),
            caller  = "modkit"
        )
    )
})

test_that("commaData() mod_type filter reduces sites", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    cd_all  <- commaData(
        files   = c(s1 = bed_file),
        colData = data.frame(sample_name = "s1", condition = "ctrl",
                              replicate = 1L, stringsAsFactors = FALSE),
        genome  = c(chr_sim = 100000L)
    )
    cd_6ma  <- commaData(
        files    = c(s1 = bed_file),
        colData  = data.frame(sample_name = "s1", condition = "ctrl",
                               replicate = 1L, stringsAsFactors = FALSE),
        genome   = c(chr_sim = 100000L),
        mod_type = "6mA"
    )
    # Only 6mA sites retained
    expect_true(nrow(cd_6ma) <= nrow(cd_all))
    expect_equal(as.character(unique(rowData(cd_6ma)$mod_type)), "6mA")
})

# ─────────────────────────────────────────────────────────────────────────────
# expected_mod_contexts constructor filter
# ─────────────────────────────────────────────────────────────────────────────

test_that("commaData: expected_mod_contexts filters to specified contexts", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    cd_all <- commaData(
        files   = c(s1 = bed_file),
        colData = data.frame(sample_name = "s1", condition = "ctrl",
                             replicate = 1L, stringsAsFactors = FALSE),
        genome  = c(chr_sim = 100000L)
    )
    expect_message(
        cd_6mA <- commaData(
            files                 = c(s1 = bed_file),
            colData               = data.frame(sample_name = "s1", condition = "ctrl",
                                               replicate = 1L, stringsAsFactors = FALSE),
            genome                = c(chr_sim = 100000L),
            expected_mod_contexts = list("6mA" = "GATC")
        ),
        regexp = "dropping"
    )
    # Only 6mA_GATC sites remain (mod_context is computed on demand)
    si <- siteInfo(cd_6mA)
    expect_true(all(si$mod_context == "6mA_GATC"))
    expect_true(nrow(cd_6mA) < nrow(cd_all))
})

test_that("commaData: expected_mod_contexts accepts multiple mod types", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    cd <- commaData(
        files                 = c(s1 = bed_file),
        colData               = data.frame(sample_name = "s1", condition = "ctrl",
                                           replicate = 1L, stringsAsFactors = FALSE),
        genome                = c(chr_sim = 100000L),
        expected_mod_contexts = list("6mA" = "GATC", "5mC" = "CCWGG")
    )
    si <- siteInfo(cd)
    expect_true(all(si$mod_context %in% c("6mA_GATC", "5mC_CCWGG")))
})

test_that("commaData: expected_mod_contexts stops if no sites remain", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    expect_error(
        commaData(
            files                 = c(s1 = bed_file),
            colData               = data.frame(sample_name = "s1", condition = "ctrl",
                                               replicate = 1L, stringsAsFactors = FALSE),
            genome                = c(chr_sim = 100000L),
            expected_mod_contexts = list("6mA" = "TTAA")  # no TTAA motif in data
        ),
        regexp = "No sites remain"
    )
})

test_that("commaData: expected_mod_contexts errors with unrecognized mod_type", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
    skip_if(bed_file == "", message = "extdata not available")

    expect_error(
        commaData(
            files                 = c(s1 = bed_file),
            colData               = data.frame(sample_name = "s1", condition = "ctrl",
                                               replicate = 1L, stringsAsFactors = FALSE),
            genome                = c(chr_sim = 100000L),
            expected_mod_contexts = list("7mX" = "GATC")
        ),
        regexp = "Unrecognized"
    )
})

# ─────────────────────────────────────────────────────────────────────────────
# mod_context is computed on demand (no longer stored in rowData)
# ─────────────────────────────────────────────────────────────────────────────

test_that("mod_context is not stored in rowData but available via siteInfo", {
    obj <- .make_minimal_commaData()
    # mod_context should not be a column in rowData/mcols
    expect_false("mod_context" %in% colnames(rowData(obj)))
    # But it should be available via siteInfo()
    si <- siteInfo(obj)
    expect_true("mod_context" %in% colnames(si))
    expect_true(all(si$mod_context == "6mA_GATC"))
})

test_that("modContexts() returns correct values without mod_context column", {
    obj <- .make_minimal_commaData()
    expect_equal(modContexts(obj), "6mA_GATC")
})

test_that("mod_type is a factor with valid levels in example data", {
    data(comma_example_data)
    mc <- mcols(rowRanges(comma_example_data))
    expect_true(is.factor(mc$mod_type))
    expect_equal(levels(mc$mod_type), c("4mC", "5mC", "6mA"))
})

test_that("validity rejects methylation beta values outside [0, 1]", {
    obj <- .make_minimal_commaData()
    assay(obj, "methylation")[1, 1] <- 1.2
    expect_error(validObject(obj), regexp = "methylation.*\\[0, 1\\]")
})

test_that("validity rejects negative or non-integer-like coverage values", {
    obj <- .make_minimal_commaData()
    assay(obj, "coverage")[1, 1] <- -1
    expect_error(validObject(obj), regexp = "coverage.*non-negative")

    obj <- .make_minimal_commaData()
    assay(obj, "coverage")[1, 1] <- 1.5
    expect_error(validObject(obj), regexp = "coverage.*integer-like")
})

test_that("validity rejects invalid raw count assay values", {
    obj <- .make_minimal_commaData()
    assay(obj, "mod_counts")[1, 1] <- -1
    expect_error(validObject(obj), regexp = "mod_counts.*non-negative")

    obj <- .make_minimal_commaData()
    assay(obj, "canonical_counts")[1, 1] <- 1.5
    expect_error(validObject(obj), regexp = "canonical_counts.*integer-like")
})

test_that("validity rejects raw count sums above coverage", {
    obj <- .make_minimal_commaData()
    assay(obj, "mod_counts")[1, 1] <- assay(obj, "coverage")[1, 1]
    assay(obj, "canonical_counts")[1, 1] <- 1L
    expect_error(validObject(obj), regexp = "mod_counts.*canonical_counts.*coverage")
})

test_that("commaData() rejects duplicate parsed site identities before merge assignment", {
    f <- tempfile(fileext = ".bed")
    row <- paste(
        "chr_sim", 99L, 100L, "a,GATC,1", 0L, "+",
        99L, 100L, "0,0,0", 20L, 50, 10L, 10L, 0L, 0L, 0L, 0L, 0L,
        sep = "\t"
    )
    writeLines(c(row, row), f)

    expect_error(
        commaData(
            files = c(s1 = f),
            colData = data.frame(
                sample_name = "s1", condition = "control", replicate = 1L,
                stringsAsFactors = FALSE
            ),
            genome = c(chr_sim = 100000L),
            caller = "modkit"
        ),
        regexp = "duplicate methylation site rows"
    )
})
