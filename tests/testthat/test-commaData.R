## Tests for the commaData S4 class, constructor, and show() method
##
## These tests use:
##  - system.file("extdata", "example_modkit.bed", package = "comma") for file-based tests
##  - Direct object construction for class/validity tests (avoids file I/O overhead)

library(testthat)
library(SummarizedExperiment)
library(S4Vectors)
library(GenomicRanges)

# ── Helper: build a minimal valid commaData without file I/O ──────────────────

.make_minimal_commaData <- function(n_sites = 5L, n_samples = 2L) {
    site_keys <- paste0("chr_sim:", seq_len(n_sites) * 100L, ":+:6mA:GATC")
    methyl    <- matrix(runif(n_sites * n_samples, 0.1, 0.9),
                        nrow = n_sites, ncol = n_samples,
                        dimnames = list(site_keys, paste0("s", seq_len(n_samples))))
    cov       <- matrix(as.integer(runif(n_sites * n_samples, 10, 50)),
                        nrow = n_sites, ncol = n_samples,
                        dimnames = list(site_keys, paste0("s", seq_len(n_samples))))
    rd <- S4Vectors::DataFrame(
        chrom       = rep("chr_sim", n_sites),
        position    = seq_len(n_sites) * 100L,
        strand      = rep("+", n_sites),
        mod_type    = rep("6mA", n_sites),
        motif       = rep("GATC", n_sites),
        mod_context = rep("6mA_GATC", n_sites),
        row.names   = site_keys
    )
    cd <- S4Vectors::DataFrame(
        sample_name = paste0("s", seq_len(n_samples)),
        condition   = rep(c("control", "treatment"), length.out = n_samples),
        replicate   = seq_len(n_samples),
        row.names   = paste0("s", seq_len(n_samples))
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = methyl, coverage = cov),
        rowData = rd,
        colData = cd
    )
    new("commaData", se, genomeInfo = c(chr_sim = 100000L),
        annotation = GenomicRanges::GRanges(),
        motifSites = GenomicRanges::GRanges())
}

# ─────────────────────────────────────────────────────────────────────────────
# Class instantiation and validity
# ─────────────────────────────────────────────────────────────────────────────

test_that("minimal commaData object is valid", {
    obj <- .make_minimal_commaData()
    expect_true(is(obj, "commaData"))
    expect_true(is(obj, "SummarizedExperiment"))
    expect_no_error(validObject(obj))
})

test_that("validity rejects missing rowData columns", {
    obj <- .make_minimal_commaData()
    # Remove required column
    rowData(obj)$mod_type <- NULL
    expect_error(validObject(obj), regexp = "mod_type")
})

test_that("validity rejects unrecognized mod_type values", {
    obj <- .make_minimal_commaData()
    rowData(obj)$mod_type <- rep("7mX", nrow(obj))
    expect_error(validObject(obj), regexp = "7mX")
})

test_that("validity rejects missing colData columns", {
    obj <- .make_minimal_commaData()
    colData(obj)$condition <- NULL
    expect_error(validObject(obj), regexp = "condition")
})

test_that("validity rejects non-integer genomeInfo", {
    obj <- .make_minimal_commaData()
    obj@genomeInfo <- list(chr_sim = 100000)   # list, not named integer
    expect_error(validObject(obj), regexp = "genomeInfo")
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
    rd$mod_context[1L] <- "6mA"  # fallback when motif is NA
    rowData(obj) <- rd
    expect_no_error(validObject(obj))
})

test_that("validity passes when all motif values are NA", {
    obj <- .make_minimal_commaData()
    rd  <- rowData(obj)
    rd$motif <- NA_character_
    rd$mod_context <- "6mA"  # fallback when motif is NA
    rowData(obj) <- rd
    expect_no_error(validObject(obj))
})

# ─────────────────────────────────────────────────────────────────────────────
# Constructor: argument validation
# ─────────────────────────────────────────────────────────────────────────────

test_that("commaData() errors on colData missing required columns", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
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
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
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
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
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
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
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
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
    skip_if(bed_file == "", message = "extdata not available")

    cd <- commaData(
        files   = c(s1 = bed_file),
        colData = data.frame(
            sample_name = "s1",
            condition   = "control",
            replicate   = 1L,
            stringsAsFactors = FALSE
        ),
        genome  = c(chr_sim = 100000L),
        caller  = "modkit"
    )

    expect_true(is(cd, "commaData"))
    expect_no_error(validObject(cd))
    expect_true(nrow(cd) > 0)
    expect_equal(ncol(cd), 1L)
    expect_true("methylation" %in% assayNames(cd))
    expect_true("coverage" %in% assayNames(cd))
})

test_that("commaData() applies min_coverage filter correctly", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
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
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
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
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
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
    expect_equal(unique(rowData(cd_6ma)$mod_type), "6mA")
})

# ─────────────────────────────────────────────────────────────────────────────
# expected_mod_contexts constructor filter
# ─────────────────────────────────────────────────────────────────────────────

test_that("commaData: expected_mod_contexts filters to specified contexts", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
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
    # Only 6mA_GATC sites remain
    expect_true(all(rowData(cd_6mA)$mod_context == "6mA_GATC"))
    expect_true(nrow(cd_6mA) < nrow(cd_all))
})

test_that("commaData: expected_mod_contexts accepts multiple mod types", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
    skip_if(bed_file == "", message = "extdata not available")

    cd <- commaData(
        files                 = c(s1 = bed_file),
        colData               = data.frame(sample_name = "s1", condition = "ctrl",
                                           replicate = 1L, stringsAsFactors = FALSE),
        genome                = c(chr_sim = 100000L),
        expected_mod_contexts = list("6mA" = "GATC", "5mC" = "CCWGG")
    )
    expect_true(all(rowData(cd)$mod_context %in% c("6mA_GATC", "5mC_CCWGG")))
})

test_that("commaData: expected_mod_contexts stops if no sites remain", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
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
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
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
# mod_context validity checks
# ─────────────────────────────────────────────────────────────────────────────

test_that("validity fails when mod_context is missing from rowData", {
    obj <- .make_minimal_commaData()
    rd  <- rowData(obj)
    rd$mod_context <- NULL
    rowData(obj) <- rd
    expect_error(validObject(obj), "mod_context")
})

test_that("validity fails when mod_context is inconsistent with mod_type + motif", {
    obj <- .make_minimal_commaData()
    rd  <- rowData(obj)
    rd$mod_context <- rep("5mC_CCWGG", nrow(rd))  # wrong context for 6mA_GATC data
    rowData(obj) <- rd
    expect_error(validObject(obj), "mod_context")
})

test_that("validity fails when mod_context contains NA", {
    obj <- .make_minimal_commaData()
    rd  <- rowData(obj)
    rd$mod_context[1L] <- NA_character_
    rowData(obj) <- rd
    expect_error(validObject(obj), "mod_context")
})
