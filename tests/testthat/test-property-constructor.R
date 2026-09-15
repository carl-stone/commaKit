# Property-based tests: commaData() constructor contract.
#
# Invariants under test:
#   1. Valid inputs produce a valid commaData object (never a half-built one).
#   2. Invalid inputs are rejected with an actionable error that names the
#      offending argument — never a deep-internals error, never silence.

library(quickcheck)

test_that("commaData() accepts valid files/colData for random sample counts", {
  for_all(
    n_samples = integer_bounded(1L, 6L, len = c(1L, 1L)),
    property = function(n_samples) {
      comps <- .qc_valid_files_coldata(n_samples)
      obj <- expect_no_error(
        suppressMessages(
          commaData(
            files = comps$files,
            colData = comps$colData,
            genome = c(chr_qc = 100000L)
          )
        )
      )
      expect_s4_class(obj, "commaData")
      expect_true(validObject(obj))
      expect_identical(ncol(obj), n_samples)
      # Genome info restricted to data chromosomes
      expect_identical(genomeSizes(obj), c(chr_qc = 100000L))
    }
  )
})

test_that("commaData() min_coverage accepts any positive integer", {
  for_all(
    min_coverage = integer_bounded(1L, 1000L, len = c(1L, 1L)),
    property = function(min_coverage) {
      comps <- .qc_valid_files_coldata(1L)
      obj <- suppressMessages(
        commaData(
          files = comps$files,
          colData = comps$colData,
          genome = c(chr_qc = 100000L),
          min_coverage = min_coverage
        )
      )
      expect_identical(
        S4Vectors::metadata(obj)$min_coverage,
        as.integer(min_coverage)
      )
    }
  )
})

test_that("commaData() rejects colData missing required columns", {
  for_all(
    property = function() {
      expect_error(
        suppressMessages(
          commaData(
            files = stats::setNames("/tmp/x.bed", "s1"),
            colData = .qc_bad_coldata_missing_col()
          )
        ),
        regexp = "colData.*missing required columns"
      )
    }
  )
})

test_that("commaData() rejects unnamed files vector", {
  for_all(
    property = function() {
      expect_error(
        commaData(
          files = .qc_invalid_files_unnamed(),
          colData = data.frame(
            sample_name = c("s1", "s2"),
            replicate = 1:2
          )
        ),
        regexp = "files must be a named character"
      )
    }
  )
})

test_that("commaData() rejects files/colData sample mismatches", {
  for_all(
    property = function() {
      coldata <- data.frame(
        sample_name = c("s1", "s2"),
        replicate = 1:2,
        stringsAsFactors = FALSE
      )
      # file name not in colData
      expect_error(
        commaData(
          files = .qc_bad_files_missing_sample(),
          colData = coldata
        ),
        regexp = "names\\(files\\).*not found in colData"
      )
      # colData sample with no file
      expect_error(
        commaData(
          files = .qc_invalid_files_extra_sample(coldata),
          colData = coldata
        ),
        regexp = "no file in 'files'"
      )
    }
  )
})

test_that("commaData() rejects invalid expected_mod_contexts shapes", {
  for_all(
    property = function() {
      coldata <- data.frame(
        sample_name = "s1",
        replicate = 1L,
        stringsAsFactors = FALSE
      )
      files <- stats::setNames("/tmp/x.bed", "s1")
      # unnamed list
      expect_error(
        commaData(
          files = files,
          colData = coldata,
          expected_mod_contexts = list("GATC")
        ),
        regexp = "expected_mod_contexts"
      )
      # bad mod_type name
      expect_error(
        commaData(
          files = files,
          colData = coldata,
          expected_mod_contexts = list("99mA" = "GATC")
        ),
        regexp = "valid mod_type"
      )
    }
  )
})

test_that("commaData() rejects invalid min_coverage values", {
  for_all(
    bad_coverage = double_(len = c(1L, 1L)),
    property = function(bad_coverage) {
      comps <- .qc_valid_files_coldata(1L)
      if (
        is.na(bad_coverage) ||
          !is.finite(bad_coverage) ||
          bad_coverage < 1 ||
          bad_coverage != floor(bad_coverage)
      ) {
        expect_error(
          commaData(
            files = comps$files,
            colData = comps$colData,
            min_coverage = bad_coverage
          ),
          regexp = "min_coverage.*single positive integer"
        )
      }
    }
  )
})
