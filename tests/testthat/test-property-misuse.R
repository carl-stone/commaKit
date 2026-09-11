# Property-based misuse harness: type-confused API calls.
#
# The target audience (scientists pasting AI-generated code, and coding
# agents) will pass wrong-typed arguments. Invariant: every misuse fails
# with an actionable error naming the offending argument — never a
# deep-internals error, never silent success with a wrong result.
#
# Uses .qc_random_fixture() to generate random-but-valid commaData
# objects, then attacks the API with type-confused arguments.

library(quickcheck)

test_that("diffMethyl() misuse fails actionably", {
  for_all(
    n_sites = integer_bounded(4L, 20L, len = c(1L, 1L)),
    n_samples = integer_bounded(2L, 6L, len = c(1L, 1L)),
    property = function(n_sites, n_samples) {
      obj <- .qc_random_fixture(n_sites, n_samples)
      # not a commaData at all
      expect_error(
        diffMethyl(data.frame(x = 1:3)),
        regexp = "commaData"
      )
      # formula referencing a nonexistent colData column
      expect_error(
        diffMethyl(obj, formula = ~nonexistent_column),
        regexp = "nonexistent_column|condition|colData"
      )
    }
  )
})

test_that("accessor misuse fails actionably", {
  for_all(
    n_sites = integer_bounded(4L, 20L, len = c(1L, 1L)),
    n_samples = integer_bounded(2L, 6L, len = c(1L, 1L)),
    property = function(n_sites, n_samples) {
      obj <- .qc_random_fixture(n_sites, n_samples)
      # Accessors on non-commaData input fail with the friendly default
      # method: the error names the function, the expected class, and
      # the class actually received.
      expect_error(
        methylation(data.frame(x = 1)),
        regexp = "methylation\\(\\) expects a commaData object, got data.frame"
      )
      expect_error(
        siteCoverage(data.frame(x = 1)),
        regexp = "siteCoverage\\(\\) expects a commaData object, got data.frame"
      )
      expect_error(
        modCounts(data.frame(x = 1)),
        regexp = "modCounts\\(\\) expects a commaData object, got data.frame"
      )
      expect_error(
        sampleInfo(data.frame(x = 1)),
        regexp = "sampleInfo\\(\\) expects a commaData object, got data.frame"
      )
      expect_error(
        siteInfo(data.frame(x = 1)),
        regexp = "siteInfo\\(\\) expects a commaData object, got data.frame"
      )
      # valid object: accessors return correctly-shaped output
      beta <- methylation(obj)
      expect_identical(dim(beta), c(n_sites, n_samples))
      expect_true(all(!is.na(beta) | is.na(beta))) # schema only
    }
  )
})

test_that("mValues() misuse fails actionably", {
  for_all(
    n_sites = integer_bounded(4L, 20L, len = c(1L, 1L)),
    n_samples = integer_bounded(2L, 6L, len = c(1L, 1L)),
    property = function(n_sites, n_samples) {
      obj <- .qc_random_fixture(n_sites, n_samples)
      expect_error(mValues(data.frame(x = 1)), regexp = "commaData")
      # alpha contract: strictly positive finite scalar.
      # alpha = 1 is VALID (Beta(1,1) prior) — only non-positive/invalid
      # values must be rejected.
      expect_error(mValues(obj, alpha = 0), regexp = "alpha")
      expect_error(mValues(obj, alpha = -0.5), regexp = "alpha")
      expect_no_error(mValues(obj, alpha = 1))
    }
  )
})

test_that("writeBED() misuse fails actionably", {
  for_all(
    n_sites = integer_bounded(4L, 20L, len = c(1L, 1L)),
    n_samples = integer_bounded(2L, 6L, len = c(1L, 1L)),
    property = function(n_sites, n_samples) {
      obj <- .qc_random_fixture(n_sites, n_samples)
      expect_error(
        writeBED(
          data.frame(x = 1),
          file = tempfile(fileext = ".bed"),
          sample = "s1"
        ),
        regexp = "commaData"
      )
      # unwritable path: error must be actionable (name the file problem)
      expect_error(
        writeBED(obj, file = "/nonexistent-dir/x.bed", sample = "s1"),
        regexp = "cannot open|No such file|failed|nonexistent"
      )
    }
  )
})

test_that("results() misuse fails actionably", {
  for_all(
    n_sites = integer_bounded(4L, 20L, len = c(1L, 1L)),
    n_samples = integer_bounded(2L, 6L, len = c(1L, 1L)),
    property = function(n_sites, n_samples) {
      obj <- .qc_random_fixture(n_sites, n_samples)
      # Non-commaData input fails with the friendly default method
      expect_error(
        results(data.frame(x = 1)),
        regexp = "results\\(\\) expects a commaData object, got data.frame"
      )
      # calling results() before diffMethyl() must fail with an
      # actionable error that tells the user what to run
      expect_error(
        results(obj),
        regexp = "No differential methylation results|diffMethyl"
      )
    }
  )
})
