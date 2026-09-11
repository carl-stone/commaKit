# Property-based tests: .parseModkit() parser contract.
#
# Invariants under test:
#   1. Valid bedMethyl rows parse to the tidy schema with beta computed
#      from authoritative counts (Nmod / Nvalid_cov).
#   2. Malformed input is rejected with an actionable error naming the
#      file and the problem — never partial silent output.

library(quickcheck)

test_that(".parseModkit() computes beta from counts for random coverage/beta", {
  for_all(
    coverage = integer_bounded(1L, 500L, len = c(1L, 1L)),
    beta = double_bounded(0, 1, len = c(1L, 1L)),
    property = function(coverage, beta) {
      path <- .qc_write_bedmethyl(
        file.path(tempdir(), paste0("qc-parse-", sample.int(2^30, 1L), ".bed")),
        .qc_make_bedmethyl_row(
          chrom = "chr_qc",
          position = 1000L,
          strand = "+",
          mod_type = "6mA",
          motif = "GATC",
          coverage = coverage,
          beta = beta
        )
      )
      df <- expect_no_warning(.parseModkit(path, "s1", min_coverage = 1L))
      expect_identical(nrow(df), 1L)
      # beta is Nmod / Nvalid_cov — the counts are authoritative
      nmod <- round(coverage * beta)
      expect_equal(df$beta, nmod / coverage)
      expect_identical(df$coverage, as.integer(coverage))
      expect_identical(df$mod_counts, as.integer(nmod))
    }
  )
})

test_that(".parseModkit() handles all valid mod codes", {
  for_all(
    mod_type = one_of(
      quickcheck::constant("6mA"),
      quickcheck::constant("5mC"),
      quickcheck::constant("4mC")
    ),
    property = function(mod_type) {
      path <- .qc_write_bedmethyl(
        file.path(tempdir(), paste0("qc-mod-", sample.int(2^30, 1L), ".bed")),
        .qc_make_bedmethyl_row(
          chrom = "chr_qc",
          position = 1000L,
          strand = "+",
          mod_type = mod_type,
          motif = "GATC",
          coverage = 20L,
          beta = 0.5
        )
      )
      df <- .parseModkit(path, "s1", min_coverage = 1L)
      expect_identical(df$mod_type, mod_type)
    }
  )
})

test_that(".parseModkit() rejects files with fewer than 18 columns", {
  for_all(
    n_cols = integer_bounded(1L, 17L, len = c(1L, 1L)),
    property = function(n_cols) {
      row <- strsplit(
        .qc_make_bedmethyl_row(
          chrom = "chr_qc",
          position = 1000L,
          strand = "+",
          mod_type = "6mA",
          motif = "GATC",
          coverage = 20L,
          beta = 0.5
        ),
        "\t"
      )[[1L]]
      path <- .qc_write_bedmethyl(
        file.path(tempdir(), paste0("qc-trunc-", sample.int(2^30, 1L), ".bed")),
        paste(utils::head(row, n_cols), collapse = "\t")
      )
      expect_error(
        .parseModkit(path, "s1"),
        regexp = "expected at least 18"
      )
    }
  )
})

test_that(".parseModkit() rejects rows with missing required fields", {
  for_all(
    property = function() {
      row <- strsplit(
        .qc_make_bedmethyl_row(
          chrom = "chr_qc",
          position = 1000L,
          strand = "+",
          mod_type = "6mA",
          motif = "GATC",
          coverage = 20L,
          beta = 0.5
        ),
        "\t"
      )[[1L]]
      # blank out Nmod (col 12)
      row[[12L]] <- ""
      path <- .qc_write_bedmethyl(
        file.path(tempdir(), paste0("qc-blank-", sample.int(2^30, 1L), ".bed")),
        paste(row, collapse = "\t")
      )
      expect_error(
        .parseModkit(path, "s1"),
        regexp = "missing required field"
      )
    }
  )
})

test_that(".parseModkit() warns on unknown mod codes and skips them", {
  for_all(
    property = function() {
      good <- .qc_make_bedmethyl_row(
        chrom = "chr_qc",
        position = 1000L,
        strand = "+",
        mod_type = "6mA",
        motif = "GATC",
        coverage = 20L,
        beta = 0.5
      )
      bad <- .qc_make_bedmethyl_row(
        chrom = "chr_qc",
        position = 2000L,
        strand = "+",
        mod_type = "6mA",
        motif = "GATC",
        coverage = 20L,
        beta = 0.5
      )
      # replace mod code with an unmapped one
      bad <- sub("\ta,GATC,1\t", "\tzzz,GATC,1\t", bad)
      path <- .qc_write_bedmethyl(
        file.path(tempdir(), paste0("qc-unk-", sample.int(2^30, 1L), ".bed")),
        c(good, bad)
      )
      df <- NULL
      expect_warning(
        df <- .parseModkit(path, "s1", min_coverage = 1L),
        regexp = "Unknown mod_code"
      )
      # only the good row survives
      expect_identical(nrow(df), 1L)
      expect_identical(df$position, 1000L)
    }
  )
})

test_that(".parseModkit() min_coverage filter drops low-coverage sites", {
  for_all(
    cov_a = integer_bounded(1L, 4L, len = c(1L, 1L)),
    cov_b = integer_bounded(5L, 50L, len = c(1L, 1L)),
    property = function(cov_a, cov_b) {
      lines <- c(
        .qc_make_bedmethyl_row("chr_qc", 1000L, "+", "6mA", "GATC", cov_a, 0.5),
        .qc_make_bedmethyl_row("chr_qc", 2000L, "+", "6mA", "GATC", cov_b, 0.5)
      )
      path <- .qc_write_bedmethyl(
        file.path(tempdir(), paste0("qc-cov-", sample.int(2^30, 1L), ".bed")),
        lines
      )
      df <- .parseModkit(path, "s1", min_coverage = 5L)
      expect_identical(df$position, 2000L)
      expect_identical(df$coverage, as.integer(cov_b))
    }
  )
})
