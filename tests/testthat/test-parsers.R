## Tests for the internal modkit bedMethyl parser: .parseModkit()

library(testthat)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

# Write a temporary modkit BED file with known content for testing
.write_tmp_modkit <- function(rows, file = tempfile(fileext = ".bed")) {
  write.table(
    rows,
    file      = file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  file
}

# A minimal valid modkit row (18 columns, real modkit bedMethyl format):
# chrom start end mod_code score strand thickStart thickEnd itemRgb
#   Nvalid_cov fraction_modified Nmod Ncanonical Nother_mod Ndelete Nfail Ndiff Nnocall
# mod_code uses compound "code,motif,position" format (e.g. "a,GATC,1")
# fraction_modified is a percentage (0-100)
.modkit_row <- function(chrom = "chr1", start = 99L, mod_code = "a,GATC,1",
                        strand = "+", cov = 20L, mod_freq = 0.9,
                        other_mod = 0L) {
  n_mod <- as.integer(round(mod_freq * cov))
  n_can <- cov - n_mod - as.integer(other_mod)
  data.frame(chrom, start, start + 1L, mod_code, cov, strand,
    start, start + 1L, "255,0,0",
    cov, mod_freq * 100,
    n_mod, n_can,
    as.integer(other_mod), 0L, 0L, 0L, 0L,
    stringsAsFactors = FALSE
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# .parseModkit() — successful parsing
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() returns correct columns", {
  f <- .write_tmp_modkit(.modkit_row())
  result <- commaKit:::.parseModkit(f, "s1")
  expect_named(result, c(
    "chrom", "position", "strand", "mod_type", "motif",
    "beta", "coverage", "mod_counts", "canonical_counts",
    "other_mod_counts"
  ))
})

test_that(".parseModkit() extracts motif from compound mod_code 'a,GATC,1'", {
  f <- .write_tmp_modkit(.modkit_row(mod_code = "a,GATC,1"))
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(result$motif, "GATC")
})

test_that(".parseModkit() extracts motif from compound mod_code 'm,CCWGG,1'", {
  f <- .write_tmp_modkit(.modkit_row(mod_code = "m,CCWGG,1"))
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(result$motif, "CCWGG")
})

test_that(".parseModkit() returns NA motif for simple mod_code without comma", {
  f <- .write_tmp_modkit(.modkit_row(mod_code = "a"))
  result <- commaKit:::.parseModkit(f, "s1")
  expect_true(is.na(result$motif))
})

test_that(".parseModkit() maps mod_code 'a' to '6mA'", {
  f <- .write_tmp_modkit(.modkit_row(mod_code = "a,GATC,1"))
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(result$mod_type, "6mA")
})

test_that(".parseModkit() maps mod_code 'm' to '5mC'", {
  f <- .write_tmp_modkit(.modkit_row(mod_code = "m,CCWGG,1"))
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(result$mod_type, "5mC")
})

test_that(".parseModkit() maps mod_code '21839' to '4mC'", {
  f <- .write_tmp_modkit(.modkit_row(mod_code = "21839,CCWGG,1"))
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(result$mod_type, "4mC")
})

test_that(".parseModkit() converts 0-based start to 1-based position", {
  f <- .write_tmp_modkit(.modkit_row(start = 99L))
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(result$position, 100L)
})

test_that(".parseModkit() preserves beta value correctly", {
  f <- .write_tmp_modkit(.modkit_row(mod_freq = 0.75))
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(result$beta, 0.75, tolerance = 1e-6)
})

test_that(".parseModkit() preserves coverage correctly", {
  f <- .write_tmp_modkit(.modkit_row(cov = 42L))
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(result$coverage, 42L)
})

test_that(".parseModkit() preserves observed count columns", {
  f <- .write_tmp_modkit(
    .modkit_row(cov = 20L, mod_freq = 0.75, other_mod = 2L)
  )
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(result$mod_counts, 15L)
  expect_equal(result$canonical_counts, 3L)
  expect_equal(result$other_mod_counts, 2L)
  expect_equal(
    result$mod_counts + result$canonical_counts + result$other_mod_counts,
    result$coverage
  )
})

test_that(".parseModkit() parses multiple mod types correctly", {
  rows <- rbind(
    .modkit_row(mod_code = "a,GATC,1", start = 99L),
    .modkit_row(mod_code = "m,CCWGG,1", start = 199L),
    .modkit_row(mod_code = "21839,CCWGG,1", start = 299L)
  )
  f <- .write_tmp_modkit(rows)
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(sort(result$mod_type), c("4mC", "5mC", "6mA"))
  expect_equal(nrow(result), 3L)
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseModkit() — min_coverage filtering
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() drops sites below min_coverage", {
  rows <- rbind(
    .modkit_row(cov = 3L, start = 99L), # below threshold
    .modkit_row(cov = 10L, start = 199L) # above threshold
  )
  f <- .write_tmp_modkit(rows)
  result <- commaKit:::.parseModkit(f, "s1", min_coverage = 5L)
  expect_equal(nrow(result), 1L)
  expect_equal(result$coverage, 10L)
})

test_that(".parseModkit() keeps all sites when min_coverage = 0", {
  rows <- rbind(
    .modkit_row(cov = 1L, start = 99L),
    .modkit_row(cov = 2L, start = 199L)
  )
  f <- .write_tmp_modkit(rows)
  result <- commaKit:::.parseModkit(f, "s1", min_coverage = 0L)
  expect_equal(nrow(result), 2L)
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseModkit() — mod_type filtering
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() filters by mod_type when specified", {
  rows <- rbind(
    .modkit_row(mod_code = "a", start = 99L),
    .modkit_row(mod_code = "m", start = 199L)
  )
  f <- .write_tmp_modkit(rows)
  result <- commaKit:::.parseModkit(f, "s1", mod_type = "6mA")
  expect_equal(nrow(result), 1L)
  expect_equal(result$mod_type, "6mA")
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseModkit() — error handling
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() errors on missing file", {
  expect_error(
    commaKit:::.parseModkit("/nonexistent/path/to/file.bed", "s1"),
    regexp = "not found"
  )
})

test_that(".parseModkit() errors on file with fewer than 18 columns", {
  f <- tempfile(fileext = ".bed")
  writeLines("chr1\t99\t100\ta,GATC,1\t255\t+", f)
  expect_error(
    commaKit:::.parseModkit(f, "s1"),
    regexp = "18"
  )
})

test_that(".parseModkit() warns on unknown mod_code and drops those rows", {
  rows <- rbind(
    .modkit_row(mod_code = "z,GATC,1", start = 99L), # unknown
    .modkit_row(mod_code = "a,GATC,1", start = 199L) # known
  )
  f <- .write_tmp_modkit(rows)
  expect_warning(
    result <- commaKit:::.parseModkit(f, "s1"),
    regexp = "z"
  )
  expect_equal(nrow(result), 1L)
  expect_equal(result$mod_type, "6mA")
})

test_that(".parseModkit() returns empty data frame for empty file", {
  f <- tempfile(fileext = ".bed")
  writeLines("", f)
  expect_message(
    result <- commaKit:::.parseModkit(f, "s1"),
    regexp = "no data"
  )
  expect_equal(nrow(result), 0L)
  expect_named(result, c(
    "chrom", "position", "strand", "mod_type", "motif",
    "beta", "coverage", "mod_counts", "canonical_counts",
    "other_mod_counts"
  ))
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseModkit() — using the package extdata file
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() parses the bundled example file without error", {
  bed_file <- system.file("extdata", "example_modkit.bed", package = "commaKit")
  skip_if(bed_file == "", message = "extdata not available")

  result <- commaKit:::.parseModkit(bed_file, "example")
  expect_true(nrow(result) > 0)
  expect_true(all(c("6mA", "5mC") %in% result$mod_type))
  expect_true(all(result$beta >= 0 & result$beta <= 1))
  expect_true(all(result$coverage > 0))
})

# ─────────────────────────────────────────────────────────────────────────────
# Production-like modkit edge cases
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() errors clearly on partial rows in mixed files", {
  complete <- as.character(unlist(.modkit_row(start = 99L), use.names = FALSE))
  partial <- complete[seq_len(13L)]
  f <- tempfile(fileext = ".bed")
  writeLines(
    c(paste(complete, collapse = "\t"), paste(partial, collapse = "\t")),
    f
  )

  expect_error(
    commaKit:::.parseModkit(f, "s1", min_coverage = 0L),
    regexp = "missing required field.*Nother_mod"
  )
})

test_that(".parseModkit() retains zero-coverage rows only when requested", {
  f <- .write_tmp_modkit(.modkit_row(cov = 0L, mod_freq = 0, start = 99L))

  default_result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(nrow(default_result), 0L)
  expect_named(default_result, c(
    "chrom", "position", "strand", "mod_type", "motif",
    "beta", "coverage", "mod_counts", "canonical_counts",
    "other_mod_counts"
  ))

  kept_result <- commaKit:::.parseModkit(f, "s1", min_coverage = 0L)
  expect_equal(nrow(kept_result), 1L)
  expect_equal(kept_result$coverage, 0L)
  expect_equal(kept_result$beta, 0)
  expect_equal(kept_result$mod_counts + kept_result$canonical_counts +
    kept_result$other_mod_counts, 0L)
})

test_that(".parseModkit() preserves unexpected motif strings without schema changes", {
  rows <- rbind(
    .modkit_row(mod_code = "a,GATC,1", start = 99L),
    .modkit_row(mod_code = "a,not-a-standard-motif,1", start = 199L),
    .modkit_row(mod_code = "m,CCWGG,1", start = 299L)
  )
  f <- .write_tmp_modkit(rows)

  result <- commaKit:::.parseModkit(f, "s1", min_coverage = 0L)
  expect_named(result, c(
    "chrom", "position", "strand", "mod_type", "motif",
    "beta", "coverage", "mod_counts", "canonical_counts",
    "other_mod_counts"
  ))
  expect_equal(nrow(result), 3L)
  expect_true("not-a-standard-motif" %in% result$motif)
})

# ─────────────────────────────────────────────────────────────────────────────
# Delimiter structure and blank-field detection (#235)
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() rejects blank required fields without shifting later values", {
  # Row with blank Nvalid_cov (column 10) but valid later fields.
  # With sep = "" the blank would be collapsed and fraction_modified
  # would shift into the Nvalid_cov position, hiding the error.
  row <- paste(c(
    "chr1", "99", "100", "a,GATC,1", "20", "+",
    "99", "100", "255,0,0",
    "", "90", "18", "2", "0", "0", "0", "0", "0"
  ), collapse = "\t")
  f <- tempfile(fileext = ".bed")
  writeLines(row, f)

  expect_error(
    commaKit:::.parseModkit(f, "s1", min_coverage = 0L),
    regexp = "missing required field.*Nvalid_cov"
  )
})

test_that(".parseModkit() rejects space-separated (non-tab) bedMethyl files", {
  row <- paste(c(
    "chr1", "99", "100", "a,GATC,1", "20", "+",
    "99", "100", "255,0,0",
    "20", "90", "18", "2", "0", "0", "0", "0", "0"
  ), collapse = " ")
  f <- tempfile(fileext = ".bed")
  writeLines(row, f)

  expect_error(
    commaKit:::.parseModkit(f, "s1"),
    regexp = "tab-separated"
  )
})

test_that(".parseModkit() preserves all-tab bedMethyl with no blank fields", {
  f <- .write_tmp_modkit(rbind(
    .modkit_row(start = 99L),
    .modkit_row(start = 199L, mod_code = "m,CCWGG,1")
  ))
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(nrow(result), 2L)
  expect_equal(sort(result$mod_type), c("5mC", "6mA"))
})

# ─────────────────────────────────────────────────────────────────────────────
# Authoritative count fields and zero-coverage drops (#237)
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() computes beta from Nmod/Nvalid_cov, not fraction_modified", {
  # Row where fraction_modified (40%) disagrees with Nmod/Nvalid_cov (7/20 = 35%)
  row <- paste(c(
    "chr1", "99", "100", "a,GATC,1", "20", "+",
    "99", "100", "255,0,0",
    "20", "40", "7", "11", "2", "0", "0", "0", "0"
  ), collapse = "\t")
  f <- tempfile(fileext = ".bed")
  writeLines(row, f)

  result <- commaKit:::.parseModkit(f, "s1", min_coverage = 0L)
  expect_equal(nrow(result), 1L)
  expect_equal(result$beta, 0.35, tolerance = 1e-6)
  expect_equal(result$mod_counts, 7L)
  expect_equal(result$coverage, 20L)
})

test_that(".parseModkit() drops zero-coverage rows with undefined fraction_modified", {
  # Row with Nvalid_cov = 0 and blank fraction_modified.
  # fraction_modified is no longer a required field, so the row passes
  # the required-fields check and is dropped by the coverage filter.
  row <- paste(c(
    "chr1", "99", "100", "a,GATC,1", "0", "+",
    "99", "100", "255,0,0",
    "0", "", "0", "0", "0", "0", "0", "0", "0"
  ), collapse = "\t")
  f <- tempfile(fileext = ".bed")
  writeLines(row, f)

  # Default min_coverage = 5 drops the zero-coverage row without error
  result <- commaKit:::.parseModkit(f, "s1")
  expect_equal(nrow(result), 0L)

  # With min_coverage = 0 the row is retained with beta = 0
  kept <- commaKit:::.parseModkit(f, "s1", min_coverage = 0L)
  expect_equal(nrow(kept), 1L)
  expect_equal(kept$coverage, 0L)
  expect_equal(kept$beta, 0)
  expect_equal(kept$mod_counts, 0L)
})
