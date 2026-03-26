## Tests for the internal parser functions:
##   .parseModkit(), .parseMegalodon(), .parseDorado()

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
                         strand = "+", cov = 20L, mod_freq = 0.9) {
    n_mod <- as.integer(round(mod_freq * cov))
    n_can <- cov - n_mod
    data.frame(chrom, start, start + 1L, mod_code, cov, strand,
               start, start + 1L, "255,0,0",
               cov, mod_freq * 100,
               n_mod, n_can,
               0L, 0L, 0L, 0L, 0L,
               stringsAsFactors = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# .parseModkit() — successful parsing
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() returns correct columns", {
    f <- .write_tmp_modkit(.modkit_row())
    result <- comma:::.parseModkit(f, "s1")
    expect_named(result, c("chrom", "position", "strand", "mod_type", "beta", "coverage"))
})

test_that(".parseModkit() maps mod_code 'a' to '6mA'", {
    f <- .write_tmp_modkit(.modkit_row(mod_code = "a,GATC,1"))
    result <- comma:::.parseModkit(f, "s1")
    expect_equal(result$mod_type, "6mA")
})

test_that(".parseModkit() maps mod_code 'm' to '5mC'", {
    f <- .write_tmp_modkit(.modkit_row(mod_code = "m,CCWGG,1"))
    result <- comma:::.parseModkit(f, "s1")
    expect_equal(result$mod_type, "5mC")
})

test_that(".parseModkit() maps mod_code '21839' to '4mC'", {
    f <- .write_tmp_modkit(.modkit_row(mod_code = "21839,CCWGG,1"))
    result <- comma:::.parseModkit(f, "s1")
    expect_equal(result$mod_type, "4mC")
})

test_that(".parseModkit() converts 0-based start to 1-based position", {
    f <- .write_tmp_modkit(.modkit_row(start = 99L))
    result <- comma:::.parseModkit(f, "s1")
    expect_equal(result$position, 100L)
})

test_that(".parseModkit() preserves beta value correctly", {
    f <- .write_tmp_modkit(.modkit_row(mod_freq = 0.75))
    result <- comma:::.parseModkit(f, "s1")
    expect_equal(result$beta, 0.75, tolerance = 1e-6)
})

test_that(".parseModkit() preserves coverage correctly", {
    f <- .write_tmp_modkit(.modkit_row(cov = 42L))
    result <- comma:::.parseModkit(f, "s1")
    expect_equal(result$coverage, 42L)
})

test_that(".parseModkit() parses multiple mod types correctly", {
    rows <- rbind(
        .modkit_row(mod_code = "a,GATC,1",     start = 99L),
        .modkit_row(mod_code = "m,CCWGG,1",    start = 199L),
        .modkit_row(mod_code = "21839,CCWGG,1", start = 299L)
    )
    f <- .write_tmp_modkit(rows)
    result <- comma:::.parseModkit(f, "s1")
    expect_equal(sort(result$mod_type), c("4mC", "5mC", "6mA"))
    expect_equal(nrow(result), 3L)
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseModkit() — min_coverage filtering
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() drops sites below min_coverage", {
    rows <- rbind(
        .modkit_row(cov = 3L,  start = 99L),   # below threshold
        .modkit_row(cov = 10L, start = 199L)   # above threshold
    )
    f <- .write_tmp_modkit(rows)
    result <- comma:::.parseModkit(f, "s1", min_coverage = 5L)
    expect_equal(nrow(result), 1L)
    expect_equal(result$coverage, 10L)
})

test_that(".parseModkit() keeps all sites when min_coverage = 0", {
    rows <- rbind(
        .modkit_row(cov = 1L, start = 99L),
        .modkit_row(cov = 2L, start = 199L)
    )
    f <- .write_tmp_modkit(rows)
    result <- comma:::.parseModkit(f, "s1", min_coverage = 0L)
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
    result <- comma:::.parseModkit(f, "s1", mod_type = "6mA")
    expect_equal(nrow(result), 1L)
    expect_equal(result$mod_type, "6mA")
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseModkit() — error handling
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() errors on missing file", {
    expect_error(
        comma:::.parseModkit("/nonexistent/path/to/file.bed", "s1"),
        regexp = "not found"
    )
})

test_that(".parseModkit() errors on file with fewer than 18 columns", {
    f <- tempfile(fileext = ".bed")
    writeLines("chr1\t99\t100\ta,GATC,1\t255\t+", f)
    expect_error(
        comma:::.parseModkit(f, "s1"),
        regexp = "18"
    )
})

test_that(".parseModkit() warns on unknown mod_code and drops those rows", {
    rows <- rbind(
        .modkit_row(mod_code = "z,GATC,1", start = 99L),   # unknown
        .modkit_row(mod_code = "a,GATC,1", start = 199L)   # known
    )
    f <- .write_tmp_modkit(rows)
    expect_warning(
        result <- comma:::.parseModkit(f, "s1"),
        regexp = "z"
    )
    expect_equal(nrow(result), 1L)
    expect_equal(result$mod_type, "6mA")
})

test_that(".parseModkit() returns empty data frame for empty file", {
    f <- tempfile(fileext = ".bed")
    writeLines("", f)
    expect_message(
        result <- comma:::.parseModkit(f, "s1"),
        regexp = "no data"
    )
    expect_equal(nrow(result), 0L)
    expect_named(result, c("chrom", "position", "strand", "mod_type", "beta", "coverage"))
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseModkit() — using the package extdata file
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseModkit() parses the bundled example file without error", {
    bed_file <- system.file("extdata", "example_modkit.bed", package = "comma")
    skip_if(bed_file == "", message = "extdata not available")

    result <- comma:::.parseModkit(bed_file, "example")
    expect_true(nrow(result) > 0)
    expect_true(all(c("6mA", "5mC") %in% result$mod_type))
    expect_true(all(result$beta >= 0 & result$beta <= 1))
    expect_true(all(result$coverage > 0))
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseDorado() — input validation (Phase 4: full implementation)
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseDorado() stops with informative error for missing file", {
    expect_error(
        comma:::.parseDorado("/nonexistent/path/to.bam", "s1"),
        regexp = "not found"
    )
})

test_that(".parseDorado() stops with informative error for non-character file", {
    expect_error(
        comma:::.parseDorado(123L, "s1"),
        regexp = "single character string"
    )
})
