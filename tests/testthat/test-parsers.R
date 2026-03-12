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

# A minimal valid modkit row (15 columns):
# chrom start end mod_code score strand cov mod_freq n_mod n_can n_other n_del n_fail n_diff n_nocall
.modkit_row <- function(chrom = "chr1", start = 99L, mod_code = "a",
                         strand = "+", cov = 20L, mod_freq = 0.9) {
    data.frame(chrom, start, start + 1L, mod_code, 255L, strand,
               cov, mod_freq,
               as.integer(mod_freq * cov), cov - as.integer(mod_freq * cov),
               0L, 0L, 0L, 0L, 0L)
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
    f <- .write_tmp_modkit(.modkit_row(mod_code = "a"))
    result <- comma:::.parseModkit(f, "s1")
    expect_equal(result$mod_type, "6mA")
})

test_that(".parseModkit() maps mod_code 'm' to '5mC'", {
    f <- .write_tmp_modkit(.modkit_row(mod_code = "m"))
    result <- comma:::.parseModkit(f, "s1")
    expect_equal(result$mod_type, "5mC")
})

test_that(".parseModkit() maps mod_code '21839' to '4mC'", {
    f <- .write_tmp_modkit(.modkit_row(mod_code = "21839"))
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
        .modkit_row(mod_code = "a", start = 99L),
        .modkit_row(mod_code = "m", start = 199L),
        .modkit_row(mod_code = "21839", start = 299L)
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

test_that(".parseModkit() errors on file with fewer than 15 columns", {
    f <- tempfile(fileext = ".bed")
    writeLines("chr1\t99\t100\ta\t255\t+", f)
    expect_error(
        comma:::.parseModkit(f, "s1"),
        regexp = "15"
    )
})

test_that(".parseModkit() warns on unknown mod_code and drops those rows", {
    rows <- rbind(
        .modkit_row(mod_code = "z", start = 99L),    # unknown
        .modkit_row(mod_code = "a", start = 199L)    # known
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
# .parseDorado() — stub error
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseDorado() stops with informative message", {
    expect_error(
        comma:::.parseDorado("any.bam", "s1"),
        regexp = "not yet implemented"
    )
    expect_error(
        comma:::.parseDorado("any.bam", "s1"),
        regexp = "modkit"
    )
})
