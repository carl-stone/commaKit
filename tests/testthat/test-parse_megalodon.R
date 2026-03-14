## Tests for the internal Megalodon parser:
##   .parseMegalodon()

library(testthat)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

# Write a temporary Megalodon per-read BED file.
# Minimum Megalodon format: chrom, start, end, read_id, score, strand, ..., mod_prob
# The parser reads col1=chrom, col2=start, col6=strand, last_col=mod_prob.
.write_tmp_megalodon <- function(rows, file = tempfile(fileext = ".bed")) {
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

# Build one per-read row (7-column minimal Megalodon BED)
.megalodon_row <- function(chrom    = "chr1",
                            start    = 99L,
                            strand   = "+",
                            read_id  = "read1",
                            mod_prob = 0.9) {
    data.frame(chrom, start, start + 1L, read_id, 255L, strand, mod_prob,
               stringsAsFactors = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# .parseMegalodon() — output structure
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseMegalodon() returns a data.frame with correct columns", {
    f <- .write_tmp_megalodon(.megalodon_row())
    result <- comma:::.parseMegalodon(f, "s1", mod_type = "6mA", min_coverage = 1L)
    expect_true(is.data.frame(result))
    expect_named(result, c("chrom", "position", "strand", "mod_type", "beta", "coverage"))
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseMegalodon() — coordinate and value handling
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseMegalodon() converts 0-based start to 1-based position", {
    f <- .write_tmp_megalodon(.megalodon_row(start = 99L))
    result <- comma:::.parseMegalodon(f, "s1", mod_type = "6mA", min_coverage = 1L)
    expect_equal(result$position, 100L)
})

test_that(".parseMegalodon() assigns mod_type from the argument", {
    f      <- .write_tmp_megalodon(.megalodon_row())
    result <- comma:::.parseMegalodon(f, "s1", mod_type = "5mC", min_coverage = 1L)
    expect_equal(result$mod_type, "5mC")
})

test_that(".parseMegalodon() preserves chrom and strand", {
    f      <- .write_tmp_megalodon(.megalodon_row(chrom = "chrX", strand = "-"))
    result <- comma:::.parseMegalodon(f, "s1", mod_type = "6mA", min_coverage = 1L)
    expect_equal(result$chrom,  "chrX")
    expect_equal(result$strand, "-")
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseMegalodon() — per-read → per-site aggregation
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseMegalodon() computes beta as mean of per-read mod_prob", {
    rows <- rbind(
        .megalodon_row(start = 99L, read_id = "r1", mod_prob = 0.8),
        .megalodon_row(start = 99L, read_id = "r2", mod_prob = 0.6)
    )
    f      <- .write_tmp_megalodon(rows)
    result <- comma:::.parseMegalodon(f, "s1", mod_type = "6mA", min_coverage = 1L)
    # mean(0.8, 0.6) = 0.7
    expect_equal(result$beta, 0.7, tolerance = 1e-6)
})

test_that(".parseMegalodon() computes coverage as read count per site", {
    rows <- rbind(
        .megalodon_row(start = 99L, read_id = "r1", mod_prob = 0.9),
        .megalodon_row(start = 99L, read_id = "r2", mod_prob = 0.8),
        .megalodon_row(start = 99L, read_id = "r3", mod_prob = 0.7)
    )
    f      <- .write_tmp_megalodon(rows)
    result <- comma:::.parseMegalodon(f, "s1", mod_type = "6mA", min_coverage = 1L)
    expect_equal(result$coverage, 3L)
})

test_that(".parseMegalodon() aggregates multiple sites independently", {
    rows <- rbind(
        .megalodon_row(start = 99L,  read_id = "r1", mod_prob = 1.0),
        .megalodon_row(start = 99L,  read_id = "r2", mod_prob = 0.0),
        .megalodon_row(start = 199L, read_id = "r3", mod_prob = 0.6),
        .megalodon_row(start = 199L, read_id = "r4", mod_prob = 0.4)
    )
    f      <- .write_tmp_megalodon(rows)
    result <- comma:::.parseMegalodon(f, "s1", mod_type = "6mA", min_coverage = 1L)
    result <- result[order(result$position), ]

    expect_equal(nrow(result), 2L)
    expect_equal(result$position,  c(100L, 200L))
    expect_equal(result$beta[1],   0.5, tolerance = 1e-6)  # (1.0 + 0.0) / 2
    expect_equal(result$beta[2],   0.5, tolerance = 1e-6)  # (0.6 + 0.4) / 2
    expect_equal(result$coverage,  c(2L, 2L))
})

test_that(".parseMegalodon() treats the same position on different strands as separate sites", {
    rows <- rbind(
        .megalodon_row(start = 99L, strand = "+", read_id = "r1", mod_prob = 0.9),
        .megalodon_row(start = 99L, strand = "-", read_id = "r2", mod_prob = 0.1)
    )
    f      <- .write_tmp_megalodon(rows)
    result <- comma:::.parseMegalodon(f, "s1", mod_type = "6mA", min_coverage = 1L)
    # Same position, different strand → two separate sites
    expect_equal(nrow(result), 2L)
    expect_setequal(result$strand, c("-", "+"))
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseMegalodon() — min_coverage filtering
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseMegalodon() drops sites below min_coverage", {
    rows <- rbind(
        # site at position 100: 2 reads — below threshold of 5
        .megalodon_row(start = 99L,  read_id = "r1", mod_prob = 0.8),
        .megalodon_row(start = 99L,  read_id = "r2", mod_prob = 0.9),
        # site at position 200: 6 reads — above threshold
        .megalodon_row(start = 199L, read_id = "r3", mod_prob = 0.7),
        .megalodon_row(start = 199L, read_id = "r4", mod_prob = 0.8),
        .megalodon_row(start = 199L, read_id = "r5", mod_prob = 0.9),
        .megalodon_row(start = 199L, read_id = "r6", mod_prob = 0.6),
        .megalodon_row(start = 199L, read_id = "r7", mod_prob = 0.7),
        .megalodon_row(start = 199L, read_id = "r8", mod_prob = 0.8)
    )
    f      <- .write_tmp_megalodon(rows)
    result <- comma:::.parseMegalodon(f, "s1", mod_type = "6mA", min_coverage = 5L)
    expect_equal(nrow(result),    1L)
    expect_equal(result$position, 200L)
})

test_that(".parseMegalodon() retains all sites when min_coverage = 0", {
    rows <- rbind(
        .megalodon_row(start = 99L,  read_id = "r1", mod_prob = 0.9),
        .megalodon_row(start = 199L, read_id = "r2", mod_prob = 0.5)
    )
    f      <- .write_tmp_megalodon(rows)
    result <- comma:::.parseMegalodon(f, "s1", mod_type = "6mA", min_coverage = 0L)
    expect_equal(nrow(result), 2L)
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseMegalodon() — mod_type defaulting
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseMegalodon() warns and defaults to '6mA' when mod_type is NULL", {
    f <- .write_tmp_megalodon(.megalodon_row())
    expect_warning(
        result <- comma:::.parseMegalodon(f, "s1", mod_type = NULL, min_coverage = 1L),
        regexp = "6mA"
    )
    expect_equal(result$mod_type, "6mA")
})

# ─────────────────────────────────────────────────────────────────────────────
# .parseMegalodon() — error handling
# ─────────────────────────────────────────────────────────────────────────────

test_that(".parseMegalodon() errors on non-existent file", {
    expect_error(
        comma:::.parseMegalodon("/nonexistent/path/file.bed", "s1", mod_type = "6mA"),
        regexp = "not found"
    )
})

test_that(".parseMegalodon() errors on file with fewer than 7 columns", {
    f <- tempfile(fileext = ".bed")
    writeLines("chr1\t99\t100\tread1\t255\t+", f)  # only 6 columns
    expect_error(
        comma:::.parseMegalodon(f, "s1", mod_type = "6mA"),
        regexp = "7"
    )
})

test_that(".parseMegalodon() errors on non-character file argument", {
    expect_error(
        comma:::.parseMegalodon(123L, "s1", mod_type = "6mA"),
        regexp = "character"
    )
})
