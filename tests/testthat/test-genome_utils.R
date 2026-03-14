## Tests for internal genome utility functions:
##   .validateGenomeInfo(), .circularIndex(), .makeSeqinfo()

library(testthat)
library(GenomicRanges)

# ─────────────────────────────────────────────────────────────────────────────
# .validateGenomeInfo() — NULL input
# ─────────────────────────────────────────────────────────────────────────────

test_that(".validateGenomeInfo() returns NULL for NULL input", {
    expect_null(comma:::.validateGenomeInfo(NULL))
})

# ─────────────────────────────────────────────────────────────────────────────
# .validateGenomeInfo() — named integer / numeric vectors
# ─────────────────────────────────────────────────────────────────────────────

test_that(".validateGenomeInfo() returns a named integer vector unchanged", {
    gi <- c(chr1 = 1000L, chr2 = 2000L)
    result <- comma:::.validateGenomeInfo(gi)
    expect_identical(result, gi)
})

test_that(".validateGenomeInfo() coerces a named numeric vector to integer", {
    gi_num <- c(chr1 = 1000, chr2 = 2000)
    result  <- comma:::.validateGenomeInfo(gi_num)
    expect_true(is.integer(result))
    expected <- as.integer(unname(gi_num))
    expect_equal(unname(result), expected)
    expect_equal(names(result), names(gi_num))
})

test_that(".validateGenomeInfo() errors on unnamed integer vector", {
    expect_error(
        comma:::.validateGenomeInfo(c(1000L, 2000L)),
        regexp = "named integer"
    )
})

test_that(".validateGenomeInfo() errors on invalid type (list)", {
    expect_error(
        comma:::.validateGenomeInfo(list(chr1 = 1000L)),
        regexp = "BSgenome"
    )
})

# ─────────────────────────────────────────────────────────────────────────────
# .validateGenomeInfo() — FASTA path
# ─────────────────────────────────────────────────────────────────────────────

test_that(".validateGenomeInfo() errors on non-existent FASTA path", {
    expect_error(
        comma:::.validateGenomeInfo("/nonexistent/path/genome.fa"),
        regexp = "not found"
    )
})

test_that(".validateGenomeInfo() reads FASTA and returns named integer vector", {
    skip_if_not_installed("Biostrings")
    fa <- tempfile(fileext = ".fa")
    writeLines(c(">chr_a", "ATCGATCG", ">chr_b", "GGGGCCCC"), fa)

    result <- comma:::.validateGenomeInfo(fa)

    expect_true(is.integer(result))
    expect_equal(names(result), c("chr_a", "chr_b"))
    expect_equal(result[["chr_a"]], 8L)
    expect_equal(result[["chr_b"]], 8L)
})

test_that(".validateGenomeInfo() FASTA result has correct sizes for unequal chromosomes", {
    skip_if_not_installed("Biostrings")
    fa <- tempfile(fileext = ".fa")
    writeLines(c(">chrA", "AAAA", ">chrB", "TTTTTTTTTT"), fa)

    result <- comma:::.validateGenomeInfo(fa)

    expect_equal(result[["chrA"]], 4L)
    expect_equal(result[["chrB"]], 10L)
})

# ─────────────────────────────────────────────────────────────────────────────
# .circularIndex() — valid (in-range) positions
# ─────────────────────────────────────────────────────────────────────────────

test_that(".circularIndex() returns positions within range unchanged", {
    expect_equal(comma:::.circularIndex(1L,   100L), 1L)
    expect_equal(comma:::.circularIndex(50L,  100L), 50L)
    expect_equal(comma:::.circularIndex(100L, 100L), 100L)
})

# ─────────────────────────────────────────────────────────────────────────────
# .circularIndex() — wrapping past the end
# ─────────────────────────────────────────────────────────────────────────────

test_that(".circularIndex() wraps positions past the genome end", {
    expect_equal(comma:::.circularIndex(101L, 100L), 1L)
    expect_equal(comma:::.circularIndex(200L, 100L), 100L)
    expect_equal(comma:::.circularIndex(201L, 100L), 1L)
})

test_that(".circularIndex() wraps position 0 to the last position", {
    # position 0 is one step before position 1 on a circular genome
    expect_equal(comma:::.circularIndex(0L, 100L), 100L)
})

# ─────────────────────────────────────────────────────────────────────────────
# .circularIndex() — vectorised input
# ─────────────────────────────────────────────────────────────────────────────

test_that(".circularIndex() is vectorised over positions", {
    result <- comma:::.circularIndex(c(1L, 100L, 101L, 200L), 100L)
    expect_equal(result, c(1L, 100L, 1L, 100L))
})

# ─────────────────────────────────────────────────────────────────────────────
# .circularIndex() — invalid genome_size
# ─────────────────────────────────────────────────────────────────────────────

test_that(".circularIndex() errors on genome_size = 0", {
    expect_error(comma:::.circularIndex(5L, 0L))
})

test_that(".circularIndex() errors on negative genome_size", {
    expect_error(comma:::.circularIndex(5L, -10L))
})

# ─────────────────────────────────────────────────────────────────────────────
# .makeSeqinfo() — NULL input
# ─────────────────────────────────────────────────────────────────────────────

test_that(".makeSeqinfo() returns NULL for NULL input", {
    expect_null(comma:::.makeSeqinfo(NULL))
})

# ─────────────────────────────────────────────────────────────────────────────
# .makeSeqinfo() — Seqinfo construction
# ─────────────────────────────────────────────────────────────────────────────

test_that(".makeSeqinfo() returns a Seqinfo object", {
    gi     <- c(chr1 = 1000L, chr2 = 2000L)
    result <- comma:::.makeSeqinfo(gi)
    expect_true(is(result, "Seqinfo"))
})

test_that(".makeSeqinfo() has correct seqnames", {
    gi     <- c(chr1 = 1000L, chr2 = 2000L)
    result <- comma:::.makeSeqinfo(gi)
    expect_equal(GenomeInfoDb::seqnames(result), c("chr1", "chr2"))
})

test_that(".makeSeqinfo() has correct seqlengths", {
    gi     <- c(chr1 = 1000L, chr2 = 2000L)
    result <- comma:::.makeSeqinfo(gi)
    expect_equal(GenomeInfoDb::seqlengths(result), c(chr1 = 1000L, chr2 = 2000L))
})

test_that(".makeSeqinfo() records genome_name when provided", {
    gi     <- c(chr1 = 1000L)
    result <- comma:::.makeSeqinfo(gi, genome_name = "test_genome")
    expect_equal(unique(GenomeInfoDb::genome(result)), "test_genome")
})

test_that(".makeSeqinfo() handles a single chromosome correctly", {
    gi     <- c(chr_sim = 100000L)
    result <- comma:::.makeSeqinfo(gi)
    expect_equal(length(result), 1L)
    expect_equal(GenomeInfoDb::seqnames(result), "chr_sim")
    expect_equal(GenomeInfoDb::seqlengths(result), c(chr_sim = 100000L))
})
