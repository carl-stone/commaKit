## Tests for internal genome utility functions:
##   .validateGenomeInfo(), .circularIndex(), .makeSeqinfo()

library(testthat)
library(GenomicRanges)

# ─────────────────────────────────────────────────────────────────────────────
# .validateGenomeInfo() — NULL input
# ─────────────────────────────────────────────────────────────────────────────

test_that(".validateGenomeInfo() returns NULL for NULL input", {
    expect_null(commaKit:::.validateGenomeInfo(NULL))
})

# ─────────────────────────────────────────────────────────────────────────────
# .validateGenomeInfo() — named integer / numeric vectors
# ─────────────────────────────────────────────────────────────────────────────

test_that(".validateGenomeInfo() returns a named integer vector unchanged", {
    gi <- c(chr1 = 1000L, chr2 = 2000L)
    result <- commaKit:::.validateGenomeInfo(gi)
    expect_identical(result, gi)
})

test_that(".validateGenomeInfo() coerces a named numeric vector to integer", {
    gi_num <- c(chr1 = 1000, chr2 = 2000)
    result  <- commaKit:::.validateGenomeInfo(gi_num)
    expect_true(is.integer(result))
    expected <- as.integer(unname(gi_num))
    expect_equal(unname(result), expected)
    expect_equal(names(result), names(gi_num))
})

test_that(".validateGenomeInfo() errors on unnamed integer vector", {
    expect_error(
        commaKit:::.validateGenomeInfo(c(1000L, 2000L)),
        regexp = "named integer"
    )
})

test_that(".validateGenomeInfo() errors on invalid type (list)", {
    expect_error(
        commaKit:::.validateGenomeInfo(list(chr1 = 1000L)),
        regexp = "BSgenome"
    )
})

# ─────────────────────────────────────────────────────────────────────────────
# .validateGenomeInfo() — FASTA path
# ─────────────────────────────────────────────────────────────────────────────

test_that(".validateGenomeInfo() errors on non-existent FASTA path", {
    expect_error(
        commaKit:::.validateGenomeInfo("/nonexistent/path/genome.fa"),
        regexp = "not found"
    )
})

test_that(".validateGenomeInfo() reads FASTA and returns named integer vector", {
    skip_if_not_installed("Biostrings")
    fa <- tempfile(fileext = ".fa")
    writeLines(c(">chr_a", "ATCGATCG", ">chr_b", "GGGGCCCC"), fa)

    result <- commaKit:::.validateGenomeInfo(fa)

    expect_true(is.integer(result))
    expect_equal(names(result), c("chr_a", "chr_b"))
    expect_equal(result[["chr_a"]], 8L)
    expect_equal(result[["chr_b"]], 8L)
})

test_that(".validateGenomeInfo() FASTA result has correct sizes for unequal chromosomes", {
    skip_if_not_installed("Biostrings")
    fa <- tempfile(fileext = ".fa")
    writeLines(c(">chrA", "AAAA", ">chrB", "TTTTTTTTTT"), fa)

    result <- commaKit:::.validateGenomeInfo(fa)

    expect_equal(result[["chrA"]], 4L)
    expect_equal(result[["chrB"]], 10L)
})

# ─────────────────────────────────────────────────────────────────────────────
# .validateGenomeInfo() — DNAStringSet input
# ─────────────────────────────────────────────────────────────────────────────

test_that(".validateGenomeInfo() accepts a named DNAStringSet", {
    skip_if_not_installed("Biostrings")
    seqs   <- Biostrings::DNAStringSet(c(chr1 = "ATCGATCG", chr2 = "GGGGCCCC"))
    result <- commaKit:::.validateGenomeInfo(seqs)
    expect_true(is.integer(result))
    expect_equal(names(result), c("chr1", "chr2"))
    expect_equal(result[["chr1"]], 8L)
    expect_equal(result[["chr2"]], 8L)
})

test_that(".validateGenomeInfo() DNAStringSet single-sequence returns correct size", {
    skip_if_not_installed("Biostrings")
    seqs   <- Biostrings::DNAStringSet(c(NC_000913 = "GATCAAAA"))
    result <- commaKit:::.validateGenomeInfo(seqs)
    expect_equal(names(result), "NC_000913")
    expect_equal(result[["NC_000913"]], 8L)
})

test_that(".validateGenomeInfo() errors on unnamed DNAStringSet", {
    skip_if_not_installed("Biostrings")
    seqs <- Biostrings::DNAStringSet(c("ATCGATCG", "GGGGCCCC"))
    expect_error(commaKit:::.validateGenomeInfo(seqs), regexp = "names")
})

# ─────────────────────────────────────────────────────────────────────────────
# .validateGenomeInfo() — DNAString input
# ─────────────────────────────────────────────────────────────────────────────

test_that(".validateGenomeInfo() errors on DNAString with helpful message", {
    skip_if_not_installed("Biostrings")
    seq <- Biostrings::DNAString("ATCGATCG")
    expect_error(
        commaKit:::.validateGenomeInfo(seq),
        regexp = "DNAString"
    )
})

# ─────────────────────────────────────────────────────────────────────────────
# .loadGenomeSequences() — shared sequence resolver
# ─────────────────────────────────────────────────────────────────────────────

test_that(".loadGenomeSequences() reads a FASTA path", {
    skip_if_not_installed("Biostrings")
    fa <- tempfile(fileext = ".fa")
    writeLines(c(">chr_a", "ATCG", ">chr_b", "GGGGCC"), fa)

    seqs <- commaKit:::.loadGenomeSequences(fa)

    expect_s4_class(seqs, "DNAStringSet")
    expect_equal(names(seqs), c("chr_a", "chr_b"))
    expect_equal(as.integer(Biostrings::width(seqs)), c(4L, 6L))
})

test_that(".loadGenomeSequences() passes through named DNAStringSet input", {
    skip_if_not_installed("Biostrings")
    seqs <- Biostrings::DNAStringSet(c(chr1 = "ATCG", chr2 = "GGCC"))

    result <- commaKit:::.loadGenomeSequences(seqs)

    expect_s4_class(result, "DNAStringSet")
    expect_equal(names(result), names(seqs))
    expect_equal(as.character(result), as.character(seqs))
})

test_that(".loadGenomeSequences() requires named sequences", {
    skip_if_not_installed("Biostrings")
    seqs <- Biostrings::DNAStringSet(c("ATCG", "GGCC"))

    expect_error(
        commaKit:::.loadGenomeSequences(seqs),
        "non-empty names"
    )
})

test_that(".validateGenomeInfo() and .loadGenomeSequences() share sequence sizes", {
    skip_if_not_installed("Biostrings")
    seqs <- Biostrings::DNAStringSet(c(chr1 = "ATCG", chr2 = "GGCCCC"))

    loaded <- commaKit:::.loadGenomeSequences(seqs)
    sizes <- commaKit:::.validateGenomeInfo(seqs)

    expect_equal(sizes, stats::setNames(
        as.integer(Biostrings::width(loaded)),
        names(loaded)
    ))
})

# ─────────────────────────────────────────────────────────────────────────────
# .circularIndex() — valid (in-range) positions
# ─────────────────────────────────────────────────────────────────────────────

test_that(".circularIndex() returns positions within range unchanged", {
    expect_equal(commaKit:::.circularIndex(1L,   100L), 1L)
    expect_equal(commaKit:::.circularIndex(50L,  100L), 50L)
    expect_equal(commaKit:::.circularIndex(100L, 100L), 100L)
})

# ─────────────────────────────────────────────────────────────────────────────
# .circularIndex() — wrapping past the end
# ─────────────────────────────────────────────────────────────────────────────

test_that(".circularIndex() wraps positions past the genome end", {
    expect_equal(commaKit:::.circularIndex(101L, 100L), 1L)
    expect_equal(commaKit:::.circularIndex(200L, 100L), 100L)
    expect_equal(commaKit:::.circularIndex(201L, 100L), 1L)
})

test_that(".circularIndex() wraps position 0 to the last position", {
    # position 0 is one step before position 1 on a circular genome
    expect_equal(commaKit:::.circularIndex(0L, 100L), 100L)
})

# ─────────────────────────────────────────────────────────────────────────────
# .circularIndex() — vectorised input
# ─────────────────────────────────────────────────────────────────────────────

test_that(".circularIndex() is vectorised over positions", {
    result <- commaKit:::.circularIndex(c(1L, 100L, 101L, 200L), 100L)
    expect_equal(result, c(1L, 100L, 1L, 100L))
})

# ─────────────────────────────────────────────────────────────────────────────
# .circularIndex() — invalid genome_size
# ─────────────────────────────────────────────────────────────────────────────

test_that(".circularIndex() errors on genome_size = 0", {
    expect_error(commaKit:::.circularIndex(5L, 0L))
})

test_that(".circularIndex() errors on negative genome_size", {
    expect_error(commaKit:::.circularIndex(5L, -10L))
})

# ─────────────────────────────────────────────────────────────────────────────
# .makeSeqinfo() — NULL input
# ─────────────────────────────────────────────────────────────────────────────

test_that(".makeSeqinfo() returns NULL for NULL input", {
    expect_null(commaKit:::.makeSeqinfo(NULL))
})

# ─────────────────────────────────────────────────────────────────────────────
# .makeSeqinfo() — Seqinfo construction
# ─────────────────────────────────────────────────────────────────────────────

test_that(".makeSeqinfo() returns a Seqinfo object", {
    gi     <- c(chr1 = 1000L, chr2 = 2000L)
    result <- commaKit:::.makeSeqinfo(gi)
    expect_true(is(result, "Seqinfo"))
})

test_that(".makeSeqinfo() has correct seqnames", {
    gi     <- c(chr1 = 1000L, chr2 = 2000L)
    result <- commaKit:::.makeSeqinfo(gi)
    expect_equal(GenomeInfoDb::seqnames(result), c("chr1", "chr2"))
})

test_that(".makeSeqinfo() has correct seqlengths", {
    gi     <- c(chr1 = 1000L, chr2 = 2000L)
    result <- commaKit:::.makeSeqinfo(gi)
    expect_equal(GenomeInfoDb::seqlengths(result), c(chr1 = 1000L, chr2 = 2000L))
})

test_that(".makeSeqinfo() defaults chromosomes to circular", {
    gi     <- c(chr1 = 1000L, chr2 = 2000L)
    result <- commaKit:::.makeSeqinfo(gi)
    expect_equal(GenomeInfoDb::isCircular(result), c(chr1 = TRUE, chr2 = TRUE))
})

test_that(".makeSeqinfo() records genome_name when provided", {
    gi     <- c(chr1 = 1000L)
    result <- commaKit:::.makeSeqinfo(gi, genome_name = "test_genome")
    expect_equal(unique(GenomeInfoDb::genome(result)), "test_genome")
})

test_that(".makeSeqinfo() handles a single chromosome correctly", {
    gi     <- c(chr_sim = 100000L)
    result <- commaKit:::.makeSeqinfo(gi)
    expect_equal(length(result), 1L)
    expect_equal(GenomeInfoDb::seqnames(result), "chr_sim")
    expect_equal(GenomeInfoDb::seqlengths(result), c(chr_sim = 100000L))
})
