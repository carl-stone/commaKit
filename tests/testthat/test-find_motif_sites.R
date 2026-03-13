## Tests for the exported findMotifSites() function.
##
## All tests require Biostrings (Suggests). Each test calls
## skip_if_not_installed("Biostrings") so the suite degrades gracefully when
## the package is absent.

library(testthat)
library(GenomicRanges)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

# Write a minimal FASTA file.
# sequences: named character vector, e.g. c(chr_test = "ATGATCG")
.write_tmp_fasta <- function(sequences, file = tempfile(fileext = ".fa")) {
    lines <- character(0)
    for (nm in names(sequences)) {
        lines <- c(lines, paste0(">", nm), sequences[[nm]])
    }
    writeLines(lines, file)
    file
}

# ─────────────────────────────────────────────────────────────────────────────
# findMotifSites() — output structure
# ─────────────────────────────────────────────────────────────────────────────

test_that("findMotifSites() returns a GRanges object", {
    skip_if_not_installed("Biostrings")
    fa     <- .write_tmp_fasta(c(chr_test = "ATCGATCGGG"))
    result <- findMotifSites(fa, "GATC")
    expect_true(is(result, "GRanges"))
})

test_that("findMotifSites() result has a 'motif' mcols column", {
    skip_if_not_installed("Biostrings")
    fa     <- .write_tmp_fasta(c(chr_test = "ATCGATCGGG"))
    result <- findMotifSites(fa, "GATC")
    expect_true("motif" %in% names(mcols(result)))
    expect_true(all(result$motif == "GATC"))
})

test_that("findMotifSites() ranges have width equal to motif length", {
    skip_if_not_installed("Biostrings")
    fa     <- .write_tmp_fasta(c(chr_test = "ATCGATCGGG"))
    result <- findMotifSites(fa, "GATC")
    expect_true(all(width(result) == nchar("GATC")))
})

# ─────────────────────────────────────────────────────────────────────────────
# findMotifSites() — palindromic motif (GATC)
#
# GATC reverse-complement is GATC, so it is a palindrome.
# The function reports hits on both strands at the same positions.
# ─────────────────────────────────────────────────────────────────────────────

test_that("findMotifSites() finds the correct position for a single palindromic motif hit", {
    skip_if_not_installed("Biostrings")
    # "ATCGATCGGG": GATC at positions 4-7
    # A(1)T(2)C(3)G(4)A(5)T(6)C(7)G(8)G(9)G(10)
    fa      <- .write_tmp_fasta(c(chr_test = "ATCGATCGGG"))
    result  <- findMotifSites(fa, "GATC")
    fwd     <- result[as.character(strand(result)) == "+"]
    expect_equal(start(fwd), 4L)
    expect_equal(end(fwd),   7L)
})

test_that("findMotifSites() reports both + and - strand hits for a palindromic motif", {
    skip_if_not_installed("Biostrings")
    fa     <- .write_tmp_fasta(c(chr_test = "ATCGATCGGG"))
    result <- findMotifSites(fa, "GATC")
    expect_true("+" %in% as.character(strand(result)))
    expect_true("-" %in% as.character(strand(result)))
})

test_that("findMotifSites() reports identical positions on both strands for a palindrome", {
    skip_if_not_installed("Biostrings")
    fa     <- .write_tmp_fasta(c(chr_test = "ATCGATCGGG"))
    result <- findMotifSites(fa, "GATC")
    fwd    <- result[as.character(strand(result)) == "+"]
    rev    <- result[as.character(strand(result)) == "-"]
    expect_equal(start(fwd), start(rev))
})

test_that("findMotifSites() finds two palindromic hits in a sequence with two occurrences", {
    skip_if_not_installed("Biostrings")
    # "GATCAGATCG": GATC at positions 1-4 and 6-9
    fa     <- .write_tmp_fasta(c(chr_test = "GATCAGATCG"))
    result <- findMotifSites(fa, "GATC")
    fwd    <- result[as.character(strand(result)) == "+"]
    expect_equal(length(fwd), 2L)
    expect_equal(sort(start(fwd)), c(1L, 6L))
})

# ─────────────────────────────────────────────────────────────────────────────
# findMotifSites() — non-palindromic motif (ACAT / RC = ATGT)
# ─────────────────────────────────────────────────────────────────────────────

test_that("findMotifSites() finds a non-palindromic motif on the + strand", {
    skip_if_not_installed("Biostrings")
    # "ACATCCCATGT": ACAT at positions 1-4
    fa     <- .write_tmp_fasta(c(chr_test = "ACATCCCATGT"))
    result <- findMotifSites(fa, "ACAT")
    fwd    <- result[as.character(strand(result)) == "+"]
    expect_equal(length(fwd), 1L)
    expect_equal(start(fwd),  1L)
})

test_that("findMotifSites() finds the reverse-complement of a non-palindromic motif on the - strand", {
    skip_if_not_installed("Biostrings")
    # "ACATCCCATGT": RC(ACAT)=ATGT found at positions 8-11 on the - strand
    fa     <- .write_tmp_fasta(c(chr_test = "ACATCCCATGT"))
    result <- findMotifSites(fa, "ACAT")
    rev    <- result[as.character(strand(result)) == "-"]
    expect_equal(length(rev), 1L)
    expect_equal(start(rev),  8L)
})

# ─────────────────────────────────────────────────────────────────────────────
# findMotifSites() — IUPAC ambiguity codes
# ─────────────────────────────────────────────────────────────────────────────

test_that("findMotifSites() matches IUPAC W (A or T) in a motif", {
    skip_if_not_installed("Biostrings")
    # CCWGG matches CCAGG (pos 1-5) and CCTGG (pos 6-10) in "CCAGGCCTGG"
    fa     <- .write_tmp_fasta(c(chr_test = "CCAGGCCTGG"))
    result <- findMotifSites(fa, "CCWGG")
    fwd    <- result[as.character(strand(result)) == "+"]
    expect_equal(length(fwd), 2L)
    expect_equal(sort(start(fwd)), c(1L, 6L))
})

# ─────────────────────────────────────────────────────────────────────────────
# findMotifSites() — multi-chromosome genome
# ─────────────────────────────────────────────────────────────────────────────

test_that("findMotifSites() searches across all chromosomes in a multi-chr FASTA", {
    skip_if_not_installed("Biostrings")
    fa     <- .write_tmp_fasta(c(chr1 = "GATCAAAA", chr2 = "AAAGATCA"))
    result <- findMotifSites(fa, "GATC")
    chroms <- as.character(seqnames(result))
    expect_true("chr1" %in% chroms)
    expect_true("chr2" %in% chroms)
})

test_that("findMotifSites() returns a zero-length GRanges when the motif is absent", {
    skip_if_not_installed("Biostrings")
    fa     <- .write_tmp_fasta(c(chr_test = "AAAAAAAAAA"))
    result <- findMotifSites(fa, "GATC")
    expect_equal(length(result), 0L)
    expect_true(is(result, "GRanges"))
})

# ─────────────────────────────────────────────────────────────────────────────
# findMotifSites() — error handling
# ─────────────────────────────────────────────────────────────────────────────

test_that("findMotifSites() errors when genome argument is missing", {
    skip_if_not_installed("Biostrings")
    expect_error(findMotifSites(motif = "GATC"), regexp = "genome")
})

test_that("findMotifSites() errors on a non-existent FASTA path", {
    skip_if_not_installed("Biostrings")
    expect_error(
        findMotifSites("/nonexistent/genome.fa", "GATC"),
        regexp = "not found"
    )
})

test_that("findMotifSites() errors on an empty motif string", {
    skip_if_not_installed("Biostrings")
    fa <- .write_tmp_fasta(c(chr_test = "ATGATCG"))
    expect_error(findMotifSites(fa, ""), regexp = "motif")
})

test_that("findMotifSites() errors on a non-character motif argument", {
    skip_if_not_installed("Biostrings")
    fa <- .write_tmp_fasta(c(chr_test = "ATGATCG"))
    expect_error(findMotifSites(fa, 42L), regexp = "motif")
})
