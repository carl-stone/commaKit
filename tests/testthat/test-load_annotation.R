## Tests for the exported loadAnnotation() function.
##
## All tests require rtracklayer (Suggests). Each test calls
## skip_if_not_installed("rtracklayer") so the suite degrades gracefully when
## the package is absent (e.g. on minimal CI images).

library(testthat)
library(GenomicRanges)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

# Write a minimal 6-column BED file (chrom, start, end, name, score, strand)
.write_tmp_bed <- function(rows, file = tempfile(fileext = ".bed")) {
    write.table(rows, file = file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    file
}

.bed_row <- function(chrom  = "chr_sim", start = 0L, end = 500L,
                      name   = "regionA", score = 0L, strand = "+") {
    data.frame(chrom, start, end, name, score, strand, stringsAsFactors = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# loadAnnotation() — GFF3: basic structure
# ─────────────────────────────────────────────────────────────────────────────

test_that("loadAnnotation() returns a GRanges from a GFF3 file", {
    skip_if_not_installed("rtracklayer")
    gff <- system.file("extdata", "example.gff3", package = "comma")
    skip_if(gff == "", message = "extdata not available")

    result <- loadAnnotation(gff)
    expect_true(is(result, "GRanges"))
})

test_that("loadAnnotation() GRanges always has 'feature_type' and 'name' mcols", {
    skip_if_not_installed("rtracklayer")
    gff <- system.file("extdata", "example.gff3", package = "comma")
    skip_if(gff == "", message = "extdata not available")

    result <- loadAnnotation(gff)
    expect_true("feature_type" %in% names(mcols(result)))
    expect_true("name"         %in% names(mcols(result)))
})

test_that("loadAnnotation() feature_type reflects the GFF3 type column", {
    skip_if_not_installed("rtracklayer")
    gff <- system.file("extdata", "example.gff3", package = "comma")
    skip_if(gff == "", message = "extdata not available")

    result <- loadAnnotation(gff)
    # example.gff3 contains gene, CDS, rRNA, tRNA features
    expect_true("gene" %in% result$feature_type)
    expect_true("CDS"  %in% result$feature_type)
    expect_true("rRNA" %in% result$feature_type)
    expect_true("tRNA" %in% result$feature_type)
})

test_that("loadAnnotation() name column uses the GFF3 Name attribute", {
    skip_if_not_installed("rtracklayer")
    gff <- system.file("extdata", "example.gff3", package = "comma")
    skip_if(gff == "", message = "extdata not available")

    result <- loadAnnotation(gff)
    # Known gene names from example.gff3
    expect_true("geneA" %in% result$name)
    expect_true("geneB" %in% result$name)
})

test_that("loadAnnotation() returns the expected feature count from example GFF3", {
    skip_if_not_installed("rtracklayer")
    gff <- system.file("extdata", "example.gff3", package = "comma")
    skip_if(gff == "", message = "extdata not available")

    result <- loadAnnotation(gff)
    # example.gff3 has 5 genes + 3 CDS + 1 rRNA + 1 tRNA = 10 features
    expect_equal(length(result), 10L)
})

# ─────────────────────────────────────────────────────────────────────────────
# loadAnnotation() — GFF3: feature_types filtering
# ─────────────────────────────────────────────────────────────────────────────

test_that("loadAnnotation() filters to a single feature_type correctly", {
    skip_if_not_installed("rtracklayer")
    gff <- system.file("extdata", "example.gff3", package = "comma")
    skip_if(gff == "", message = "extdata not available")

    result <- loadAnnotation(gff, feature_types = "gene")
    # example.gff3 has exactly 5 gene features
    expect_equal(length(result), 5L)
    expect_true(all(result$feature_type == "gene"))
})

test_that("loadAnnotation() filters to multiple feature_types correctly", {
    skip_if_not_installed("rtracklayer")
    gff <- system.file("extdata", "example.gff3", package = "comma")
    skip_if(gff == "", message = "extdata not available")

    result <- loadAnnotation(gff, feature_types = c("gene", "CDS"))
    expect_true(all(result$feature_type %in% c("gene", "CDS")))
    expect_false("rRNA" %in% result$feature_type)
    expect_false("tRNA" %in% result$feature_type)
})

test_that("loadAnnotation() warns and returns empty GRanges when feature_types has no matches", {
    skip_if_not_installed("rtracklayer")
    gff <- system.file("extdata", "example.gff3", package = "comma")
    skip_if(gff == "", message = "extdata not available")

    expect_warning(
        result <- loadAnnotation(gff, feature_types = "exon"),
        regexp = "exon"
    )
    expect_equal(length(result), 0L)
})

test_that("loadAnnotation() with NULL feature_types returns all features", {
    skip_if_not_installed("rtracklayer")
    gff <- system.file("extdata", "example.gff3", package = "comma")
    skip_if(gff == "", message = "extdata not available")

    result_null <- loadAnnotation(gff, feature_types = NULL)
    result_all  <- loadAnnotation(gff)
    expect_equal(length(result_null), length(result_all))
})

# ─────────────────────────────────────────────────────────────────────────────
# loadAnnotation() — BED format
# ─────────────────────────────────────────────────────────────────────────────

test_that("loadAnnotation() reads a BED file and returns GRanges", {
    skip_if_not_installed("rtracklayer")
    rows <- rbind(
        .bed_row("chr_sim",  0L,   500L, "regionA"),
        .bed_row("chr_sim", 600L, 1200L, "regionB")
    )
    f      <- .write_tmp_bed(rows)
    result <- loadAnnotation(f)
    expect_true(is(result, "GRanges"))
    expect_equal(length(result), 2L)
})

test_that("loadAnnotation() assigns feature_type 'region' for BED features", {
    skip_if_not_installed("rtracklayer")
    rows   <- .bed_row("chr_sim", 0L, 500L, "regionA")
    f      <- .write_tmp_bed(rows)
    result <- loadAnnotation(f)
    expect_true(all(result$feature_type == "region"))
})

test_that("loadAnnotation() preserves the BED name column", {
    skip_if_not_installed("rtracklayer")
    rows   <- .bed_row("chr_sim", 0L, 500L, "mySpecialRegion")
    f      <- .write_tmp_bed(rows)
    result <- loadAnnotation(f)
    expect_true("mySpecialRegion" %in% result$name)
})

# ─────────────────────────────────────────────────────────────────────────────
# loadAnnotation() — error handling
# ─────────────────────────────────────────────────────────────────────────────

test_that("loadAnnotation() errors on a non-existent file path", {
    skip_if_not_installed("rtracklayer")
    expect_error(
        loadAnnotation("/nonexistent/path/annotation.gff3"),
        regexp = "not found"
    )
})

test_that("loadAnnotation() errors on a non-character file argument", {
    skip_if_not_installed("rtracklayer")
    expect_error(
        loadAnnotation(123L),
        regexp = "character string"
    )
})

test_that("loadAnnotation() errors when a vector of paths is supplied", {
    skip_if_not_installed("rtracklayer")
    expect_error(
        loadAnnotation(c("file1.gff3", "file2.gff3")),
        regexp = "character string"
    )
})
