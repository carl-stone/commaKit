# Tests for .parseDorado() and its internal helpers — Phase 4

# ─── .cigarToRefPos() ────────────────────────────────────────────────────────

test_that("cigarToRefPos: simple all-match CIGAR", {
    # 5M: 5 match operations from ref position 100
    result <- comma:::.cigarToRefPos("5M", ref_start = 100L, seq_bases = "ACGTA")
    expect_equal(result, 100L:104L)
})

test_that("cigarToRefPos: CIGAR with deletion", {
    # 3M2D2M: 3 matches, 2 deleted ref bases, 2 more matches
    # read: ACGTA (5 bases), ref positions: 100,101,102, skip 103+104, 105,106
    result <- comma:::.cigarToRefPos("3M2D2M", ref_start = 100L, seq_bases = "ACGTA")
    expect_equal(result, c(100L, 101L, 102L, 105L, 106L))
})

test_that("cigarToRefPos: CIGAR with insertion", {
    # 3M2I2M: 3 matches, 2 inserted read bases (no ref pos), 2 more matches
    # read: ACGTAA (7 bases), ref positions: 100,101,102, NA,NA, 103,104
    result <- comma:::.cigarToRefPos("3M2I2M", ref_start = 100L, seq_bases = "ACGTAAG")
    expect_equal(result, c(100L, 101L, 102L, NA_integer_, NA_integer_, 103L, 104L))
})

test_that("cigarToRefPos: CIGAR with soft clip", {
    # 2S3M: 2 soft-clipped bases (NA ref pos), 3 matches
    result <- comma:::.cigarToRefPos("2S3M", ref_start = 50L, seq_bases = "XXACG")
    expect_equal(result, c(NA_integer_, NA_integer_, 50L, 51L, 52L))
})

test_that("cigarToRefPos: returns NULL for empty cigar", {
    result <- comma:::.cigarToRefPos("", ref_start = 1L, seq_bases = "A")
    expect_null(result)
})

# ─── .parseMmTag() ────────────────────────────────────────────────────────────

test_that("parseMmTag: parses single-mod-type MM tag", {
    # MM: "A+a?,0" means first A in read is modified (6mA)
    # ML: 200 → probability 200/255 ≈ 0.78 → is_mod = TRUE
    result <- comma:::.parseMmTag(
        mm_tag   = "A+a?,0",
        ml_tag   = as.raw(200L),
        seq_bases = "ACGT"
    )
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 1L)
    expect_equal(result$mod_type, "6mA")
    expect_equal(result$read_pos, 1L)  # first A at position 1
    expect_true(result$is_mod)
})

test_that("parseMmTag: ML probability <= 0.5 → is_mod = FALSE", {
    # ML: 100 → 100/255 ≈ 0.39 < 0.5 → not modified
    result <- comma:::.parseMmTag(
        mm_tag    = "A+a?,0",
        ml_tag    = as.raw(100L),
        seq_bases = "ACGT"
    )
    expect_false(result$is_mod)
})

test_that("parseMmTag: delta offset positions multiple modifications", {
    # MM: "A+a?,0,1" — two 6mA modifications
    # delta 0 → first A; delta 1 → skip 1 A, so third A overall (2nd A skipped)
    # Sequence "AACGTA": A at pos 1,2,6; delta 0 → pos1, delta 1 → pos6
    result <- comma:::.parseMmTag(
        mm_tag    = "A+a?,0,1",
        ml_tag    = as.raw(c(200L, 210L)),
        seq_bases = "AACGTA"
    )
    expect_equal(nrow(result), 2L)
    expect_equal(result$read_pos, c(1L, 6L))
})

test_that("parseMmTag: returns NULL for null tags", {
    expect_null(comma:::.parseMmTag(NULL, as.raw(200L), "ACGT"))
    expect_null(comma:::.parseMmTag("A+a?,0", NULL, "ACGT"))
})

test_that("parseMmTag: skips unknown mod_code", {
    # mod_code "z" is not in .MODKIT_CODE_MAP
    result <- comma:::.parseMmTag(
        mm_tag    = "A+z?,0",
        ml_tag    = as.raw(200L),
        seq_bases = "ACGT"
    )
    expect_null(result)
})

# ─── .parseDorado() error handling ────────────────────────────────────────────

test_that("parseDorado: error on missing file", {
    expect_error(
        comma:::.parseDorado("/nonexistent/path.bam", "sample1"),
        "not found"
    )
})

test_that("parseDorado: error on non-character file argument", {
    expect_error(
        comma:::.parseDorado(123L, "sample1"),
        "must be a single character string"
    )
})

test_that("parseDorado: error on length-2 file argument", {
    expect_error(
        comma:::.parseDorado(c("a.bam", "b.bam"), "sample1"),
        "must be a single character string"
    )
})

# ─── .cigarToRefPos() additional edge cases ───────────────────────────────────

test_that("cigarToRefPos: hard clip (H) does not consume read positions", {
    # 2H5M: 2 hard-clipped bases (not in seq), then 5 matches
    # Hard clips are NOT in the sequence string, so seq has 5 bases
    result <- comma:::.cigarToRefPos("2H5M", ref_start = 10L, seq_bases = "ACGTA")
    # All 5 read positions map to ref 10..14
    expect_equal(result, 10L:14L)
})

test_that("cigarToRefPos: all soft-clip returns all NA", {
    # 5S: entire read is soft-clipped; no ref positions assigned
    result <- comma:::.cigarToRefPos("5S", ref_start = 1L, seq_bases = "ACGTA")
    expect_equal(result, rep(NA_integer_, 5L))
})

test_that("cigarToRefPos: mixed soft-clip and match assigns NA then ref positions", {
    # 3S2M: 3 soft-clipped (NA), then 2 matches starting at ref 50
    result <- comma:::.cigarToRefPos("3S2M", ref_start = 50L, seq_bases = "XXXAC")
    expect_equal(result, c(NA_integer_, NA_integer_, NA_integer_, 50L, 51L))
})

test_that("cigarToRefPos: N (skipped region) advances reference like deletion", {
    # 2M3N2M: 2 matches, skip 3 ref bases (intron), 2 more matches
    result <- comma:::.cigarToRefPos("2M3N2M", ref_start = 100L, seq_bases = "ACGT")
    expect_equal(result, c(100L, 101L, 105L, 106L))
})

# ─── .parseMmTag() additional tests ──────────────────────────────────────────

test_that("parseMmTag: 5mC (mod_code 'm') is recognised", {
    # MM: "C+m?,0" — first C in read is potentially 5-methylcytosine
    result <- comma:::.parseMmTag(
        mm_tag    = "C+m?,0",
        ml_tag    = as.raw(200L),
        seq_bases = "ACGT"
    )
    expect_false(is.null(result))
    expect_equal(result$mod_type[1], "5mC")
})

test_that("parseMmTag: multi-mod-type MM tag returns rows for each type", {
    # MM: "A+a?,0;C+m?,0" — one 6mA and one 5mC in same read
    # seq: "ACGT" — A at pos 1, C at pos 2
    result <- comma:::.parseMmTag(
        mm_tag    = "A+a?,0;C+m?,0",
        ml_tag    = as.raw(c(200L, 200L)),
        seq_bases = "ACGT"
    )
    expect_equal(nrow(result), 2L)
    expect_true("6mA" %in% result$mod_type)
    expect_true("5mC" %in% result$mod_type)
})

test_that("parseMmTag: zero modifications encoded (empty delta list) returns NULL", {
    # MM: "A+a?" with no deltas means zero modifications for this type
    result <- comma:::.parseMmTag(
        mm_tag    = "A+a?",
        ml_tag    = as.raw(integer(0)),
        seq_bases = "ACGT"
    )
    # No deltas parsed → NULL
    expect_null(result)
})

test_that("parseMmTag: ML at exactly boundary (128/255 ≈ 0.502) is_mod = TRUE", {
    # 128/255 ≈ 0.502 > 0.5 → is_mod TRUE
    result <- comma:::.parseMmTag(
        mm_tag    = "A+a?,0",
        ml_tag    = as.raw(128L),
        seq_bases = "ACGT"
    )
    expect_true(result$is_mod[1])
})

test_that("parseMmTag: ML at exactly 127/255 ≈ 0.498 is_mod = FALSE", {
    # 127/255 ≈ 0.498 ≤ 0.5 → is_mod FALSE
    result <- comma:::.parseMmTag(
        mm_tag    = "A+a?,0",
        ml_tag    = as.raw(127L),
        seq_bases = "ACGT"
    )
    expect_false(result$is_mod[1])
})
