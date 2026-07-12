# Tests for .parseDorado() and its internal helpers — Phase 4

.make_dorado_test_bam <- function(records) {
  skip_if_not_installed("Rsamtools")

  sam_file <- tempfile(fileext = ".sam")
  bam_prefix <- tempfile()

  writeLines(c(
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:1000",
    records
  ), sam_file)

  bam_file <- Rsamtools::asBam(
    sam_file,
    destination = bam_prefix,
    overwrite = TRUE
  )
  Rsamtools::indexBam(bam_file)

  bam_file
}

.make_sam_record <- function(qname,
                             pos,
                             cigar,
                             seq,
                             mm,
                             ml,
                             flag = 0L) {
  paste(
    qname,
    flag,
    "chr1",
    pos,
    60L,
    cigar,
    "*",
    0L,
    0L,
    seq,
    "*",
    paste0("MM:Z:", mm),
    paste0("ML:B:C,", paste(ml, collapse = ",")),
    sep = "\t"
  )
}

# ─── .cigarToRefPos() ────────────────────────────────────────────────────────

test_that("cigarToRefPos: simple all-match CIGAR", {
  # 5M: 5 match operations from ref position 100
  result <- commaKit:::.cigarToRefPos(
    "5M",
    ref_start = 100L,
    seq_bases = "ACGTA"
  )
  expect_equal(result, 100L:104L)
})

test_that("cigarToRefPos: CIGAR with deletion", {
  # 3M2D2M: 3 matches, 2 deleted ref bases, 2 more matches
  # read: ACGTA (5 bases), ref positions: 100,101,102, skip 103+104, 105,106
  result <- commaKit:::.cigarToRefPos(
    "3M2D2M",
    ref_start = 100L,
    seq_bases = "ACGTA"
  )
  expect_equal(result, c(100L, 101L, 102L, 105L, 106L))
})

test_that("cigarToRefPos: CIGAR with insertion", {
  # 3M2I2M: 3 matches, 2 inserted read bases (no ref pos), 2 more matches
  # read: ACGTAA (7 bases), ref positions: 100,101,102, NA,NA, 103,104
  result <- commaKit:::.cigarToRefPos(
    "3M2I2M",
    ref_start = 100L,
    seq_bases = "ACGTAAG"
  )
  expect_equal(
    result,
    c(100L, 101L, 102L, NA_integer_, NA_integer_, 103L, 104L)
  )
})

test_that("cigarToRefPos: CIGAR with soft clip", {
  # 2S3M: 2 soft-clipped bases (NA ref pos), 3 matches
  result <- commaKit:::.cigarToRefPos(
    "2S3M",
    ref_start = 50L,
    seq_bases = "XXACG"
  )
  expect_equal(result, c(NA_integer_, NA_integer_, 50L, 51L, 52L))
})

test_that("cigarToRefPos: returns NULL for empty cigar", {
  result <- commaKit:::.cigarToRefPos("", ref_start = 1L, seq_bases = "A")
  expect_null(result)
})

test_that("cigarToRefPos: malformed CIGAR with trailing junk returns NULL", {
  result <- commaKit:::.cigarToRefPos(
    "5Mbogus",
    ref_start = 1L,
    seq_bases = "ACGTA"
  )
  expect_null(result)
})

test_that("cigarToRefPos: overlong CIGAR stops after exhausting read", {
  result <- commaKit:::.cigarToRefPos(
    "5M1M",
    ref_start = 10L,
    seq_bases = "ACGTA"
  )
  expect_equal(result, 10L:14L)
})

test_that("cigarToRefPos: matches GenomicAlignments reference/query ranges", {
  skip_if_not_installed("GenomicAlignments")
  skip_if_not_installed("IRanges")

  cases <- list(
    list(cigar = "5M", ref_start = 100L, seq_bases = "ACGTA"),
    list(cigar = "3M2D2M", ref_start = 100L, seq_bases = "ACGTA"),
    list(cigar = "3M2I2M", ref_start = 100L, seq_bases = "ACGTAAG"),
    list(cigar = "2S3M", ref_start = 50L, seq_bases = "XXACG"),
    list(cigar = "2M3N2M", ref_start = 100L, seq_bases = "ACGT")
  )

  for (case in cases) {
    observed <- commaKit:::.cigarToRefPos(
      case$cigar,
      ref_start = case$ref_start,
      seq_bases = case$seq_bases
    )

    query_ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(
      case$cigar,
      ops = c("M", "=", "X"),
      drop.empty.ranges = FALSE
    )[[1L]]
    reference_ranges <- GenomicAlignments::cigarRangesAlongReferenceSpace(
      case$cigar,
      pos = case$ref_start,
      ops = c("M", "=", "X"),
      drop.empty.ranges = FALSE
    )[[1L]]

    expected <- rep(NA_integer_, nchar(case$seq_bases))
    for (i in seq_along(query_ranges)) {
      query_pos <- seq(
        IRanges::start(query_ranges)[[i]],
        IRanges::end(query_ranges)[[i]]
      )
      reference_pos <- seq(
        IRanges::start(reference_ranges)[[i]],
        IRanges::end(reference_ranges)[[i]]
      )
      expected[query_pos] <- reference_pos[seq_along(query_pos)]
    }

    expect_equal(observed, expected, info = case$cigar)
  }
})

# ─── .parseMmTag() ────────────────────────────────────────────────────────────

test_that("parseMmTag: parses single-mod-type MM tag", {
  # MM: "A+a?,0" means first A in read is modified (6mA)
  # ML: 200 → probability 200/255 ≈ 0.78 → is_mod = TRUE
  result <- commaKit:::.parseMmTag(
    mm_tag = "A+a?,0",
    ml_tag = as.raw(200L),
    seq_bases = "ACGT"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_equal(result$mod_type, "6mA")
  expect_equal(result$read_pos, 1L) # first A at position 1
  expect_true(result$is_mod)
})

test_that("parseMmTag: ML probability <= 0.5 → is_mod = FALSE", {
  # ML: 100 → 100/255 ≈ 0.39 < 0.5 → not modified
  result <- commaKit:::.parseMmTag(
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
  result <- commaKit:::.parseMmTag(
    mm_tag    = "A+a?,0,1",
    ml_tag    = as.raw(c(200L, 210L)),
    seq_bases = "AACGTA"
  )
  expect_equal(nrow(result), 2L)
  expect_equal(result$read_pos, c(1L, 6L))
})

test_that("parseMmTag: returns NULL for null tags", {
  expect_null(commaKit:::.parseMmTag(NULL, as.raw(200L), "ACGT"))
  expect_null(commaKit:::.parseMmTag("A+a?,0", NULL, "ACGT"))
})

test_that("parseMmTag: skips unknown mod_code", {
  # mod_code "z" is not in .MODKIT_CODE_MAP
  result <- commaKit:::.parseMmTag(
    mm_tag    = "A+z?,0",
    ml_tag    = as.raw(200L),
    seq_bases = "ACGT"
  )
  expect_null(result)
})

# ─── .parseDorado() error handling ────────────────────────────────────────────

test_that("parseDorado: error on missing file", {
  expect_error(
    commaKit:::.parseDorado("/nonexistent/path.bam", "sample1"),
    "not found"
  )
})

test_that("parseDorado: error on non-character file argument", {
  expect_error(
    commaKit:::.parseDorado(123L, "sample1"),
    "must be a single character string"
  )
})

test_that("parseDorado: error on length-2 file argument", {
  expect_error(
    commaKit:::.parseDorado(c("a.bam", "b.bam"), "sample1"),
    "must be a single character string"
  )
})

test_that("parseDorado: synthetic BAM reaches site aggregation path", {
  bam_file <- .make_dorado_test_bam(c(
    .make_sam_record(
      qname = "read_mod",
      pos = 100L,
      cigar = "4M",
      seq = "ACGT",
      mm = "A+a?,0;C+m?,0;",
      ml = c(200L, 100L)
    ),
    .make_sam_record(
      qname = "read_canonical",
      pos = 100L,
      cigar = "4M",
      seq = "ACGT",
      mm = "A+a?,0;",
      ml = 100L
    )
  ))

  result <- commaKit:::.parseDorado(
    bam_file,
    sample_name = "sample1",
    min_coverage = 1L
  )
  result <- result[order(result$position, result$mod_type), ]

  expect_s3_class(result, "data.frame")
  expect_named(result, c(
    "chrom", "position", "strand", "mod_type", "motif", "beta", "coverage",
    "mod_counts", "canonical_counts", "other_mod_counts"
  ))
  expect_equal(result$chrom, c("chr1", "chr1"))
  expect_equal(result$position, c(100L, 101L))
  expect_equal(result$strand, c("+", "+"))
  expect_equal(result$mod_type, c("6mA", "5mC"))
  expect_true(all(is.na(result$motif)))
  expect_equal(result$beta, c(0.5, 0))
  expect_equal(result$coverage, c(2L, 1L))
  expect_equal(result$mod_counts, c(1L, 0L))
  expect_equal(result$canonical_counts, c(1L, 1L))
  expect_true(all(is.na(result$other_mod_counts)))
})

test_that("parseDorado: malformed MM/ML reads are skipped without recycling", {
  bam_file <- .make_dorado_test_bam(c(
    .make_sam_record(
      qname = "read_valid",
      pos = 100L,
      cigar = "4M",
      seq = "ACGT",
      mm = "A+a?,0;",
      ml = 200L
    ),
    .make_sam_record(
      qname = "read_truncated",
      pos = 100L,
      cigar = "4M",
      seq = "AAGT",
      mm = "A+a?,0,0;",
      ml = 200L
    )
  ))

  result <- commaKit:::.parseDorado(
    bam_file,
    sample_name = "sample1",
    min_coverage = 1L
  )

  expect_equal(nrow(result), 1L)
  expect_equal(result$position, 100L)
  expect_equal(result$coverage, 1L)
  expect_equal(result$mod_counts, 1L)
  expect_equal(result$canonical_counts, 0L)
})

test_that("parseDorado: CIGAR-unmapped modification calls are dropped", {
  bam_file <- .make_dorado_test_bam(c(
    .make_sam_record(
      qname = "read_inserted_call",
      pos = 100L,
      cigar = "1M1I2M",
      seq = "AAAA",
      mm = "A+a?,1;",
      ml = 200L
    ),
    .make_sam_record(
      qname = "read_on_reference",
      pos = 100L,
      cigar = "4M",
      seq = "AAAA",
      mm = "A+a?,0;",
      ml = 200L
    )
  ))

  result <- commaKit:::.parseDorado(
    bam_file,
    sample_name = "sample1",
    min_coverage = 1L
  )

  expect_equal(nrow(result), 1L)
  expect_equal(result$position, 100L)
  expect_equal(result$coverage, 1L)
  expect_equal(result$mod_counts, 1L)
})

# ─── .cigarToRefPos() additional edge cases ───────────────────────────────────

test_that("cigarToRefPos: hard clip (H) does not consume read positions", {
  # 2H5M: 2 hard-clipped bases (not in seq), then 5 matches
  # Hard clips are NOT in the sequence string, so seq has 5 bases
  result <- commaKit:::.cigarToRefPos(
    "2H5M",
    ref_start = 10L,
    seq_bases = "ACGTA"
  )
  # All 5 read positions map to ref 10..14
  expect_equal(result, 10L:14L)
})

test_that("cigarToRefPos: all soft-clip returns all NA", {
  # 5S: entire read is soft-clipped; no ref positions assigned
  result <- commaKit:::.cigarToRefPos("5S", ref_start = 1L, seq_bases = "ACGTA")
  expect_equal(result, rep(NA_integer_, 5L))
})

test_that(
  "cigarToRefPos: mixed soft-clip and match assigns NA then ref positions",
  {
    # 3S2M: 3 soft-clipped (NA), then 2 matches starting at ref 50
    result <- commaKit:::.cigarToRefPos(
      "3S2M",
      ref_start = 50L,
      seq_bases = "XXXAC"
    )
    expect_equal(result, c(NA_integer_, NA_integer_, NA_integer_, 50L, 51L))
  }
)

test_that(
  "cigarToRefPos: N (skipped region) advances reference like deletion",
  {
    # 2M3N2M: 2 matches, skip 3 ref bases (intron), 2 more matches
    result <- commaKit:::.cigarToRefPos(
      "2M3N2M",
      ref_start = 100L,
      seq_bases = "ACGT"
    )
    expect_equal(result, c(100L, 101L, 105L, 106L))
  }
)

# ─── .parseMmTag() additional tests ──────────────────────────────────────────

test_that("parseMmTag: 5mC (mod_code 'm') is recognised", {
  # MM: "C+m?,0" — first C in read is potentially 5-methylcytosine
  result <- commaKit:::.parseMmTag(
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
  result <- commaKit:::.parseMmTag(
    mm_tag    = "A+a?,0;C+m?,0",
    ml_tag    = as.raw(c(200L, 200L)),
    seq_bases = "ACGT"
  )
  expect_equal(nrow(result), 2L)
  expect_true("6mA" %in% result$mod_type)
  expect_true("5mC" %in% result$mod_type)
})

test_that(
  "parseMmTag: zero modifications encoded (empty delta list) returns NULL",
  {
    # MM: "A+a?" with no deltas means zero modifications for this type
    result <- commaKit:::.parseMmTag(
      mm_tag    = "A+a?",
      ml_tag    = as.raw(integer(0)),
      seq_bases = "ACGT"
    )
    # No deltas parsed → NULL
    expect_null(result)
  }
)

test_that(
  "parseMmTag: ML at exactly boundary (128/255 ≈ 0.502) is_mod = TRUE",
  {
    # 128/255 ≈ 0.502 > 0.5 → is_mod TRUE
    result <- commaKit:::.parseMmTag(
      mm_tag    = "A+a?,0",
      ml_tag    = as.raw(128L),
      seq_bases = "ACGT"
    )
    expect_true(result$is_mod[1])
  }
)

test_that("parseMmTag: ML at exactly 127/255 ≈ 0.498 is_mod = FALSE", {
  # 127/255 ≈ 0.498 ≤ 0.5 → is_mod FALSE
  result <- commaKit:::.parseMmTag(
    mm_tag    = "A+a?,0",
    ml_tag    = as.raw(127L),
    seq_bases = "ACGT"
  )
  expect_false(result$is_mod[1])
})

# ─────────────────────────────────────────────────────────────────────────────
# Production-like Dorado MM/ML/CIGAR edge cases
# ─────────────────────────────────────────────────────────────────────────────

test_that("cigarToRefPos: malformed CIGAR returns NULL rather than partial map", {
  expect_null(commaKit:::.cigarToRefPos(
    "3M2Z2M",
    ref_start = 100L,
    seq_bases = "ACGTACG"
  ))
})

test_that("Dorado MM/ML helper path drops modification calls on insertions", {
  ref_positions <- commaKit:::.cigarToRefPos(
    "1M1I2M",
    ref_start = 100L,
    seq_bases = "AAAA"
  )
  calls <- commaKit:::.parseMmTag(
    mm_tag = "A+a?,1",
    ml_tag = as.raw(200L),
    seq_bases = "AAAA"
  )

  expect_equal(calls$read_pos, 2L)
  expect_true(is.na(ref_positions[calls$read_pos]))
})

test_that("parseMmTag: truncated ML array returns NULL without recycling", {
  result <- commaKit:::.parseMmTag(
    mm_tag = "A+a?,0,0",
    ml_tag = as.raw(200L),
    seq_bases = "AAAA"
  )
  expect_null(result)
})
