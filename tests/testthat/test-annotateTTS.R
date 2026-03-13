# annotateTTS is an internal (non-exported) function.
# Access it via CoMMA:::annotateTTS.

test_that("annotatedOnly = FALSE: sense strand site annotated using TU Right boundary", {
  # tuA: Right=2000, Strand="+". Window size 200.
  # Position 1900 is within [2000-200, 2000+200] = [1800, 2200]. ✓
  # RelPos = 1900 - 2000 = -100
  methyl <- data.frame(Position = 1900L, beta = 0.5, stringsAsFactors = FALSE)
  result <- CoMMA:::annotateTTS(methyl, fixture_tu_df, location = "Position",
                                size = 200, long = FALSE, annotatedOnly = FALSE)
  expect_true("RelPos_+1" %in% names(result))
  expect_equal(result$`RelPos_+1`, -100L)
})

test_that("annotatedOnly = FALSE: antisense strand site annotated using TU Left boundary", {
  # tuB: Left=2500, Strand="-". Window size 200.
  # Position 2600 is within [2500-200, 2500+200] = [2300, 2700]. ✓
  # RelPos = 2500 - 2600 = -100
  methyl <- data.frame(Position = 2600L, beta = 0.5, stringsAsFactors = FALSE)
  result <- CoMMA:::annotateTTS(methyl, fixture_tu_df, location = "Position",
                                size = 200, long = FALSE, annotatedOnly = FALSE)
  expect_true("RelPos_-1" %in% names(result))
  expect_equal(result$`RelPos_-1`, -100L)
})

test_that("annotatedOnly = FALSE: position outside all TTS windows gets NoTTS = 'X'", {
  methyl <- data.frame(Position = 10L, beta = 0.5, stringsAsFactors = FALSE)
  result <- CoMMA:::annotateTTS(methyl, fixture_tu_df, location = "Position",
                                size = 200, long = FALSE, annotatedOnly = FALSE)
  expect_equal(result$NoTTS, "X")
})

test_that("annotatedOnly = TRUE: uses terminator features, not TU boundaries", {
  # fixture_term_df contains a Rho-Independent-Terminators at Left=2050, Right=2150, Strand="+"
  # Position 2100 is within [2050-200, 2050+200] = [1850, 2250]. ✓
  # RelPos = 2100 - 2050 = 50
  methyl <- data.frame(Position = 2100L, beta = 0.5, stringsAsFactors = FALSE)
  result <- CoMMA:::annotateTTS(methyl, fixture_term_df, location = "Position",
                                size = 200, long = FALSE, annotatedOnly = TRUE)
  expect_true("RelPos_+1" %in% names(result))
  expect_equal(result$`RelPos_+1`, 50L)
})

test_that("annotatedOnly = TRUE: TU rows in meta_df are ignored", {
  # Position 1900 is near tuA TTS (Right=2000) but annotatedOnly=TRUE should
  # only look at terminator features, not Transcription-Units.
  methyl <- data.frame(Position = 1900L, beta = 0.5, stringsAsFactors = FALSE)
  result <- CoMMA:::annotateTTS(methyl, fixture_term_df, location = "Position",
                                size = 50, long = FALSE, annotatedOnly = TRUE)
  # Only termA at [2050,2150] is considered; 1900 is outside its window → NoTTS
  expect_equal(result$NoTTS, "X")
})

test_that("long = TRUE returns RelPos column with no NA values", {
  methyl <- data.frame(Position = 1900L, beta = 0.5, stringsAsFactors = FALSE)
  result <- CoMMA:::annotateTTS(methyl, fixture_tu_df, location = "Position",
                                size = 200, long = TRUE, annotatedOnly = FALSE)
  expect_true("RelPos" %in% names(result))
  expect_false(any(is.na(result$RelPos)))
})

test_that("long = FALSE preserves input row count", {
  result <- CoMMA:::annotateTTS(fixture_methyl_df, fixture_tu_df, location = "Position",
                                size = 500, long = FALSE, annotatedOnly = FALSE)
  expect_equal(nrow(result), nrow(fixture_methyl_df))
})
