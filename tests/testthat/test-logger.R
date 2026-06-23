library(testthat)

.write_logger_modkit <- function(file = tempfile(fileext = ".bed")) {
  rows <- data.frame(
    chrom = "chr_sim",
    start = 99L,
    end = 100L,
    mod_code = "a,GATC,1",
    score = 20L,
    strand = "+",
    thickStart = 99L,
    thickEnd = 100L,
    itemRgb = "255,0,0",
    Nvalid_cov = 20L,
    fraction_modified = 50,
    Nmod = 10L,
    Ncanonical = 7L,
    Nother_mod = 3L,
    Ndelete = 0L,
    Nfail = 0L,
    Ndiff = 0L,
    Nnocall = 0L,
    stringsAsFactors = FALSE
  )
  write.table(rows,
    file = file, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  file
}

test_that("structured logger is silent unless enabled", {
  out <- character()
  con <- textConnection("out", "w", local = TRUE)
  old <- options(commaKit.log = FALSE, commaKit.log.connection = con)

  logged <- commaKit:::.comma_log_event(
    "suppressed_event",
    component = "test",
    .time = as.POSIXct("2026-01-02 03:04:05", tz = "UTC")
  )

  close(con)
  options(old)
  expect_false(logged)
  expect_equal(out, character())
})

test_that("structured logger emits machine-readable JSON fields", {
  out <- character()
  con <- textConnection("out", "w", local = TRUE)
  old <- options(commaKit.log = TRUE, commaKit.log.connection = con)

  logged <- commaKit:::.comma_log_event(
    "sample_parse_finished",
    component = "parser",
    sample_name = "s1",
    site_count = 2L,
    labels = c("6mA", "5mC"),
    note = "quoted \"value\"",
    ok = TRUE,
    .time = as.POSIXct("2026-01-02 03:04:05", tz = "UTC")
  )

  close(con)
  options(old)
  expect_true(logged)
  expect_length(out, 1L)
  expect_match(out, '"timestamp":"2026-01-02T03:04:05.000Z"', fixed = TRUE)
  expect_match(out, '"level":"info"', fixed = TRUE)
  expect_match(out, '"event":"sample_parse_finished"', fixed = TRUE)
  expect_match(out, '"component":"parser"', fixed = TRUE)
  expect_match(out, '"sample_name":"s1"', fixed = TRUE)
  expect_match(out, '"site_count":2', fixed = TRUE)
  expect_match(out, '"labels":["6mA","5mC"]', fixed = TRUE)
  expect_match(out, '"note":"quoted \\"value\\""', fixed = TRUE)
  expect_match(out, '"ok":true', fixed = TRUE)
})

test_that("commaData constructor records structured lifecycle events", {
  out <- character()
  con <- textConnection("out", "w", local = TRUE)
  old <- options(commaKit.log = TRUE, commaKit.log.connection = con)
  bed_file <- .write_logger_modkit()

  obj <- suppressMessages(commaData(
    files = c(s1 = bed_file),
    colData = data.frame(sample_name = "s1", replicate = 1L),
    genome = c(chr_sim = 1000L),
    caller = "modkit",
    min_coverage = 1L
  ))

  close(con)
  options(old)
  log_text <- paste(out, collapse = "\n")
  expect_true(is(obj, "commaData"))
  expect_match(log_text, '"event":"commaData_construct_started"', fixed = TRUE)
  expect_match(log_text, '"event":"sample_parse_started"', fixed = TRUE)
  expect_match(log_text, '"event":"sample_parse_finished"', fixed = TRUE)
  expect_match(log_text, '"event":"commaData_construct_finished"', fixed = TRUE)
  expect_match(log_text, '"sample_count":1', fixed = TRUE)
  expect_match(log_text, '"site_count":1', fixed = TRUE)
  expect_match(log_text, '"caller":"modkit"', fixed = TRUE)
})
