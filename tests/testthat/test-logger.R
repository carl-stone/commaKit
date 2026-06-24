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
  on.exit({
    try(if (isOpen(con)) close(con), silent = TRUE)
    options(old)
  })

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
  on.exit({
    try(if (isOpen(con)) close(con), silent = TRUE)
    options(old)
  })

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

test_that("structured logger default stream is separate from messages", {
  old <- options(commaKit.log = TRUE)
  options(commaKit.log.connection = NULL)
  on.exit(options(old))

  messages <- NULL
  logged <- FALSE
  output <- capture.output(
    messages <- capture.output(
      {
        logged <- commaKit:::.comma_log_event(
          "default_stream",
          component = "test",
          .time = as.POSIXct("2026-01-02 03:04:05", tz = "UTC")
        )
        message("human-readable progress")
      },
      type = "message"
    ),
    type = "output"
  )

  expect_true(logged)
  expect_match(output, '"event":"default_stream"', fixed = TRUE)
  expect_equal(messages, "human-readable progress")
})

test_that("commaData constructor records structured lifecycle events", {
  out <- character()
  con <- textConnection("out", "w", local = TRUE)
  old <- options(commaKit.log = TRUE, commaKit.log.connection = con)
  bed_file <- .write_logger_modkit()
  on.exit({
    try(if (isOpen(con)) close(con), silent = TRUE)
    options(old)
    unlink(bed_file)
  })

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

test_that("structured logger uses JSON-safe numeric and string encoding", {
  out <- character()
  con <- textConnection("out", "w", local = TRUE)
  old <- options(
    commaKit.log = TRUE,
    commaKit.log.connection = con,
    OutDec = ","
  )
  on.exit({
    try(if (isOpen(con)) close(con), silent = TRUE)
    options(old)
  })

  logged <- commaKit:::.comma_log_event(
    "encoding_check",
    component = "test",
    ratio = 1.25,
    note = paste0("back", "\b", "form", "\f", "unit", "\001"),
    .time = as.POSIXct("2026-01-02 03:04:05", tz = "UTC")
  )

  close(con)
  options(old)
  expect_true(logged)
  expect_length(out, 1L)
  expect_match(out, '"ratio":1.25', fixed = TRUE)
  expect_match(out, '"note":"back\\bform\\funit\\u0001"', fixed = TRUE)
})

test_that("structured logger failures do not interrupt callers", {
  out <- character()
  con <- textConnection("out", "w", local = TRUE)
  old <- options(commaKit.log = TRUE, commaKit.log.connection = con)
  on.exit({
    try(if (isOpen(con)) close(con), silent = TRUE)
    options(old)
  })

  logged <- expect_no_error(commaKit:::.comma_log_event(
    "bad_fields",
    level = "info",
    component = NULL,
    list(unnamed = TRUE)
  ))

  close(con)
  options(old)
  expect_false(logged)
  expect_equal(out, character())
})

test_that("contextual error tracking captures breadcrumbs and user context", {
  captured <- NULL
  old <- options(
    commaKit.error_tracking = TRUE,
    commaKit.error_tracking.breadcrumbs = list(),
    commaKit.error_tracking.user = list(
      id = "agent-1",
      email = "agent@example.test"
    ),
    commaKit.error_tracking.release = "commaKit@test",
    commaKit.error_tracking.environment = "test",
    commaKit.error_tracking.reporter = function(event) {
      captured <<- event
      TRUE
    },
    commaKit.sentry.dsn = ""
  )
  on.exit(options(old))

  commaKit:::.comma_log_event(
    "sample_parse_started",
    component = "commaData",
    caller = "modkit",
    sample_name = "s1",
    .time = as.POSIXct("2026-01-02 03:04:05", tz = "UTC")
  )
  tracked <- commaKit:::.comma_track_error(
    simpleError("synthetic failure"),
    component = "commaData",
    operation = "sample_parse",
    caller = "modkit",
    sample_name = "s1",
    .time = as.POSIXct("2026-01-02 03:04:06", tz = "UTC")
  )

  expect_true(tracked)
  expect_type(captured, "list")
  expect_equal(captured$logger, "commaKit")
  expect_equal(captured$release, "commaKit@test")
  expect_equal(captured$environment, "test")
  expect_equal(captured$user$id, "agent-1")
  expect_equal(captured$tags$component, "commaData")
  expect_equal(captured$tags$operation, "sample_parse")
  expect_equal(captured$extra$caller, "modkit")
  expect_equal(captured$extra$sample_name, "s1")
  expect_equal(
    captured$exception$values[[1L]]$value,
    "synthetic failure"
  )
  expect_true(length(captured$exception$values[[1L]]$stacktrace$frames) > 0L)
  messages <- vapply(
    captured$breadcrumbs$values,
    `[[`,
    character(1),
    "message"
  )
  expect_true("sample_parse_started" %in% messages)
})

test_that("commaData parse errors are reported with sample context", {
  captured <- NULL
  missing_file <- tempfile(fileext = ".bed")
  old <- options(
    commaKit.error_tracking = TRUE,
    commaKit.error_tracking.breadcrumbs = list(),
    commaKit.error_tracking.reporter = function(event) {
      captured <<- event
      TRUE
    },
    commaKit.sentry.dsn = ""
  )
  on.exit(options(old))

  expect_error(
    suppressMessages(commaData(
      files = c(s1 = missing_file),
      colData = data.frame(sample_name = "s1", replicate = 1L),
      genome = c(chr_sim = 1000L),
      caller = "modkit",
      min_coverage = 1L
    )),
    "file not found|not found"
  )

  expect_type(captured, "list")
  expect_equal(captured$transaction, "sample_parse")
  expect_equal(captured$extra$caller, "modkit")
  expect_equal(captured$extra$sample_name, "s1")
  expect_equal(captured$extra$file, missing_file)
  messages <- vapply(
    captured$breadcrumbs$values,
    `[[`,
    character(1),
    "message"
  )
  expect_true("commaData_construct_started" %in% messages)
  expect_true("sample_parse_started" %in% messages)
})

test_that("Sentry DSN is mapped to the envelope endpoint", {
  endpoint <- commaKit:::.comma_sentry_envelope_url(
    "https://public@example.sentry.io/12345"
  )

  expect_equal(endpoint, "https://example.sentry.io/api/12345/envelope/")
})
