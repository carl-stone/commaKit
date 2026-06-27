readiness_env <- new.env(parent = baseenv())
readiness_script <- system.file(
  "scripts",
  "readiness-checks.R",
  package = "commaKit"
)
if (!nzchar(readiness_script)) {
  readiness_script <- testthat::test_path(
    "..",
    "..",
    "inst",
    "scripts",
    "readiness-checks.R"
  )
}
sys.source(readiness_script, envir = readiness_env)

test_that("large-file readiness check enforces byte limit", {
  root <- tempfile("readiness-root-")
  dir.create(root)
  writeLines(strrep("x", 12), file.path(root, "too-large.txt"))

  expect_error(
    readiness_env$readiness_check_large_files(
      files = "too-large.txt",
      root = root,
      max_bytes = NA_real_
    ),
    "single positive numeric"
  )

  expect_error(
    readiness_env$readiness_check_large_files(
      files = "too-large.txt",
      root = root,
      max_bytes = 10
    ),
    "Files exceed"
  )
})

test_that("readiness_run can limit checks to explicit files", {
  root <- tempfile("readiness-root-")
  dir.create(root)
  marker <- paste0("TO", "DO")
  writeLines(
    paste(marker, "unrelated scratch", sep = ": "),
    file.path(root, "scratch.R")
  )
  writeLines("ok <- TRUE", file.path(root, "staged.R"))
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)

  expect_true(readiness_env$readiness_run(
    "debt-markers",
    files = "staged.R",
    root = root
  ))
})

test_that("debt-marker readiness check requires issue references", {
  root <- tempfile("readiness-root-")
  dir.create(root)
  marker <- paste0("TO", "DO")
  writeLines(
    paste(marker, "explain this later", sep = ": "),
    file.path(root, "bad.R")
  )
  writeLines(paste0(marker, "(#123): tracked work"), file.path(root, "good.R"))

  expect_error(
    readiness_env$readiness_check_debt_markers(
      files = c("bad.R", "good.R"),
      root = root
    ),
    "Debt markers must reference"
  )

  expect_true(readiness_env$readiness_check_debt_markers(
    files = "good.R",
    root = root
  ))
})

test_that("AGENTS link readiness check validates local links", {
  root <- tempfile("readiness-root-")
  dir.create(file.path(root, "docs"), recursive = TRUE)
  writeLines("ok", file.path(root, "docs", "guide.md"))
  writeLines(
    c("[Guide](docs/guide.md)", "[External](https://example.com/docs)"),
    file.path(root, "AGENTS.md")
  )

  expect_true(readiness_env$readiness_check_agents_links(
    files = "AGENTS.md",
    root = root
  ))

  writeLines("[Missing](docs/missing.md)", file.path(root, "AGENTS.md"))
  expect_error(
    readiness_env$readiness_check_agents_links(
      files = "AGENTS.md",
      root = root
    ),
    "missing local links"
  )
})

test_that("complexity readiness check enforces branch thresholds", {
  root <- tempfile("readiness-root-")
  dir.create(file.path(root, "R"), recursive = TRUE)
  writeLines(
    c(
      "simple <- function(x) {",
      "  if (x > 1) x else 0",
      "}",
      "complex <- function(x) {",
      "  if (x > 1) x <- x + 1",
      "  if (x > 2) x <- x + 1",
      "  if (x > 3) x <- x + 1",
      "  x",
      "}"
    ),
    file.path(root, "R", "quality.R")
  )

  expect_true(readiness_env$readiness_check_complexity(
    files = "R/quality.R",
    root = root,
    max_complexity = 4L
  ))
  expect_error(
    readiness_env$readiness_check_complexity(
      files = "R/quality.R",
      root = root,
      max_complexity = 3L
    ),
    "cyclomatic complexity"
  )
})

test_that("dead-code readiness check flags unused internal helpers", {
  root <- tempfile("readiness-root-")
  dir.create(file.path(root, "R"), recursive = TRUE)
  writeLines(
    c(
      ".used_helper <- function(x) x + 1",
      ".unused_helper <- function(x) x - 1",
      "public <- function(x) .used_helper(x)"
    ),
    file.path(root, "R", "helpers.R")
  )

  expect_error(
    readiness_env$readiness_check_dead_code(
      files = "R/helpers.R",
      root = root
    ),
    ".unused_helper",
    fixed = TRUE
  )
})

test_that("duplicate-code readiness check compares normalized windows", {
  root <- tempfile("readiness-root-")
  dir.create(file.path(root, "R"), recursive = TRUE)
  repeated <- c(
    "alpha <- 1",
    "beta <- 2",
    "gamma <- alpha + beta"
  )
  writeLines(
    c(
      "one <- function() {",
      repeated,
      "}",
      "two <- function() {",
      repeated,
      "}"
    ),
    file.path(root, "R", "duplicate.R")
  )

  expect_error(
    readiness_env$readiness_check_duplicate_code(
      files = "R/duplicate.R",
      root = root,
      window_size = 3L
    ),
    "Duplicate code windows"
  )
})
