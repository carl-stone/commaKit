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
      max_bytes = 10
    ),
    "Files exceed"
  )
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
