test_that("R Markdown checks ignore sources without Markdown outputs", {
  skip_if_not_installed("rmarkdown")

  script <- testthat::test_path("..", "..", "dev", "precommit.R")
  vignette <- testthat::test_path(
    "..",
    "..",
    "vignettes",
    "import-troubleshooting.Rmd"
  )
  skip_if_not(file.exists(script), "pre-commit script is unavailable")
  skip_if_not(file.exists(vignette), "source vignette is unavailable")

  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", script, "rmarkdown", vignette),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }

  expect_equal(status, 0L, info = paste(output, collapse = "\n"))
})
