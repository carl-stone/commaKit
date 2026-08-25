.commaDataVignettePath <- function() {
  vignette_name <- "understanding-commaData.Rmd"
  candidates <- c(
    file.path("vignettes", vignette_name),
    testthat::test_path(
      "..",
      "..",
      "vignettes",
      "understanding-commaData.Rmd"
    ),
    system.file("doc", vignette_name, package = "commaKit", mustWork = FALSE)
  )
  candidates[file.exists(candidates)][1]
}

test_that("commaData vignette knits with package-provided data", {
  vignette_file <- .commaDataVignettePath()
  skip_if(
    is.na(vignette_file),
    "understanding-commaData vignette source is not available in this build"
  )

  output_file <- tempfile(fileext = ".md")
  expect_no_error(
    knitr::knit(
      input = vignette_file,
      output = output_file,
      quiet = TRUE,
      envir = new.env(parent = globalenv())
    )
  )

  expect_true(file.exists(output_file))
})
