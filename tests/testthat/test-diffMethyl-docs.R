.diffMethylDocPath <- function(...) {
  relative_path <- file.path(...)
  source_dir <- Sys.getenv("COMMAKIT_SOURCE_DIR", unset = "")
  candidates <- c(
    if (nzchar(source_dir)) {
      file.path(source_dir, relative_path)
    } else {
      character()
    },
    testthat::test_path("..", "..", relative_path),
    testthat::test_path("..", "..", "00_pkg_src", "commaKit", relative_path)
  )
  found <- candidates[file.exists(candidates)]
  if (length(found)) found[1] else NA_character_
}

.renderDiffMethylRd <- function(path) {
  rendered_path <- tempfile(fileext = ".txt")
  on.exit(unlink(rendered_path), add = TRUE)

  tools::Rd2txt(
    tools::parse_Rd(path),
    out = rendered_path
  )
  paste(readLines(rendered_path, warn = FALSE), collapse = "\n")
}

test_that(
  "docs-lint: diffMethyl preserves method and limitation guidance",
  {
    source_path <- .diffMethylDocPath("R", "diffMethyl.R")
    skip_if(
      is.na(source_path),
      "diffMethyl source is unavailable in this build"
    )
    source <- paste(readLines(source_path, warn = FALSE), collapse = "\n")

    # These are high-value user guidance, not implementation contracts.
    expect_match(source, "Default method")
    expect_match(source, "quasi_f")
    expect_match(source, "general-purpose alternative")
    expect_match(source, "Small-sample note")
    expect_match(source, "extremely low statistical power")
    expect_match(source, "exploratory only")
  }
)

test_that("rendered docs: diffMethyl Rd keeps method guidance", {
  rd_path <- .diffMethylDocPath("man", "diffMethyl.Rd")
  skip_if(
    is.na(rd_path),
    "diffMethyl Rd is unavailable in this build"
  )
  rendered <- .renderDiffMethylRd(rd_path)

  expect_match(rendered, "Default method")
  expect_match(rendered, "methylkit")
  expect_match(rendered, "quasi_f")
  expect_match(rendered, "general-purpose alternative")
  expect_match(rendered, "Small-sample note")
  expect_match(rendered, "exploratory only")
})

test_that("rendered vignette: diffMethyl keeps method guidance", {
  vignette_path <- .diffMethylDocPath(
    "vignettes",
    "getting-started.Rmd"
  )
  skip_if(
    is.na(vignette_path),
    "getting-started vignette source is unavailable in this build"
  )
  skip_if_not_installed("knitr")

  output_path <- tempfile(fileext = ".md")
  on.exit(unlink(output_path), add = TRUE)
  expect_no_error(
    knitr::knit(
      input = vignette_path,
      output = output_path,
      quiet = TRUE,
      envir = new.env(parent = globalenv())
    )
  )
  rendered <- paste(readLines(output_path, warn = FALSE), collapse = "\n")

  expect_match(rendered, 'method = "methylkit"')
  expect_match(rendered, "default")
  expect_match(rendered, 'method = "quasi_f"')
  expect_match(rendered, "general-purpose alternative")
  expect_match(rendered, 'method = "limma"')
  expect_match(rendered, "M-values")
})
