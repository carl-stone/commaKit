.sourceDocPath <- function(...) {
    rel <- file.path(...)
    candidates <- c(
        testthat::test_path("..", "..", rel),
        testthat::test_path("..", "..", "00_pkg_src", "commaKit", rel),
        testthat::test_path("..", "..", "..", "00_pkg_src", "commaKit", rel)
    )
    found <- candidates[file.exists(candidates)]
    if (length(found)) found[1] else NA_character_
}

test_that("public docs use commaKit load and install instructions", {
    readme <- .sourceDocPath("README.Rmd")
    if (is.na(readme)) {
        readme <- .sourceDocPath("README.md")
    }

    docs <- c(
        readme,
        .sourceDocPath("vignettes", "getting-started.Rmd"),
        .sourceDocPath("vignettes", "understanding-commaData.Rmd"),
        .sourceDocPath("vignettes", "multiple-modification-types.Rmd"),
        .sourceDocPath("vignettes", "import-troubleshooting.Rmd")
    )

    docs <- unname(docs)
    present <- file.exists(docs)
    expect_true(all(present))

    text <- paste(unlist(lapply(docs[present], readLines, warn = FALSE)), collapse = "\n")

    expect_match(text, 'devtools::install_github\\("carl-stone/commaKit"\\)')
    expect_match(text, "library\\(commaKit\\)")
    expect_no_match(text, "library\\(comma\\)")
    expect_no_match(text, 'package = "comma"')
    expect_no_match(text, 'BiocManager::install\\("comma"\\)')
    expect_no_match(text, "carl-stone/comma($|[^K])")
})
