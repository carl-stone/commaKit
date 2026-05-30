.commaDataVignettePath <- function() {
    vignette_name <- "understanding-commaData.Rmd"
    candidates <- c(
        file.path("vignettes", vignette_name),
        testthat::test_path(
            "..", "..", "vignettes", "understanding-commaData.Rmd"
        ),
        system.file("doc", vignette_name, package = "comma", mustWork = FALSE)
    )
    candidates[file.exists(candidates)][1]
}

test_that("commaData vignette source has expected teaching contract", {
    vignette_file <- .commaDataVignettePath()

    expect_true(file.exists(vignette_file))

    source <- readLines(vignette_file, warn = FALSE)
    source_text <- paste(source, collapse = "\n")

    expect_match(source_text, "VignetteIndexEntry\\{Understanding the commaData object\\}")
    expect_match(source_text, "VignetteEngine\\{knitr::rmarkdown\\}")
    expect_match(source_text, "data\\(comma_example_data\\)")
    expect_match(source_text, "mod_context.*modification \\+ context")
    expect_match(source_text, "modTypes\\(comma_example_data\\)")
    expect_match(source_text, "motifs\\(comma_example_data\\)")
    expect_match(source_text, "modContexts\\(comma_example_data\\)")
    expect_match(source_text, "siteCoverage\\(comma_example_data\\)")

    expect_false(grepl("modType(", source_text, fixed = TRUE))
    expect_false(grepl("modContext(", source_text, fixed = TRUE))
    expect_false(grepl("coverage(comma_example_data)", source_text, fixed = TRUE))
})

test_that("commaData vignette knits with package-provided data", {
    vignette_file <- .commaDataVignettePath()
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

    rendered <- paste(readLines(output_file, warn = FALSE), collapse = "\n")
    expect_match(rendered, "Understanding the commaData object")
    expect_match(rendered, "6mA_GATC|5mC_CCWGG")
    expect_match(rendered, "dm_delta_beta")
})
