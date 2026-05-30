test_that("import troubleshooting vignette covers common data import failures", {
    vignette_path <- testthat::test_path("..", "..", "vignettes", "import-troubleshooting.Rmd")
    testthat::skip_if_not(file.exists(vignette_path), "source vignettes are unavailable")

    text <- readLines(vignette_path, warn = FALSE)

    required_patterns <- c(
        "Modkit BED Format Problems",
        "Choosing The Caller",
        "Genome Size Problems",
        "Annotation File Problems",
        "expected at least 18",
        "caller = \"modkit\"",
        "caller = \"dorado\"",
        "caller = \"megalodon\"",
        "BSgenomeObject\\$NC_000913",
        "rtracklayer"
    )

    for (pattern in required_patterns) {
        expect_true(any(grepl(pattern, text)), info = pattern)
    }
})
