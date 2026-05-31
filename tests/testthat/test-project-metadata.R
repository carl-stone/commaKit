test_that("package metadata uses the commaKit identity", {
    desc <- packageDescription("commaKit")

    expect_identical(desc$Package, "commaKit")
    expect_identical(desc$Version, "0.2.0")
    expect_match(desc$Title, "Comparative Microbial Methylomics Analysis Kit")
    expect_match(desc$URL, "github\\.com/carl-stone/commaKit")
    expect_match(desc$URL, "carl-stone\\.github\\.io/commaKit")
    expect_match(desc$BugReports, "github\\.com/carl-stone/commaKit/issues")
    expect_no_match(paste(desc$URL, desc$BugReports), "github\\.com/carl-stone/comma($|[^K])")
})
