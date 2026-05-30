test_that("package metadata uses the current repository identity", {
    desc <- packageDescription("comma")

    expect_identical(desc$Version, "0.2.0")
    expect_match(desc$URL, "github\\.com/carl-stone/comma")
    expect_match(desc$URL, "carl-stone\\.github\\.io/comma")
    expect_match(desc$BugReports, "github\\.com/carl-stone/comma/issues")
    expect_no_match(paste(desc$URL, desc$BugReports), "CoMMA|commaKit")
})
