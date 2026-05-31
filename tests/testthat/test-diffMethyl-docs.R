test_that("diffMethyl method docs preserve default and practical guidance", {
    source_path <- c(
        test_path("..", "..", "R", "diffMethyl.R"),
        test_path("..", "..", "00_pkg_src", "commaKit", "R", "diffMethyl.R")
    )
    source <- paste(readLines(source_path[file.exists(source_path)][1]), collapse = "\n")

    vignette_path <- c(
        test_path("..", "..", "vignettes", "getting-started.Rmd"),
        test_path("..", "..", "00_pkg_src", "commaKit", "vignettes", "getting-started.Rmd")
    )
    vignette <- paste(
        readLines(vignette_path[file.exists(vignette_path)][1]),
        collapse = "\n"
    )

    expect_match(source, 'method\\s*=\\s*c\\("methylkit", "limma", "quasi_f"\\)')
    expect_match(source, 'Default method \\(\\\\code\\{method = "methylkit"\\}\\)')
    expect_match(source, 'quasi_f.*general-purpose alternative')

    expect_match(vignette, 'keeps `method = "methylkit"` as the default')
    expect_match(vignette, '`method = "quasi_f"` is a good general-purpose alternative')
    expect_match(vignette, '`method = "limma"` uses limma')
})
