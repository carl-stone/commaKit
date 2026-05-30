test_that("R version hint matches renv metadata", {
    project_root <- test_path("..", "..")
    skip_if_not(
        file.exists(file.path(project_root, ".R-version")),
        ".R-version is a development file excluded from package builds"
    )

    r_version <- readLines(file.path(project_root, ".R-version"), warn = FALSE)
    settings <- readLines(file.path(project_root, "renv", "settings.json"), warn = FALSE)

    expect_identical(r_version, "4.5.3")
    expect_true(any(grepl('"r.version": "4.5.3"', settings, fixed = TRUE)))
})

test_that("run-r-4.5 creates an overlay Rscript that uses overlay R", {
    project_root <- test_path("..", "..")
    skip_if_not(
        file.exists(file.path(project_root, "dev", "run-r-4.5.sh")),
        "development wrapper is excluded from package builds"
    )

    fake_home <- file.path(tempdir(), "fake-r-home")
    overlay_home <- file.path(tempdir(), paste0("comma-r-4.5-home.", sample.int(1e6, 1)))
    unlink(c(fake_home, overlay_home), recursive = TRUE, force = TRUE)
    dir.create(file.path(fake_home, "bin", "exec"), recursive = TRUE)
    dir.create(file.path(fake_home, "etc"), recursive = TRUE)
    dir.create(file.path(fake_home, "share"), recursive = TRUE)
    dir.create(file.path(fake_home, "include"), recursive = TRUE)
    dir.create(file.path(fake_home, "doc"), recursive = TRUE)

    fake_r <- c(
        "#!/bin/sh",
        "echo \"$0\"",
        "echo \"$@\""
    )
    writeLines(fake_r, file.path(fake_home, "bin", "exec", "R"))
    writeLines(fake_r, file.path(fake_home, "bin", "R"))
    writeLines(c("#!/bin/sh", "echo wrong-rscript"), file.path(fake_home, "bin", "Rscript"))
    writeLines(
        c(
            'LIBR = -L"/Library/Frameworks/R.framework/Resources/lib" -lR',
            'R_HOME = /Library/Frameworks/R.framework/Resources'
        ),
        file.path(fake_home, "etc", "Makeconf")
    )
    Sys.chmod(c(
        file.path(fake_home, "bin", "exec", "R"),
        file.path(fake_home, "bin", "R"),
        file.path(fake_home, "bin", "Rscript")
    ), "0755")

    script <- file.path(project_root, "dev", "run-r-4.5.sh")
    result <- system2(
        script,
        c("-e", shQuote("cat('ok')")),
        env = c(
            paste0("COMMA_R45_HOME=", fake_home),
            paste0("COMMA_R45_OVERLAY_HOME=", overlay_home)
        ),
        stdout = TRUE,
        stderr = TRUE
    )
    expect_true(any(grepl("comma-r-4.5-home.", result, fixed = TRUE)))
    expect_true(any(grepl("/bin/exec/R", result, fixed = TRUE)))

    rscript_result <- system2(
        file.path(overlay_home, "bin", "Rscript"),
        c("-e", shQuote("cat('ok')")),
        stdout = TRUE,
        stderr = TRUE
    )

    expect_true(any(grepl("comma-r-4.5-home.", rscript_result, fixed = TRUE)))
    expect_true(any(grepl("/bin/R", rscript_result, fixed = TRUE)))
    expect_true(any(grepl("--no-echo --no-restore -e", rscript_result, fixed = TRUE)))

    rscript_expr_args <- system2(
        file.path(overlay_home, "bin", "Rscript"),
        c("-e", shQuote("cat(commandArgs(TRUE))"), "vignettes/getting-started.Rmd"),
        stdout = TRUE,
        stderr = TRUE
    )
    expect_true(any(grepl("--args vignettes/getting-started.Rmd", rscript_expr_args, fixed = TRUE)))
    expect_false(any(grepl("--file=vignettes/getting-started.Rmd", rscript_expr_args, fixed = TRUE)))
})
