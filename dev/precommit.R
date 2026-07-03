#!/usr/bin/env Rscript

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    message(
      "Package '", package, "' is required for this pre-commit hook. ",
      "Install project development dependencies before retrying."
    )
    quit(status = 1, save = "no")
  }
}

is_r_file <- function(path) {
  grepl("\\.[Rr](md)?$", path) || basename(path) == ".Rprofile"
}

is_lintr_file <- function(path) {
  basename(path) == ".lintr"
}

existing_r_files <- function(files) {
  files[file.exists(files) & vapply(files, is_r_file, logical(1))]
}

changed_lintr_config <- function(files) {
  any(file.exists(files) & vapply(files, is_lintr_file, logical(1)))
}

lint_files <- function(files) {
  lints <- lapply(files, lintr::lint)
  structure(do.call(c, lints), class = c("lints", "list"))
}

check_style <- function(files) {
  require_package("styler")

  files <- existing_r_files(files)
  result <- tryCatch(
    {
      if (length(files) == 0) {
        styler::style_pkg(dry = "fail")
      } else {
        styler::style_file(files, dry = "fail")
      }
      TRUE
    },
    error = function(e) e
  )

  if (inherits(result, "error")) {
    message("R style check failed.")
    message("Run this command, stage the formatting changes, and retry:")
    message("  Rscript -e 'styler::style_pkg()'")
    message("")
    message("Original styler error:")
    message(conditionMessage(result))
    quit(status = 1, save = "no")
  }

  invisible(TRUE)
}


existing_rmd_files <- function(files) {
  rmd_files <- files[file.exists(files) & grepl("\\.[Rr]md$", files)]
  md_files <- files[
    file.exists(files) & grepl("\\.md$", files, ignore.case = TRUE)
  ]
  rmd_siblings <- sub("\\.md$", ".Rmd", md_files, ignore.case = TRUE)
  unique(c(rmd_files, rmd_siblings[file.exists(rmd_siblings)]))
}

rendered_md_path <- function(rmd_file) {
  sub("\\.[Rr]md$", ".md", rmd_file)
}

check_rmarkdown_rendered <- function(files) {
  require_package("rmarkdown")

  requested_files <- files
  files <- existing_rmd_files(files)
  if (length(files) == 0L && length(requested_files) == 0L) {
    files <- existing_rmd_files(Sys.glob(c("*.Rmd", "vignettes/*.Rmd")))
  }

  output_files <- vapply(files, rendered_md_path, character(1))
  missing <- files[!file.exists(output_files)]
  files <- files[file.exists(output_files)]
  if (length(files) == 0L && length(missing) == 0L) {
    message("No R Markdown files with Markdown outputs to check.")
    return(invisible(TRUE))
  }

  stale <- if (length(missing) > 0L) {
    paste0(missing, " -> ", rendered_md_path(missing), " (missing)")
  } else {
    character()
  }
  for (file in files) {
    output_file <- rendered_md_path(file)
    output_dir <- tempfile("rmd-render-")
    dir.create(output_dir)
    on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

    rendered <- rmarkdown::render(
      file,
      output_format = "github_document",
      output_dir = output_dir,
      quiet = TRUE,
      envir = new.env(parent = globalenv())
    )
    expected <- readLines(output_file, warn = FALSE)
    actual <- readLines(rendered, warn = FALSE)
    if (!identical(expected, actual)) {
      stale <- c(stale, paste0(file, " -> ", output_file))
    }
  }

  if (length(stale) > 0L) {
    message("R Markdown outputs are stale.")
    message(
      "Render the stale source file(s), stage the generated output, ",
      "and retry."
    )
    message("For README changes, run:")
    message("  Rscript -e 'rmarkdown::render(\"README.Rmd\")'")
    message("")
    message("Stale pairs:")
    message("- ", paste(stale, collapse = "\n- "))
    quit(status = 1, save = "no")
  }

  message("R Markdown outputs are up to date.")
  invisible(TRUE)
}

check_lint <- function(files) {
  require_package("lintr")

  if (changed_lintr_config(files)) {
    lints <- lintr::lint_package()
  } else {
    files <- existing_r_files(files)
    if (length(files) == 0) {
      message("No changed R files to lint.")
      return(invisible())
    }

    lints <- lint_files(files)
  }

  if (length(lints) == 0) {
    message("No lint issues found.")
  }

  print(lints)

  if (length(lints) > 0) {
    quit(status = 1, save = "no")
  }
}

check_tests <- function() {
  require_package("devtools")
  devtools::test()
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  message(
    "Usage: Rscript dev/precommit.R ",
    "<style|lint|rmarkdown|test> [files...]"
  )
  quit(status = 1, save = "no")
}

command <- args[[1]]
files <- args[-1]

switch(command,
  style = check_style(files),
  lint = check_lint(files),
  rmarkdown = check_rmarkdown_rendered(files),
  test = check_tests(),
  {
    message("Unknown pre-commit command: ", command)
    quit(status = 1, save = "no")
  }
)
