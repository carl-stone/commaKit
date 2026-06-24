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
  if (length(files) == 0) {
    styler::style_pkg(dry = "fail")
  } else {
    styler::style_file(files, dry = "fail")
  }
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
  message("Usage: Rscript dev/precommit.R <style|lint|test> [files...]")
  quit(status = 1, save = "no")
}

command <- args[[1]]
files <- args[-1]

switch(command,
  style = check_style(files),
  lint = check_lint(files),
  test = check_tests(),
  {
    message("Unknown pre-commit command: ", command)
    quit(status = 1, save = "no")
  }
)
