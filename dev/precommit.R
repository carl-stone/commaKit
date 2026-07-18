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

  requested_files <- files
  files <- existing_r_files(files)
  result <- tryCatch(
    {
      if (length(files) == 0 && length(requested_files) > 0L) {
        message("No changed R files to style.")
      } else if (length(files) == 0) {
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

is_roxygen_relevant_file <- function(path) {
  grepl("^R/.*\\.[Rr]$", path) ||
    grepl("^man/.*\\.Rd$", path) ||
    basename(path) %in% c("DESCRIPTION", "NAMESPACE")
}

package_files <- function() {
  files <- system2("git", c("ls-files"), stdout = TRUE)
  files[nzchar(files) & file.exists(files)]
}

copy_package_for_document <- function() {
  source_files <- package_files()
  document_outputs <- grepl("^man/.*\\.Rd$", source_files)
  source_files <- source_files[!document_outputs]

  pkg_dir <- file.path(tempdir(), paste0("commakit-document-", Sys.getpid()))
  if (dir.exists(pkg_dir)) {
    unlink(pkg_dir, recursive = TRUE)
  }
  dir.create(pkg_dir, recursive = TRUE)

  for (file in source_files) {
    dest <- file.path(pkg_dir, file)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file, dest, overwrite = TRUE)) {
      stop("Could not copy package file to temporary document check: ", file)
    }
  }

  pkg_dir
}

read_file_raw <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readBin(path, what = "raw", n = file.info(path)$size)
}

generated_document_paths <- function(root) {
  man_dir <- file.path(root, "man")
  rd_files <- if (dir.exists(man_dir)) {
    file.path("man", list.files(man_dir, pattern = "\\.Rd$"))
  } else {
    character()
  }
  c("NAMESPACE", rd_files)
}

check_roxygen_documented <- function(files) {
  require_package("devtools")

  requested_files <- files
  has_requested_files <- length(requested_files) > 0L
  has_roxygen_files <- any(vapply(
    requested_files,
    is_roxygen_relevant_file,
    logical(1)
  ))
  if (has_requested_files && !has_roxygen_files) {
    message("No roxygen source or generated documentation files to check.")
    return(invisible(TRUE))
  }

  pkg_dir <- copy_package_for_document()
  on.exit(unlink(pkg_dir, recursive = TRUE), add = TRUE)

  result <- tryCatch(
    {
      devtools::document(pkg = pkg_dir, quiet = TRUE)
      TRUE
    },
    error = function(e) e
  )

  if (inherits(result, "error")) {
    message("Could not run devtools::document() in a temporary package copy.")
    message("Original error:")
    message(conditionMessage(result))
    quit(status = 1, save = "no")
  }

  expected <- generated_document_paths(pkg_dir)
  actual <- generated_document_paths(".")
  paths <- sort(unique(c(expected, actual)))
  stale <- character()

  for (path in paths) {
    expected_raw <- read_file_raw(file.path(pkg_dir, path))
    actual_raw <- read_file_raw(path)
    if (!identical(expected_raw, actual_raw)) {
      stale <- c(stale, path)
    }
  }

  if (length(stale) > 0L) {
    message("Roxygen-generated documentation is stale.")
    message("Run this command, stage the generated changes, and retry:")
    message("  Rscript -e 'devtools::document()'")
    message("")
    message("Stale generated file(s):")
    message("- ", paste(stale, collapse = "\n- "))
    quit(status = 1, save = "no")
  }

  message("Roxygen-generated documentation is up to date.")
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

trim_trailing_whitespace <- function(lines) {
  sub("[[:space:]]+$", "", lines)
}

normalize_pipe_table <- function(line) {
  if (!grepl("^\\|.*\\|$", line)) {
    return(line)
  }
  contents <- sub("^\\||\\|$", "", line)
  cells <- trimws(strsplit(contents, "\\|", fixed = FALSE)[[1]])
  cells[grepl("^:?-+:?$", cells)] <- "---"
  paste0("| ", paste(cells, collapse = " | "), " |")
}

normalize_rendered_markdown <- function(lines) {
  vapply(trim_trailing_whitespace(lines), normalize_pipe_table, character(1))
}

report_render_diff <- function(expected_path, actual_path) {
  diff <- suppressWarnings(system2(
    "diff",
    c(
      "-u", "--label", expected_path, expected_path,
      "--label", actual_path, actual_path
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  if (length(diff) > 0L) {
    message("Rendered-output diff:")
    message(paste(diff, collapse = "\n"))
  }
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
    render_dir <- tempfile("rmd-render-")
    dir.create(render_dir)
    on.exit(unlink(render_dir, recursive = TRUE), add = TRUE)
    render_input <- file.path(render_dir, basename(file))
    if (!file.copy(file, render_input, overwrite = TRUE)) {
      stop("Could not copy R Markdown source to temporary render directory.")
    }

    rendered <- rmarkdown::render(
      render_input,
      output_format = "github_document",
      knit_root_dir = normalizePath(dirname(file), mustWork = TRUE),
      quiet = TRUE,
      envir = new.env(parent = globalenv())
    )
    expected <- normalize_rendered_markdown(
      readLines(output_file, warn = FALSE)
    )
    actual <- normalize_rendered_markdown(readLines(rendered, warn = FALSE))
    if (!identical(expected, actual)) {
      stale <- c(stale, paste0(file, " -> ", output_file))
      report_render_diff(output_file, rendered)
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
    "<style|lint|roxygen|rmarkdown|test> [files...]"
  )
  quit(status = 1, save = "no")
}

command <- args[[1]]
files <- args[-1]

switch(command,
  style = check_style(files),
  lint = check_lint(files),
  roxygen = check_roxygen_documented(files),
  rmarkdown = check_rmarkdown_rendered(files),
  test = check_tests(),
  {
    message("Unknown pre-commit command: ", command)
    quit(status = 1, save = "no")
  }
)
