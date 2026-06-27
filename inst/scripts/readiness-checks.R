#!/usr/bin/env Rscript

readiness_repo_root <- function() {
  root <- suppressWarnings(system2(
    "git",
    c("rev-parse", "--show-toplevel"),
    stdout = TRUE,
    stderr = FALSE
  ))
  if (length(root) == 1L && nzchar(root)) {
    return(normalizePath(root, winslash = "/", mustWork = TRUE))
  }

  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  if (nzchar(workspace) && dir.exists(file.path(workspace, ".git"))) {
    return(normalizePath(workspace, winslash = "/", mustWork = TRUE))
  }

  if (dir.exists(".git")) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }

  stop("Could not determine git repository root.", call. = FALSE)
}

readiness_git_files <- function(root = readiness_repo_root()) {
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  files <- system2(
    "git",
    c("ls-files", "--cached", "--others", "--exclude-standard"),
    stdout = TRUE
  )
  files[nzchar(files)]
}

readiness_abs_path <- function(root, files) {
  normalizePath(file.path(root, files), winslash = "/", mustWork = FALSE)
}

readiness_is_text_file <- function(path) {
  info <- file.info(path)
  if (is.na(info$size) || info$size == 0L || info$size > 1024^2) {
    return(FALSE)
  }
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  sample <- readBin(con, "raw", n = min(info$size, 4096L))
  !any(sample == as.raw(0L))
}

readiness_format_bytes <- function(bytes) {
  units <- c("B", "KB", "MB", "GB")
  value <- as.numeric(bytes)
  unit <- 1L
  while (value >= 1024 && unit < length(units)) {
    value <- value / 1024
    unit <- unit + 1L
  }
  sprintf("%.1f %s", value, units[[unit]])
}

readiness_check_large_files <- function(
  files = readiness_git_files(),
  root = readiness_repo_root(),
  max_bytes = as.numeric(Sys.getenv(
    "COMMAKIT_MAX_FILE_BYTES",
    unset = 5 * 1024^2
  ))
) {
  paths <- readiness_abs_path(root, files)
  sizes <- file.info(paths)$size
  too_large <- !is.na(sizes) & sizes > max_bytes
  if (any(too_large)) {
    details <- paste0(
      files[too_large],
      " (",
      vapply(sizes[too_large], readiness_format_bytes, character(1)),
      ")"
    )
    stop(
      "Files exceed ",
      readiness_format_bytes(max_bytes),
      " limit:\n- ",
      paste(details, collapse = "\n- "),
      call. = FALSE
    )
  }
  message("Large-file check passed for ", length(files), " files.")
  invisible(TRUE)
}

readiness_text_candidates <- function(files) {
  grepl(
    "\\.(R|Rmd|md|ya?ml|json|toml|txt|sh|Dockerfile)$|(^|/)Dockerfile$",
    files,
    ignore.case = TRUE
  )
}

readiness_check_debt_markers <- function(
  files = readiness_git_files(),
  root = readiness_repo_root()
) {
  files <- files[readiness_text_candidates(files)]
  marker_words <- c(paste0("TO", "DO"), paste0("FIX", "ME"))
  marker_re <- paste0("\\b(", paste(marker_words, collapse = "|"), ")\\b")
  issue_re <- paste(
    "(#[0-9]+",
    "[A-Z][A-Z0-9]+-[0-9]+",
    "https://github\\.com/[^[:space:]]+/issues/[0-9]+)",
    sep = "|"
  )
  violations <- character()
  for (file in files) {
    path <- file.path(root, file)
    if (!file.exists(path) || !readiness_is_text_file(path)) next
    lines <- readLines(path, warn = FALSE)
    hits <- grepl(marker_re, lines, ignore.case = TRUE) &
      !grepl(issue_re, lines)
    if (any(hits)) {
      violations <- c(
        violations,
        paste0(file, ":", which(hits), ": ", trimws(lines[hits]))
      )
    }
  }
  if (length(violations) > 0L) {
    stop(
      "Debt markers must reference an issue or ticket:\n- ",
      paste(violations, collapse = "\n- "),
      call. = FALSE
    )
  }
  message("Debt-marker check passed for ", length(files), " text files.")
  invisible(TRUE)
}

readiness_markdown_links <- function(path) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  matches <- gregexec("\\[[^]]+\\]\\(([^)[:space:]]+)", text, perl = TRUE)
  captured <- regmatches(text, matches)[[1L]]
  if (length(captured) == 0L) {
    return(character())
  }
  if (is.matrix(captured)) {
    return(captured[2L, ])
  }
  sub("^\\[[^]]+\\]\\(([^)[:space:]]+).*$", "\\1", captured)
}

readiness_is_local_link <- function(link) {
  !grepl("^(#|https?://|mailto:)", link, ignore.case = TRUE)
}

readiness_check_agents_links <- function(
  files = readiness_git_files(),
  root = readiness_repo_root()
) {
  agents_files <- files[basename(files) == "AGENTS.md"]
  missing <- character()
  for (file in agents_files) {
    path <- file.path(root, file)
    links <- readiness_markdown_links(path)
    links <- links[readiness_is_local_link(links)]
    links <- sub("#.*$", "", links)
    links <- links[nzchar(links)]
    for (link in links) {
      target <- normalizePath(
        file.path(dirname(path), link),
        winslash = "/",
        mustWork = FALSE
      )
      if (!file.exists(target)) {
        missing <- c(missing, paste0(file, " -> ", link))
      }
    }
  }
  if (length(missing) > 0L) {
    stop(
      "AGENTS.md contains missing local links:\n- ",
      paste(missing, collapse = "\n- "),
      call. = FALSE
    )
  }
  message("AGENTS.md link check passed for ", length(agents_files), " files.")
  invisible(TRUE)
}

readiness_run <- function(checks) {
  available <- c("large-files", "debt-markers", "agents-links")
  if (identical(checks, "all")) {
    checks <- available
  }
  unknown <- setdiff(checks, available)
  if (length(unknown) > 0L) {
    stop("Unknown readiness check: ", paste(unknown, collapse = ", "))
  }

  for (check in checks) {
    switch(check,
      "large-files" = readiness_check_large_files(),
      "debt-markers" = readiness_check_debt_markers(),
      "agents-links" = readiness_check_agents_links()
    )
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L) {
    args <- "all"
  }
  tryCatch(
    readiness_run(args),
    error = function(error) {
      message(conditionMessage(error))
      quit(status = 1L, save = "no")
    }
  )
}
