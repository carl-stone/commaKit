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

readiness_tracked_files <- function(root = readiness_repo_root()) {
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  files <- suppressWarnings(system2(
    "git",
    c("ls-files", "--cached"),
    stdout = TRUE,
    stderr = FALSE
  ))
  if (!is.null(attr(files, "status"))) {
    return(character())
  }
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
  if (length(max_bytes) != 1L || is.na(max_bytes) || max_bytes <= 0) {
    stop(
      "'max_bytes' must be a single positive numeric file-size limit.",
      call. = FALSE
    )
  }
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
  matches <- gregexpr("\\[[^]]+\\]\\([^)[:space:]]+\\)", text, perl = TRUE)
  captured <- regmatches(text, matches)[[1L]]
  if (length(captured) == 0L || identical(captured, character(0))) {
    return(character())
  }
  sub("^\\[[^]]+\\]\\(([^)[:space:]]+)\\)$", "\\1", captured)
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

readiness_r_quality_files <- function(files) {
  files[
    grepl("\\.[Rr]$", files) &
      grepl("^(R|inst/scripts)/", files) &
      !grepl("^(renv|docs|man|droid-wiki|\\.git)/", files)
  ]
}

readiness_descendant_ids <- function(parse_data, parent_id) {
  ids <- parent_id
  current <- parent_id
  repeat {
    children <- parse_data$id[parse_data$parent %in% current]
    children <- setdiff(children, ids)
    if (length(children) == 0L) {
      break
    }
    ids <- c(ids, children)
    current <- children
  }
  ids
}

readiness_function_name <- function(parse_data, function_row) {
  siblings <- parse_data[
    parse_data$parent == function_row$parent &
      parse_data$line1 <= function_row$line1, ,
    drop = FALSE
  ]
  symbols <- siblings[siblings$token == "SYMBOL", , drop = FALSE]
  if (nrow(symbols) > 0L) {
    return(utils::tail(symbols$text, 1L))
  }
  paste0("<anonymous>:", function_row$line1)
}

readiness_function_complexity <- function(parse_data, function_row) {
  ids <- readiness_descendant_ids(parse_data, function_row$parent)
  body <- parse_data[parse_data$id %in% ids, , drop = FALSE]
  branch_tokens <- c("IF", "FOR", "WHILE", "REPEAT", "AND", "OR")
  branch_calls <- c("ifelse", "tryCatch", "withCallingHandlers")
  1L +
    sum(body$token %in% branch_tokens) +
    sum(body$token == "SYMBOL_FUNCTION_CALL" & body$text %in% branch_calls)
}

readiness_check_complexity <- function(
  files = readiness_git_files(),
  root = readiness_repo_root(),
  max_complexity = as.integer(Sys.getenv(
    "COMMAKIT_MAX_CYCLOMATIC_COMPLEXITY",
    unset = "50"
  ))
) {
  files <- readiness_r_quality_files(files)
  violations <- character()
  max_seen <- 0L
  for (file in files) {
    path <- file.path(root, file)
    if (!file.exists(path)) next
    parsed <- tryCatch(parse(path, keep.source = TRUE), error = identity)
    if (inherits(parsed, "error")) {
      stop("Could not parse ", file, ": ", conditionMessage(parsed),
        call. = FALSE
      )
    }
    parse_data <- utils::getParseData(parsed)
    if (is.null(parse_data)) next
    functions <- parse_data[parse_data$token == "FUNCTION", , drop = FALSE]
    for (i in seq_len(nrow(functions))) {
      complexity <- readiness_function_complexity(parse_data, functions[i, ])
      max_seen <- max(max_seen, complexity)
      if (complexity > max_complexity) {
        violations <- c(
          violations,
          paste0(
            file,
            ":",
            functions$line1[[i]],
            " ",
            readiness_function_name(parse_data, functions[i, ]),
            " complexity ",
            complexity
          )
        )
      }
    }
  }
  if (length(violations) > 0L) {
    stop(
      "Functions exceed cyclomatic complexity limit ",
      max_complexity,
      ":\n- ",
      paste(violations, collapse = "\n- "),
      call. = FALSE
    )
  }
  message(
    "Complexity check passed for ",
    length(files),
    " R files; maximum observed complexity was ",
    max_seen,
    "."
  )
  invisible(TRUE)
}

readiness_internal_definitions <- function(files, root) {
  defs <- data.frame(
    name = character(),
    file = character(),
    line = integer(),
    stringsAsFactors = FALSE
  )
  definition_re <- paste0(
    "^\\s*(\\.[A-Za-z][A-Za-z0-9._]*)\\s*",
    "(<-|=)\\s*function\\s*\\(.*$"
  )
  for (file in files) {
    path <- file.path(root, file)
    if (!file.exists(path) || !readiness_is_text_file(path)) next
    lines <- readLines(path, warn = FALSE)
    hits <- grep(definition_re, lines, perl = TRUE)
    if (length(hits) == 0L) next
    defs <- rbind(
      defs,
      data.frame(
        name = sub(definition_re, "\\1", lines[hits], perl = TRUE),
        file = file,
        line = hits,
        stringsAsFactors = FALSE
      )
    )
  }
  defs
}

readiness_check_dead_code <- function(
  files = readiness_git_files(),
  root = readiness_repo_root()
) {
  r_files <- readiness_r_quality_files(files)
  searchable <- readiness_tracked_files(root)
  if (length(searchable) == 0L) {
    searchable <- files
  }
  searchable <- searchable[readiness_text_candidates(searchable)]
  defs <- readiness_internal_definitions(r_files, root)
  unused <- character()
  for (i in seq_len(nrow(defs))) {
    reference_count <- 0L
    name_re <- paste0(
      "(?<![A-Za-z0-9._])", gsub("\\.", "\\\\.", defs$name[[i]]),
      "(?![A-Za-z0-9._])"
    )
    for (file in searchable) {
      path <- file.path(root, file)
      if (!file.exists(path) || !readiness_is_text_file(path)) next
      lines <- readLines(path, warn = FALSE)
      reference_count <- reference_count +
        sum(grepl(name_re, lines, perl = TRUE))
    }
    if (reference_count <= 1L) {
      unused <- c(
        unused,
        paste0(defs$file[[i]], ":", defs$line[[i]], " ", defs$name[[i]])
      )
    }
  }
  if (length(unused) > 0L) {
    stop(
      "Internal helper definitions appear unused:\n- ",
      paste(unused, collapse = "\n- "),
      call. = FALSE
    )
  }
  message(
    "Dead internal-code check passed for ",
    nrow(defs),
    " internal helpers."
  )
  invisible(TRUE)
}

readiness_normalize_code_line <- function(line) {
  line <- trimws(line)
  line <- sub("\\s+#.*$", "", line)
  line <- gsub("\"([^\"\\\\]|\\\\.)*\"", "\"\"", line, perl = TRUE)
  line <- gsub("'([^'\\\\]|\\\\.)*'", "''", line, perl = TRUE)
  line <- gsub("[[:space:]]+", " ", line)
  trimws(line)
}

readiness_duplicate_windows <- function(lines, window_size) {
  if (length(lines) < window_size) {
    return(character())
  }
  vapply(
    seq_len(length(lines) - window_size + 1L),
    function(i) paste(lines[seq.int(i, i + window_size - 1L)], collapse = "\n"),
    character(1)
  )
}

readiness_check_duplicate_code <- function(
  files = readiness_git_files(),
  root = readiness_repo_root(),
  window_size = as.integer(Sys.getenv(
    "COMMAKIT_DUPLICATE_CODE_WINDOW",
    unset = "30"
  ))
) {
  files <- readiness_r_quality_files(files)
  files <- files[grepl("^R/", files)]
  windows <- list()
  for (file in files) {
    path <- file.path(root, file)
    if (!file.exists(path) || !readiness_is_text_file(path)) next
    lines <- readLines(path, warn = FALSE)
    lines <- vapply(lines, readiness_normalize_code_line, character(1))
    lines <- lines[nzchar(lines) & !lines %in% c("{", "}", "},", ")", "),")]
    chunks <- readiness_duplicate_windows(lines, window_size)
    if (length(chunks) == 0L) next
    windows <- c(
      windows,
      stats::setNames(
        as.list(seq_along(chunks)),
        paste(file, seq_along(chunks), chunks, sep = "\001")
      )
    )
  }
  if (length(windows) == 0L) {
    message("Duplicate-code check passed; no windows to compare.")
    return(invisible(TRUE))
  }
  keys <- sub("^[^\001]+\001[^\001]+\001", "", names(windows))
  duplicates <- names(which(table(keys) > 1L))
  if (length(duplicates) > 0L) {
    locations <- vapply(
      utils::head(duplicates, 20L),
      function(key) {
        hits <- names(windows)[keys == key]
        paste(sub("\001.*$", "", hits), collapse = ", ")
      },
      character(1)
    )
    stop(
      "Duplicate code windows of ",
      window_size,
      " normalized lines found:\n- ",
      paste(locations, collapse = "\n- "),
      call. = FALSE
    )
  }
  message(
    "Duplicate-code check passed across ",
    length(files),
    " R files with ",
    window_size,
    "-line windows."
  )
  invisible(TRUE)
}

readiness_check_code_quality <- function(
  files = readiness_git_files(),
  root = readiness_repo_root()
) {
  readiness_check_complexity(files = files, root = root)
  readiness_check_dead_code(files = files, root = root)
  readiness_check_duplicate_code(files = files, root = root)
  invisible(TRUE)
}

readiness_run <- function(checks, files = NULL, root = NULL) {
  available <- c(
    "large-files",
    "debt-markers",
    "agents-links",
    "complexity",
    "dead-code",
    "duplicate-code",
    "code-quality"
  )
  if (identical(checks, "all")) {
    checks <- c("large-files", "debt-markers", "agents-links", "code-quality")
  }
  unknown <- setdiff(checks, available)
  if (length(unknown) > 0L) {
    stop("Unknown readiness check: ", paste(unknown, collapse = ", "))
  }
  if (is.null(files)) {
    files_root <- if (is.null(root)) readiness_repo_root() else root
    files <- readiness_git_files(root = files_root)
  }
  if (is.null(root)) {
    root <- readiness_repo_root()
  }

  for (check in checks) {
    switch(check,
      "large-files" = readiness_check_large_files(files = files, root = root),
      "debt-markers" = readiness_check_debt_markers(files = files, root = root),
      "agents-links" = readiness_check_agents_links(files = files, root = root),
      "complexity" = readiness_check_complexity(files = files, root = root),
      "dead-code" = readiness_check_dead_code(files = files, root = root),
      "duplicate-code" = readiness_check_duplicate_code(
        files = files,
        root = root
      ),
      "code-quality" = readiness_check_code_quality(files = files, root = root)
    )
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L) {
    args <- "all"
  }
  checks <- args[[1L]]
  files <- args[-1L]
  if (length(files) == 0L) {
    files <- NULL
  }
  tryCatch(
    readiness_run(checks, files = files),
    error = function(error) {
      message(conditionMessage(error))
      quit(status = 1L, save = "no")
    }
  )
}
