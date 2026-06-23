.comma_log_enabled <- function() {
  opt <- getOption("commaKit.log", FALSE)
  env <- Sys.getenv("COMMAKIT_LOG", unset = "")
  enabled_values <- c("1", "true", "yes", "json", "structured")

  isTRUE(opt) ||
    tolower(as.character(opt)[1L]) %in% enabled_values ||
    tolower(env) %in% enabled_values
}

.comma_log_connection <- function() {
  getOption("commaKit.log.connection", stderr())
}

.comma_log_event <- function(event,
                             level = "info",
                             component = NULL,
                             ...,
                             .time = Sys.time(),
                             .stream = .comma_log_connection()) {
  if (!.comma_log_enabled()) {
    return(invisible(FALSE))
  }

  fields <- list(
    timestamp = format(
      as.POSIXct(.time, tz = "UTC"),
      "%Y-%m-%dT%H:%M:%OS3Z",
      tz = "UTC"
    ),
    level = level,
    event = event,
    package = "commaKit"
  )
  if (!is.null(component)) {
    fields$component <- component
  }
  fields <- c(fields, list(...))
  fields <- fields[!vapply(fields, is.null, logical(1))]

  cat(.comma_json_object(fields), "\n", sep = "", file = .stream)
  invisible(TRUE)
}

.comma_json_object <- function(fields) {
  if (is.null(names(fields)) ||
    any(is.na(names(fields))) ||
    any(names(fields) == "")) {
    stop("Structured log fields must be named.", call. = FALSE)
  }

  entries <- paste(
    vapply(
      names(fields),
      function(name) {
        paste0(.comma_json_string(name), ":", .comma_json_value(fields[[name]]))
      },
      character(1)
    ),
    collapse = ","
  )
  paste0("{", entries, "}")
}

.comma_json_value <- function(value) {
  if (is.null(value) || length(value) == 0L || all(is.na(value))) {
    return("null")
  }

  if (inherits(value, c("POSIXct", "POSIXlt"))) {
    return(.comma_json_string(format(
      as.POSIXct(value[1L], tz = "UTC"),
      "%Y-%m-%dT%H:%M:%OS3Z",
      tz = "UTC"
    )))
  }

  if (inherits(value, "Date")) {
    return(.comma_json_string(format(value[1L], "%Y-%m-%d")))
  }

  if (is.factor(value)) {
    value <- as.character(value)
  }

  if (is.list(value)) {
    if (!is.null(names(value)) && all(nzchar(names(value)))) {
      return(.comma_json_object(value))
    }
    return(paste0(
      "[",
      paste(vapply(value, .comma_json_value, character(1)), collapse = ","),
      "]"
    ))
  }

  if (length(value) > 1L) {
    return(paste0(
      "[",
      paste(vapply(as.list(value), .comma_json_value, character(1)),
        collapse = ","
      ),
      "]"
    ))
  }

  if (is.logical(value)) {
    return(if (isTRUE(value)) "true" else "false")
  }

  if (is.numeric(value)) {
    if (!is.finite(value)) {
      return("null")
    }
    return(format(value, scientific = FALSE, trim = TRUE))
  }

  .comma_json_string(as.character(value))
}

.comma_json_string <- function(value) {
  chars <- strsplit(enc2utf8(value), "", fixed = TRUE)
  escaped <- vapply(
    chars,
    function(characters) {
      paste(
        vapply(
          characters,
          function(character) {
            switch(character,
              "\\" = "\\\\",
              "\"" = "\\\"",
              "\n" = "\\n",
              "\r" = "\\r",
              "\t" = "\\t",
              character
            )
          },
          character(1)
        ),
        collapse = ""
      )
    },
    character(1)
  )
  paste0("\"", escaped, "\"")
}
