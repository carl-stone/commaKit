.comma_log_enabled <- function() {
  opt <- getOption("commaKit.log", FALSE)
  env <- Sys.getenv("COMMAKIT_LOG", unset = "")
  enabled_values <- c("1", "true", "yes", "json", "structured")

  isTRUE(opt) ||
    tolower(as.character(opt)[1L]) %in% enabled_values ||
    tolower(env) %in% enabled_values
}

.comma_log_connection <- function() {
  getOption("commaKit.log.connection", stdout())
}

.comma_log_event <- function(event,
                             level = "info",
                             component = NULL,
                             ...,
                             .time = Sys.time(),
                             .stream = .comma_log_connection()) {
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
  fields <- .comma_redact_fields(fields)

  .comma_add_breadcrumb(
    message = event,
    category = if (is.null(component)) "commaKit" else component,
    level = level,
    data = fields,
    .time = .time
  )

  if (!.comma_log_enabled()) {
    return(invisible(FALSE))
  }

  ok <- tryCatch(
    {
      cat(.comma_json_object(fields), "\n", sep = "", file = .stream)
      TRUE
    },
    error = function(e) FALSE
  )
  invisible(ok)
}

.comma_error_tracking_enabled <- function() {
  opt <- getOption("commaKit.error_tracking", FALSE)
  env <- Sys.getenv("COMMAKIT_ERROR_TRACKING", unset = "")
  enabled_values <- c("1", "true", "yes", "sentry", "structured")

  isTRUE(opt) ||
    tolower(as.character(opt)[1L]) %in% enabled_values ||
    tolower(env) %in% enabled_values ||
    nzchar(.comma_sentry_dsn())
}

.comma_sentry_dsn <- function() {
  opt <- getOption("commaKit.sentry.dsn", "")
  env <- Sys.getenv("COMMAKIT_SENTRY_DSN", unset = "")
  sentry_env <- Sys.getenv("SENTRY_DSN", unset = "")
  dsn <- as.character(c(opt, env, sentry_env))
  dsn <- dsn[!is.na(dsn) & nzchar(dsn)]
  if (length(dsn) == 0L) {
    return("")
  }
  as.character(dsn[1L])
}

.comma_add_breadcrumb <- function(message,
                                  category = "commaKit",
                                  level = "info",
                                  data = list(),
                                  .time = Sys.time()) {
  if (!.comma_error_tracking_enabled()) {
    return(invisible(FALSE))
  }

  breadcrumb <- list(
    timestamp = format(
      as.POSIXct(.time, tz = "UTC"),
      "%Y-%m-%dT%H:%M:%OS3Z",
      tz = "UTC"
    ),
    type = "default",
    category = category,
    level = level,
    message = message,
    data = .comma_sanitize_context(data)
  )
  max_breadcrumbs <- getOption("commaKit.error_tracking.max_breadcrumbs", 50L)
  max_breadcrumbs <- suppressWarnings(as.integer(max_breadcrumbs)[1L])
  if (is.na(max_breadcrumbs) || max_breadcrumbs < 1L) {
    max_breadcrumbs <- 50L
  }
  breadcrumbs <- getOption("commaKit.error_tracking.breadcrumbs", list())
  breadcrumbs <- c(breadcrumbs, list(breadcrumb))
  if (length(breadcrumbs) > max_breadcrumbs) {
    breadcrumbs <- breadcrumbs[
      seq.int(length(breadcrumbs) - max_breadcrumbs + 1L, length(breadcrumbs))
    ]
  }
  options(commaKit.error_tracking.breadcrumbs = breadcrumbs)
  invisible(TRUE)
}

.comma_track_error <- function(condition,
                               component,
                               operation,
                               ...,
                               .time = Sys.time(),
                               .calls = sys.calls()) {
  context <- .comma_sanitize_context(list(...))
  .comma_log_event(
    "error_captured",
    level = "error",
    component = component,
    operation = operation,
    error_class = class(condition),
    error_message = conditionMessage(condition),
    context = context,
    .time = .time
  )

  if (!.comma_error_tracking_enabled()) {
    return(invisible(FALSE))
  }

  event <- .comma_error_event(
    condition = condition,
    component = component,
    operation = operation,
    context = context,
    .time = .time,
    .calls = .calls
  )
  ok <- .comma_send_error_event(event)
  invisible(ok)
}

.comma_error_event <- function(condition,
                               component,
                               operation,
                               context = list(),
                               .time = Sys.time(),
                               .calls = sys.calls()) {
  list(
    event_id = .comma_event_id(),
    timestamp = format(
      as.POSIXct(.time, tz = "UTC"),
      "%Y-%m-%dT%H:%M:%OS3Z",
      tz = "UTC"
    ),
    platform = "r",
    logger = "commaKit",
    level = "error",
    release = .comma_release(),
    environment = .comma_environment(),
    transaction = operation,
    tags = list(
      package = "commaKit",
      component = component,
      operation = operation
    ),
    user = .comma_user_context(),
    contexts = list(
      runtime = list(name = "R", version = as.character(getRversion())),
      trace = list(operation = operation)
    ),
    extra = context,
    breadcrumbs = list(
      values = getOption("commaKit.error_tracking.breadcrumbs", list())
    ),
    exception = list(
      values = list(list(
        type = class(condition)[1L],
        value = conditionMessage(condition),
        stacktrace = list(frames = .comma_stack_frames(.calls))
      ))
    )
  )
}

.comma_send_error_event <- function(event) {
  stored <- .comma_store_error_event(event)
  reporter <- getOption("commaKit.error_tracking.reporter", NULL)
  if (is.function(reporter)) {
    return(tryCatch(isTRUE(reporter(event)), error = function(e) FALSE))
  }

  dsn <- .comma_sentry_dsn()
  if (!nzchar(dsn) || !requireNamespace("curl", quietly = TRUE)) {
    return(stored)
  }

  endpoint <- .comma_sentry_envelope_url(dsn)
  if (!nzchar(endpoint)) {
    return(stored)
  }

  envelope <- paste(
    .comma_json_object(list(
      event_id = event$event_id,
      dsn = dsn,
      sent_at = event$timestamp
    )),
    .comma_json_object(list(type = "event")),
    .comma_json_object(event),
    sep = "\n"
  )
  tryCatch(
    {
      response <- curl::curl_fetch_memory(
        endpoint,
        handle = curl::new_handle(
          post = TRUE,
          postfields = charToRaw(envelope),
          timeout = .comma_error_tracking_timeout(),
          connecttimeout = .comma_connect_timeout(),
          httpheader = c("Content-Type" = "application/x-sentry-envelope")
        )
      )
      response$status_code >= 200L && response$status_code < 300L
    },
    error = function(e) FALSE
  )
}

.comma_store_error_event <- function(event) {
  max_events <- getOption("commaKit.error_tracking.max_events", 25L)
  max_events <- suppressWarnings(as.integer(max_events)[1L])
  if (is.na(max_events) || max_events < 1L) {
    max_events <- 25L
  }
  events <- getOption("commaKit.error_tracking.events", list())
  if (!is.list(events)) {
    events <- list()
  }
  events <- c(events, list(event))
  if (length(events) > max_events) {
    events <- events[seq.int(length(events) - max_events + 1L, length(events))]
  }
  options(commaKit.error_tracking.events = events)
  invisible(TRUE)
}

.comma_error_tracking_timeout <- function() {
  timeout <- getOption(
    "commaKit.error_tracking.timeout",
    Sys.getenv("COMMAKIT_ERROR_TRACKING_TIMEOUT", unset = "2")
  )
  timeout <- suppressWarnings(as.numeric(timeout)[1L])
  if (is.na(timeout) || timeout <= 0) {
    return(2)
  }
  timeout
}

.comma_connect_timeout <- function() {
  timeout <- getOption(
    "commaKit.error_tracking.connect_timeout",
    Sys.getenv("COMMAKIT_ERROR_TRACKING_CONNECT_TIMEOUT", unset = "1")
  )
  timeout <- suppressWarnings(as.numeric(timeout)[1L])
  if (is.na(timeout) || timeout <= 0) {
    return(1)
  }
  timeout
}

.comma_sentry_envelope_url <- function(dsn) {
  parsed <- utils::URLdecode(dsn)
  match <- regexec(
    "^([^:]+)://([^@]+)@([^/]+)(/.*)?/([^/]+)$",
    parsed
  )
  parts <- regmatches(parsed, match)[[1L]]
  if (length(parts) != 6L) {
    return("")
  }
  protocol <- parts[2L]
  host <- parts[4L]
  path <- parts[5L]
  project_id <- parts[6L]
  if (is.na(path) || !nzchar(path)) {
    path <- ""
  }
  paste0(protocol, "://", host, path, "/api/", project_id, "/envelope/")
}

.comma_stack_frames <- function(calls) {
  if (length(calls) == 0L) {
    return(list())
  }
  calls <- tail(calls, 40L)
  frames <- lapply(calls, function(call) {
    frame <- list(
      "function" = .comma_call_function(call),
      context_line = paste(deparse(call), collapse = " ")
    )
    srcref <- attr(call, "srcref", exact = TRUE)
    if (!is.null(srcref)) {
      srcfile <- attr(srcref, "srcfile", exact = TRUE)
      frame$filename <- if (!is.null(srcfile)) srcfile$filename else NULL
      frame$abs_path <- frame$filename
      frame$lineno <- as.integer(srcref[1L])
      frame$colno <- as.integer(srcref[5L])
    }
    frame[!vapply(frame, is.null, logical(1))]
  })
  frames
}

.comma_call_function <- function(call) {
  if (!is.call(call) || length(call) == 0L) {
    return("<unknown>")
  }
  fun <- call[[1L]]
  paste(deparse(fun), collapse = "::")
}

.comma_user_context <- function() {
  user <- getOption("commaKit.error_tracking.user", NULL)
  if (is.null(user)) {
    user <- list(
      id = Sys.getenv("COMMAKIT_USER_ID", unset = ""),
      email = Sys.getenv("COMMAKIT_USER_EMAIL", unset = ""),
      username = Sys.getenv("COMMAKIT_USER_NAME", unset = "")
    )
  }
  if (!is.list(user)) {
    user <- list(id = as.character(user)[1L])
  }
  user <- .comma_sanitize_context(user)
  present <- vapply(
    user,
    function(value) {
      length(value) > 0L && nzchar(as.character(value)[1L])
    },
    logical(1)
  )
  user[present]
}

.comma_environment <- function() {
  env <- getOption(
    "commaKit.error_tracking.environment",
    Sys.getenv("COMMAKIT_ENVIRONMENT", unset = "")
  )
  env <- as.character(env)[1L]
  if (is.na(env) || !nzchar(env)) {
    env <- Sys.getenv("SENTRY_ENVIRONMENT", unset = "production")
  }
  env <- as.character(env)[1L]
  if (is.na(env) || !nzchar(env)) {
    return("production")
  }
  env
}

.comma_release <- function() {
  release <- getOption(
    "commaKit.error_tracking.release",
    Sys.getenv("SENTRY_RELEASE", unset = "")
  )
  release <- as.character(release)[1L]
  if (!is.na(release) && nzchar(release)) {
    return(release)
  }
  version <- tryCatch(
    as.character(utils::packageVersion("commaKit")),
    error = function(e) "unknown"
  )
  paste0("commaKit@", version)
}

.comma_event_id <- function() {
  counter <- getOption("commaKit.error_tracking.event_counter", 0L)
  counter <- suppressWarnings(as.integer(counter)[1L])
  if (is.na(counter)) {
    counter <- 0L
  }
  counter <- counter + 1L
  options(commaKit.error_tracking.event_counter = counter)

  seed <- paste(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    Sys.getpid(),
    counter,
    sep = ":"
  )
  ints <- utf8ToInt(seed)
  hashes <- vapply(
    0:3,
    function(offset) {
      sum((ints + offset) * seq_along(ints)) %% 4294967296
    },
    numeric(1)
  )
  paste(vapply(hashes, .comma_hex8, character(1)), collapse = "")
}

.comma_hex8 <- function(value) {
  digits <- strsplit("0123456789abcdef", "", fixed = TRUE)[[1L]]
  value <- floor(as.numeric(value))
  chars <- character(8L)
  for (i in 8:1) {
    chars[[i]] <- digits[[value %% 16L + 1L]]
    value <- floor(value / 16L)
  }
  paste(chars, collapse = "")
}

.comma_sanitize_context <- function(context) {
  if (!is.list(context)) {
    return(as.character(context)[1L])
  }
  field_names <- names(context)
  if (is.null(field_names)) {
    field_names <- rep("", length(context))
  }
  cleaned <- Map(function(value, name) {
    if (.comma_sensitive_field(name)) {
      return("[REDACTED]")
    }
    if (is.null(value) || length(value) == 0L || all(is.na(value))) {
      return(NULL)
    }
    if (is.list(value)) {
      return(.comma_sanitize_context(value))
    }
    if (inherits(value, c("POSIXct", "POSIXlt", "Date"))) {
      return(value[1L])
    }
    if (length(value) > 20L) {
      value <- c(value[seq_len(20L)], "...")
    }
    as.character(value)
  }, context, field_names)
  cleaned[!vapply(cleaned, is.null, logical(1))]
}

.comma_redact_fields <- function(fields) {
  if (!is.list(fields)) {
    return(fields)
  }
  field_names <- names(fields)
  if (is.null(field_names)) {
    field_names <- rep("", length(fields))
  }
  redacted <- Map(function(value, name) {
    if (.comma_sensitive_field(name)) {
      return("[REDACTED]")
    }
    if (is.list(value)) {
      return(.comma_redact_fields(value))
    }
    value
  }, fields, field_names)
  if (!is.null(names(fields))) {
    names(redacted) <- names(fields)
  }
  redacted
}

.comma_sensitive_field <- function(name) {
  if (is.null(name) || is.na(name) || !nzchar(name)) {
    return(FALSE)
  }
  pattern <- paste(
    "password",
    "passwd",
    "secret",
    "token",
    "api[_-]?key",
    "access[_-]?key",
    "private[_-]?key",
    "dsn",
    "authorization",
    "cookie",
    "email",
    sep = "|"
  )
  grepl(
    pattern,
    name,
    ignore.case = TRUE
  )
}

.comma_json_object <- function(fields) {
  invalid_names <- is.null(names(fields)) ||
    any(is.na(names(fields))) ||
    any(names(fields) == "")
  if (invalid_names) {
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
    return(format(value,
      scientific = FALSE, trim = TRUE,
      decimal.mark = "."
    ))
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
            code <- utf8ToInt(character)
            if (code < 32L) {
              return(switch(character,
                "\b" = "\\b",
                "\f" = "\\f",
                "\n" = "\\n",
                "\r" = "\\r",
                "\t" = "\\t",
                sprintf("\\u%04x", code)
              ))
            }
            switch(character,
              "\\" = "\\\\",
              "\"" = "\\\"",
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
