#' @importFrom utils read.table
#' @importFrom stats aggregate
NULL

#' Parse a Megalodon per-read methylation file into a tidy per-site data frame
#'
#' Reads a Megalodon per-read modification output file, aggregates per-read
#' calls to per-site beta values and coverage, and returns a tidy data frame
#' compatible with the \code{\link{commaData}} constructor. This is an internal
#' function called when \code{caller = "megalodon"}.
#'
#' @details
#' Megalodon per-read output (\code{modified_bases.5mC.bed} or similar) has
#' the format:
#' \preformatted{
#'   chrom  start  end  read_id  score  strand  ...  mod_prob
#' }
#' where \code{mod_prob} is the per-read probability of modification. This
#' function aggregates across reads at each site by computing:
#' \itemize{
#'   \item \code{beta} = mean of per-read probabilities at each site
#'   \item \code{coverage} = number of reads overlapping each site
#' }
#' Sites with \code{coverage < min_coverage} are dropped.
#'
#' @param file Character string. Path to the Megalodon per-read BED file.
#' @param sample_name Character string. Sample name (used in messages).
#' @param mod_type Character string. Modification type to assign to all sites
#'   (e.g., \code{"6mA"}). Megalodon files are modification-type-specific,
#'   so the type cannot be auto-detected from the file alone and must be
#'   supplied explicitly.
#' @param min_coverage Integer. Minimum read depth. Default \code{5}.
#'
#' @return A \code{data.frame} with columns: \code{chrom}, \code{position}
#'   (1-based), \code{strand}, \code{mod_type}, \code{motif} (always
#'   \code{NA} because Megalodon files do not encode motif context),
#'   \code{beta},
#'   \code{coverage}, \code{mod_counts}, \code{canonical_counts}, and
#'   \code{other_mod_counts}. Count columns are \code{NA} because this parser
#'   aggregates probabilities rather than observed calls.
#'
#' @keywords internal
.parseMegalodon <- function(file,
                            sample_name,
                            mod_type = NULL,
                            min_coverage = 5L) {
  if (!is.character(file) || length(file) != 1) {
    stop("file must be a single character string path")
  }
  if (!file.exists(file)) {
    stop("Megalodon file not found: ", file)
  }
  min_coverage <- as.integer(min_coverage)

  invalid_mod_type <- is.null(mod_type)
  invalid_mod_type <- invalid_mod_type || !is.character(mod_type)
  invalid_mod_type <- invalid_mod_type || length(mod_type) != 1L
  invalid_mod_type <- invalid_mod_type || is.na(mod_type)
  if (invalid_mod_type) {
    stop(
      "'mod_type' must be explicitly supplied as a single modification ",
      "type for Megalodon files. Megalodon output does not encode the ",
      "modification type in the file."
    )
  }
  mod_type_errors <- .checkModTypeValues(mod_type)
  if (length(mod_type_errors) > 0L) {
    stop(paste(mod_type_errors, collapse = "\n"))
  }

  # Read file.
  raw <- tryCatch(
    read.table(
      file,
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE,
      comment.char = "#",
      fill = TRUE
    ),
    error = function(e) {
      stop("Failed to read Megalodon file '", file, "': ", e$message)
    }
  )

  if (nrow(raw) == 0L) {
    return(.emptyParseResult())
  }

  if (ncol(raw) < 7L) {
    stop(
      "Megalodon file '", file, "' has ", ncol(raw), " columns; ",
      "expected at least 7 (chrom, start, end, read_id, score, ",
      "strand, mod_prob)."
    )
  }

  # Standard Megalodon per-read BED columns (minimum 7)
  # Col 1=chrom, 2=start, 3=end, 4=read_id, 5=score, 6=strand, last=mod_prob
  chrom <- as.character(raw[[1]])
  start <- .parseMegalodonNumericField(
    raw[[2]],
    file = file,
    field = "start coordinate",
    integer = TRUE,
    min = 0
  )
  .parseMegalodonNumericField(
    raw[[3]],
    file = file,
    field = "end coordinate",
    integer = TRUE,
    min = 0
  )
  strand <- as.character(raw[[6]])
  mod_prob <- .parseMegalodonNumericField(
    raw[[ncol(raw)]],
    file = file,
    field = "mod_prob",
    min = 0,
    max = 1
  )

  # Aggregate per-read calls to per-site values.
  # Group by genomic position (chrom, position, strand) and compute
  # per-site beta (mean of mod_prob) and coverage (count of reads).
  position <- start + 1L

  site_df <- data.frame(
    chrom = chrom,
    position = position,
    strand = strand,
    mod_type = mod_type,
    mod_prob = mod_prob,
    stringsAsFactors = FALSE
  )

  # Summarize each keyed site once so beta and coverage cannot drift by order.
  site_summary <- stats::aggregate(
    mod_prob ~ chrom + position + strand + mod_type,
    data = site_df,
    FUN = function(x) {
      c(beta = mean(x), coverage = length(x))
    }
  )
  summary_metrics <- site_summary$mod_prob

  result <- data.frame(
    chrom = site_summary$chrom,
    position = site_summary$position,
    strand = site_summary$strand,
    mod_type = site_summary$mod_type,
    motif = NA_character_,
    beta = summary_metrics[, "beta"],
    coverage = as.integer(summary_metrics[, "coverage"]),
    mod_counts = NA_integer_,
    canonical_counts = NA_integer_,
    other_mod_counts = NA_integer_,
    stringsAsFactors = FALSE
  )

  # Apply min_coverage filter.
  result <- result[result$coverage >= min_coverage, , drop = FALSE]
  rownames(result) <- NULL
  result
}

.parseMegalodonNumericField <- function(values,
                                        file,
                                        field,
                                        integer = FALSE,
                                        min = -Inf,
                                        max = Inf) {
  values_chr <- trimws(as.character(values))
  value_is_nan <- if (is.numeric(values)) {
    is.nan(values)
  } else {
    rep(FALSE, length(values_chr))
  }
  missing <- (is.na(values) & !value_is_nan) | is.na(values_chr) |
    !nzchar(values_chr)

  if (any(missing)) {
    row_idx <- which(missing)[[1L]]
    stop(
      "Megalodon file '", file, "' has missing ", field,
      " in row ", row_idx, "."
    )
  }

  parsed <- suppressWarnings(as.numeric(values_chr))
  not_finite <- !is.finite(parsed)
  if (any(not_finite)) {
    row_idx <- which(not_finite)[[1L]]
    stop(
      "Megalodon file '", file, "' has malformed ", field,
      " in row ", row_idx, ": '", values_chr[[row_idx]],
      "'. Expected a finite numeric value."
    )
  }

  if (integer) {
    not_integer <- parsed != floor(parsed)
    if (any(not_integer)) {
      row_idx <- which(not_integer)[[1L]]
      stop(
        "Megalodon file '", file, "' has malformed ", field,
        " in row ", row_idx, ": '", values_chr[[row_idx]],
        "'. Expected an integer value."
      )
    }
    parsed_integer <- suppressWarnings(as.integer(parsed))
    overflow <- is.na(parsed_integer) & !is.na(parsed)
    if (any(overflow)) {
      row_idx <- which(overflow)[[1L]]
      stop(
        "Megalodon file '", file, "' has out-of-range ", field,
        " in row ", row_idx, ": '", values_chr[[row_idx]],
        "'. Expected a value between ", min, " and ", max, "."
      )
    }
    parsed <- parsed_integer
  }

  out_of_range <- parsed < min | parsed > max
  if (any(out_of_range, na.rm = TRUE)) {
    row_idx <- which(out_of_range)[[1L]]
    stop(
      "Megalodon file '", file, "' has out-of-range ", field,
      " in row ", row_idx, ": '", values_chr[[row_idx]],
      "'. Expected a value between ", min, " and ", max, "."
    )
  }

  parsed
}
