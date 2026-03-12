NULL

# ─── mod_code → mod_type mapping ─────────────────────────────────────────────

.MODKIT_CODE_MAP <- c(
    "a"     = "6mA",
    "m"     = "5mC",
    "21839" = "4mC"
)

# modkit pileup BED column names (15 columns)
.MODKIT_COLS <- c(
    "chrom", "start", "end", "mod_code", "score", "strand",
    "coverage", "mod_frequency", "n_mod", "n_canonical",
    "n_other_mod", "n_delete", "n_fail", "n_diff", "n_no_call"
)

#' Parse a modkit pileup BED file into a tidy per-site data frame
#'
#' Reads a single-sample modkit \code{pileup} output file and returns a tidy
#' data frame of per-site methylation values. This is an internal function
#' called by \code{\link{commaData}}.
#'
#' @param file Character string. Path to the modkit pileup BED file.
#' @param sample_name Character string. Name for this sample (used in messages
#'   only; not added to the returned data frame).
#' @param mod_type Character vector or \code{NULL}. If provided, only sites
#'   with a matching modification type (e.g., \code{"6mA"}, \code{"5mC"},
#'   \code{"4mC"}) are returned. \code{NULL} retains all modification types.
#' @param min_coverage Integer. Minimum read depth required to retain a site.
#'   Sites with \code{coverage < min_coverage} are dropped. Default \code{5}.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{\code{chrom}}{Chromosome name (character).}
#'     \item{\code{position}}{1-based genomic position (integer).}
#'     \item{\code{strand}}{Strand, \code{"+"} or \code{"-"} (character).}
#'     \item{\code{mod_type}}{Modification type: \code{"6mA"}, \code{"5mC"},
#'       or \code{"4mC"} (character).}
#'     \item{\code{beta}}{Proportion of reads called methylated, range 0–1
#'       (numeric).}
#'     \item{\code{coverage}}{Total read depth at this site (integer).}
#'   }
#'
#' @keywords internal
.parseModkit <- function(file, sample_name, mod_type = NULL, min_coverage = 5L) {
    # ── Validate inputs ─────────────────────────────────────────────────────
    if (!is.character(file) || length(file) != 1) {
        stop("file must be a single character string path")
    }
    if (!file.exists(file)) {
        stop("modkit BED file not found: ", file)
    }
    min_coverage <- as.integer(min_coverage)

    # ── Read file ───────────────────────────────────────────────────────────
    raw <- tryCatch(
        read.table(
            file,
            header          = FALSE,
            sep             = "\t",
            stringsAsFactors = FALSE,
            comment.char    = "#",
            fill            = TRUE
        ),
        error = function(e) {
            if (grepl("no lines available in input", e$message, fixed = TRUE)) {
                return(NULL)  # empty file
            }
            stop("Failed to read modkit BED file '", file, "': ", e$message)
        }
    )

    if (is.null(raw) || nrow(raw) == 0L) {
        message("Note: modkit BED file '", file, "' contains no data rows")
        return(.emptyModkitResult())
    }

    if (ncol(raw) < 15L) {
        stop(
            "modkit BED file '", file, "' has ", ncol(raw), " columns; ",
            "expected at least 15 (modkit pileup format). ",
            "Check that the file is a modkit pileup output."
        )
    }
    # Use only the first 15 columns
    raw <- raw[, seq_len(15L), drop = FALSE]
    colnames(raw) <- .MODKIT_COLS

    # ── Map mod_code → mod_type ─────────────────────────────────────────────
    raw$mod_code <- as.character(raw$mod_code)
    mapped        <- .MODKIT_CODE_MAP[raw$mod_code]
    unknown_codes <- unique(raw$mod_code[is.na(mapped)])
    if (length(unknown_codes) > 0) {
        warning(
            "Unknown mod_code values in '", file, "' (skipped): ",
            paste(unknown_codes, collapse = ", "),
            ". Known codes: ", paste(names(.MODKIT_CODE_MAP), collapse = ", ")
        )
    }
    raw$mod_type_mapped <- mapped

    # Drop rows with unrecognized mod_code
    raw <- raw[!is.na(raw$mod_type_mapped), , drop = FALSE]

    # ── Apply min_coverage filter ───────────────────────────────────────────
    raw$coverage <- as.integer(raw$coverage)
    raw <- raw[raw$coverage >= min_coverage, , drop = FALSE]

    # ── Apply mod_type filter ───────────────────────────────────────────────
    if (!is.null(mod_type)) {
        raw <- raw[raw$mod_type_mapped %in% mod_type, , drop = FALSE]
    }

    if (nrow(raw) == 0L) {
        return(.emptyModkitResult())
    }

    # ── Build result ────────────────────────────────────────────────────────
    data.frame(
        chrom    = as.character(raw$chrom),
        position = as.integer(raw$start) + 1L,  # BED is 0-based → 1-based
        strand   = as.character(raw$strand),
        mod_type = raw$mod_type_mapped,
        beta     = as.numeric(raw$mod_frequency),
        coverage = raw$coverage,
        stringsAsFactors = FALSE
    )
}

#' Empty modkit parse result (zero-row data frame with correct schema)
#' @keywords internal
.emptyModkitResult <- function() {
    data.frame(
        chrom    = character(0),
        position = integer(0),
        strand   = character(0),
        mod_type = character(0),
        beta     = numeric(0),
        coverage = integer(0),
        stringsAsFactors = FALSE
    )
}
