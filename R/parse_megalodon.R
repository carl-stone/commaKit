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
#'   \code{NA} — Megalodon files do not encode motif context), \code{beta},
#'   \code{coverage}, \code{mod_counts}, and \code{canonical_counts}. Count
#'   columns are \code{NA} because this parser aggregates probabilities rather
#'   than observed calls.
#'
#' @keywords internal
.parseMegalodon <- function(file, sample_name, mod_type = NULL, min_coverage = 5L) {
    if (!is.character(file) || length(file) != 1) {
        stop("file must be a single character string path")
    }
    if (!file.exists(file)) {
        stop("Megalodon file not found: ", file)
    }
    min_coverage <- as.integer(min_coverage)

    if (is.null(mod_type) || !is.character(mod_type) || length(mod_type) != 1L ||
            is.na(mod_type)) {
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

    # ── Read file ───────────────────────────────────────────────────────────
    raw <- tryCatch(
        read.table(
            file,
            header           = FALSE,
            sep              = "\t",
            stringsAsFactors = FALSE,
            comment.char     = "#",
            fill             = TRUE
        ),
        error = function(e) stop("Failed to read Megalodon file '", file, "': ", e$message)
    )

    if (nrow(raw) == 0L) {
        return(.emptyParseResult())
    }

    if (ncol(raw) < 7L) {
        stop(
            "Megalodon file '", file, "' has ", ncol(raw), " columns; ",
            "expected at least 7 (chrom, start, end, read_id, score, strand, mod_prob)."
        )
    }

    # Standard Megalodon per-read BED columns (minimum 7)
    # Col 1=chrom, 2=start, 3=end, 4=read_id, 5=score, 6=strand, last=mod_prob
    chrom    <- as.character(raw[[1]])
    start    <- as.integer(raw[[2]])
    strand   <- as.character(raw[[6]])
    mod_prob <- as.numeric(raw[[ncol(raw)]])

    # ── Aggregate per-read → per-site ───────────────────────────────────────
    # Group by genomic position (chrom, position, strand) and compute
    # per-site beta (mean of mod_prob) and coverage (count of reads).
    position <- start + 1L

    site_df <- data.frame(
        chrom    = chrom,
        position = position,
        strand   = strand,
        mod_prob = mod_prob,
        stringsAsFactors = FALSE
    )

    # Use aggregate() for base-R dedup — no string keys needed
    agg_beta <- stats::aggregate(mod_prob ~ chrom + position + strand,
                           data = site_df,
                           FUN = mean, na.rm = TRUE)
    agg_cov  <- stats::aggregate(mod_prob ~ chrom + position + strand,
                           data = site_df,
                           FUN = length)

    # Merge the two aggregates (same grouping columns, same row order)
    result <- data.frame(
        chrom    = agg_beta$chrom,
        position = agg_beta$position,
        strand   = agg_beta$strand,
        mod_type = mod_type,
        motif    = NA_character_,
        beta     = agg_beta$mod_prob,
        coverage = agg_cov$mod_prob,
        mod_counts = NA_integer_,
        canonical_counts = NA_integer_,
        stringsAsFactors = FALSE
    )

    # ── Apply min_coverage filter ───────────────────────────────────────────
    result <- result[result$coverage >= min_coverage, , drop = FALSE]
    rownames(result) <- NULL
    result
}
