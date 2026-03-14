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
#' @param mod_type Character string or \code{NULL}. Modification type to assign
#'   to all sites (e.g., \code{"6mA"}). Megalodon files are modification-
#'   type-specific, so the type cannot be auto-detected from the file alone.
#'   If \code{NULL}, defaults to \code{"6mA"} with a warning.
#' @param min_coverage Integer. Minimum read depth. Default \code{5}.
#'
#' @return A \code{data.frame} with columns: \code{chrom}, \code{position}
#'   (1-based), \code{strand}, \code{mod_type}, \code{beta}, \code{coverage}.
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

    if (is.null(mod_type)) {
        warning(
            "mod_type not specified for Megalodon file '", file, "'. ",
            "Defaulting to '6mA'. Set mod_type explicitly in commaData()."
        )
        mod_type <- "6mA"
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
        return(.emptyModkitResult())
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
    # Site key: chrom:position:strand (1-based)
    position <- start + 1L
    site_key <- paste(chrom, position, strand, sep = ":")

    beta_vals <- tapply(mod_prob, site_key, mean,  na.rm = TRUE)
    cov_vals  <- tapply(mod_prob, site_key, length)

    all_keys   <- names(beta_vals)
    key_parts  <- strsplit(all_keys, ":", fixed = TRUE)
    out_chrom  <- vapply(key_parts, `[[`, character(1), 1L)
    out_pos    <- as.integer(vapply(key_parts, `[[`, character(1), 2L))
    out_strand <- vapply(key_parts, `[[`, character(1), 3L)

    result <- data.frame(
        chrom    = out_chrom,
        position = out_pos,
        strand   = out_strand,
        mod_type = mod_type,
        beta     = as.numeric(beta_vals),
        coverage = as.integer(cov_vals),
        stringsAsFactors = FALSE
    )

    # ── Apply min_coverage filter ───────────────────────────────────────────
    result <- result[result$coverage >= min_coverage, , drop = FALSE]
    rownames(result) <- NULL
    result
}
