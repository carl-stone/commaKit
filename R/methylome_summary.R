#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData assay
NULL

#' Summarize per-sample methylation and coverage distributions
#'
#' Computes per-sample summary statistics for methylation beta values and
#' sequencing coverage in a \code{\link{commaData}} object. Returns a tidy
#' \code{data.frame} suitable for direct use with \pkg{ggplot2} or for
#' tabular reporting.
#'
#' @param object A \code{\link{commaData}} object.
#' @param mod_type Character string or \code{NULL}. If provided, only sites
#'   of the specified modification type (e.g., \code{"6mA"}) are included in
#'   the summary. If \code{NULL} (default), all modification types are
#'   summarized together.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   matching sequence context motif(s) are included (e.g., \code{"GATC"}).
#'   If \code{NULL} (default), all motifs are included.
#' @param mod_context Character vector or \code{NULL}. If provided, only sites
#'   with a matching modification context are included (e.g.,
#'   \code{"6mA_GATC"}). Applied after any \code{mod_type} and \code{motif}
#'   filters. Use \code{\link{modContexts}} to see available values.
#'
#' @return A \code{data.frame} with one row per sample, containing:
#'   \describe{
#'     \item{\code{sample_name}}{Sample identifier.}
#'     \item{\code{condition}}{Experimental condition, from
#'       \code{sampleInfo(object)$condition}.}
#'     \item{\code{mod_type}}{The modification type summarized
#'       (\code{"all"} if \code{mod_type = NULL}).}
#'     \item{\code{n_sites}}{Total number of sites considered.}
#'     \item{\code{n_covered}}{Number of sites with non-\code{NA} methylation
#'       in this sample (i.e., sites above the coverage threshold).}
#'     \item{\code{mean_beta}}{Mean beta value across covered sites.}
#'     \item{\code{median_beta}}{Median beta value across covered sites.}
#'     \item{\code{sd_beta}}{Standard deviation of beta values across covered
#'       sites.}
#'     \item{\code{frac_methylated}}{Fraction of covered sites with
#'       \eqn{\beta > 0.5} (broadly methylated).}
#'     \item{\code{mean_coverage}}{Mean sequencing depth across all sites
#'       (including sites below the \code{min_coverage} threshold, which have
#'       coverage stored as 0 or their raw depth).}
#'     \item{\code{median_coverage}}{Median sequencing depth.}
#'   }
#'
#' @examples
#' data(comma_example_data)
#' ms <- methylomeSummary(comma_example_data)
#' ms
#'
#' # Summarize only 6mA sites
#' ms_6mA <- methylomeSummary(comma_example_data, mod_type = "6mA")
#' ms_6mA[, c("sample_name", "condition", "mean_beta", "n_covered")]
#'
#' @seealso \code{\link{methylation}}, \code{\link{coverage}},
#'   \code{\link{sampleInfo}}
#'
#' @export
methylomeSummary <- function(object, mod_type = NULL, motif = NULL,
                             mod_context = NULL) {
    # ── Input validation ──────────────────────────────────────────────────────
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }

    # ── Filter by mod_type ────────────────────────────────────────────────────
    mt_label <- if (is.null(mod_type)) "all" else mod_type

    if (!is.null(mod_type)) {
        available <- modTypes(object)
        if (!mod_type %in% available) {
            stop(
                "'mod_type' = '", mod_type, "' not found in object. ",
                "Available types: ", paste(available, collapse = ", ")
            )
        }
        object <- subset(object, mod_type = mod_type)
    }

    # ── Filter by motif ───────────────────────────────────────────────────────
    if (!is.null(motif)) {
        available_m <- motifs(object)
        bad_m <- setdiff(motif, available_m)
        if (length(bad_m) > 0L) {
            stop(
                "'motif' value(s) not found in object: ",
                paste(bad_m, collapse = ", "),
                ". Available: ", paste(available_m, collapse = ", ")
            )
        }
        object <- subset(object, motif = motif)
    }

    # ── Filter by mod_context ─────────────────────────────────────────────────
    if (!is.null(mod_context)) {
        available_mc <- modContexts(object)
        bad_mc <- setdiff(mod_context, available_mc)
        if (length(bad_mc) > 0L) {
            stop(
                "'mod_context' value(s) not found in object: ",
                paste(bad_mc, collapse = ", "),
                ". Available: ", paste(available_mc, collapse = ", ")
            )
        }
        object <- subset(object, mod_context = mod_context)
    }

    methyl_mat <- methylation(object)
    cov_mat    <- coverage(object)
    si         <- sampleInfo(object)
    sample_nms <- colnames(methyl_mat)
    n_sites    <- nrow(methyl_mat)

    # ── Per-sample statistics ─────────────────────────────────────────────────
    rows <- lapply(sample_nms, function(samp) {
        betas <- methyl_mat[, samp]
        covs  <- cov_mat[, samp]

        covered   <- !is.na(betas)
        n_covered <- sum(covered)
        b_cov     <- betas[covered]
        c_all     <- as.numeric(covs)

        data.frame(
            sample_name      = samp,
            condition        = si$condition[si$sample_name == samp],
            mod_type         = mt_label,
            n_sites          = n_sites,
            n_covered        = n_covered,
            mean_beta        = if (n_covered > 0) mean(b_cov)    else NA_real_,
            median_beta      = if (n_covered > 0) stats::median(b_cov) else NA_real_,
            sd_beta          = if (n_covered > 1) stats::sd(b_cov)     else NA_real_,
            frac_methylated  = if (n_covered > 0) mean(b_cov > 0.5)    else NA_real_,
            mean_coverage    = mean(c_all, na.rm = TRUE),
            median_coverage  = stats::median(c_all, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, rows)
}
