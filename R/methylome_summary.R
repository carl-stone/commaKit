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
#' @param mod_type Character vector or \code{NULL}. If provided, only sites
#'   of the specified modification type(s) (e.g., \code{"6mA"},
#'   \code{c("6mA", "5mC")}) are included in the summary. If \code{NULL}
#'   (default), all modification types are summarized together.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   matching sequence context motif(s) are included (e.g., \code{"GATC"}).
#'   If \code{NULL} (default), all motifs are included.
#' @param mod_context Character vector or \code{NULL}. If provided, only sites
#'   with a matching modification context are included (e.g.,
#'   \code{"6mA:GATC"}). Applied after any \code{mod_type} and \code{motif}
#'   filters. Use \code{\link{modContexts}} to see available values.
#'
#' @details
#' Methylation beta summaries (\code{mean_beta}, \code{median_beta},
#' \code{sd_beta}, and \code{frac_methylated}) are computed over covered sites:
#' sites with non-\code{NA} beta values after the object's coverage threshold
#' has been applied. Coverage summaries (\code{mean_coverage} and
#' \code{median_coverage}) are computed over non-missing coverage values for
#' sites retained after filtering, including sites whose beta value is
#' \code{NA} in a sample because the site did not meet the coverage threshold.
#'
#' @return A \code{data.frame} with one row per sample, containing:
#'   \describe{
#'     \item{\code{sample_name}}{Sample identifier.}
#'     \item{\code{condition}}{Experimental condition, from
#'       the optional \code{condition} column in \code{sampleInfo(object)}, or
#'       \code{NA} when that metadata is absent.}
#'     \item{\code{mod_type}}{The modification type summarized
#'       (\code{"all"} if \code{mod_type = NULL}).}
#'     \item{\code{n_sites}}{Total number of sites considered after filters.}
#'     \item{\code{n_covered}}{Number of sites with non-\code{NA} methylation
#'       in this sample (i.e., sites above the coverage threshold); this is the
#'       denominator for beta summaries.}
#'     \item{\code{mean_beta}}{Mean beta value across covered sites.}
#'     \item{\code{median_beta}}{Median beta value across covered sites.}
#'     \item{\code{sd_beta}}{Standard deviation of beta values across covered
#'       sites.}
#'     \item{\code{frac_methylated}}{Fraction of covered sites with
#'       \eqn{\beta > 0.5} (broadly methylated).}
#'     \item{\code{mean_coverage}}{Mean sequencing depth across non-missing
#'       coverage values for retained sites, including sites below the
#'       \code{min_coverage} threshold when coverage is available.}
#'     \item{\code{median_coverage}}{Median sequencing depth across non-missing
#'       coverage values for retained sites, including sites below the
#'       \code{min_coverage} threshold when coverage is available.}
#'     \item{\code{caller}}{Methylation caller that produced the data
#'       (e.g., \code{"modkit"}), or \code{NA} if not stored.}
#'     \item{\code{min_coverage}}{Minimum coverage threshold applied at
#'       construction, or \code{NA} if not stored.}
#'   }
#'
#' @examples
#' data(comma_example_data)
#' ms <- methylomeSummary(comma_example_data)
#' ms
#' # Beta summaries use n_covered; coverage summaries use non-missing coverage.
#' ms[, c("sample_name", "n_sites", "n_covered", "mean_beta", "mean_coverage")]
#'
#' # Summarize only 6mA sites
#' ms_6mA <- methylomeSummary(comma_example_data, mod_type = "6mA")
#' ms_6mA[, c("sample_name", "condition", "mean_beta", "n_covered")]
#'
#' @seealso \code{\link{methylation}}, \code{\link{siteCoverage}},
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
    mt_label <- if (is.null(mod_type)) "all" else paste(mod_type, collapse = ",")

    object <- .applySiteFilters(
        object,
        mod_type = mod_type,
        motif = motif,
        mod_context = mod_context,
        caller = "methylomeSummary()"
    )

    methyl_mat <- methylation(object)
    cov_mat    <- siteCoverage(object)
    si         <- sampleInfo(object)
    sample_nms <- colnames(methyl_mat)
    n_sites    <- nrow(methyl_mat)
    obj_caller <- caller(object)
    obj_min_cov <- minCoverage(object)

    # ── Per-sample statistics ─────────────────────────────────────────────────
    rows <- lapply(sample_nms, function(samp) {
        betas <- methyl_mat[, samp]
        covs  <- cov_mat[, samp]
        condition_value <- if ("condition" %in% colnames(si)) {
            si$condition[si$sample_name == samp]
        } else {
            NA_character_
        }

        covered   <- !is.na(betas)
        n_covered <- sum(covered)
        b_cov     <- betas[covered]
        c_all     <- as.numeric(covs)

        data.frame(
            sample_name      = samp,
            condition        = condition_value,
            mod_type         = mt_label,
            n_sites          = n_sites,
            n_covered        = n_covered,
            mean_beta        = if (n_covered > 0) mean(b_cov)    else NA_real_,
            median_beta      = if (n_covered > 0) stats::median(b_cov) else NA_real_,
            sd_beta          = if (n_covered > 1) stats::sd(b_cov)     else NA_real_,
            frac_methylated  = if (n_covered > 0) mean(b_cov > 0.5)    else NA_real_,
            mean_coverage    = mean(c_all, na.rm = TRUE),
            median_coverage  = stats::median(c_all, na.rm = TRUE),
            caller           = obj_caller,
            min_coverage     = obj_min_cov,
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, rows)
}
