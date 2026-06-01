#' Compute M-values from a commaData object
#'
#' Converts per-site beta values (methylation fractions) and read depths into
#' M-values using a pseudocount-offset logit transformation. M-values are
#' variance-stabilized relative to beta values and are better suited for
#' distance-based analyses such as PCA or hierarchical clustering.
#'
#' @param object A \code{\link{commaData}} object containing methylation beta
#'   values and coverage (read depth) assays.
#' @param alpha Positive numeric scalar. Pseudocount added to both the
#'   methylated and unmethylated read counts before log-transformation. Prevents
#'   infinite values at beta = 0 or beta = 1 and corresponds to a symmetric
#'   Beta(alpha, alpha) prior on the methylation fraction. Default \code{0.5}.
#' @param mod_type Character vector of modification types to include
#'   (e.g., \code{"6mA"}, \code{c("6mA", "5mC")}). If \code{NULL} (default),
#'   M-values are computed for all sites in \code{object}.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   matching sequence context motif(s) are included. If \code{NULL} (default),
#'   all motifs are included.
#' @param mod_context Character vector or \code{NULL}. If provided, only sites
#'   with a matching modification context are included (e.g.,
#'   \code{"6mA_GATC"}). Applied after any \code{mod_type} and \code{motif}
#'   filters.
#'
#' @details
#' The M-value for a site in one sample is computed as:
#'
#' \deqn{M = \log_2\!\left(\frac{M_{\mathrm{reads}} + \alpha}{U_{\mathrm{reads}} + \alpha}\right)}
#'
#' where \eqn{M_{\mathrm{reads}}} and \eqn{U_{\mathrm{reads}}} are observed
#' modified and canonical/unmodified read counts when count assays are
#' available. For legacy or probability-only objects, counts are reconstructed
#' from \code{round(beta * coverage)} as a documented fallback. The
#' pseudocount offset is \eqn{\alpha}.
#'
#' Sites with zero coverage or \code{NA} beta values are returned as \code{NA}.
#' The pseudocount \code{alpha} must be strictly positive to avoid \code{-Inf}
#' or \code{NaN} values in the output.
#'
#' @return A numeric matrix of M-values with the same dimensions and
#'   \code{dimnames} as \code{methylation(object)} (or the subset of rows
#'   matching \code{mod_type} if specified). Positive values indicate
#'   hypermethylation; negative values indicate hypomethylation; zero corresponds
#'   to a beta value of approximately 0.5.
#'
#' @examples
#' data(comma_example_data)
#'
#' # Compute M-values for all modification types
#' m <- mValues(comma_example_data)
#' dim(m)          # same as dim(methylation(comma_example_data))
#' range(m, na.rm = TRUE)
#'
#' # Only 6mA sites
#' m6 <- mValues(comma_example_data, mod_type = "6mA")
#'
#' # Use a smaller pseudocount
#' m_tight <- mValues(comma_example_data, alpha = 0.1)
#'
#' @seealso \code{\link{methylation}}, \code{\link{siteCoverage}},
#'   \code{\link{modCounts}}, \code{\link{canonicalCounts}},
#'   \code{\link{plot_pca}}
#'
#' @export
mValues <- function(object, alpha = 0.5, mod_type = NULL, motif = NULL,
                    mod_context = NULL) {

    ## --- Input validation ---------------------------------------------------
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
            alpha <= 0) {
        stop("'alpha' must be a single positive finite number.")
    }

    ## --- Optional site filters ----------------------------------------------
    object <- .applySiteFilters(
        object,
        mod_type = mod_type,
        motif = motif,
        mod_context = mod_context,
        caller = "mValues()"
    )

    ## --- Compute M-values ---------------------------------------------------
    beta_mat <- methylation(object)   # sites × samples, values in [0, 1]
    cov_mat  <- siteCoverage(object)      # sites × samples, non-negative integers
    count_mats <- .resolveCountMatrices(
        beta_mat,
        cov_mat,
        mod_counts_mat = .optionalAssay(object, "mod_counts"),
        canonical_counts_mat = .optionalAssay(object, "canonical_counts")
    )

    ## Clamp to the physically possible range so malformed/manual objects cannot
    ## produce negative unmethylated counts and NaN M-values.
    m_reads <- count_mats$modified
    m_reads <- pmax(0, pmin(m_reads, cov_mat))
    u_reads <- pmax(0, count_mats$unmodified)
    dim(m_reads) <- dim(beta_mat)
    dim(u_reads) <- dim(beta_mat)
    dimnames(m_reads) <- dimnames(beta_mat)
    dimnames(u_reads) <- dimnames(beta_mat)

    ## Sites with coverage == 0 or filtered beta values return NA rather than a
    ## pseudocount-only log ratio. Count assays may preserve raw observations
    ## even when the beta layer has been filtered by min_coverage.
    missing_value <- is.na(beta_mat) | is.na(cov_mat) | cov_mat == 0L
    m_reads[missing_value] <- NA_real_
    u_reads[missing_value] <- NA_real_

    ## M-value = log2((M + alpha) / (U + alpha))
    log2((m_reads + alpha) / (u_reads + alpha))
}
