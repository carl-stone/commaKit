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
#' @param mod_type Character string specifying a single modification type
#'   (e.g., \code{"6mA"}, \code{"5mC"}). If \code{NULL} (default), M-values
#'   are computed for all sites in \code{object}.
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
#' where \eqn{M_{\mathrm{reads}} = \mathrm{round}(\beta \times \mathrm{coverage})}
#' is the estimated number of methylated reads, \eqn{U_{\mathrm{reads}} =
#' \mathrm{coverage} - M_{\mathrm{reads}}} is the estimated number of
#' unmethylated reads, and \eqn{\alpha} is the pseudocount offset.
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
#' @seealso \code{\link{methylation}}, \code{\link[GenomicRanges]{coverage}},
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

    ## --- Optional mod_type filter -------------------------------------------
    if (!is.null(mod_type)) {
        available <- modTypes(object)
        if (!mod_type %in% available) {
            stop("'mod_type' = '", mod_type, "' not found in object. ",
                 "Available types: ", paste(available, collapse = ", "), ".")
        }
        object <- subset(object, mod_type = mod_type)
    }

    ## --- Optional motif filter ----------------------------------------------
    if (!is.null(motif)) {
        available_m <- motifs(object)
        bad_m <- setdiff(motif, available_m)
        if (length(bad_m) > 0L) {
            stop("'motif' value(s) not found in object: ",
                 paste(bad_m, collapse = ", "),
                 ". Available: ", paste(available_m, collapse = ", "), ".")
        }
        object <- subset(object, motif = motif)
    }

    ## --- Optional mod_context filter ----------------------------------------
    if (!is.null(mod_context)) {
        available_mc <- modContexts(object)
        bad_mc <- setdiff(mod_context, available_mc)
        if (length(bad_mc) > 0L) {
            stop("'mod_context' value(s) not found in object: ",
                 paste(bad_mc, collapse = ", "),
                 ". Available: ", paste(available_mc, collapse = ", "), ".")
        }
        object <- subset(object, mod_context = mod_context)
    }

    ## --- Compute M-values ---------------------------------------------------
    beta_mat <- methylation(object)   # sites × samples, values in [0, 1]
    cov_mat  <- coverage(object)      # sites × samples, non-negative integers

    ## Estimated methylated read counts (round to nearest integer)
    m_reads <- round(beta_mat * cov_mat)
    u_reads <- cov_mat - m_reads

    ## Sites with coverage == 0: set m_reads and u_reads to NA so that the
    ## log ratio returns NA rather than log2(alpha / alpha) = 0 (misleading).
    m_reads[cov_mat == 0L] <- NA_real_
    u_reads[cov_mat == 0L] <- NA_real_

    ## M-value = log2((M + alpha) / (U + alpha))
    log2((m_reads + alpha) / (u_reads + alpha))
}
