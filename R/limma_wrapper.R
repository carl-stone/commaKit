#' @importFrom stats model.matrix relevel
NULL

# ─── limma eBayes wrapper ─────────────────────────────────────────────────────

#' Per-site moderated t-test via limma eBayes for differential methylation
#'
#' An internal wrapper that uses \pkg{limma}'s empirical Bayes moderated
#' t-test to identify differentially methylated sites. Called by
#' \code{\link{diffMethyl}} when \code{method = "limma"}.
#'
#' \pkg{limma} must be installed (it is listed in \code{Suggests}). If it
#' is not available, this function stops with an informative message.
#'
#' @details
#' Beta values are first transformed to M-values:
#' \deqn{M = \log_2\!\left(\frac{n_{\mathrm{mod}} + \alpha}{n_{\mathrm{unmod}} + \alpha}\right)}
#' where \eqn{\alpha} is a pseudocount (default 0.5). M-values are
#' approximately normally distributed and homoscedastic, making OLS
#' appropriate.
#'
#' A linear model is fitted across all sites simultaneously with
#' \code{\link[limma]{lmFit}}. \code{\link[limma]{eBayes}} then estimates an
#' empirical Bayes prior on the residual variance across all sites and computes
#' a moderated posterior variance per site — shrinking the noisy per-site
#' estimate toward the genome-wide mean. P-values are derived from a moderated
#' t-statistic with posterior degrees of freedom
#' \eqn{d_0 + df_{\mathrm{residual}}}.
#'
#' Only sites where all samples have non-NA M-values (i.e., non-zero coverage
#' after \code{min_coverage} thresholding) are passed to limma. Sites with any
#' \code{NA} retain \code{NA} in all result columns.
#'
#' Effect sizes (\code{delta_beta}) and per-group means are reported on the
#' original beta (0–1) scale for interpretability, not back-transformed from
#' M-value coefficients.
#'
#' @param methyl_mat Numeric matrix (sites × samples) of beta values.
#'   \code{NA} indicates below-coverage sites.
#' @param coverage_mat Integer matrix (sites × samples) of read depths.
#' @param mod_counts_mat Optional integer matrix of observed modified-read
#'   counts. If supplied, these counts are preferred over reconstructing from
#'   beta values.
#' @param canonical_counts_mat Optional integer matrix of observed
#'   canonical-read counts.
#' @param other_mod_counts_mat Optional integer matrix of observed non-target
#'   modified-read counts. When present with \code{canonical_counts_mat}, these
#'   counts are included in the non-target denominator so it matches coverage.
#' @param site_df Data frame with columns \code{chrom}, \code{position},
#'   \code{strand}, \code{mod_type}, \code{motif} — one row per site.
#' @param coldata \code{data.frame} with at least one column matching the
#'   RHS variable in \code{formula} (typically \code{condition}).
#' @param formula One-sided formula specifying the design (e.g.,
#'   \code{~ condition}).
#' @param alpha Positive numeric pseudocount added to modified and unmodified
#'   read counts before log-transformation. Default \code{0.5}.
#'
#' @return A \code{data.frame} with one row per site (same row order as
#'   \code{methyl_mat}), containing:
#'   \describe{
#'     \item{\code{pvalue}}{Moderated t-test p-value from \code{eBayes}.
#'       \code{NA} for sites with any missing data.}
#'     \item{\code{delta_beta}}{Effect size (treatment mean beta minus
#'       reference mean beta) on the 0–1 scale. \code{NA} where group means
#'       cannot be computed.}
#'     \item{\code{mean_beta_<level>}}{One column per condition level
#'       containing the per-group observed mean beta value.}
#'   }
#'
#' @keywords internal
.runLimma <- function(methyl_mat, coverage_mat, site_df, coldata, formula, alpha = 0.5,
                      ref_level = NULL, design_info = NULL,
                      mod_counts_mat = NULL, canonical_counts_mat = NULL,
                      other_mod_counts_mat = NULL) {
    # ── Dependency check ──────────────────────────────────────────────────────
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop(
            "Package 'limma' is required for method = \"limma\".\n",
            "Install it with: BiocManager::install(\"limma\")\n",
            "Alternatively, use method = \"methylkit\" if methylKit is available."
        )
    }

    # ── Validate alpha ────────────────────────────────────────────────────────
    if (!is.numeric(alpha) || length(alpha) != 1L ||
            !is.finite(alpha) || alpha <= 0) {
        stop("'alpha' must be a single positive finite number.")
    }

    # ── Resolve two-level design and group statistics ─────────────────────────
    if (is.null(design_info)) {
        design_info <- .resolveDiffMethylDesign(coldata, formula, ref_level = ref_level)
    }
    primary_var <- design_info$primary_var
    ref_level   <- design_info$ref_level
    treat_level <- design_info$treat_level
    cond_levels <- design_info$cond_levels
    cond        <- design_info$cond

    n_sites <- nrow(methyl_mat)

    group_stats    <- .computeDiffMethylGroupStats(methyl_mat, design_info)
    group_means    <- group_stats$group_means
    delta_beta_vec <- group_stats$delta_beta
    count_mats <- .resolveCountMatrices(
        methyl_mat,
        coverage_mat,
        mod_counts_mat = mod_counts_mat,
        canonical_counts_mat = canonical_counts_mat,
        other_mod_counts_mat = other_mod_counts_mat
    )

    # ── Compute M-value matrix ────────────────────────────────────────────────
    n_mod   <- count_mats$modified
    n_unmod <- count_mats$unmodified
    # Clamp to [0, coverage] to guard against floating-point edge cases
    n_mod   <- pmax(0, pmin(n_mod,   coverage_mat))
    n_unmod <- pmax(0, n_unmod)
    M_mat   <- log2((n_mod + alpha) / (n_unmod + alpha))
    dim(M_mat)      <- dim(coverage_mat)
    dimnames(M_mat) <- dimnames(coverage_mat)
    # Sites with zero or NA coverage → NA
    M_mat[is.na(methyl_mat) | is.na(coverage_mat) | coverage_mat == 0L] <- NA_real_

    # ── Identify complete-case sites ──────────────────────────────────────────
    complete_sites <- which(apply(!is.na(M_mat), 1L, all))
    pvalue_vec     <- rep(NA_real_, n_sites)

    if (length(complete_sites) < 2L) {
        # Not enough sites with complete data to estimate the eBayes prior
        result <- data.frame(
            pvalue     = pvalue_vec,
            delta_beta = delta_beta_vec,
            stringsAsFactors = FALSE
        )
        for (lv in cond_levels) {
            result[[paste0("mean_beta_", lv)]] <- group_means[, lv]
        }
        return(result)
    }

    M_complete <- M_mat[complete_sites, , drop = FALSE]

    # ── Build design matrix ───────────────────────────────────────────────────
    # Relevel the primary variable so model.matrix() encodes contrasts against
    # ref_level, regardless of whether the original column was a factor or char.
    coldata[[primary_var]] <- relevel(
        factor(coldata[[primary_var]]),
        ref = ref_level
    )
    design <- stats::model.matrix(formula, data = coldata)

    # ── Fit linear model + eBayes ─────────────────────────────────────────────
    fit <- limma::lmFit(M_complete, design)
    fit <- limma::eBayes(fit)

    # ── Extract p-values for the contrast coefficient ─────────────────────────
    coef_names    <- colnames(design)
    contrast_cols <- grep(primary_var, coef_names, value = TRUE)
    if (length(contrast_cols) == 0L) {
        stop(
            "Could not find a coefficient for '", primary_var,
            "' in the design matrix. ",
            "Available coefficients: ", paste(coef_names, collapse = ", ")
        )
    }
    # Take the last matching coefficient (mirrors .betaBinomialTest() behaviour)
    contrast_col <- contrast_cols[[length(contrast_cols)]]

    # fit$p.value is (complete sites) × (coefficients)
    pvalue_vec[complete_sites] <- fit$p.value[, contrast_col]

    # ── Assemble result ───────────────────────────────────────────────────────
    result <- data.frame(
        pvalue     = pvalue_vec,
        delta_beta = delta_beta_vec,
        stringsAsFactors = FALSE
    )
    for (lv in cond_levels) {
        result[[paste0("mean_beta_", lv)]] <- group_means[, lv]
    }

    result
}
