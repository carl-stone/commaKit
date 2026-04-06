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
.runLimma <- function(methyl_mat, coverage_mat, coldata, formula, alpha = 0.5,
                      ref_level = NULL) {
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

    # ── Parse formula ─────────────────────────────────────────────────────────
    rhs_vars <- all.vars(formula)
    if (length(rhs_vars) == 0L) {
        stop("'formula' must contain at least one RHS variable (e.g., ~ condition).")
    }
    primary_var <- rhs_vars[[1L]]

    if (!primary_var %in% colnames(coldata)) {
        stop(
            "Variable '", primary_var, "' from formula not found in sample metadata. ",
            "Available columns: ", paste(colnames(coldata), collapse = ", ")
        )
    }

    cond        <- as.character(coldata[[primary_var]])
    all_levels  <- sort(unique(cond))

    if (length(all_levels) < 2L) {
        stop(
            "Differential methylation requires at least 2 distinct levels of '",
            primary_var, "'. Found only: '", all_levels[[1L]], "'."
        )
    }

    # Use provided ref_level, or fall back to alphabetically first
    if (is.null(ref_level)) {
        ref_level <- all_levels[[1L]]
    }
    cond_levels <- c(ref_level, setdiff(all_levels, ref_level))
    treat_level <- cond_levels[[2L]]

    n_sites <- nrow(methyl_mat)

    # ── Pre-compute group means (beta scale, vectorised) ─────────────────────
    group_idx <- lapply(cond_levels, function(lv) which(cond == lv))
    names(group_idx) <- cond_levels

    group_means <- vapply(cond_levels, function(lv) {
        idx <- group_idx[[lv]]
        if (length(idx) == 1L) {
            methyl_mat[, idx]
        } else {
            rowMeans(methyl_mat[, idx, drop = FALSE], na.rm = TRUE)
        }
    }, numeric(n_sites))
    if (is.null(dim(group_means))) {
        group_means <- matrix(group_means, nrow = 1L,
                              dimnames = list(rownames(methyl_mat), cond_levels))
    }
    group_means[is.nan(group_means)] <- NA_real_

    delta_beta_vec <- group_means[, treat_level] - group_means[, ref_level]

    # ── Compute M-value matrix ────────────────────────────────────────────────
    n_mod   <- round(methyl_mat * coverage_mat)
    n_unmod <- coverage_mat - n_mod
    # Clamp to [0, coverage] to guard against floating-point edge cases
    n_mod   <- pmax(0, pmin(n_mod,   coverage_mat))
    n_unmod <- pmax(0, n_unmod)
    M_mat   <- log2((n_mod + alpha) / (n_unmod + alpha))
    dim(M_mat)      <- dim(coverage_mat)
    dimnames(M_mat) <- dimnames(coverage_mat)
    # Sites with zero or NA coverage → NA
    M_mat[is.na(coverage_mat) | coverage_mat == 0L] <- NA_real_

    # ── Identify complete-case sites ──────────────────────────────────────────
    complete_sites <- which(apply(!is.na(M_mat), 1L, all))
    pvalue_vec     <- rep(NA_real_, n_sites)

    if (length(complete_sites) < 2L) {
        # Not enough sites with complete data to estimate the eBayes prior
        result <- data.frame(
            pvalue     = pvalue_vec,
            delta_beta = delta_beta_vec,
            row.names  = rownames(methyl_mat),
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
        row.names  = rownames(methyl_mat),
        stringsAsFactors = FALSE
    )
    for (lv in cond_levels) {
        result[[paste0("mean_beta_", lv)]] <- group_means[, lv]
    }

    result
}
