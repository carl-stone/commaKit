#' @importFrom stats glm quasibinomial coef
NULL

# ─── Beta-binomial (quasibinomial GLM) test ───────────────────────────────────

#' Per-site quasibinomial GLM for differential methylation
#'
#' Fits a per-site quasibinomial generalized linear model to test for
#' differential methylation between conditions. This function is an internal
#' engine called by \code{\link{diffMethyl}} when \code{method = "beta_binomial"}.
#'
#' The model treats methylation counts as overdispersed binomial data:
#' for each site, the number of modified reads is modelled as
#' \eqn{\text{Binomial}(n_{\text{coverage}}, p)} with a quasibinomial
#' family to account for overdispersion. The formula is fitted site-by-site
#' using \code{\link[stats]{glm}}.
#'
#' @details
#' **Small-sample caveat:** With few replicates (e.g., 2 control + 1 treatment),
#' residual degrees of freedom will be very small (df = 1 in this case). The
#' quasibinomial dispersion estimate and the resulting p-values remain
#' mathematically valid but will have very low power. Interpret results from
#' n < 3 per group with caution.
#'
#' Sites are skipped (result is \code{NA}) if:
#' \itemize{
#'   \item All samples have \code{NA} beta values (low coverage).
#'   \item Fewer than 2 non-NA samples are present.
#'   \item The GLM fails to converge or the contrast coefficient cannot be
#'     extracted (e.g., singular fit, only one condition level present after
#'     NA removal).
#' }
#'
#' @param methyl_mat Numeric matrix (sites × samples) of beta values.
#'   \code{NA} indicates below-coverage sites.
#' @param coverage_mat Integer matrix (sites × samples) of read depths.
#' @param coldata \code{data.frame} with at least one column matching the
#'   RHS variable in \code{formula} (typically \code{condition}).
#' @param formula One-sided formula specifying the design (e.g.,
#'   \code{~ condition}).
#'
#' @return A \code{data.frame} with one row per site (same row order as
#'   \code{methyl_mat}), containing:
#'   \describe{
#'     \item{\code{pvalue}}{Raw p-value from the Wald test on the condition
#'       coefficient. \code{NA} for untestable sites.}
#'     \item{\code{delta_beta}}{Estimated effect size: mean beta in the
#'       treatment group minus mean beta in the reference (control) group.
#'       Uses the model matrix reference level as the baseline. \code{NA}
#'       for untestable sites.}
#'     \item{\code{mean_beta_<level>}}{One column per unique condition level,
#'       named \code{mean_beta_<level>} (e.g., \code{mean_beta_control},
#'       \code{mean_beta_treatment}). Contains the per-group observed mean
#'       beta value computed directly from non-NA data.}
#'   }
#'
#' @keywords internal
.betaBinomialTest <- function(methyl_mat, coverage_mat, coldata, formula,
                              ref_level = NULL) {
    # ── Parse formula ─────────────────────────────────────────────────────────
    rhs_vars <- all.vars(formula)
    if (length(rhs_vars) == 0L) {
        stop("'formula' must contain at least one RHS variable (e.g., ~ condition).")
    }
    # Use the first variable for the primary contrast
    primary_var <- rhs_vars[[1L]]

    if (!primary_var %in% colnames(coldata)) {
        stop(
            "Variable '", primary_var, "' from formula not found in sample metadata. ",
            "Available columns: ", paste(colnames(coldata), collapse = ", ")
        )
    }

    cond        <- as.character(coldata[[primary_var]])
    all_levels  <- unique(cond)

    if (length(all_levels) < 2L) {
        stop(
            "Differential methylation requires at least 2 distinct levels of '",
            primary_var, "'. Found only: '", all_levels, "'."
        )
    }

    # Use provided ref_level, or fall back to alphabetically first
    if (is.null(ref_level)) {
        ref_level <- sort(all_levels)[[1L]]
    }
    cond_levels <- c(ref_level, setdiff(sort(all_levels), ref_level))
    treat_level <- cond_levels[[2L]]

    n_sites   <- nrow(methyl_mat)
    n_samples <- ncol(methyl_mat)

    # ── Pre-compute mean beta per condition (vectorised) ─────────────────────
    # For each condition level, column indices belonging to that group
    group_idx <- lapply(cond_levels, function(lv) which(cond == lv))
    names(group_idx) <- cond_levels

    # per-group mean beta: site × n_levels matrix
    group_means <- vapply(cond_levels, function(lv) {
        idx <- group_idx[[lv]]
        if (length(idx) == 1L) {
            methyl_mat[, idx]
        } else {
            rowMeans(methyl_mat[, idx, drop = FALSE], na.rm = TRUE)
        }
    }, numeric(n_sites))
    # Ensure result is always a matrix (vapply returns a vector when n_sites == 1)
    if (is.null(dim(group_means))) {
        group_means <- matrix(group_means, nrow = 1L,
                              dimnames = list(rownames(methyl_mat), cond_levels))
    }
    # group_means is n_sites × n_levels; NaN where all samples are NA
    group_means[is.nan(group_means)] <- NA_real_

    # delta_beta = treat - ref (if exactly 2 groups); or treat1 - ref for >2
    delta_beta_vec <- group_means[, treat_level] - group_means[, ref_level]

    # ── Per-site GLM ─────────────────────────────────────────────────────────
    pvalue_vec <- numeric(n_sites)
    pvalue_vec[] <- NA_real_

    for (i in seq_len(n_sites)) {
        beta_i <- methyl_mat[i, ]
        cov_i  <- coverage_mat[i, ]

        # Require at least 2 non-NA samples with positive coverage
        ok <- !is.na(beta_i) & !is.na(cov_i) & cov_i > 0L
        if (sum(ok) < 2L) next

        # Require at least 2 distinct condition levels among non-NA samples
        cond_ok <- cond[ok]
        if (length(unique(cond_ok)) < 2L) next

        n_mod   <- round(beta_i[ok] * cov_i[ok])
        n_unmod <- cov_i[ok] - n_mod

        # Clamp to [0, coverage] to guard against floating-point edge cases
        n_mod   <- pmax(0L, pmin(n_mod, cov_i[ok]))
        n_unmod <- pmax(0L, n_unmod)

        df_glm <- data.frame(
            n_mod   = n_mod,
            n_unmod = n_unmod,
            stringsAsFactors = FALSE
        )
        # Set factor levels so GLM encodes contrasts against ref_level
        df_glm[[primary_var]] <- factor(
            cond_ok,
            levels = c(ref_level, setdiff(unique(cond_ok), ref_level))
        )

        fit <- tryCatch(
            glm(
                cbind(n_mod, n_unmod) ~ .,
                data   = df_glm,
                family = quasibinomial()
            ),
            error   = function(e) NULL,
            warning = function(w) {
                tryCatch(
                    suppressWarnings(glm(
                        cbind(n_mod, n_unmod) ~ .,
                        data   = df_glm,
                        family = quasibinomial()
                    )),
                    error = function(e2) NULL
                )
            }
        )

        if (is.null(fit)) next
        if (fit$df.residual < 1L) next

        cs <- tryCatch(coef(summary(fit)), error = function(e) NULL)
        if (is.null(cs)) next

        # Find the row corresponding to the contrast term (condition level)
        row_nm <- rownames(cs)
        contrast_row <- grep(primary_var, row_nm, value = TRUE)
        if (length(contrast_row) == 0L) next

        # Take the last matching row (highest non-reference level by default)
        p_val <- cs[contrast_row[[length(contrast_row)]], "Pr(>|t|)"]
        if (!is.na(p_val)) {
            pvalue_vec[[i]] <- p_val
        }
    }

    # ── Assemble result ───────────────────────────────────────────────────────
    result <- data.frame(
        pvalue     = pvalue_vec,
        delta_beta = delta_beta_vec,
        row.names  = rownames(methyl_mat),
        stringsAsFactors = FALSE
    )

    # Add per-group mean beta columns
    for (lv in sort(cond_levels)) {
        col_nm <- paste0("mean_beta_", lv)
        result[[col_nm]] <- group_means[, lv]
    }

    result
}
