#' @importFrom stats glm quasibinomial coef pt
NULL

# ─── Quasi-likelihood F-test (EB on quasibinomial dispersions) ────────────────

#' Per-site quasi-likelihood F-test for differential methylation
#'
#' An internal wrapper that combines the per-site quasibinomial GLM of
#' \code{.betaBinomialTest} with empirical Bayes shrinkage of the per-site
#' dispersion estimates. Called by \code{\link{diffMethyl}} when
#' \code{method = "quasi_f"}.
#'
#' \pkg{limma} must be installed (it is listed in \code{Suggests}).
#'
#' @details
#' The method runs in three passes:
#'
#' \strong{Pass 1 — per-site GLM.}
#' The same quasibinomial model as \code{.betaBinomialTest} is fitted at each
#' site:
#' \deqn{\mathrm{glm}(\mathrm{cbind}(n_{\mathrm{mod}},\, n_{\mathrm{unmod}})
#'   \sim \mathrm{condition},\; \mathrm{family} = \mathrm{quasibinomial}())}
#' For each site \eqn{j}, three quantities are collected:
#' \itemize{
#'   \item \eqn{\hat\phi_j = }\code{fit\$dispersion} — Pearson chi-squared
#'     dispersion estimate
#'   \item \eqn{df_j = }\code{fit\$df.residual} — residual degrees of freedom
#'   \item \eqn{\tilde{t}_j^{(0)} = t_j \times \sqrt{\hat\phi_j}} — the
#'     "unscaled" t-statistic (independent of \eqn{\hat\phi_j}), where
#'     \eqn{t_j} is the Wald t-statistic from \code{coef(summary(fit))}
#' }
#'
#' \strong{Pass 2 — empirical Bayes dispersion shrinkage.}
#' \code{\link[limma]{squeezeVar}} pools the \eqn{\{\hat\phi_j\}} estimates
#' across all testable sites, fits a log-normal prior, and returns posterior
#' dispersion estimates \eqn{\{\tilde\phi_j\}} and a prior degrees-of-freedom
#' scalar \eqn{d_0}.
#'
#' \strong{Pass 3 — moderated test statistic.}
#' The posterior t-statistic and p-value for each site are:
#' \deqn{\tilde{t}_j = \frac{\tilde{t}_j^{(0)}}{\sqrt{\tilde\phi_j}},
#'   \quad p_j = 2\,P(T \leq -|\tilde{t}_j|),\;
#'   T \sim t(d_0 + df_j)}
#' The additional \eqn{d_0} degrees of freedom are the power gain over the
#' unadjusted quasibinomial test.
#'
#' This procedure is methodologically equivalent to the quasi-likelihood F-test
#' of \pkg{edgeR} (\code{glmQLFTest}), adapted for methylation proportions
#' (quasibinomial) rather than RNA-seq counts (quasi-negative-binomial).
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
#'     \item{\code{pvalue}}{Moderated quasi-F p-value. \code{NA} for
#'       untestable sites.}
#'     \item{\code{delta_beta}}{Effect size (treatment mean beta minus
#'       reference mean beta) on the 0–1 scale.}
#'     \item{\code{mean_beta_<level>}}{One column per condition level
#'       containing the per-group observed mean beta value.}
#'   }
#'
#' @keywords internal
.runQuasiF <- function(methyl_mat, coverage_mat, coldata, formula,
                       ref_level = NULL) {
    # ── Dependency check ──────────────────────────────────────────────────────
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop(
            "Package 'limma' is required for method = \"quasi_f\".\n",
            "Install it with: BiocManager::install(\"limma\")\n",
            "Alternatively, use method = \"methylkit\" if methylKit is available."
        )
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

    n_sites   <- nrow(methyl_mat)
    n_samples <- ncol(methyl_mat)

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

    # ── Pass 1: per-site GLM — collect dispersion and unscaled t ─────────────
    phi_vec    <- rep(NA_real_, n_sites)
    df_vec     <- rep(NA_integer_, n_sites)
    t_unscaled <- rep(NA_real_, n_sites)

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

        # Clamp to [0, coverage]
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

        if (is.null(fit) || fit$df.residual < 1L) next

        # summary(fit) computes the Pearson dispersion estimate for quasibinomial;
        # fit$dispersion is NULL for quasi families — must use summary(fit)$dispersion
        sm <- tryCatch(summary(fit), error = function(e) NULL)
        if (is.null(sm)) next

        phi_hat <- sm$dispersion
        if (is.null(phi_hat) || length(phi_hat) != 1L ||
                is.na(phi_hat) || phi_hat <= 0) next

        cs <- sm$coefficients
        if (is.null(cs)) next

        row_nm       <- rownames(cs)
        contrast_row <- grep(primary_var, row_nm, value = TRUE)
        if (length(contrast_row) == 0L) next

        cr <- contrast_row[[length(contrast_row)]]

        # t_j is the Wald t-statistic from the GLM (uses phi_hat in SE)
        # t_unscaled_j = t_j × sqrt(phi_hat) = beta_hat / unscaled_SE
        # This is independent of phi_hat and is what we carry forward.
        t_j <- cs[cr, "t value"]
        if (is.na(t_j)) next

        phi_vec[i]    <- phi_hat
        df_vec[i]     <- fit$df.residual
        t_unscaled[i] <- t_j * sqrt(phi_hat)
    }

    # ── Pass 2: EB shrinkage on dispersions via limma::squeezeVar ────────────
    testable   <- !is.na(phi_vec) & !is.na(t_unscaled)
    pvalue_vec <- rep(NA_real_, n_sites)

    if (sum(testable) < 2L) {
        # Not enough sites to estimate the EB prior; return all NA
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

    squeezed <- limma::squeezeVar(phi_vec[testable], df_vec[testable])
    phi_post <- squeezed$var.post    # posterior dispersions (same length as phi_vec[testable])
    df_prior <- squeezed$df.prior    # scalar: estimated prior degrees of freedom

    # Guard against non-finite df_prior (degenerate case: all dispersions identical)
    if (is.na(df_prior) || !is.finite(df_prior)) {
        df_prior <- 0
    }

    # ── Pass 3: recompute t-stats and p-values with posterior dispersion ──────
    t_post  <- t_unscaled[testable] / sqrt(phi_post)
    df_post <- df_prior + df_vec[testable]
    p_post  <- 2 * pt(-abs(t_post), df = df_post)

    pvalue_vec[testable] <- p_post

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
