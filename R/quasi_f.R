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
.runQuasiF <- function(methyl_mat, coverage_mat, site_df, coldata, formula,
                       ref_level = NULL, design_info = NULL,
                       mod_counts_mat = NULL, canonical_counts_mat = NULL,
                       other_mod_counts_mat = NULL) {
    # ── Dependency check ──────────────────────────────────────────────────────
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop(
            "Package 'limma' is required for method = \"quasi_f\".\n",
            "Install it with: BiocManager::install(\"limma\")\n",
            "Alternatively, use method = \"methylkit\" if methylKit is available."
        )
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

        n_mod   <- count_mats$modified[i, ok]
        n_unmod <- count_mats$unmodified[i, ok]

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
    had_data <- rowSums(!is.na(methyl_mat) & !is.na(coverage_mat) & coverage_mat > 0L) >= 2L
    failed_with_data <- had_data & (is.na(phi_vec) | is.na(t_unscaled))
    if (sum(had_data) > 0L && sum(failed_with_data) / sum(had_data) > 0.5) {
        warning(
            "quasi_f: GLM fitting failed for ", sum(failed_with_data),
            " of ", sum(had_data), " sites with sufficient observed data. ",
            "These sites retain p = NA.",
            call. = FALSE
        )
    }

    testable   <- !is.na(phi_vec) & !is.na(t_unscaled)
    pvalue_vec <- rep(NA_real_, n_sites)

    if (sum(testable) < 2L) {
        # Not enough sites to estimate the EB prior; return all NA
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
        stringsAsFactors = FALSE
    )
    for (lv in cond_levels) {
        result[[paste0("mean_beta_", lv)]] <- group_means[, lv]
    }

    result
}
