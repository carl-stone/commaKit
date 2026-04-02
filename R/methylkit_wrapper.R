NULL

# ─── methylKit wrapper ────────────────────────────────────────────────────────

#' Run differential methylation via methylKit
#'
#' An internal wrapper that uses \pkg{methylKit} to test for differential
#' methylation. Called by \code{\link{diffMethyl}} when
#' \code{method = "methylkit"}.
#'
#' \pkg{methylKit} must be installed (it is listed in \code{Suggests}). If it
#' is not available, this function stops with an informative message.
#'
#' @details
#' The wrapper converts the methylation and coverage matrices from a
#' \code{\link{commaData}} object into the format expected by
#' \code{methylKit::methylRawList}, runs \code{methylKit::unite()} and
#' \code{methylKit::calculateDiffMeth()}, and returns results in the same
#' standardised format as \code{.betaBinomialTest()}.
#'
#' Only the first RHS variable of \code{formula} is used as the grouping
#' variable. Complex formulas with interactions or batch terms are not currently
#' supported by this wrapper.
#'
#' @param methyl_mat Numeric matrix (sites × samples) of beta values.
#' @param coverage_mat Integer matrix (sites × samples) of read depths.
#' @param coldata \code{data.frame} with at least one column matching the
#'   RHS variable in \code{formula}.
#' @param formula One-sided formula (e.g., \code{~ condition}).
#'
#' @return A \code{data.frame} with the same columns as \code{.betaBinomialTest()}:
#'   \code{pvalue}, \code{delta_beta}, and one \code{mean_beta_<level>} column
#'   per condition level. Row names are site keys.
#'
#' @keywords internal
.runMethylKit <- function(methyl_mat, coverage_mat, coldata, formula,
                          ref_level = NULL) {
    # ── Dependency check ──────────────────────────────────────────────────────
    if (!requireNamespace("methylKit", quietly = TRUE)) {
        stop(
            "Package 'methylKit' is required for method = \"methylkit\".\n",
            "Install it with: BiocManager::install(\"methylKit\")\n",
            "Alternatively, use method = \"beta_binomial\" (no extra packages needed)."
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

    # methylKit requires integer treatment codes: 0 = reference, 1 = treatment
    # Use provided ref_level, or fall back to alphabetically first
    if (is.null(ref_level)) {
        ref_level <- all_levels[[1L]]
    }
    cond_levels <- c(ref_level, setdiff(all_levels, ref_level))
    treat_level <- cond_levels[[2L]]
    treatment   <- as.integer(cond != ref_level)  # 0/1

    # ── Parse site keys into components ───────────────────────────────────────
    site_keys <- rownames(methyl_mat)
    # site_key format: "chrom:position:strand:mod_type"
    key_parts <- strsplit(site_keys, ":", fixed = TRUE)
    chroms    <- vapply(key_parts, `[[`, character(1), 1L)
    positions <- as.integer(vapply(key_parts, `[[`, character(1), 2L))
    strands   <- vapply(key_parts, `[[`, character(1), 3L)
    n_sites   <- length(site_keys)

    # ── Build methylKit objects per sample ─────────────────────────────────────
    if (requireNamespace("methylKit", quietly = TRUE)) {
      # Construct methylRaw objects directly — bypasses methRead() which only
      # accepts file paths, not in-memory data frames
      sample_list <- lapply(seq_len(ncol(methyl_mat)), function(j) {
        beta_j <- methyl_mat[, j]
        cov_j  <- as.integer(coverage_mat[, j])

        # Replace NA with 0 coverage for methylKit (it handles 0-coverage sites
        # via unite() with min.per.group)
        cov_j[is.na(cov_j)]  <- 0L
        cov_j[is.na(beta_j)] <- 0L

        n_meth <- as.integer(round(beta_j * cov_j))
        n_meth[is.na(n_meth)] <- 0L
        n_meth <- pmax(0L, pmin(n_meth, cov_j))

        df <- data.frame(
          chr      = chroms,
          start    = positions,
          end      = positions,
          strand   = strands,
          coverage = cov_j,
          numCs    = n_meth,
          numTs    = cov_j - n_meth,
          stringsAsFactors = FALSE
        )

        # methylRaw extends data.frame; pass df as first positional arg
        methods::new("methylRaw", df,
                     sample.id  = colnames(methyl_mat)[[j]],
                     assembly   = "custom",
                     context    = "none",
                     resolution = "base")
      })

      mk_list <- methods::new("methylRawList",
                              .Data     = sample_list,
                              treatment = as.integer(treatment))
    }

    # ── Unite and test ────────────────────────────────────────────────────────
    if (requireNamespace("methylKit", quietly = TRUE)) {
      mk_united <- tryCatch(
        methylKit::unite(mk_list, destrand = FALSE, min.per.group = 1L),
        error = function(e) {
          stop("methylKit::unite() failed: ", e$message)
        }
      )

      mk_diff <- tryCatch(
        suppressWarnings(
          methylKit::calculateDiffMeth(mk_united, weighted.mean = FALSE)
        ),
        error = function(e) {
          stop("methylKit::calculateDiffMeth() failed: ", e$message)
        }
      )
    }

    # ── Extract results and standardise format ────────────────────────────────
    diff_df  <- as.data.frame(mk_diff)
    # methylKit result columns: chr, start, end, strand, pvalue, qvalue, meth.diff
    # meth.diff = methylation % difference (treatment - control), in percent

    # Build a lookup from position (as character) back to site key
    mk_key <- paste0(diff_df$chr, ":", diff_df$start, ":", diff_df$strand)

    # Pre-compute group means from our matrices (same logic as beta_binomial)
    group_idx  <- lapply(cond_levels, function(lv) which(cond == lv))
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
    pvalue_vec     <- rep(NA_real_, n_sites)

    # Match methylKit results back to original site order by key
    our_key_for_match <- paste0(chroms, ":", positions, ":", strands)
    mk_match <- match(our_key_for_match, mk_key)
    found    <- !is.na(mk_match)
    if (any(found)) {
        pvalue_vec[found] <- diff_df$pvalue[mk_match[found]]
    }

    # ── Assemble result ───────────────────────────────────────────────────────
    result <- data.frame(
        pvalue     = pvalue_vec,
        delta_beta = delta_beta_vec,
        row.names  = site_keys,
        stringsAsFactors = FALSE
    )
    for (lv in cond_levels) {
        result[[paste0("mean_beta_", lv)]] <- group_means[, lv]
    }

    result
}
