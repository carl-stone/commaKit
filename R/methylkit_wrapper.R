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
#'   RHS variable in \code{formula}.
#' @param formula One-sided formula (e.g., \code{~ condition}).
#'
#' @return A \code{data.frame} with the same columns as \code{.betaBinomialTest()}:
#'   \code{pvalue}, \code{delta_beta}, and one \code{mean_beta_<level>} column
#'   per condition level. Row names are site keys.
#'
#' @keywords internal
.runMethylKit <- function(methyl_mat, coverage_mat, site_df, coldata, formula,
                          ref_level = NULL, design_info = NULL,
                          mod_counts_mat = NULL, canonical_counts_mat = NULL,
                          other_mod_counts_mat = NULL) {
    # ── Dependency check ──────────────────────────────────────────────────────
    if (!requireNamespace("methylKit", quietly = TRUE)) {
        stop(
            "Package 'methylKit' is required for method = \"methylkit\".\n",
            "Install it with: BiocManager::install(\"methylKit\")\n",
            "Install it with: BiocManager::install(\"methylKit\"), or use method = \"quasi_f\" or \"limma\"."
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
    treatment   <- as.integer(cond == treat_level)  # 0 = reference, 1 = treatment
    count_mats <- .resolveCountMatrices(
        methyl_mat,
        coverage_mat,
        mod_counts_mat = mod_counts_mat,
        canonical_counts_mat = canonical_counts_mat,
        other_mod_counts_mat = other_mod_counts_mat
    )

    message(
        "methylKit: comparing '", treat_level, "' (treatment) vs '",
        ref_level, "' (reference/control)"
    )

    # ── Get site positions from site_df ───────────────────────────────────────
    chroms    <- site_df$chrom
    positions <- site_df$position
    strands   <- site_df$strand
    n_sites   <- nrow(site_df)

    # ── Filter zero-variance sites ───────────────────────────────────────────
    # methylKit's calculateDiffMeth() crashes when any site has identical
    # counts across all samples (all methylated or all unmethylated). The
    # internal mgcv::uniquecombs() returns a vector instead of a matrix, and
    # methylKit's reshape logic hardcodes ncol = 4 (2-sample assumption),
    # producing garbage that cascades into a split() error. Filter these sites
    # out before building methylRaw objects. These sites are untestable, so
    # they retain p = NA and are excluded from multiple-testing correction.
    all_meth  <- apply(methyl_mat, 1L, function(x) all(x == 1, na.rm = TRUE))
    all_unmeth <- apply(methyl_mat, 1L, function(x) all(x == 0, na.rm = TRUE))
    skip_idx  <- which(all_meth | all_unmeth)
    n_skipped <- length(skip_idx)

    if (n_skipped > 0L) {
        message(
            "methylKit: ", n_skipped, " site(s) with all samples at ",
            "beta = 1 or beta = 0 excluded (zero within-group variance; ",
            "methylKit GLM cannot estimate). Assigned p = NA."
        )
    }

    keep_idx <- if (n_skipped > 0L) setdiff(seq_len(n_sites), skip_idx) else seq_len(n_sites)

    # If all sites are zero-variance, return early with p = NA
    if (length(keep_idx) == 0L) {
        group_stats <- .computeDiffMethylGroupStats(methyl_mat, design_info)
        result <- data.frame(
            pvalue     = rep(NA_real_, n_sites),
            delta_beta = group_stats$delta_beta,
            stringsAsFactors = FALSE
        )
        for (lv in cond_levels) {
            result[[paste0("mean_beta_", lv)]] <- group_stats$group_means[, lv]
        }
        return(result)
    }

    # ── Build methylKit objects per sample ─────────────────────────────────────
    if (requireNamespace("methylKit", quietly = TRUE)) {
      # Construct methylRaw objects directly — bypasses methRead() which only
      # accepts file paths, not in-memory data frames
      sample_list <- lapply(seq_len(ncol(methyl_mat)), function(j) {
        beta_j <- methyl_mat[keep_idx, j]
        n_meth <- as.integer(count_mats$modified[keep_idx, j])
        n_unmeth <- as.integer(count_mats$unmodified[keep_idx, j])
        cov_j <- n_meth + n_unmeth

        # Replace NA with 0 coverage for methylKit (it handles 0-coverage sites
        # via unite() with min.per.group)
        cov_j[is.na(cov_j)]  <- 0L
        n_meth[is.na(n_meth)] <- 0L
        n_unmeth[is.na(n_unmeth)] <- 0L
        cov_j[is.na(beta_j)] <- 0L
        n_meth <- pmax(0L, pmin(n_meth, cov_j))
        n_unmeth <- pmax(0L, cov_j - n_meth)

        df <- data.frame(
          chr      = chroms[keep_idx],
          start    = positions[keep_idx],
          end      = positions[keep_idx],
          strand   = strands[keep_idx],
          coverage = cov_j,
          numCs    = n_meth,
          numTs    = n_unmeth,
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

      # methylKit's unite() does not filter sites where ALL samples have
      # zero coverage (e.g. when every sample had coverage below min_coverage
      # and was set to 0).  glm.fit inside logReg crashes with
      # "object of type 'closure' is not subsettable" when all weights are 0.
      # Remove such sites from the united object before calling
      # calculateDiffMeth; they retain p = NA via the untestable-site path above.
      united_df  <- methylKit::getData(mk_united)
      cov_cols   <- seq(5L, ncol(united_df), by = 3L)
      site_total_cov <- rowSums(as.matrix(united_df[, cov_cols, drop = FALSE]))
      keep_united <- which(site_total_cov > 0L)
      n_dropped   <- nrow(united_df) - length(keep_united)
      if (n_dropped > 0L) {
        message(
          "methylKit: ", n_dropped, " site(s) with zero total coverage across ",
          "all samples removed before testing (p = NA assigned)."
        )
        mk_united <- mk_united[keep_united, ]
      }

      mk_warn_counts <- list()
      mk_diff <- tryCatch(
        withCallingHandlers(
          suppressMessages(   # suppress only the "group: 0/1" message we've replaced
            methylKit::calculateDiffMeth(mk_united, weighted.mean = FALSE)
          ),
          warning = function(w) {
            key <- trimws(conditionMessage(w))
            n_prev <- mk_warn_counts[[key]]
            mk_warn_counts[[key]] <<- (if (is.null(n_prev)) 0L else n_prev) + 1L
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          stop("methylKit::calculateDiffMeth() failed: ", e$message)
        }
      )
      .emitMethylKitWarnings(mk_warn_counts)
    }

    # ── Extract results and standardise format ────────────────────────────────
    diff_df  <- as.data.frame(mk_diff)
    # methylKit result columns: chr, start, end, strand, pvalue, qvalue, meth.diff
    # meth.diff = methylation % difference (treatment - control), in percent

    # Build a lookup from position (as character) back to site key
    mk_key <- paste0(diff_df$chr, ":", diff_df$start, ":", diff_df$strand)

    # Pre-compute group means from our matrices (mirroring quasi_f/limma wrappers)
    group_stats    <- .computeDiffMethylGroupStats(methyl_mat, design_info)
    group_means    <- group_stats$group_means
    delta_beta_vec <- group_stats$delta_beta
    # Initialise p = NA for all sites; untestable sites stay NA and are excluded
    # from multiple-testing correction, matching the limma and quasi_f backends.
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
        stringsAsFactors = FALSE
    )
    for (lv in cond_levels) {
        result[[paste0("mean_beta_", lv)]] <- group_means[, lv]
    }

    result
}

# ─── Warning translation helper ───────────────────────────────────────────────

#' Translate and emit collected methylKit warnings as informative summaries
#'
#' Called after \code{methylKit::calculateDiffMeth()} to emit one summarized
#' warning per unique warning type rather than per-site spam. Unknown warning
#' types are re-emitted verbatim so nothing is silently lost.
#'
#' @param warn_counts Named list mapping warning message text to integer counts.
#' @keywords internal
.emitMethylKitWarnings <- function(warn_counts) {
    if (length(warn_counts) == 0L) return(invisible(NULL))
    for (msg in names(warn_counts)) {
        n <- warn_counts[[msg]]
        if (grepl("fitted probabilities numerically 0 or 1", msg)) {
            warning(
                "methylKit: ", n, " site(s) had fitted probabilities at 0 or 1. ",
                "This typically occurs when methylation fractions are at boundary ",
                "values (fully methylated or unmethylated), common in 6mA datasets ",
                "with few replicates. GLM estimates may be unreliable for these sites. ",
                "Consider method = \"limma\" if many sites are affected.",
                call. = FALSE
            )
        } else if (grepl("algorithm did not converge", msg)) {
            warning(
                "methylKit: GLM failed to converge for ", n, " site(s). ",
                "Results at these sites may be unreliable.",
                call. = FALSE
            )
        } else {
            # Unknown warning type — re-emit verbatim so nothing is silently lost
            warning("methylKit (", n, "x): ", msg, call. = FALSE)
        }
    }
}
