#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData colData assay
#' @importFrom S4Vectors metadata
NULL

# ─── diffMethyl() ─────────────────────────────────────────────────────────────

#' Identify differentially methylated sites between conditions
#'
#' The main function for differential methylation analysis in \pkg{comma}.
#' Analogous to \code{DESeq2::DESeq()}, \code{diffMethyl()} accepts a
#' \code{\link{commaData}} object, fits a statistical model to each methylation
#' site, and returns the same object enriched with per-site test results stored
#' as new columns in \code{rowData}.
#'
#' @details
#' \strong{Statistical model (\code{method = "beta_binomial"}):}
#' A per-site quasibinomial GLM is fitted using \code{\link[stats]{glm}}:
#' \deqn{
#'   \text{Binomial}(n_{\text{mod}},\ n_{\text{total}}) \sim \text{condition}
#' }
#' where \eqn{n_{\text{mod}} = \text{round}(\beta \times \text{coverage})} and
#' \eqn{n_{\text{total}} = \text{coverage}}. The quasibinomial family accounts
#' for overdispersion. P-values are extracted from the Wald t-test on the
#' contrast coefficient. This method requires no additional packages beyond
#' base R.
#'
#' \strong{Alternative model (\code{method = "limma"}):}
#' Beta values are transformed to M-values via
#' \eqn{M = \log_2((n_{\mathrm{mod}} + \alpha) / (n_{\mathrm{unmod}} + \alpha))},
#' then \code{\link[limma]{lmFit}} fits an OLS model per site and
#' \code{\link[limma]{eBayes}} applies empirical Bayes variance shrinkage —
#' borrowing information across all sites to stabilize the per-site variance
#' estimate. This gives substantially more power than \code{"beta_binomial"}
#' when replicates are few (n < 3 per group). Requires \pkg{limma}
#' (\code{BiocManager::install("limma")}). Effect sizes are reported on the
#' original beta scale.
#'
#' \strong{Alternative model (\code{method = "methylkit"}):}
#' Wraps \code{methylKit::calculateDiffMeth()}, which uses logistic regression.
#' Requires \pkg{methylKit} to be installed
#' (\code{BiocManager::install("methylKit")}). Returns results in the same
#' standardised format.
#'
#' \strong{Multiple mod types:} When \code{mod_type = NULL} (default), all
#' modification types present in the object are tested independently and
#' results are combined. Sites of a mod type that is not being tested receive
#' \code{NA} in all \code{dm_*} columns.
#'
#' \strong{Small-sample note:} Differential methylation testing with very few
#' replicates (e.g., n = 1 per group) is mathematically possible but has
#' extremely low statistical power. Treat such results as exploratory only.
#'
#' \strong{Result columns added to \code{rowData}:}
#' \describe{
#'   \item{\code{dm_pvalue}}{Raw p-value from the GLM Wald test.}
#'   \item{\code{dm_padj}}{Adjusted p-value (Benjamini-Hochberg by default).}
#'   \item{\code{dm_delta_beta}}{Effect size: mean methylation in the treatment
#'     group minus mean methylation in the reference (control) group. Positive
#'     values indicate hypermethylation in treatment.}
#'   \item{\code{dm_mean_beta_<condition>}}{One column per condition level
#'     (named after the actual condition values), containing the per-group
#'     mean beta value for each site.}
#' }
#' Analysis parameters and result column names are stored in
#' \code{metadata(object)$diffMethyl_params} and
#' \code{metadata(object)$diffMethyl_result_cols}.
#'
#' @param object A \code{\link{commaData}} object with at least two samples in
#'   distinct conditions.
#' @param formula A one-sided formula specifying the design. The RHS variable
#'   must match a column in \code{sampleInfo(object)} (e.g., \code{~ condition}).
#'   Default is \code{~ condition}.
#' @param method Character string selecting the statistical backend.
#'   \code{"beta_binomial"} (default) uses a quasibinomial GLM via base R.
#'   \code{"limma"} applies empirical Bayes variance shrinkage via
#'   \code{\link[limma]{eBayes}} on M-value-transformed data; recommended when
#'   replicates are few (n < 3 per group). Requires \pkg{limma}.
#'   \code{"methylkit"} wraps \code{methylKit::calculateDiffMeth()}, requiring
#'   \pkg{methylKit} to be installed.
#' @param alpha Positive numeric pseudocount used to compute M-values when
#'   \code{method = "limma"}:
#'   \eqn{M = \log_2((n_{\mathrm{mod}} + \alpha) / (n_{\mathrm{unmod}} + \alpha))}.
#'   Default \code{0.5} (a Beta(0.5, 0.5) prior). Ignored for other methods.
#' @param mod_type Character vector or \code{NULL}. Modification type(s) to
#'   test (e.g., \code{"6mA"}, \code{c("6mA", "5mC")}). If \code{NULL}
#'   (default), all modification types present in \code{object} are tested.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   matching sequence context motif(s) are tested (e.g., \code{"GATC"}).
#'   Uses \code{\link{motifs}} to validate the requested values. If \code{NULL}
#'   (default), all motifs (including \code{NA}) are included.
#' @param min_coverage Integer. Minimum per-sample read depth required to
#'   include a site in testing. Sites where any sample has coverage below
#'   this threshold are treated as \code{NA} in that sample. Default \code{5L}.
#' @param p_adjust_method Character string. Multiple testing correction method,
#'   passed to \code{\link[stats]{p.adjust}}. Default \code{"BH"}
#'   (Benjamini-Hochberg). Other options: \code{"bonferroni"}, \code{"holm"},
#'   \code{"BY"}, \code{"none"}.
#' @param ... Additional arguments (reserved for future use).
#'
#' @return The input \code{commaData} object with additional columns in
#'   \code{rowData}: \code{dm_pvalue}, \code{dm_padj}, \code{dm_delta_beta},
#'   and one \code{dm_mean_beta_<condition>} column per condition level. The
#'   \code{metadata} slot is updated with analysis parameters and result
#'   column names.
#'
#' @seealso \code{\link{results}} to extract the test results as a tidy
#'   \code{data.frame}; \code{\link{filterResults}} to filter by significance
#'   thresholds.
#'
#' @examples
#' data(comma_example_data)
#' # Test for differential 6mA methylation between conditions
#' dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
#'
#' # How many sites are significant?
#' rd <- as.data.frame(SummarizedExperiment::rowData(dm))
#' sum(rd$dm_padj < 0.05 & abs(rd$dm_delta_beta) >= 0.2, na.rm = TRUE)
#'
#' @export
diffMethyl <- function(
    object,
    formula         = ~ condition,
    method          = c("beta_binomial", "methylkit", "limma"),
    mod_type        = NULL,
    motif           = NULL,
    min_coverage    = 5L,
    alpha           = 0.5,
    p_adjust_method = "BH",
    ...
) {
    # ── Input validation ──────────────────────────────────────────────────────
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    if (!inherits(formula, "formula")) {
        stop("'formula' must be a formula object (e.g., ~ condition).")
    }
    method      <- match.arg(method)
    min_coverage <- as.integer(min_coverage)

    if (method == "methylkit" && !requireNamespace("methylKit", quietly = TRUE)) {
        stop(
            "Package 'methylKit' is required for method = \"methylkit\".\n",
            "Install it with: BiocManager::install(\"methylKit\")"
        )
    }

    if (method == "limma" && !requireNamespace("limma", quietly = TRUE)) {
        stop(
            "Package 'limma' is required for method = \"limma\".\n",
            "Install it with: BiocManager::install(\"limma\")"
        )
    }

    if (method == "limma") {
        if (!is.numeric(alpha) || length(alpha) != 1L ||
                !is.finite(alpha) || alpha <= 0) {
            stop("'alpha' must be a single positive finite number.")
        }
    }

    if (!is.null(mod_type)) {
        available <- modTypes(object)
        bad <- setdiff(mod_type, available)
        if (length(bad) > 0L) {
            stop(
                "mod_type value(s) not found in object: ",
                paste(bad, collapse = ", "),
                ". Available: ", paste(available, collapse = ", ")
            )
        }
    }

    if (!is.null(motif)) {
        available_m <- motifs(object)
        bad_m <- setdiff(motif, available_m)
        if (length(bad_m) > 0L) {
            stop(
                "motif value(s) not found in object: ",
                paste(bad_m, collapse = ", "),
                ". Available: ", paste(available_m, collapse = ", ")
            )
        }
        object <- subset(object, motif = motif)
    }

    # Validate formula variable exists in colData
    rhs_vars <- all.vars(formula)
    if (length(rhs_vars) == 0L) {
        stop("'formula' must contain at least one RHS variable (e.g., ~ condition).")
    }
    cd <- as.data.frame(colData(object))
    primary_var <- rhs_vars[[1L]]
    if (!primary_var %in% colnames(cd)) {
        stop(
            "Variable '", primary_var, "' in formula not found in sampleInfo. ",
            "Available columns: ", paste(colnames(cd), collapse = ", ")
        )
    }
    if (length(unique(cd[[primary_var]])) < 2L) {
        stop(
            "Column '", primary_var, "' must have at least 2 distinct values. ",
            "Found: '", unique(cd[[primary_var]]), "'"
        )
    }

    # ── Determine which mod types to test ─────────────────────────────────────
    test_types <- if (!is.null(mod_type)) mod_type else modTypes(object)
    cond_levels <- sort(unique(cd[[primary_var]]))

    # ── Extract full matrices ─────────────────────────────────────────────────
    methyl_full  <- methylation(object)
    cov_full     <- coverage(object)
    rd_full      <- as.data.frame(rowData(object))
    n_sites_all  <- nrow(methyl_full)

    # ── Initialise result columns (all NA) ───────────────────────────────────
    pvalue_all     <- rep(NA_real_, n_sites_all)
    delta_beta_all <- rep(NA_real_, n_sites_all)
    mean_beta_cols <- lapply(cond_levels, function(lv) rep(NA_real_, n_sites_all))
    names(mean_beta_cols) <- paste0("dm_mean_beta_", cond_levels)

    # ── Test each mod type independently ─────────────────────────────────────
    for (mt in test_types) {
        site_idx <- which(rd_full$mod_type == mt)
        if (length(site_idx) == 0L) next

        methyl_sub <- methyl_full[site_idx, , drop = FALSE]
        cov_sub    <- cov_full[site_idx, , drop = FALSE]

        # Apply min_coverage: set beta to NA where coverage < threshold
        low_cov <- !is.na(cov_sub) & cov_sub < min_coverage
        methyl_sub[low_cov] <- NA_real_

        # Dispatch to statistical backend
        res_sub <- tryCatch(
            if (method == "beta_binomial") {
                .betaBinomialTest(methyl_sub, cov_sub, cd, formula)
            } else if (method == "limma") {
                .runLimma(methyl_sub, cov_sub, cd, formula, alpha = alpha)
            } else {
                .runMethylKit(methyl_sub, cov_sub, cd, formula)
            },
            error = function(e) {
                warning(
                    "diffMethyl() failed for mod_type = '", mt, "': ",
                    e$message, ". Skipping this mod type."
                )
                NULL
            }
        )

        if (is.null(res_sub)) next

        # Write back to full-object vectors
        pvalue_all[site_idx]     <- res_sub$pvalue
        delta_beta_all[site_idx] <- res_sub$delta_beta

        for (lv in cond_levels) {
            col_nm <- paste0("mean_beta_", lv)
            dm_col  <- paste0("dm_mean_beta_", lv)
            if (col_nm %in% colnames(res_sub)) {
                mean_beta_cols[[dm_col]][site_idx] <- res_sub[[col_nm]]
            }
        }
    }

    # ── Apply genome-wide multiple testing correction ─────────────────────────
    padj_all <- .applyMultipleTesting(pvalue_all, method = p_adjust_method)

    # ── Update rowData ────────────────────────────────────────────────────────
    rd_new <- rowData(object)
    rd_new$dm_pvalue     <- pvalue_all
    rd_new$dm_padj       <- padj_all
    rd_new$dm_delta_beta <- delta_beta_all
    for (col_nm in names(mean_beta_cols)) {
        rd_new[[col_nm]] <- mean_beta_cols[[col_nm]]
    }

    result_cols <- c("dm_pvalue", "dm_padj", "dm_delta_beta",
                     names(mean_beta_cols))

    # ── Build returned object ─────────────────────────────────────────────────
    # Reconstruct commaData with updated rowData
    se_new <- SummarizedExperiment::SummarizedExperiment(
        assays  = SummarizedExperiment::assays(object),
        rowData = rd_new,
        colData = colData(object)
    )
    out <- new("commaData", se_new,
               genomeInfo = object@genomeInfo,
               annotation = object@annotation,
               motifSites = object@motifSites)

    # Copy existing metadata, then add diffMethyl entries
    S4Vectors::metadata(out) <- S4Vectors::metadata(object)
    S4Vectors::metadata(out)$diffMethyl_result_cols <- result_cols
    S4Vectors::metadata(out)$diffMethyl_params <- list(
        formula         = deparse(formula),
        method          = method,
        mod_type        = mod_type,
        p_adjust_method = p_adjust_method,
        min_coverage    = min_coverage,
        alpha           = alpha,
        timestamp       = Sys.time()
    )

    out
}
