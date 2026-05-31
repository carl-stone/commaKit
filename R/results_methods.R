#' @importFrom methods setGeneric setMethod is
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
NULL

# ─── results() ────────────────────────────────────────────────────────────────

#' Extract differential methylation results as a tidy data frame
#'
#' Retrieves the per-site differential methylation statistics added by
#' \code{\link{diffMethyl}} and returns them as a tidy \code{data.frame}
#' suitable for downstream analysis and plotting.
#'
#' @param object A \code{\link{commaData}} object on which
#'   \code{\link{diffMethyl}} has been run.
#' @param mod_type Character vector or \code{NULL}. If provided, only sites
#'   of the specified modification type(s) are returned. If \code{NULL} (default),
#'   results for all modification types are returned.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   matching sequence context motif(s) are returned. If \code{NULL} (default),
#'   all motifs are returned.
#' @param mod_context Character vector or \code{NULL}. If provided, only sites
#'   with a matching modification context (e.g., \code{"6mA_GATC"}) are
#'   returned. Use \code{\link{modContexts}(object)} to see available values.
#'   Applied in addition to any \code{mod_type} or \code{motif} filters.
#' @param ... Ignored (for S4 generic compatibility).
#'
#' @return A \code{data.frame} with one row per methylation site, containing:
#'   \describe{
#'     \item{\code{chrom}}{Chromosome name.}
#'     \item{\code{position}}{1-based genomic position.}
#'     \item{\code{strand}}{Strand (\code{"+"} or \code{"-"}).}
#'     \item{\code{mod_type}}{Modification type (e.g., \code{"6mA"}).}
#'     \item{\code{dm_pvalue}}{Raw p-value from the statistical test.}
#'     \item{\code{dm_padj}}{Adjusted p-value (Benjamini-Hochberg by default).}
#'     \item{\code{dm_delta_beta}}{Effect size: mean methylation in the
#'       treatment group minus mean methylation in the reference group.}
#'     \item{\code{dm_mean_beta_<condition>}}{One column per condition level
#'       with per-group mean beta values.}
#'   }
#'   Any other annotation columns present in \code{rowData(object)} (e.g.,
#'   from \code{\link{annotateSites}}) are also included.
#'
#' @seealso \code{\link{diffMethyl}}, \code{\link{filterResults}}
#'
#' @examples
#' data(comma_example_data)
#' dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
#' res <- results(dm)
#' head(res[order(res$dm_padj), ])
#'
#' @export
setGeneric("results", function(object, ...) standardGeneric("results"))

#' @rdname results
setMethod("results", "commaData", function(object, mod_type = NULL, motif = NULL,
                                           mod_context = NULL, ...) {
    # ── Check diffMethyl has been run ─────────────────────────────────────────
    md <- metadata(object)
    if (is.null(md$diffMethyl_result_cols)) {
        stop(
            "No differential methylation results found in this commaData object.\n",
            "run diffMethyl() first:\n",
            "  dm <- diffMethyl(object, formula = ~ condition)"
        )
    }

    object <- .applySiteFilters(
        object,
        mod_type = mod_type,
        motif = motif,
        mod_context = mod_context,
        caller = "results()"
    )

    as.data.frame(siteInfo(object))
})

# ─── filterResults() ──────────────────────────────────────────────────────────

#' Filter differential methylation results by significance thresholds
#'
#' A convenience wrapper around \code{\link{results}} that filters sites by
#' adjusted p-value and absolute effect size.
#'
#' @param object A \code{\link{commaData}} object on which
#'   \code{\link{diffMethyl}} has been run.
#' @param padj Numeric. Maximum adjusted p-value threshold (inclusive).
#'   Default \code{0.05}.
#' @param delta_beta Numeric. Minimum absolute effect size threshold
#'   (\eqn{|\Delta\beta|}) (inclusive). Default \code{0.1}. Set to \code{0} to
#'   disable filtering on effect size.
#' @param mod_type Character vector or \code{NULL}. Passed to
#'   \code{\link{results}} for optional modification type filtering.
#' @param motif Character vector or \code{NULL}. Passed to
#'   \code{\link{results}} for optional sequence context motif filtering.
#' @param mod_context Character vector or \code{NULL}. Passed to
#'   \code{\link{results}} for optional modification context filtering
#'   (e.g., \code{"6mA_GATC"}).
#' @param ... Ignored.
#'
#' @return A \code{data.frame} (same format as \code{\link{results}}) containing
#'   only sites where \code{dm_padj <= padj} \strong{and}
#'   \code{abs(dm_delta_beta) >= delta_beta}. Sites with \code{NA} values in
#'   either column are excluded.
#'
#' @seealso \code{\link{diffMethyl}}, \code{\link{results}}
#'
#' @examples
#' data(comma_example_data)
#' dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
#' sig <- filterResults(dm, padj = 0.05, delta_beta = 0.2)
#' nrow(sig)
#'
#' @export
setGeneric("filterResults",
           function(object, ...) standardGeneric("filterResults"))

#' @rdname filterResults
setMethod("filterResults", "commaData",
          function(object, padj = 0.05, delta_beta = 0.1,
                   mod_type = NULL, motif = NULL, mod_context = NULL, ...) {
    res <- results(object, mod_type = mod_type, motif = motif,
                   mod_context = mod_context)

    if (!"dm_padj" %in% colnames(res)) {
        stop(
            "Column 'dm_padj' not found in results. ",
            "Ensure diffMethyl() has completed successfully."
        )
    }
    if (!"dm_delta_beta" %in% colnames(res)) {
        stop(
            "Column 'dm_delta_beta' not found in results. ",
            "Ensure diffMethyl() has completed successfully."
        )
    }

    if (!is.numeric(padj) || length(padj) != 1L ||
            is.na(padj) || !is.finite(padj) || padj < 0) {
        stop("'padj' must be a single non-NA, non-negative finite number.")
    }
    if (!is.numeric(delta_beta) || length(delta_beta) != 1L ||
            is.na(delta_beta) || !is.finite(delta_beta) || delta_beta < 0) {
        stop("'delta_beta' must be a single non-NA, non-negative finite number.")
    }

    keep <- !is.na(res$dm_padj) &
            !is.na(res$dm_delta_beta) &
            res$dm_padj <= padj &
            abs(res$dm_delta_beta) >= delta_beta

    res[keep, , drop = FALSE]
})
