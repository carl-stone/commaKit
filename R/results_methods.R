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
#' @param mod_type Character string or \code{NULL}. If provided, only sites
#'   of the specified modification type are returned. If \code{NULL} (default),
#'   results for all modification types are returned.
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
setMethod("results", "commaData", function(object, mod_type = NULL, ...) {
    # ── Check diffMethyl has been run ─────────────────────────────────────────
    md <- metadata(object)
    if (is.null(md$diffMethyl_result_cols)) {
        stop(
            "No differential methylation results found in this commaData object.\n",
            "run diffMethyl() first:\n",
            "  dm <- diffMethyl(object, formula = ~ condition)"
        )
    }

    # ── Extract rowData as data.frame ─────────────────────────────────────────
    rd <- as.data.frame(rowData(object))

    # ── Optional mod_type filter ──────────────────────────────────────────────
    if (!is.null(mod_type)) {
        available <- unique(rd$mod_type)
        bad <- setdiff(mod_type, available)
        if (length(bad) > 0L) {
            stop(
                "'mod_type' value(s) not found: ",
                paste(bad, collapse = ", "),
                ". Available: ", paste(sort(available), collapse = ", ")
            )
        }
        rd <- rd[rd$mod_type %in% mod_type, , drop = FALSE]
    }

    rd
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
#' @param mod_type Character string or \code{NULL}. Passed to
#'   \code{\link{results}} for optional modification type filtering.
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
          function(object, padj = 0.05, delta_beta = 0.1, mod_type = NULL, ...) {
    res <- results(object, mod_type = mod_type)

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

    keep <- !is.na(res$dm_padj) &
            !is.na(res$dm_delta_beta) &
            res$dm_padj <= padj &
            abs(res$dm_delta_beta) >= delta_beta

    res[keep, , drop = FALSE]
})
