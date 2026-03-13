NULL

# ─── Multiple testing correction ──────────────────────────────────────────────

#' Apply multiple testing correction to a vector of p-values
#'
#' A thin wrapper around \code{\link[stats]{p.adjust}} that operates on a named
#' numeric vector of p-values and returns adjusted values in the same order.
#' Used internally by \code{\link{diffMethyl}}.
#'
#' @param pvalues Named numeric vector of raw p-values. \code{NA} values are
#'   passed through (i.e., sites that could not be tested remain \code{NA}
#'   after adjustment).
#' @param method Character string. Correction method passed to
#'   \code{\link[stats]{p.adjust}}. Default \code{"BH"} (Benjamini-Hochberg
#'   FDR control). Other options include \code{"bonferroni"}, \code{"holm"},
#'   \code{"BY"}, \code{"fdr"}, and \code{"none"}.
#'
#' @return Named numeric vector of adjusted p-values, same length and names as
#'   \code{pvalues}.
#'
#' @keywords internal
.applyMultipleTesting <- function(pvalues, method = "BH") {
    stats::p.adjust(pvalues, method = method)
}
