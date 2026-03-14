#' Volcano plot for differential methylation results
#'
#' Produces a volcano plot from differential methylation results computed by
#' \code{\link{diffMethyl}()} and extracted with \code{\link{results}()}.
#' Sites are colored by significance category based on user-supplied adjusted
#' p-value and delta-beta thresholds.
#'
#' @param results A \code{data.frame} returned by \code{\link{results}()},
#'   containing at minimum the columns \code{dm_delta_beta} (numeric, effect
#'   size as beta difference treatment minus control) and \code{dm_padj}
#'   (numeric, BH-adjusted p-value in [0, 1]).
#' @param delta_beta_threshold Numeric scalar in (0, 1). Sites with
#'   \code{|dm_delta_beta| >= delta_beta_threshold} AND
#'   \code{dm_padj <= padj_threshold} are considered significant. Default
#'   \code{0.2}.
#' @param padj_threshold Numeric scalar in (0, 1). Adjusted p-value cutoff for
#'   significance. Default \code{0.05}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The x-axis shows
#'   \code{dm_delta_beta} (effect size), the y-axis shows
#'   \code{-log10(dm_padj)} (significance). Sites are colored as
#'   \code{"Hypermethylated"} (positive delta-beta, significant),
#'   \code{"Hypomethylated"} (negative delta-beta, significant), or
#'   \code{"Not significant"}. Dashed lines mark the threshold boundaries.
#'   Sites with \code{NA} adjusted p-value are excluded from the plot.
#'
#' @examples
#' data(comma_example_data)
#' cd_dm <- diffMethyl(comma_example_data, ~ condition)
#' res <- results(cd_dm)
#' plot_volcano(res)
#'
#' # Custom thresholds
#' plot_volcano(res, delta_beta_threshold = 0.3, padj_threshold = 0.01)
#'
#' @seealso \code{\link{diffMethyl}}, \code{\link{results}},
#'   \code{\link{filterResults}}
#'
#' @export
plot_volcano <- function(results,
                         delta_beta_threshold = 0.2,
                         padj_threshold = 0.05) {

    ## --- Input validation ---------------------------------------------------
    if (!is.data.frame(results)) {
        stop("'results' must be a data.frame returned by results().")
    }
    required_cols <- c("dm_delta_beta", "dm_padj")
    missing_cols <- setdiff(required_cols, colnames(results))
    if (length(missing_cols) > 0) {
        stop("'results' is missing required column(s): ",
             paste(missing_cols, collapse = ", "), ". ",
             "Ensure 'results' was produced by results() after diffMethyl().")
    }
    if (!is.numeric(delta_beta_threshold) || length(delta_beta_threshold) != 1L ||
        is.na(delta_beta_threshold) || delta_beta_threshold <= 0 ||
        delta_beta_threshold >= 1) {
        stop("'delta_beta_threshold' must be a single numeric value in (0, 1).")
    }
    if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
        is.na(padj_threshold) || padj_threshold <= 0 || padj_threshold >= 1) {
        stop("'padj_threshold' must be a single numeric value in (0, 1).")
    }

    ## --- Build plot data ----------------------------------------------------
    df <- results[, c("dm_delta_beta", "dm_padj"), drop = FALSE]

    ## Exclude rows with NA padj
    df <- df[!is.na(df$dm_padj), , drop = FALSE]

    if (nrow(df) == 0L) {
        stop("No rows with non-NA 'dm_padj' values found in 'results'.")
    }

    ## Compute y-axis: -log10(padj); clamp padj=0 to min positive double to
    ## avoid -Inf (can happen with very extreme p-values)
    padj_clamped <- pmax(df$dm_padj, .Machine$double.xmin)
    df$neg_log10_padj <- -log10(padj_clamped)

    ## Assign significance category
    is_sig_p    <- df$dm_padj <= padj_threshold
    is_sig_db_pos <- df$dm_delta_beta >= delta_beta_threshold
    is_sig_db_neg <- df$dm_delta_beta <= -delta_beta_threshold

    df$significance <- "Not significant"
    df$significance[is_sig_p & is_sig_db_pos] <- "Hypermethylated"
    df$significance[is_sig_p & is_sig_db_neg] <- "Hypomethylated"
    df$significance <- factor(df$significance,
                              levels = c("Hypermethylated",
                                         "Hypomethylated",
                                         "Not significant"))

    ## --- Build ggplot -------------------------------------------------------
    sig_colors <- c(
        "Hypermethylated" = "#d73027",
        "Hypomethylated"  = "#4575b4",
        "Not significant" = "grey60"
    )

    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = .data[["dm_delta_beta"]],
            y = .data[["neg_log10_padj"]],
            color = .data[["significance"]]
        )
    ) +
        ggplot2::geom_point(alpha = 0.6, size = 1.2) +
        ggplot2::geom_vline(
            xintercept = c(-delta_beta_threshold, delta_beta_threshold),
            linetype = "dashed", color = "grey40", linewidth = 0.5
        ) +
        ggplot2::geom_hline(
            yintercept = -log10(padj_threshold),
            linetype = "dashed", color = "grey40", linewidth = 0.5
        ) +
        ggplot2::scale_color_manual(
            values = sig_colors,
            drop = FALSE
        ) +
        ggplot2::labs(
            x     = "\u0394\u03b2 (treatment \u2212 control)",
            y     = "-log\u2081\u2080(adjusted p-value)",
            title = "Differential Methylation Volcano Plot",
            color = NULL
        ) +
        ggplot2::theme_bw()

    p
}
