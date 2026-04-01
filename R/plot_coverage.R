#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline facet_wrap
#'   scale_x_log10 labs theme_bw waiver
NULL

#' Plot coverage depth distribution
#'
#' Produces a histogram of sequencing depth (coverage) across sites for each
#' sample in a \code{\link{commaData}} object. Useful for QC to assess whether
#' coverage is sufficient and consistent across samples.
#'
#' @param object A \code{\link{commaData}} object.
#' @param mod_type Character string specifying a single modification type
#'   (e.g., \code{"6mA"}, \code{"5mC"}). If \code{NULL} (default), all sites
#'   from all modification types are included.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   matching sequence context motif(s) are included (e.g., \code{"GATC"}).
#'   If \code{NULL} (default), all motifs are included.
#' @param per_sample Logical. If \code{TRUE} (default), the plot is faceted by
#'   sample, producing one histogram panel per sample. If \code{FALSE}, all
#'   samples are overlaid on a single plot with per-sample colors.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The x-axis shows coverage
#'   depth on a log10 scale; the y-axis shows the number of sites at each
#'   depth. A vertical dashed line marks the median coverage per sample (when
#'   \code{per_sample = TRUE}) or across all samples (when
#'   \code{per_sample = FALSE}). Sites with \code{NA} coverage are silently
#'   excluded.
#'
#' @examples
#' data(comma_example_data)
#' plot_coverage(comma_example_data)
#'
#' # Overlay all samples on one plot
#' plot_coverage(comma_example_data, per_sample = FALSE)
#'
#' # One modification type only
#' plot_coverage(comma_example_data, mod_type = "6mA")
#'
#' @seealso \code{\link{coverageDepth}}, \code{\link{varianceByDepth}},
#'   \code{\link{plot_methylation_distribution}}
#'
#' @export
plot_coverage <- function(object,
                          mod_type    = NULL,
                          motif       = NULL,
                          mod_context = NULL,
                          per_sample  = TRUE) {

    ## --- Input validation ---------------------------------------------------
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    if (!is.logical(per_sample) || length(per_sample) != 1L || is.na(per_sample)) {
        stop("'per_sample' must be TRUE or FALSE.")
    }

    ## --- Optional mod_type filter -------------------------------------------
    if (!is.null(mod_type)) {
        available <- modTypes(object)
        if (!mod_type %in% available) {
            stop("'mod_type' = '", mod_type, "' not found in object. ",
                 "Available types: ", paste(available, collapse = ", "), ".")
        }
        object <- subset(object, mod_type = mod_type)
    }
    if (!is.null(motif)) {
        object <- subset(object, motif = motif)
    }
    if (!is.null(mod_context)) {
        object <- subset(object, mod_context = mod_context)
    }

    ## --- Extract data -------------------------------------------------------
    cov_mat    <- coverage(object)
    si         <- sampleInfo(object)
    sample_nms <- colnames(cov_mat)
    n_sites    <- nrow(cov_mat)
    n_samples  <- length(sample_nms)

    ## Reshape to long data.frame
    df <- data.frame(
        depth       = as.vector(cov_mat),
        sample_name = rep(sample_nms, each = n_sites),
        stringsAsFactors = FALSE
    )

    ## Drop NA coverage values
    df <- df[!is.na(df$depth), , drop = FALSE]

    if (nrow(df) == 0L) {
        stop("No non-NA coverage values found after filtering.")
    }

    ## Join condition from sampleInfo
    si_sub <- si[, c("sample_name", "condition"), drop = FALSE]
    df <- merge(df, si_sub, by = "sample_name", all.x = TRUE)

    ## Compute median coverage per sample for vlines
    med_per_samp <- tapply(df$depth, df$sample_name, stats::median)
    med_df <- data.frame(
        sample_name = names(med_per_samp),
        median_depth = as.numeric(med_per_samp),
        stringsAsFactors = FALSE
    )

    ## --- Build ggplot -------------------------------------------------------
    if (per_sample) {
        p <- ggplot2::ggplot(
            df,
            ggplot2::aes(x = .data[["depth"]], fill = .data[["sample_name"]])
        ) +
            ggplot2::geom_histogram(alpha = 0.8, bins = 30, color = "white") +
            ggplot2::geom_vline(
                data = med_df,
                ggplot2::aes(xintercept = .data[["median_depth"]]),
                linetype = "dashed", color = "grey30", linewidth = 0.6
            ) +
            ggplot2::facet_wrap("sample_name") +
            ggplot2::labs(fill = "Sample")
    } else {
        overall_median <- stats::median(df$depth)
        p <- ggplot2::ggplot(
            df,
            ggplot2::aes(
                x     = .data[["depth"]],
                fill  = .data[["sample_name"]],
                color = .data[["sample_name"]]
            )
        ) +
            ggplot2::geom_histogram(
                alpha = 0.4, bins = 30, position = "identity"
            ) +
            ggplot2::geom_vline(
                xintercept = overall_median,
                linetype = "dashed", color = "grey30", linewidth = 0.6
            ) +
            ggplot2::labs(fill = "Sample", color = "Sample")
    }

    p <- p +
        ggplot2::scale_x_log10(
            labels = scales_comma_label()
        ) +
        ggplot2::labs(
            x     = expression("Coverage depth (reads, " * log[10] * " scale)"),
            y     = "Number of sites",
            title = "Coverage Depth Distribution"
        ) +
        ggplot2::theme_bw()

    p
}

## Internal helper: format axis labels with commas if scales is available,
## otherwise fall back to default.
scales_comma_label <- function() {
    if (requireNamespace("scales", quietly = TRUE)) {
        scales::comma
    } else {
        ggplot2::waiver()
    }
}
