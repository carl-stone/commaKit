#' @importFrom ggplot2 ggplot aes geom_density scale_x_continuous facet_wrap
#'   labs theme_bw
NULL

#' Plot methylation beta value distributions
#'
#' Produces a density plot of methylation beta values (0--1) for each sample
#' in a \code{\link{commaData}} object. Useful for QC and for comparing
#' methylation level distributions across samples and modification types.
#'
#' @param object A \code{\link{commaData}} object.
#' @param mod_type Character string specifying a single modification type to
#'   plot (e.g., \code{"6mA"}, \code{"5mC"}). If \code{NULL} (default), all
#'   modification types are included and the plot is faceted by
#'   \code{mod_type}.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   matching sequence context motif(s) are included (e.g., \code{"GATC"}).
#'   If \code{NULL} (default), all motifs are included.
#' @param per_sample Logical. If \code{TRUE} (default), a separate density
#'   curve is drawn for each sample. If \code{FALSE}, a single aggregate
#'   density curve is drawn per modification type.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The x-axis shows beta
#'   values (0 = unmethylated, 1 = fully methylated); the y-axis shows kernel
#'   density. When \code{per_sample = TRUE}, curves are colored by
#'   \code{sample_name}. When multiple modification types are present (and
#'   \code{mod_type = NULL}), the plot is faceted by \code{mod_type}. Sites
#'   with \code{NA} beta values (below coverage threshold) are silently
#'   excluded.
#'
#' @examples
#' data(comma_example_data)
#' plot_methylation_distribution(comma_example_data)
#'
#' # One modification type only
#' plot_methylation_distribution(comma_example_data, mod_type = "6mA")
#'
#' # Aggregate across samples
#' plot_methylation_distribution(comma_example_data, per_sample = FALSE)
#'
#' @seealso \code{\link{methylomeSummary}}, \code{\link{plot_coverage}}
#'
#' @export
plot_methylation_distribution <- function(object,
                                          mod_type = NULL,
                                          motif    = NULL,
                                          per_sample = TRUE) {

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

    ## --- Extract data -------------------------------------------------------
    methyl_mat  <- methylation(object)
    si          <- sampleInfo(object)
    site_info   <- siteInfo(object)
    sample_nms  <- colnames(methyl_mat)
    n_sites     <- nrow(methyl_mat)
    n_samples   <- length(sample_nms)

    ## Reshape to long data.frame (vectorized)
    df <- data.frame(
        beta        = as.vector(methyl_mat),
        sample_name = rep(sample_nms, each = n_sites),
        mod_type    = rep(site_info$mod_type, times = n_samples),
        stringsAsFactors = FALSE
    )

    ## Drop NA beta values (sites below min_coverage)
    df <- df[!is.na(df$beta), , drop = FALSE]

    if (nrow(df) == 0L) {
        stop("No non-NA methylation values found after filtering. ",
             "Check coverage thresholds.")
    }

    ## Join condition from sampleInfo
    si_sub <- si[, c("sample_name", "condition"), drop = FALSE]
    df <- merge(df, si_sub, by = "sample_name", all.x = TRUE)

    ## --- Build ggplot -------------------------------------------------------
    multi_mod <- length(unique(df$mod_type)) > 1L

    if (per_sample) {
        p <- ggplot2::ggplot(
            df,
            ggplot2::aes(
                x    = .data[["beta"]],
                color = .data[["sample_name"]],
                fill  = .data[["sample_name"]]
            )
        ) +
            ggplot2::geom_density(alpha = 0.3) +
            ggplot2::labs(color = "Sample", fill = "Sample")
    } else {
        p <- ggplot2::ggplot(
            df,
            ggplot2::aes(x = .data[["beta"]])
        ) +
            ggplot2::geom_density(fill = "steelblue", alpha = 0.4, color = "steelblue4")
    }

    p <- p +
        ggplot2::scale_x_continuous(
            limits = c(0, 1),
            expand = c(0.01, 0.01),
            name   = "Methylation"
        ) +
        ggplot2::labs(
            y     = "Density",
            title = "Methylation Beta Distribution"
        ) +
        ggplot2::theme_bw()

    if (multi_mod) {
        p <- p + ggplot2::facet_wrap("mod_type")
    }

    p
}
