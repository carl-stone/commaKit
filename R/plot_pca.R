#' @importFrom ggplot2 ggplot aes geom_point geom_text scale_shape_discrete
#'   labs theme_bw
NULL

#' PCA of methylation profiles
#'
#' Performs principal component analysis (PCA) on per-sample methylation
#' profiles and plots PC1 vs PC2. Useful for sample-level QC, detecting
#' outliers, and assessing whether biological conditions separate in
#' methylation space.
#'
#' @param object A \code{\link{commaData}} object.
#' @param mod_type Character string specifying a single modification type
#'   (e.g., \code{"6mA"}, \code{"5mC"}). If \code{NULL} (default), all sites
#'   from all modification types are used.
#' @param color_by Character string naming a column in \code{sampleInfo(object)}
#'   to use for point color. Default \code{"condition"}.
#' @param shape_by Character string naming a column in \code{sampleInfo(object)}
#'   to use for point shape. If \code{NULL} (default), all points use the same
#'   shape.
#'
#' @details
#' Beta values are first converted to M-values via \code{\link{mValues}}
#' (using \code{alpha = 0.5}) before PCA. M-values are variance-stabilized
#' relative to raw beta values, making distance-based analyses more reliable
#' especially when many sites are near 0 or 1. Sites with any \code{NA}
#' M-values across samples (including sites with zero coverage) are removed
#' to ensure a complete data matrix. PCA is computed via \code{stats::prcomp}
#' with centering (\code{center = TRUE}) and without scaling
#' (\code{scale. = FALSE}). A warning is issued if fewer than three samples
#' are present.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. PC1 and PC2 are shown on
#'   the x- and y-axes, respectively, with percentage of variance explained
#'   shown in the axis labels. Each point represents one sample and is labeled
#'   with its \code{sample_name}. Points are colored by \code{color_by}.
#'
#' @examples
#' data(comma_example_data)
#' plot_pca(comma_example_data)
#'
#' # Color by condition, shape by replicate
#' plot_pca(comma_example_data, color_by = "condition")
#'
#' # Only 6mA sites
#' plot_pca(comma_example_data, mod_type = "6mA")
#'
#' @seealso \code{\link{methylomeSummary}}, \code{\link{plot_methylation_distribution}}
#'
#' @export
plot_pca <- function(object,
                     mod_type  = NULL,
                     color_by  = "condition",
                     shape_by  = NULL) {

    ## --- Input validation ---------------------------------------------------
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
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

    ## --- Validate color_by / shape_by ---------------------------------------
    si <- sampleInfo(object)
    if (!color_by %in% colnames(si)) {
        stop("'color_by' column '", color_by, "' not found in sampleInfo(object). ",
             "Available columns: ", paste(colnames(si), collapse = ", "), ".")
    }
    if (!is.null(shape_by) && !shape_by %in% colnames(si)) {
        stop("'shape_by' column '", shape_by, "' not found in sampleInfo(object). ",
             "Available columns: ", paste(colnames(si), collapse = ", "), ".")
    }

    ## --- Build complete-case matrix (M-value transformed) -------------------
    methyl_mat <- mValues(object)     # variance-stabilized M-values
    n_samples  <- ncol(methyl_mat)

    if (n_samples < 2L) {
        stop("PCA requires at least 2 samples; object has ", n_samples, ".")
    }
    if (n_samples < 3L) {
        warning("Fewer than 3 samples available; PCA results may not be meaningful.")
    }

    ## Keep only sites with no NA M-values across all samples
    complete_sites <- which(rowSums(is.na(methyl_mat)) == 0L)
    if (length(complete_sites) < 2L) {
        stop("Fewer than 2 sites have non-NA M-values across all samples. ",
             "Cannot compute PCA. Try reducing 'min_coverage' in commaData() or ",
             "using a larger dataset.")
    }

    mat <- t(methyl_mat[complete_sites, , drop = FALSE])  # samples x sites

    ## --- Compute PCA --------------------------------------------------------
    pca     <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
    pct_var <- round(summary(pca)$importance[2L, seq_len(min(2L, ncol(pca$x)))] * 100,
                     digits = 1L)

    ## Build scores data.frame – use only however many PCs are available (≤ 2)
    n_pcs     <- min(2L, ncol(pca$x))
    scores_df <- as.data.frame(pca$x[, seq_len(n_pcs), drop = FALSE])
    colnames(scores_df)[seq_len(n_pcs)] <- c("PC1", "PC2")[seq_len(n_pcs)]
    ## Add a dummy PC2 column of zeros when only 1 PC was computed so that the
    ## ggplot2 aes mapping to "PC2" does not fail.
    if (n_pcs < 2L) {
        scores_df[["PC2"]] <- 0
    }
    scores_df$sample_name <- rownames(scores_df)

    ## Join sampleInfo columns
    scores_df <- merge(scores_df, si, by = "sample_name", all.x = TRUE)

    ## --- Build ggplot -------------------------------------------------------
    x_label <- if (length(pct_var) >= 1L) {
        paste0("PC1 (", pct_var[1L], "% variance)")
    } else "PC1"
    y_label <- if (length(pct_var) >= 2L) {
        paste0("PC2 (", pct_var[2L], "% variance)")
    } else "PC2"

    ## Build aes dynamically to support optional shape_by
    if (!is.null(shape_by)) {
        p <- ggplot2::ggplot(
            scores_df,
            ggplot2::aes(
                x     = .data[["PC1"]],
                y     = .data[["PC2"]],
                color = .data[[color_by]],
                shape = .data[[shape_by]],
                label = .data[["sample_name"]]
            )
        ) +
            ggplot2::geom_point(size = 3.5) +
            ggplot2::scale_shape_discrete(name = shape_by)
    } else {
        p <- ggplot2::ggplot(
            scores_df,
            ggplot2::aes(
                x     = .data[["PC1"]],
                y     = .data[["PC2"]],
                color = .data[[color_by]],
                label = .data[["sample_name"]]
            )
        ) +
            ggplot2::geom_point(size = 3.5)
    }

    p <- p +
        ggplot2::geom_text(
            vjust = -0.9, size = 3,
            show.legend = FALSE
        ) +
        ggplot2::labs(
            x     = x_label,
            y     = y_label,
            title = "PCA of Methylation Profiles",
            color = color_by
        ) +
        ggplot2::theme_bw()

    p
}
