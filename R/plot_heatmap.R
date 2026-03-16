#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 scale_x_discrete
#'   scale_y_discrete labs theme_minimal theme element_text element_blank
NULL

#' Heatmap of top differentially methylated sites
#'
#' Produces a heatmap showing methylation beta values for the top
#' differentially methylated sites across all samples. Sites are ranked by
#' adjusted p-value and ordered vertically by effect size (delta beta) to
#' reveal condition-specific patterns.
#'
#' @param results A \code{data.frame} returned by \code{\link{results}()},
#'   containing at minimum the columns \code{chrom}, \code{position},
#'   \code{strand}, \code{mod_type}, \code{dm_padj}, and
#'   \code{dm_delta_beta}.
#' @param object A \code{\link{commaData}} object that was used to produce
#'   \code{results}. Used to extract the methylation matrix for selected sites.
#' @param n_sites Positive integer. The number of top sites (ranked by
#'   ascending \code{dm_padj}) to include in the heatmap. If fewer significant
#'   (non-\code{NA} padj) sites exist, all are shown. Default \code{50}.
#' @param annotation_cols Character vector naming columns from
#'   \code{sampleInfo(object)} to display as a colored annotation bar above
#'   the heatmap. If \code{NULL} (default), the \code{condition} column is
#'   used.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The x-axis shows sample
#'   names; the y-axis shows the selected sites (y-axis labels suppressed for
#'   readability). The fill color encodes methylation beta (blue = 0,
#'   white = 0.5, red = 1). \code{NA} values are shown in light grey. An
#'   annotation strip below shows sample-level metadata encoded by color.
#'
#' @examples
#' data(comma_example_data)
#' cd_dm <- diffMethyl(comma_example_data, ~ condition)
#' res   <- results(cd_dm)
#' plot_heatmap(res, cd_dm)
#'
#' # Show only top 20 sites
#' plot_heatmap(res, cd_dm, n_sites = 20)
#'
#' @seealso \code{\link{diffMethyl}}, \code{\link{results}},
#'   \code{\link{plot_volcano}}
#'
#' @export
plot_heatmap <- function(results,
                         object,
                         n_sites        = 50L,
                         annotation_cols = NULL) {

    ## --- Input validation ---------------------------------------------------
    if (!is.data.frame(results)) {
        stop("'results' must be a data.frame returned by results().")
    }
    required_cols <- c("chrom", "position", "strand", "mod_type",
                       "dm_padj", "dm_delta_beta")
    missing_cols <- setdiff(required_cols, colnames(results))
    if (length(missing_cols) > 0L) {
        stop("'results' is missing required column(s): ",
             paste(missing_cols, collapse = ", "), ". ",
             "Ensure 'results' was produced by results() after diffMethyl().")
    }
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    n_sites <- as.integer(n_sites)
    if (is.na(n_sites) || n_sites < 1L) {
        stop("'n_sites' must be a positive integer.")
    }

    ## --- Select top n_sites by padj -----------------------------------------
    res_nonNA <- results[!is.na(results$dm_padj), , drop = FALSE]
    if (nrow(res_nonNA) == 0L) {
        stop("No rows with non-NA 'dm_padj' in 'results'. ",
             "Run diffMethyl() first.")
    }
    res_sorted <- res_nonNA[order(res_nonNA$dm_padj), , drop = FALSE]
    n_use <- min(n_sites, nrow(res_sorted))
    top_res <- res_sorted[seq_len(n_use), , drop = FALSE]

    ## --- Match to object rowData --------------------------------------------
    make_key <- function(chrom, position, strand, mod_type) {
        paste(chrom, position, strand, mod_type, sep = ":")
    }

    rd        <- siteInfo(object)
    obj_keys  <- make_key(rd$chrom, rd$position, rd$strand, rd$mod_type)
    res_keys  <- make_key(top_res$chrom, top_res$position,
                          top_res$strand, top_res$mod_type)

    row_idx <- match(res_keys, obj_keys)
    missing_keys <- sum(is.na(row_idx))
    if (missing_keys == n_use) {
        stop("None of the selected sites were found in 'object'. ",
             "Ensure 'results' was produced from the same commaData object.")
    }
    if (missing_keys > 0L) {
        warning(missing_keys, " site(s) from 'results' could not be matched in 'object' and will be excluded.")
        top_res  <- top_res[!is.na(row_idx), , drop = FALSE]
        row_idx  <- row_idx[!is.na(row_idx)]
    }

    ## --- Extract methylation submatrix -------------------------------------
    methyl_mat <- methylation(object)[row_idx, , drop = FALSE]
    sample_nms <- colnames(methyl_mat)

    ## Order sites by dm_delta_beta for a meaningful visual grouping
    delta_order <- order(top_res$dm_delta_beta)
    methyl_mat  <- methyl_mat[delta_order, , drop = FALSE]
    top_res     <- top_res[delta_order, , drop = FALSE]

    ## Build site keys for y-axis ordering
    site_keys <- make_key(top_res$chrom, top_res$position,
                          top_res$strand, top_res$mod_type)
    n_final <- length(site_keys)

    ## --- Reshape to long data.frame ----------------------------------------
    df <- data.frame(
        site_key    = rep(site_keys, times = length(sample_nms)),
        sample_name = rep(sample_nms, each  = n_final),
        beta        = as.vector(methyl_mat),
        stringsAsFactors = FALSE
    )
    df$site_key <- factor(df$site_key, levels = site_keys)

    ## --- Sample annotation strip -------------------------------------------
    si <- sampleInfo(object)

    if (is.null(annotation_cols)) {
        annotation_cols <- "condition"
    }
    annotation_cols <- intersect(annotation_cols, colnames(si))

    ## --- Build ggplot -------------------------------------------------------
    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
            x    = .data[["sample_name"]],
            y    = .data[["site_key"]],
            fill = .data[["beta"]]
        )
    ) +
        ggplot2::geom_tile(color = "white", linewidth = 0.1) +
        ggplot2::scale_fill_gradient2(
            low      = "#4575b4",
            mid      = "white",
            high     = "#d73027",
            midpoint = 0.5,
            limits   = c(0, 1),
            na.value = "grey85",
            name     = "Methylation"
        ) +
        ggplot2::scale_y_discrete(limits = site_keys) +
        ggplot2::scale_x_discrete(limits = sample_nms) +
        ggplot2::labs(
            x     = NULL,
            y     = paste0("Top ", n_final, " differential sites\n(ordered by delta methylation)"),
            title = "Differential Methylation Heatmap"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y  = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            panel.grid   = ggplot2::element_blank(),
            axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1)
        )

    ## If annotation_cols provided, append an annotation strip below
    if (length(annotation_cols) > 0L) {
        ## Build long annotation data.frame
        annot_rows <- lapply(annotation_cols, function(col) {
            data.frame(
                sample_name   = si$sample_name,
                annot_value   = as.character(si[[col]]),
                annot_var     = col,
                stringsAsFactors = FALSE
            )
        })
        annot_df <- do.call(rbind, annot_rows)
        annot_df$sample_name <- factor(annot_df$sample_name, levels = sample_nms)

        p_annot <- ggplot2::ggplot(
            annot_df,
            ggplot2::aes(
                x    = .data[["sample_name"]],
                y    = .data[["annot_var"]],
                fill = .data[["annot_value"]]
            )
        ) +
            ggplot2::geom_tile(color = "white", linewidth = 0.3) +
            ggplot2::scale_x_discrete(limits = sample_nms) +
            ggplot2::labs(x = NULL, y = NULL, fill = "Sample\nannotation") +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
                panel.grid   = ggplot2::element_blank()
            )

        if (requireNamespace("patchwork", quietly = TRUE)) {
            p <- patchwork::wrap_plots(
                p_annot, p,
                ncol    = 1L,
                heights = c(length(annotation_cols), 10)
            )
        } else {
            message("Install the 'patchwork' package to display sample annotation strips. ",
                    "Returning heatmap without annotation strip.")
        }
    }

    p
}
