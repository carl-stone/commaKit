#' @importFrom ggplot2 ggplot aes geom_line geom_vline scale_x_continuous
#'   scale_y_continuous labs theme_bw
NULL

#' Metagene plot of methylation across genomic features
#'
#' Computes average methylation beta values at normalized positions across
#' genomic features (e.g., genes from TSS to TTS) and plots the smoothed
#' profile. Useful for assessing whether methylation is enriched at particular
#' positions within a feature class.
#'
#' @param object A \code{\link{commaData}} object with a non-empty
#'   \code{annotation} slot or a user-supplied \code{features} GRanges.
#' @param feature Character string specifying the feature type to use as the
#'   reference. Must match a value in the \code{feature_type} metadata column
#'   of the annotation. Default \code{"gene"}.
#' @param mod_type Character string specifying a single modification type
#'   (e.g., \code{"6mA"}, \code{"5mC"}). If \code{NULL} (default), all
#'   modification types are used.
#' @param n_bins Positive integer. Number of equal-width bins to divide the
#'   normalized feature position \eqn{[0, 1]} into. Default \code{50}.
#'
#' @details
#' Internally calls \code{\link{annotateSites}(type = "metagene")} to compute
#' normalized positions (0 = TSS, 1 = TTS) for each methylation site that
#' overlaps a feature of the requested type. Sites that do not overlap any
#' feature are excluded from the plot. The mean beta value is then computed
#' within each position bin for each sample.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The x-axis shows normalized
#'   position within the feature (0 = TSS, 0.5 = midpoint, 1 = TTS); the
#'   y-axis shows mean methylation (beta). One line is drawn per sample,
#'   colored by sample name. Dashed vertical lines mark the TSS (0) and
#'   TTS (1).
#'
#' @examples
#' data(comma_example_data)
#' plot_metagene(comma_example_data, feature = "gene")
#'
#' # Only 6mA sites
#' plot_metagene(comma_example_data, feature = "gene", mod_type = "6mA")
#'
#' @seealso \code{\link{annotateSites}}, \code{\link{plot_genome_track}}
#'
#' @export
plot_metagene <- function(object,
                          feature  = "gene",
                          mod_type = NULL,
                          motif    = NULL,
                          n_bins   = 50L) {

    ## --- Input validation ---------------------------------------------------
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    if (!is.character(feature) || length(feature) != 1L) {
        stop("'feature' must be a single character string (e.g., 'gene').")
    }
    n_bins <- as.integer(n_bins)
    if (is.na(n_bins) || n_bins < 2L) {
        stop("'n_bins' must be an integer >= 2.")
    }

    ## --- Validate annotation ------------------------------------------------
    annot_gr <- annotation(object)
    if (is.null(annot_gr) || length(annot_gr) == 0L) {
        stop("annotation(object) is empty. ",
             "Provide a commaData object with annotation features, or annotate first with loadAnnotation().")
    }

    ## Filter to requested feature type
    if (!"feature_type" %in% colnames(S4Vectors::mcols(annot_gr))) {
        stop("annotation(object) does not have a 'feature_type' metadata column. ",
             "Use loadAnnotation() to build the annotation GRanges.")
    }
    feat_gr <- annot_gr[annot_gr$feature_type == feature]
    if (length(feat_gr) == 0L) {
        available_types <- unique(as.character(annot_gr$feature_type))
        stop("No features of type '", feature, "' found in annotation(object). ",
             "Available feature types: ", paste(available_types, collapse = ", "), ".")
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

    ## --- Run metagene annotation -------------------------------------------
    annotated <- annotateSites(object, features = feat_gr, type = "metagene")
    rd        <- as.data.frame(SummarizedExperiment::rowData(annotated))

    ## metagene_positions is a NumericList (one list element per site)
    pos_list <- rd$metagene_positions
    site_lengths <- lengths(pos_list)

    ## Keep only sites with at least one metagene position
    has_overlap <- site_lengths > 0L
    if (!any(has_overlap)) {
        stop("No methylation sites overlap any '", feature, "' features. ",
             "The metagene plot cannot be produced.")
    }

    ## Build expanded data.frame: one row per (site × feature) pair
    ## site_idx is 1-based into the filtered object
    site_idx     <- rep(which(has_overlap), times = site_lengths[has_overlap])
    metagene_pos <- unlist(pos_list[has_overlap], use.names = FALSE)

    ## --- Extract beta values for overlapping sites -------------------------
    methyl_mat <- methylation(object)
    sample_nms <- colnames(methyl_mat)
    n_samples  <- length(sample_nms)
    si         <- sampleInfo(object)

    ## Build per-sample long data.frame
    rows <- lapply(sample_nms, function(samp) {
        betas <- methyl_mat[site_idx, samp]
        data.frame(
            sample_name  = samp,
            metagene_pos = metagene_pos,
            beta         = betas,
            stringsAsFactors = FALSE
        )
    })
    long_df <- do.call(rbind, rows)

    ## Drop NAs
    long_df <- long_df[!is.na(long_df$beta), , drop = FALSE]

    if (nrow(long_df) == 0L) {
        stop("All overlapping sites have NA beta values. Cannot produce metagene plot.")
    }

    ## --- Bin and summarize --------------------------------------------------
    breaks    <- seq(0, 1, length.out = n_bins + 1L)
    bin_idx   <- findInterval(long_df$metagene_pos, breaks, rightmost.closed = TRUE)
    ## clamp edge cases to valid bin range
    bin_idx   <- pmax(1L, pmin(bin_idx, n_bins))
    bin_centers <- (breaks[-length(breaks)] + breaks[-1L]) / 2

    long_df$bin_center <- bin_centers[bin_idx]

    ## Compute mean beta per (sample, bin_center)
    mean_vals <- tapply(
        long_df$beta,
        list(long_df$sample_name, long_df$bin_center),
        mean, na.rm = TRUE
    )

    summary_rows <- lapply(sample_nms, function(samp) {
        bin_labels <- colnames(mean_vals)
        vals <- mean_vals[samp, ]
        data.frame(
            sample_name = samp,
            bin_center  = as.numeric(bin_labels),
            mean_beta   = as.numeric(vals),
            stringsAsFactors = FALSE
        )
    })
    summary_df <- do.call(rbind, summary_rows)
    summary_df <- summary_df[!is.na(summary_df$mean_beta), , drop = FALSE]

    ## Ensure consistent sample ordering
    summary_df$sample_name <- factor(summary_df$sample_name, levels = sample_nms)

    ## --- Build ggplot -------------------------------------------------------
    p <- ggplot2::ggplot(
        summary_df,
        ggplot2::aes(
            x     = .data[["bin_center"]],
            y     = .data[["mean_beta"]],
            color = .data[["sample_name"]]
        )
    ) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_vline(
            xintercept = c(0, 1),
            linetype = "dashed", color = "grey50", linewidth = 0.5
        ) +
        ggplot2::scale_x_continuous(
            breaks = c(0, 0.5, 1),
            labels = c("TSS", "Middle", "TTS"),
            name   = paste0("Relative position within '", feature, "'")
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, 1),
            name   = "Mean methylation"
        ) +
        ggplot2::labs(
            title = paste0("Metagene: ", feature),
            color = "Sample"
        ) +
        ggplot2::theme_bw()

    p
}
