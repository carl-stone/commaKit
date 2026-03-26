#' @importFrom ggplot2 ggplot aes geom_point geom_rect scale_color_manual
#'   scale_x_continuous scale_y_continuous facet_wrap labs theme_bw theme
#'   element_text element_blank
NULL

#' Genome browser–style methylation track plot
#'
#' Plots methylation beta values for individual sites along a chromosome
#' region in a genome browser–style layout, with one panel per sample.
#' Optionally overlays genomic feature annotations as colored rectangles
#' in a separate track below the methylation data.
#'
#' @param object A \code{\link{commaData}} object.
#' @param chromosome Character string. The chromosome (sequence name) to plot.
#'   Must be present in \code{names(genome(object))}.
#' @param start Integer or \code{NULL}. Start position of the region to display
#'   (1-based, inclusive). If \code{NULL}, the plot begins at position 1.
#' @param end Integer or \code{NULL}. End position of the region to display
#'   (1-based, inclusive). If \code{NULL}, the plot extends to the end of the
#'   chromosome.
#' @param mod_type Character string specifying a single modification type to
#'   display (e.g., \code{"6mA"}, \code{"5mC"}). If \code{NULL} (default),
#'   all modification types are shown, colored differently.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   matching sequence context motif(s) are displayed (e.g., \code{"GATC"}).
#'   If \code{NULL} (default), all motifs are included.
#' @param annotation A \code{GRanges} object of genomic features to display in
#'   the annotation track, \code{NULL} (default, uses \code{annotation(object)}
#'   if available), or \code{FALSE} to suppress the annotation track entirely.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The methylation track shows
#'   individual sites as points (x = genomic position, y = beta value 0--1)
#'   colored by modification type, faceted by sample. If annotation features
#'   are present on the selected chromosome (and \code{annotation} is not
#'   \code{FALSE}), they are displayed as colored rectangles in a separate
#'   annotation panel below.
#'
#' @examples
#' data(comma_example_data)
#' plot_genome_track(comma_example_data, chromosome = "chr_sim")
#'
#' # Restrict to a region
#' plot_genome_track(comma_example_data, chromosome = "chr_sim",
#'                   start = 1000, end = 50000)
#'
#' # One modification type, no annotation
#' plot_genome_track(comma_example_data, chromosome = "chr_sim",
#'                   mod_type = "6mA", annotation = FALSE)
#'
#' @seealso \code{\link{annotateSites}}, \code{\link{plot_metagene}}
#'
#' @export
plot_genome_track <- function(object,
                              chromosome,
                              start      = NULL,
                              end        = NULL,
                              mod_type   = NULL,
                              motif      = NULL,
                              annotation = NULL) {

    ## --- Input validation ---------------------------------------------------
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    if (!is.character(chromosome) || length(chromosome) != 1L) {
        stop("'chromosome' must be a single character string.")
    }
    genome_info <- genome(object)
    if (!is.null(genome_info) && length(genome_info) > 0L &&
        !chromosome %in% names(genome_info)) {
        stop("'chromosome' = '", chromosome, "' not found in genome(object). ",
             "Available chromosomes: ", paste(names(genome_info), collapse = ", "), ".")
    }
    if (!is.null(start) && (!is.numeric(start) || length(start) != 1L || start < 1)) {
        stop("'start' must be a single positive integer or NULL.")
    }
    if (!is.null(end) && (!is.numeric(end) || length(end) != 1L || end < 1)) {
        stop("'end' must be a single positive integer or NULL.")
    }
    if (!is.null(start) && !is.null(end) && end < start) {
        stop("'end' must be >= 'start'.")
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

    ## --- Extract and filter site data ---------------------------------------
    si_df      <- siteInfo(object)
    methyl_mat <- methylation(object)
    sample_nms <- colnames(methyl_mat)
    n_samples  <- length(sample_nms)
    n_sites    <- nrow(methyl_mat)

    ## Filter to chromosome
    chr_idx <- which(si_df$chrom == chromosome)
    if (length(chr_idx) == 0L) {
        stop("No methylation sites found on chromosome '", chromosome, "'.")
    }

    ## Apply start/end range
    pos_vals <- si_df$position[chr_idx]
    if (!is.null(start)) {
        keep <- pos_vals >= start
        chr_idx <- chr_idx[keep]
        pos_vals <- pos_vals[keep]
    }
    if (!is.null(end)) {
        keep <- pos_vals <= end
        chr_idx <- chr_idx[keep]
        pos_vals <- pos_vals[keep]
    }

    if (length(chr_idx) == 0L) {
        stop("No methylation sites found in the specified region of '", chromosome, "'.")
    }

    ## Build long data.frame: one row per (site, sample)
    site_sub   <- si_df[chr_idx, , drop = FALSE]
    methyl_sub <- methyl_mat[chr_idx, , drop = FALSE]

    df <- data.frame(
        position    = rep(site_sub$position, times = n_samples),
        beta        = as.vector(methyl_sub),
        sample_name = rep(sample_nms, each = length(chr_idx)),
        mod_type    = rep(site_sub$mod_type, times = n_samples),
        panel       = "Methylation",
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$beta), , drop = FALSE]

    ## --- Resolve annotation -------------------------------------------------
    use_annotation <- FALSE
    annot_gr <- NULL

    if (isFALSE(annotation)) {
        use_annotation <- FALSE
    } else if (is.null(annotation)) {
        annot_gr <- annotation(object)
        use_annotation <- !is.null(annot_gr) && length(annot_gr) > 0L
    } else if (inherits(annotation, "GRanges")) {
        annot_gr <- annotation
        use_annotation <- length(annot_gr) > 0L
    } else {
        stop("'annotation' must be NULL, FALSE, or a GRanges object.")
    }

    ## Filter annotation to this chromosome
    if (use_annotation) {
        chr_annot <- annot_gr[as.character(GenomicRanges::seqnames(annot_gr)) == chromosome]
        if (!is.null(start)) {
            chr_annot <- chr_annot[GenomicRanges::end(chr_annot) >= start]
        }
        if (!is.null(end)) {
            chr_annot <- chr_annot[GenomicRanges::start(chr_annot) <= end]
        }
        use_annotation <- length(chr_annot) > 0L
        if (use_annotation) annot_gr <- chr_annot
    }

    ## --- Color palette for mod_type -----------------------------------------
    mod_colors <- c("6mA" = "#e41a1c", "5mC" = "#377eb8", "4mC" = "#4daf4a")
    present_mods <- unique(df$mod_type)
    missing_mods <- setdiff(present_mods, names(mod_colors))
    for (m in missing_mods) {
        mod_colors[m] <- "grey50"
    }

    ## --- Build methylation track --------------------------------------------
    ## sample_name ordering consistent with sampleInfo order
    df$sample_name <- factor(df$sample_name, levels = sample_nms)

    region_start <- if (!is.null(start)) start else min(df$position)
    region_end   <- if (!is.null(end))   end   else max(df$position)

    p_meth <- ggplot2::ggplot(
        df,
        ggplot2::aes(
            x     = .data[["position"]],
            y     = .data[["beta"]],
            color = .data[["mod_type"]]
        )
    ) +
        ggplot2::geom_point(alpha = 0.7, size = 0.9) +
        ggplot2::scale_color_manual(
            values = mod_colors,
            name   = "Mod type"
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, 1),
            name   = "Methylation (beta)"
        ) +
        ggplot2::scale_x_continuous(
            limits = c(region_start, region_end),
            labels = function(x) format(x, big.mark = ",", scientific = FALSE)
        ) +
        ggplot2::facet_wrap("sample_name", ncol = 1L) +
        ggplot2::labs(
            x     = paste0("Position on ", chromosome),
            title = paste0("Genome track: ", chromosome)
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(strip.text = ggplot2::element_text(size = 9))

    ## If no annotation, return the methylation plot directly
    if (!use_annotation) {
        return(p_meth)
    }

    ## --- Build annotation track as a secondary ggplot -----------------------
    ## Represent features as rectangles; use feature_type for color
    feat_df <- data.frame(
        feat_start = GenomicRanges::start(annot_gr),
        feat_end   = GenomicRanges::end(annot_gr),
        feature_type = if ("feature_type" %in% colnames(S4Vectors::mcols(annot_gr))) {
            as.character(annot_gr$feature_type)
        } else {
            rep("feature", length(annot_gr))
        },
        feat_name = if ("name" %in% colnames(S4Vectors::mcols(annot_gr))) {
            as.character(annot_gr$name)
        } else {
            rep("", length(annot_gr))
        },
        stringsAsFactors = FALSE
    )

    ## Assign y positions (stagger overlapping features if needed — simple 1-row)
    feat_df$ymin <- 0.1
    feat_df$ymax <- 0.9

    p_annot <- ggplot2::ggplot(
        feat_df,
        ggplot2::aes(
            xmin = .data[["feat_start"]],
            xmax = .data[["feat_end"]],
            ymin = .data[["ymin"]],
            ymax = .data[["ymax"]],
            fill = .data[["feature_type"]]
        )
    ) +
        ggplot2::geom_rect(color = "grey30", linewidth = 0.3, alpha = 0.7) +
        ggplot2::scale_x_continuous(
            limits = c(region_start, region_end),
            labels = function(x) format(x, big.mark = ",", scientific = FALSE)
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1), breaks = NULL) +
        ggplot2::labs(
            x    = paste0("Position on ", chromosome),
            y    = "Annotation",
            fill = "Feature type"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.y  = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            panel.grid   = ggplot2::element_blank()
        )

    ## Combine tracks using patchwork if available; otherwise return meth plot
    ## with a message about annotation.
    if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- patchwork::wrap_plots(
            p_meth, p_annot,
            ncol   = 1L,
            heights = c(3, 1)
        )
        return(combined)
    } else {
        message("Install the 'patchwork' package to display the annotation track below the methylation track. ",
                "Returning methylation track only.")
        return(p_meth)
    }
}
