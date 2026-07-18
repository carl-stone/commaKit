#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_vline
#'   scale_x_continuous scale_y_continuous facet_wrap labs theme_bw
#' @importFrom GenomicRanges GRanges findOverlaps mcols
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom SummarizedExperiment rowData
#' @importFrom methods is
NULL

#' TSS-centered methylation profile
#'
#' Plots individual methylation sites at their absolute base-pair position
#' relative to the nearest transcription start site (TSS), showing the raw
#' spatial distribution of methylation around gene starts. Optionally colours
#' sites by regulatory feature overlap (e.g., sigma factor binding sites,
#' \eqn{-10}/\eqn{-35} elements, transcription factor binding sites).
#'
#' @param object A \code{\link{commaData}} object with a non-empty
#'   \code{annotation} slot.
#' @param feature_type Character string specifying the feature type whose start
#'   coordinate defines the TSS (e.g., \code{"gene"}, \code{"CDS"}). Must
#'   match a value in the \code{feature_type} metadata column of
#'   \code{annotation(object)}. Default \code{"gene"}.
#' @param window Positive integer. Half-window size in base pairs around the
#'   TSS. Sites beyond \code{window} bp upstream or downstream of every TSS
#'   are excluded. Default \code{500L}.
#' @param regulatory_feature_types Character vector or \code{NULL}. Feature
#'   types from \code{annotation(object)} used to label sites when
#'   \code{color_by = "regulatory_element"} (e.g.,
#'   \code{c("sigma_binding", "promoter_-10", "promoter_-35",
#'   "TF_binding")}). Required when \code{color_by = "regulatory_element"};
#'   ignored otherwise.
#' @param mod_type Character vector or \code{NULL}. If provided, only sites
#'   of the specified modification type(s) (e.g., \code{"6mA"},
#'   \code{c("6mA", "5mC")}) are included.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   the specified sequence motif(s) are included.
#' @param mod_context Character vector or \code{NULL}. If provided, only sites
#'   whose \code{mod_context} rowData column matches one of the supplied values
#'   are included. If \code{NULL} (default), all modification contexts are
#'   included.
#' @param color_by Character string controlling the colour aesthetic:
#'   \describe{
#'     \item{\code{"sample"}}{One colour per sample (default).}
#'     \item{\code{"regulatory_element"}}{Sites coloured by the first
#'       regulatory feature they overlap (from \code{regulatory_feature_types});
#'       sites with no overlap are labelled \code{"None"} and shown in grey.
#'       Requires \code{regulatory_feature_types}.}
#'     \item{\code{"mod_type"}}{One colour per modification type.}
#'     \item{\code{"mod_context"}}{One colour per modification context.}
#'     \item{\code{"none"}}{No colour mapping; all points drawn in the default
#'       colour. Useful when \code{facet_by} provides sufficient grouping.}
#'   }
#' @param facet_by Character string controlling optional faceting:
#'   \code{"none"} (default), \code{"sample"}, \code{"mod_type"}, or
#'   \code{"mod_context"}.
#' @param alpha Numeric in \eqn{(0, 1]}. Point transparency. Default
#'   \code{0.4}.
#' @param show_smooth Logical. If \code{TRUE}, a ggplot2 loess smoothing line
#'   is overlaid per colour group within each facet. Sparse or unstable groups
#'   follow ggplot2's native smoothing behavior. Default \code{FALSE}.
#' @param smooth_span Numeric in \eqn{(0, 1]}. Loess span parameter passed
#'   to ggplot2's smoothing statistic. Default \code{0.3}.
#'
#' @details
#' Unlike \code{\link{plot_metagene}}, which normalises positions to
#' \eqn{[0, 1]} and averages across sites, \code{plot_tss_profile} shows
#' every individual site at its exact signed base-pair distance from the
#' nearest TSS (negative = upstream, positive = downstream). This preserves
#' the absolute spacing of promoter elements relative to the start of
#' transcription.
#'
#' Internally calls \code{\link{annotateSites}(keep = "proximity")} to
#' compute signed TSS distances. Strand awareness follows the same
#' convention: for \code{+} strand features, position 0 is the lowest
#' coordinate; for \code{-} strand features, position 0 is the highest
#' coordinate (the biological TSS). When a site is within \code{window} bp
#' of multiple TSS features, only the nearest (by absolute distance) is used.
#'
#' If \code{color_by = "regulatory_element"} and none of the specified
#' \code{regulatory_feature_types} are found in \code{annotation(object)},
#' the function issues a message and falls back to \code{color_by = "sample"}
#' so the plot still renders.
#'
#' Colour mappings use ggplot2's default discrete scale. This keeps colour
#' assignment and smoothing grouping in one ggplot2-native contract; colours
#' are not lightened or darkened separately for points and smooths. When
#' \code{show_smooth = TRUE}, \code{geom_smooth()} fits a loess curve for each
#' mapped colour group in each facet. Sparse or numerically unstable groups
#' are handled by ggplot2 and may produce a warning or no rendered curve.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The x-axis shows signed
#'   position relative to the TSS in base pairs; the y-axis shows methylation
#'   beta value (0--1). A dashed vertical line marks the TSS at x = 0.
#'
#' @examples
#' data(comma_example_data)
#' plot_tss_profile(comma_example_data, feature_type = "gene", window = 500L)
#'
#' # Colour by modification type, facet by sample
#' plot_tss_profile(comma_example_data,
#'   feature_type = "gene",
#'   color_by = "mod_type", facet_by = "sample"
#' )
#'
#' # Overlay loess smooth
#' plot_tss_profile(comma_example_data,
#'   feature_type = "gene",
#'   show_smooth = TRUE
#' )
#'
#' @seealso \code{\link{annotateSites}}, \code{\link{plot_metagene}},
#'   \code{\link{loadAnnotation}}
#'
#' @export
plot_tss_profile <- function(object,
                             feature_type = "gene",
                             window = 500L,
                             regulatory_feature_types = NULL,
                             mod_type = NULL,
                             motif = NULL,
                             mod_context = NULL,
                             color_by = c(
                               "sample",
                               "regulatory_element",
                               "mod_type",
                               "mod_context",
                               "none"
                             ),
                             facet_by = c(
                               "none", "sample",
                               "mod_type",
                               "mod_context"
                             ),
                             alpha = 0.4,
                             show_smooth = FALSE,
                             smooth_span = 0.3) {
  ## ── A. Input validation ──────────────────────────────────────────────────
  if (!is(object, "commaData")) {
    stop("'object' must be a commaData object.")
  }
  if (!is.character(feature_type) || length(feature_type) != 1L ||
      is.na(feature_type)) {
    stop("'feature_type' must be a single non-NA character string.")
  }
  window <- as.integer(window)
  if (is.na(window) || window < 1L) {
    stop("'window' must be a positive integer (base pairs).")
  }
  color_by <- match.arg(color_by)
  facet_by <- match.arg(facet_by)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha > 1) {
    stop("'alpha' must be a single numeric value in (0, 1].")
  }
  if (!is.logical(show_smooth) || length(show_smooth) != 1L ||
      is.na(show_smooth)) {
    stop("'show_smooth' must be TRUE or FALSE.")
  }
  if (!is.numeric(smooth_span) || length(smooth_span) != 1L ||
      is.na(smooth_span) || smooth_span <= 0 || smooth_span > 1) {
    stop("'smooth_span' must be a single numeric value in (0, 1].")
  }
  if (color_by == "regulatory_element" && is.null(regulatory_feature_types)) {
    stop(
      "'color_by = \"regulatory_element\"' requires ",
      "'regulatory_feature_types' to be specified."
    )
  }

  ## ── B. Annotation validation ─────────────────────────────────────────────
  annot_gr <- annotation(object)
  if (is.null(annot_gr) || length(annot_gr) == 0L) {
    stop(
      "annotation(object) is empty. ",
      "Provide annotation features via loadAnnotation()."
    )
  }
  if (!"feature_type" %in% colnames(S4Vectors::mcols(annot_gr))) {
    stop(
      "annotation(object) does not have a 'feature_type' metadata ",
      "column. Use loadAnnotation() to build the annotation GRanges."
    )
  }
  tss_gr <- annot_gr[annot_gr$feature_type == feature_type]
  if (length(tss_gr) == 0L) {
    available_types <- unique(as.character(annot_gr$feature_type))
    stop(
      "No features of type '", feature_type,
      "' found in annotation(object). ",
      "Available feature types: ",
      paste(available_types, collapse = ", "), "."
    )
  }

  ## Narrow features to 1-bp TSS points so that rel_positions from
  ## annotateSites(type = "proximity") measures distance from the TSS itself
  ## rather than from the nearest feature boundary.
  feat_strand <- as.character(GenomicRanges::strand(tss_gr))
  tss_pos <- ifelse(feat_strand == "-",
    GenomicRanges::end(tss_gr),
    GenomicRanges::start(tss_gr)
  )
  tss_1bp <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(tss_gr),
    ranges   = IRanges::IRanges(start = tss_pos, width = 1L),
    strand   = GenomicRanges::strand(tss_gr)
  )
  GenomicRanges::mcols(tss_1bp) <- GenomicRanges::mcols(tss_gr)
  tss_gr <- tss_1bp

  ## ── C. Filter by mod_type / motif / mod_context ──────────────────────────
  object <- .applySiteFilters(
    object,
    mod_type = mod_type,
    motif = motif,
    mod_context = mod_context,
    caller = "plot_tss_profile()"
  )

  ## ── D. Proximity annotation ───────────────────────────────────────────────
  ## annotateSites(keep = "proximity") on 1-bp TSS GRanges adds rel_position
  ## (IntegerList): negative = upstream of TSS, positive = downstream.
  annotated <- annotateSites(object,
    features = tss_gr, keep = "proximity",
    window = window
  )
  rd <- as.data.frame(siteInfo(annotated))

  ## ── E. Nearest TSS offset per site (vectorized) ───────────────────────────
  pos_list <- rd$rel_position
  has_nearby <- lengths(pos_list) > 0L

  nearest_rel_pos <- rep(NA_integer_, nrow(rd))
  nearest_rel_pos[has_nearby] <- vapply(
    pos_list[has_nearby],
    function(v) v[which.min(abs(v))],
    integer(1L)
  )
  ## Exclude sites whose nearest TSS is beyond the requested window.
  ## annotateSites proximity finds features within 'window' bp of the site
  ## boundary, but a site inside a gene body can have a large signed distance
  ## to the TSS. Filter those out here.
  nearest_rel_pos[abs(nearest_rel_pos) > window] <- NA_integer_

  keep <- !is.na(nearest_rel_pos)
  if (!any(keep)) {
    stop(
      "No methylation sites found within ", window, " bp of any '",
      feature_type, "' TSS. Consider increasing 'window'."
    )
  }
  annotated_sub <- annotated[keep, ]
  rel_pos_vals <- nearest_rel_pos[keep]
  rd_sub <- as.data.frame(siteInfo(annotated_sub))

  ## ── F. Build long data.frame (sites × samples) ────────────────────────────
  methyl_mat <- methylation(annotated_sub)
  sample_nms <- colnames(methyl_mat)
  n_sites_sub <- nrow(methyl_mat)
  n_samples <- length(sample_nms)

  df <- data.frame(
    rel_pos = rep(rel_pos_vals, times = n_samples),
    beta = as.vector(methyl_mat),
    sample_name = rep(sample_nms, each = n_sites_sub),
    mod_type_col = rep(rd_sub$mod_type, times = n_samples),
    mod_context = rep(rd_sub$mod_context, times = n_samples),
    chrom = rep(rd_sub$chrom, times = n_samples),
    position = rep(rd_sub$position, times = n_samples),
    stringsAsFactors = FALSE
  )
  ## Rename mod_type_col → mod_type (avoid clash with parameter name)
  names(df)[names(df) == "mod_type_col"] <- "mod_type"

  df_keep <- !is.na(df$beta)
  df <- df[df_keep, , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("All sites within the window have NA beta values.")
  }

  ## ── G. Regulatory element coloring ───────────────────────────────────────
  if (color_by == "regulatory_element") {
    reg_gr <- annot_gr[annot_gr$feature_type %in% regulatory_feature_types]

    if (length(reg_gr) == 0L) {
      available_types <- unique(as.character(annot_gr$feature_type))
      message(
        "None of the 'regulatory_feature_types' (",
        paste(regulatory_feature_types, collapse = ", "),
        ") found in annotation(object). ",
        "Available feature types: ",
        paste(available_types, collapse = ", "),
        ". Falling back to color_by = 'sample'."
      )
      color_by <- "sample"
    } else {
      sites_gr <- GenomicRanges::GRanges(
        seqnames = rd_sub$chrom,
        ranges   = IRanges::IRanges(start = rd_sub$position, width = 1L),
        strand   = rd_sub$strand
      )
      reg_hits <- GenomicRanges::findOverlaps(sites_gr, reg_gr,
        ignore.strand = TRUE
      )
      reg_label <- rep("None", nrow(rd_sub))
      if (length(reg_hits) > 0L) {
        q_idx <- S4Vectors::queryHits(reg_hits)
        s_idx <- S4Vectors::subjectHits(reg_hits)
        ## First match per site (annotation order)
        first_match <- !duplicated(q_idx)
        reg_label[q_idx[first_match]] <-
          as.character(reg_gr$feature_type[s_idx[first_match]])
      }
      ## Factor: regulatory types alphabetically, then "None" last
      reg_levels <- c(
        sort(unique(reg_label[reg_label != "None"])),
        "None"
      )
      reg_label <- factor(reg_label, levels = reg_levels)
      ## Expand across samples
      df$regulatory_element <- rep(reg_label, times = n_samples)[df_keep]
    }
  }

  ## ── H. Map color variable name ────────────────────────────────────────────
  color_var <- switch(color_by,
    sample             = "sample_name",
    regulatory_element = "regulatory_element",
    mod_type           = "mod_type",
    mod_context        = "mod_context",
    none               = NULL
  )

  ## ── I. Build ggplot ───────────────────────────────────────────────────────
  base_aes <- if (!is.null(color_var)) {
    ggplot2::aes(
      x     = .data[["rel_pos"]],
      y     = .data[["beta"]],
      color = .data[[color_var]]
    )
  } else {
    ggplot2::aes(
      x = .data[["rel_pos"]],
      y = .data[["beta"]]
    )
  }
  p <- ggplot2::ggplot(df, base_aes) +
    ggplot2::geom_point(
      alpha = alpha, size = 0.8,
      color = if (is.null(color_var)) "black" else NULL
    ) +
    ggplot2::geom_vline(
      xintercept = 0L,
      linetype   = "dashed",
      color      = "grey40",
      linewidth  = 0.5
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-window, window),
      name = paste0(
        "Position relative to TSS (bp)\n[",
        feature_type, "]"
      ),
      labels = function(x) {
        ifelse(x > 0, paste0("+", x),
          ifelse(x == 0, "TSS", as.character(x))
        )
      }
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      name   = "Methylation (beta)"
    ) +
    ggplot2::labs(
      title = paste0(
        "TSS methylation profile: ", feature_type,
        " (\u00b1", window, " bp)"
      ),
      color = if (!is.null(color_var)) {
        switch(color_by,
          sample             = "Sample",
          regulatory_element = "Regulatory element",
          mod_type           = "Modification type",
          mod_context        = "Modification context"
        )
      } else {
        NULL
      }
    ) +
    ggplot2::theme_bw()

  ## ── J. Optional loess smooth overlay ─────────────────────────────────────
  if (show_smooth) {
    p <- p + ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      span = smooth_span,
      se = FALSE,
      na.rm = TRUE,
      linewidth = 1.0
    )
  }

  ## ── K. Faceting ──────────────────────────────────────────────────────────
  if (facet_by == "sample") {
    p <- p + ggplot2::facet_wrap("sample_name", ncol = 1L)
  } else if (facet_by == "mod_type") {
    p <- p + ggplot2::facet_wrap("mod_type")
  } else if (facet_by == "mod_context") {
    p <- p + ggplot2::facet_wrap("mod_context")
  }

  p
}
