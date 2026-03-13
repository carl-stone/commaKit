#' @importFrom GenomicRanges GRanges findOverlaps distanceToNearest start end
#'   width strand
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame queryHits subjectHits
#' @importFrom SummarizedExperiment rowData
#' @importFrom methods is
NULL

#' Annotate methylation sites relative to genomic features
#'
#' Assigns genomic feature annotations to methylation sites stored in a
#' \code{\link{commaData}} object using vectorized
#' \code{\link[GenomicRanges]{findOverlaps}} queries — no nested for-loops.
#' Three annotation modes are available: \code{"overlap"} assigns feature
#' identity to each site, \code{"proximity"} reports the distance and
#' signed offset to the nearest feature boundary, and \code{"metagene"}
#' reports the fractional position within the overlapping feature.
#'
#' @param object A \code{\link{commaData}} object.
#' @param features A \code{\link[GenomicRanges]{GRanges}} of genomic features
#'   to annotate against. If \code{NULL} (default), the annotation stored in
#'   \code{object} via \code{\link{annotation}(object)} is used. Must have
#'   mcols columns named by \code{feature_col} and \code{name_col}.
#' @param type Character string specifying the annotation mode. One of:
#'   \describe{
#'     \item{\code{"overlap"}}{(default) Each site is assigned the feature type
#'       and name of the overlapping feature. Sites that overlap no feature
#'       receive \code{"intergenic"}.}
#'     \item{\code{"proximity"}}{Each site is assigned the nearest feature, its
#'       distance, and a signed relative position (negative = upstream of the
#'       feature start; positive = downstream).}
#'     \item{\code{"metagene"}}{Each site that overlaps a feature is assigned
#'       a fractional position within that feature (0 = feature start, 1 =
#'       feature end). Strand-aware: for \code{"-"} strand features, 0 is at
#'       the feature end (highest coordinate) and 1 is at the feature start
#'       (lowest coordinate). Non-overlapping sites receive \code{NA}.}
#'   }
#' @param feature_col Character string. Name of the \code{mcols} column in
#'   \code{features} that contains the feature type (e.g.,
#'   \code{"feature_type"}). Default: \code{"feature_type"}.
#' @param name_col Character string. Name of the \code{mcols} column in
#'   \code{features} that contains the feature name (e.g., \code{"name"}).
#'   Default: \code{"name"}.
#' @param window Integer. Window size in base pairs for
#'   \code{type = "proximity"}. Only features within this distance are
#'   considered. Default: \code{500L}.
#'
#' @return A \code{\link{commaData}} object identical to \code{object} except
#'   that \code{rowData} has been extended with new annotation columns:
#'   \describe{
#'     \item{For \code{type = "overlap"}:}{\code{feature_type} and
#'       \code{feature_name} columns.}
#'     \item{For \code{type = "proximity"}:}{\code{nearest_feature},
#'       \code{distance_to_feature}, and \code{rel_pos} columns.}
#'     \item{For \code{type = "metagene"}:}{\code{metagene_feature} and
#'       \code{metagene_pos} columns.}
#'   }
#'
#' @examples
#' data(comma_example_data)
#' # Overlap annotation using built-in annotation
#' annotated <- annotateSites(comma_example_data)
#' head(siteInfo(annotated)[, c("chrom", "position", "feature_type", "feature_name")])
#'
#' # Metagene annotation
#' mg <- annotateSites(comma_example_data, type = "metagene")
#' head(siteInfo(mg)[, c("position", "metagene_pos")])
#'
#' @export
annotateSites <- function(object,
                          features    = NULL,
                          type        = c("overlap", "proximity", "metagene"),
                          feature_col = "feature_type",
                          name_col    = "name",
                          window      = 500L) {
    # ── Input validation ──────────────────────────────────────────────────────
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    type <- match.arg(type)

    # Resolve feature set
    if (is.null(features)) {
        features <- annotation(object)
    }
    if (!is(features, "GRanges")) {
        stop("'features' must be a GRanges object (or NULL to use annotation(object)).")
    }
    if (length(features) == 0) {
        stop(
            "No features available for annotation. ",
            "Provide a non-empty 'features' GRanges, or supply annotation when ",
            "constructing the commaData object."
        )
    }

    # Validate required mcols columns
    if (!feature_col %in% names(GenomicRanges::mcols(features))) {
        stop(
            "Column '", feature_col, "' not found in mcols(features). ",
            "Available columns: ",
            paste(names(GenomicRanges::mcols(features)), collapse = ", ")
        )
    }
    if (!name_col %in% names(GenomicRanges::mcols(features))) {
        stop(
            "Column '", name_col, "' not found in mcols(features). ",
            "Available columns: ",
            paste(names(GenomicRanges::mcols(features)), collapse = ", ")
        )
    }

    # ── Build GRanges for sites ───────────────────────────────────────────────
    rd        <- rowData(object)
    sites_gr  <- GenomicRanges::GRanges(
        seqnames = rd$chrom,
        ranges   = IRanges::IRanges(start = rd$position, width = 1L),
        strand   = rd$strand
    )

    # ── Dispatch annotation mode ──────────────────────────────────────────────
    if (type == "overlap") {
        rd <- .annotateSites_overlap(rd, sites_gr, features, feature_col, name_col)
    } else if (type == "proximity") {
        rd <- .annotateSites_proximity(rd, sites_gr, features, feature_col, name_col, window)
    } else {
        rd <- .annotateSites_metagene(rd, sites_gr, features, feature_col, name_col)
    }

    # ── Return updated commaData ──────────────────────────────────────────────
    rowData(object) <- rd
    object
}

# ── Internal: overlap annotation ─────────────────────────────────────────────

.annotateSites_overlap <- function(rd, sites_gr, features, feature_col, name_col) {
    hits <- GenomicRanges::findOverlaps(sites_gr, features, ignore.strand = TRUE)

    # For sites with multiple overlapping features, keep the first match
    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)
    first_hit <- !duplicated(q_idx)
    q_idx <- q_idx[first_hit]
    s_idx <- s_idx[first_hit]

    feat_types  <- rep("intergenic", nrow(rd))
    feat_names  <- rep("intergenic", nrow(rd))

    feat_types[q_idx] <- as.character(GenomicRanges::mcols(features)[[feature_col]][s_idx])
    feat_names[q_idx] <- as.character(GenomicRanges::mcols(features)[[name_col]][s_idx])

    rd$feature_type <- feat_types
    rd$feature_name <- feat_names
    rd
}

# ── Internal: proximity annotation ───────────────────────────────────────────

.annotateSites_proximity <- function(rd, sites_gr, features, feature_col, name_col, window) {
    # Expand window around each site
    sites_expanded <- GenomicRanges::resize(sites_gr, width = 2L * window + 1L, fix = "center")
    # Clip to positive coordinates
    GenomicRanges::start(sites_expanded) <- pmax(GenomicRanges::start(sites_expanded), 1L)

    hits <- GenomicRanges::distanceToNearest(sites_gr, features, ignore.strand = TRUE)

    n_sites <- nrow(rd)
    nearest_feat <- rep(NA_character_, n_sites)
    dist_to_feat <- rep(NA_integer_,   n_sites)
    rel_pos      <- rep(NA_integer_,   n_sites)

    if (length(hits) > 0) {
        q_idx <- S4Vectors::queryHits(hits)
        s_idx <- S4Vectors::subjectHits(hits)
        dists <- S4Vectors::mcols(hits)$distance

        # Filter to within window
        in_window <- dists <= window
        q_idx <- q_idx[in_window]
        s_idx <- s_idx[in_window]
        dists <- dists[in_window]

        if (length(q_idx) > 0) {
            feat_starts <- GenomicRanges::start(features)[s_idx]
            feat_ends   <- GenomicRanges::end(features)[s_idx]
            feat_strand <- as.character(GenomicRanges::strand(features))[s_idx]
            site_pos    <- rd$position[q_idx]

            # TSS is start for + strand, end for - strand
            tss <- ifelse(feat_strand == "-", feat_ends, feat_starts)
            # Signed distance from TSS: positive = downstream (into gene)
            signed_dist <- ifelse(
                feat_strand == "-",
                tss - site_pos,
                site_pos - tss
            )

            nearest_feat[q_idx] <- as.character(GenomicRanges::mcols(features)[[name_col]][s_idx])
            dist_to_feat[q_idx] <- as.integer(dists)
            rel_pos[q_idx]      <- as.integer(signed_dist)
        }
    }

    rd$nearest_feature      <- nearest_feat
    rd$distance_to_feature  <- dist_to_feat
    rd$rel_pos              <- rel_pos
    rd
}

# ── Internal: metagene annotation ────────────────────────────────────────────

.annotateSites_metagene <- function(rd, sites_gr, features, feature_col, name_col) {
    hits <- GenomicRanges::findOverlaps(sites_gr, features, ignore.strand = TRUE)

    n_sites       <- nrow(rd)
    metagene_feat <- rep(NA_character_, n_sites)
    metagene_pos  <- rep(NA_real_,      n_sites)

    if (length(hits) > 0) {
        q_idx <- S4Vectors::queryHits(hits)
        s_idx <- S4Vectors::subjectHits(hits)

        # For sites in multiple features, keep the first match only
        first_hit <- !duplicated(q_idx)
        q_idx <- q_idx[first_hit]
        s_idx <- s_idx[first_hit]

        feat_starts  <- GenomicRanges::start(features)[s_idx]
        feat_ends    <- GenomicRanges::end(features)[s_idx]
        feat_widths  <- GenomicRanges::width(features)[s_idx]
        feat_strand  <- as.character(GenomicRanges::strand(features))[s_idx]
        site_pos     <- rd$position[q_idx]

        # Fractional position: 0 = feature start (TSS), 1 = feature end (TTS)
        # For - strand: 0 = high coordinate (TTS on genome), 1 = low coordinate (TSS on genome)
        frac <- ifelse(
            feat_strand == "-",
            (feat_ends - site_pos) / (feat_widths - 1L),
            (site_pos - feat_starts) / (feat_widths - 1L)
        )
        # Clamp to [0, 1] in case of rounding or single-bp features
        frac <- pmax(0, pmin(1, frac))

        metagene_feat[q_idx] <- as.character(GenomicRanges::mcols(features)[[name_col]][s_idx])
        metagene_pos[q_idx]  <- frac
    }

    rd$metagene_feature <- metagene_feat
    rd$metagene_pos     <- metagene_pos
    rd
}
