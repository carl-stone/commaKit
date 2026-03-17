#' @importFrom GenomicRanges GRanges findOverlaps start end width strand resize mcols
#' @importFrom IRanges IRanges CharacterList IntegerList NumericList
#' @importFrom S4Vectors DataFrame queryHits subjectHits splitAsList
#' @importFrom SummarizedExperiment rowData "rowData<-"
#' @importFrom methods is
NULL

#' Annotate methylation sites relative to genomic features
#'
#' Assigns genomic feature annotations to methylation sites stored in a
#' \code{\link{commaData}} object using
#' \code{\link[GenomicRanges]{findOverlaps}}.
#' Three annotation modes are available: \code{"overlap"} assigns all
#' overlapping feature identities to each site, \code{"proximity"} reports
#' all features within a distance window and their signed offsets, and
#' \code{"metagene"} reports fractional positions within every overlapping
#' feature.
#'
#' All three modes return every matching feature per site. Results are stored as
#' \code{\link[IRanges]{CharacterList}},
#' \code{\link[IRanges]{IntegerList}}, or \code{NumericList} columns in
#' \code{rowData}.
#' Sites with no overlapping/nearby features receive length-0 list elements;
#' test for them with \code{lengths(col) == 0}.
#'
#' @param object A \code{\link{commaData}} object.
#' @param features A \code{\link[GenomicRanges]{GRanges}} of genomic features
#'   to annotate against. If \code{NULL} (default), the annotation stored in
#'   \code{object} via \code{\link[BiocGenerics]{annotation}(object)} is used.
#'   Must have mcols columns named by \code{feature_col} and \code{name_col}.
#' @param type Character string specifying the annotation mode. One of:
#'   \describe{
#'     \item{\code{"overlap"}}{(default) Each site is assigned all
#'       overlapping feature types and names. Sites that overlap no feature
#'       receive length-0 \code{CharacterList} elements.}
#'     \item{\code{"proximity"}}{Each site is assigned all features
#'       within \code{window} bp: their names, absolute distances, and signed
#'       relative positions (negative = upstream; positive = downstream of the
#'       feature TSS). Sites with no nearby features receive length-0 elements.}
#'     \item{\code{"metagene"}}{Each site that overlaps a feature is assigned
#'       a fractional position within that feature (0 = feature start, 1 =
#'       feature end) for every overlapping feature. Strand-aware: for
#'       \code{"-"} strand features, 0 is at the feature end (highest
#'       coordinate) and 1 is at the feature start (lowest coordinate).
#'       Non-overlapping sites receive length-0 elements.}
#'   }
#' @param feature_col Character string. Name of the \code{mcols} column in
#'   \code{features} that contains the feature type (e.g.,
#'   \code{"feature_type"}). Default: \code{"feature_type"}.
#' @param name_col Character string. Name of the \code{mcols} column in
#'   \code{features} that contains the feature name (e.g., \code{"name"}).
#'   Default: \code{"name"}.
#' @param window Integer. Window size in base pairs for
#'   \code{type = "proximity"}. All features within this distance are
#'   returned. Default: \code{500L}.
#'
#' @return A \code{\link{commaData}} object identical to \code{object} except
#'   that \code{rowData} has been extended with new list-valued annotation
#'   columns:
#'   \describe{
#'     \item{For \code{type = "overlap"}:}{\code{feature_types}
#'       (\code{CharacterList}) and \code{feature_names} (\code{CharacterList})
#'       — all overlapping feature types and names per site.
#'       Intergenic sites: \code{lengths(feature_types) == 0}.}
#'     \item{For \code{type = "proximity"}:}{\code{nearby_features}
#'       (\code{CharacterList}), \code{distances_to_features}
#'       (\code{IntegerList}), and \code{rel_positions} (\code{IntegerList})
#'       — all features within \code{window} bp. Sites with none:
#'       \code{lengths(nearby_features) == 0}.}
#'     \item{For \code{type = "metagene"}:}{\code{metagene_features}
#'       (\code{CharacterList}) and \code{metagene_positions}
#'       (\code{NumericList}) — all overlapping feature names and their
#'       fractional positions in \eqn{[0, 1]}. Non-overlapping sites:
#'       \code{lengths(metagene_features) == 0}.}
#'   }
#'
#' @examples
#' data(comma_example_data)
#' # Overlap annotation using built-in annotation
#' annotated <- annotateSites(comma_example_data)
#' si <- siteInfo(annotated)
#' # All overlapping feature types for the first site:
#' si$feature_types[[1]]
#' # Number of sites that overlap at least one feature:
#' sum(lengths(si$feature_types) > 0)
#' # Intergenic sites:
#' sum(lengths(si$feature_types) == 0)
#'
#' # Metagene annotation
#' mg <- annotateSites(comma_example_data, type = "metagene")
#' si_mg <- siteInfo(mg)
#' # Metagene positions for the first overlapping site:
#' si_mg$metagene_positions[[which(lengths(si_mg$metagene_positions) > 0)[1]]]
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
# Retains ALL overlapping features per site as CharacterList columns.
# This is a deliberate design choice for bacterial genomes where features
# are highly overlapping. Do NOT revert to single-match (!duplicated) behavior.

.annotateSites_overlap <- function(rd, sites_gr, features, feature_col, name_col) {
    n_sites <- nrow(rd)
    hits    <- GenomicRanges::findOverlaps(sites_gr, features, ignore.strand = TRUE)

    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)

    feat_types_raw <- as.character(GenomicRanges::mcols(features)[[feature_col]][s_idx])
    feat_names_raw <- as.character(GenomicRanges::mcols(features)[[name_col]][s_idx])

    # splitAsList groups values by site index; sites with no hit get character(0)
    site_factor <- factor(q_idx, levels = seq_len(n_sites))
    rd$feature_types <- S4Vectors::splitAsList(feat_types_raw, site_factor)
    rd$feature_names <- S4Vectors::splitAsList(feat_names_raw, site_factor)
    rd
}

# ── Internal: proximity annotation ───────────────────────────────────────────
# Returns ALL features within 'window' bp per site as IntegerList/CharacterList.
# This is a deliberate design choice. Do NOT revert to distanceToNearest()
# (single-match) behavior.

.annotateSites_proximity <- function(rd, sites_gr, features, feature_col, name_col, window) {
    n_sites <- nrow(rd)
    window  <- as.integer(window)

    # Expand a search window around each site; clip start to 1
    sites_expanded <- GenomicRanges::resize(sites_gr, width = 2L * window + 1L, fix = "center")
    GenomicRanges::start(sites_expanded) <- pmax(GenomicRanges::start(sites_expanded), 1L)

    hits  <- GenomicRanges::findOverlaps(sites_expanded, features, ignore.strand = TRUE)

    site_factor <- factor(integer(0), levels = seq_len(n_sites))  # empty default
    empty_cl    <- S4Vectors::splitAsList(character(0), site_factor)
    empty_il    <- S4Vectors::splitAsList(integer(0),   site_factor)

    if (length(hits) == 0L) {
        rd$nearby_features      <- empty_cl
        rd$distances_to_features <- empty_il
        rd$rel_positions        <- empty_il
        return(rd)
    }

    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)

    feat_starts <- GenomicRanges::start(features)[s_idx]
    feat_ends   <- GenomicRanges::end(features)[s_idx]
    feat_strand <- as.character(GenomicRanges::strand(features))[s_idx]
    feat_names  <- as.character(GenomicRanges::mcols(features)[[name_col]][s_idx])
    site_pos    <- rd$position[q_idx]

    # Absolute distance: 0 if site is inside feature, else distance to nearest edge
    dist_raw <- as.integer(pmax(0L, pmax(feat_starts - site_pos, site_pos - feat_ends)))

    # Signed position relative to TSS (start for + strand, end for - strand)
    tss         <- ifelse(feat_strand == "-", feat_ends, feat_starts)
    signed_dist <- as.integer(ifelse(feat_strand == "-", tss - site_pos, site_pos - tss))

    site_factor <- factor(q_idx, levels = seq_len(n_sites))
    rd$nearby_features       <- S4Vectors::splitAsList(feat_names,  site_factor)
    rd$distances_to_features <- S4Vectors::splitAsList(dist_raw,    site_factor)
    rd$rel_positions         <- S4Vectors::splitAsList(signed_dist, site_factor)
    rd
}

# ── Internal: metagene annotation ────────────────────────────────────────────
# Returns fractional position within ALL overlapping features per site.
# This is a deliberate design choice. Do NOT revert to single-match behavior.

.annotateSites_metagene <- function(rd, sites_gr, features, feature_col, name_col) {
    n_sites <- nrow(rd)
    hits    <- GenomicRanges::findOverlaps(sites_gr, features, ignore.strand = TRUE)

    site_factor <- factor(integer(0), levels = seq_len(n_sites))
    empty_cl    <- S4Vectors::splitAsList(character(0), site_factor)
    empty_nl    <- S4Vectors::splitAsList(numeric(0),   site_factor)

    if (length(hits) == 0L) {
        rd$metagene_features  <- empty_cl
        rd$metagene_positions <- empty_nl
        return(rd)
    }

    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)

    feat_starts <- GenomicRanges::start(features)[s_idx]
    feat_ends   <- GenomicRanges::end(features)[s_idx]
    feat_widths <- GenomicRanges::width(features)[s_idx]
    feat_strand <- as.character(GenomicRanges::strand(features))[s_idx]
    feat_names  <- as.character(GenomicRanges::mcols(features)[[name_col]][s_idx])
    site_pos    <- rd$position[q_idx]

    # Fractional position [0, 1]: 0 = TSS, 1 = TTS
    # For - strand: 0 = high coordinate (TTS on genome), 1 = low coordinate (TSS on genome)
    # For 1-bp features feat_widths - 1L == 0; use pmax(..., 1L) to avoid
    # division by zero producing NaN (which pmax/pmin do not clamp).
    denom <- pmax(feat_widths - 1L, 1L)
    frac <- ifelse(
        feat_strand == "-",
        (feat_ends - site_pos) / denom,
        (site_pos - feat_starts) / denom
    )
    # Clamp to [0, 1] in case of rounding
    frac <- pmax(0, pmin(1, frac))

    site_factor <- factor(q_idx, levels = seq_len(n_sites))
    rd$metagene_features  <- S4Vectors::splitAsList(feat_names, site_factor)
    rd$metagene_positions <- S4Vectors::splitAsList(frac,       site_factor)
    rd
}
