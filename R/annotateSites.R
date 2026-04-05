#' @importFrom GenomicRanges GRanges findOverlaps start end width strand resize mcols
#' @importFrom IRanges IRanges CharacterList IntegerList NumericList
#' @importFrom S4Vectors DataFrame queryHits subjectHits splitAsList mendoapply
#' @importFrom SummarizedExperiment rowData "rowData<-"
#' @importFrom methods is
NULL

#' Annotate methylation sites relative to genomic features
#'
#' Assigns genomic feature annotations to methylation sites stored in a
#' \code{\link{commaData}} object using
#' \code{\link[GenomicRanges]{findOverlaps}}.
#'
#' The function searches within \code{window} bp of each site and returns
#' every feature found. Four parallel list columns are added to
#' \code{rowData}: feature types, feature names, signed relative position,
#' and fractional position within the feature.
#'
#' Signed \code{rel_position} convention (strand-aware):
#' \itemize{
#'   \item \strong{0} — site is \emph{inside} the feature.
#'   \item \strong{negative} — site is upstream (before the feature in the
#'     direction of transcription).
#'   \item \strong{positive} — site is downstream (after the feature).
#' }
#'
#' \code{frac_position} is in \eqn{[0, 1]} when the site is inside a feature
#' (TSS = 0, TTS = 1, strand-aware) and \code{NA} when the site is outside.
#'
#' All list columns are strictly parallel: element \eqn{j} in
#' \code{feature_types[[i]]} corresponds to element \eqn{j} in
#' \code{feature_names[[i]]}, \code{rel_position[[i]]}, and
#' \code{frac_position[[i]]}. This invariant is preserved by any metadata
#' columns added via \code{metadata_cols}.
#'
#' @param object A \code{\link{commaData}} object.
#' @param features A \code{\link[GenomicRanges]{GRanges}} of genomic features
#'   to annotate against. If \code{NULL} (default), the annotation stored in
#'   \code{object} via \code{\link[BiocGenerics]{annotation}(object)} is used.
#'   Must have mcols columns named by \code{feature_col} and \code{name_col}.
#' @param feature_col Character string. Name of the \code{mcols} column in
#'   \code{features} that contains the feature type. Default:
#'   \code{"feature_type"}.
#' @param name_col Character string. Name of the \code{mcols} column in
#'   \code{features} that contains the feature name. Default: \code{"name"}.
#' @param window Integer. Search window in base pairs. All features within
#'   this distance of each site are returned. Default: \code{50L}.
#' @param keep Character string controlling which output columns are retained:
#'   \describe{
#'     \item{\code{"all"}}{(default) Return all four columns for all
#'       associations — including sites within the window but not inside a
#'       feature.}
#'     \item{\code{"overlap"}}{Subset each site's associations to features
#'       where \code{rel_position == 0} (inside the feature). Drop the
#'       \code{rel_position} and \code{frac_position} columns.}
#'     \item{\code{"proximity"}}{Keep all associations. Drop the
#'       \code{frac_position} column.}
#'     \item{\code{"metagene"}}{Subset to \code{rel_position == 0}. Drop the
#'       \code{rel_position} column. Retain \code{frac_position}.}
#'   }
#' @param metadata_cols Character vector or \code{NULL}. Names of additional
#'   \code{mcols} columns in \code{features} to pass through as parallel list
#'   columns in \code{rowData}. Each column \code{X} is stored as
#'   \code{X_values} (\code{CharacterList}). Default: \code{NULL}.
#'
#' @return A \code{\link{commaData}} object identical to \code{object} except
#'   that \code{rowData} has been extended with annotation list columns.
#'   With \code{keep = "all"} (default):
#'   \describe{
#'     \item{\code{feature_types}}{CharacterList. Feature type for each
#'       association per site.}
#'     \item{\code{feature_names}}{CharacterList. Feature name for each
#'       association per site.}
#'     \item{\code{rel_position}}{IntegerList. Signed relative position
#'       (bp): 0 inside, negative upstream, positive downstream.}
#'     \item{\code{frac_position}}{NumericList. Fractional position in
#'       \eqn{[0, 1]} inside features; \code{NA} outside.}
#'   }
#'   Intergenic sites (no features within \code{window}) receive length-0
#'   list elements in all columns.
#'
#' @examples
#' data(comma_example_data)
#' # Default: unified output with all four columns
#' annotated <- annotateSites(comma_example_data)
#' si <- siteInfo(annotated)
#' # Sites inside at least one feature (rel_position == 0):
#' sum(sapply(as.list(si$rel_position), function(x) any(x == 0L)))
#' # Fractional position of first inside site:
#' inside_idx <- which(lengths(si$frac_position) > 0)[1]
#' si$frac_position[[inside_idx]]
#'
#' # Backward-compatible overlap output (no rel_position / frac_position):
#' ann_ov <- annotateSites(comma_example_data, keep = "overlap")
#' siteInfo(ann_ov)$feature_names[[1]]
#'
#' @export
annotateSites <- function(object,
                          features     = NULL,
                          feature_col  = "feature_type",
                          name_col     = "name",
                          window       = 50L,
                          keep         = c("all", "overlap", "proximity", "metagene"),
                          metadata_cols = NULL) {
    # ── Input validation ──────────────────────────────────────────────────────
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    keep <- match.arg(keep)

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
    feat_mcols <- names(GenomicRanges::mcols(features))
    if (!feature_col %in% feat_mcols) {
        stop(
            "Column '", feature_col, "' not found in mcols(features). ",
            "Available columns: ", paste(feat_mcols, collapse = ", ")
        )
    }
    if (!name_col %in% feat_mcols) {
        stop(
            "Column '", name_col, "' not found in mcols(features). ",
            "Available columns: ", paste(feat_mcols, collapse = ", ")
        )
    }
    if (!is.null(metadata_cols)) {
        missing_meta <- setdiff(metadata_cols, feat_mcols)
        if (length(missing_meta) > 0L) {
            stop(
                "metadata_cols not found in mcols(features): ",
                paste(missing_meta, collapse = ", ")
            )
        }
    }

    # ── Build GRanges for sites ───────────────────────────────────────────────
    rd       <- rowData(object)
    sites_gr <- GenomicRanges::GRanges(
        seqnames = rd$chrom,
        ranges   = IRanges::IRanges(start = rd$position, width = 1L),
        strand   = rd$strand
    )

    # ── Unified annotation ────────────────────────────────────────────────────
    rd <- .annotateSites_unified(rd, sites_gr, features, feature_col, name_col,
                                 window, metadata_cols)

    # ── Post-annotation keep filter ───────────────────────────────────────────
    meta_out_cols <- if (!is.null(metadata_cols)) paste0(metadata_cols, "_values") else character(0)
    rd <- .applyKeepFilter(rd, keep, meta_out_cols)

    # ── Return updated commaData ──────────────────────────────────────────────
    rowData(object) <- rd
    object
}

# ── Internal: unified annotation backend ──────────────────────────────────────
# Computes all four parallel list columns (feature_types, feature_names,
# rel_position, frac_position) from a single findOverlaps call using the
# proximity-style expanded window. All columns are strictly parallel:
# element j in feature_types[[i]] corresponds to element j in all other
# list columns for site i. This invariant must be preserved throughout.

.annotateSites_unified <- function(rd, sites_gr, features, feature_col, name_col,
                                   window, metadata_cols) {
    n_sites <- nrow(rd)
    window  <- as.integer(window)

    # Empty list factories for the no-hit early return
    site_factor0 <- factor(integer(0), levels = seq_len(n_sites))
    empty_cl     <- S4Vectors::splitAsList(character(0), site_factor0)
    empty_il     <- S4Vectors::splitAsList(integer(0),   site_factor0)
    empty_nl     <- S4Vectors::splitAsList(numeric(0),   site_factor0)

    # Expand a search window around each site; clip start to 1
    sites_expanded <- GenomicRanges::resize(sites_gr, width = 2L * window + 1L, fix = "center")
    GenomicRanges::start(sites_expanded) <- pmax(GenomicRanges::start(sites_expanded), 1L)

    hits <- GenomicRanges::findOverlaps(sites_expanded, features, ignore.strand = TRUE)

    if (length(hits) == 0L) {
        rd$feature_types <- empty_cl
        rd$feature_names <- empty_cl
        rd$rel_position  <- empty_il
        rd$frac_position <- empty_nl
        if (!is.null(metadata_cols)) {
            for (col in metadata_cols) {
                rd[[paste0(col, "_values")]] <- empty_cl
            }
        }
        return(rd)
    }

    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)

    feat_starts <- GenomicRanges::start(features)[s_idx]
    feat_ends   <- GenomicRanges::end(features)[s_idx]
    feat_widths <- GenomicRanges::width(features)[s_idx]
    feat_strand <- as.character(GenomicRanges::strand(features))[s_idx]
    feat_types  <- as.character(GenomicRanges::mcols(features)[[feature_col]][s_idx])
    feat_names  <- as.character(GenomicRanges::mcols(features)[[name_col]][s_idx])
    site_pos    <- rd$position[q_idx]

    # ── rel_position: signed, strand-aware ────────────────────────────────────
    # 0 when site is inside the feature.
    # For + strand (and *): negative = upstream (lower coord), positive = downstream.
    # For - strand: negative = upstream (higher coord, biological), positive = downstream.
    inside     <- site_pos >= feat_starts & site_pos <= feat_ends
    pos_signed <- ifelse(inside, 0L,
                      ifelse(site_pos < feat_starts,
                             site_pos - feat_starts,   # negative: upstream
                             site_pos - feat_ends))    # positive: downstream
    neg_signed <- ifelse(inside, 0L,
                      ifelse(site_pos > feat_ends,
                             feat_ends - site_pos,     # negative: upstream on -
                             feat_starts - site_pos))  # positive: downstream on -
    rel_pos <- as.integer(ifelse(feat_strand == "-", neg_signed, pos_signed))

    # ── frac_position: [0, 1] inside features, NA outside ────────────────────
    # TSS = 0, TTS = 1 (strand-aware).
    # For 1-bp features: pmax(..., 1L) avoids division by zero.
    denom    <- pmax(feat_widths - 1L, 1L)
    frac_raw <- ifelse(
        feat_strand == "-",
        (feat_ends - site_pos) / denom,      # TSS at feat_end on - strand
        (site_pos - feat_starts) / denom     # TSS at feat_start on + strand
    )
    frac_raw <- pmax(0, pmin(1, frac_raw))   # clamp [0, 1]
    frac_pos <- ifelse(inside, frac_raw, NA_real_)

    # ── Split into parallel list columns ──────────────────────────────────────
    site_factor      <- factor(q_idx, levels = seq_len(n_sites))
    rd$feature_types <- S4Vectors::splitAsList(feat_types, site_factor)
    rd$feature_names <- S4Vectors::splitAsList(feat_names, site_factor)
    rd$rel_position  <- S4Vectors::splitAsList(rel_pos,    site_factor)
    rd$frac_position <- S4Vectors::splitAsList(frac_pos,   site_factor)

    # ── Optional additional metadata columns ──────────────────────────────────
    if (!is.null(metadata_cols)) {
        for (col in metadata_cols) {
            meta_vals <- as.character(GenomicRanges::mcols(features)[[col]][s_idx])
            rd[[paste0(col, "_values")]] <- S4Vectors::splitAsList(meta_vals, site_factor)
        }
    }

    rd
}

# ── Internal: post-annotation keep filter ─────────────────────────────────────
# Applies the 'keep' filter to rowData after unified annotation.
# Uses S4Vectors::mendoapply() to preserve typed List subclasses (CharacterList,
# NumericList) while subsetting all list columns in parallel.

.applyKeepFilter <- function(rd, keep, meta_out_cols) {
    if (keep == "all") return(rd)

    if (keep %in% c("overlap", "metagene")) {
        # Subset all list columns to indices where rel_position == 0
        keep_idx <- S4Vectors::mendoapply(function(rp) which(rp == 0L), rd$rel_position)

        rd$feature_types <- S4Vectors::mendoapply(function(x, idx) x[idx], rd$feature_types, keep_idx)
        rd$feature_names <- S4Vectors::mendoapply(function(x, idx) x[idx], rd$feature_names, keep_idx)

        if (keep == "metagene") {
            # Keep frac_position (now NA-free since all are inside)
            rd$frac_position <- S4Vectors::mendoapply(function(x, idx) x[idx], rd$frac_position, keep_idx)
        }

        # Subset metadata columns too
        for (col in meta_out_cols) {
            rd[[col]] <- S4Vectors::mendoapply(function(x, idx) x[idx], rd[[col]], keep_idx)
        }

        # Drop position columns per keep mode
        rd$rel_position <- NULL
        if (keep == "overlap") {
            rd$frac_position <- NULL
        }

    } else if (keep == "proximity") {
        # Keep all associations, just drop frac_position
        rd$frac_position <- NULL
    }

    rd
}
