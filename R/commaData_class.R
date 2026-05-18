#' @importFrom methods new setClass setGeneric setMethod setValidity validObject is isVirtualClass
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames rowData colData rowRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqinfo seqlengths
NULL

# ─── Allowed values ──────────────────────────────────────────────────────────

.VALID_MOD_TYPES <- c("6mA", "5mC", "4mC")

# ─── Class definition ────────────────────────────────────────────────────────

#' commaData: the central data object for the comma package
#'
#' \code{commaData} is an S4 class that extends
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} to store
#' genome-wide bacterial methylation data from Oxford Nanopore sequencing.
#' It is the central object accepted and returned by all \code{comma}
#' analysis functions.
#'
#' @details
#' Genomic annotation and motif site positions are stored in
#' \code{metadata(object)} rather than as dedicated slots. Use
#' \code{\link[BiocGenerics]{annotation}(object)} and
#' \code{\link{motifSites}(object)} to access them.
#'
#' @details
#' Genome size information is stored in the \code{Seqinfo} attached to
#' \code{rowRanges(object)}, accessible via \code{seqlengths(object)} or
#' \code{seqinfo(object)}. The \code{genome()} accessor returns the same
#' named integer vector for backward compatibility.
#'
#' @details
#' The class stores methylation data in two assay matrices (accessible via
#' \code{\link[SummarizedExperiment]{assay}}):
#' \describe{
#'   \item{\code{"methylation"}}{Beta values (proportion of reads called
#'     methylated, range 0-1). Sites with coverage below the
#'     \code{min_coverage} threshold are stored as \code{NA}.}
#'   \item{\code{"coverage"}}{Integer read depth at each site.}
#' }
#'
#' Genomic positions are stored in
#' \code{\link[SummarizedExperiment]{rowRanges}(object)}, a
#' \code{\link[GenomicRanges]{GRanges}} with one 1-bp range per methylation
#' site. Per-site metadata is in the \code{mcols} of this GRanges and
#' includes at minimum: \code{mod_type} and \code{motif}. The \code{motif}
#' column stores the sequence context of each site (e.g., \code{"GATC"} or
#' \code{"CCWGG"}) as extracted from the modkit \code{mod_code} field. It is
#' \code{NA} for Dorado and Megalodon callers. The \code{mod_context} is
#' computed on demand from \code{mod_type} and \code{motif} (e.g.,
#' \code{"6mA_GATC"}, \code{"5mC_CCWGG"}), or just \code{mod_type} when
#' motif is unavailable (e.g., \code{"6mA"} for Dorado/Megalodon data). Use
#' \code{\link{modContexts}(object)} or \code{\link{siteInfo}(object)} to
#' retrieve it. All analyses default to running independently per
#' \code{mod_context} group to prevent spurious mixing of biologically distinct
#' methylation events.
#'
#' For convenience, \code{\link{siteInfo}(object)} returns a flat
#' \code{DataFrame} combining the genomic coordinates (chrom, position,
#' strand) with the mcols columns.
#'
#' Per-sample metadata is in \code{colData(object)} and includes at minimum:
#' \code{sample_name}, \code{condition}, \code{replicate}.
#'
#' @return An object of class \code{commaData}. Use
#'   \code{\link{commaData}} to construct instances.
#'
#' @seealso \code{\link{commaData}} for the constructor,
#'   \code{\link{methylation}}, \code{\link[GenomicRanges]{coverage}},
#'   \code{\link{sampleInfo}}, \code{\link{siteInfo}},
#'   \code{\link{modTypes}}, \code{\link{modContexts}},
#'   \code{\link[BiocGenerics]{annotation}} for accessors.
#'
#' @name commaData-class
#' @exportClass commaData
setClass(
    "commaData",
    contains = "RangedSummarizedExperiment"
)

# ─── Validity ────────────────────────────────────────────────────────────────

setValidity("commaData", function(object) {
    errors <- character(0)

    # ── rowRanges required mcols ────────────────────────────────────────────
    required_mcol_cols <- c("mod_type", "motif")
    rr <- rowRanges(object)
    mc <- GenomicRanges::mcols(rr)
    missing_cols <- setdiff(required_mcol_cols, colnames(mc))
    if (length(missing_cols) > 0) {
        errors <- c(errors, paste0(
            "rowRanges mcols is missing required columns: ",
            paste(missing_cols, collapse = ", ")
        ))
    }

    # ── rowRanges must be 1-bp ranges (one per methylation site) ──────────
    if (length(rr) > 0L && !is(rr, "GRangesList")) {
        widths <- GenomicRanges::width(rr)
        if (any(widths != 1L)) {
            n_bad <- sum(widths != 1L)
            errors <- c(errors, paste0(
                "rowRanges must contain 1-bp ranges (one per site), ",
                "but ", n_bad, " range(s) have width != 1. ",
                "Downstream code treats each row as a single position."
            ))
        }
    }

    # ── mod_type allowed values ─────────────────────────────────────────────
    if ("mod_type" %in% colnames(mc)) {
        bad_types <- setdiff(unique(mc$mod_type), .VALID_MOD_TYPES)
        if (length(bad_types) > 0) {
            errors <- c(errors, paste0(
                "rowRanges mcols$mod_type contains unrecognized values: ",
                paste(bad_types, collapse = ", "),
                ". Allowed values: ",
                paste(.VALID_MOD_TYPES, collapse = ", ")
            ))
        }
    }

    # ── motif column type ───────────────────────────────────────────────────
    if ("motif" %in% colnames(mc)) {
        if (!is.character(mc$motif) && !all(is.na(mc$motif))) {
            errors <- c(errors, "rowRanges mcols$motif must be a character vector (NA allowed)")
        }
    }



    # ── colData required columns ────────────────────────────────────────────
    required_col_cols <- c("sample_name", "condition", "replicate")
    cd <- colData(object)
    missing_cols2 <- setdiff(required_col_cols, colnames(cd))
    if (length(missing_cols2) > 0) {
        errors <- c(errors, paste0(
            "colData is missing required columns: ",
            paste(missing_cols2, collapse = ", ")
        ))
    }

    # ── assay names ─────────────────────────────────────────────────────────
    expected_assays <- c("methylation", "coverage")
    present_assays  <- assayNames(object)
    missing_assays  <- setdiff(expected_assays, present_assays)
    if (length(missing_assays) > 0) {
        errors <- c(errors, paste0(
            "Missing required assays: ",
            paste(missing_assays, collapse = ", ")
        ))
    }

    if (length(errors) == 0) TRUE else errors
})

# ─── show() ──────────────────────────────────────────────────────────────────

#' @importFrom methods show
setMethod("show", "commaData", function(object) {
    n_sites   <- nrow(object)
    n_samples <- ncol(object)

    cat("class: commaData\n")
    cat("sites:", n_sites, "| samples:", n_samples, "\n")

    # mod types, motifs, and contexts
    rd <- rowData(object)  # for RSE, rowData() returns mcols(rowRanges())
    if ("mod_type" %in% colnames(rd) && n_sites > 0) {
        mt <- sort(unique(rd$mod_type))
        cat("mod types:", paste(mt, collapse = ", "), "\n")
    }
    if ("motif" %in% colnames(rd) && n_sites > 0) {
        m <- sort(unique(rd$motif[!is.na(rd$motif)]))
        cat("motifs:", if (length(m) == 0L) "not available" else paste(m, collapse = ", "), "\n")
    }
    if (n_sites > 0) {
        mc <- sort(unique(modContexts(object)))
        cat("mod contexts:", paste(mc, collapse = ", "), "\n")
    }

    # conditions
    cd <- colData(object)
    if ("condition" %in% colnames(cd) && n_samples > 0) {
        cond <- sort(unique(cd$condition))
        cat("conditions:", paste(cond, collapse = ", "), "\n")
    }

    # genome info (from Seqinfo)
    sl <- GenomeInfoDb::seqlengths(object)
    if (length(sl) > 0 && !all(is.na(sl))) {
        total_bp <- sum(sl, na.rm = TRUE)
        cat("genome:", length(sl), ifelse(length(sl) == 1, "chromosome", "chromosomes"),
            paste0("(", format(total_bp, big.mark = ","), " bp total)"), "\n")
    } else {
        cat("genome: not provided\n")
    }

    # annotation / motif sites
    n_ann <- length(annotation(object))
    n_mot <- length(motifSites(object))
    cat("annotation:", if (n_ann == 0) "none" else paste(n_ann, "features"), "\n")
    cat("motif sites:", if (n_mot == 0) "none" else paste(format(n_mot, big.mark = ","), "sites"), "\n")
})
