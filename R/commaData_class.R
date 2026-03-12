#' @importFrom methods new setClass setGeneric setMethod setValidity validObject
#'   is isVirtualClass
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#'   rowData colData
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
NULL

# ─── Allowed values ──────────────────────────────────────────────────────────

.VALID_MOD_TYPES <- c("6mA", "5mC", "4mC")

# ─── Class definition ────────────────────────────────────────────────────────

#' commaData: the central data object for the comma package
#'
#' \code{commaData} is an S4 class that extends
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} to store
#' genome-wide bacterial methylation data from Oxford Nanopore sequencing.
#' It is the central object accepted and returned by all \code{comma}
#' analysis functions.
#'
#' @slot genomeInfo Named integer vector of chromosome sizes
#'   (chromosome name → length in bp).
#' @slot annotation \code{\link[GenomicRanges]{GRanges}} of genomic features
#'   loaded from a GFF3 or BED file. May be an empty \code{GRanges} if no
#'   annotation was provided.
#' @slot motifSites \code{\link[GenomicRanges]{GRanges}} of all instances of
#'   the user-specified sequence motif in the genome (e.g., all GATC sites).
#'   May be an empty \code{GRanges} if no motif was specified.
#'
#' @details
#' The class stores methylation data in two assay matrices (accessible via
#' \code{\link[SummarizedExperiment]{assay}}):
#' \describe{
#'   \item{\code{"methylation"}}{Beta values (proportion of reads called
#'     methylated, range 0–1). Sites with coverage below the
#'     \code{min_coverage} threshold are stored as \code{NA}.}
#'   \item{\code{"coverage"}}{Integer read depth at each site.}
#' }
#'
#' Per-site metadata is in \code{rowData(object)} and includes at minimum:
#' \code{chrom}, \code{position}, \code{strand}, \code{mod_type}.
#'
#' Per-sample metadata is in \code{colData(object)} and includes at minimum:
#' \code{sample_name}, \code{condition}, \code{replicate}.
#'
#' @seealso \code{\link{commaData}} for the constructor,
#'   \code{\link{methylation}}, \code{\link{coverage}},
#'   \code{\link{sampleInfo}}, \code{\link{siteInfo}},
#'   \code{\link{modTypes}}, \code{\link{genome}},
#'   \code{\link{annotation}} for accessors.
#'
#' @name commaData-class
#' @exportClass commaData
setClass(
    "commaData",
    contains = "SummarizedExperiment",
    representation(
        genomeInfo  = "ANY",   # named integer vector or NULL
        annotation  = "GRanges",
        motifSites  = "GRanges"
    ),
    prototype(
        genomeInfo  = NULL,
        annotation  = GenomicRanges::GRanges(),
        motifSites  = GenomicRanges::GRanges()
    )
)

# ─── Validity ────────────────────────────────────────────────────────────────

setValidity("commaData", function(object) {
    errors <- character(0)

    # ── rowData required columns ────────────────────────────────────────────
    required_row_cols <- c("chrom", "position", "strand", "mod_type")
    rd <- rowData(object)
    missing_cols <- setdiff(required_row_cols, colnames(rd))
    if (length(missing_cols) > 0) {
        errors <- c(errors, paste0(
            "rowData is missing required columns: ",
            paste(missing_cols, collapse = ", ")
        ))
    }

    # ── mod_type allowed values ─────────────────────────────────────────────
    if ("mod_type" %in% colnames(rd)) {
        bad_types <- setdiff(unique(rd$mod_type), .VALID_MOD_TYPES)
        if (length(bad_types) > 0) {
            errors <- c(errors, paste0(
                "rowData$mod_type contains unrecognized values: ",
                paste(bad_types, collapse = ", "),
                ". Allowed values: ",
                paste(.VALID_MOD_TYPES, collapse = ", ")
            ))
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

    # ── genomeInfo type ─────────────────────────────────────────────────────
    gi <- object@genomeInfo
    if (!is.null(gi)) {
        if (!is.integer(gi) || is.null(names(gi))) {
            errors <- c(errors,
                "genomeInfo must be a named integer vector or NULL"
            )
        }
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

    # mod types
    rd <- rowData(object)
    if ("mod_type" %in% colnames(rd) && n_sites > 0) {
        mt <- sort(unique(rd$mod_type))
        cat("mod types:", paste(mt, collapse = ", "), "\n")
    }

    # conditions
    cd <- colData(object)
    if ("condition" %in% colnames(cd) && n_samples > 0) {
        cond <- sort(unique(cd$condition))
        cat("conditions:", paste(cond, collapse = ", "), "\n")
    }

    # genome info
    gi <- object@genomeInfo
    if (!is.null(gi) && length(gi) > 0) {
        total_bp <- sum(gi)
        cat("genome:", length(gi), ifelse(length(gi) == 1, "chromosome", "chromosomes"),
            paste0("(", format(total_bp, big.mark = ","), " bp total)"), "\n")
    } else {
        cat("genome: not provided\n")
    }

    # annotation / motif sites
    n_ann <- length(object@annotation)
    n_mot <- length(object@motifSites)
    cat("annotation:", if (n_ann == 0) "none" else paste(n_ann, "features"), "\n")
    cat("motif sites:", if (n_mot == 0) "none" else paste(format(n_mot, big.mark = ","), "sites"), "\n")
})
