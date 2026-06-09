#' @importFrom methods new setClass setGeneric setMethod setValidity validObject is isVirtualClass
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames rowData colData rowRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqinfo seqlengths
NULL

# ─── Allowed values ──────────────────────────────────────────────────────────

.VALID_MOD_TYPES <- c("6mA", "5mC", "4mC")

# ─── Class definition ────────────────────────────────────────────────────────

#' commaData: the central data object for the commaKit package
#'
#' \code{commaData} is an S4 class that extends
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} to store
#' genome-wide bacterial methylation data from Oxford Nanopore sequencing.
#' It is the central object accepted and returned by all commaKit
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
#' \code{seqinfo(object)}. Use \code{\link{genomeSizes}} for chromosome sizes; \code{genome()}
#' is retained only as a backward-compatible size-vector method.
#'
#' @details
#' The class stores methylation data in assay matrices (accessible via
#' \code{\link[SummarizedExperiment]{assay}}):
#' \describe{
#'   \item{\code{"methylation"}}{Beta values (proportion of reads called
#'     methylated, range 0-1). Sites with coverage below the
#'     \code{min_coverage} threshold are stored as \code{NA}.}
#'   \item{\code{"coverage"}}{Integer read depth at each site.}
#'   \item{\code{"mod_counts"}}{Observed reads called as the target
#'     modification, when available from the caller.}
#'   \item{\code{"canonical_counts"}}{Observed reads called canonical or
#'     unmodified, when available from the caller.}
#' }
#'
#' Genomic positions are stored in
#' \code{\link[SummarizedExperiment]{rowRanges}(object)}, a
#' \code{\link[GenomicRanges]{GRanges}} with one 1-bp range per methylation
#' site. Per-site metadata is in the \code{mcols} of this GRanges and
#' includes at minimum: \code{mod_type} and \code{motif}. The \code{mod_type}
#' column is a factor with levels \code{c("4mC", "5mC", "6mA")}, enforcing
#' valid values at the data structure level. The \code{motif}
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
#' The methylation caller and minimum coverage threshold are stored in
#' \code{metadata(object)} and accessible via \code{\link{caller}(object)}
#' and \code{\link{minCoverage}(object)}.
#'
#' Assay-layer provenance and defaults are stored in
#' \code{metadata(object)$assay_provenance} and
#' \code{metadata(object)$assay_defaults}. Use \code{\link{assayLayers}} for a
#' tabular summary.
#'
#' Differential methylation result layers are stored in
#' \code{metadata(object)$diffMethyl_results} with provenance in
#' \code{metadata(object)$diffMethyl_result_layers}. Use
#' \code{\link{resultLayers}} to list named result runs.
#'
#' @return An object of class \code{commaData}. Use
#'   \code{\link{commaData}} to construct instances.
#'
#' @seealso \code{\link{commaData}} for the constructor,
#'   \code{\link{methylation}}, \code{\link{siteCoverage}},
#'   \code{\link{modCounts}}, \code{\link{canonicalCounts}},
#'   \code{\link{assayLayers}}, \code{\link{assayProvenance}},
#'   \code{\link{resultLayers}},
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

    # ── mod_type must be a factor with valid levels ─────────────────────────
    if ("mod_type" %in% colnames(mc)) {
        if (!is.factor(mc$mod_type)) {
            errors <- c(errors, paste0(
                "rowRanges mcols$mod_type must be a factor. ",
                "Use factor(mod_type, levels = .VALID_MOD_TYPES) when constructing."
            ))
        } else {
            # Check levels and values using shared helper
            mt_errors <- .checkModTypeValues(
                values = as.character(unique(mc$mod_type)),
                levels = levels(mc$mod_type)
            )
            errors <- c(errors, mt_errors)
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


    # ── assay value invariants ──────────────────────────────────────────────
    if ("methylation" %in% assayNames(object)) {
        methyl <- SummarizedExperiment::assay(object, "methylation")
        methyl_vals <- methyl[!is.na(methyl)]
        if (length(methyl_vals) > 0L &&
                any(!is.finite(methyl_vals) | methyl_vals < 0 | methyl_vals > 1)) {
            errors <- c(errors,
                "assay 'methylation' must contain beta values in [0, 1] or NA"
            )
        }
    }
    if ("coverage" %in% assayNames(object)) {
        cov <- SummarizedExperiment::assay(object, "coverage")
        cov_vals <- cov[!is.na(cov)]
        if (length(cov_vals) > 0L) {
            bad_cov <- !is.finite(cov_vals) | cov_vals < 0 |
                       abs(cov_vals - round(cov_vals)) > sqrt(.Machine$double.eps)
            if (any(bad_cov)) {
                errors <- c(errors,
                    "assay 'coverage' must contain non-negative integer-like values or NA"
                )
            }
        }
    }
    for (assay_name in c("mod_counts", "canonical_counts")) {
        if (assay_name %in% assayNames(object)) {
            counts <- SummarizedExperiment::assay(object, assay_name)
            count_vals <- counts[!is.na(counts)]
            if (length(count_vals) > 0L) {
                bad_counts <- !is.finite(count_vals) | count_vals < 0 |
                    abs(count_vals - round(count_vals)) > sqrt(.Machine$double.eps)
                if (any(bad_counts)) {
                    errors <- c(errors, paste0(
                        "assay '", assay_name,
                        "' must contain non-negative integer-like values or NA"
                    ))
                }
            }
        }
    }
    if (all(c("mod_counts", "canonical_counts", "coverage") %in% assayNames(object))) {
        mod_counts <- SummarizedExperiment::assay(object, "mod_counts")
        canonical_counts <- SummarizedExperiment::assay(object, "canonical_counts")
        cov <- SummarizedExperiment::assay(object, "coverage")
        comparable <- !is.na(mod_counts) & !is.na(canonical_counts) & !is.na(cov)
        too_large <- comparable & (mod_counts + canonical_counts > cov)
        if (any(too_large)) {
            errors <- c(errors,
                "assays 'mod_counts' + 'canonical_counts' must be <= assay 'coverage'"
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
    cat("sites: ", n_sites, " | samples: ", n_samples, "\n", sep = "")

    # mod types, motifs, and contexts
    rd <- rowData(object)  # for RSE, rowData() returns mcols(rowRanges())
    if ("mod_type" %in% colnames(rd) && n_sites > 0) {
        mt <- sort(unique(rd$mod_type))
        cat("mod types: ", paste(mt, collapse = ", "), "\n", sep = "")
    }
    if ("motif" %in% colnames(rd) && n_sites > 0) {
        m <- sort(unique(rd$motif[!is.na(rd$motif)]))
        cat("motifs: ",
            if (length(m) == 0L) "not available" else paste(m, collapse = ", "),
            "\n",
            sep = "")
    }
    if (n_sites > 0) {
        mc <- sort(unique(modContexts(object)))
        cat("mod contexts: ", paste(mc, collapse = ", "), "\n", sep = "")
    }

    # conditions
    cd <- colData(object)
    if ("condition" %in% colnames(cd) && n_samples > 0) {
        cond <- sort(unique(cd$condition))
        cat("conditions: ", paste(cond, collapse = ", "), "\n", sep = "")
    }

    # genome info (from Seqinfo)
    sl <- GenomeInfoDb::seqlengths(object)
    if (length(sl) > 0 && !all(is.na(sl))) {
        total_bp <- sum(sl, na.rm = TRUE)
        cat("genome: ",
            length(sl), " ",
            ifelse(length(sl) == 1, "chromosome", "chromosomes"), " ",
            paste0("(", format(total_bp, big.mark = ","), " bp total)"),
            "\n",
            sep = "")
    } else {
        cat("genome: not provided\n")
    }

    # annotation / motif sites
    n_ann <- length(annotation(object))
    n_mot <- length(motifSites(object))
    cat("annotation: ", if (n_ann == 0) "none" else paste(n_ann, "features"), "\n", sep = "")
    cat("motif sites: ",
        if (n_mot == 0) "none" else paste(format(n_mot, big.mark = ","), "sites"),
        "\n",
        sep = "")

    # caller and min_coverage
    cl <- S4Vectors::metadata(object)$caller
    mc <- S4Vectors::metadata(object)$min_coverage
    if (!is.null(cl)) cat("caller: ", cl, "\n", sep = "")
    if (!is.null(mc)) cat("min_coverage: ", mc, "\n", sep = "")
})
