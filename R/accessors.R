#' @importFrom methods setGeneric setMethod callNextMethod
#' @importFrom SummarizedExperiment assay assayNames rowData "rowData<-" colData rowRanges
#' @importFrom BiocGenerics annotation start strand
#' @importFrom IRanges coverage
#' @importFrom GenomeInfoDb genome seqnames seqlengths seqinfo
#' @importFrom GenomicRanges mcols
NULL

# ─── methylation() ───────────────────────────────────────────────────────────

#' Accessor for the methylation (beta value) matrix
#'
#' Retrieves the sites × samples matrix of methylation beta values from a
#' \code{\link{commaData}} object. Values are in the range 0–1 (proportion
#' of reads called methylated). Sites below the \code{min_coverage} threshold
#' set at object creation are \code{NA}.
#'
#' @param object A \code{commaData} object.
#'
#' @return A numeric matrix with rows corresponding to methylation sites and
#'   columns corresponding to samples. Rownames are site keys
#'   (\code{"chrom:position:strand:mod_type"}); column names are sample names.
#'
#' @seealso \code{\link{siteCoverage}}, \code{\link{siteInfo}},
#'   \code{\link{sampleInfo}}
#'
#' @examples
#' data(comma_example_data)
#' m <- methylation(comma_example_data)
#' dim(m)
#' head(m)
#'
#' @export
setGeneric("methylation", function(object) standardGeneric("methylation"))

#' @rdname methylation
setMethod("methylation", "commaData", function(object) {
    assay(object, "methylation")
})

# ─── siteCoverage() ─────────────────────────────────────────────────────────

#' Accessor for the sequencing coverage (read depth) matrix
#'
#' Retrieves the sites × samples matrix of read depth from a
#' \code{\link{commaData}} object. This package-specific accessor avoids
#' overloading \code{IRanges::coverage()}, whose conventional Bioconductor
#' meaning is genomic/Rle coverage computation.
#'
#' @param object A \code{commaData} object.
#'
#' @return An integer matrix with rows corresponding to methylation sites and
#'   columns corresponding to samples.
#'
#' @seealso \code{\link{methylation}}, \code{\link{siteInfo}}
#'
#' @examples
#' data(comma_example_data)
#' cov <- siteCoverage(comma_example_data)
#' summary(as.vector(cov))
#'
#' @export
setGeneric("siteCoverage", function(object) standardGeneric("siteCoverage"))

#' @rdname siteCoverage
setMethod("siteCoverage", "commaData", function(object) {
    assay(object, "coverage")
})

# ─── modCounts() ─────────────────────────────────────────────────────────────

#' Accessor for observed modified-read counts
#'
#' Retrieves the sites × samples matrix of observed reads called as the target
#' modification. This assay is available for callers that report count-like
#' methylation evidence directly, such as modkit pileup. Older objects or
#' probability-only callers may contain \code{NA} values.
#'
#' @param object A \code{commaData} object.
#'
#' @return An integer matrix with rows corresponding to methylation sites and
#'   columns corresponding to samples.
#'
#' @seealso \code{\link{canonicalCounts}}, \code{\link{siteCoverage}},
#'   \code{\link{methylation}}
#'
#' @examples
#' data(comma_example_data)
#' mod <- modCounts(comma_example_data)
#' dim(mod)
#'
#' @export
setGeneric("modCounts", function(object) standardGeneric("modCounts"))

#' @rdname modCounts
setMethod("modCounts", "commaData", function(object) {
    if (!"mod_counts" %in% assayNames(object)) {
        stop(
            "assay 'mod_counts' not found. This object may have been created ",
            "before raw count assays were added."
        )
    }
    assay(object, "mod_counts")
})

# ─── canonicalCounts() ───────────────────────────────────────────────────────

#' Accessor for observed canonical-read counts
#'
#' Retrieves the sites × samples matrix of observed reads called as canonical
#' or unmodified for the site. For callers that cannot provide a true
#' canonical-count decomposition, values may be \code{NA}; consult
#' \code{\link{assayProvenance}} for source details.
#'
#' @param object A \code{commaData} object.
#'
#' @return An integer matrix with rows corresponding to methylation sites and
#'   columns corresponding to samples.
#'
#' @seealso \code{\link{modCounts}}, \code{\link{siteCoverage}},
#'   \code{\link{methylation}}
#'
#' @examples
#' data(comma_example_data)
#' canonical <- canonicalCounts(comma_example_data)
#' dim(canonical)
#'
#' @export
setGeneric("canonicalCounts", function(object) standardGeneric("canonicalCounts"))

#' @rdname canonicalCounts
setMethod("canonicalCounts", "commaData", function(object) {
    if (!"canonical_counts" %in% assayNames(object)) {
        stop(
            "assay 'canonical_counts' not found. This object may have been ",
            "created before raw count assays were added."
        )
    }
    assay(object, "canonical_counts")
})

# ─── assayProvenance() ───────────────────────────────────────────────────────

#' Accessor for assay provenance metadata
#'
#' Returns metadata describing how assay matrices were produced. The constructor
#' records whether count assays were observed directly from caller output,
#' unavailable, or reconstructed for synthetic/test objects.
#'
#' @param object A \code{commaData} object.
#'
#' @return A named \code{list}. Returns an empty list for legacy objects without
#'   assay provenance metadata.
#'
#' @seealso \code{\link{modCounts}}, \code{\link{canonicalCounts}},
#'   \code{\link{methylation}}
#'
#' @examples
#' data(comma_example_data)
#' assayProvenance(comma_example_data)
#'
#' @export
setGeneric("assayProvenance", function(object) standardGeneric("assayProvenance"))

#' @rdname assayProvenance
setMethod("assayProvenance", "commaData", function(object) {
    provenance <- S4Vectors::metadata(object)$assay_provenance
    if (is.null(provenance)) list() else provenance
})

# ─── coverage() compatibility wrapper ────────────────────────────────────────

#' Deprecated coverage accessor for commaData objects
#'
#' \code{coverage(commaData)} is deprecated because \code{coverage()} is an
#' established IRanges/GenomicRanges generic for computing genomic coverage, not
#' for retrieving an assay matrix. Use \code{\link{siteCoverage}} instead.
#'
#' @param x A \code{commaData} object.
#' @param shift,width,weight,... Inherited from \code{IRanges::coverage}. These
#'   arguments are not meaningful for the commaData assay accessor and must be
#'   left at their defaults.
#'
#' @return An integer matrix with rows corresponding to methylation sites and
#'   columns corresponding to samples.
#'
#' @export
setMethod("coverage", "commaData", function(x, shift = 0L, width = NULL, weight = 1L, ...) {
    dots <- list(...)
    if (!identical(shift, 0L) || !is.null(width) || !identical(weight, 1L) ||
            length(dots) > 0L) {
        stop(
            "coverage() for commaData objects is deprecated and does not ",
            "support IRanges::coverage arguments such as 'shift', 'width', ",
            "or 'weight'. Use siteCoverage(object) to retrieve the coverage ",
            "assay matrix."
        )
    }
    .Deprecated("siteCoverage")
    siteCoverage(x)
})

# ─── sampleInfo() ────────────────────────────────────────────────────────────

#' Accessor for per-sample metadata
#'
#' Returns the per-sample metadata table from a \code{\link{commaData}} object.
#' Equivalent to \code{colData(object)} but returns a plain \code{data.frame}
#' for ease of use.
#'
#' @param object A \code{commaData} object.
#'
#' @return A \code{data.frame} with one row per sample. Always contains columns
#'   \code{sample_name}, \code{condition}, and \code{replicate}. May contain
#'   additional columns such as \code{caller} and \code{file_path}.
#'
#' @seealso \code{\link{siteInfo}}, \code{\link{modTypes}}
#'
#' @examples
#' data(comma_example_data)
#' sampleInfo(comma_example_data)
#'
#' @export
setGeneric("sampleInfo", function(object) standardGeneric("sampleInfo"))

#' @rdname sampleInfo
setMethod("sampleInfo", "commaData", function(object) {
    as.data.frame(colData(object))
})

# ─── siteInfo() ──────────────────────────────────────────────────────────────

#' Accessor for per-site metadata
#'
#' Returns the per-site metadata table from a \code{\link{commaData}} object.
#' Reconstructs a flat \code{DataFrame} from the \code{rowRanges()} GRanges,
#' combining genomic coordinates (chrom, position, strand) with the mcols
#' columns (mod_type, motif, mod_context, plus any annotation/result columns).
#' This provides a backward-compatible interface to the pre-Schema-v2
#' \code{rowData()} layout.
#'
#' @param object A \code{commaData} object.
#'
#' @return A \code{\link[S4Vectors]{DataFrame}} with one row per methylation site.
#'   Always contains columns \code{chrom}, \code{position}, \code{strand},
#'   \code{mod_type}, \code{motif} (the sequence context; \code{NA} for
#'   Dorado/Megalodon callers), \code{mod_context} (the composite
#'   modification context, e.g., \code{"6mA_GATC"}), and \code{site_key}
#'   (a human-readable label with fixed
#'   \code{"chrom:position:strand:mod_type:motif"} fields, e.g.,
#'   \code{"chr1:512:+:6mA:GATC"}; computed on demand, not used for
#'   internal matching). May contain additional
#'   annotation columns added by \code{\link[=annotateSites]{annotateSites()}}
#'   or result columns from \code{\link{diffMethyl}()}.
#'
#' @seealso \code{\link{methylation}}, \code{\link{modTypes}}
#'
#' @examples
#' data(comma_example_data)
#' head(siteInfo(comma_example_data))
#'
#' @export
setGeneric("siteInfo", function(object) standardGeneric("siteInfo"))

#' @rdname siteInfo
setMethod("siteInfo", "commaData", function(object) {
    rr <- rowRanges(object)
    mc <- GenomicRanges::mcols(rr)
    df <- S4Vectors::DataFrame(
        chrom       = as.character(GenomeInfoDb::seqnames(rr)),
        position    = BiocGenerics::start(rr),
        strand      = as.character(BiocGenerics::strand(rr)),
        mc,
        row.names   = NULL
    )
    # Add computed mod_context column if not already present
    if (!"mod_context" %in% colnames(df) &&
        "mod_type" %in% colnames(mc) && "motif" %in% colnames(mc)) {
        df$mod_context <- .computeModContext(mc$mod_type, mc$motif)
    }
    # Add computed site_key column for human readability (not used for matching).
    if ("mod_type" %in% colnames(mc) && "motif" %in% colnames(mc)) {
        df$site_key <- paste(df$chrom, df$position, df$strand,
                             as.character(df$mod_type), df$motif, sep = ":")
    }
    df
})

# ─── modTypes() ──────────────────────────────────────────────────────────────

#' Return the modification types present in a commaData object
#'
#' Returns the unique methylation modification types stored in a
#' \code{\link{commaData}} object.
#'
#' @param object A \code{commaData} object.
#'
#' @return A character vector of modification types present in
#'   \code{rowData(object)$mod_type} (e.g., \code{c("6mA", "5mC")}).
#'
#' @examples
#' data(comma_example_data)
#' modTypes(comma_example_data)
#'
#' @export
setGeneric("modTypes", function(object) standardGeneric("modTypes"))

#' @rdname modTypes
setMethod("modTypes", "commaData", function(object) {
    sort(unique(as.character(rowData(object)$mod_type)))
})

# ─── motifs() ────────────────────────────────────────────────────────────────

#' Accessor for sequence context motifs present in a commaData object
#'
#' Returns the sorted unique motif strings stored in
#' \code{rowData(object)$motif}. \code{NA} values (sites from Dorado or
#' Megalodon callers where motif context is unavailable) are excluded from
#' the result.
#'
#' @param object A \code{commaData} object.
#'
#' @return A sorted character vector of unique non-\code{NA} motif strings
#'   (e.g., \code{c("CCWGG", "GATC")}). Returns \code{character(0)} if all
#'   motif values are \code{NA} (e.g., Dorado-only data).
#'
#' @examples
#' data(comma_example_data)
#' motifs(comma_example_data)
#'
#' @export
setGeneric("motifs", function(object) standardGeneric("motifs"))

#' @rdname motifs
setMethod("motifs", "commaData", function(object) {
    all_m <- rowData(object)$motif
    sort(unique(all_m[!is.na(all_m)]))
})

# ─── modContexts() ───────────────────────────────────────────────────────────

#' Return the modification contexts present in a commaData object
#'
#' Returns the unique modification contexts stored in a
#' \code{\link{commaData}} object. A \code{mod_context} is a composite string
#' combining modification type and sequence motif:
#' \code{paste(mod_type, motif, sep = "_")} when motif information is available
#' (e.g., \code{"6mA_GATC"}, \code{"5mC_CCWGG"}), or just \code{mod_type} for
#' callers that do not provide per-site motif context (e.g., \code{"6mA"} for
#' Dorado or Megalodon data).
#'
#' All differential methylation analyses run independently per
#' \code{mod_context} group by default, preventing spurious pooling of
#' biologically distinct methylation events (e.g., 6mA at GATC motifs from
#' Dam methyltransferase versus any cytosine methylation detected at GATC
#' positions, which is likely artefactual).
#'
#' @param object A \code{commaData} object.
#'
#' @return A sorted character vector of unique \code{mod_context} strings
#'   present in \code{rowData(object)$mod_context}
#'   (e.g., \code{c("5mC_CCWGG", "6mA_GATC")}).
#'
#' @seealso \code{\link{modTypes}}, \code{\link{motifs}}, \code{\link{subset}}
#'
#' @examples
#' data(comma_example_data)
#' modContexts(comma_example_data)
#'
#' @export
setGeneric("modContexts", function(object) standardGeneric("modContexts"))

#' @rdname modContexts
setMethod("modContexts", "commaData", function(object) {
    mc <- GenomicRanges::mcols(rowRanges(object))
    sort(unique(.computeModContext(mc$mod_type, mc$motif)))
})

# ─── genome() ────────────────────────────────────────────────────────────────

#' Accessor for genome size information
#'
#' Returns the chromosome sizes stored in a \code{\link{commaData}} object.
#' Genome size information is stored in the \code{Seqinfo} attached to
#' \code{rowRanges(object)}. This accessor returns \code{seqlengths(object)}
#' for backward compatibility.
#'
#' @param x A \code{commaData} object.
#'
#' @return A named integer vector of chromosome sizes
#'   (chromosome name -> length in bp), or \code{NULL} if no genome information
#'   was provided at construction.
#'
#' @examples
#' data(comma_example_data)
#' genome(comma_example_data)
#'
#' @export
setMethod("genome", "commaData", function(x) {
    sl <- GenomeInfoDb::seqlengths(x)
    if (length(sl) == 0 || all(is.na(sl))) NULL else sl
})

# ─── annotation() ────────────────────────────────────────────────────────────

#' Accessor for genomic feature annotation
#'
#' Returns the \code{\link[GenomicRanges]{GRanges}} of genomic features stored
#' in a \code{\link{commaData}} object. This is the annotation loaded from a
#' GFF3 or BED file at construction time.
#'
#' @param object A \code{commaData} object.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object. May be empty (length
#'   0) if no annotation was provided when creating the object.
#'
#' @examples
#' data(comma_example_data)
#' annotation(comma_example_data)
#'
#' @export
setMethod("annotation", "commaData", function(object) {
    md <- S4Vectors::metadata(object)
    if (is.null(md$annotation)) GenomicRanges::GRanges() else md$annotation
})

# ─── motifSites() ────────────────────────────────────────────────────────────

#' Accessor for motif site positions
#'
#' Returns the \code{\link[GenomicRanges]{GRanges}} of all instances of the
#' user-specified sequence motif in the genome, as computed by
#' \code{\link{findMotifSites}} during object construction.
#'
#' @param object A \code{commaData} object.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object. May be empty (length
#'   0) if no motif was specified at construction.
#'
#' @examples
#' data(comma_example_data)
#' motifSites(comma_example_data)
#'
#' @export
setGeneric("motifSites", function(object) standardGeneric("motifSites"))

#' @rdname motifSites
setMethod("motifSites", "commaData", function(object) {
    md <- S4Vectors::metadata(object)
    if (is.null(md$motifSites)) GenomicRanges::GRanges() else md$motifSites
})

# ─── [ subsetting ────────────────────────────────────────────────────────────

#' Subset a commaData object by sites and/or samples
#'
#' Standard bracket-based subsetting. Rows correspond to methylation sites;
#' columns correspond to samples. The resulting object is a valid
#' \code{commaData} with all assays, rowData, colData, and custom slots
#' updated consistently.
#'
#' @param x A \code{commaData} object.
#' @param i Row (site) index: integer, logical, or character vector.
#' @param j Column (sample) index: integer, logical, or character vector.
#' @param ... Not used.
#' @param drop Ignored (required by generic).
#'
#' @return A \code{commaData} object with the selected sites and samples.
#'
#' @examples
#' data(comma_example_data)
#' # First 50 sites, all samples
#' sub <- comma_example_data[1:50, ]
#' dim(sub)
#'
#' @export
setMethod("[", "commaData", function(x, i, j, ..., drop = FALSE) {
    # Delegate to SummarizedExperiment's [ method, then re-wrap
    se_sub <- callNextMethod()

    # metadata is automatically preserved by RSE subsetting
    new("commaData", se_sub)
})

# ─── filterSites() ───────────────────────────────────────────────────────────

#' Filter a commaData object by condition, modification type, or chromosome
#'
#' A convenience function for filtering a \code{\link{commaData}} object by
#' common criteria. For arbitrary index-based subsetting, use \code{[}. This
#' package-specific name avoids exporting a broad \code{subset()} generic that
#' masks \code{base::subset()}.
#'
#' @param x A \code{commaData} object.
#' @param mod_type Character vector or \code{NULL}. If provided, only sites
#'   with a matching modification type are kept (e.g., \code{"6mA"}).
#' @param condition Character vector or \code{NULL}. If provided, only samples
#'   matching the specified condition(s) are kept.
#' @param chrom Character vector or \code{NULL}. If provided, only sites on
#'   the specified chromosome(s) are kept.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   a matching sequence context motif are kept (e.g., \code{"GATC"}). Sites
#'   with \code{NA} motif values are excluded when this filter is active.
#'   Use \code{\link{motifs}} to see which motifs are present.
#' @param mod_context Character vector or \code{NULL}. If provided, only sites
#'   with a matching modification context are kept (e.g.,
#'   \code{"6mA_GATC"}, \code{"5mC_CCWGG"}). A \code{mod_context} value is
#'   \code{paste(mod_type, motif, sep = "_")} when motif is available, or just
#'   \code{mod_type} for Dorado/Megalodon data. Use \code{\link{modContexts}}
#'   to see which contexts are present. When provided, this filter is applied in
#'   addition to (ANDed with) any \code{mod_type} or \code{motif} filters.
#' @param ... Ignored.
#'
#' @return A \code{commaData} object containing only the selected sites and
#'   samples.
#'
#' @examples
#' data(comma_example_data)
#' # Only 6mA sites
#' six_ma <- filterSites(comma_example_data, mod_type = "6mA")
#' modTypes(six_ma)
#'
#' # Only GATC-context sites
#' gatc <- filterSites(comma_example_data, motif = "GATC")
#' nrow(gatc)
#'
#' # Filter by mod_context (equivalent to the above for modkit data)
#' gatc2 <- filterSites(comma_example_data, mod_context = "6mA_GATC")
#' nrow(gatc2)
#'
#' @export
filterSites <- function(x, mod_type = NULL, condition = NULL, chrom = NULL,
                        motif = NULL, mod_context = NULL, ...) {
    if (!is(x, "commaData")) {
        stop("'x' must be a commaData object.")
    }

    rr <- rowRanges(x)
    mc <- GenomicRanges::mcols(rr)
    cd <- colData(x)

    # Site filter
    site_keep <- rep(TRUE, nrow(x))
    if (!is.null(mod_type)) {
        site_keep <- site_keep & (mc$mod_type %in% mod_type)
    }
    if (!is.null(chrom)) {
        site_keep <- site_keep & (as.character(GenomeInfoDb::seqnames(rr)) %in% chrom)
    }
    if (!is.null(motif)) {
        site_keep <- site_keep & (!is.na(mc$motif)) & (mc$motif %in% motif)
    }
    if (!is.null(mod_context)) {
        # Compute mod_context on demand for filtering
        computed_ctx <- .computeModContext(mc$mod_type, mc$motif)
        site_keep <- site_keep & (computed_ctx %in% mod_context)
    }

    # Sample filter
    samp_keep <- rep(TRUE, ncol(x))
    if (!is.null(condition)) {
        samp_keep <- samp_keep & (cd$condition %in% condition)
    }

    x[site_keep, samp_keep]
}

#' Deprecated subset method for commaData objects
#'
#' \code{subset(commaData)} is deprecated to avoid package-level masking of
#' \code{base::subset()}. Use \code{\link{filterSites}} for common
#' commaData filters or \code{[} for index-based subsetting.
#'
#' @inheritParams filterSites
#' @return A \code{commaData} object containing only the selected sites and
#'   samples.
#' @export
subset.commaData <- function(x, mod_type = NULL, condition = NULL, chrom = NULL,
                             motif = NULL, mod_context = NULL, ...) {
    warning(
        "subset.commaData() is deprecated; use filterSites() instead.",
        call. = FALSE
    )
    filterSites(x, mod_type = mod_type, condition = condition, chrom = chrom,
                motif = motif, mod_context = mod_context, ...)
}

# SummarizedExperiment/S4Vectors promote subset() to an S4 generic when
# attached, so register the same compatibility path for S4 dispatch.
setMethod("subset", "commaData", function(x, ...) {
    subset.commaData(x, ...)
})

# ─── caller() ────────────────────────────────────────────────────────────────

#' Accessor for the methylation caller
#'
#' Returns the name of the methylation caller that produced the data
#' (e.g., \code{"modkit"}, \code{"megalodon"}, or \code{"dorado"}).
#' The caller is stored in \code{metadata(object)} at construction time.
#'
#' @param object A \code{commaData} object.
#'
#' @return A character string naming the caller, or \code{NA} if not stored
#'   (e.g., objects created before caller storage was implemented).
#'
#' @examples
#' data(comma_example_data)
#' caller(comma_example_data)
#'
#' @export
setGeneric("caller", function(object) standardGeneric("caller"))

#' @rdname caller
setMethod("caller", "commaData", function(object) {
    md <- S4Vectors::metadata(object)
    if (is.null(md$caller)) NA_character_ else md$caller
})

# ─── minCoverage() ───────────────────────────────────────────────────────────

#' Accessor for the minimum coverage threshold
#'
#' Returns the minimum read depth threshold that was applied at construction
#' time. Sites with coverage below this threshold have their beta value set
#' to \code{NA}.
#'
#' @param object A \code{commaData} object.
#'
#' @return An integer (the minimum coverage threshold), or \code{NA_integer_}
#'   if not stored (e.g., objects created before min_coverage storage was
#'   implemented).
#'
#' @examples
#' data(comma_example_data)
#' minCoverage(comma_example_data)
#'
#' @export
setGeneric("minCoverage", function(object) standardGeneric("minCoverage"))

#' @rdname minCoverage
setMethod("minCoverage", "commaData", function(object) {
    md <- S4Vectors::metadata(object)
    if (is.null(md$min_coverage)) NA_integer_ else md$min_coverage
})

# ─── .computeModContext() ──────────────────────────────────────────────────

# Internal helper: compute mod_context from mod_type and motif vectors.
# Returns "mod_type_motif" when motif is known, or just "mod_type" when
# motif is NA (e.g., Dorado/Megalodon callers).
.computeModContext <- function(mod_type, motif) {
    # Convert factor to character (paste/ifelse on factor return integer codes)
    mod_type <- as.character(mod_type)
    ifelse(is.na(motif), mod_type, paste(mod_type, motif, sep = "_"))
}

# ─── .validateModType() ───────────────────────────────────────────────────

# Internal helper: validate that all requested mod_type values exist in the
# object. Called by 13 functions that accept a mod_type filter parameter.
# Stops with an informative error listing the invalid values and available
# types. Returns invisibly if all values are valid.
.validateModType <- function(mod_type, object) {
    available <- modTypes(object)
    bad <- setdiff(mod_type, available)
    if (length(bad) > 0L) {
        stop("'mod_type' value(s) not found in object: ",
             paste(bad, collapse = ", "),
             ". Available types: ", paste(available, collapse = ", "), ".")
    }
    invisible(NULL)
}

# Centralized validator for site-level filters used by exported helpers.
.validateSiteFilterValues <- function(arg, values, available, object_label = "object") {
    bad <- setdiff(values, available)
    if (length(bad) > 0L) {
        available_label <- if (length(available) > 0L) {
            paste(available, collapse = ", ")
        } else {
            "<none>"
        }
        stop(
            "'", arg, "' value(s) not found in ", object_label, ": ",
            paste(bad, collapse = ", "),
            ". Available: ", available_label, "."
        )
    }
    invisible(NULL)
}

.siteFilterLabel <- function(arg, values) {
    paste0(arg, " = '", paste(values, collapse = "', '"), "'")
}

.stopEmptySiteFilter <- function(filters, caller = NULL) {
    where <- if (is.null(caller)) "" else paste0(" in ", caller)
    stop(
        "No sites remain after applying site filters", where, ": ",
        paste(filters, collapse = "; "), "."
    )
}

# Apply mod_type, motif, and mod_context filters with consistent validation and
# empty-result handling. Filters are applied sequentially so motif/mod_context
# values are validated against the sites that remain after earlier filters.
.applySiteFilters <- function(object, mod_type = NULL, motif = NULL,
                              mod_context = NULL, stop_on_empty = TRUE,
                              caller = NULL) {
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }

    filters <- character(0)
    if (!is.null(mod_type)) {
        .validateModType(mod_type, object)
        object <- filterSites(object, mod_type = mod_type)
        filters <- c(filters, .siteFilterLabel("mod_type", mod_type))
        if (stop_on_empty && nrow(object) == 0L) {
            .stopEmptySiteFilter(filters, caller = caller)
        }
    }

    if (!is.null(motif)) {
        .validateSiteFilterValues("motif", motif, motifs(object))
        object <- filterSites(object, motif = motif)
        filters <- c(filters, .siteFilterLabel("motif", motif))
        if (stop_on_empty && nrow(object) == 0L) {
            .stopEmptySiteFilter(filters, caller = caller)
        }
    }

    if (!is.null(mod_context)) {
        .validateSiteFilterValues("mod_context", mod_context, modContexts(object))
        object <- filterSites(object, mod_context = mod_context)
        filters <- c(filters, .siteFilterLabel("mod_context", mod_context))
        if (stop_on_empty && nrow(object) == 0L) {
            .stopEmptySiteFilter(filters, caller = caller)
        }
    }

    object
}

# ─── .checkModTypeValues() ────────────────────────────────────────────────

# Internal helper: check that values and/or factor levels are valid mod_types.
# Returns a character vector of error messages (empty if valid). Used by the
# commaData validity method and by .validateModType(). This is the single
# source of truth for what constitutes a valid mod_type value.
.checkModTypeValues <- function(values, levels = NULL) {
    errors <- character(0)
    # Check factor levels if provided
    if (!is.null(levels)) {
        bad_levels <- setdiff(levels, .VALID_MOD_TYPES)
        if (length(bad_levels) > 0L) {
            errors <- c(errors, paste0(
                "rowRanges mcols$mod_type has invalid factor levels: ",
                paste(bad_levels, collapse = ", "),
                ". Allowed levels: ",
                paste(.VALID_MOD_TYPES, collapse = ", ")
            ))
        }
    }
    # Check actual values
    bad_vals <- setdiff(values, .VALID_MOD_TYPES)
    if (length(bad_vals) > 0L) {
        errors <- c(errors, paste0(
            "rowRanges mcols$mod_type contains unrecognized values: ",
            paste(bad_vals, collapse = ", "),
            ". Allowed values: ",
            paste(.VALID_MOD_TYPES, collapse = ", ")
        ))
    }
    errors
}
