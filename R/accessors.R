#' @importFrom methods setGeneric setMethod callNextMethod
#' @importFrom SummarizedExperiment assay rowData "rowData<-" colData rowRanges
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
#' @seealso \code{\link[=coverage,commaData-method]{coverage}}, \code{\link{siteInfo}},
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

# ─── coverage() ──────────────────────────────────────────────────────────────

#' Accessor for the sequencing coverage (read depth) matrix
#'
#' Retrieves the sites × samples matrix of read depth from a
#' \code{\link{commaData}} object.
#'
#' @param x A \code{commaData} object.
#' @param shift Not used; inherited from the \code{IRanges::coverage} generic.
#' @param width Not used; inherited from the \code{IRanges::coverage} generic.
#' @param weight Not used; inherited from the \code{IRanges::coverage} generic.
#' @param ... Not used.
#'
#' @return An integer matrix with rows corresponding to methylation sites and
#'   columns corresponding to samples.
#'
#' @seealso \code{\link{methylation}}, \code{\link{siteInfo}}
#'
#' @examples
#' data(comma_example_data)
#' cov <- coverage(comma_example_data)
#' summary(as.vector(cov))
#'
#' @export
setMethod("coverage", "commaData", function(x, shift = 0L, width = NULL, weight = 1L, ...) {
    assay(x, "coverage")
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
#'   Dorado/Megalodon callers), and \code{mod_context} (the composite
#'   modification context, e.g., \code{"6mA_GATC"}). May contain additional
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
        row.names   = names(rr)
    )
    # Add computed mod_context column if not already present
    if (!"mod_context" %in% colnames(df) &&
        "mod_type" %in% colnames(mc) && "motif" %in% colnames(mc)) {
        df$mod_context <- .computeModContext(mc$mod_type, mc$motif)
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
    sort(unique(rowData(object)$mod_type))
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

# ─── subset() ────────────────────────────────────────────────────────────────

#' Subset a commaData object by condition, modification type, or chromosome
#'
#' A convenience function for filtering a \code{\link{commaData}} object by
#' common criteria. For arbitrary index-based subsetting, use \code{[}.
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
#' six_ma <- subset(comma_example_data, mod_type = "6mA")
#' modTypes(six_ma)
#'
#' # Only GATC-context sites
#' gatc <- subset(comma_example_data, motif = "GATC")
#' nrow(gatc)
#'
#' # Filter by mod_context (equivalent to the above for modkit data)
#' gatc2 <- subset(comma_example_data, mod_context = "6mA_GATC")
#' nrow(gatc2)
#'
#' @export
setGeneric("subset", function(x, ...) standardGeneric("subset"))

#' @rdname subset
setMethod("subset", "commaData", function(x, mod_type = NULL,
                                           condition = NULL,
                                           chrom = NULL,
                                           motif = NULL,
                                           mod_context = NULL, ...) {
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
    ifelse(is.na(motif), mod_type, paste(mod_type, motif, sep = "_"))
}
