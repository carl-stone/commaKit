#' @importFrom methods setGeneric setMethod
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom BiocGenerics annotation
#' @importFrom IRanges coverage
#' @importFrom GenomeInfoDb genome
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
#' @seealso \code{\link{coverage}}, \code{\link{siteInfo}},
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
#' @param object A \code{commaData} object.
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
#' Equivalent to \code{rowData(object)} but returns a plain \code{data.frame}.
#'
#' @param object A \code{commaData} object.
#'
#' @return A \code{data.frame} with one row per methylation site. Always
#'   contains columns \code{chrom}, \code{position}, \code{strand}, and
#'   \code{mod_type}. May contain additional annotation columns added by
#'   \code{\link[=annotateSites]{annotateSites()}}.
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
    as.data.frame(rowData(object))
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

# ─── genome() ────────────────────────────────────────────────────────────────

#' Accessor for genome size information
#'
#' Returns the chromosome sizes stored in a \code{\link{commaData}} object.
#'
#' @param object A \code{commaData} object.
#'
#' @return A named integer vector of chromosome sizes
#'   (chromosome name → length in bp), or \code{NULL} if no genome information
#'   was provided at construction.
#'
#' @examples
#' data(comma_example_data)
#' genome(comma_example_data)
#'
#' @export
setMethod("genome", "commaData", function(x) {
    x@genomeInfo
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
    object@annotation
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
    object@motifSites
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

    new("commaData",
        se_sub,
        genomeInfo = x@genomeInfo,
        annotation = x@annotation,
        motifSites = x@motifSites
    )
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
#' @export
setGeneric("subset", function(x, ...) standardGeneric("subset"))

#' @rdname subset
setMethod("subset", "commaData", function(x, mod_type = NULL,
                                           condition = NULL,
                                           chrom = NULL, ...) {
    rd <- rowData(x)
    cd <- colData(x)

    # Site filter
    site_keep <- rep(TRUE, nrow(x))
    if (!is.null(mod_type)) {
        site_keep <- site_keep & (rd$mod_type %in% mod_type)
    }
    if (!is.null(chrom)) {
        site_keep <- site_keep & (rd$chrom %in% chrom)
    }

    # Sample filter
    samp_keep <- rep(TRUE, ncol(x))
    if (!is.null(condition)) {
        samp_keep <- samp_keep & (cd$condition %in% condition)
    }

    x[site_keep, samp_keep]
})
