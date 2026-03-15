#' @importFrom GenomicRanges GRanges mcols mcols<-
NULL

#' Load genomic feature annotations from a GFF3 or BED file
#'
#' Reads a GFF3 or BED annotation file and returns a
#' \code{\link[GenomicRanges]{GRanges}} object with standardized metadata
#' columns. The result can be passed directly to the \code{annotation} argument
#' of \code{\link{commaData}}.
#'
#' @param file Character string. Path to a GFF3 (\code{.gff}, \code{.gff3},
#'   \code{.gff.gz}, \code{.gff3.gz}) or BED (\code{.bed}) file.
#' @param feature_types Character vector or \code{NULL}. If provided, only
#'   features with a matching \code{type} (GFF3) or \code{name} (BED) are
#'   retained. Common GFF3 types include \code{"gene"}, \code{"CDS"},
#'   \code{"rRNA"}, \code{"tRNA"}. \code{NULL} retains all features.
#' @param ... Additional arguments passed to \code{rtracklayer::import()}.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object. The \code{mcols}
#'   always include:
#'   \describe{
#'     \item{\code{feature_type}}{Character. The feature type (from GFF3
#'       \code{type} column, or \code{"region"} for BED).}
#'     \item{\code{name}}{Character. Feature name or identifier (from GFF3
#'       \code{Name} or \code{ID} attribute, or BED \code{name} column).}
#'   }
#'   Additional metadata columns from the source file are preserved.
#'
#' @details
#' This function requires the \code{rtracklayer} package (Bioconductor):
#' \preformatted{
#'   BiocManager::install("rtracklayer")
#' }
#' Both NCBI-style GFF3 (with \code{gene_biotype}, \code{product} attributes)
#' and Ensembl-style GFF3 are supported.
#'
#' @examples
#' # Load the bundled example GFF3 annotation
#' gff_file <- system.file("extdata", "example.gff3", package = "comma")
#' if (requireNamespace("rtracklayer", quietly = TRUE)) {
#'   ann <- loadAnnotation(gff_file)
#'   ann <- loadAnnotation(gff_file, feature_types = "gene")
#' }
#'
#' \donttest{
#' # Load only genes and CDS from your own file
#' ann <- loadAnnotation("my_genome.gff3", feature_types = c("gene", "CDS"))
#' }
#'
#' @export
loadAnnotation <- function(file, feature_types = NULL, ...) {
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        stop(
            "Package 'rtracklayer' is required for loadAnnotation(). ",
            "Install it with: BiocManager::install('rtracklayer')"
        )
    }

    if (!is.character(file) || length(file) != 1) {
        stop("file must be a single character string path")
    }
    if (!file.exists(file)) {
        stop("Annotation file not found: ", file)
    }

    # ── Detect format from extension ────────────────────────────────────────
    ext <- .annotationFileExt(file)

    gr <- tryCatch(
        rtracklayer::import(file, ...),
        error = function(e) stop("Failed to read annotation file '", file, "': ", e$message)
    )

    if (length(gr) == 0L) {
        message("Note: annotation file '", file, "' contains no features")
        return(gr)
    }

    # ── Standardize mcols ───────────────────────────────────────────────────
    gr <- .standardizeAnnotationMcols(gr, ext)

    # ── Filter by feature_types ─────────────────────────────────────────────
    if (!is.null(feature_types)) {
        keep <- GenomicRanges::mcols(gr)$feature_type %in% feature_types
        gr   <- gr[keep]
        if (length(gr) == 0L) {
            warning(
                "No features of the requested type(s) found in '", file, "': ",
                paste(feature_types, collapse = ", ")
            )
        }
    }

    gr
}

#' Detect annotation file format from file extension
#' @return Character string, either \code{"gff"} or \code{"bed"}.
#' @keywords internal
.annotationFileExt <- function(file) {
    # Strip compression suffix first
    base <- sub("\\.gz$|\\.bz2$|\\.xz$", "", file, ignore.case = TRUE)
    ext  <- tolower(tools::file_ext(base))
    if (ext %in% c("gff", "gff3")) return("gff")
    if (ext == "bed") return("bed")
    # Default: try GFF
    warning("Unrecognized file extension for annotation file '", file,
            "'. Attempting to read as GFF3.")
    "gff"
}

#' Standardize annotation GRanges mcols to always have feature_type and name
#' @return A \code{GRanges} object with standardized \code{feature_type} and
#'   \code{name} metadata columns.
#' @keywords internal
.standardizeAnnotationMcols <- function(gr, ext) {
    mc <- GenomicRanges::mcols(gr)

    if (ext == "gff") {
        # GFF3: 'type' column → feature_type
        if ("type" %in% colnames(mc)) {
            mc$feature_type <- as.character(mc$type)
        } else {
            mc$feature_type <- NA_character_
        }
        # Name: prefer Name attribute, fall back to ID
        if ("Name" %in% colnames(mc)) {
            mc$name <- as.character(mc$Name)
        } else if ("ID" %in% colnames(mc)) {
            mc$name <- as.character(mc$ID)
        } else {
            mc$name <- NA_character_
        }
    } else {
        # BED
        mc$feature_type <- "region"
        if ("name" %in% colnames(mc)) {
            mc$name <- as.character(mc$name)
        } else {
            mc$name <- NA_character_
        }
    }

    GenomicRanges::mcols(gr) <- mc
    gr
}
