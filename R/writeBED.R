#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData assay
NULL

#' Export methylation data as a BED file
#'
#' Writes per-site methylation beta values for a single sample from a
#' \code{\link{commaData}} object to a 9-column BED file (BED9 format),
#' suitable for visualisation in IGV, UCSC Genome Browser, or other genome
#' browsers that support the \code{itemRGB} field.
#'
#' @param object A \code{\link{commaData}} object.
#' @param file Character string. Path to the output BED file. The file will be
#'   created or overwritten.
#' @param sample Character string. Name of the sample to export. Must match a
#'   column name in \code{methylation(object)}.
#' @param mod_type Character string or \code{NULL}. If provided, only sites of
#'   the specified modification type are written (e.g., \code{"6mA"}). If
#'   \code{NULL} (default), all sites are written.
#' @param rgb_scale Logical. If \code{TRUE} (default), an \code{itemRGB}
#'   column is added based on the methylation score using a blue-to-red
#'   gradient (low = blue, high = red). If \code{FALSE}, \code{itemRGB} is
#'   set to \code{"0,0,0"} for all sites.
#' @param track_name Character string. Name shown in the genome browser track
#'   header. Defaults to \code{sample}.
#' @param track_description Character string. Description shown in the track
#'   header. Defaults to \code{"methylation beta values"}.
#'
#' @return Invisibly returns the path to the written file (\code{file}).
#'
#' @details
#' The BED score field (column 5) contains the beta value multiplied by 1000
#' and rounded to the nearest integer (range 0–1000), which is the standard
#' convention for methylation BED files. The \code{itemRGB} colour gradient
#' transitions as follows: score ≤ 200 = blue (0,0,255), score ≤ 400 =
#' blue-purple (83,0,172), score ≤ 600 = purple (167,0,85), score ≤ 800 =
#' red-purple (222,0,28), score ≤ 1000 = red (250,0,0).
#'
#' Sites with \code{NA} methylation (below the coverage threshold) are
#' excluded from the output.
#'
#' @examples
#' \dontrun{
#' data(comma_example_data)
#' writeBED(comma_example_data,
#'          file    = "ctrl_1_methylation.bed",
#'          sample  = "ctrl_1",
#'          mod_type = "6mA")
#' }
#'
#' @seealso \code{\link{methylation}}, \code{\link{siteInfo}}
#'
#' @export
writeBED <- function(object,
                     file,
                     sample,
                     mod_type          = NULL,
                     rgb_scale         = TRUE,
                     track_name        = sample,
                     track_description = "methylation beta values") {
    # ── Input validation ──────────────────────────────────────────────────────
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    if (missing(file) || !is.character(file) || length(file) != 1) {
        stop("'file' must be a single character string specifying the output path.")
    }
    if (missing(sample) || !is.character(sample) || length(sample) != 1) {
        stop("'sample' must be a single character string matching a sample in the object.")
    }

    available_samples <- colnames(methylation(object))
    if (!sample %in% available_samples) {
        stop(
            "'sample' = '", sample, "' not found in object. ",
            "Available samples: ", paste(available_samples, collapse = ", ")
        )
    }

    # ── Filter by mod_type ────────────────────────────────────────────────────
    if (!is.null(mod_type)) {
        object <- subset(object, mod_type = mod_type)
        if (nrow(object) == 0) {
            stop("No sites remain after filtering for mod_type = '", mod_type, "'.")
        }
    }

    # ── Extract data ──────────────────────────────────────────────────────────
    rd    <- as.data.frame(rowData(object))
    betas <- methylation(object)[, sample]

    # Exclude NA sites (below coverage threshold)
    keep  <- !is.na(betas)
    rd    <- rd[keep, , drop = FALSE]
    betas <- betas[keep]

    if (nrow(rd) == 0) {
        warning("No sites with non-NA methylation for sample '", sample, "'. ",
                "Writing empty BED file.")
        cat(paste0('track name="', track_name, '" description="', track_description,
                   '" itemRgb="On"\n'),
            file = file)
        return(invisible(file))
    }

    # ── Build BED rows ────────────────────────────────────────────────────────
    # BED is 0-based; our positions are 1-based
    score <- as.integer(round(betas * 1000))

    if (rgb_scale) {
        rgb <- rep("250,0,0", length(score))
        rgb[score <= 800] <- "222,0,28"
        rgb[score <= 600] <- "167,0,85"
        rgb[score <= 400] <- "83,0,172"
        rgb[score <= 200] <- "0,0,255"
    } else {
        rgb <- rep("0,0,0", length(score))
    }

    bed_df <- data.frame(
        chrom       = rd$chrom,
        chromStart  = rd$position - 1L,     # 0-based
        chromEnd    = rd$position,           # half-open
        name        = rd$position,
        score       = score,
        strand      = rd$strand,
        thickStart  = rd$position - 1L,
        thickEnd    = rd$position,
        itemRGB     = rgb,
        stringsAsFactors = FALSE
    )

    # ── Write output ──────────────────────────────────────────────────────────
    track_line <- paste0('track name="', track_name, '" description="',
                         track_description, '" itemRgb="On"')
    cat(track_line, "\n", file = file, sep = "")

    utils::write.table(
        bed_df,
        file      = file,
        sep       = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote     = FALSE,
        append    = TRUE
    )

    invisible(file)
}
