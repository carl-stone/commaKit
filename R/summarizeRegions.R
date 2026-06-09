#' Summarize methylation counts over genomic regions
#'
#' Aggregates site-level count evidence over user-supplied genomic regions.
#' This is a descriptive regional summary, not a DMR caller, smoothing method,
#' or region-level statistical test.
#'
#' @param object A \code{\link{commaData}} object with count evidence assays.
#' @param regions A \code{\link[GenomicRanges]{GRanges}} object defining
#'   regions to summarize.
#' @param min_sites Integer. Minimum number of overlapping sites with usable
#'   count evidence required before \code{region_methylation} is reported.
#'   Regions with fewer usable sites return \code{NA} methylation for that
#'   sample. Default: \code{1L}.
#' @param mod_type,motif,mod_context Optional site filters. If provided, only
#'   sites matching these modification annotations are summarized.
#'
#' @return A tidy \code{data.frame} with one row per region × sample containing
#'   region coordinates/metadata, \code{sample_name}, \code{n_sites}, summed
#'   count evidence, and \code{region_methylation} computed as
#'   \code{sum(mod_counts) / sum(valid_coverage)}.
#'
#' @details
#' The primary regional methylation statistic is count-based:
#' \deqn{region\_methylation = sum(mod\_counts) / sum(valid\_coverage)}
#' where \code{valid_coverage} is the \code{coverage} assay for sites with
#' non-missing modified counts and coverage. For modkit-style inputs this
#' corresponds to \code{Nvalid_cov = Nmod + Nother_mod + Ncanonical}.
#'
#' \code{summarizeRegions()} intentionally does not average beta values by
#' default. Beta-only/reconstructed-count import paths should record provenance
#' and message when counts are reconstructed; regional summaries consume the
#' count assays present in the \code{commaData} object.
#'
#' @examples
#' data(comma_example_data)
#' regions <- GenomicRanges::GRanges(
#'     seqnames = "chr_sim",
#'     ranges = IRanges::IRanges(start = 1, end = 5000)
#' )
#' summarizeRegions(comma_example_data, regions, mod_type = "6mA")
#'
#' @export
summarizeRegions <- function(object, regions, min_sites = 1L,
                             mod_type = NULL, motif = NULL,
                             mod_context = NULL) {
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    if (!is(regions, "GRanges")) {
        stop("'regions' must be a GRanges object.")
    }
    if (!is.numeric(min_sites) || length(min_sites) != 1L ||
            is.na(min_sites) || min_sites < 0) {
        stop("'min_sites' must be a single non-negative integer.")
    }
    min_sites <- as.integer(min_sites)

    assay_names <- SummarizedExperiment::assayNames(object)
    missing_assays <- setdiff(c("mod_counts", "coverage"), assay_names)
    if (length(missing_assays) > 0L) {
        stop(
            "summarizeRegions() requires count evidence assays: ",
            paste(missing_assays, collapse = ", "), "."
        )
    }

    rr <- SummarizedExperiment::rowRanges(object)
    rd <- SummarizedExperiment::rowData(object)
    keep <- rep(TRUE, nrow(object))

    if (!is.null(mod_type)) {
        .validateModType(mod_type, object)
        keep <- keep & rd$mod_type %in% mod_type
    }
    if (!is.null(motif)) {
        keep <- keep & rd$motif %in% motif
    }
    if (!is.null(mod_context)) {
        keep <- keep & .computeModContext(rd$mod_type, rd$motif) %in% mod_context
    }

    site_idx <- which(keep)
    overlaps <- GenomicRanges::findOverlaps(regions, rr[site_idx])
    overlap_region <- S4Vectors::queryHits(overlaps)
    overlap_site <- site_idx[S4Vectors::subjectHits(overlaps)]

    mod_counts <- SummarizedExperiment::assay(object, "mod_counts")
    valid_coverage <- SummarizedExperiment::assay(object, "coverage")
    canonical_counts <- if ("canonical_counts" %in% assay_names) {
        SummarizedExperiment::assay(object, "canonical_counts")
    } else {
        NULL
    }
    other_mod_counts <- if ("other_mod_counts" %in% assay_names) {
        SummarizedExperiment::assay(object, "other_mod_counts")
    } else {
        NULL
    }

    cd <- as.data.frame(SummarizedExperiment::colData(object))
    sample_names <- if ("sample_name" %in% colnames(cd)) {
        as.character(cd$sample_name)
    } else {
        colnames(object)
    }
    if (is.null(sample_names)) {
        sample_names <- paste0("sample_", seq_len(ncol(object)))
    }

    region_df <- .regionSummaryFrame(regions)
    rows <- vector("list", length(regions) * ncol(object))
    row_i <- 0L
    for (region_i in seq_along(regions)) {
        sites <- overlap_site[overlap_region == region_i]
        for (sample_i in seq_len(ncol(object))) {
            row_i <- row_i + 1L
            rows[[row_i]] <- cbind(
                region_df[region_i, , drop = FALSE],
                .summarizeRegionSample(
                    sites = sites,
                    sample_i = sample_i,
                    sample_name = sample_names[[sample_i]],
                    mod_counts = mod_counts,
                    valid_coverage = valid_coverage,
                    canonical_counts = canonical_counts,
                    other_mod_counts = other_mod_counts,
                    min_sites = min_sites
                ),
                stringsAsFactors = FALSE
            )
        }
    }

    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

.regionSummaryFrame <- function(regions) {
    region_id <- names(regions)
    if (is.null(region_id)) {
        region_id <- paste0("region_", seq_along(regions))
    }
    missing_id <- is.na(region_id) | !nzchar(region_id)
    region_id[missing_id] <- paste0("region_", which(missing_id))

    out <- data.frame(
        region_id = region_id,
        seqnames = as.character(GenomicRanges::seqnames(regions)),
        start = GenomicRanges::start(regions),
        end = GenomicRanges::end(regions),
        width = GenomicRanges::width(regions),
        strand = as.character(GenomicRanges::strand(regions)),
        stringsAsFactors = FALSE
    )
    region_mcols <- as.data.frame(GenomicRanges::mcols(regions))
    if (ncol(region_mcols) > 0L) {
        colnames(region_mcols) <- paste0("region_", colnames(region_mcols))
        out <- cbind(out, region_mcols, stringsAsFactors = FALSE)
    }
    out
}

.summarizeRegionSample <- function(sites, sample_i, sample_name,
                                   mod_counts, valid_coverage,
                                   canonical_counts = NULL,
                                   other_mod_counts = NULL,
                                   min_sites = 1L) {
    if (length(sites) == 0L) {
        valid <- logical(0)
        mod_vals <- cov_vals <- numeric(0)
    } else {
        mod_vals <- mod_counts[sites, sample_i]
        cov_vals <- valid_coverage[sites, sample_i]
        valid <- !is.na(mod_vals) & !is.na(cov_vals)
        mod_vals <- mod_vals[valid]
        cov_vals <- cov_vals[valid]
    }

    n_sites <- length(mod_vals)
    total_mod <- sum(mod_vals, na.rm = TRUE)
    total_cov <- sum(cov_vals, na.rm = TRUE)
    region_methylation <- if (n_sites >= min_sites && total_cov > 0) {
        total_mod / total_cov
    } else {
        NA_real_
    }

    out <- data.frame(
        sample_name = sample_name,
        n_sites = n_sites,
        total_mod_counts = total_mod,
        total_valid_coverage = total_cov,
        region_methylation = region_methylation,
        stringsAsFactors = FALSE
    )

    if (!is.null(canonical_counts)) {
        canon_vals <- if (length(sites) == 0L) numeric(0) else canonical_counts[sites, sample_i][valid]
        out$total_canonical_counts <- .sumOrNA(canon_vals)
    }
    if (!is.null(other_mod_counts)) {
        other_vals <- if (length(sites) == 0L) numeric(0) else other_mod_counts[sites, sample_i][valid]
        out$total_other_mod_counts <- .sumOrNA(other_vals)
    }
    out
}

.sumOrNA <- function(x) {
    if (length(x) == 0L) {
        return(0)
    }
    if (all(is.na(x))) {
        return(NA_real_)
    }
    sum(x, na.rm = TRUE)
}
