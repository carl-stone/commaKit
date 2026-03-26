#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData assay
NULL

#' Windowed sequencing depth across the genome
#'
#' Bins the genome into non-overlapping windows and computes the average (or
#' median) sequencing depth in each window for each sample. Returns a tidy
#' \code{data.frame} suitable for plotting or downstream QC analysis.
#'
#' @param object A \code{\link{commaData}} object.
#' @param window Positive integer. Window size in base pairs.
#' @param method Character string. Aggregation method within each window.
#'   One of \code{"mean"} (default) or \code{"median"}.
#' @param log2_transform Logical. If \code{TRUE}, the depth values are
#'   log2-transformed (using \eqn{log2(depth + 1)} to handle zeros). Default:
#'   \code{FALSE}.
#'
#' @return A \code{data.frame} with one row per (chromosome window, sample),
#'   containing:
#'   \describe{
#'     \item{\code{chrom}}{Chromosome name.}
#'     \item{\code{window_start}}{First base of the window (1-based).}
#'     \item{\code{window_end}}{Last base of the window (1-based).}
#'     \item{\code{sample_name}}{Sample identifier.}
#'     \item{\code{depth}}{Mean or median sequencing depth in the window.}
#'     \item{\code{log2_depth}}{Log2-transformed depth (only present if
#'       \code{log2_transform = TRUE}).}
#'   }
#'
#' @details
#' Depth is computed only at positions with observed methylation sites. Windows
#' with no sites have \code{depth = NA}.
#'
#' If genome size information is stored in \code{genome(object)}, windows are
#' sized to fit the chromosomes exactly (the last window may be smaller than
#' \code{window}). If genome information is absent, only the range spanned by
#' observed sites is covered.
#'
#' @examples
#' data(comma_example_data)
#' cd <- coverageDepth(comma_example_data, window = 10000L)
#' head(cd)
#'
#' @seealso \code{\link{varianceByDepth}}, \code{\link{methylomeSummary}}
#'
#' @export
coverageDepth <- function(object,
                          window,
                          method         = c("mean", "median"),
                          log2_transform = FALSE) {
    # ── Input validation ──────────────────────────────────────────────────────
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    method <- match.arg(method)
    if (missing(window) || !is.numeric(window) || length(window) != 1 || window < 1) {
        stop("'window' must be a positive integer specifying window size in bp.")
    }
    window <- as.integer(window)

    rd         <- as.data.frame(rowData(object))
    cov_mat    <- coverage(object)
    sample_nms <- colnames(cov_mat)
    agg_fn     <- if (method == "mean") mean else stats::median

    # Determine chromosome sizes
    genome_info <- genome(object)
    chroms <- if (!is.null(genome_info)) names(genome_info) else unique(rd$chrom)

    result_list <- vector("list", length(chroms))

    for (ci in seq_along(chroms)) {
        chr      <- chroms[ci]
        chr_size <- if (!is.null(genome_info)) genome_info[[chr]] else max(rd$position[rd$chrom == chr])

        # Define window breakpoints
        win_starts <- seq(1L, as.integer(chr_size), by = window)
        win_ends   <- pmin(win_starts + window - 1L, as.integer(chr_size))

        chr_idx  <- which(rd$chrom == chr)
        chr_pos  <- rd$position[chr_idx]

        sample_dfs <- vector("list", length(sample_nms))

        for (si in seq_along(sample_nms)) {
            samp     <- sample_nms[si]
            chr_cov  <- as.numeric(cov_mat[chr_idx, samp])

            depths <- vapply(seq_along(win_starts), function(wi) {
                in_win <- chr_pos >= win_starts[wi] & chr_pos <= win_ends[wi]
                if (!any(in_win)) return(NA_real_)
                agg_fn(chr_cov[in_win], na.rm = TRUE)
            }, numeric(1))

            df <- data.frame(
                chrom        = chr,
                window_start = win_starts,
                window_end   = win_ends,
                sample_name  = samp,
                depth        = depths,
                stringsAsFactors = FALSE
            )
            if (log2_transform) {
                df$log2_depth <- log2(df$depth + 1)
            }
            sample_dfs[[si]] <- df
        }

        result_list[[ci]] <- do.call(rbind, sample_dfs)
    }

    do.call(rbind, result_list)
}


#' Methylation variance as a function of sequencing depth
#'
#' For each coverage level (or coverage bin), computes the variance of
#' methylation beta values across sites at that depth. This is useful for
#' diagnosing whether low-coverage sites have inflated methylation variance
#' and for setting appropriate coverage thresholds.
#'
#' @param object A \code{\link{commaData}} object.
#' @param coverage_bins Integer vector specifying the coverage levels to
#'   include. If \code{NULL} (default), all unique coverage levels observed
#'   across all samples are used. Useful to pass \code{5:30} to focus on a
#'   specific depth range.
#' @param mod_type Character string or \code{NULL}. If provided, only sites
#'   of the specified modification type are included. Default: \code{NULL}
#'   (all types).
#'
#' @return A \code{data.frame} with one row per (coverage level, sample),
#'   containing:
#'   \describe{
#'     \item{\code{coverage}}{Sequencing depth (integer).}
#'     \item{\code{sample_name}}{Sample identifier.}
#'     \item{\code{variance}}{Variance of beta values at sites with exactly
#'       this coverage level. \code{NA} if fewer than 2 sites are at this
#'       level.}
#'     \item{\code{n_sites}}{Number of sites at this coverage level.}
#'   }
#'
#' @examples
#' data(comma_example_data)
#' vd <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
#' head(vd)
#'
#' @seealso \code{\link{coverageDepth}}, \code{\link{methylomeSummary}}
#'
#' @export
varianceByDepth <- function(object,
                            coverage_bins = NULL,
                            mod_type      = NULL,
                            motif         = NULL) {
    # ── Input validation ──────────────────────────────────────────────────────
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }

    if (!is.null(mod_type)) {
        available <- modTypes(object)
        if (!mod_type %in% available) {
            stop(
                "'mod_type' = '", mod_type, "' not found in object. ",
                "Available types: ", paste(available, collapse = ", ")
            )
        }
        object <- subset(object, mod_type = mod_type)
    }

    if (!is.null(motif)) {
        available_m <- motifs(object)
        bad_m <- setdiff(motif, available_m)
        if (length(bad_m) > 0L) {
            stop(
                "'motif' value(s) not found in object: ",
                paste(bad_m, collapse = ", "),
                ". Available: ", paste(available_m, collapse = ", ")
            )
        }
        object <- subset(object, motif = motif)
    }

    methyl_mat <- methylation(object)
    cov_mat    <- coverage(object)
    sample_nms <- colnames(methyl_mat)

    result_list <- vector("list", length(sample_nms))

    for (si in seq_along(sample_nms)) {
        samp    <- sample_nms[si]
        betas   <- methyl_mat[, samp]
        covs    <- as.integer(cov_mat[, samp])

        # Restrict to specified coverage levels
        levels <- if (!is.null(coverage_bins)) {
            as.integer(coverage_bins)
        } else {
            sort(unique(covs[!is.na(covs)]))
        }

        rows <- lapply(levels, function(lv) {
            idx <- which(covs == lv & !is.na(betas))
            n   <- length(idx)
            data.frame(
                coverage    = lv,
                sample_name = samp,
                variance    = if (n >= 2) stats::var(betas[idx]) else NA_real_,
                n_sites     = n,
                stringsAsFactors = FALSE
            )
        })

        result_list[[si]] <- do.call(rbind, rows)
    }

    do.call(rbind, result_list)
}
