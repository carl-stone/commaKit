#' @importFrom zoo rollapply
#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData assay
NULL

#' Sliding window methylation summary along the genome
#'
#' Computes a per-position sliding window statistic (median or mean) of
#' methylation beta values for each sample in a \code{\link{commaData}}
#' object. The genome size for each chromosome is read from
#' \code{genome(object)}, so no organism-specific values are ever hardcoded.
#' Circular genome wrap-around is supported via
#' \code{.circularIndex()}.
#'
#' @param object A \code{\link{commaData}} object. Must have genome size
#'   information in \code{genome(object)} (i.e., it must have been constructed
#'   with a \code{genome} argument).
#' @param window Positive integer. Window size in base pairs. The smoothed
#'   value at position \eqn{p} is computed from positions
#'   \eqn{[p - \lfloor w/2 \rfloor,\; p + \lfloor w/2 \rfloor]}.
#' @param stat Character string. Summary statistic to apply within each window.
#'   One of \code{"median"} (default) or \code{"mean"}.
#' @param mod_type Character string or \code{NULL}. If provided, only sites of
#'   the specified modification type (e.g., \code{"6mA"}) are included in the
#'   smoothing. If \code{NULL} (default), all sites are used.
#' @param circular Logical. If \code{TRUE} (default), positions at the ends of
#'   each chromosome are wrapped around so that the window at position 1 can
#'   draw from positions near the chromosome end, and vice versa. Appropriate
#'   for circular bacterial chromosomes.
#'
#' @return A \code{data.frame} with one row per (chromosome, position, sample)
#'   combination, containing:
#'   \describe{
#'     \item{\code{chrom}}{Chromosome name (character).}
#'     \item{\code{position}}{Genomic position, 1-based (integer).}
#'     \item{\code{sample_name}}{Sample identifier (character).}
#'     \item{\code{window_stat}}{Smoothed beta value (numeric). The column is
#'       named \code{window_median} or \code{window_mean} depending on
#'       \code{stat}.}
#'   }
#'
#' @details
#' Because most positions in a genome have no methylation site, the beta value
#' vector for each chromosome-sample pair is sparse (mostly \code{NA}).
#' \code{zoo::rollapply} is called with \code{na.rm = TRUE} so that windows
#' spanning regions with no data still produce a value where at least one site
#' is present.
#'
#' Positions where every site in the window is \code{NA} (i.e., no coverage)
#' remain \code{NA} in the output.
#'
#' @seealso \code{\link{methylation}}, \code{\link[GenomeInfoDb]{genome}},
#'   \code{\link{methylomeSummary}}
#'
#' @examples
#' data(comma_example_data)
#' \donttest{
#' sw <- slidingWindow(comma_example_data, window = 5000L)
#' head(sw)
#' # Filter to one sample
#' sw_ctrl1 <- sw[sw$sample_name == "ctrl_1", ]
#' }
#'
#' @export
slidingWindow <- function(object,
                          window,
                          stat     = c("median", "mean"),
                          mod_type = NULL,
                          circular = TRUE) {
    # ── Input validation ──────────────────────────────────────────────────────
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    stat <- match.arg(stat)

    if (missing(window) || !is.numeric(window) || length(window) != 1 ||
        window < 1) {
        stop("'window' must be a positive integer specifying window size in bp.")
    }
    window <- as.integer(window)

    genome_info <- genome(object)
    if (is.null(genome_info) || length(genome_info) == 0) {
        stop(
            "genome(object) is NULL or empty. ",
            "Provide genome size information when constructing the commaData object."
        )
    }

    min_chr <- min(genome_info)
    if (window > min_chr) {
        stop(
            "'window' (", window, " bp) exceeds the smallest chromosome size (",
            min_chr, " bp). Reduce window size."
        )
    }

    # ── Filter by mod_type if requested ──────────────────────────────────────
    if (!is.null(mod_type)) {
        object <- subset(object, mod_type = mod_type)
        if (nrow(object) == 0) {
            stop("No sites remain after filtering for mod_type = '", mod_type, "'.")
        }
    }

    rd          <- as.data.frame(rowData(object))
    methyl_mat  <- methylation(object)
    sample_nms  <- colnames(methyl_mat)
    stat_colnm  <- paste0("window_", stat)
    stat_fn     <- if (stat == "median") stats::median else base::mean

    result_list <- vector("list", length(genome_info))

    for (ci in seq_along(genome_info)) {
        chr      <- names(genome_info)[ci]
        chr_size <- genome_info[ci]

        # Sites on this chromosome
        chr_idx  <- which(rd$chrom == chr)

        # Per-sample smoothing
        sample_dfs <- vector("list", length(sample_nms))

        for (si in seq_along(sample_nms)) {
            samp <- sample_nms[si]

            # Build full-length beta vector (1 to chr_size), NA where no site
            beta_vec <- rep(NA_real_, chr_size)
            if (length(chr_idx) > 0) {
                positions <- rd$position[chr_idx]
                # Clip positions to valid range
                valid     <- positions >= 1L & positions <= chr_size
                beta_vec[positions[valid]] <- methyl_mat[chr_idx[valid], samp]
            }

            # Circular padding: prepend tail and append head
            if (circular) {
                half    <- as.integer(floor(window / 2))
                padded  <- c(beta_vec[(chr_size - half + 1L):chr_size],
                              beta_vec,
                              beta_vec[1L:half])
                smoothed_padded <- zoo::rollapply(
                    padded,
                    width   = window,
                    FUN     = function(x) stat_fn(x, na.rm = TRUE),
                    partial = TRUE,
                    align   = "center",
                    fill    = NA_real_
                )
                smoothed <- smoothed_padded[(half + 1L):(half + chr_size)]
            } else {
                smoothed <- zoo::rollapply(
                    beta_vec,
                    width   = window,
                    FUN     = function(x) stat_fn(x, na.rm = TRUE),
                    partial = TRUE,
                    align   = "center",
                    fill    = NA_real_
                )
            }

            # NaN → NA (produced when all values in window are NA and na.rm=TRUE)
            smoothed[is.nan(smoothed)] <- NA_real_

            sample_dfs[[si]] <- data.frame(
                chrom       = chr,
                position    = seq_len(chr_size),
                sample_name = samp,
                stringsAsFactors = FALSE
            )
            sample_dfs[[si]][[stat_colnm]] <- as.numeric(smoothed)
        }

        result_list[[ci]] <- do.call(rbind, sample_dfs)
    }

    do.call(rbind, result_list)
}
