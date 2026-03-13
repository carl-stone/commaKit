NULL

# ─── Dorado BAM parser ────────────────────────────────────────────────────────

#' Parse a Dorado BAM file with MM/ML base modification tags
#'
#' Reads a Dorado-aligned BAM file containing MM (base modification) and
#' ML (modification likelihood) tags, and aggregates per-read modification
#' probabilities into per-site beta values. This is an internal function called
#' by \code{\link{commaData}} when \code{caller = "dorado"}.
#'
#' @details
#' The function reads the BAM using \code{\link[Rsamtools]{scanBam}} and
#' processes the MM/ML tags from each aligned read. The MM tag encodes which
#' bases carry modifications and their positions within the read sequence
#' (as inter-base offsets); the parallel ML tag provides modification
#' probabilities (0–255 scaled to 0–1).
#'
#' A base is called modified when its ML probability is \eqn{> 0.5}.
#' Per-site statistics are aggregated by counting modified reads
#' (\eqn{n_{\text{mod}}}) and total reads at each genomic position
#' (\eqn{n_{\text{total}}}). The beta value is
#' \eqn{\beta = n_{\text{mod}} / n_{\text{total}}}.
#'
#' CIGAR operations are used to map read positions to reference coordinates.
#' Soft-clipped (\code{S}) and hard-clipped (\code{H}) bases are excluded;
#' insertions (\code{I}) consume read positions without advancing the reference.
#'
#' \strong{Recommended workflow:} For most users, it is simpler to first run
#' \code{modkit pileup} on the Dorado BAM, then load the resulting BED file
#' with \code{caller = "modkit"}. Direct BAM parsing is provided for users who
#' prefer to avoid the modkit step.
#'
#' @param file Character string. Path to the Dorado-aligned BAM file. An
#'   accompanying \code{.bai} index file must exist (required by
#'   \code{Rsamtools}).
#' @param sample_name Character string. Sample name (used in messages only).
#' @param mod_type Character vector or \code{NULL}. If provided, only sites
#'   with a matching modification type are returned. \code{NULL} retains all
#'   types.
#' @param min_coverage Integer. Minimum read depth required to retain a site.
#'   Default \code{5L}.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{\code{chrom}}{Chromosome name (character).}
#'     \item{\code{position}}{1-based genomic position (integer).}
#'     \item{\code{strand}}{Strand, \code{"+"} or \code{"-"} (character).}
#'     \item{\code{mod_type}}{Modification type: \code{"6mA"}, \code{"5mC"},
#'       or \code{"4mC"} (character).}
#'     \item{\code{beta}}{Proportion of reads called methylated, range 0–1
#'       (numeric).}
#'     \item{\code{coverage}}{Total read depth at this site (integer).}
#'   }
#'
#' @keywords internal
.parseDorado <- function(file, sample_name, mod_type = NULL, min_coverage = 5L) {
    # ── Validate inputs ───────────────────────────────────────────────────────
    if (!is.character(file) || length(file) != 1L) {
        stop("'file' must be a single character string path.")
    }
    if (!file.exists(file)) {
        stop("Dorado BAM file not found: ", file)
    }
    min_coverage <- as.integer(min_coverage)

    # ── Read BAM ──────────────────────────────────────────────────────────────
    bam_file <- Rsamtools::BamFile(file)

    scan_param <- Rsamtools::ScanBamParam(
        what = c("rname", "strand", "pos", "cigar", "seq", "flag"),
        tag  = c("MM", "ML")
    )

    reads <- tryCatch(
        Rsamtools::scanBam(bam_file, param = scan_param)[[1L]],
        error = function(e) {
            stop("Failed to read BAM file '", file, "': ", e$message)
        }
    )

    n_reads <- length(reads$pos)
    if (n_reads == 0L) {
        message("Note: BAM file '", file, "' contains no aligned reads.")
        return(.emptyModkitResult())
    }

    # ── Parse modifications from each read ────────────────────────────────────
    # Collect per-site modification calls as a list of data.frames
    site_records <- vector("list", n_reads)

    for (i in seq_len(n_reads)) {
        mm_tag <- reads$tag$MM[[i]]
        ml_tag <- reads$tag$ML[[i]]

        # Skip reads without modification tags
        if (is.null(mm_tag) || is.null(ml_tag) || nchar(mm_tag) == 0L) next

        pos_ref   <- reads$pos[[i]]    # 1-based leftmost mapping position
        cigar_str <- reads$cigar[[i]]
        seq_bases <- as.character(reads$seq[[i]])
        strand    <- as.character(reads$strand[[i]])
        chrom     <- as.character(reads$rname[[i]])
        flag      <- reads$flag[[i]]

        if (is.na(pos_ref) || is.na(cigar_str) || is.null(seq_bases)) next
        if (nchar(seq_bases) == 0L) next

        # Parse CIGAR → read-to-reference position map
        ref_positions <- .cigarToRefPos(cigar_str, pos_ref, seq_bases)
        if (is.null(ref_positions)) next

        # Parse MM tag → list of modification calls (read position + mod_type)
        mod_calls <- .parseMmTag(mm_tag, ml_tag, seq_bases)
        if (is.null(mod_calls) || nrow(mod_calls) == 0L) next

        # Map read positions to reference positions
        valid <- mod_calls$read_pos %in% seq_along(ref_positions)
        if (!any(valid)) next
        mod_calls <- mod_calls[valid, , drop = FALSE]

        ref_pos_for_mod <- ref_positions[mod_calls$read_pos]
        on_ref          <- !is.na(ref_pos_for_mod)
        if (!any(on_ref)) next

        site_records[[i]] <- data.frame(
            chrom     = chrom,
            position  = ref_pos_for_mod[on_ref],
            strand    = strand,
            mod_type  = mod_calls$mod_type[on_ref],
            is_mod    = mod_calls$is_mod[on_ref],
            stringsAsFactors = FALSE
        )
    }

    # ── Aggregate per site ────────────────────────────────────────────────────
    all_records <- do.call(rbind, site_records[!vapply(site_records, is.null, logical(1))])

    if (is.null(all_records) || nrow(all_records) == 0L) {
        message("Note: No modification records found in BAM '", file, "'.")
        return(.emptyModkitResult())
    }

    # Group by site key
    all_records$key <- paste0(
        all_records$chrom, ":",
        all_records$position, ":",
        all_records$strand, ":",
        all_records$mod_type
    )

    agg_list <- tapply(seq_len(nrow(all_records)), all_records$key, function(idx) {
        sub <- all_records[idx, , drop = FALSE]
        data.frame(
            chrom    = sub$chrom[[1L]],
            position = sub$position[[1L]],
            strand   = sub$strand[[1L]],
            mod_type = sub$mod_type[[1L]],
            coverage = nrow(sub),
            n_mod    = sum(sub$is_mod),
            stringsAsFactors = FALSE
        )
    }, simplify = FALSE)

    agg_df <- do.call(rbind, agg_list)
    agg_df$beta <- agg_df$n_mod / agg_df$coverage

    # ── Apply filters ─────────────────────────────────────────────────────────
    agg_df <- agg_df[agg_df$coverage >= min_coverage, , drop = FALSE]

    if (!is.null(mod_type)) {
        agg_df <- agg_df[agg_df$mod_type %in% mod_type, , drop = FALSE]
    }

    if (nrow(agg_df) == 0L) {
        return(.emptyModkitResult())
    }

    # ── Return standard format ────────────────────────────────────────────────
    data.frame(
        chrom    = as.character(agg_df$chrom),
        position = as.integer(agg_df$position),
        strand   = as.character(agg_df$strand),
        mod_type = as.character(agg_df$mod_type),
        beta     = as.numeric(agg_df$beta),
        coverage = as.integer(agg_df$coverage),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

# ─── Helper: CIGAR → reference position map ───────────────────────────────────

#' Map read positions to reference (genomic) positions via CIGAR string
#'
#' Returns a vector of length equal to the read sequence, where element \code{i}
#' is the 1-based reference position corresponding to read base \code{i}, or
#' \code{NA} for soft-clipped or inserted bases.
#'
#' @param cigar_str Character string. CIGAR string (e.g., \code{"5M2I3M"}).
#' @param ref_start Integer. 1-based leftmost reference position of the
#'   alignment.
#' @param seq_bases Character string of read sequence bases.
#'
#' @return Integer vector of length \code{nchar(seq_bases)}, or \code{NULL} on
#'   parse failure.
#'
#' @keywords internal
.cigarToRefPos <- function(cigar_str, ref_start, seq_bases) {
    read_len <- nchar(seq_bases)
    if (read_len == 0L || is.na(cigar_str)) return(NULL)

    # Parse CIGAR into (op, len) pairs
    matches <- gregexpr("([0-9]+)([MIDNSHP=X])", cigar_str, perl = TRUE)
    ops_raw <- regmatches(cigar_str, matches)[[1L]]
    if (length(ops_raw) == 0L) return(NULL)

    ops <- sub("^([0-9]+)([MIDNSHP=X])$", "\\2", ops_raw)
    lens <- as.integer(sub("^([0-9]+)([MIDNSHP=X])$", "\\1", ops_raw))

    # Map each read position to its reference position
    ref_pos_map <- integer(read_len)
    ref_pos_map[] <- NA_integer_

    ref_cur  <- ref_start  # current reference position (1-based)
    read_cur <- 1L         # current read position (1-based)

    for (k in seq_along(ops)) {
        op  <- ops[[k]]
        len <- lens[[k]]

        if (op %in% c("M", "=", "X")) {
            # Match/mismatch: advances both read and reference
            read_end <- read_cur + len - 1L
            if (read_end > read_len) read_end <- read_len
            n_use <- read_end - read_cur + 1L
            ref_pos_map[read_cur:read_end] <- seq(ref_cur, by = 1L, length.out = n_use)
            ref_cur  <- ref_cur  + len
            read_cur <- read_cur + len

        } else if (op == "I") {
            # Insertion: advances read only; ref positions remain NA
            read_cur <- read_cur + len

        } else if (op %in% c("D", "N")) {
            # Deletion/skip: advances reference only
            ref_cur <- ref_cur + len

        } else if (op == "S") {
            # Soft clip: read bases present but not aligned; NA ref pos
            read_cur <- read_cur + len

        } else if (op == "H") {
            # Hard clip: bases not in seq; skip
            next

        }
        # P (padding) skipped

        if (read_cur > read_len + 1L) break
    }

    ref_pos_map
}

# ─── Helper: Parse MM/ML tags ─────────────────────────────────────────────────

#' Parse MM and ML BAM tags into per-base modification calls
#'
#' Interprets the MM (base modification) and ML (modification likelihood) tags
#' according to the SAM specification and returns a data frame of read positions
#' with their modification type and whether they are called modified.
#'
#' @param mm_tag Character string. The MM tag value from the BAM record,
#'   e.g., \code{"A+a?,0,1,3;C+m?,5;"}.
#' @param ml_tag Integer vector. The ML tag values (0–255), parallel to the
#'   modifications listed in \code{mm_tag}.
#' @param seq_bases Character string. Read sequence.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{\code{read_pos}}{1-based position in the read.}
#'     \item{\code{mod_type}}{Modification type string (e.g., \code{"6mA"}).}
#'     \item{\code{is_mod}}{Logical; \code{TRUE} when ML probability > 0.5.}
#'   }
#'   Returns \code{NULL} on parse failure.
#'
#' @keywords internal
.parseMmTag <- function(mm_tag, ml_tag, seq_bases) {
    if (is.null(mm_tag) || is.null(ml_tag)) return(NULL)

    seq_vec <- strsplit(seq_bases, "")[[1L]]
    n_bases <- length(seq_vec)

    # Split MM tag on ";" — each block describes one modification type
    blocks <- strsplit(mm_tag, ";")[[1L]]
    blocks <- blocks[nzchar(blocks)]
    if (length(blocks) == 0L) return(NULL)

    result_list <- vector("list", length(blocks))
    ml_offset   <- 0L  # running index into ml_tag

    for (b in seq_along(blocks)) {
        block <- trimws(blocks[[b]])
        # Format: "<base_type>+<mod_code>[?,.]<delta>,<delta>,..."
        # e.g., "A+a?,0,1,3"  or  "C+m?,5"
        header_match <- regexpr("^([ACGT])([+-])([^,?. ]+)[?,.]?", block, perl = TRUE)
        if (header_match == -1L) next

        header_len  <- attr(header_match, "match.length")
        header_str  <- substr(block, 1L, header_len)
        rest        <- substr(block, header_len + 1L, nchar(block))

        # Extract components from header
        base_type <- substr(header_str, 1L, 1L)
        mod_code  <- sub("^[ACGT][+-]([^,?. ]+)[?,.]?$", "\\1", header_str)

        # Map mod_code → mod_type using the same map as modkit parser
        mt <- .MODKIT_CODE_MAP[[mod_code]]
        if (is.na(mt) || is.null(mt)) next  # unknown modification; skip

        # Parse delta-position array (comma-separated integers)
        if (nchar(rest) == 0L) next
        delta_str <- trimws(rest)
        if (nchar(delta_str) == 0L) next

        deltas <- tryCatch(
            as.integer(strsplit(delta_str, ",")[[1L]]),
            error = function(e) NULL
        )
        if (is.null(deltas) || length(deltas) == 0L) next

        n_mods_this_block <- length(deltas)

        # Locate modified bases in read sequence using delta-position encoding
        # Deltas are inter-base offsets: number of matching-type bases to skip
        # between consecutive modified positions
        base_positions <- which(seq_vec == base_type)
        if (length(base_positions) == 0L) {
            ml_offset <- ml_offset + n_mods_this_block
            next
        }

        read_positions <- integer(n_mods_this_block)
        base_cursor    <- 0L  # index into base_positions (0-based for cumsum)

        valid_block <- TRUE
        for (d in seq_len(n_mods_this_block)) {
            base_cursor <- base_cursor + deltas[[d]] + 1L
            if (base_cursor > length(base_positions)) {
                valid_block <- FALSE
                break
            }
            read_positions[[d]] <- base_positions[[base_cursor]]
        }

        if (!valid_block) {
            ml_offset <- ml_offset + n_mods_this_block
            next
        }

        # Retrieve corresponding ML probabilities
        ml_idx <- ml_offset + seq_len(n_mods_this_block)
        if (max(ml_idx) > length(ml_tag)) {
            ml_offset <- ml_offset + n_mods_this_block
            next
        }
        probs  <- as.integer(ml_tag[ml_idx]) / 255.0
        is_mod <- probs > 0.5

        result_list[[b]] <- data.frame(
            read_pos = read_positions,
            mod_type = mt,
            is_mod   = is_mod,
            stringsAsFactors = FALSE
        )

        ml_offset <- ml_offset + n_mods_this_block
    }

    result_list <- result_list[!vapply(result_list, is.null, logical(1))]
    if (length(result_list) == 0L) return(NULL)
    do.call(rbind, result_list)
}
