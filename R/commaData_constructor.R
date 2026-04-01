#' @importFrom methods new validObject
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
NULL

#' Create a commaData object from methylation calling output files
#'
#' Constructor for the \code{\link{commaData-class}} S4 class. Parses one or
#' more methylation calling output files (modkit, Megalodon, or Dorado), merges
#' them into a sites × samples matrix representation, and optionally loads
#' genomic annotation and motif site positions.
#'
#' @param files Named character vector mapping sample names to file paths.
#'   Names must match \code{colData$sample_name}. Example:
#'   \code{c(ctrl_1 = "/path/to/ctrl_1.bed", treat_1 = "/path/to/treat_1.bed")}.
#' @param colData A \code{data.frame} with one row per sample. Must contain
#'   columns \code{sample_name}, \code{condition}, and \code{replicate}.
#'   Additional columns (e.g., \code{file_path}, \code{batch}) are preserved.
#' @param genome Genome size information: a named integer vector of chromosome
#'   sizes (e.g., \code{c(NC_000913 = 4641652L)}), a path to a FASTA file, a
#'   \code{DNAStringSet} (Biostrings), or a \code{BSgenome} object. For
#'   single-chromosome genomes pass the \code{BSgenome} object directly or a
#'   named integer vector — do not index into the BSgenome with \code{$}
#'   (e.g., \code{BSgenome.Ecoli.NCBI.20080805$NC_000913}) as that yields a
#'   \code{DNAString} which has no chromosome name and cannot be used. Set to
#'   \code{NULL} to omit genome information (not recommended). When a
#'   multi-sequence source is provided, genomeInfo is automatically restricted
#'   to chromosomes present in the data.
#' @param annotation Optional. Path to a GFF3 or BED annotation file, or a
#'   pre-loaded \code{\link[GenomicRanges]{GRanges}} object. If \code{NULL},
#'   the annotation slot is left empty.
#' @param mod_type Optional character vector specifying which modification
#'   types to retain (e.g., \code{"6mA"} or \code{c("6mA", "5mC")}). If
#'   \code{NULL}, all modification types detected in the files are kept.
#' @param motif Optional character string. A DNA sequence motif (e.g.,
#'   \code{"GATC"}) to locate in the genome using \code{\link{findMotifSites}}.
#'   The results are stored in the \code{motifSites} slot as a genome-wide
#'   \code{GRanges} of all motif instances. Requires \code{genome} to be a
#'   FASTA path or \code{BSgenome} object (not a named integer vector). If
#'   \code{NULL}, the \code{motifSites} slot is left empty. \emph{Note:} this
#'   argument is distinct from \code{rowData(object)$motif}, which stores the
#'   per-site sequence context extracted automatically from the modkit
#'   \code{mod_code} field (e.g., \code{"a,GATC,1"} → \code{motif = "GATC"})
#'   and is \code{NA} for Dorado and Megalodon callers.
#' @param expected_mod_contexts Named list or \code{NULL}. If provided,
#'   specifies which modification type / sequence motif combinations to retain.
#'   Names must be modification type strings (e.g., \code{"6mA"}, \code{"5mC"}).
#'   Values are character vectors of motif strings (e.g., \code{c("GATC",
#'   "ACCACC")}). Sites whose \code{mod_context}
#'   (\code{paste(mod_type, motif, sep = "_")}) does not match any name–value
#'   pair are dropped before the object is assembled. A message is emitted
#'   reporting the number of sites dropped per modification type. Use
#'   \code{NULL} (default) to retain all sites.
#'   Example: \code{list("6mA" = "GATC", "5mC" = c("CCWGG", "CCGG"))}.
#'   \emph{Note:} for Dorado/Megalodon callers where \code{motif} is
#'   \code{NA}, the \code{mod_context} falls back to just \code{mod_type}
#'   (e.g., \code{"6mA"}), so those sites are only retained if you include
#'   \code{NA} in the motif vector for that type
#'   (e.g., \code{list("6mA" = NA)}).
#' @param min_coverage Integer. Minimum read depth to include a site. Sites
#'   present in a sample with coverage below this threshold have their beta
#'   value set to \code{NA}. Sites absent from a sample entirely are also
#'   \code{NA}. Default \code{5}.
#' @param caller Character string specifying the methylation caller that
#'   produced the input files. One of \code{"modkit"} (default),
#'   \code{"megalodon"}, or \code{"dorado"}.
#'
#' @return A valid \code{\link{commaData}} object.
#'
#' @details
#' The constructor uses a parse-then-merge strategy:
#' \enumerate{
#'   \item Each file is parsed independently using the appropriate parser.
#'   \item Sites are identified by a 5-part key:
#'     \code{"chrom:position:strand:mod_type:motif"} (motif is \code{"NA"}
#'     for Dorado and Megalodon callers).
#'   \item The union of all sites across all samples is taken.
#'   \item Beta values and coverage are arranged into sites × samples matrices,
#'     with \code{NA} for samples that do not cover a given site.
#'   \item Sites where coverage is below \code{min_coverage} in a sample have
#'     their beta value set to \code{NA} (but coverage is preserved).
#' }
#'
#' @examples
#' \dontrun{
#' # Load two modkit BED files
#' cd <- commaData(
#'   files   = c(
#'     ctrl_1  = "ctrl_1_modkit.bed",
#'     treat_1 = "treat_1_modkit.bed"
#'   ),
#'   colData = data.frame(
#'     sample_name = c("ctrl_1", "treat_1"),
#'     condition   = c("control", "treatment"),
#'     replicate   = c(1L, 1L)
#'   ),
#'   genome    = c(chr1 = 4641652L),
#'   annotation = "MG1655.gff3",
#'   caller    = "modkit"
#' )
#' cd
#' }
#'
#' @seealso \code{\link{commaData-class}}, \code{\link{methylation}},
#'   \code{\link[=coverage,commaData-method]{coverage}}, \code{\link{sampleInfo}}, \code{\link{siteInfo}},
#'   \code{\link{modTypes}}, \code{\link{loadAnnotation}},
#'   \code{\link{findMotifSites}}
#'
#' @export
commaData <- function(files,
                      colData,
                      genome               = NULL,
                      annotation           = NULL,
                      mod_type             = NULL,
                      motif                = NULL,
                      expected_mod_contexts = NULL,
                      min_coverage         = 5L,
                      caller               = "modkit") {

    min_coverage <- as.integer(min_coverage)
    caller       <- match.arg(caller, c("modkit", "megalodon", "dorado"))

    # ── Validate expected_mod_contexts ───────────────────────────────────────
    if (!is.null(expected_mod_contexts)) {
        if (!is.list(expected_mod_contexts) ||
                is.null(names(expected_mod_contexts)) ||
                any(names(expected_mod_contexts) == "")) {
            stop(
                "'expected_mod_contexts' must be a named list mapping modification ",
                "type strings to character vectors of motif strings ",
                "(e.g., list(\"6mA\" = \"GATC\", \"5mC\" = \"CCWGG\")) or NULL."
            )
        }
        bad_types <- setdiff(names(expected_mod_contexts), .VALID_MOD_TYPES)
        if (length(bad_types) > 0L) {
            stop(
                "Names in 'expected_mod_contexts' must be valid mod_type values. ",
                "Unrecognized: ", paste(bad_types, collapse = ", "),
                ". Allowed: ", paste(.VALID_MOD_TYPES, collapse = ", ")
            )
        }
    }

    # ── Validate colData ────────────────────────────────────────────────────
    if (!is.data.frame(colData)) {
        stop("colData must be a data.frame")
    }
    required_cd_cols <- c("sample_name", "condition", "replicate")
    missing_cd <- setdiff(required_cd_cols, colnames(colData))
    if (length(missing_cd) > 0) {
        stop(
            "colData is missing required columns: ",
            paste(missing_cd, collapse = ", ")
        )
    }

    # ── Validate files ──────────────────────────────────────────────────────
    if (!is.character(files) || is.null(names(files))) {
        stop(
            "files must be a named character vector mapping sample names to ",
            "file paths (e.g., c(ctrl_1 = '/path/to/ctrl_1.bed'))"
        )
    }
    missing_samples <- setdiff(names(files), colData$sample_name)
    if (length(missing_samples) > 0) {
        stop(
            "names(files) contains sample names not found in colData$sample_name: ",
            paste(missing_samples, collapse = ", ")
        )
    }
    extra_samples <- setdiff(colData$sample_name, names(files))
    if (length(extra_samples) > 0) {
        stop(
            "colData$sample_name contains samples with no file in 'files': ",
            paste(extra_samples, collapse = ", ")
        )
    }

    # ── Select parser ───────────────────────────────────────────────────────
    parser_fn <- switch(caller,
        modkit   = .parseModkit,
        megalodon = .parseMegalodon,
        dorado   = .parseDorado
    )

    # ── Parse each sample ───────────────────────────────────────────────────
    sample_names <- colData$sample_name
    parsed_list  <- vector("list", length(sample_names))
    names(parsed_list) <- sample_names

    for (sn in sample_names) {
        message("Parsing ", caller, " file for sample '", sn, "'...")
        parsed_list[[sn]] <- parser_fn(
            file        = files[sn],
            sample_name = sn,
            mod_type    = mod_type,
            min_coverage = 1L   # apply min_coverage AFTER merging (see below)
        )
    }

    # ── Build site universe ─────────────────────────────────────────────────
    all_sites <- unique(do.call(rbind, lapply(parsed_list, function(df) {
        if (nrow(df) == 0L) return(NULL)
        df[, c("chrom", "position", "strand", "mod_type", "motif"), drop = FALSE]
    })))

    if (is.null(all_sites) || nrow(all_sites) == 0L) {
        stop(
            "No methylation sites passed filtering across all samples. ",
            "Check your files, min_coverage, and mod_type arguments."
        )
    }

    # Stable sort: chrom, position, strand, mod_type, motif (NA last)
    ord      <- order(all_sites$chrom, all_sites$position,
                      all_sites$strand, all_sites$mod_type,
                      all_sites$motif)
    all_sites <- all_sites[ord, , drop = FALSE]
    rownames(all_sites) <- NULL

    # ── Compute mod_context ──────────────────────────────────────────────────
    # "6mA_GATC" when motif is known; falls back to "6mA" for NA-motif callers.
    all_sites$mod_context <- ifelse(
        is.na(all_sites$motif),
        all_sites$mod_type,
        paste(all_sites$mod_type, all_sites$motif, sep = "_")
    )

    # ── Apply expected_mod_contexts filter ───────────────────────────────────
    if (!is.null(expected_mod_contexts)) {
        # Build the set of allowed mod_context strings from the named list.
        # NA motif values in the list produce a fallback context (just mod_type).
        allowed_contexts <- character(0L)
        for (mt in names(expected_mod_contexts)) {
            motifs_for_mt <- expected_mod_contexts[[mt]]
            na_motifs  <- is.na(motifs_for_mt)
            str_motifs <- motifs_for_mt[!na_motifs]
            if (length(str_motifs) > 0L) {
                allowed_contexts <- c(allowed_contexts,
                                      paste(mt, str_motifs, sep = "_"))
            }
            if (any(na_motifs)) {
                allowed_contexts <- c(allowed_contexts, mt)
            }
        }
        allowed_contexts <- unique(allowed_contexts)

        drop_mask <- !(all_sites$mod_context %in% allowed_contexts)
        if (any(drop_mask)) {
            dropped <- all_sites[drop_mask, , drop = FALSE]
            for (mt in unique(dropped$mod_type)) {
                n_drop <- sum(dropped$mod_type == mt)
                message(
                    "expected_mod_contexts: dropping ", n_drop,
                    " site(s) with mod_type='", mt,
                    "' not in expected contexts."
                )
            }
            all_sites <- all_sites[!drop_mask, , drop = FALSE]
            rownames(all_sites) <- NULL
        }
        if (nrow(all_sites) == 0L) {
            stop(
                "No sites remain after applying 'expected_mod_contexts' filter. ",
                "Check that the named list matches the mod_type and motif values ",
                "present in your data."
            )
        }
    }

    site_keys <- paste(all_sites$chrom, all_sites$position,
                       all_sites$strand, all_sites$mod_type,
                       all_sites$motif, sep = ":")
    n_sites   <- length(site_keys)

    # ── Build matrices ──────────────────────────────────────────────────────
    n_samples    <- length(sample_names)
    methyl_mat   <- matrix(NA_real_,    nrow = n_sites, ncol = n_samples,
                           dimnames = list(site_keys, sample_names))
    coverage_mat <- matrix(NA_integer_, nrow = n_sites, ncol = n_samples,
                           dimnames = list(site_keys, sample_names))

    for (sn in sample_names) {
        df <- parsed_list[[sn]]
        if (nrow(df) == 0L) next

        df_keys <- paste(df$chrom, df$position, df$strand, df$mod_type, df$motif, sep = ":")
        idx      <- match(df_keys, site_keys)
        valid    <- !is.na(idx)

        methyl_mat[idx[valid], sn]   <- df$beta[valid]
        coverage_mat[idx[valid], sn] <- df$coverage[valid]
    }

    # ── Apply min_coverage: set beta NA where coverage < threshold ──────────
    below_threshold <- !is.na(coverage_mat) & coverage_mat < min_coverage
    methyl_mat[below_threshold] <- NA_real_

    # ── Build rowData ───────────────────────────────────────────────────────
    row_df <- S4Vectors::DataFrame(
        chrom       = all_sites$chrom,
        position    = all_sites$position,
        strand      = all_sites$strand,
        mod_type    = all_sites$mod_type,
        motif       = all_sites$motif,
        mod_context = all_sites$mod_context,
        row.names   = site_keys
    )

    # ── Build colData ───────────────────────────────────────────────────────
    # Reorder colData to match sample_names order in files
    cd_ordered <- as.data.frame(
        colData[match(sample_names, colData$sample_name), , drop = FALSE]
    )
    rownames(cd_ordered) <- cd_ordered$sample_name
    col_df <- S4Vectors::DataFrame(cd_ordered)

    # ── Genome info ─────────────────────────────────────────────────────────
    genome_info <- .validateGenomeInfo(genome)

    # Restrict genomeInfo to chromosomes actually present in the data
    if (!is.null(genome_info)) {
        data_chroms  <- unique(all_sites$chrom)
        extra_chroms <- setdiff(names(genome_info), data_chroms)
        if (length(extra_chroms) > 0L) {
            message(
                "Dropping ", length(extra_chroms), " chromosome(s) from genomeInfo ",
                "not present in data: ",
                paste(extra_chroms, collapse = ", ")
            )
            genome_info <- genome_info[names(genome_info) %in% data_chroms]
        }
    }

    # ── Annotation ──────────────────────────────────────────────────────────
    ann_gr <- GenomicRanges::GRanges()
    if (!is.null(annotation)) {
        if (is(annotation, "GRanges")) {
            ann_gr <- annotation
        } else if (is.character(annotation)) {
            ann_gr <- loadAnnotation(annotation)
        } else {
            stop("annotation must be a GRanges object or a file path string")
        }
    }

    # ── Motif sites ─────────────────────────────────────────────────────────
    motif_gr <- GenomicRanges::GRanges()
    if (!is.null(motif)) {
        if (is.null(genome) || is.integer(genome) || is.numeric(genome)) {
            warning(
                "motif specified but genome is a named integer vector (not a ",
                "FASTA/BSgenome). findMotifSites() requires sequence data. ",
                "Provide a FASTA path or BSgenome object to locate motif sites."
            )
        } else {
            message("Finding motif sites for '", motif, "'...")
            motif_gr <- findMotifSites(genome = genome, motif = motif)
        }
    }

    # ── Assemble SummarizedExperiment ────────────────────────────────────────
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = methyl_mat, coverage = coverage_mat),
        rowData = row_df,
        colData = col_df
    )

    # ── Construct commaData ─────────────────────────────────────────────────
    obj <- new("commaData",
               se,
               genomeInfo = genome_info,
               annotation = ann_gr,
               motifSites = motif_gr)

    validObject(obj)
    obj
}
