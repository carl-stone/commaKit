#' @importFrom methods new validObject
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo seqinfo
NULL

.validateMinCoverage <- function(min_coverage) {
  min_coverage_error <- "'min_coverage' must be a single positive integer."
  if (!is.numeric(min_coverage)) stop(min_coverage_error)
  if (length(min_coverage) != 1L) stop(min_coverage_error)
  if (is.na(min_coverage)) stop(min_coverage_error)
  if (!is.finite(min_coverage)) stop(min_coverage_error)
  if (min_coverage < 1L) stop(min_coverage_error)
  if (min_coverage != floor(min_coverage)) stop(min_coverage_error)
  if (min_coverage > .Machine$integer.max) stop(min_coverage_error)
  as.integer(min_coverage)
}

#' Create a commaData object from methylation calling output files
#'
#' Constructor for the \code{\link{commaData-class}} S4 class. Parses one or
#' more modkit pileup bedMethyl files, merges them into a sites × samples matrix
#' representation, and optionally loads genomic annotation and motif positions.
#'
#' @param files Named character vector mapping sample names to file paths.
#'   Names must match \code{colData$sample_name}. Example:
#'   \code{c(ctrl_1 = "/path/to/ctrl_1.bed", treat_1 = "/path/to/treat_1.bed")}.
#' @param colData A \code{data.frame} with one row per sample. Must contain
#'   columns \code{sample_name} and \code{replicate}. A \code{condition}
#'   column is optional and is used by design-aware functions such as
#'   \code{\link{diffMethyl}} when requested in their formula.
#'   Additional columns (e.g., \code{file_path}, \code{batch}) are preserved.
#' @param genome Genome size information: a named integer vector of chromosome
#'   sizes (e.g., \code{c(NC_000913 = 4641652L)}), a path to a FASTA file, a
#'   \code{DNAStringSet} (Biostrings), or a \code{BSgenome} object. For
#'   single-chromosome genomes pass the \code{BSgenome} object directly or a
#'   named integer vector — do not index into the BSgenome with \code{$}
#'   (e.g., \code{BSgenome.Ecoli.NCBI.20080805$NC_000913}) as that yields a
#'   \code{DNAString} which has no chromosome name and cannot be used. Set to
#'   \code{NULL} to omit genome information (not recommended). When a
#'   multi-sequence source is provided, genome info is automatically restricted
#'   to chromosomes present in the data. Chromosomes are assumed circular by
#'   default when genome information is attached, matching the package's
#'   bacterial genome default; edit \code{seqinfo(rowRanges(object))} after
#'   construction if a chromosome should be treated as linear.
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
#'   FASTA path, \code{DNAStringSet}, or \code{BSgenome} object (not a named
#'   integer vector). If \code{NULL}, the \code{motifSites} slot is left empty.
#'   \emph{Note:} this
#'   argument is distinct from \code{rowData(object)$motif}, which stores the
#'   per-site sequence context extracted automatically from the modkit
#'   \code{mod_code} field (e.g., \code{"a,GATC,1"} → \code{motif = "GATC"}).
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
#' @param min_coverage Integer. Minimum read depth to include a site. Sites
#'   present in a sample with coverage below this threshold have their beta
#'   value set to \code{NA}. Sites absent from a sample entirely are also
#'   \code{NA}. Default \code{5}.
#' @return A valid \code{\link{commaData}} object.
#'
#' @details
#' The constructor uses a parse-then-merge strategy:
#' \enumerate{
#'   \item Each file is parsed independently using the appropriate parser.
#'   \item Sites are identified by their genomic position (chromosome, position,
#'     strand) plus modification type and motif context.
#'   \item The union of all sites across all samples is taken, using
#'     \code{findOverlaps()} for alignment.
#'   \item Beta values and coverage are arranged into sites \eqn{\times} samples
#'     matrices, with \code{NA} for samples that do not cover a given site.
#'   \item Observed modified, canonical, and non-target modified read counts
#'     are preserved as \code{mod_counts}, \code{canonical_counts}, and
#'     \code{other_mod_counts} assays.
#'   \item Assay-layer provenance and default roles are recorded in
#'     \code{metadata(object)$assay_provenance} and
#'     \code{metadata(object)$assay_defaults}.
#'   \item Sites where coverage is below \code{min_coverage} in a sample have
#'     their beta value set to \code{NA} (but coverage is preserved).
#' }
#'
#' @examples
#' # Construct a commaData object from built-in example data
#' data(comma_example_data)
#' comma_example_data
#'
#' \dontrun{
#' # Load two modkit BED files (requires user-provided files)
#' cd <- commaData(
#'   files = c(
#'     ctrl_1  = "ctrl_1_modkit.bed",
#'     treat_1 = "treat_1_modkit.bed"
#'   ),
#'   colData = data.frame(
#'     sample_name = c("ctrl_1", "treat_1"),
#'     condition   = c("control", "treatment"),
#'     replicate   = c(1L, 1L)
#'   ),
#'   genome = c(chr1 = 4641652L),
#'   annotation = "MG1655.gff3"
#' )
#' cd
#' }
#'
#' @seealso \code{\link{commaData-class}}, \code{\link{methylation}},
#'   \code{\link{siteCoverage}}, \code{\link{modCounts}},
#'   \code{\link{canonicalCounts}}, \code{\link{assayLayers}},
#'   \code{\link{assayProvenance}},
#'   \code{\link{sampleInfo}}, \code{\link{siteInfo}},
#'   \code{\link{modTypes}}, \code{\link{loadAnnotation}},
#'   \code{\link{findMotifSites}}
#'
#' @export
commaData <- function(files,
                      colData,
                      genome = NULL,
                      annotation = NULL,
                      mod_type = NULL,
                      motif = NULL,
                      expected_mod_contexts = NULL,
                      min_coverage = 5L) {
  min_coverage <- .validateMinCoverage(min_coverage)

  # ── Validate expected_mod_contexts ───────────────────────────────────────
  if (!is.null(expected_mod_contexts)) {
    invalid_expected_mod_contexts <- !is.list(expected_mod_contexts)
    invalid_expected_mod_contexts <- invalid_expected_mod_contexts ||
      is.null(names(expected_mod_contexts))
    invalid_expected_mod_contexts <- invalid_expected_mod_contexts ||
      any(names(expected_mod_contexts) == "")
    if (invalid_expected_mod_contexts) {
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
  required_cd_cols <- c("sample_name", "replicate")
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

  # ── Parse each sample ───────────────────────────────────────────────────
  sample_names <- colData$sample_name
  parsed_list <- vector("list", length(sample_names))
  names(parsed_list) <- sample_names

  for (sn in sample_names) {
    message("Parsing modkit bedMethyl file for sample '", sn, "'...")
    parsed_list[[sn]] <- .parseModkit(
      file = files[sn],
      sample_name = sn,
      mod_type = mod_type,
      min_coverage = 1L # apply min_coverage AFTER merging (see below)
    )
  }

  # ── Build site universe ─────────────────────────────────────────────────
  all_sites <- unique(do.call(rbind, lapply(parsed_list, function(df) {
    if (nrow(df) == 0L) {
      return(NULL)
    }
    df[, c("chrom", "position", "strand", "mod_type", "motif"), drop = FALSE]
  })))

  if (is.null(all_sites) || nrow(all_sites) == 0L) {
    stop(
      "No methylation sites passed filtering across all samples. ",
      "Check your files, min_coverage, and mod_type arguments."
    )
  }

  # Stable sort: chrom, position, strand, mod_type, motif (NA last)
  ord <- order(
    all_sites$chrom, all_sites$position,
    all_sites$strand, all_sites$mod_type,
    all_sites$motif
  )
  all_sites <- all_sites[ord, , drop = FALSE]
  rownames(all_sites) <- NULL

  # ── Compute mod_context on the fly for filtering ────────────────────────
  # "6mA_GATC" when motif is known; falls back to "6mA" when motif is NA.
  site_ctx <- ifelse(
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
      na_motifs <- is.na(motifs_for_mt)
      str_motifs <- motifs_for_mt[!na_motifs]
      if (length(str_motifs) > 0L) {
        allowed_contexts <- c(
          allowed_contexts,
          paste(mt, str_motifs, sep = "_")
        )
      }
      if (any(na_motifs)) {
        allowed_contexts <- c(allowed_contexts, mt)
      }
    }
    allowed_contexts <- unique(allowed_contexts)

    drop_mask <- !(site_ctx %in% allowed_contexts)
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

  n_sites <- nrow(all_sites)

  # ── Build matrices ──────────────────────────────────────────────────────
  n_samples <- length(sample_names)
  methyl_mat <- matrix(NA_real_,
    nrow = n_sites, ncol = n_samples,
    dimnames = list(NULL, sample_names)
  )
  coverage_mat <- matrix(NA_integer_,
    nrow = n_sites, ncol = n_samples,
    dimnames = list(NULL, sample_names)
  )
  mod_counts_mat <- matrix(NA_integer_,
    nrow = n_sites, ncol = n_samples,
    dimnames = list(NULL, sample_names)
  )
  canonical_counts_mat <- matrix(NA_integer_,
    nrow = n_sites, ncol = n_samples,
    dimnames = list(NULL, sample_names)
  )
  other_mod_counts_mat <- matrix(NA_integer_,
    nrow = n_sites, ncol = n_samples,
    dimnames = list(NULL, sample_names)
  )

  # ── Build rowRanges (GRanges) for findOverlaps merge ────────────────────
  site_gr <- GenomicRanges::GRanges(
    seqnames = all_sites$chrom,
    ranges = IRanges::IRanges(start = all_sites$position, width = 1L),
    strand = all_sites$strand,
    mod_type = factor(all_sites$mod_type, levels = .VALID_MOD_TYPES),
    motif = all_sites$motif
  )

  for (sn in sample_names) {
    df <- parsed_list[[sn]]
    if (nrow(df) == 0L) next

    # Parsed rows must be unique for the site identity used by the merge.
    # Otherwise later matrix assignment would silently let the last duplicate
    # win.
    identity_cols <- c("chrom", "position", "strand", "mod_type", "motif")
    df_identity <- df[, identity_cols, drop = FALSE]
    if (anyDuplicated(df_identity)) {
      dup_mask <- duplicated(df_identity) |
        duplicated(df_identity, fromLast = TRUE)
      dup <- df_identity[dup_mask, , drop = FALSE]
      stop(
        "Parser returned duplicate methylation site rows for sample '", sn,
        "'. Duplicate chrom/position/strand/mod_type/motif entries ",
        "must be aggregated before commaData() can merge samples. ",
        "First duplicate: ", paste(dup[1L, ], collapse = ":")
      )
    }

    # Align parsed sites to the site universe using findOverlaps()
    df_gr <- GenomicRanges::GRanges(
      seqnames = df$chrom,
      ranges   = IRanges::IRanges(start = df$position, width = 1L),
      strand   = df$strand
    )
    GenomicRanges::mcols(df_gr)$mod_type <- df$mod_type
    GenomicRanges::mcols(df_gr)$motif <- df$motif

    hits <- GenomicRanges::findOverlaps(df_gr, site_gr, type = "equal")
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)

    # Verify mod_type and motif match (same position can have multiple types)
    mc_q <- GenomicRanges::mcols(df_gr)[qh, c("mod_type", "motif")]
    mc_s <- GenomicRanges::mcols(site_gr)[sh, c("mod_type", "motif")]
    type_match <- mc_q$mod_type == mc_s$mod_type
    both_missing_motif <- is.na(mc_q$motif) & is.na(mc_s$motif)
    both_known_matching_motif <- !is.na(mc_q$motif) &
      !is.na(mc_s$motif) &
      mc_q$motif == mc_s$motif
    motif_match <- both_missing_motif | both_known_matching_motif
    valid <- type_match & motif_match

    if (anyDuplicated(qh[valid])) {
      stop(
        "Ambiguous site merge for sample '", sn,
        "': at least one parsed row matched multiple site-universe rows."
      )
    }

    idx <- rep(NA_integer_, nrow(df))
    idx[qh[valid]] <- sh[valid]

    valid_idx <- !is.na(idx)
    methyl_mat[idx[valid_idx], sn] <- df$beta[valid_idx]
    coverage_mat[idx[valid_idx], sn] <- df$coverage[valid_idx]
    if ("mod_counts" %in% colnames(df)) {
      mod_counts_mat[idx[valid_idx], sn] <- df$mod_counts[valid_idx]
    }
    if ("canonical_counts" %in% colnames(df)) {
      canonical_counts_mat[idx[valid_idx], sn] <- df$canonical_counts[valid_idx]
    }
    if ("other_mod_counts" %in% colnames(df)) {
      other_mod_counts_mat[idx[valid_idx], sn] <- df$other_mod_counts[valid_idx]
    }
  }

  # ── Apply min_coverage: set beta NA where coverage < threshold ──────────
  below_threshold <- !is.na(coverage_mat) & coverage_mat < min_coverage
  methyl_mat[below_threshold] <- NA_real_

  # ── Genome info ─────────────────────────────────────────────────────────
  genome_info <- .validateGenomeInfo(genome)

  # Restrict genome info to chromosomes actually present in the data, while
  # failing clearly when methylation files contain chromosomes absent from the
  # supplied genome. A later Seqinfo assignment would fail too, but with a
  # lower-level message that is hard to connect to the input files.
  if (!is.null(genome_info)) {
    data_chroms <- unique(all_sites$chrom)
    missing_chroms <- setdiff(data_chroms, names(genome_info))
    if (length(missing_chroms) > 0L) {
      stop(
        "genome is missing chromosome(s) present in methylation data: ",
        paste(missing_chroms, collapse = ", "), ". ",
        "Chromosome names in 'genome' must match the imported files."
      )
    }

    extra_chroms <- setdiff(names(genome_info), data_chroms)
    if (length(extra_chroms) > 0L) {
      message(
        "Dropping ", length(extra_chroms), " chromosome(s) from genome info ",
        "not present in data: ",
        paste(extra_chroms, collapse = ", ")
      )
      genome_info <- genome_info[names(genome_info) %in% data_chroms]
    }
  }

  # ── Attach Seqinfo to rowRanges ──────────────────────────────────────
  if (!is.null(genome_info)) {
    GenomeInfoDb::seqinfo(site_gr) <- .makeSeqinfo(genome_info)
  }

  # ── Build colData ───────────────────────────────────────────────────────
  # Reorder colData to match sample_names order in files
  cd_ordered <- as.data.frame(
    colData[match(sample_names, colData$sample_name), , drop = FALSE]
  )
  rownames(cd_ordered) <- cd_ordered$sample_name
  col_df <- S4Vectors::DataFrame(cd_ordered)

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

  # ── Assemble RangedSummarizedExperiment ─────────────────────────────────
  rse <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      methylation = methyl_mat,
      coverage = coverage_mat,
      mod_counts = mod_counts_mat,
      canonical_counts = canonical_counts_mat,
      other_mod_counts = other_mod_counts_mat
    ),
    rowRanges = site_gr,
    colData = col_df
  )

  # ── Construct commaData ─────────────────────────────────────────────────
  obj <- new("commaData", rse)

  # Store annotation and motifSites in metadata (not as S4 slots)
  S4Vectors::metadata(obj)$annotation <- ann_gr
  S4Vectors::metadata(obj)$motifSites <- motif_gr

  # Store construction threshold in metadata for reproducibility
  S4Vectors::metadata(obj)$min_coverage <- min_coverage
  count_provenance <- "observed_modkit_pileup"
  S4Vectors::metadata(obj)$assay_defaults <- list(
    methylation = "methylation",
    coverage = "coverage",
    mod_counts = "mod_counts",
    canonical_counts = "canonical_counts",
    other_mod_counts = "other_mod_counts"
  )
  S4Vectors::metadata(obj)$assay_provenance <- list(
    methylation = .makeAssayLayerRecord(
      type = "filtered_beta",
      source = "modkit",
      role = "methylation",
      parent_assays = c("coverage"),
      method = "modkit_beta_filter",
      params = list(min_coverage = min_coverage, filtered_assay = TRUE),
      default_for = "methylation"
    ),
    coverage = .makeAssayLayerRecord(
      type = "observed_total_coverage",
      source = "modkit",
      role = "coverage",
      method = "modkit_coverage",
      default_for = "coverage"
    ),
    mod_counts = .makeAssayLayerRecord(
      type = "observed_counts",
      source = count_provenance,
      role = "mod_counts",
      parent_assays = "coverage",
      method = count_provenance,
      default_for = "mod_counts"
    ),
    canonical_counts = .makeAssayLayerRecord(
      type = "observed_counts",
      source = count_provenance,
      role = "canonical_counts",
      parent_assays = "coverage",
      method = count_provenance,
      default_for = "canonical_counts"
    ),
    other_mod_counts = .makeAssayLayerRecord(
      type = "observed_counts",
      source = count_provenance,
      role = "other_mod_counts",
      parent_assays = "coverage",
      method = count_provenance,
      default_for = "other_mod_counts"
    )
  )

  validObject(obj)
  obj
}
