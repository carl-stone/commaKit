#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
NULL

# ─── Constants ────────────────────────────────────────────────────────────────

# Feature types that are large regions; overlap_only defaults to TRUE for these.
.COMMA_REGION_FEATURE_TYPES <- c(
    "gene", "CDS", "mRNA", "tRNA", "rRNA", "ncRNA",
    "operon", "repeat_region", "prophage", "region", "insertion_sequence"
)

# Map EcoCyc sigma factor names to E. coli gene symbols.
# Users can override via custom TERM2GENE for other organisms.
.SIGMA_FACTOR_GENE_MAP <- c(
    Sigma70 = "rpoD", Sigma24 = "rpoE", Sigma32 = "rpoH",
    Sigma28 = "fliA", Sigma38 = "rpoS", Sigma54 = "rpoN",
    Sigma19 = "fecI"
)

# ─── Existing internal helpers (kept for backward compatibility) ───────────────

# Map differential methylation results from sites to genes.
#
# Explodes the CharacterList/list annotation column so each site × gene
# combination becomes its own row.  Intergenic sites (length-0 entries) are
# silently excluded.
#
# @param res_df  data.frame returned by results()
# @param gene_col  name of the list-valued column containing gene IDs
# @return data.frame with columns gene_id, site_key, dm_padj, dm_delta_beta
.siteToGeneMap <- function(res_df, gene_col) {
    if (!gene_col %in% colnames(res_df)) {
        stop(
            "Column '", gene_col, "' not found in results.\n",
            "Run annotateSites() before enrichMethylation(), e.g.:\n",
            "  object <- annotateSites(object, keep = \"overlap\")"
        )
    }

    gene_lists <- res_df[[gene_col]]
    if (is(gene_lists, "CharacterList")) {
        gene_lists <- as.list(gene_lists)
    }

    lens <- lengths(gene_lists)
    has_gene <- lens > 0L

    if (!any(has_gene)) {
        warning(
            "No sites have gene annotations in column '", gene_col, "'. ",
            "Returning empty gene map."
        )
        return(data.frame(
            gene_id       = character(),
            site_key      = character(),
            dm_padj       = numeric(),
            dm_delta_beta = numeric(),
            stringsAsFactors = FALSE
        ))
    }

    res_sub        <- res_df[has_gene, , drop = FALSE]
    gene_lists_sub <- gene_lists[has_gene]
    lens_sub       <- lens[has_gene]
    site_keys      <- paste(res_sub$chrom, res_sub$position,
                            res_sub$strand, sep = ":")

    result <- data.frame(
        gene_id       = unlist(gene_lists_sub, use.names = FALSE),
        site_key      = rep(site_keys, lens_sub),
        dm_padj       = rep(res_sub$dm_padj, lens_sub),
        dm_delta_beta = rep(res_sub$dm_delta_beta, lens_sub),
        stringsAsFactors = FALSE
    )
    result[!is.na(result$gene_id), , drop = FALSE]
}


# Compute a per-gene methylation score for GSEA ranking.
#
# @param site_gene_df  data.frame from .siteToGeneMap() or .extractGeneRoles()
# @param score_metric  "combined" | "padj" | "delta_beta"
# @param agg           how to aggregate site scores per gene: "max" | "mean"
# @return  named numeric vector sorted decreasing, suitable for clusterProfiler::GSEA()
.computeGeneScores <- function(site_gene_df, score_metric = "combined",
                               agg = "max") {
    valid <- !is.na(site_gene_df$dm_padj) & !is.na(site_gene_df$dm_delta_beta)
    df <- site_gene_df[valid, , drop = FALSE]

    if (nrow(df) == 0L) {
        return(setNames(numeric(0L), character(0L)))
    }

    eps <- .Machine$double.eps
    df$site_score <- switch(
        score_metric,
        combined   = -log10(pmax(df$dm_padj, eps)) * sign(df$dm_delta_beta),
        padj       = -log10(pmax(df$dm_padj, eps)),
        delta_beta = df$dm_delta_beta,
        stop("'score_metric' must be one of: \"combined\", \"padj\", \"delta_beta\"")
    )

    genes <- unique(df$gene_id)
    genes <- genes[!is.na(genes)]

    if (length(genes) == 0L) {
        return(setNames(numeric(0L), character(0L)))
    }

    gene_scores <- vapply(genes, function(g) {
        s <- df$site_score[df$gene_id == g]
        s <- s[!is.na(s)]
        if (length(s) == 0L) return(NA_real_)
        if (agg == "max") s[which.max(abs(s))] else mean(s)
    }, FUN.VALUE = numeric(1L))
    names(gene_scores) <- genes
    gene_scores <- gene_scores[!is.na(gene_scores)]

    sort(gene_scores, decreasing = TRUE)
}

# ─── New internal helpers for regulatory role extraction ──────────────────────

# Parse target gene names from feature names based on feature type.
#
# Always returns a list of character vectors (one element per input name,
# each element may have length > 1 for operonic names).
# NA indicates the name could not be parsed.
#
# @param feature_names  character vector of raw feature names
# @param feature_type   single string
# @param tu_values      optional character vector: transcription_unit GFF3
#                       attribute values, parallel to feature_names (preferred
#                       source for protein_binding_site and RNA_binding_site)
# @return list of character vectors
.parseTargetGenes <- function(feature_names, feature_type, tu_values = NULL) {
    n <- length(feature_names)

    switch(feature_type,
        gene = , CDS = , tRNA = , rRNA = , ncRNA = {
            # Identity: feature name is already the gene symbol
            as.list(feature_names)
        },

        promoter = ,
        minus_10_signal = ,
        minus_35_signal = ,
        transcription_factor_binding_site = {
            # Strip promoter number suffix (p, p1, p2, ...) then split on "-"
            # e.g. "geneA-geneBp1" -> c("geneA","geneB"); "aaeRp" -> "aaeR"
            stripped <- sub("p\\d*$", "", feature_names, perl = TRUE)
            strsplit(stripped, "-", fixed = TRUE)
        },

        protein_binding_site = {
            result <- vector("list", n)
            for (i in seq_len(n)) {
                nm <- feature_names[i]
                tu <- if (!is.null(tu_values)) tu_values[i] else NA_character_

                # Preferred: transcription_unit attribute (explicit EcoCyc field)
                if (!is.na(tu) && nzchar(tu)) {
                    # Split TUs by " // " separator; use each TU name as target
                    tus <- trimws(strsplit(tu, " // ", fixed = TRUE)[[1]])
                    result[[i]] <- tus[nzchar(tus)]
                } else if (grepl("\\bof\\s+\\S+$", nm, perl = TRUE)) {
                    # Fallback A: "... of {promoter_name}" - parse as promoter
                    prom <- sub(".*\\bof\\s+(\\S+)$", "\\1", nm, perl = TRUE)
                    stripped <- sub("p\\d*$", "", prom, perl = TRUE)
                    result[[i]] <- strsplit(stripped, "-", fixed = TRUE)[[1]]
                } else {
                    # Fallback B: no target info available
                    result[[i]] <- NA_character_
                }
            }
            result
        },

        RNA_binding_site = {
            result <- vector("list", n)
            for (i in seq_len(n)) {
                nm <- feature_names[i]
                tu <- if (!is.null(tu_values)) tu_values[i] else NA_character_

                # Preferred: transcription_unit attribute
                if (!is.na(tu) && nzchar(tu)) {
                    tus <- trimws(strsplit(tu, " // ", fixed = TRUE)[[1]])
                    result[[i]] <- tus[nzchar(tus)]
                } else if (grepl("\\bregulating\\s+", nm, perl = TRUE)) {
                    # Fallback: "... regulating {gene}" - split on "-" for operons
                    target_str <- sub(".*\\bregulating\\s+", "", nm, perl = TRUE)
                    result[[i]] <- strsplit(target_str, "-", fixed = TRUE)[[1]]
                } else {
                    result[[i]] <- NA_character_
                }
            }
            result
        },

        terminator = {
            # "{gene} terminator" -> "{gene}"
            stripped <- sub("\\s+terminator$", "", feature_names, perl = TRUE)
            as.list(stripped)
        },

        {
            # Default: identity (insertion_sequence, repeat_region, etc.)
            as.list(feature_names)
        }
    )
}

# Parse regulator gene names from feature names + optional subtype values.
#
# Always returns a list of character vectors. NA for types with no regulator
# or where the regulator cannot be parsed.
#
# @param feature_names   character vector of raw feature names
# @param feature_type    single string
# @param subtype_values  optional character vector: feature_subtype column
#                        (sigma factor identity for TFBSs)
# @return list of character vectors
.parseRegulatorGenes <- function(feature_names, feature_type,
                                  subtype_values = NULL) {
    n <- length(feature_names)
    na_list <- as.list(rep(NA_character_, n))

    switch(feature_type,
        transcription_factor_binding_site = ,
        minus_10_signal = ,
        minus_35_signal = {
            # Regulator: sigma factor gene from feature_subtype_values
            if (is.null(subtype_values)) {
                return(na_list)
            }
            result <- vector("list", n)
            for (i in seq_len(n)) {
                sigma <- subtype_values[i]
                gene  <- .SIGMA_FACTOR_GENE_MAP[sigma]
                result[[i]] <- if (is.na(gene)) NA_character_ else unname(gene)
            }
            result
        },

        protein_binding_site = {
            # Extract binding protein name before " binding site" or
            # " DNA-binding-site"; lowercase first character only
            result <- vector("list", n)
            pat <- "^([A-Za-z][A-Za-z0-9]+)[- ](binding|DNA-binding)"
            for (i in seq_len(n)) {
                nm <- feature_names[i]
                m  <- regmatches(nm, regexpr(pat, nm, perl = TRUE))
                if (length(m) > 0L && nzchar(m)) {
                    protein <- sub("[- ](binding|DNA-binding).*", "", m, perl = TRUE)
                    gene <- paste0(tolower(substring(protein, 1, 1)),
                                   substring(protein, 2))
                    result[[i]] <- gene
                } else {
                    result[[i]] <- NA_character_
                }
            }
            result
        },

        RNA_binding_site = {
            # Extract first word; if it matches a bacterial gene-name pattern,
            # lowercase first char and use as regulator gene.
            # Gene name pattern: 3-6 chars, first letter uppercase, rest
            # lowercase with optional trailing digit/uppercase (e.g. ArcZ,
            # AcnB, DsrA, RyhB). Riboswitches (e.g. "adenosylcobalamin") won't
            # match because they're too long or all-lowercase.
            result <- vector("list", n)
            gene_pat <- "^[A-Z][a-z]{2,5}[0-9A-Z]?$"
            for (i in seq_len(n)) {
                nm         <- feature_names[i]
                first_word <- sub("^(\\S+).*", "\\1", nm, perl = TRUE)
                if (grepl(gene_pat, first_word, perl = TRUE)) {
                    gene <- paste0(tolower(substring(first_word, 1, 1)),
                                   substring(first_word, 2))
                    result[[i]] <- gene
                } else {
                    result[[i]] <- NA_character_
                }
            }
            result
        },

        {
            # No regulator for this feature type (gene, terminator, IS, etc.)
            na_list
        }
    )
}

# Build a site-gene role map for a single feature type.
#
# Filters sites to those with the given feature type, applies the overlap_only
# filter if requested, parses feature names into target and/or regulator genes,
# and returns a data.frame with columns:
#   gene_id, role, role_type, site_key, dm_padj, dm_delta_beta
#
# @param res_df           data.frame from results()
# @param ft               single feature type string (or NULL for all)
# @param gene_col         column in res_df with feature_names CharacterList
# @param overlap_only     logical: restrict to rel_position==0 hits
# @param rel_position_col name of the rel_position column
# @return data.frame or NULL if no sites match
.extractGeneRoles <- function(res_df, ft, gene_col,
                               overlap_only = FALSE,
                               rel_position_col = "rel_position") {
    if (!gene_col %in% colnames(res_df)) {
        stop("Column '", gene_col, "' not found. Run annotateSites() first.")
    }

    ft_col <- "feature_types"
    if (!ft_col %in% colnames(res_df)) {
        stop("Column 'feature_types' not found. Run annotateSites() first.")
    }

    # ── Step 1: Find per-site indices matching ft ──────────────────────────────
    ft_lists <- if (is(res_df[[ft_col]], "CharacterList")) {
        as.list(res_df[[ft_col]])
    } else {
        res_df[[ft_col]]
    }

    type_idx <- lapply(ft_lists, function(x) {
        if (length(x) == 0L || all(is.na(x))) integer(0L) else which(x == ft)
    })
    has_match <- lengths(type_idx) > 0L

    if (!any(has_match)) return(NULL)

    res_sub      <- res_df[has_match, , drop = FALSE]
    type_idx_sub <- type_idx[has_match]

    # ── Step 2: Optional overlap_only filter ──────────────────────────────────
    if (isTRUE(overlap_only) && rel_position_col %in% colnames(res_df)) {
        rp_all <- if (is(res_sub[[rel_position_col]], "IntegerList")) {
            as.list(res_sub[[rel_position_col]])
        } else {
            res_sub[[rel_position_col]]
        }
        # Keep only ft-indices where rel_position == 0
        type_idx_sub <- mapply(function(tidx, rp) {
            if (length(tidx) == 0L) return(integer(0L))
            inside_ft <- vapply(tidx, function(j) {
                length(rp) >= j && !is.na(rp[j]) && rp[j] == 0L
            }, logical(1L))
            tidx[inside_ft]
        }, type_idx_sub, rp_all, SIMPLIFY = FALSE)

        still_match <- lengths(type_idx_sub) > 0L
        if (!any(still_match)) return(NULL)
        res_sub      <- res_sub[still_match, , drop = FALSE]
        type_idx_sub <- type_idx_sub[still_match]
    }

    # ── Step 3: Extract raw feature names for this ft ─────────────────────────
    gn_lists <- if (is(res_sub[[gene_col]], "CharacterList")) {
        as.list(res_sub[[gene_col]])
    } else {
        res_sub[[gene_col]]
    }
    raw_per_site <- mapply(function(nms, idx) as.character(nms[idx]),
                           gn_lists, type_idx_sub, SIMPLIFY = FALSE)

    # ── Step 4: Extract optional metadata columns ──────────────────────────────
    sub_col     <- "feature_subtype_values"
    tu_col      <- "transcription_unit_values"
    has_subtype <- sub_col %in% colnames(res_sub)
    has_tu      <- tu_col  %in% colnames(res_sub)

    subtype_per_site <- if (has_subtype) {
        st_lists <- if (is(res_sub[[sub_col]], "CharacterList")) {
            as.list(res_sub[[sub_col]])
        } else {
            res_sub[[sub_col]]
        }
        mapply(function(st, idx) as.character(st[idx]),
               st_lists, type_idx_sub, SIMPLIFY = FALSE)
    } else NULL

    tu_per_site <- if (has_tu) {
        tu_lists <- if (is(res_sub[[tu_col]], "CharacterList")) {
            as.list(res_sub[[tu_col]])
        } else {
            res_sub[[tu_col]]
        }
        mapply(function(tu, idx) as.character(tu[idx]),
               tu_lists, type_idx_sub, SIMPLIFY = FALSE)
    } else NULL

    # ── Step 5: Parse target and regulator genes ───────────────────────────────
    site_keys <- paste(res_sub$chrom, res_sub$position, res_sub$strand, sep = ":")
    n_sites   <- nrow(res_sub)

    rows <- vector("list", n_sites)

    for (i in seq_len(n_sites)) {
        nms      <- raw_per_site[[i]]
        subtypes <- if (!is.null(subtype_per_site)) subtype_per_site[[i]] else NULL
        tus      <- if (!is.null(tu_per_site))      tu_per_site[[i]]      else NULL
        sk       <- site_keys[i]
        padj     <- res_sub$dm_padj[i]
        db       <- res_sub$dm_delta_beta[i]

        if (length(nms) == 0L) {
            rows[[i]] <- NULL
            next
        }

        # Target genes (one entry per name, possibly expanded by operonic split)
        tgt_lists  <- .parseTargetGenes(nms, ft, tus)
        tgt_genes  <- unlist(tgt_lists, use.names = FALSE)
        tgt_genes  <- tgt_genes[!is.na(tgt_genes) & nzchar(tgt_genes)]

        # Regulator genes
        reg_lists  <- .parseRegulatorGenes(nms, ft, subtypes)
        reg_genes  <- unlist(reg_lists, use.names = FALSE)
        reg_genes  <- reg_genes[!is.na(reg_genes) & nzchar(reg_genes)]
        reg_genes  <- unique(reg_genes)  # deduplicate (e.g. autoregulatory sites)

        # role_type for this feature_type
        reg_type <- switch(ft,
            transcription_factor_binding_site = ,
            minus_10_signal = ,
            minus_35_signal = "sigma_factor",
            protein_binding_site = "TF_protein",
            RNA_binding_site     = "RNA_regulator",
            NA_character_
        )

        site_rows <- list()
        if (length(tgt_genes) > 0L) {
            tgt_df <- data.frame(
                gene_id       = unique(tgt_genes),
                role          = "target",
                role_type     = NA_character_,
                site_key      = sk,
                dm_padj       = padj,
                dm_delta_beta = db,
                stringsAsFactors = FALSE
            )
            site_rows <- c(site_rows, list(tgt_df))
        }
        if (length(reg_genes) > 0L) {
            reg_df <- data.frame(
                gene_id       = unique(reg_genes),
                role          = "regulator",
                role_type     = reg_type,
                site_key      = sk,
                dm_padj       = padj,
                dm_delta_beta = db,
                stringsAsFactors = FALSE
            )
            site_rows <- c(site_rows, list(reg_df))
        }
        rows[[i]] <- if (length(site_rows) > 0L) do.call(rbind, site_rows) else NULL
    }

    combined <- do.call(rbind, rows)
    if (is.null(combined) || nrow(combined) == 0L) return(NULL)
    combined
}

# Run ORA and/or GSEA for a prepared site-gene data.frame.
#
# This is the refactored ORA+GSEA block from enrichMethylation(), now callable
# once per feature_type × gene_role combination.
#
# @param sg           data.frame with columns gene_id, site_key, dm_padj,
#                     dm_delta_beta
# @param universe     character vector of universe gene IDs
# @param method, OrgDb, keyType, ont, organism, TERM2GENE, TERM2NAME,
#   padj_threshold, delta_beta_threshold, score_metric, gene_score_agg,
#   pvalueCutoff, qvalueCutoff, minGSSize, maxGSSize — same as enrichMethylation()
# @return list(go = ..., kegg = ...)
.runEnrichmentForGeneMap <- function(sg, universe, method, OrgDb, keyType, ont,
                                      organism, TERM2GENE, TERM2NAME,
                                      padj_threshold, delta_beta_threshold,
                                      score_metric, gene_score_agg,
                                      pvalueCutoff, qvalueCutoff,
                                      minGSSize, maxGSSize) {
    both <- length(method) > 1L

    go_ora   <- NULL
    go_gsea  <- NULL
    keg_ora  <- NULL
    keg_gsea <- NULL

    # ── ORA ───────────────────────────────────────────────────────────────────
    if ("ora" %in% method) {
        sig_mask <- !is.na(sg$dm_padj) & !is.na(sg$dm_delta_beta) &
                    sg$dm_padj <= padj_threshold &
                    abs(sg$dm_delta_beta) >= delta_beta_threshold

        sig_genes <- unique(sg$gene_id[sig_mask])

        if (length(sig_genes) == 0L) {
            warning(
                "No significantly differentially methylated genes found ",
                "(padj <= ", padj_threshold, " and |delta_beta| >= ",
                delta_beta_threshold, "). ORA will not be run."
            )
        } else {
            if (!is.null(TERM2GENE)) {
                go_ora <- clusterProfiler::enricher(
                    gene         = sig_genes,
                    universe     = universe,
                    TERM2GENE    = TERM2GENE,
                    TERM2NAME    = TERM2NAME,
                    pvalueCutoff = pvalueCutoff,
                    qvalueCutoff = qvalueCutoff,
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize
                )
            } else if (!is.null(OrgDb)) {
                go_ora <- clusterProfiler::enrichGO(
                    gene         = sig_genes,
                    OrgDb        = OrgDb,
                    keyType      = keyType,
                    ont          = ont,
                    universe     = universe,
                    pvalueCutoff = pvalueCutoff,
                    qvalueCutoff = qvalueCutoff,
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize,
                    readable     = FALSE
                )
            }

            if (!is.null(organism)) {
                keg_ora <- clusterProfiler::enrichKEGG(
                    gene         = sig_genes,
                    organism     = organism,
                    keyType      = keyType,
                    universe     = universe,
                    pvalueCutoff = pvalueCutoff,
                    qvalueCutoff = qvalueCutoff,
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize
                )
            }
        }
    }

    # ── GSEA ──────────────────────────────────────────────────────────────────
    if ("gsea" %in% method) {
        gene_scores <- .computeGeneScores(sg, score_metric, gene_score_agg)

        if (length(gene_scores) == 0L) {
            warning("No valid gene scores computed. GSEA will not be run.")
        } else {
            if (!is.null(TERM2GENE)) {
                go_gsea <- clusterProfiler::GSEA(
                    geneList     = gene_scores,
                    TERM2GENE    = TERM2GENE,
                    TERM2NAME    = TERM2NAME,
                    pvalueCutoff = pvalueCutoff,
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize
                )
            } else if (!is.null(OrgDb)) {
                go_gsea <- clusterProfiler::gseGO(
                    geneList     = gene_scores,
                    OrgDb        = OrgDb,
                    keyType      = keyType,
                    ont          = ont,
                    pvalueCutoff = pvalueCutoff,
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize
                )
            }

            if (!is.null(organism)) {
                keg_gsea <- clusterProfiler::gseKEGG(
                    geneList     = gene_scores,
                    organism     = organism,
                    keyType      = keyType,
                    pvalueCutoff = pvalueCutoff,
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize
                )
            }
        }
    }

    if (both) {
        list(go   = list(ora = go_ora,  gsea = go_gsea),
             kegg = list(ora = keg_ora, gsea = keg_gsea))
    } else if ("ora" %in% method) {
        list(go = go_ora,  kegg = keg_ora)
    } else {
        list(go = go_gsea, kegg = keg_gsea)
    }
}

# ─── enrichMethylation() ──────────────────────────────────────────────────────

#' Gene set enrichment analysis of differential methylation results
#'
#' Maps per-site differential methylation statistics up to the gene level and
#' runs GO and/or KEGG enrichment analysis using \pkg{clusterProfiler}.
#' Supports both over-representation analysis (ORA) and gene set enrichment
#' analysis (GSEA).  Works with model organisms via standard Bioconductor
#' databases, or with any organism via a custom \code{TERM2GENE} mapping.
#'
#' Different feature types ask different biological questions, and the
#' appropriate comparison universe changes accordingly.  Use \code{gene_role}
#' to specify which perspective to test:
#' \describe{
#'   \item{\code{"target"}}{(default) Which target genes have methylated
#'     sites in/near this feature type?  Universe = all genes with site
#'     coverage.}
#'   \item{\code{"regulator"}}{Which regulators (sigma factors, TFs, RNA
#'     regulators) have methylated binding sites?  Universe = all regulators
#'     of the same type found in the annotation for this feature type.}
#'   \item{\code{"both"}}{Run target and regulator enrichments separately,
#'     returning a named sub-list \code{list(target=..., regulator=...)}.}
#' }
#'
#' @section Prerequisites:
#' Before calling \code{enrichMethylation()}, you must:
#' \enumerate{
#'   \item Run \code{\link{diffMethyl}} to compute per-site \code{dm_padj}
#'     and \code{dm_delta_beta} values.
#'   \item Run \code{\link{annotateSites}} to assign feature identifiers to
#'     each site.  For regulatory feature types, pass
#'     \code{metadata_cols = c("feature_subtype", "transcription_unit")}
#'     to capture sigma factor identity and target gene information.
#' }
#'
#' @section Gene-to-pathway mapping:
#' Supply at least one of the following:
#' \describe{
#'   \item{\code{TERM2GENE}}{A two-column \code{data.frame} with columns
#'     \code{term} (pathway/GO term ID) and \code{gene} (gene identifier
#'     matching values in \code{gene_col}).  Used for the GO slot when
#'     provided; takes precedence over \code{OrgDb}.}
#'   \item{\code{OrgDb}}{A Bioconductor \code{OrgDb} annotation object (e.g.,
#'     \code{org.EcK12.eg.db}) for GO enrichment.  Ignored when
#'     \code{TERM2GENE} is supplied.}
#'   \item{\code{organism}}{A KEGG organism code (e.g., \code{"eco"} for
#'     \emph{E. coli} K-12).  Requires internet access; see
#'     \code{\link[clusterProfiler]{enrichKEGG}}.}
#' }
#'
#' @section Feature-type-specific parsing:
#' Gene names are extracted from the annotation using feature-type-aware rules:
#' \describe{
#'   \item{\code{"gene"}, CDS, etc.}{Feature name is used directly.}
#'   \item{\code{"promoter"}, \code{"minus_10_signal"},
#'     \code{"minus_35_signal"}, \code{"transcription_factor_binding_site"}}{
#'     Strip trailing \code{p} or \code{p\eqn{N}} suffix, then split on
#'     \code{"-"} for operonic promoters.}
#'   \item{\code{"protein_binding_site"}}{Target gene from
#'     \code{transcription_unit} attribute (if available via
#'     \code{metadata_cols}) or from \emph{"... of \{promoter\}"} pattern.
#'     Regulator gene from binding-protein name.}
#'   \item{\code{"RNA_binding_site"}}{Target from \code{transcription_unit}
#'     or \emph{"regulating \{gene\}"} pattern.  Regulator from first word if
#'     it matches a gene-name pattern.}
#'   \item{\code{"terminator"}}{Strip \emph{" terminator"} suffix.}
#' }
#'
#' @section GSEA ranking:
#' When \code{method} includes \code{"gsea"}, a per-gene score is computed
#' from the site-level statistics and genes are ranked in decreasing order.
#' Three metrics are available via \code{score_metric}:
#' \describe{
#'   \item{\code{"combined"}}{(default) \eqn{-\log_{10}(\text{padj}) \times
#'     \text{sign}(\Delta\beta)}; captures both significance and direction.}
#'   \item{\code{"padj"}}{\eqn{-\log_{10}(\text{padj})}; significance only,
#'     direction-agnostic.}
#'   \item{\code{"delta_beta"}}{\eqn{\Delta\beta} directly; effect-size
#'     ranking.}
#' }
#' When multiple sites map to the same gene, the score is aggregated across
#' sites using \code{gene_score_agg}: \code{"max"} (default; selects the site
#' with the largest absolute score, preserving its sign) or \code{"mean"}.
#'
#' @param object A \code{\link{commaData}} object on which
#'   \code{\link{diffMethyl}} and \code{\link{annotateSites}} have been run,
#'   \strong{or} a \code{data.frame} produced by \code{\link{results}()}.
#' @param method Character vector; one or both of \code{"ora"} (over-
#'   representation analysis) and \code{"gsea"} (gene set enrichment
#'   analysis).  Default \code{"ora"}.
#' @param OrgDb A Bioconductor \code{OrgDb} annotation object for GO
#'   enrichment (e.g., \code{org.EcK12.eg.db}).  \code{NULL} skips GO
#'   analysis unless \code{TERM2GENE} is provided.
#' @param keyType Character string; key type for gene IDs in \code{OrgDb}
#'   (e.g., \code{"SYMBOL"}, \code{"ENTREZID"}).  Default \code{"SYMBOL"}.
#' @param ont Character string; GO ontology.  One of \code{"BP"},
#'   \code{"MF"}, \code{"CC"}, or \code{"ALL"}.  Default \code{"BP"}.
#' @param organism Character string; KEGG organism code (e.g., \code{"eco"}).
#'   \code{NULL} skips KEGG.  Requires internet access.
#' @param TERM2GENE A two-column \code{data.frame} with columns \code{term}
#'   and \code{gene}.  Takes precedence over \code{OrgDb}.
#' @param TERM2NAME Optional two-column \code{data.frame} with columns
#'   \code{term} and \code{name}.
#' @param gene_col Character string; the \code{rowData} column containing gene
#'   identifiers per site.  Default \code{"feature_names"}.
#' @param feature_type Character vector or \code{NULL}. Feature type(s) to
#'   analyse.  When more than one type is given, each is run separately and
#'   results are returned as a named list.  Default \code{"gene"}.  Set
#'   \code{NULL} to include all annotated sites.
#' @param gene_role Character string; which role to test.  One of
#'   \code{"target"} (default), \code{"regulator"}, or \code{"both"}.
#'   See Details.
#' @param overlap_only Logical, \code{NULL}, or a named logical vector.
#'   When \code{TRUE}, only sites where \code{rel_position == 0} (inside the
#'   feature) contribute.  \code{NULL} (default) auto-detects: \code{TRUE}
#'   for large region features (gene, CDS, etc.), \code{FALSE} for small
#'   regulatory features (promoter, binding sites, etc.).
#' @param padj_threshold Numeric; ORA significance threshold. Default
#'   \code{0.05}.
#' @param delta_beta_threshold Numeric; minimum \eqn{|\Delta\beta|} for ORA.
#'   Default \code{0.1}.
#' @param score_metric Character string; GSEA ranking metric.  One of
#'   \code{"combined"} (default), \code{"padj"}, or \code{"delta_beta"}.
#' @param gene_score_agg Character string; aggregation across sites per gene.
#'   Either \code{"max"} (default) or \code{"mean"}.
#' @param mod_type Character string or \code{NULL}; modification-type filter
#'   passed to \code{\link{results}}.  Ignored when \code{object} is a
#'   \code{data.frame}.
#' @param mod_context Character string or \code{NULL}; mod-context filter
#'   passed to \code{\link{results}}.  Ignored for \code{data.frame} input.
#' @param pvalueCutoff Numeric; p-value cutoff for \pkg{clusterProfiler}.
#'   Default \code{0.05}.
#' @param qvalueCutoff Numeric; q-value cutoff (ORA only). Default
#'   \code{0.2}.
#' @param minGSSize Integer; minimum gene set size. Default \code{10}.
#' @param maxGSSize Integer; maximum gene set size. Default \code{500}.
#'
#' @return When a single \code{feature_type} and \code{gene_role != "both"}
#'   are requested: a named list \code{list(go = ..., kegg = ...)}
#'   (backward-compatible).  Each element is a \pkg{clusterProfiler}
#'   \code{enrichResult} (ORA) or \code{gseaResult} (GSEA), or \code{NULL}.
#'   When both ORA and GSEA are requested, each element is a nested list
#'   \code{list(ora = ..., gsea = ...)}.  When multiple \code{feature_type}
#'   values are given, a named list per feature type is returned.  When
#'   \code{gene_role = "both"}, each feature type's result is a list
#'   \code{list(target = ..., regulator = ...)}.
#'
#' @seealso \code{\link{diffMethyl}}, \code{\link{annotateSites}},
#'   \code{\link{results}}, \code{\link{filterResults}}
#'
#' @examples
#' \donttest{
#' # Requires clusterProfiler and a custom TERM2GENE mapping
#' if (requireNamespace("clusterProfiler", quietly = TRUE)) {
#'   data(comma_example_data)
#'   dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
#'   ann <- annotateSites(dm, annotation(comma_example_data), keep = "overlap")
#'
#'   # Custom TERM2GENE (works without network access or OrgDb)
#'   fake_t2g <- data.frame(
#'     term = c("PATH:01", "PATH:01", "PATH:02"),
#'     gene = c("geneA",  "geneB",   "geneC")
#'   )
#'   res <- enrichMethylation(ann, TERM2GENE = fake_t2g, method = c("ora", "gsea"))
#'   str(res, max.level = 2)
#' }
#' }
#'
#' @export
enrichMethylation <- function(object,
                               method               = "ora",
                               OrgDb                = NULL,
                               keyType              = "SYMBOL",
                               ont                  = "BP",
                               organism             = NULL,
                               TERM2GENE            = NULL,
                               TERM2NAME            = NULL,
                               gene_col             = "feature_names",
                               feature_type         = "gene",
                               gene_role            = c("target", "regulator", "both"),
                               overlap_only         = NULL,
                               padj_threshold       = 0.05,
                               delta_beta_threshold = 0.1,
                               score_metric         = "combined",
                               gene_score_agg       = "max",
                               mod_type             = NULL,
                               mod_context          = NULL,
                               pvalueCutoff         = 0.05,
                               qvalueCutoff         = 0.2,
                               minGSSize            = 10L,
                               maxGSSize            = 500L) {

    # ── Validate inputs ───────────────────────────────────────────────────────
    method         <- match.arg(method, choices = c("ora", "gsea"), several.ok = TRUE)
    gene_role      <- match.arg(gene_role)
    ont            <- match.arg(ont,            c("BP", "MF", "CC", "ALL"))
    score_metric   <- match.arg(score_metric,   c("combined", "padj", "delta_beta"))
    gene_score_agg <- match.arg(gene_score_agg, c("max", "mean"))

    # ── Validate object type first (before dependency checks) ─────────────────
    if (!is(object, "commaData") && !is.data.frame(object)) {
        stop("'object' must be a commaData object or a data.frame from results().")
    }

    if (is.null(OrgDb) && is.null(organism) && is.null(TERM2GENE)) {
        stop(
            "No gene-to-term mapping supplied. Provide at least one of:\n",
            "  OrgDb     -- Bioconductor OrgDb object for GO (e.g., org.EcK12.eg.db)\n",
            "  organism  -- KEGG organism code (e.g., \"eco\")\n",
            "  TERM2GENE -- custom data.frame(term, gene)"
        )
    }

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop(
            "Package 'clusterProfiler' is required.\n",
            "Install it with: BiocManager::install(\"clusterProfiler\")"
        )
    }

    # ── Extract results data.frame ────────────────────────────────────────────
    if (is(object, "commaData")) {
        res_df <- results(object, mod_type = mod_type, mod_context = mod_context)
    } else {
        res_df <- object
        if (!is.null(mod_type) || !is.null(mod_context)) {
            warning(
                "'object' is a data.frame; 'mod_type' and 'mod_context' filters ",
                "are ignored. Pre-filter the data.frame before calling enrichMethylation()."
            )
        }
    }

    if (!"dm_padj" %in% colnames(res_df)) {
        stop(
            "No differential methylation results found.\n",
            "Run diffMethyl() first:\n",
            "  dm <- diffMethyl(object, formula = ~ condition)"
        )
    }

    # ── Build the per-feature-type loop ───────────────────────────────────────
    ft_loop   <- if (is.null(feature_type)) list(NULL) else as.list(feature_type)
    is_single <- length(ft_loop) == 1L

    results_by_ft <- lapply(ft_loop, function(ft) {
        # Determine effective overlap_only for this feature type
        if (is.null(ft)) {
            eff_overlap_only <- FALSE
        } else if (is.null(overlap_only)) {
            eff_overlap_only <- ft %in% .COMMA_REGION_FEATURE_TYPES
        } else if (is.logical(overlap_only) && !is.null(names(overlap_only))) {
            val <- overlap_only[ft]
            eff_overlap_only <- if (is.na(val)) ft %in% .COMMA_REGION_FEATURE_TYPES else val
        } else {
            eff_overlap_only <- as.logical(overlap_only)[1L]
        }

        # Build site-gene role map
        if (is.null(ft)) {
            # NULL feature_type: use legacy siteToGeneMap (no role information)
            sg_all <- .siteToGeneMap(res_df, gene_col)
            if (is.null(sg_all) || nrow(sg_all) == 0L) {
                warning("Gene map is empty — returning NULL.")
                return(list(go = NULL, kegg = NULL))
            }
            universe <- unique(sg_all$gene_id[!is.na(sg_all$gene_id)])
            return(.runEnrichmentForGeneMap(sg_all, universe, method,
                                            OrgDb, keyType, ont, organism,
                                            TERM2GENE, TERM2NAME,
                                            padj_threshold, delta_beta_threshold,
                                            score_metric, gene_score_agg,
                                            pvalueCutoff, qvalueCutoff,
                                            minGSSize, maxGSSize))
        }

        # New code path: use .extractGeneRoles()
        if (!"feature_types" %in% colnames(res_df)) {
            warning(
                "'feature_types' column not found. Run annotateSites() first. ",
                "Returning NULL for feature_type = '", ft, "'."
            )
            return(list(go = NULL, kegg = NULL))
        }

        role_table <- .extractGeneRoles(res_df, ft, gene_col, eff_overlap_only)

        if (is.null(role_table)) {
            warning("No sites with feature_type '", ft, "'. Returning NULL.")
            return(list(go = NULL, kegg = NULL))
        }

        # Helper to run enrichment for one role subset
        .enrich_for_role <- function(rt, role_name) {
            sg <- rt[rt$role == role_name & !is.na(rt$gene_id), , drop = FALSE]
            if (nrow(sg) == 0L) {
                warning("No ", role_name, " genes found for feature_type '",
                        ft, "'. Returning NULL.")
                return(list(go = NULL, kegg = NULL))
            }
            if (role_name == "target") {
                # Universe = all target genes (standard background)
                universe <- unique(sg$gene_id)
            } else {
                # Universe = all regulator genes of this type in the annotation
                universe <- unique(rt$gene_id[rt$role == "regulator" &
                                                  !is.na(rt$gene_id)])
            }
            .runEnrichmentForGeneMap(sg, universe, method,
                                     OrgDb, keyType, ont, organism,
                                     TERM2GENE, TERM2NAME,
                                     padj_threshold, delta_beta_threshold,
                                     score_metric, gene_score_agg,
                                     pvalueCutoff, qvalueCutoff,
                                     minGSSize, maxGSSize)
        }

        if (gene_role == "both") {
            list(
                target    = .enrich_for_role(role_table, "target"),
                regulator = .enrich_for_role(role_table, "regulator")
            )
        } else {
            .enrich_for_role(role_table, gene_role)
        }
    })

    # Name the results by feature_type
    if (!is.null(feature_type)) {
        names(results_by_ft) <- feature_type
    }

    # Backward compatibility: single feature_type → unwrap
    if (is_single) results_by_ft[[1L]] else results_by_ft
}
