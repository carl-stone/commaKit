#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
#' @importFrom stats setNames
NULL

# --- Constants ----------------------------------------------------------------

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

# --- Existing internal helpers (kept for backward compatibility) ---------------

# Map differential methylation results from sites to genes.
#
# Explodes the CharacterList/list annotation column so each site x gene
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

# --- New internal helpers for regulatory role extraction ----------------------

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

    # -- Step 1: Find per-site indices matching ft ------------------------------
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

    # -- Step 2: Optional overlap_only filter ----------------------------------
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

    # -- Step 3: Extract raw feature names for this ft -------------------------
    gn_lists <- if (is(res_sub[[gene_col]], "CharacterList")) {
        as.list(res_sub[[gene_col]])
    } else {
        res_sub[[gene_col]]
    }
    raw_per_site <- mapply(function(nms, idx) as.character(nms[idx]),
                           gn_lists, type_idx_sub, SIMPLIFY = FALSE)

    # -- Step 4: Extract optional metadata columns ------------------------------
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

    # -- Step 5: Parse target and regulator genes -------------------------------
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
# once per feature_type x gene_role combination.
#
# @param sg           data.frame with columns gene_id, site_key, dm_padj,
#                     dm_delta_beta
# @param universe     character vector of universe gene IDs
# @param method, OrgDb, keyType, ont, organism, TERM2GENE, TERM2NAME,
#   kegg_term2gene, kegg_term2name,
#   padj_threshold, delta_beta_threshold, score_metric, gene_score_agg,
#   pvalueCutoff, qvalueCutoff, minGSSize, maxGSSize -- same as enrichMethylation()
# @return list(go = ..., kegg = ...)
.runEnrichmentForGeneMap <- function(sg, universe, method, OrgDb, keyType, ont,
                                      organism, TERM2GENE, TERM2NAME,
                                      kegg_term2gene, kegg_term2name,
                                      padj_threshold, delta_beta_threshold,
                                      score_metric, gene_score_agg,
                                      pvalueCutoff, qvalueCutoff,
                                      minGSSize, maxGSSize) {
    both <- length(method) > 1L

    go_ora   <- NULL
    go_gsea  <- NULL
    keg_ora  <- NULL
    keg_gsea <- NULL

    # -- ORA -------------------------------------------------------------------
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

            if (!is.null(kegg_term2gene)) {
                keg_t2n <- if (!is.null(kegg_term2name) && nrow(kegg_term2name) > 0L)
                               kegg_term2name else NULL
                keg_ora <- clusterProfiler::enricher(
                    gene         = sig_genes,
                    universe     = universe,
                    TERM2GENE    = kegg_term2gene,
                    TERM2NAME    = keg_t2n,
                    pvalueCutoff = pvalueCutoff,
                    qvalueCutoff = qvalueCutoff,
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize
                )
            } else if (!is.null(organism)) {
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

    # -- GSEA ------------------------------------------------------------------
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

            if (!is.null(kegg_term2gene)) {
                keg_t2n <- if (!is.null(kegg_term2name) && nrow(kegg_term2name) > 0L)
                               kegg_term2name else NULL
                keg_gsea <- clusterProfiler::GSEA(
                    geneList     = gene_scores,
                    TERM2GENE    = kegg_term2gene,
                    TERM2NAME    = keg_t2n,
                    pvalueCutoff = pvalueCutoff,
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize
                )
            } else if (!is.null(organism)) {
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

# --- buildKEGGTermGene() -------------------------------------------------------

#' Build a KEGG term-to-gene mapping for use with enrichMethylation()
#'
#' Fetches all KEGG pathway-gene associations for an organism using only two
#' API calls (\code{keggLink} and \code{keggList}), then returns the result as
#' a pair of data frames that can be passed directly to
#' \code{\link{enrichMethylation}()} via the \code{kegg_term2gene} and
#' \code{kegg_term2name} arguments.  Optionally caches the result to an RDS
#' file so that subsequent calls are fully offline.
#'
#' This function exists because \code{\link[clusterProfiler]{enrichKEGG}} and
#' \code{\link[clusterProfiler]{gseKEGG}} fire one HTTP request per pathway
#' when fetching gene lists, which quickly exceeds the KEGG API rate limit for
#' organisms with many pathways.  \code{buildKEGGTermGene} retrieves the same
#' data in two bulk calls — one for gene-pathway links and one for pathway
#' names.
#'
#' @param organism Character string; KEGG organism code (e.g., \code{"eco"} for
#'   \emph{Escherichia coli} K-12, \code{"hsa"} for \emph{Homo sapiens}).
#'   Browse organism codes at
#'   \url{https://www.genome.jp/kegg/catalog/org_list.html}.
#' @param file Character string or \code{NULL}.  Path to an RDS cache file.
#'   \itemize{
#'     \item If \code{file} already exists, the cached object is loaded and
#'       returned immediately — no network access occurs.
#'     \item If \code{file} does not exist, KEGG is queried and the result is
#'       saved to \code{file}.
#'     \item \code{NULL} (default) disables caching.
#'   }
#'   A warning is issued when a cache file is older than 90 days, suggesting
#'   the user refresh it to pick up pathway database updates.
#' @param strip_prefix Logical; if \code{TRUE} (default), the organism prefix
#'   is stripped from gene IDs (e.g., \code{"eco:b0001"} becomes
#'   \code{"b0001"}) and the \code{"path:"} prefix is stripped from pathway IDs
#'   (e.g., \code{"path:eco00010"} becomes \code{"eco00010"}).  The trailing
#'   organism qualifier is also stripped from pathway names (e.g., \code{"...
#'   - Escherichia coli K-12"} becomes \code{"..."}).
#' @param id_map A \code{data.frame} with columns \code{symbol} and
#'   \code{kegg_id} as returned by \code{\link{buildKEGGGeneIDMap}()}.
#'   When provided, the \code{gene} column of \code{term2gene} is translated
#'   from KEGG-internal IDs (e.g. b-numbers) to gene symbols before returning,
#'   so that genes match the identifiers in your annotation data.  KEGG IDs
#'   with no matching symbol are preserved as-is; no pathway genes are dropped.
#'   \code{NULL} (default) leaves identifiers unchanged.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{term2gene}}{A \code{data.frame} with columns \code{term}
#'       (KEGG pathway ID) and \code{gene} (gene ID matching the identifiers in
#'       your annotation data).}
#'     \item{\code{term2name}}{A \code{data.frame} with columns \code{term}
#'       (KEGG pathway ID) and \code{name} (human-readable pathway
#'       description).}
#'   }
#'   Pass these to \code{\link{enrichMethylation}()} as
#'   \code{kegg_term2gene = result$term2gene, kegg_term2name = result$term2name}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("KEGGREST", quietly = TRUE)) {
#'   # Fetch once and cache to disk
#'   kegg <- buildKEGGTermGene("eco", file = "eco_kegg.rds")
#'
#'   # Load from cache on subsequent calls (no network)
#'   kegg <- buildKEGGTermGene("eco", file = "eco_kegg.rds")
#'
#'   # Translate b-numbers to symbols with an id_map:
#'   # id_map <- buildKEGGGeneIDMap("eco", OrgDb = org.EcK12.eg.db)
#'   # kegg   <- buildKEGGTermGene("eco", file = "eco_kegg.rds", id_map = id_map)
#'   # res    <- enrichMethylation(obj,
#'   #             kegg_term2gene = kegg$term2gene,
#'   #             kegg_term2name = kegg$term2name)
#' }
#' }
#'
#' @seealso \code{\link{buildKEGGGeneIDMap}}, \code{\link{enrichMethylation}}
#' @export
buildKEGGTermGene <- function(organism, file = NULL, strip_prefix = TRUE,
                               id_map = NULL) {
    if (!is.character(organism) || length(organism) != 1L || nchar(organism) == 0L) {
        stop("'organism' must be a non-empty character string (e.g., \"eco\").")
    }

    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop(
            "Package 'KEGGREST' is required for buildKEGGTermGene().\n",
            "Install it with: BiocManager::install(\"KEGGREST\")"
        )
    }

    # -- Load from cache if it exists ------------------------------------------
    if (!is.null(file) && file.exists(file)) {
        age_days <- as.numeric(difftime(Sys.time(), file.mtime(file), units = "days"))
        if (age_days > 90) {
            warning(
                sprintf("KEGG cache file '%s' is %.0f days old. ", file, age_days),
                "Consider refreshing it by deleting the file and re-running ",
                "buildKEGGTermGene()."
            )
        }
        message(sprintf("Loading KEGG data from cache: %s", file))
        return(readRDS(file))
    }

    message(sprintf("Fetching KEGG pathway data for organism '%s' ...", organism))

    # -- API call 1: gene-to-pathway links (one bulk request) ------------------
    links <- tryCatch(
        KEGGREST::keggLink("pathway", organism),
        error = function(e) {
            stop(
                sprintf("KEGG API call failed for organism '%s'.\n", organism),
                "Check the organism code at:\n",
                "  https://www.genome.jp/kegg/catalog/org_list.html\n",
                "Original error: ", conditionMessage(e)
            )
        }
    )

    if (length(links) == 0L) {
        stop(sprintf(
            "No KEGG pathway associations returned for organism '%s'. ",
            organism
        ), "Verify the organism code.")
    }

    # -- API call 2: pathway descriptions (one bulk request) -------------------
    path_list <- tryCatch(
        KEGGREST::keggList("pathway", organism),
        error = function(e) {
            warning(
                "Could not fetch KEGG pathway names: ", conditionMessage(e),
                "\nterm2name will be empty."
            )
            character(0)
        }
    )

    # -- Reshape results -------------------------------------------------------
    gene_ids    <- names(links)   # e.g. "eco:b0001"
    pathway_ids <- unname(links)  # e.g. "path:eco00010"

    if (strip_prefix) {
        gene_ids    <- sub(paste0("^", organism, ":"), "", gene_ids)
        pathway_ids <- sub("^path:", "", pathway_ids)
    }

    term2gene <- data.frame(
        term = pathway_ids,
        gene = gene_ids,
        stringsAsFactors = FALSE
    )
    term2gene <- term2gene[!duplicated(term2gene), , drop = FALSE]
    rownames(term2gene) <- NULL

    # -- Optionally translate KEGG IDs to gene symbols -------------------------
    if (!is.null(id_map)) {
        .validateKEGGIDMap(id_map)
        mapped <- id_map$symbol[match(term2gene$gene, id_map$kegg_id)]
        term2gene$gene <- ifelse(is.na(mapped), term2gene$gene, mapped)
    }

    if (length(path_list) > 0L) {
        pw_ids   <- names(path_list)   # e.g. "path:eco00010"
        pw_names <- unname(path_list)  # e.g. "Glycolysis / Gluconeogenesis - Escherichia coli K-12"
        if (strip_prefix) {
            pw_ids   <- sub("^path:", "", pw_ids)
            pw_names <- vapply(strsplit(pw_names, " - "), function(parts) {
                if (length(parts) <= 1L) parts
                else paste(parts[-length(parts)], collapse = " - ")
            }, character(1L))
        }
        term2name <- data.frame(
            term = pw_ids,
            name = pw_names,
            stringsAsFactors = FALSE
        )
    } else {
        term2name <- data.frame(
            term = character(0),
            name = character(0),
            stringsAsFactors = FALSE
        )
    }

    result <- list(term2gene = term2gene, term2name = term2name)

    # -- Save to cache ---------------------------------------------------------
    if (!is.null(file)) {
        saveRDS(result, file)
        message(sprintf("KEGG data cached to: %s", file))
    }

    message(sprintf(
        "Done. %d gene-pathway associations across %d pathways.",
        nrow(term2gene), nrow(term2name)
    ))

    result
}

# --- buildKEGGGeneIDMap() -----------------------------------------------------

# Validate that id_map is a data.frame with columns 'symbol' and 'kegg_id'.
.validateKEGGIDMap <- function(id_map) {
    if (!is.data.frame(id_map) ||
            !all(c("symbol", "kegg_id") %in% colnames(id_map))) {
        stop(
            "'id_map' must be a data.frame with columns 'symbol' and 'kegg_id'.\n",
            "Use buildKEGGGeneIDMap() to create it."
        )
    }
}

#' Build a KEGG gene ID map for symbol translation
#'
#' Fetches the complete NCBI Gene ID \eqn{\leftrightarrow} KEGG gene ID
#' correspondence for an organism in a single API call
#' (\code{KEGGREST::keggConv}), then joins it against a gene symbol table
#' supplied either as a Bioconductor \code{OrgDb} object or as a plain
#' \code{data.frame}.  The result maps each gene symbol to the organism's
#' KEGG-internal gene identifier (e.g. b-numbers for \emph{E. coli} K-12),
#' and can be passed to \code{\link{buildKEGGTermGene}()} to translate the
#' \code{gene} column of \code{term2gene} from KEGG IDs to symbols that match
#' your annotation data.
#'
#' @section Why this is needed:
#' \code{\link{buildKEGGTermGene}()} returns pathway-gene associations where
#' genes are identified by KEGG-internal IDs (b-numbers for \emph{E. coli},
#' Entrez IDs for human, etc.).  Annotation data produced by
#' \code{\link{loadAnnotation}()} uses gene symbols from the GFF3
#' \code{Name=} attribute.  Without translation, no genes will match during
#' enrichment.
#'
#' @section Input modes:
#' Exactly one of \code{OrgDb} or \code{entrez2symbol} must be supplied:
#' \describe{
#'   \item{\code{OrgDb}}{A Bioconductor \code{OrgDb} annotation object (e.g.,
#'     \code{org.EcK12.eg.db}).  NCBI Gene IDs (\code{id_col}) and gene symbols
#'     (\code{keys_col}) are extracted via \code{AnnotationDbi::select()}.
#'     Requires internet access only for the single \code{keggConv} call.}
#'   \item{\code{entrez2symbol}}{A two-column \code{data.frame} with columns
#'     \code{entrez_id} (character NCBI Gene IDs) and \code{symbol} (gene
#'     symbols).  Use this when no \code{OrgDb} package is available, e.g.
#'     when you have downloaded a mapping from NCBI's \code{gene_info} file.}
#' }
#'
#' @param organism Character string; KEGG organism code (e.g., \code{"eco"}
#'   for \emph{Escherichia coli} K-12).  See
#'   \url{https://www.genome.jp/kegg/catalog/org_list.html}.
#' @param OrgDb A Bioconductor \code{OrgDb} object for automatic symbol
#'   lookup.  \code{NULL} (default) if \code{entrez2symbol} is provided
#'   instead.
#' @param entrez2symbol A \code{data.frame} with columns \code{entrez_id} and
#'   \code{symbol}.  \code{NULL} (default) if \code{OrgDb} is provided.
#' @param keys_col Character string; the \code{OrgDb} column containing gene
#'   symbols.  Default \code{"SYMBOL"}.  Ignored when \code{entrez2symbol} is
#'   provided.
#' @param id_col Character string; the \code{OrgDb} column containing NCBI
#'   Gene IDs.  Default \code{"ENTREZID"}.  Ignored when \code{entrez2symbol}
#'   is provided.
#' @param file Character string or \code{NULL}.  Path to an RDS cache file.
#'   \itemize{
#'     \item If \code{file} already exists, it is loaded and returned without
#'       any API calls.
#'     \item If \code{file} does not exist, the mapping is built and saved.
#'     \item \code{NULL} (default) disables caching.
#'   }
#'   A warning is issued when the cache is older than 90 days.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{\code{symbol}}{Gene symbol matching identifiers in your annotation
#'       data (e.g. \code{"lacZ"}, \code{"rpoD"}).}
#'     \item{\code{kegg_id}}{Corresponding KEGG gene identifier after stripping
#'       the organism prefix (e.g. \code{"b0344"} for \emph{E. coli} K-12).}
#'   }
#'   Pass this to \code{\link{buildKEGGTermGene}()} as \code{id_map = ...} to
#'   translate the \code{gene} column of \code{term2gene} from KEGG IDs to
#'   symbols.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("KEGGREST", quietly = TRUE) &&
#'     requireNamespace("org.EcK12.eg.db", quietly = TRUE)) {
#'   id_map <- buildKEGGGeneIDMap("eco",
#'                                OrgDb = org.EcK12.eg.db,
#'                                file  = "eco_id_map.rds")
#'
#'   # Apply when building pathway map:
#'   kegg <- buildKEGGTermGene("eco", file = "eco_kegg.rds", id_map = id_map)
#'   # kegg$term2gene$gene now contains symbols instead of b-numbers
#' }
#'
#' # Manual table alternative:
#' if (requireNamespace("KEGGREST", quietly = TRUE)) {
#'   ent2sym <- data.frame(
#'     entrez_id = c("945076", "945803"),
#'     symbol    = c("lacZ",   "lacY"),
#'     stringsAsFactors = FALSE
#'   )
#'   id_map <- buildKEGGGeneIDMap("eco", entrez2symbol = ent2sym)
#' }
#' }
#'
#' @seealso \code{\link{buildKEGGTermGene}}, \code{\link{enrichMethylation}}
#' @export
buildKEGGGeneIDMap <- function(organism,
                                OrgDb         = NULL,
                                entrez2symbol = NULL,
                                keys_col      = "SYMBOL",
                                id_col        = "ENTREZID",
                                file          = NULL) {
    # -- Input validation ------------------------------------------------------
    if (!is.character(organism) || length(organism) != 1L || nchar(organism) == 0L) {
        stop("'organism' must be a non-empty character string (e.g., \"eco\").")
    }
    if (is.null(OrgDb) && is.null(entrez2symbol)) {
        stop(
            "Provide at least one of:\n",
            "  OrgDb         -- a Bioconductor OrgDb object\n",
            "  entrez2symbol -- a data.frame(entrez_id, symbol)"
        )
    }
    if (!is.null(entrez2symbol)) {
        if (!is.data.frame(entrez2symbol) ||
                !all(c("entrez_id", "symbol") %in% colnames(entrez2symbol))) {
            stop(
                "'entrez2symbol' must be a data.frame with columns ",
                "'entrez_id' and 'symbol'."
            )
        }
    }

    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop(
            "Package 'KEGGREST' is required.\n",
            "Install it with: BiocManager::install(\"KEGGREST\")"
        )
    }

    # -- Load from cache -------------------------------------------------------
    if (!is.null(file) && file.exists(file)) {
        age_days <- as.numeric(difftime(Sys.time(), file.mtime(file), units = "days"))
        if (age_days > 90) {
            warning(
                sprintf("KEGG ID map cache '%s' is %.0f days old. ", file, age_days),
                "Consider refreshing it by deleting the file and re-running ",
                "buildKEGGGeneIDMap()."
            )
        }
        message(sprintf("Loading KEGG ID map from cache: %s", file))
        return(readRDS(file))
    }

    message(sprintf(
        "Fetching KEGG gene ID map for organism '%s' ...", organism
    ))

    # -- API call: NCBI Gene ID <-> KEGG gene ID (one bulk request) ------------
    conv <- tryCatch(
        KEGGREST::keggConv(organism, "ncbi-geneid"),
        error = function(e) {
            stop(
                sprintf("KEGG API call failed for organism '%s'.\n", organism),
                "Check the organism code at:\n",
                "  https://www.genome.jp/kegg/catalog/org_list.html\n",
                "Original error: ", conditionMessage(e)
            )
        }
    )

    if (length(conv) == 0L) {
        stop(sprintf(
            "keggConv returned no entries for organism '%s'. ",
            organism
        ), "Verify the organism code.")
    }

    # conv: names = "ncbi-geneid:945076", values = "eco:b0344"
    entrez_ids <- sub("^ncbi-geneid:", "", names(conv))
    kegg_ids   <- sub(paste0("^", organism, ":"), "", unname(conv))

    kegg_conv_df <- data.frame(
        entrez_id = entrez_ids,
        kegg_id   = kegg_ids,
        stringsAsFactors = FALSE
    )

    # -- Get ENTREZID <-> SYMBOL map -------------------------------------------
    if (!is.null(entrez2symbol)) {
        sym_df <- entrez2symbol[, c("entrez_id", "symbol"), drop = FALSE]
        sym_df$entrez_id <- as.character(sym_df$entrez_id)
    } else {
        # OrgDb path
        if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
            stop(
                "Package 'AnnotationDbi' is required when using OrgDb.\n",
                "Install it with: BiocManager::install(\"AnnotationDbi\")"
            )
        }
        avail_cols <- AnnotationDbi::columns(OrgDb)
        for (col in c(id_col, keys_col)) {
            if (!col %in% avail_cols) {
                stop(sprintf(
                    "Column '%s' not found in OrgDb. Available columns:\n  %s",
                    col, paste(avail_cols, collapse = ", ")
                ))
            }
        }
        sym_raw <- AnnotationDbi::select(
            OrgDb,
            keys    = AnnotationDbi::keys(OrgDb, keytype = id_col),
            columns = c(id_col, keys_col),
            keytype = id_col
        )
        sym_df <- data.frame(
            entrez_id = as.character(sym_raw[[id_col]]),
            symbol    = as.character(sym_raw[[keys_col]]),
            stringsAsFactors = FALSE
        )
    }

    # Drop NA rows from symbol table
    sym_df <- sym_df[!is.na(sym_df$entrez_id) & !is.na(sym_df$symbol), ,
                     drop = FALSE]

    # -- Join kegg_conv_df with sym_df on entrez_id ----------------------------
    matched_symbol <- sym_df$symbol[match(kegg_conv_df$entrez_id,
                                          sym_df$entrez_id)]
    result_df <- data.frame(
        symbol  = matched_symbol,
        kegg_id = kegg_conv_df$kegg_id,
        stringsAsFactors = FALSE
    )
    result_df <- result_df[!is.na(result_df$symbol), , drop = FALSE]
    result_df <- result_df[!duplicated(result_df), , drop = FALSE]
    rownames(result_df) <- NULL

    if (nrow(result_df) == 0L) {
        warning(
            "No symbol matches found after joining KEGG and ",
            if (is.null(entrez2symbol)) "OrgDb" else "entrez2symbol",
            " data. Check that gene identifiers overlap."
        )
    }

    # -- Save to cache ---------------------------------------------------------
    if (!is.null(file)) {
        saveRDS(result_df, file)
        message(sprintf("KEGG ID map cached to: %s", file))
    }

    message(sprintf(
        "Done. %d gene symbols mapped to KEGG IDs.",
        nrow(result_df)
    ))

    result_df
}

# --- enrichMethylation() ------------------------------------------------------

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
#'   \item{\code{kegg_term2gene}}{(Recommended for KEGG) A two-column
#'     \code{data.frame} pre-built by \code{\link{buildKEGGTermGene}()}.
#'     No network access is required at analysis time; results appear in the
#'     \code{$kegg} slot.}
#'   \item{\code{TERM2GENE}}{A two-column \code{data.frame} for custom GO
#'     enrichment.  Results appear in the \code{$go} slot.}
#'   \item{\code{OrgDb}}{A Bioconductor \code{OrgDb} annotation object (e.g.,
#'     \code{org.EcK12.eg.db}) for GO enrichment.  Ignored when
#'     \code{TERM2GENE} is supplied.}
#'   \item{\code{organism}}{A KEGG organism code (e.g., \code{"eco"}).
#'     Makes one live HTTP request per KEGG pathway; may exceed the API rate
#'     limit.  Ignored when \code{kegg_term2gene} is provided.}
#' }
#'
#' @section Feature-type-specific parsing:
#' Gene names are extracted from the annotation using feature-type-aware rules:
#' \describe{
#'   \item{\code{"gene"}, CDS, etc.}{Feature name is used directly.}
#'   \item{\code{"promoter"}, \code{"minus_10_signal"},
#'     \code{"minus_35_signal"}, \code{"transcription_factor_binding_site"}}{
#'     Strip trailing \code{p} or \code{p}\eqn{N} suffix, then split on
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
#'   \code{NULL} skips KEGG.  Requires internet access and fires one HTTP
#'   request per KEGG pathway, which can exceed the API rate limit for
#'   organisms with many pathways.  Prefer \code{kegg_term2gene} (built once
#'   via \code{\link{buildKEGGTermGene}()}) for reliable KEGG analysis.
#' @param TERM2GENE A two-column \code{data.frame} with columns \code{term}
#'   and \code{gene}.  When provided, takes precedence over \code{OrgDb} for
#'   GO analysis; results appear in the \code{$go} slot.  To supply a
#'   pre-built KEGG mapping, use \code{kegg_term2gene} instead so that results
#'   land in the \code{$kegg} slot.
#' @param TERM2NAME Optional two-column \code{data.frame} with columns
#'   \code{term} and \code{name}.  Paired with \code{TERM2GENE} for GO
#'   enrichment.
#' @param kegg_term2gene A two-column \code{data.frame} with columns
#'   \code{term} and \code{gene} containing pre-fetched KEGG pathway-gene
#'   associations.  Build it once per organism with
#'   \code{\link{buildKEGGTermGene}()} and reuse across analyses — no live
#'   KEGG API calls are made.  Results appear in the \code{$kegg} slot,
#'   consistent with the \code{organism} path.  Takes precedence over
#'   \code{organism} when both are supplied.
#' @param kegg_term2name Optional two-column \code{data.frame} with columns
#'   \code{term} and \code{name} containing KEGG pathway descriptions.
#'   Returned by \code{\link{buildKEGGTermGene}()} alongside
#'   \code{kegg_term2gene}.
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
                               kegg_term2gene       = NULL,
                               kegg_term2name       = NULL,
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

    # -- Validate inputs -------------------------------------------------------
    method         <- match.arg(method, choices = c("ora", "gsea"), several.ok = TRUE)
    gene_role      <- match.arg(gene_role)
    ont            <- match.arg(ont,            c("BP", "MF", "CC", "ALL"))
    score_metric   <- match.arg(score_metric,   c("combined", "padj", "delta_beta"))
    gene_score_agg <- match.arg(gene_score_agg, c("max", "mean"))

    # -- Validate object type first (before dependency checks) -----------------
    if (!is(object, "commaData") && !is.data.frame(object)) {
        stop("'object' must be a commaData object or a data.frame from results().")
    }

    if (is.null(OrgDb) && is.null(organism) && is.null(TERM2GENE) &&
            is.null(kegg_term2gene)) {
        stop(
            "No gene-to-term mapping supplied. Provide at least one of:\n",
            "  kegg_term2gene -- pre-built KEGG mapping from buildKEGGTermGene()\n",
            "  OrgDb          -- Bioconductor OrgDb object for GO (e.g., org.EcK12.eg.db)\n",
            "  organism       -- KEGG organism code (e.g., \"eco\") [live API]\n",
            "  TERM2GENE      -- custom data.frame(term, gene)"
        )
    }

    if (!is.null(kegg_term2gene)) {
        if (!is.data.frame(kegg_term2gene) ||
                !all(c("term", "gene") %in% colnames(kegg_term2gene))) {
            stop(
                "'kegg_term2gene' must be a data.frame with columns 'term' and 'gene'.\n",
                "Use buildKEGGTermGene() to create it."
            )
        }
    }
    if (!is.null(kegg_term2name)) {
        if (!is.data.frame(kegg_term2name) ||
                !all(c("term", "name") %in% colnames(kegg_term2name))) {
            stop(
                "'kegg_term2name' must be a data.frame with columns 'term' and 'name'."
            )
        }
    }

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop(
            "Package 'clusterProfiler' is required.\n",
            "Install it with: BiocManager::install(\"clusterProfiler\")"
        )
    }

    # -- Extract results data.frame --------------------------------------------
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

    # -- Build the per-feature-type loop ---------------------------------------
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
                warning("Gene map is empty -- returning NULL.")
                return(list(go = NULL, kegg = NULL))
            }
            universe <- unique(sg_all$gene_id[!is.na(sg_all$gene_id)])
            return(.runEnrichmentForGeneMap(sg_all, universe, method,
                                            OrgDb, keyType, ont, organism,
                                            TERM2GENE, TERM2NAME,
                                            kegg_term2gene, kegg_term2name,
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
                                     kegg_term2gene, kegg_term2name,
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

    # Backward compatibility: single feature_type -> unwrap
    if (is_single) results_by_ft[[1L]] else results_by_ft
}
