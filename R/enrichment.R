#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
NULL

# ─── Internal helpers ─────────────────────────────────────────────────────────

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
            "  object <- annotateSites(object, type = \"overlap\")"
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

    data.frame(
        gene_id       = unlist(gene_lists_sub, use.names = FALSE),
        site_key      = rep(site_keys, lens_sub),
        dm_padj       = rep(res_sub$dm_padj, lens_sub),
        dm_delta_beta = rep(res_sub$dm_delta_beta, lens_sub),
        stringsAsFactors = FALSE
    )
}


# Compute a per-gene methylation score for GSEA ranking.
#
# @param site_gene_df  data.frame from .siteToGeneMap()
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
    gene_scores <- vapply(genes, function(g) {
        s <- df$site_score[df$gene_id == g]
        if (agg == "max") s[which.max(abs(s))] else mean(s)
    }, FUN.VALUE = numeric(1L))
    names(gene_scores) <- genes

    sort(gene_scores, decreasing = TRUE)
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
#' @section Prerequisites:
#' Before calling \code{enrichMethylation()}, you must:
#' \enumerate{
#'   \item Run \code{\link{diffMethyl}} to compute per-site \code{dm_padj}
#'     and \code{dm_delta_beta} values.
#'   \item Run \code{\link{annotateSites}} (with \code{type = "overlap"}) to
#'     assign gene identifiers to each site.  The resulting
#'     \code{feature_names} column (a \code{CharacterList}) is used by default.
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
#'   \code{\link{diffMethyl}} \strong{and} \code{\link{annotateSites}} have
#'   been run.
#' @param method Character vector; one or both of \code{"ora"} (over-
#'   representation analysis) and \code{"gsea"} (gene set enrichment
#'   analysis).  Default \code{"ora"}.
#' @param OrgDb A Bioconductor \code{OrgDb} annotation object for GO
#'   enrichment (e.g., \code{org.EcK12.eg.db}).  If \code{NULL} and no
#'   \code{TERM2GENE} is provided, GO analysis is skipped.
#' @param keyType Character string specifying the key type used for gene IDs
#'   in \code{OrgDb} (e.g., \code{"SYMBOL"}, \code{"ENTREZID"}).
#'   Must match the identifiers stored in \code{gene_col}.  Default
#'   \code{"SYMBOL"}.
#' @param ont Character string specifying the GO ontology to use when
#'   \code{OrgDb} is provided.  One of \code{"BP"}, \code{"MF"}, \code{"CC"},
#'   or \code{"ALL"}.  Default \code{"BP"}.
#' @param organism Character string; KEGG organism code (e.g., \code{"eco"}).
#'   If \code{NULL}, KEGG analysis is skipped.  Requires internet access.
#' @param TERM2GENE A two-column \code{data.frame} with columns \code{term}
#'   and \code{gene} mapping pathway/term identifiers to gene identifiers.
#'   When provided, this is used for GO analysis via
#'   \code{\link[clusterProfiler]{enricher}} (ORA) or
#'   \code{\link[clusterProfiler]{GSEA}} and takes precedence over
#'   \code{OrgDb}.  Ignored for KEGG (use \code{organism} for KEGG).
#' @param TERM2NAME An optional two-column \code{data.frame} with columns
#'   \code{term} and \code{name} providing human-readable descriptions of
#'   term IDs in \code{TERM2GENE}.
#' @param gene_col Character string; the \code{rowData} column containing gene
#'   identifiers per site (a \code{CharacterList} or \code{list} column added
#'   by \code{\link{annotateSites}}).  Default \code{"feature_names"}.
#' @param padj_threshold Numeric; adjusted p-value threshold for classifying
#'   a site as differentially methylated in ORA.  Default \code{0.05}.
#' @param delta_beta_threshold Numeric; minimum absolute effect size
#'   (\eqn{|\Delta\beta|}) for classifying a site as differentially methylated
#'   in ORA.  Default \code{0.1}.
#' @param score_metric Character string controlling the per-site scoring
#'   function used for GSEA ranking.  One of \code{"combined"} (default),
#'   \code{"padj"}, or \code{"delta_beta"}.  See Details.
#' @param gene_score_agg Character string; how to aggregate per-site scores to
#'   a per-gene score when multiple sites map to the same gene.  Either
#'   \code{"max"} (default; largest absolute score, sign preserved) or
#'   \code{"mean"}.
#' @param mod_type Character string or \code{NULL}; passed to
#'   \code{\link{results}} for modification-type filtering before enrichment.
#' @param mod_context Character string or \code{NULL}; passed to
#'   \code{\link{results}} for modification-context filtering.
#' @param pvalueCutoff Numeric; p-value cutoff passed to
#'   \pkg{clusterProfiler}.  Default \code{0.05}.
#' @param qvalueCutoff Numeric; q-value cutoff passed to
#'   \pkg{clusterProfiler} (ORA only).  Default \code{0.2}.
#' @param minGSSize Integer; minimum gene set size passed to
#'   \pkg{clusterProfiler}.  Default \code{10}.
#' @param maxGSSize Integer; maximum gene set size passed to
#'   \pkg{clusterProfiler}.  Default \code{500}.
#'
#' @return A named list with elements \code{$go} and \code{$kegg} (either or
#'   both may be \code{NULL} if the corresponding analysis was not requested or
#'   not possible).  When a single \code{method} is requested, each element is
#'   a \pkg{clusterProfiler} \code{enrichResult} (ORA) or \code{gseaResult}
#'   (GSEA), or \code{NULL}.  When both \code{"ora"} and \code{"gsea"} are
#'   requested, each element is itself a named list with \code{$ora} and
#'   \code{$gsea} slots.
#'
#' @seealso \code{\link{diffMethyl}}, \code{\link{annotateSites}},
#'   \code{\link{results}}, \code{\link{filterResults}}
#'
#' @examples
#' \donttest{
#' # Requires clusterProfiler and a custom TERM2GENE mapping
#' if (requireNamespace("clusterProfiler", quietly = TRUE)) {
#'   data(comma_example_data)
#'   dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
#'   ann <- annotateSites(dm, annotation(comma_example_data), type = "overlap")
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

    # ── Validate inputs before loading optional packages ──────────────────────
    method <- match.arg(method, choices = c("ora", "gsea"), several.ok = TRUE)

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

    ont         <- match.arg(ont,          c("BP", "MF", "CC", "ALL"))
    score_metric    <- match.arg(score_metric, c("combined", "padj", "delta_beta"))
    gene_score_agg  <- match.arg(gene_score_agg, c("max", "mean"))

    # ── Extract results ───────────────────────────────────────────────────────
    res_df <- results(object, mod_type = mod_type, mod_context = mod_context)

    if (!"dm_padj" %in% colnames(res_df)) {
        stop(
            "No differential methylation results found.\n",
            "Run diffMethyl() first:\n",
            "  dm <- diffMethyl(object, formula = ~ condition)"
        )
    }

    # ── Site-to-gene mapping ──────────────────────────────────────────────────
    sg <- .siteToGeneMap(res_df, gene_col)

    if (nrow(sg) == 0L) {
        warning("Gene map is empty — returning NULL for all analyses.")
        return(list(go = NULL, kegg = NULL))
    }

    # ── Run analyses ──────────────────────────────────────────────────────────
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

        sig_genes      <- unique(sg$gene_id[sig_mask])
        universe_genes <- unique(sg$gene_id)

        if (length(sig_genes) == 0L) {
            warning(
                "No significantly differentially methylated genes found ",
                "(padj <= ", padj_threshold, " and |delta_beta| >= ",
                delta_beta_threshold, ").\n",
                "ORA will not be run. Consider relaxing the thresholds."
            )
        } else {
            # GO ORA
            if (!is.null(TERM2GENE)) {
                go_ora <- clusterProfiler::enricher(
                    gene         = sig_genes,
                    universe     = universe_genes,
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
                    universe     = universe_genes,
                    pvalueCutoff = pvalueCutoff,
                    qvalueCutoff = qvalueCutoff,
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize,
                    readable     = FALSE
                )
            }

            # KEGG ORA
            if (!is.null(organism)) {
                keg_ora <- clusterProfiler::enrichKEGG(
                    gene         = sig_genes,
                    organism     = organism,
                    keyType      = keyType,
                    universe     = universe_genes,
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
            # GO GSEA
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

            # KEGG GSEA
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

    # ── Assemble return value ─────────────────────────────────────────────────
    if (both) {
        list(
            go   = list(ora = go_ora,  gsea = go_gsea),
            kegg = list(ora = keg_ora, gsea = keg_gsea)
        )
    } else if ("ora" %in% method) {
        list(go = go_ora, kegg = keg_ora)
    } else {
        list(go = go_gsea, kegg = keg_gsea)
    }
}
