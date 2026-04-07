---
paths:
  - "R/**/*.R"
  - "data-raw/**/*.R"
---

# R Source — Architecture, Dependencies, API

## commaData S4 Class

Every analysis function accepts a `commaData` object (modeled on DESeq2's `DESeqDataSet`):

```
commaData
├── methylation    # sites × samples matrix of beta values (0–1)
├── coverage       # sites × samples matrix of read depth
├── rowData        # per-site: chrom, position, strand, motif, mod_type, mod_context
├── colData        # per-sample: sample_name, condition, replicate, caller, file_path
├── genomeInfo     # chromosome names and sizes (named integer vector or NULL)
├── annotation     # GRanges of genomic features (from GFF3/BED)
├── motifSites     # GRanges of all motif instances in genome
└── metadata       # list: package version, creation date, user fields
```

`rowData` includes `mod_type` (`"6mA"`, `"5mC"`, `"4mC"`) and `mod_context` (e.g. `"6mA_GATC"`, `"5mC_CCWGG"`) as first-class columns. `mod_context` is required in all objects (v0.8.0+); old objects without it will fail `validObject()`.

### Constructor

```r
commaData(
  files,                  # named character vector: sample_name → file_path
  colData,                # data frame: sample_name, condition, replicate (minimum)
  genome,                 # BSgenome, FASTA path, or named integer vector of chr sizes
  annotation,             # GFF3 path or GRanges (optional)
  mod_type,               # "6mA", "5mC", "4mC", or NULL to auto-detect
  motif,                  # regex motif string (e.g., "GATC") or NULL
  min_coverage,           # integer, default 5
  caller,                 # "dorado", "modkit", "megalodon"
  expected_mod_contexts   # named list, e.g. list("6mA"=c("GATC"), "5mC"="CCWGG")
                          # drops sites with unexpected mod_type × motif combinations
)
```

### annotateSites() list-column design (do not revert)

`annotateSites()` stores **all** overlapping/nearby features per site as `CharacterList`/`IntegerList`/`NumericList` columns in `rowData` — not just the first or closest. This reflects the highly overlapping nature of bacterial genome annotations (genes, promoters, TF binding sites). Do NOT revert to single-match (`!duplicated()` or `distanceToNearest()`) behavior. Intergenic sites receive length-0 list elements; test with `lengths(col) == 0`.

Four parallel list columns are always computed by the unified backend (`.annotateSites_unified()`):
- `feature_types`: GFF3 type per association (`CharacterList`)
- `feature_names`: feature name per association (`CharacterList`)
- `rel_position`: signed distance from feature; 0 = inside, negative = upstream, positive = downstream, strand-aware (`IntegerList`)
- `frac_position`: normalized position [0,1] inside feature (TSS=0); `NA` outside (`NumericList`)

Additional metadata columns can be passed through via `metadata_cols` — stored as `{col}_values` CharacterList columns, parallel to the above.

`keep` post-filter:
- `"all"` — all four columns
- `"overlap"` — subset to `rel_position==0`; drops `rel_position`/`frac_position`
- `"proximity"` — keep all associations; drops `frac_position`
- `"metagene"` — subset to `rel_position==0`; keeps `frac_position`

### loadAnnotation() — feature_subtype preservation

EcoCyc GFF3 files encode sigma factor identity as a `feature_type=Sigma70` attribute in GFF3 column 9. `.standardizeAnnotationMcols()` now saves this to `feature_subtype` **before** overwriting `feature_type` with the GFF3 type column. Result: `feature_subtype` in mcols holds sigma factor/binding-protein identity (e.g., `"Sigma70"`, `"DnaA-binding site"`); `feature_type` holds the SO type (e.g., `"transcription_factor_binding_site"`). The `transcription_unit` attribute (for binding sites encoding regulated operons) is preserved as-is by `rtracklayer::import()`.

---

## Dependencies

### Hard imports (`Imports` in DESCRIPTION)

| Package | Purpose |
|---|---|
| `GenomicRanges` | Core genomic interval arithmetic; `findOverlaps()` for annotation |
| `GenomeInfoDb` | Chromosome/genome metadata management |
| `SummarizedExperiment` | Base class infrastructure for `commaData` |
| `IRanges` | Range operations (via GenomicRanges) |
| `S4Vectors` | DataFrame and other S4 infrastructure |
| `BiocGenerics` | Bioconductor generic methods |
| `Rsamtools` | BAM parsing for Dorado (MM/ML tags) |
| `zoo` | Rolling window in `slidingWindow()` |
| `ggplot2` | All visualization |
| `methods` | S4 class system |
| `stats` | `prcomp`, `glm`, `loess`, etc. |
| `utils` | Base R utilities |

`dplyr` and `tidyr` are **not** in `Imports`. If needed, import individually (never `tidyverse`).

### Soft dependencies (`Suggests` in DESCRIPTION)

| Package | Purpose |
|---|---|
| `BSgenome` | Genome sequence for `findMotifSites()` |
| `clusterProfiler` | GO/KEGG enrichment in `enrichMethylation()` |
| `KEGGREST` | Bulk KEGG API calls in `buildKEGGTermGene()` and `buildKEGGGeneIDMap()` |
| `AnnotationDbi` | OrgDb symbol lookup in `buildKEGGGeneIDMap()` (Bioconductor core) |
| `Biostrings` | Sequence pattern matching for motif search |
| `BiocStyle` | Vignette styling |
| `limma` | eBayes for `diffMethyl(method="limma"\|"quasi_f")` |
| `methylKit` | Alternative DM backend for `diffMethyl(method="methylkit")` |
| `patchwork` | Multi-panel plots in `plot_genome_track()`, `plot_heatmap()` |
| `rtracklayer` | GFF3 import via `import()` |
| `testthat` | Testing framework (edition 3) |
| `knitr` / `rmarkdown` | Vignette rendering |
| `ggrepel` | Available for volcano labels (not currently used) |
| `ComplexHeatmap` | Available as heatmap backend (not currently used) |
| `scales` | Axis/color scale helpers |

---

## Input Format Reference

### modkit BED (primary format)

15-column modkit `pileup` output:
```
chrom, start, end, mod_code, score, strand, coverage, mod_frequency,
n_mod, n_canonical, n_other_mod, n_delete, n_fail, n_diff, n_no_call
```

`mod_code` mapping: `a` → 6mA, `m` → 5mC, `21839` → 4mC

`mod_frequency` = beta value (0–1); `coverage` = total read depth.
Coordinates are **0-based** — parser converts: `position = start + 1`.

### Dorado BAM

MM/ML tags parsed via `Rsamtools::scanBam()`. CIGAR-decoded to reference positions. Aggregates per-read calls to per-site beta values. Handles 6mA, 5mC, 4mC in one BAM. Invoked by `caller = "dorado"`. Slower than modkit; prefer running `modkit pileup` first.

### Megalodon (legacy)

Per-read aggregation to per-site beta values. `mod_type` must be provided explicitly (cannot be inferred).

---

## Exported API Quick Reference

> **Keep this list current.** Add or remove entries here whenever an exported function is added, renamed, or removed.

```r
# S4 class + constructor
commaData(files, colData, genome, annotation, mod_type, motif, min_coverage, caller,
          expected_mod_contexts)

# Accessors
methylation(object)      # → sites × samples beta matrix
coverage(object)         # → sites × samples integer matrix
sampleInfo(object)       # → per-sample DataFrame
siteInfo(object)         # → per-site DataFrame (chrom, position, strand, mod_type, mod_context, ...)
modTypes(object)         # → character vector of modification types
modContexts(object)      # → sorted unique mod_context strings (e.g. "6mA_GATC")
motifs(object)           # → character vector of unique sequence context motifs
genome(object)           # → named integer vector of chromosome sizes
annotation(object)       # → GRanges of genomic features
motifSites(object)       # → GRanges of motif instances

# Subsetting
object[sites, samples]
subset(object, mod_type, mod_context, condition, chrom)

# Utilities
loadAnnotation(file, feature_types)   # GFF3/BED → GRanges
findMotifSites(genome, motif)         # genome + motif → GRanges

# Analysis
annotateSites(object, features, feature_col, name_col, window, keep, metadata_cols)
# keep = "all"|"overlap"|"proximity"|"metagene" (default "all")
# always produces: feature_types, feature_names, rel_position, frac_position
# metadata_cols: GFF3 attribute names to pass through as {col}_values CharacterList columns
slidingWindow(object, window, stat, mod_context, ...)
methylomeSummary(object, mod_type, mod_context)
coverageDepth(object, window, method, ...)
varianceByDepth(object, coverage_bins)
writeBED(object, file, sample, mod_context, ...)
mValues(object, alpha, mod_type, mod_context)

# Differential methylation
diffMethyl(object, formula, method, mod_type, mod_context, min_coverage, alpha, p_adjust_method)
# method: "beta_binomial" | "quasi_f" | "limma" | "methylkit"
results(object, mod_type, mod_context)
filterResults(object, padj, delta_beta, ...)

# Enrichment (requires annotateSites() + diffMethyl() first)
buildKEGGGeneIDMap(organism, OrgDb, entrez2symbol, keys_col, id_col, file)
# organism: KEGG org code; OrgDb or entrez2symbol required (at least one)
# Returns: data.frame(symbol, kegg_id) — gene symbol <-> KEGG ID mapping
# Uses 1 API call (keggConv); cache to avoid repeated network use.
# Pass to buildKEGGTermGene(id_map=...) to translate b-numbers → symbols.

buildKEGGTermGene(organism, file, strip_prefix, id_map)
# organism: KEGG org code (e.g. "eco"); file: optional RDS cache path
# id_map: optional data.frame(symbol, kegg_id) from buildKEGGGeneIDMap()
#         when provided, translates gene column from KEGG IDs to symbols
# Returns: list(term2gene = data.frame(term, gene),
#               term2name = data.frame(term, name))
# Uses 2 API calls total (keggLink + keggList); cache to avoid repeated network use.
# Pass result to enrichMethylation() via kegg_term2gene/kegg_term2name.

enrichMethylation(object, method, OrgDb, keyType, ont, organism,
                  TERM2GENE, TERM2NAME,
                  kegg_term2gene, kegg_term2name,
                  gene_col, feature_type, gene_role, overlap_only,
                  padj_threshold, delta_beta_threshold, score_metric,
                  gene_score_agg, mod_type, mod_context,
                  pvalueCutoff, qvalueCutoff, minGSSize, maxGSSize)
# object: commaData OR data.frame from results()
# feature_type: character vector (run separately per type); NULL = all
# gene_role: "target"|"regulator"|"both" (see below)
# kegg_term2gene: pre-built KEGG mapping (avoids live API rate limits)
#   → routes results to $kegg slot via enricher()/GSEA()
# organism: live KEGG API path (one req/pathway; may hit rate limit)
#   → ignored when kegg_term2gene is provided
# Returns: list(go=..., kegg=...)   [single feature_type, gene_role != "both"]
#          list(gene=list(go,kegg), promoter=...) [multiple feature_types]
#          list(target=..., regulator=...) [gene_role="both"]
#
# gene_role semantics:
#   "target"    → universe = all target genes in data
#   "regulator" → universe = only regulators of that type in annotation
#   "both"      → runs both separately, returns named sub-list
#
# .COMMA_REGION_FEATURE_TYPES: c("gene","CDS","mRNA","tRNA","rRNA","ncRNA",
#   "operon","repeat_region","prophage","region","insertion_sequence")
#   overlap_only auto-defaults to TRUE for these types
#
# .SIGMA_FACTOR_GENE_MAP: Sigma70→rpoD, Sigma24→rpoE, Sigma32→rpoH,
#   Sigma28→fliA, Sigma38→rpoS, Sigma54→rpoN, Sigma19→fecI

# Visualization (all return ggplot/patchwork; all accept mod_context filter)
plot_methylation_distribution(object, mod_type, mod_context, per_sample)
plot_genome_track(object, chromosome, start, end, mod_type, mod_context)
plot_metagene(object, feature, mod_type, mod_context, window)
plot_volcano(results_df, delta_beta_threshold, padj_threshold)
plot_heatmap(object, result_df, n_sites, annotation_cols)
plot_pca(object, mod_type, mod_context, color_by, shape_by, return_data)
plot_coverage(object, mod_context, per_sample)
plot_tss_profile(object, feature_type, window, mod_type, mod_context, motif,
                 color_by, facet_by, alpha, show_smooth, smooth_span,
                 regulatory_feature_types)
```
