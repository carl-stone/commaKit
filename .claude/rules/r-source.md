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
annotateSites(object, features, type, ...)    # type = "overlap"|"proximity"|"metagene"
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
enrichMethylation(object, method, OrgDb, keyType, ont, organism, TERM2GENE, TERM2NAME,
                  gene_col, padj_threshold, delta_beta_threshold, score_metric,
                  gene_score_agg, mod_type, mod_context,
                  pvalueCutoff, qvalueCutoff, minGSSize, maxGSSize)
# Returns list(go=..., kegg=...)

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
