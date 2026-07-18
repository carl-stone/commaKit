---
type: Architecture Contract
title: commaKit Core Architecture
description: Current data structures, user-facing API groups, and behavioral contracts agents must preserve.
resource: R
tags: [architecture, r, bioconductor, commaData, mod_context, assays, results]
timestamp: 2026-06-15T00:00:00Z
status: current
owner: Carl Stone
---

# Summary

commaKit analyzes bacterial DNA methylation from Oxford Nanopore data. The package is modification-type agnostic and centers on the `commaData` S4 class.

# commaData

`commaData` extends `RangedSummarizedExperiment`.

Rows are methylation sites. Columns are samples. Genomic sites live in `rowRanges(object)` as 1-bp `GRanges`; per-site metadata lives in `mcols(rowRanges(object))`.

Required per-site metadata:

| Field | Meaning |
|---|---|
| `mod_type` | Chemical modification, currently valid values `6mA`, `5mC`, `4mC`. |
| `motif` | Sequence context, such as `GATC` or `CCWGG`; may be `NA` for callers without motif context. |

Computed accessors:

| Accessor | Meaning |
|---|---|
| `modContexts(object)` | Unique `mod_context` values. Uses `paste(mod_type, motif, sep = "_")` when motif exists, otherwise `mod_type`. |
| `siteInfo(object)` | Flat `DataFrame` with `chrom`, `position`, `strand`, metadata columns, computed `mod_context`, and computed `site_key`. |
| `sampleInfo(object)` | Per-sample metadata from `colData(object)`. |
| `genomeSizes(object)` | Chromosome sizes from `Seqinfo`. |

`site_key` is a human-readable display label with colon-separated fields, for example `chr1:512:+:6mA:GATC`. It is not an alignment key.

# Assays

Core assays:

| Assay | Accessor | Meaning |
|---|---|---|
| `methylation` | `methylation(object)` | Beta values in [0, 1], with `NA` for below-threshold sites. |
| `coverage` | `siteCoverage(object)` | Non-negative integer-like read depth. |
| `mod_counts` | `modCounts(object)` | Reads called as the target modification, when available. |
| `canonical_counts` | `canonicalCounts(object)` | Reads called canonical/unmodified, when available. |
| `other_mod_counts` | `otherModCounts(object)` | Reads called as non-target modifications, when available. |

`coverage(object)` is deprecated for `commaData`; use `siteCoverage(object)`.

Assay layer metadata lives in the named
`metadata(object)$assay_provenance` record list and named
`metadata(object)$assay_defaults` role map. An explicit `default_for` in a
provenance record may fill a missing role-map entry when loading older objects.
Assay names alone do not infer layer metadata or defaults; absent metadata is
reported as `NA`. Use `assayLayers(object)` and `assayProvenance(object)`
rather than ad hoc metadata inspection in user-facing code.

# Result Layers

`diffMethyl()` stores named result layers in `metadata(object)$diffMethyl_results`, with provenance in `metadata(object)$diffMethyl_result_layers` and an active/default result in `metadata(object)$diffMethyl_default_result`.

The supported result shape is the three named registries above: result tables,
layer records, and the explicit active layer name. The active result is
mirrored into bare `dm_*` row metadata columns for backward compatibility. The
older `diffMethyl_result_cols` and `diffMethyl_params` fields, inferred
default names, and row metadata without a named layer are not used to
reconstruct results. Use `resultLayers(object)` to inspect named runs.
`results()` and `filterResults()` can retrieve named result layers.

# Pipeline

1. `commaData()` constructs an object from modkit, Dorado, or Megalodon output.
2. `annotateSites()` annotates sites to genomic features and stores all associations as list-columns.
3. `diffMethyl()` tests differential methylation by `mod_context`.
4. `results()` and `filterResults()` extract and filter result layers.
5. `enrichMethylation()` runs ORA/GSEA over target/regulator genes.
6. QC and diagnostics include `methylomeSummary()`, `coverageDepth()`, `varianceByDepth()`, `mValues()`, and plots.
7. Regional/export utilities include `slidingWindow()`, `summarizeRegions()`, and `writeBED()`.

# Parser Contracts

- modkit pileup BED is the primary path.
- Dorado BAM reads MM/ML tags and can contain 6mA, 5mC, and 4mC in one file.
- Megalodon is legacy and requires explicit `mod_type`.
- Constructor merging aligns by `GRanges` overlaps plus `mod_type`/motif checks, not by row names or string keys.

# Annotation Contract

`annotateSites()` stores all feature associations per site in list-columns:

- `feature_types`
- `feature_names`
- `rel_position`
- `frac_position`

Metadata passthrough columns become `{column}_values` list-columns. Intergenic/no-hit sites use length-0 list elements. Do not replace this with `distanceToNearest()` or first-hit behavior.

# Differential Methylation Contract

- Loop by `mod_context`, not `mod_type`.
- Report effect sizes on the beta scale.
- Adjust p-values genome-wide across all contexts tested in that call.
- Keep `dm_padj` as commaKit's backend-independent adjusted p-value. methylKit
  q-values are backend-specific evidence and, when present, live in
  `dm_methylkit_qvalue` only on methylKit result layers.
- The default backend is `methylkit` for compatibility with established methylKit workflows. `quasi_f` is the preferred package-native alternative when users want count-aware empirical-Bayes shrinkage and genome-wide multiple-testing correction handled inside commaKit.
- Formula support is one-sided with 2+ levels; multi-level formulas must fail clearly.

# User-Facing API Groups

Core and accessors:

```r
commaData()
methylation(); siteCoverage(); modCounts(); canonicalCounts(); otherModCounts()
sampleInfo(); siteInfo(); modTypes(); modContexts(); motifs()
genomeSizes(); annotation(); motifSites(); caller(); minCoverage()
assayLayers(); assayProvenance(); resultLayers()
filterSites()
```

Analysis, enrichment, and export:

```r
annotateSites(); diffMethyl(); results(); filterResults()
buildKEGGGeneIDMap(); buildKEGGTermGene(); enrichMethylation()
methylomeSummary(); coverageDepth(); varianceByDepth()
mValues(); slidingWindow(); summarizeRegions(); writeBED()
loadAnnotation(); findMotifSites()
```

Plots:

```r
plot_coverage(); plot_methylation_distribution(); plot_pca()
plot_genome_track(); plot_metagene(); plot_tss_profile()
plot_volcano(); plot_heatmap()
```

# See Also

- [Design Decisions](design-decisions.md)
- [Known Issues](known-issues.md)
- [Project Status](project-status.md)
