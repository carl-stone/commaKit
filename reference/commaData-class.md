# commaData: the central data object for the commaKit package

`commaData` is an S4 class that extends
[`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
to store genome-wide bacterial methylation data from Oxford Nanopore
sequencing. It is the central object accepted and returned by all
commaKit analysis functions.

## Value

An object of class `commaData`. Use
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
to construct instances.

## Details

Genomic annotation and motif site positions are stored in
`metadata(object)` rather than as dedicated slots. Use
[`annotation`](https://rdrr.io/pkg/BiocGenerics/man/annotation.html)`(object)`
and
[`motifSites`](https://carl-stone.github.io/commaKit/reference/motifSites.md)`(object)`
to access them.

Genome size information is stored in the `Seqinfo` attached to
`rowRanges(object)`, accessible via `seqlengths(object)` or
`seqinfo(object)`. Use
[`genomeSizes`](https://carl-stone.github.io/commaKit/reference/genomeSizes.md)
for chromosome sizes; `genome()` is retained only as a
backward-compatible size-vector method.

The class stores methylation data in assay matrices (accessible via
[`assay`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)):

- `"methylation"`:

  Beta values (proportion of reads called methylated, range 0-1). Sites
  with coverage below the `min_coverage` threshold are stored as `NA`.

- `"coverage"`:

  Integer read depth at each site.

- `"mod_counts"`:

  Observed reads called as the target modification by modkit.

- `"canonical_counts"`:

  Observed reads called canonical or unmodified by modkit.

- `"other_mod_counts"`:

  Observed reads called as non-target modifications by modkit. For
  modkit pileup, `coverage` is the valid denominator
  `mod_counts + canonical_counts + other_mod_counts`.

Genomic positions are stored in `rowRanges(object)`, a
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
with one 1-bp range per methylation site. Per-site metadata is in the
`mcols` of this GRanges and includes at minimum: `mod_type` and `motif`.
The `mod_type` column is a factor with levels `c("4mC", "5mC", "6mA")`,
enforcing valid values at the data structure level. The `motif` column
stores the sequence context of each site (e.g., `"GATC"` or `"CCWGG"`)
as extracted from the modkit `mod_code` field. The `mod_context` is
computed on demand from `mod_type` and `motif` (e.g., `"6mA_GATC"`,
`"5mC_CCWGG"`), or just `mod_type` when motif is unavailable. Use
[`modContexts`](https://carl-stone.github.io/commaKit/reference/modContexts.md)`(object)`
or
[`siteInfo`](https://carl-stone.github.io/commaKit/reference/siteInfo.md)`(object)`
to retrieve it. All analyses default to running independently per
`mod_context` group to prevent spurious mixing of biologically distinct
methylation events.

For convenience,
[`siteInfo`](https://carl-stone.github.io/commaKit/reference/siteInfo.md)`(object)`
returns a flat `DataFrame` combining the genomic coordinates (chrom,
position, strand) with the mcols columns.

Per-sample metadata is in `colData(object)` and includes at minimum:
`sample_name` and `replicate`. A `condition` column is optional
container metadata; functions that need a grouping/design variable (such
as
[`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md))
validate that requirement locally.

The minimum coverage threshold is stored in `metadata(object)` and
accessible via
[`minCoverage`](https://carl-stone.github.io/commaKit/reference/minCoverage.md)`(object)`.

Assay-layer provenance and defaults are stored in
`metadata(object)$assay_provenance` and
`metadata(object)$assay_defaults`. Use
[`assayLayers`](https://carl-stone.github.io/commaKit/reference/assayLayers.md)
for a tabular summary.

Differential methylation result layers are stored in
`metadata(object)$diffMethyl_results` with provenance in
`metadata(object)$diffMethyl_result_layers`. Use
[`resultLayers`](https://carl-stone.github.io/commaKit/reference/resultLayers.md)
to list named result runs.

## See also

[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
for the constructor,
[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md),
[`siteCoverage`](https://carl-stone.github.io/commaKit/reference/siteCoverage.md),
[`modCounts`](https://carl-stone.github.io/commaKit/reference/modCounts.md),
[`canonicalCounts`](https://carl-stone.github.io/commaKit/reference/canonicalCounts.md),
[`assayLayers`](https://carl-stone.github.io/commaKit/reference/assayLayers.md),
[`assayProvenance`](https://carl-stone.github.io/commaKit/reference/assayProvenance.md),
[`resultLayers`](https://carl-stone.github.io/commaKit/reference/resultLayers.md),
[`sampleInfo`](https://carl-stone.github.io/commaKit/reference/sampleInfo.md),
[`siteInfo`](https://carl-stone.github.io/commaKit/reference/siteInfo.md),
[`modTypes`](https://carl-stone.github.io/commaKit/reference/modTypes.md),
[`modContexts`](https://carl-stone.github.io/commaKit/reference/modContexts.md),
[`annotation`](https://rdrr.io/pkg/BiocGenerics/man/annotation.html) for
accessors.
