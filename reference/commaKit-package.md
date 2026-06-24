# commaKit: Comparative Microbial Methylomics Analysis Kit

The commaKit package provides a complete toolkit for genome-wide
analysis of bacterial DNA methylation from Oxford Nanopore sequencing
data. It supports the three major modification types (6mA, 5mC, 4mC),
handles input from modkit pileup (primary), Dorado BAM, and Megalodon
(legacy) callers, and provides a unified data container, annotation
utilities, differential methylation testing, and a full set of
visualization functions.

## Value

No return value. This page provides package-level documentation. See
individual function pages for return values.

## Details

commaKit was previously developed under the package names comma and
CoMMA.

## Main workflow

1.  **Load data** —
    [`commaData()`](https://carl-stone.github.io/commaKit/reference/commaData.md)
    constructs the central `commaData` S4 object from per-sample
    methylation files.

2.  **QC** —
    [`methylomeSummary()`](https://carl-stone.github.io/commaKit/reference/methylomeSummary.md),
    [`coverageDepth()`](https://carl-stone.github.io/commaKit/reference/coverageDepth.md),
    [`plot_coverage()`](https://carl-stone.github.io/commaKit/reference/plot_coverage.md),
    [`plot_methylation_distribution()`](https://carl-stone.github.io/commaKit/reference/plot_methylation_distribution.md),
    and
    [`plot_pca()`](https://carl-stone.github.io/commaKit/reference/plot_pca.md)
    provide sample-level quality assessment.

3.  **Annotate** —
    [`loadAnnotation()`](https://carl-stone.github.io/commaKit/reference/loadAnnotation.md)
    imports a GFF3 or BED file;
    [`annotateSites()`](https://carl-stone.github.io/commaKit/reference/annotateSites.md)
    maps methylation sites to genomic features.

4.  **Visualize** —
    [`plot_genome_track()`](https://carl-stone.github.io/commaKit/reference/plot_genome_track.md)
    and
    [`plot_metagene()`](https://carl-stone.github.io/commaKit/reference/plot_metagene.md)
    show methylation in a genomic context.

5.  **Differential methylation** —
    [`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
    tests each site;
    [`results()`](https://carl-stone.github.io/commaKit/reference/results.md)
    and
    [`filterResults()`](https://carl-stone.github.io/commaKit/reference/filterResults.md)
    extract and filter the results table;
    [`plot_volcano()`](https://carl-stone.github.io/commaKit/reference/plot_volcano.md)
    and
    [`plot_heatmap()`](https://carl-stone.github.io/commaKit/reference/plot_heatmap.md)
    visualize the findings.

## Key classes and constructors

- [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md):

  The central S4 data container, extending `SummarizedExperiment`.
  Stores methylation (beta), coverage, modified-count, and
  canonical-count matrices, per-site and per-sample metadata, genome
  information, genomic annotation, and motif site locations.

- [`assayLayers`](https://carl-stone.github.io/commaKit/reference/assayLayers.md):

  Tabular registry of assay layers, provenance, and default roles.

- [`resultLayers`](https://carl-stone.github.io/commaKit/reference/resultLayers.md):

  Tabular registry of named differential methylation result runs.

## Package options

commaKit records structured breadcrumbs for parsing and construction
events when contextual error tracking is enabled. Set
`options(commaKit.error_tracking = TRUE)` or
`COMMAKIT_ERROR_TRACKING=true` to capture recent error events in
`options("commaKit.error_tracking.events")`. Set
`options(commaKit.error_tracking.reporter = function(event) ...)` to
receive event payloads in process, or set
`options(commaKit.sentry.dsn = "...")` / `SENTRY_DSN` to send captured
errors to Sentry when the suggested curl package is installed. Optional
context can be supplied through `commaKit.error_tracking.user`,
`commaKit.error_tracking.release`,
`commaKit.error_tracking.environment`, and the timeout options
`commaKit.error_tracking.timeout` and
`commaKit.error_tracking.connect_timeout`.

## References

The modkit pileup format is documented at
<https://nanoporetech.github.io/modkit/>.

## See also

Useful links:

- <https://github.com/carl-stone/commaKit>

- <https://carl-stone.github.io/commaKit/>

- Report bugs at <https://github.com/carl-stone/commaKit/issues>

## Author

**Maintainer**: Carl Stone <carl.j.stone@vanderbilt.edu>
([ORCID](https://orcid.org/0000-0003-0232-5223))

Authors:

- Carl Stone <carl.j.stone@vanderbilt.edu>
  ([ORCID](https://orcid.org/0000-0003-0232-5223))
