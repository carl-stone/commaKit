# Core Architecture

commaKit analyzes bacterial DNA methylation from Oxford Nanopore modkit pileup
output. The package is modification-type agnostic and centers on the
`commaData` S4 class.

## Terminology

- **modkit bedMethyl input** is modkit pileup bedMethyl accepted by the parser,
  not generic BED or an arbitrary TSV.
- **modkit count fields** (`Nvalid_cov`, `Nmod`, `Ncanonical`, and
  `Nother_mod`) are authoritative for coverage and methylation calculations;
  `fraction_modified` is derived.
- A **site** is a 1-bp genomic position after caller evidence has been mapped
  and aggregated. A **read position** is an offset within one sequencing read.
- A **modified-base call** becomes site evidence only after it maps to a
  reference site.
- `dm_padj` is commaKit's adjusted p-value. Backend-specific statistics retain
  backend-specific names and meanings.

## Data model

`commaData` extends `RangedSummarizedExperiment`:

- rows are 1-bp genomic methylation sites stored in `rowRanges()`;
- columns are samples described by `colData()`;
- `mod_type` is the chemical modification;
- `motif` is sequence context;
- `mod_context` combines modification and motif and is derived on demand.

`site_key` is a display label, not an alignment key. Site alignment uses
`GenomicRanges::findOverlaps()` plus explicit modification and motif matching.

## Assays and results

Core assays are `methylation`, `coverage`, `mod_counts`, `canonical_counts`,
and `other_mod_counts`. Use `siteCoverage()` for package-specific coverage;
`coverage()` is retained only for compatibility.

Derived assay layers are named and retain provenance. Differential methylation
runs are stored as named result layers, with the active result mirrored into
`dm_*` row metadata for compatibility.

## Annotation

`annotateSites()` stores every relevant association in list-columns. A site may
overlap several biological features; the implementation must not collapse this
to a single nearest or first match.

## Differential methylation

`diffMethyl()` must:

- test each `mod_context` independently;
- report effect sizes on the beta scale;
- adjust p-values over the complete testing family from the call;
- reserve `dm_padj` for commaKit's backend-independent adjusted p-value;
- keep backend-specific statistics in backend-specific columns.

Genome sizes come from `Seqinfo`/`seqlengths()`, never hardcoded organisms.
