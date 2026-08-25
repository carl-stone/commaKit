# Design Decisions

This file records rationale that is easy to lose by reading code alone.

## Modification context defines the testing strata

Modification types can occur in biologically distinct motifs. Differential
methylation therefore operates by `mod_context`, not merely by `mod_type`.
Multiple-testing correction still covers all sites tested in the call.

## Effect sizes are reported on the beta scale

Backends may fit transformed values internally, but `dm_delta_beta` represents
the difference in methylation proportion. This keeps results interpretable as
percentage-point changes.

## Annotation preserves multiple associations

Bacterial genomes are densely annotated. A methylation site may overlap a gene,
promoter, binding site, and operon. `annotateSites()` therefore uses list-columns
rather than selecting one association.

## Genomic ranges are authoritative

`commaData` extends `RangedSummarizedExperiment`, genomic coordinates live in
`rowRanges()`, and site matching uses genomic overlaps. Row names and formatted
keys are not reliable alignment mechanisms.

## Derived values are not duplicated unnecessarily

`mod_context` and `site_key` are computed from authoritative fields. This avoids
stale derived metadata when modification or motif information changes.

## Raw and derived analyses can coexist

Raw assays remain stable. Derived assay layers and differential-methylation
result layers are named and carry provenance so several transformations or
analyses can coexist without overwriting the underlying observations.

## Adjusted p-values have one package-level meaning

`dm_padj` is computed by commaKit over the complete testing family. Statistics
whose definitions belong to a backend, such as methylKit q-values, remain in
backend-specific columns.

## Enrichment distinguishes biological roles

Target genes and regulator genes answer different biological questions and use
different background universes. `enrichMethylation()` keeps those roles
separate.

## Network-backed annotation should be bulk and cacheable

KEGG mappings are built with bulk requests and can be cached locally. Repeated
per-pathway requests are slow and vulnerable to rate limits.
