# comma 0.2.0

## Major new features

* Added `commaData` S4 class (extends `SummarizedExperiment`) — the central data
  object for all comma analyses. Stores methylation beta values and coverage as
  assay matrices with per-site and per-sample metadata.

* Added `commaData()` constructor: parses one or more modkit pileup BED files
  (or Megalodon files) and merges them into a `commaData` object with consistent
  site × sample matrices. Supports multi-sample, multi-modification-type data.

* Added modkit parser (`.parseModkit()`): reads modkit `pileup` BED output
  (primary ONT methylation calling format). Supports all three modification
  types: 6mA (`a`), 5mC (`m`), and 4mC (`21839`).

* Added Megalodon parser (`.parseMegalodon()`): reads legacy Megalodon per-read
  BED output for backward compatibility with existing datasets.

* Added Dorado BAM parser stub (`.parseDorado()`): not yet implemented;
  provides a clear error message directing users to use `modkit pileup` first.

* Added `findMotifSites()`: locates all instances of a DNA sequence motif (e.g.,
  GATC, CCWGG) in a genome (BSgenome object or FASTA file) and returns a
  `GRanges` of all positions on both strands.

* Added `loadAnnotation()`: reads GFF3 or BED annotation files and returns a
  `GRanges` with standardized `feature_type` and `name` metadata columns.

* Added accessor functions: `methylation()`, `coverage()`, `sampleInfo()`,
  `siteInfo()`, `modTypes()`, `genome()`, `annotation()`, `motifSites()`.

* Added subsetting: `[` for sites/samples and `subset()` for filtering by
  `mod_type`, `condition`, or `chrom`.

* Added internal genome utilities: `.validateGenomeInfo()`, `.circularIndex()`,
  `.makeSeqinfo()`.

## Data changes

* **Removed** all 9 MG1655-specific bundled datasets (`data/*.rda`). These were
  specific to *E. coli* K-12 MG1655 and not usable with other organisms.

* **Added** `comma_example_data`: a synthetic `commaData` object with a
  simulated 100 kb genome, 3 samples (2 conditions), and two modification types
  (6mA and 5mC). Used in vignettes and tests. Run
  `data-raw/create_example_data.R` to regenerate.

* **Added** `inst/extdata/example_modkit.bed`: small modkit pileup BED for
  constructor tests and examples.

* **Added** `inst/extdata/example.gff3`: matching GFF3 annotation for the
  example genome.

## Dependency changes

* **Added** to `Imports`: `GenomicRanges`, `SummarizedExperiment`, `IRanges`,
  `S4Vectors`, `BiocGenerics`, `Rsamtools`, `zoo`

* **Added** to `Suggests`: `BSgenome`, `Biostrings`, `rtracklayer`

* **Removed** from `Imports`: `forcats`

## Package changes

* Package renamed from `CoMMA` to `comma` (Comparative Methylomics for
  Microbial Analysis). The new name reflects the broader scope — all
  modification types (6mA, 5mC, 4mC), not just adenine methylation.

* Version bumped to `0.2.0`.

* `R (>= 4.1.0)` is now the minimum R version requirement.

* Added `biocViews` field in preparation for Bioconductor submission.

## Notes

* The existing functions `annotateMethylSites()`, `annotateTSS()`, and
  `methylRollingMedian()` are preserved but unchanged. They do not yet accept
  `commaData` objects and retain their known limitations (nested for-loops,
  hardcoded genome size). These will be replaced in Phase 3 (v0.3.0).

---

# comma 0.1.0

* Initial release as `CoMMA` (Comparison of Microbial Methylated Adenines).
* Exported functions: `annotateMethylSites()`, `annotateTSS()`,
  `methylRollingMedian()`.
* Data: 9 MG1655-specific `.rda` datasets.
