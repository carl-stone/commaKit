# comma 0.3.0

## Major new features

* **`annotateSites()`** — new vectorized annotation function that replaces
  `annotateMethylSites()`, `annotateTSS()`, and `annotateTTS()`. Accepts a
  `commaData` object and uses `GenomicRanges::findOverlaps()` internally —
  no nested for-loops. Supports three annotation modes:
  - `type = "overlap"`: assigns feature type and name to each site;
    non-overlapping sites are labelled `"intergenic"`.
  - `type = "proximity"`: reports the nearest feature, its distance (in bp),
    and a signed relative position relative to the feature's TSS.
  - `type = "metagene"`: reports a fractional position (0–1) within the
    overlapping feature, strand-aware.
  Returns an updated `commaData` with new columns in `rowData`.

* **`slidingWindow()`** — new generalised sliding window function that
  replaces `methylRollingMedian()` and `methylRollingMean()`. Accepts a
  `commaData` object; genome size is always read from `genome(object)` — no
  hardcoded values. Supports `stat = "median"` or `stat = "mean"`, per-sample
  output, `mod_type` filtering, and circular genome wrap-around. Returns a
  tidy `data.frame`.

* **`methylomeSummary()`** — new function that computes per-sample summary
  statistics (mean/median/SD of beta values, fraction methylated, mean/median
  coverage). Returns a tidy `data.frame` suitable for `ggplot2` or tabular
  reporting. Supports `mod_type` filtering.

* **`coverageDepth()`** — new function that bins the genome into
  non-overlapping windows and computes mean or median sequencing depth per
  window per sample. Replaces `calculateMethylSiteDepth()`. Returns a tidy
  `data.frame`.

* **`varianceByDepth()`** — new function that computes per-sample methylation
  variance as a function of sequencing depth. Replaces `varByCoverage()`.
  Accepts a `commaData` object; no hardcoded column names. Returns a tidy
  `data.frame`.

* **`writeBED()`** — rewritten from scratch. Accepts a `commaData` object,
  an output path, and a sample name. Removes all hardcoded developer paths
  and chromosome names. Writes a standard BED9 file with an optional
  blue-to-red `itemRGB` methylation scale.

## Breaking changes

* **`annotateMethylSites()`**, **`annotateTSS()`**, and **`annotateTTS()`**
  have been removed. Use `annotateSites()` instead.

* **`methylRollingMedian()`** and **`methylRollingMean()`** have been
  removed. Use `slidingWindow()` instead.

* **`calculateMethylSiteDepth()`** and **`varByCoverage()`** have been
  removed (they were not exported). Use `coverageDepth()` and
  `varianceByDepth()` instead.

## Package changes

* Root-level clutter removed: `functions.R`, `testscript.R`,
  `WT_6mA_Mg.txt`, `WT_6mA_all_callers.txt`, `all_site_annotations.txt`,
  and `all_site_annotations_60p.txt` have been deleted.

* `methylKitGATC.R` (historical analysis script) moved to
  `inst/scripts/methylKitGATC_historical.R` with a descriptive header
  comment. It is not part of the package.

* Version bumped to `0.3.0`.

## Dependency changes

* **Added** `GenomeInfoDb` to `Imports` (was transitively available but is
  now explicitly declared).

---

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
