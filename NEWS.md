# comma 0.5.0

## Major new features

* **`plot_methylation_distribution()`** — beta value density plot per sample,
  coloured by sample name, and faceted by modification type when multiple types
  are present. Useful for QC and comparing methylation level distributions.

* **`plot_genome_track()`** — genome browser-style scatter plot of methylation
  beta values vs. genomic position. Supports positional windowing via `start`/
  `end` arguments, `mod_type` filtering, and optional feature annotation
  rectangles from `annotation(object)`.

* **`plot_metagene()`** — average methylation profile across a class of genomic
  features, normalized to fractional position [0 = TSS, 1 = TTS]. Uses
  `annotateSites(type = "metagene")` internally.

* **`plot_volcano()`** — volcano plot for differential methylation results.
  Accepts the `data.frame` output of `results()`. Colors points as
  Hypermethylated, Hypomethylated, or Not significant based on user-supplied
  `padj_threshold` and `delta_beta_threshold`.

* **`plot_heatmap()`** — `geom_tile` heatmap of the top differentially
  methylated sites (ranked by adjusted p-value) across all samples. Uses
  `ggplot2` only; no `ComplexHeatmap` dependency required.

* **`plot_pca()`** — PCA of per-sample methylation profiles using
  `stats::prcomp()`. Points colored and optionally shaped by any column in
  `sampleInfo(object)`.

* **`plot_coverage()`** — histogram of sequencing depth per site, per sample.
  Useful for coverage QC before differential methylation testing.

## Documentation

* Added **two vignettes** required for Bioconductor submission:
  - `"Getting Started with comma"` — end-to-end workflow using
    `comma_example_data`.
  - `"Working with Multiple Modification Types"` — joint 6mA + 5mC analysis.

* Added **package-level documentation** (`?comma`) describing the overall
  workflow.

## Dependency changes

* Added `BiocStyle`, `ComplexHeatmap`, `ggrepel`, and `patchwork` to
  `Suggests`. None are required for core functions; `patchwork` enables
  combined annotation tracks in `plot_genome_track()` and `plot_heatmap()`.

* Added `VignetteBuilder: knitr` to DESCRIPTION.

* Added `Visualization` to `biocViews`.

# comma 0.4.0

## Major new features

* **`diffMethyl()`** — the core differential methylation function. Accepts a
  `commaData` object and a design formula, tests each methylation site for
  differential methylation between conditions, and returns the object enriched
  with per-site statistics in `rowData`. Supports two statistical backends:
  - `method = "beta_binomial"` (default): per-site quasibinomial GLM via base
    R `stats::glm()`. No extra packages required.
  - `method = "methylkit"`: wraps `methylKit::calculateDiffMeth()` for users
    who prefer the methylKit pipeline. Requires `methylKit` to be installed.

  New `rowData` columns: `dm_pvalue`, `dm_padj` (Benjamini-Hochberg by
  default), `dm_delta_beta` (effect size: treatment minus control mean beta),
  and one `dm_mean_beta_<condition>` column per condition level. Analysis
  parameters are stored in `metadata(object)$diffMethyl_params`.

* **`results()`** — S4 method to extract differential methylation results from
  a `commaData` object as a tidy `data.frame`. Supports optional `mod_type`
  filtering.

* **`filterResults()`** — S4 method to filter the results table by adjusted
  p-value (`padj`) and absolute effect size (`delta_beta`) thresholds.

* **`.parseDorado()` — full Dorado BAM parser** replacing the previous stub.
  Reads MM/ML base modification tags from Dorado-aligned BAM files using
  `Rsamtools::scanBam()`, parses modification positions from CIGAR-decoded read
  coordinates, and aggregates per-read calls into per-site beta values.
  Supports 6mA, 5mC, and 4mC in the same BAM file. Invoked automatically by
  `commaData(..., caller = "dorado")`.

## Dependency changes

* **Added** `methylKit` to `Suggests`. Required only for
  `diffMethyl(..., method = "methylkit")`.

---

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
