# CLAUDE.md — AI Assistant Guide for `comma`

> This file is the primary reference for AI coding agents working on this repository. Read it fully before making any changes.

---

## 1. Project Overview

**Package name:** `comma` (was `CoMMA`; rename complete in package internals)
**Full name:** Comparative Methylomics for Microbial Analysis
**Author:** Carl Stone, Vanderbilt University (carl.j.stone@vanderbilt.edu)
**Current version:** 0.5.0
**License:** MIT
**Target:** Bioconductor submission at v1.0.0

`comma` is an R package for analyzing bacterial DNA methylation from Oxford Nanopore sequencing data. It characterizes genome-wide methylation patterns, annotates methylation sites relative to genomic features, identifies differentially methylated sites between conditions, and visualizes results. Phases 1–5 of the architectural refactor are complete; Bioconductor submission prep is the current priority.

### Scientific scope

`comma` is **modification-type agnostic**. It was originally conceived for N6-methyladenine (6mA, Dam methyltransferase, GATC motif in *E. coli*), but 4mC and 5mC have significant and understudied roles even in model bacteria, and Dorado now detects all three modification types simultaneously from a single sequencing run. Every data structure, function signature, and analysis module must treat methylation type as a first-class parameter, not an assumption. The package must work equally well for 6mA, 4mC (various bacterial MTases), 5mC (Dcm and analogs), and any other modification type that nanopore callers produce.

---

## 2. Current Repository State

### What exists now (v0.5.0 — Phases 1–5 complete)

```
CoMMA/
├── R/
│   ├── comma-package.R          # ✅ Package-level ?comma documentation (~49 lines)
│   ├── commaData_class.R        # ✅ S4 class def (~174 lines), show(), validity()
│   ├── commaData_constructor.R  # ✅ commaData() constructor (~258 lines)
│   ├── accessors.R              # ✅ methylation(), coverage(), sampleInfo(), etc. (~298 lines)
│   ├── parse_modkit.R           # ✅ PRIMARY input format parser (~143 lines)
│   ├── parse_dorado.R           # ✅ Full Dorado BAM parser with MM/ML tags (~385 lines)
│   ├── parse_megalodon.R        # ✅ Backward compatibility parser (~112 lines)
│   ├── load_annotation.R        # ✅ GFF3/BED → GRanges (~142 lines)
│   ├── find_motif_sites.R       # ✅ FASTA + motif regex → GRanges (~135 lines)
│   ├── genome_utils.R           # ✅ Genome validation, circular arithmetic (~105 lines)
│   ├── comma_example_data.R     # ✅ Roxygen docs for synthetic dataset (~57 lines)
│   ├── annotateSites.R          # ✅ Vectorized annotation using findOverlaps() (~278 lines)
│   ├── sliding_window.R         # ✅ Generalized slidingWindow() (~176 lines)
│   ├── methylome_summary.R      # ✅ Per-sample QC stats (~106 lines)
│   ├── coverage_analysis.R      # ✅ coverageDepth() + varianceByDepth() (~205 lines)
│   ├── writeBED.R               # ✅ Fully generalized BED9 export (~149 lines)
│   ├── diffMethyl.R             # ✅ Main differential methylation interface (~256 lines)
│   ├── beta_binomial.R          # ✅ Quasibinomial GLM per-site engine (~198 lines)
│   ├── methylkit_wrapper.R      # ✅ methylKit alternative method wrapper (~185 lines)
│   ├── multiple_testing.R       # ✅ BH/FDR correction utility (~25 lines)
│   ├── results_methods.R        # ✅ results() + filterResults() S4 methods (~139 lines)
│   ├── plot_distribution.R      # ✅ Phase 5: plot_methylation_distribution() (~126 lines)
│   ├── plot_genome_track.R      # ✅ Phase 5: plot_genome_track() (~269 lines)
│   ├── plot_metagene.R          # ✅ Phase 5: plot_metagene() (~195 lines)
│   ├── plot_volcano.R           # ✅ Phase 5: plot_volcano() (~134 lines)
│   ├── plot_heatmap.R           # ✅ Phase 5: plot_heatmap() (~207 lines)
│   ├── m_values.R               # ✅ mValues() — M-value transformation from beta + coverage
│   ├── plot_pca.R               # ✅ Phase 5: plot_pca() with M-value transform + return_data (~196 lines)
│   └── plot_coverage.R          # ✅ Phase 5: plot_coverage() (~147 lines)
├── vignettes/
│   ├── getting-started.Rmd           # ✅ End-to-end workflow (~214 lines)
│   └── multiple-modification-types.Rmd  # ✅ Joint 6mA + 5mC analysis (~173 lines)
├── data/
│   └── comma_example_data.rda   # ✅ Synthetic commaData (300 sites, 3 samples, chr_sim 100kb)
├── data-raw/
│   └── create_example_data.R    # ✅ Generates comma_example_data (set.seed(42))
├── inst/
│   ├── extdata/
│   │   ├── example_modkit.bed   # ✅ 20-site modkit pileup example (chr_sim, mixed 6mA+5mC)
│   │   └── example.gff3         # ✅ 5-gene GFF3 annotation (chr_sim)
│   └── scripts/
│       └── methylKitGATC_historical.R  # ✅ Historical reference (moved from root in Phase 3)
├── tests/testthat/
│   ├── test-commaData.R               # ✅ ~20 tests: S4 class, validity, constructor, show()
│   ├── test-parsers.R                 # ✅ ~15 tests: .parseModkit(), coverage filter
│   ├── test-accessors.R               # ✅ ~20 tests: all accessor methods, subsetting
│   ├── test-genome_utils.R            # ✅ Tests for .validateGenomeInfo(), .circularIndex(), .makeSeqinfo()
│   ├── test-load_annotation.R         # ✅ Tests for loadAnnotation() GFF3/BED parsing
│   ├── test-find_motif_sites.R        # ✅ Tests for findMotifSites()
│   ├── test-parse_megalodon.R         # ✅ Tests for .parseMegalodon()
│   ├── test-annotateSites.R           # ✅ ~20 tests for annotateSites() (overlap/proximity/metagene)
│   ├── test-slidingWindow.R           # ✅ ~15 tests for slidingWindow()
│   ├── test-methylomeSummary.R        # ✅ ~11 tests for methylomeSummary() (incl. all-NA sample)
│   ├── test-coverageAnalysis.R        # ✅ ~8 tests for coverageDepth() and varianceByDepth()
│   ├── test-writeBED.R                # ✅ ~20 tests for writeBED() (file creation, BED format, RGB, errors)
│   ├── test-diffMethyl.R              # ✅ ~30 tests for diffMethyl(), .betaBinomialTest(), .applyMultipleTesting(), ground-truth recovery
│   ├── test-results.R                 # ✅ ~23 tests for results() and filterResults() (incl. boundary cases)
│   ├── test-parse_dorado.R            # ✅ ~21 tests for .cigarToRefPos(), .parseMmTag(), .parseDorado()
│   ├── test-plot_distribution.R       # ✅ Phase 5: plot_methylation_distribution() tests (~159 lines)
│   ├── test-plot_genome_track.R       # ✅ Phase 5: plot_genome_track() tests (~145 lines)
│   ├── test-plot_metagene.R           # ✅ Phase 5: plot_metagene() tests (~120 lines)
│   ├── test-plot_volcano.R            # ✅ Phase 5: plot_volcano() tests (~103 lines)
│   ├── test-plot_heatmap.R            # ✅ Phase 5: plot_heatmap() tests (~129 lines)
│   ├── test-m_values.R                # ✅ mValues() tests (~24 tests: formula, NA propagation, alpha, mod_type)
│   ├── test-plot_pca.R                # ✅ Phase 5: plot_pca() tests, including return_data (~22 tests)
│   └── test-plot_coverage.R           # ✅ Phase 5: plot_coverage() tests (~117 lines)
├── man/                          # Roxygen2-generated docs (all current)
├── .github/workflows/
│   ├── r.yml                     # rcmdcheck on push/PR (R 3.6.3 + 4.1.1, macOS-latest)
│   └── render-rmarkdown.yaml     # Auto-renders .Rmd on push
├── DESCRIPTION                   # v0.5.0; 12 Imports, 11 Suggests; R >= 4.1.0
├── NAMESPACE                     # commaData class + all Phase 1–5 exports
├── NEWS.md                       # v0.5.0, v0.4.0, v0.3.0, v0.2.0, and v0.1.0 entries
├── README.md / README.Rmd        # Reflects v0.3.0 — needs update for v0.5.0
└── CLAUDE.md                     # ← THIS FILE (AI assistant guide)
```

**Note:** All root-level legacy files have been removed or moved:
- Deleted: `functions.R`, `testscript.R`, `WT_6mA_Mg.txt`, `WT_6mA_all_callers.txt`, `all_site_annotations.txt`, `all_site_annotations_60p.txt`
- Moved: `methylKitGATC.R` → `inst/scripts/methylKitGATC_historical.R`
- Deleted: `comma_pm.md` — all essential content migrated into this file

**Note:** There is **no** `tests/testthat/helper-fixtures.R` — fixtures are defined inline within each test file, or `comma_example_data` is used directly via `data(comma_example_data)`.

### Implemented in v0.2.0 (Phase 1 & 2)

- **`commaData` S4 class** — extends `SummarizedExperiment`; slots: `genomeInfo`, `annotation`, `motifSites`; full `validity()` and `show()` methods
- **`commaData()` constructor** — dispatches to parser by `caller` arg; merges multi-sample matrices using site key (`chrom:position:strand:mod_type`); applies `min_coverage` thresholding
- **Modkit parser** (`.parseModkit`) — reads 15-column modkit `pileup` BED; maps mod codes (`a`→6mA, `m`→5mC, `21839`→4mC); 0-based→1-based conversion
- **Megalodon parser** (`.parseMegalodon`) — per-read aggregation to per-site beta values; explicit `mod_type` required
- **Accessor S4 methods** — `methylation()`, `coverage()`, `sampleInfo()`, `siteInfo()`, `modTypes()`, `genome()`, `annotation()`, `motifSites()`, `[`, `subset()`
- **`loadAnnotation()`** — GFF3/BED → GRanges with standardized feature_type/name columns
- **`findMotifSites()`** — BSgenome or FASTA + motif regex → GRanges (both strands, IUPAC support)
- **Genome utilities** — `.validateGenomeInfo()`, `.circularIndex()`, `.makeSeqinfo()`
- **`comma_example_data`** — synthetic commaData: 300 sites (200×6mA, 100×5mC), 3 samples, chr_sim (100 kb), differential ground truth in `rowData$is_diff`
- **Example files** — `inst/extdata/example_modkit.bed` (20 sites: 10×6mA, 5×5mC) and `inst/extdata/example.gff3` (5 genes on chr_sim)

### Added in v0.3.0 (Phase 3)

- **`annotateSites()`** — vectorized annotation using `GenomicRanges::findOverlaps()`; three modes (overlap, proximity, metagene); replaces deleted `annotateMethylSites()`, `annotateTSS()`, `annotateTTS()`
- **`slidingWindow()`** — generalized genome-wide smoothing; stat = "median" | "mean"; genome size always from `genomeInfo`, never hardcoded; uses `zoo::rollapply()`; replaces deleted `methylRollingMedian()`, `methylRollingMean()`
- **`methylomeSummary()`** — per-sample distribution stats (mean/median/SD beta, n_sites, frac_methylated, coverage stats); returns tidy data.frame for ggplot2
- **`coverageDepth()`** — windowed sequencing depth across genome; optional log2 transform; replaces `calculateMethylSiteDepth()`
- **`varianceByDepth()`** — methylation variance stratified by coverage level; replaces `varByCoverage()`
- **`writeBED()`** — fully rewritten; accepts `commaData`, output path, sample name; writes BED9 format with itemRGB methylation scale; no hardcoded paths

### Added in v0.4.0 (Phase 4)

- **`diffMethyl()`** — main differential methylation function (modeled on DESeq2's `DESeq()`); accepts `commaData` + formula; returns enriched `commaData` with `dm_pvalue`, `dm_padj`, `dm_delta_beta`, `dm_mean_beta_<cond>` in `rowData`; supports `method = "beta_binomial"` (default) and `method = "methylkit"` (requires methylKit)
- **`results()`** — S4 method to extract diff methylation table as tidy `data.frame`
- **`filterResults()`** — S4 method to filter results by padj and delta_beta thresholds
- **`.parseDorado()`** — full Dorado BAM parser; reads MM/ML tags via `Rsamtools::scanBam()`, CIGAR-decodes read positions, aggregates to per-site beta values; handles 6mA, 5mC, and 4mC in one BAM
- **`beta_binomial.R`** — internal per-site quasibinomial GLM engine
- **`methylkit_wrapper.R`** — internal methylKit dispatch wrapper
- **`multiple_testing.R`** — internal BH/FDR correction utility

### Added in v0.5.0 (Phase 5) — Current version

- **`plot_methylation_distribution()`** — beta value density plot per sample, coloured by sample name, faceted by modification type; QC and distribution comparison
- **`plot_genome_track()`** — genome browser-style scatter plot of methylation beta values vs. genomic position; supports positional windowing (`start`/`end`), `mod_type` filtering, and optional feature annotation rectangles from `annotation(object)`; uses `patchwork` for combined tracks
- **`plot_metagene()`** — average methylation profile across a class of genomic features, normalized to fractional position [0 = TSS, 1 = TTS]; uses `annotateSites(type = "metagene")` internally
- **`plot_volcano()`** — volcano plot for differential methylation results; accepts output of `results()`; colors points as Hypermethylated / Hypomethylated / Not significant
- **`plot_heatmap()`** — `geom_tile` heatmap of top differentially methylated sites ranked by adjusted p-value; ggplot2 only (no ComplexHeatmap dependency required)
- **`plot_pca()`** — PCA of per-sample methylation profiles using `stats::prcomp()`; points colored and optionally shaped by any column in `sampleInfo(object)`
- **`plot_coverage()`** — histogram of sequencing depth per site, per sample; coverage QC
- **Vignettes** — `vignettes/getting-started.Rmd` (end-to-end workflow) and `vignettes/multiple-modification-types.Rmd` (joint 6mA + 5mC analysis)
- **`comma-package.R`** — package-level `?comma` documentation page (Bioconductor requirement)
- **Tests for all plot functions** — every `plot_*()` function has a dedicated test file

### Breaking changes in v0.3.0

The following functions were **removed** (use their replacements):

| Removed | Replacement |
|---|---|
| `annotateMethylSites()` | `annotateSites(type = "overlap")` |
| `annotateTSS()` | `annotateSites(type = "proximity")` |
| `annotateTTS()` | `annotateSites(type = "proximity")` |
| `methylRollingMedian()` | `slidingWindow(stat = "median")` |
| `methylRollingMean()` | `slidingWindow(stat = "mean")` |
| `calculateMethylSiteDepth()` | `coverageDepth()` |
| `varByCoverage()` | `varianceByDepth()` |

---

## 3. Target Architecture

### Central data object: `commaData` S4 class

Every analysis function accepts a `commaData` object. Modeled on DESeq2's `DESeqDataSet`.

```r
commaData
├── methylation    # sites × samples matrix of beta values (0-1)
├── coverage       # sites × samples matrix of read depth
├── rowData        # per-site: chrom, position, strand, motif, mod_type
├── colData        # per-sample: sample_name, condition, replicate, caller, file_path
├── genomeInfo     # chromosome names and sizes (named integer vector or NULL)
├── annotation     # GRanges of genomic features (from GFF3/BED)
├── motifSites     # GRanges of all motif instances in genome
└── metadata       # list: package version, creation date, user fields
```

Key: `rowData` includes `mod_type` (`"6mA"`, `"5mC"`, `"4mC"`) as a first-class column. A single object can hold multiple modification types.

### Constructor

```r
commaData(
  files,         # named character vector: sample_name → file_path
  colData,       # data frame: sample_name, condition, replicate (minimum)
  genome,        # BSgenome, FASTA path, or named integer vector of chr sizes
  annotation,    # GFF3 path or GRanges (optional)
  mod_type,      # "6mA", "5mC", "4mC", or NULL to auto-detect
  motif,         # regex motif string (e.g., "GATC") or NULL
  min_coverage,  # integer, default 5
  caller         # "dorado", "modkit", "megalodon"
)
```

---

## 4. Naming Conventions

Follow these strictly. Inconsistency here will make the package feel unprofessional.

| Category | Convention | Examples |
|---|---|---|
| S4 class | `camelCase`, lowercase first letter | `commaData` |
| Constructor | Same as class name | `commaData()` |
| Analysis functions | `verbNoun()` camelCase | `annotateSites()`, `diffMethyl()`, `findMotifSites()` |
| Plot functions | `plot_noun()` snake_case with `plot_` prefix | `plot_volcano()`, `plot_metagene()` |
| Internal functions | `.` prefix | `.parseBetaValues()`, `.circularIndex()` |
| Arguments | `snake_case` throughout | `mod_type`, `min_coverage`, `position_col` |
| Test files | `test-functionName.R` | `test-annotateSites.R` |

---

## 5. Development Phases

| Phase | Version | Key deliverable | Status |
|---|---|---|---|
| 1 — Data Infrastructure | 0.2.0 | `commaData` S4 class + modkit parser + constructor | ✅ Complete |
| 2 — Genome Generalization | 0.2.0 | Remove all hardcoded MG1655 assumptions | ✅ Complete |
| 3 — Refactor Functions | 0.3.0 | Vectorized annotation, sliding window, coverage analysis, cleanup | ✅ Complete |
| 4 — Differential Methylation | 0.4.0 | `diffMethyl()` with beta-binomial model | ✅ Complete |
| 5 — Visualization & Release | 0.5.0 | All `plot_*()` functions, vignettes, package docs | ✅ Complete |
| Bioconductor submission | 1.0.0 | `BiocCheck` passing, full docs, DOI | ⏳ Next |

**All five phases are complete.** Bioconductor submission prep is the current priority.

---

## 6. Key Implementation Rules

### Always

- Every exported function must accept a `commaData` object as its primary input
- Use `GenomicRanges::findOverlaps()` for any genomic interval overlap — never nested for-loops
- Return tidy dataframes (or updated `commaData`) suitable for direct use with ggplot2
- All `plot_*()` functions return a `ggplot` object (or `patchwork` composite), not a rendered image
- Use `mod_type` as a parameter or infer it from `commaData@rowData$mod_type`
- Treat genome size as a parameter from `commaData@genomeInfo`, never hardcode
- Document every exported function with full roxygen2: `@param`, `@return`, `@examples`
- Write tests for every exported function using `testthat`
- **Preserve all annotated features**: `annotateSites()` stores ALL overlapping/nearby features per site as `CharacterList`/`IntegerList`/`NumericList` columns in `rowData`, not just the first or closest. This is a deliberate design decision reflecting the highly overlapping nature of bacterial genome annotations (genes, promoters, TF binding sites). Do NOT revert to single-match (`!duplicated()` or `distanceToNearest()`) behavior. Intergenic/non-overlapping sites receive length-0 list elements; test with `lengths(col) == 0`.

### Never

- Never hardcode genome size, chromosome names, or organism-specific values
- Never use nested R for-loops over genomic positions (use `GenomicRanges` instead)
- Never hardcode file paths
- Never import `tidyverse` as a package dependency — import `dplyr`, `tidyr` individually (Bioconductor requirement)
- Never write stub documentation like "A dataframe." or "A string." for `@param`/`@return`
- Never add features outside the current phase's scope without discussing first

### Performance

Any function that touches genomic positions must be vectorized:
- Use `GenomicRanges::findOverlaps()` for overlap queries
- Use `zoo::rollapply()` for sliding windows
- Use matrix operations, not element-wise R loops
- If R-level vectorization is insufficient, consider Rcpp

---

## 7. Dependencies

### Hard imports (`Imports` in DESCRIPTION — all currently declared)

| Package | Purpose |
|---|---|
| `GenomicRanges` | Core genomic interval arithmetic; `findOverlaps()` for annotation |
| `GenomeInfoDb` | Chromosome/genome metadata management |
| `SummarizedExperiment` | Base class infrastructure for `commaData` |
| `IRanges` | Range operations (via GenomicRanges) |
| `S4Vectors` | DataFrame and other S4 infrastructure used by SummarizedExperiment |
| `BiocGenerics` | Bioconductor generic methods |
| `Rsamtools` | BAM file parsing for Dorado input (`.parseDorado()` MM/ML tag parser) |
| `zoo` | Rolling window operations in `slidingWindow()` |
| `ggplot2` | All visualization |
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |
| `methods` | S4 class system (base R, but must be declared) |

### Soft dependencies (`Suggests` in DESCRIPTION — all currently declared)

| Package | Purpose |
|---|---|
| `BSgenome` | Genome sequence access for `findMotifSites()` |
| `Biostrings` | Sequence pattern matching for motif search |
| `BiocStyle` | Vignette styling (Bioconductor standard) |
| `ComplexHeatmap` | Available as alternative heatmap backend (not currently used — `plot_heatmap()` uses ggplot2) |
| `ggrepel` | Available for volcano plot labels (not currently used — `plot_volcano()` uses ggplot2 directly) |
| `methylKit` | Alternative differential methylation backend for `diffMethyl(..., method = "methylkit")` |
| `patchwork` | Multi-panel plot assembly in `plot_genome_track()` and `plot_heatmap()` |
| `rtracklayer` | GFF3 import via `import()` |
| `testthat` | Testing framework (edition 3) |
| `knitr` | R markdown processing for vignettes |
| `rmarkdown` | Vignette rendering |

---

## 8. Input Format Reference

### modkit BED (primary target format)

modkit `pileup` output columns (15 total):
```
chrom, start, end, mod_code, score, strand, coverage, mod_frequency,
n_mod, n_canonical, n_other_mod, n_delete, n_fail, n_diff, n_no_call
```

`mod_code` values to handle:
- `a` = 6mA (N6-methyladenine)
- `m` = 5mC (5-methylcytosine)
- `21839` = 4mC (N4-methylcytosine)

`mod_frequency` is the beta value (0–1); `coverage` is total read depth.
Coordinates are **0-based** — the parser converts to 1-based by computing `position = start + 1`.

### Dorado BAM

MM/ML tags in BAM format. Full implementation in `parse_dorado.R`. Reads MM/ML modification tags via `Rsamtools::scanBam()`, maps read positions to reference coordinates via CIGAR decoding, and aggregates per-read calls into per-site beta values. Supports 6mA, 5mC, and 4mC in one BAM. Invoked by `commaData(..., caller = "dorado")`.

The recommended workflow remains running `modkit pileup` first, then loading with `caller = "modkit"`, as direct BAM parsing is slower.

### Megalodon (backward compatibility)

Legacy format from earlier analysis. See `inst/scripts/methylKitGATC_historical.R` for historical parsing context. `mod_type` must be provided explicitly (cannot be inferred from file).

---

## 9. Testing

### Framework

`testthat` (edition 3). Tests live in `tests/testthat/`.

### Test fixture

Use `comma_example_data` — a synthetic `commaData` object created in Phase 1 (script: `data-raw/create_example_data.R`, `set.seed(42)`):
- **300 sites**: 200 × 6mA, 100 × 5mC
- **3 samples**: ctrl_1, ctrl_2, treat_1
- **2 conditions**: control (n=2), treatment (n=1)
- **Genome**: chr_sim, 100 kb
- **Ground truth**: ~30 of 200 6mA sites are differentially methylated (control ~0.9, treatment ~0.25); marked in `rowData$is_diff`
- **Annotation**: 5 simulated genes (GRanges)

### Current test files

| File | Coverage |
|---|---|
| `test-commaData.R` | S4 class validity, constructor, bad inputs, show() — ~20 tests |
| `test-parsers.R` | Modkit column mapping, mod codes, coverage filter — ~15 tests |
| `test-accessors.R` | Matrix shape, value ranges, multi-mod-type, subsetting — ~20 tests |
| `test-genome_utils.R` | .validateGenomeInfo, .circularIndex, .makeSeqinfo |
| `test-load_annotation.R` | GFF3/BED parsing, feature_type filtering |
| `test-find_motif_sites.R` | Motif search, both strands, palindromic motifs |
| `test-parse_megalodon.R` | .parseMegalodon aggregation, mod_type requirement |
| `test-annotateSites.R` | overlap/proximity/metagene modes, edge cases — ~20 tests |
| `test-slidingWindow.R` | stat modes, circular wrap, genome-size inference — ~15 tests |
| `test-methylomeSummary.R` | per-sample stats, mod_type filtering, all-NA sample column — ~11 tests |
| `test-coverageAnalysis.R` | coverageDepth() windowing, varianceByDepth() bins — ~8 tests |
| `test-writeBED.R` | file creation, track header, 9-col BED, RGB color bands, NA exclusion, mod_type filtering, errors — ~20 tests |
| `test-diffMethyl.R` | diffMethyl() basic, statistical correctness, mod_type/min_coverage/p_adjust, errors, .applyMultipleTesting() contract, .betaBinomialTest() edge cases, ground-truth recovery on comma_example_data — ~30 tests |
| `test-results.R` | results() and filterResults(): output shape, filtering, thresholds, boundary conditions, errors — ~23 tests |
| `test-parse_dorado.R` | .cigarToRefPos() (H/S clips, N skips, mixed), .parseMmTag() (5mC, multi-mod, ML boundary values 127/255 and 128/255), .parseDorado() error handling — ~21 tests |
| `test-plot_distribution.R` | plot_methylation_distribution() returns ggplot, mod_type filtering, per_sample |
| `test-plot_genome_track.R` | plot_genome_track() returns ggplot/patchwork, windowing, annotation |
| `test-plot_metagene.R` | plot_metagene() returns ggplot, feature normalization |
| `test-plot_volcano.R` | plot_volcano() returns ggplot, thresholds, coloring |
| `test-plot_heatmap.R` | plot_heatmap() returns ggplot, top-N sites, sample annotation |
| `test-m_values.R` | mValues(): formula correctness, NA/zero-coverage propagation, alpha validation, mod_type filter — ~24 tests |
| `test-plot_pca.R` | plot_pca(): returns ggplot, color_by/shape_by, return_data data.frame + percentVar attr — ~22 tests |
| `test-plot_coverage.R` | plot_coverage() returns ggplot, per_sample mode |

### Required coverage

Every exported function needs tests for:
- Valid input → correct output
- Invalid input → informative error message (not a cryptic R error)
- Edge cases (empty data, NA handling, single sample, etc.)

### Running tests

```r
devtools::test()
# or for a specific file:
testthat::test_file("tests/testthat/test-annotateSites.R")
```

---

## 10. Documentation

### roxygen2

All documentation is generated via roxygen2. Run `devtools::document()` to rebuild `man/` and `NAMESPACE`.

**Every exported function needs:**
```r
#' Short title (one line)
#'
#' Longer description paragraph explaining what the function does
#' and why you would use it.
#'
#' @param object A \code{commaData} object.
#' @param mod_type Character string specifying the modification type
#'   (e.g., \code{"6mA"}, \code{"5mC"}). If \code{NULL}, uses all
#'   types present in \code{object}.
#' @return Description of what is returned and its structure.
#' @examples
#' # Example that runs without error
#' data(comma_example_data)
#' result <- functionName(comma_example_data)
#' @export
```

Stub documentation (`"A dataframe."`, `"A string."`) is not acceptable.

### Package-level documentation

`R/comma-package.R` provides the `?comma` package documentation page. It describes the five-step workflow: Load → QC → Annotate → Visualize → Differential methylation. This satisfies the Bioconductor requirement.

### Vignettes

Two vignettes in `vignettes/`:
- **`getting-started.Rmd`** (~214 lines) — end-to-end workflow using `comma_example_data`: construct → characterize → diff methylation → visualize
- **`multiple-modification-types.Rmd`** (~173 lines) — joint 6mA + 5mC analysis; demonstrates subsetting by `mod_type`, comparing patterns, visualizing both simultaneously

---

## 11. Bioconductor Requirements

The package targets Bioconductor submission at v1.0.0.

| Requirement | Status |
|---|---|
| Individual package imports (`dplyr`, not `tidyverse`) | ✅ Done |
| S4 classes with proper `validity()` methods | ✅ Done |
| `show()` methods for all S4 classes | ✅ Done |
| Package-level `?comma` documentation page | ✅ Done |
| At least two vignettes | ✅ Done |
| `NEWS.md` with version history | ✅ Done |
| `biocViews` declared | ✅ Done (Sequencing, Epigenetics, Coverage, DifferentialMethylation, GenomeAnnotation, DataImport, Visualization) |
| `R CMD check --as-cran` zero errors/warnings | ⏳ Verify |
| `BiocCheck::BiocCheck()` zero errors | ⏳ Run and fix |
| Bundled data < 5 MB total | ⏳ Verify |
| Zenodo DOI | ⏳ Register before submission |
| Version bumped to 1.0.0 | ⏳ Pending |

---

## 12. What Has Been Cleaned Up

| File | Action | Status |
|---|---|---|
| `functions.R` (root) | Deleted — superseded by `R/` implementations | ✅ Done |
| `testscript.R` (root) | Deleted — not needed | ✅ Done |
| `methylKitGATC.R` (root) | Moved to `inst/scripts/methylKitGATC_historical.R` | ✅ Done |
| `WT_6mA_Mg.txt` (root) | Deleted — legacy data file | ✅ Done |
| `WT_6mA_all_callers.txt` (root) | Deleted — legacy data file | ✅ Done |
| `all_site_annotations.txt` (root) | Deleted — legacy data file | ✅ Done |
| `all_site_annotations_60p.txt` (root) | Deleted — legacy data file | ✅ Done |
| `data/*.rda` (MG1655 files) | Removed — replaced by `comma_example_data` | ✅ Done |
| `writeBED.R` | Rewritten — generalized, no hardcoded paths | ✅ Done |
| `parse_dorado.R` stub | Replaced with full Dorado BAM parser (MM/ML tags) | ✅ Done |
| `comma_pm.md` | Deleted — all essential content migrated to CLAUDE.md | ✅ Done |

---

## 13. Git and CI/CD

### Branches

- Stable: `master` — do not push here directly; work through PRs
- Tagged: `0.2.0` — snapshot of Phase 1+2 complete state
- Development branches follow `claude/<description>-<id>` naming pattern for AI-initiated work

### Commit style

Use descriptive, imperative commit messages:
```
Add commaData S4 class definition and show() method
Fix circular genome arithmetic to use genomeInfo slot
Replace nested for-loops in annotateSites with findOverlaps
Implement Phase 5 (v0.5.0): visualization functions and vignettes
```

### CI/CD

GitHub Actions (`.github/workflows/r.yml`) runs `rcmdcheck` on push/PR against R 3.6.3 and 4.1.1 on macOS-latest. Keep the package passing `R CMD check` throughout development.

`.github/workflows/render-rmarkdown.yaml` auto-renders `.Rmd` files on push (for README.Rmd → README.md).

---

## 14. R Environment Setup in This Environment

### When you need R

**Not every task requires R to be running.** You can safely work without invoking R for:
- Writing or editing R source code (`.R` files, `CLAUDE.md`, `DESCRIPTION`, `NEWS.md`, etc.)
- Writing documentation, vignettes, or roxygen2 comments
- Reviewing logic, refactoring, or adding new functions

**You DO need a working R environment when:**
- Running the test suite (`devtools::test()`)
- Checking the package (`devtools::check()`)
- Rebuilding documentation (`devtools::document()`)
- Verifying that new code actually runs without errors

### What is already installed

This is an Ubuntu/Debian Linux environment. R 4.3.3 and all required package dependencies for `comma` are **pre-installed**:

- **R itself:** `/usr/bin/R` (r-base 4.3.3 via `apt`)
- **Bioconductor core:** `GenomicRanges`, `IRanges`, `SummarizedExperiment`, `S4Vectors`, `GenomeInfoDb`, `Rsamtools`, `Biostrings`, `BSgenome`, `rtracklayer`, `BiocGenerics` — installed via `apt` (`r-bioc-*`)
- **CRAN packages:** `zoo`, `ggplot2`, `dplyr`, `tidyr`, `devtools`, `testthat`, `knitr`, `rmarkdown`, `ggrepel`, `patchwork` — installed in `/usr/local/lib/R/site-library/`
- **`BiocManager`** is available for installing additional Bioconductor packages if needed

Verify R is available with:
```bash
R --version
```

### Installing missing packages (if needed)

`sudo` is available without a password in this environment. Prefer `apt` for system-level packages (faster, no compilation).

```bash
# Bioconductor packages
sudo apt install -y r-bioc-<pkgname>   # e.g., r-bioc-genomicranges

# CRAN packages
sudo apt install -y r-cran-<pkgname>   # e.g., r-cran-zoo
```

Fallback — install from R directly (writes to `/usr/local/lib/R/site-library/`):
```r
install.packages("pkgname")
BiocManager::install("pkgname")
```

### Running package checks and tests

```bash
# From repo root
Rscript -e "devtools::test()"
Rscript -e "devtools::check()"
Rscript -e "devtools::document()"
```

**Do not give up because R is not on PATH or packages appear missing.** R 4.3.3 is at `/usr/bin/R` and all `comma` dependencies are pre-installed. Install any missing package with `sudo apt install r-bioc-<name>` or `install.packages()`.

---

## 15. Out of Scope for v1.0

Do not implement these without explicit discussion:

- Multi-species comparative methylomics
- Integration with transcriptomics (RNA-seq correlation)
- Motif discovery
- Phage/plasmid methylation analysis
- Shiny interactive browser
- Python or command-line interface
- Genome browser track export beyond BED (bigWig, etc.)

---

## 16. Quick Reference

### Install for development

```r
devtools::install()   # or devtools::load_all()
```

### Run checks

```r
devtools::check()        # Full R CMD check
devtools::test()         # Tests only
devtools::document()     # Rebuild docs from roxygen2
BiocCheck::BiocCheck()   # Bioconductor-specific checks
```

### Currently exported API (v0.5.0)

```r
# S4 class (Phase 1)
commaData(files, colData, genome, annotation, mod_type, motif, min_coverage, caller)

# Accessors (Phase 1) — all accept a commaData object
methylation(object)      # → sites × samples beta matrix
coverage(object)         # → sites × samples integer matrix
sampleInfo(object)       # → per-sample DataFrame
siteInfo(object)         # → per-site DataFrame (chrom, position, strand, mod_type, ...)
modTypes(object)         # → character vector of modification types present
genome(object)           # → named integer vector of chromosome sizes
annotation(object)       # → GRanges of genomic features
motifSites(object)       # → GRanges of motif instances

# Subsetting
object[sites, samples]   # numeric/logical index
subset(object, ...)      # subset by mod_type, condition, chrom

# Utilities (Phase 1)
loadAnnotation(file, feature_types)   # GFF3/BED → GRanges
findMotifSites(genome, motif)         # genome + motif → GRanges

# Analysis functions (Phase 3)
annotateSites(object, features, type, ...)    # type = "overlap"|"proximity"|"metagene"
slidingWindow(object, window, stat, ...)      # stat = "median"|"mean"
methylomeSummary(object, mod_type)            # per-sample QC stats → tidy data.frame
coverageDepth(object, window, method, ...)    # windowed sequencing depth → tidy data.frame
varianceByDepth(object, coverage_bins)        # methylation variance by depth → tidy data.frame
writeBED(object, file, sample, ...)           # write BED9 output file
mValues(object, alpha, mod_type)              # beta values → M-values matrix (variance-stabilized)

# Differential methylation (Phase 4)
diffMethyl(object, formula, method, mod_type, min_coverage, p_adjust_method)
                                              # → commaData with dm_* results in rowData
results(object, mod_type)                     # → tidy data.frame of diff methylation results
filterResults(object, padj, delta_beta, ...)  # → filtered data.frame

# Visualization (Phase 5) — all return a ggplot (or patchwork) object
plot_methylation_distribution(object, mod_type, per_sample)
plot_genome_track(object, chromosome, start, end, mod_type)
plot_metagene(object, feature, mod_type, window)
plot_volcano(results_df, delta_beta_threshold, padj_threshold)
plot_heatmap(object, result_df, n_sites, annotation_cols)
plot_pca(object, mod_type, color_by, shape_by, return_data)
plot_coverage(object, per_sample)
```

---

## 17. Bioconductor Submission Checklist (Next Steps)

These are the remaining tasks before v1.0.0 submission:

1. **Update README** — README.md/README.Rmd still reflects v0.3.0; update to showcase Phase 5 visualization functions and the complete five-step workflow
2. **Run `BiocCheck::BiocCheck()`** — address all errors and warnings
3. **Verify `R CMD check --as-cran`** — must pass with zero errors and zero warnings
4. **Register DOI** — create a Zenodo release and register a DOI before submission
5. **Verify data size** — bundled data must be < 5 MB total (`data/` + `inst/extdata/`)
6. **Bump version to 1.0.0** — update DESCRIPTION and add NEWS.md entry
7. **Submit to Bioconductor** — follow instructions at https://contributions.bioconductor.org/

---

*Last updated: March 2026 (v0.5.0 — all five phases complete; 7 visualization functions, 2 vignettes, and package-level docs added in Phase 5; Bioconductor submission is the current priority; Section 14 documents R environment setup for agents)*
