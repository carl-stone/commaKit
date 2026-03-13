# CLAUDE.md — AI Assistant Guide for `comma`

> This file is the primary reference for AI coding agents working on this repository. Read it fully before making any changes. The authoritative project management document is `comma_pm.md` — read that too before implementing anything significant.

---

## 1. Project Overview

**Package name:** `comma` (was `CoMMA`; rename complete in package internals)
**Full name:** Comparative Methylomics for Microbial Analysis
**Author:** Carl Stone, Vanderbilt University (carl.j.stone@vanderbilt.edu)
**Current version:** 0.4.0
**License:** MIT
**Target:** Bioconductor submission at v1.0.0

`comma` is an R package for analyzing bacterial DNA methylation from Oxford Nanopore sequencing data. It characterizes genome-wide methylation patterns, annotates methylation sites relative to genomic features, and identifies differentially methylated sites between conditions. Phases 1–4 of the architectural refactor are complete — see `comma_pm.md` for the full design specification and roadmap.

---

## 2. Current Repository State

### What exists now (v0.4.0 — Phases 1, 2, 3 & 4 complete)

```
CoMMA/
├── R/
│   ├── commaData_class.R        # ✅ S4 class def (~175 lines), show(), validity()
│   ├── commaData_constructor.R  # ✅ commaData() constructor (~270 lines)
│   ├── accessors.R              # ✅ methylation(), coverage(), sampleInfo(), etc. (~150 lines)
│   ├── parse_modkit.R           # ✅ PRIMARY input format parser (~175 lines)
│   ├── parse_dorado.R           # ✅ Phase 4: full Dorado BAM parser with MM/ML tags (~220 lines)
│   ├── parse_megalodon.R        # ✅ backward compatibility parser
│   ├── load_annotation.R        # ✅ GFF3/BED → GRanges (~143 lines)
│   ├── find_motif_sites.R       # ✅ FASTA + motif regex → GRanges (~130 lines)
│   ├── genome_utils.R           # ✅ genome validation, circular arithmetic (~104 lines)
│   ├── comma_example_data.R     # ✅ roxygen docs for synthetic dataset
│   ├── annotateSites.R          # ✅ Phase 3: vectorized annotation using findOverlaps() (~252 lines)
│   ├── sliding_window.R         # ✅ Phase 3: generalized slidingWindow() (~168 lines)
│   ├── methylome_summary.R      # ✅ Phase 3: per-sample QC stats (~106 lines)
│   ├── coverage_analysis.R      # ✅ Phase 3: coverageDepth() + varianceByDepth() (~205 lines)
│   ├── writeBED.R               # ✅ Phase 3: rewritten, fully generalized (~191 lines)
│   ├── diffMethyl.R             # ✅ Phase 4: main differential methylation interface (~220 lines)
│   ├── beta_binomial.R          # ✅ Phase 4: quasibinomial GLM per-site engine (~165 lines)
│   ├── methylkit_wrapper.R      # ✅ Phase 4: methylKit alternative method wrapper (~165 lines)
│   ├── multiple_testing.R       # ✅ Phase 4: BH/FDR correction utility (~30 lines)
│   └── results_methods.R        # ✅ Phase 4: results() + filterResults() S4 methods (~140 lines)
├── data/
│   └── comma_example_data.rda   # ✅ synthetic commaData (300 sites, 3 samples, chr_sim 100kb)
├── data-raw/
│   └── create_example_data.R    # ✅ generates comma_example_data (set.seed(42))
├── inst/
│   ├── extdata/
│   │   ├── example_modkit.bed   # ✅ 20-site modkit pileup example (chr_sim, mixed 6mA+5mC)
│   │   └── example.gff3         # ✅ 5-gene GFF3 annotation (chr_sim)
│   └── scripts/
│       └── methylKitGATC_historical.R  # ✅ Moved from root in Phase 3 (historical reference)
├── tests/testthat/
│   ├── helper-fixtures.R              # ✅ shared test fixtures (minimal 10kb genome, 3 sites)
│   ├── test-commaData.R               # ✅ ~20 tests: S4 class, validity, constructor, show()
│   ├── test-parsers.R                 # ✅ ~15 tests: .parseModkit(), .parseDorado() (now full impl), coverage filter
│   ├── test-accessors.R               # ✅ ~20 tests: all accessor methods, subsetting
│   ├── test-genome_utils.R            # ✅ tests for .validateGenomeInfo(), .circularIndex(), .makeSeqinfo()
│   ├── test-load_annotation.R         # ✅ tests for loadAnnotation() GFF3/BED parsing
│   ├── test-find_motif_sites.R        # ✅ tests for findMotifSites()
│   ├── test-parse_megalodon.R         # ✅ tests for .parseMegalodon()
│   ├── test-annotateSites.R           # ✅ ~20 tests for annotateSites() (overlap/proximity/metagene)
│   ├── test-slidingWindow.R           # ✅ ~15 tests for slidingWindow()
│   ├── test-methylomeSummary.R        # ✅ ~10 tests for methylomeSummary()
│   ├── test-coverageAnalysis.R        # ✅ ~8 tests for coverageDepth() and varianceByDepth()
│   ├── test-diffMethyl.R              # ✅ Phase 4: ~25 tests for diffMethyl()
│   ├── test-results.R                 # ✅ Phase 4: ~15 tests for results() and filterResults()
│   └── test-parse_dorado.R            # ✅ Phase 4: ~10 tests for .parseDorado() helpers
├── man/                          # Roxygen2-generated docs (all current)
├── .github/workflows/
│   ├── r.yml                     # rcmdcheck on push/PR (R 3.6.3 + 4.1.1, macOS-latest)
│   └── render-rmarkdown.yaml     # auto-renders .Rmd on push
├── DESCRIPTION                   # v0.4.0; 12 Imports, 7 Suggests; R >= 4.1.0
├── NAMESPACE                     # commaData class + all Phase 1–4 exports
├── NEWS.md                       # v0.4.0, v0.3.0, v0.2.0, and v0.1.0 entries
├── README.md / README.Rmd        # ✅ Updated for v0.3.0 with Phase 1/2/3 examples
├── CLAUDE.md                     # ← THIS FILE (AI assistant guide)
└── comma_pm.md                   # ← PROJECT MANAGEMENT DOCUMENT (read this)
```

**Note:** All root-level legacy files have been removed or moved in Phase 3:
- Deleted: `functions.R`, `testscript.R`, `WT_6mA_Mg.txt`, `WT_6mA_all_callers.txt`, `all_site_annotations.txt`, `all_site_annotations_60p.txt`
- Moved: `methylKitGATC.R` → `inst/scripts/methylKitGATC_historical.R`

### Implemented in v0.2.0 (Phase 1 & 2)

- **`commaData` S4 class** — extends `SummarizedExperiment`; slots: `genomeInfo`, `annotation`, `motifSites`; full `validity()` and `show()` methods
- **`commaData()` constructor** — dispatches to parser by `caller` arg; merges multi-sample matrices using site key (`chrom:position:strand:mod_type`); applies `min_coverage` thresholding
- **Modkit parser** (`.parseModkit`) — reads 15-column modkit `pileup` BED; maps mod codes (`a`→6mA, `m`→5mC, `21839`→4mC); 0-based→1-based conversion
- **Megalodon parser** (`.parseMegalodon`) — per-read aggregation to per-site beta values; explicit `mod_type` required
- **Dorado parser** (`.parseDorado`) — intentional stub; fails with helpful error recommending `modkit pileup`
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

### Added in v0.4.0 (Phase 4)

- **`diffMethyl()`** — main differential methylation function (modeled on DESeq2's `DESeq()`); accepts `commaData` + formula; returns enriched `commaData` with `dm_pvalue`, `dm_padj`, `dm_delta_beta`, `dm_mean_beta_<cond>` in `rowData`; supports `method = "beta_binomial"` (default, no extra deps) and `method = "methylkit"` (requires methylKit)
- **`results()`** — S4 method to extract diff methylation table as tidy `data.frame`
- **`filterResults()`** — S4 method to filter results by padj and delta_beta thresholds
- **`.parseDorado()`** — full Dorado BAM parser replacing the stub; reads MM/ML tags via `Rsamtools::scanBam()`, CIGAR-decodes read positions, aggregates to per-site beta values; handles 6mA, 5mC, and 4mC in one BAM
- **`beta_binomial.R`** — internal per-site quasibinomial GLM engine
- **`methylkit_wrapper.R`** — internal methylKit dispatch wrapper
- **`multiple_testing.R`** — internal BH/FDR correction utility

### Remaining issues (to fix in Phase 5)

1. **No vignettes** — required for Bioconductor submission; planned for Phase 5
2. **No visualization functions** — all `plot_*()` functions are Phase 5

---

## 3. Target Architecture (Where We're Going)

The full design is in `comma_pm.md`. This is the condensed version for quick reference.

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

### Target file structure (✅ = exists, ⚠️ = partial, ⏳ = pending)

```
R/
├── commaData_class.R         # ✅ S4 class def, show(), validity()
├── commaData_constructor.R   # ✅ commaData() constructor
├── accessors.R               # ✅ methylation(), coverage(), sampleInfo(), etc.
├── parse_modkit.R            # ✅ PRIMARY input format
├── parse_dorado.R            # ✅ Phase 4: full Dorado BAM parser (MM/ML tags)
├── parse_megalodon.R         # ✅ Backward compatibility
├── load_annotation.R         # ✅ GFF3/BED → GRanges
├── find_motif_sites.R        # ✅ FASTA + motif regex → GRanges
├── genome_utils.R            # ✅ Circular arithmetic, seqinfo helpers
├── annotateSites.R           # ✅ Phase 3: vectorized annotation
├── sliding_window.R          # ✅ Phase 3: generalized smoothing
├── methylome_summary.R       # ✅ Phase 3: per-sample QC stats
├── coverage_analysis.R       # ✅ Phase 3: depth windowing, variance
├── writeBED.R                # ✅ Phase 3: generalized BED export
├── diffMethyl.R              # ✅ Phase 4: main differential methylation interface
├── beta_binomial.R           # ✅ Phase 4: quasibinomial GLM per-site engine
├── methylkit_wrapper.R       # ✅ Phase 4: methylKit as alternative method
├── multiple_testing.R        # ✅ Phase 4: BH / q-value correction
├── results_methods.R         # ✅ Phase 4: results(), filterResults()
├── plot_distribution.R       # ⏳ Phase 5
├── plot_genome_track.R       # ⏳ Phase 5
├── plot_metagene.R           # ⏳ Phase 5
├── plot_volcano.R            # ⏳ Phase 5
├── plot_heatmap.R            # ⏳ Phase 5
├── plot_pca.R                # ⏳ Phase 5
└── plot_coverage.R           # ⏳ Phase 5
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

See `comma_pm.md` Section 4 for full task lists. Work sequentially — each phase depends on the previous.

| Phase | Version | Key deliverable | Status |
|---|---|---|---|
| 1 — Data Infrastructure | 0.2.0 | `commaData` S4 class + modkit parser + constructor | ✅ Complete |
| 2 — Genome Generalization | 0.2.0 | Remove all hardcoded MG1655 assumptions | ✅ Complete |
| 3 — Refactor Functions | 0.3.0 | Vectorized annotation, sliding window, coverage analysis, cleanup | ✅ Complete |
| 4 — Differential Methylation | 0.4.0 | `diffMethyl()` with beta-binomial model | ✅ Complete |
| 5 — Visualization & Release | 0.5.0 | All `plot_*()` functions, real tests, vignettes | ⏳ Next |
| Bioconductor submission | 1.0.0 | `BiocCheck` passing, full docs | ⏳ Pending |

**Phases 1, 2, and 3 are complete.** Phase 4 is the current priority. Do not skip ahead to Phase 5.

---

## 6. Key Implementation Rules

### Always

- Every exported function must accept a `commaData` object as its primary input
- Use `GenomicRanges::findOverlaps()` for any genomic interval overlap — never nested for-loops
- Return tidy dataframes (or updated `commaData`) suitable for direct use with ggplot2
- All `plot_*()` functions return a `ggplot` object, not a rendered image
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
| `Rsamtools` | BAM file parsing for Dorado input (stub; full parser Phase 4) |
| `zoo` | Rolling window operations in `slidingWindow()` |
| `ggplot2` | All visualization (Phase 5) |
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |
| `methods` | S4 class system (base R, but must be declared) |

### Soft dependencies (`Suggests` in DESCRIPTION — all currently declared)

| Package | Purpose |
|---|---|
| `BSgenome` | Genome sequence access for `findMotifSites()` |
| `Biostrings` | Sequence pattern matching for motif search |
| `rtracklayer` | GFF3 import via `import()` |
| `testthat` | Testing framework (edition 3) |
| `knitr` | R markdown processing for vignettes |
| `rmarkdown` | Vignette rendering |

> **Note:** `ComplexHeatmap` and `ggrepel` are in the future plan but not yet declared in DESCRIPTION. Add them when Phase 5 visualization functions are implemented. `methylKit` is planned as a Suggests for the Phase 4 methylKit wrapper but is not yet declared — add it then.

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

MM/ML tags in BAM format. Currently a stub (`parse_dorado.R`) that errors with a helpful message recommending `modkit pileup`. Full implementation is Phase 4. Use `Rsamtools` for reading.

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

Also available: `tests/testthat/helper-fixtures.R` — minimal shared fixtures (10kb genome, 3 sites, 2 features) for fast unit tests that don't need the full example dataset.

### Current test files

| File | Coverage | Tests |
|---|---|---|
| `test-commaData.R` | S4 class validity, constructor, bad inputs, show() | ~20 |
| `test-parsers.R` | Modkit column mapping, mod codes, coverage filter, Dorado stub | ~15 |
| `test-accessors.R` | Matrix shape, value ranges, multi-mod-type, subsetting | ~20 |
| `test-genome_utils.R` | .validateGenomeInfo, .circularIndex, .makeSeqinfo | ~5 |
| `test-load_annotation.R` | GFF3/BED parsing, feature_type filtering | ~5 |
| `test-find_motif_sites.R` | Motif search, both strands, palindromic motifs | ~5 |
| `test-parse_megalodon.R` | .parseMegalodon aggregation, mod_type requirement | ~5 |
| `test-annotateSites.R` | overlap/proximity/metagene modes, edge cases | ~20 |
| `test-slidingWindow.R` | stat modes, circular wrap, genome-size inference | ~15 |
| `test-methylomeSummary.R` | per-sample stats, mod_type filtering | ~10 |
| `test-coverageAnalysis.R` | coverageDepth() windowing, varianceByDepth() bins | ~8 |

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

Write a `?comma` package documentation page explaining the overall workflow. This is a Bioconductor requirement — planned for Phase 5.

---

## 11. Bioconductor Requirements

The package targets Bioconductor submission at v1.0.0. Follow these requirements throughout development:

- Individual package imports only (`dplyr`, not `tidyverse`)
- `S4` classes with proper `validity()` methods
- `show()` methods for all S4 classes
- Package passes `R CMD check --as-cran` with zero errors and zero warnings
- `BiocCheck::BiocCheck()` passes with zero errors
- Bundled data < 5 MB total
- At least two vignettes
- `NEWS.md` with version history

---

## 12. What Has Been Cleaned Up

All Phase 3 cleanup tasks are complete:

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
| `parse_dorado.R` stub | Replaced with full Dorado BAM parser (MM/ML tags) in Phase 4 | ✅ Done |

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
Implement Phase 3 (v0.3.0): vectorized annotation, sliding window, coverage analysis
```

### CI/CD

GitHub Actions (`.github/workflows/r.yml`) runs `rcmdcheck` on push/PR against R 3.6.3 and 4.1.1 on macOS-latest. Keep the package passing `R CMD check` throughout development.

`.github/workflows/render-rmarkdown.yaml` auto-renders `.Rmd` files on push (for README.Rmd → README.md).

---

## 14. Out of Scope for v1.0

Do not implement these without explicit discussion:

- Multi-species comparative methylomics
- Integration with transcriptomics (RNA-seq correlation)
- Motif discovery
- Phage/plasmid methylation analysis
- Shiny interactive browser
- Python or command-line interface
- Genome browser track export beyond BED (bigWig, etc.)

---

## 15. Quick Reference

### Install for development

```r
# From repo root
devtools::install()
# or
devtools::load_all()
```

### Run checks

```r
devtools::check()        # Full R CMD check
devtools::test()         # Tests only
devtools::document()     # Rebuild docs from roxygen2
BiocCheck::BiocCheck()   # Bioconductor-specific checks (Phase 5)
```

### Currently exported API (v0.4.0)

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
annotateSites(object, features, type, ...)    # vectorized annotation; type = "overlap"|"proximity"|"metagene"
slidingWindow(object, window, stat, ...)      # genome-wide smoothing; stat = "median"|"mean"
methylomeSummary(object, mod_type)            # per-sample QC stats → tidy data.frame
coverageDepth(object, window, method, ...)    # windowed sequencing depth → tidy data.frame
varianceByDepth(object, coverage_bins)        # methylation variance by depth → tidy data.frame
writeBED(object, file, sample, ...)           # write BED9 output file

# Differential methylation (Phase 4)
diffMethyl(object, formula, method, mod_type, min_coverage, p_adjust_method)
                                              # → commaData with dm_* results in rowData
results(object, mod_type)                     # → tidy data.frame of diff methylation results
filterResults(object, padj, delta_beta, ...)  # → filtered data.frame
```

---

## 16. Phase 5 Implementation Guide (Next Phase)

Phase 4 is complete. This section describes Phase 5.

### Goal: Visualization, Tests, Vignettes, and Bioconductor Release (v0.5.0)

Phase 5 completes the user-facing experience and prepares for Bioconductor submission.

### Priority order for Phase 5

1. **Visualization functions** — all `plot_*()` functions return `ggplot` objects:
   - `plot_methylation_distribution(object, mod_type, per_sample)` — beta density/ECDF per sample
   - `plot_genome_track(object, chromosome, start, end, mod_type)` — genome browser style
   - `plot_volcano(results_df, delta_beta_threshold, padj_threshold)` — volcano plot
   - `plot_heatmap(object, result_df, n_sites)` — heatmap of top diff sites
   - `plot_metagene(object, feature, mod_type, window)` — metagene plot
   - `plot_pca(object, mod_type, color_by)` — PCA on methylation profiles
   - `plot_coverage(object, per_sample)` — coverage distribution QC

2. **Vignettes** — required for Bioconductor:
   - "Getting Started with comma" — end-to-end workflow using `comma_example_data`
   - "Working with Multiple Modification Types" — 6mA + 5mC joint analysis

3. **Package-level documentation** — `?comma` page explaining overall workflow

4. **Bioconductor submission prep** — `BiocCheck::BiocCheck()`, `R CMD check --as-cran`

---

*Last updated: March 2026 (v0.4.0 — Phases 1, 2, 3 & 4 complete; Phase 5 is next priority)*
*See `comma_pm.md` for complete design rationale and implementation details.*
