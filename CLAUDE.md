# CLAUDE.md — AI Assistant Guide for `comma`

> This file is the primary reference for AI coding agents working on this repository. Read it fully before making any changes. The authoritative project management document is `comma_pm.md` — read that too before implementing anything significant.

---

## 1. Project Overview

**Package name:** `comma` (was `CoMMA`; rename is pending)
**Full name:** Comparative Methylomics for Microbial Analysis
**Author:** Carl Stone, Vanderbilt University (carl.j.stone@vanderbilt.edu)
**Current version:** 0.1.0
**License:** MIT
**Target:** Bioconductor submission at v1.0.0

`comma` is an R package for analyzing bacterial DNA methylation from Oxford Nanopore sequencing data. It characterizes genome-wide methylation patterns, annotates methylation sites relative to genomic features, and identifies differentially methylated sites between conditions. The package is undergoing a major architectural refactor — see `comma_pm.md` for the full design specification.

---

## 2. Current Repository State

### What exists now (v0.1.0)

```
CoMMA/
├── R/
│   ├── annotateMethylSites.R    # Exported: overlap annotation (SLOW — nested loops)
│   ├── annotateTSS.R            # Exported: TSS proximity annotation
│   ├── methylRollingMedian.R    # Exported: sliding window median
│   ├── annotateTTS.R            # Not exported: TTS annotation
│   ├── methylRollingMean.R      # Not exported: sliding window mean
│   ├── calculateMethylSiteDepth.R # Not exported: windowed depth
│   ├── varByCoverage.R          # Not exported: variance by coverage
│   └── writeBED.R               # BROKEN: hardcoded developer paths
├── data/                        # 9 MG1655-specific .rda files (~9.5 MB)
├── tests/testthat/
│   └── test-annotateMethylSites.R  # Placeholder only (2*2=4)
├── man/                         # Roxygen2-generated docs (3 functions)
├── .github/workflows/
│   ├── r.yml                    # rcmdcheck on push/PR
│   └── render-rmarkdown.yaml
├── DESCRIPTION
├── NAMESPACE
├── README.md / README.Rmd
├── comma_pm.md                  # ← PROJECT MANAGEMENT DOCUMENT (read this)
├── functions.R                  # Root-level scratchpad — NOT part of package
├── methylKitGATC.R              # 513-line historical analysis script — NOT packaged
└── testscript.R                 # Root-level scratch — NOT part of package
```

### Known issues to never replicate

1. **No data object** — functions accept raw dataframes with assumed column names (`Position`, `beta`, `Strand`, `Left`, `Right`, `Type`, `Site`) and no validation
2. **Hardcoded genome size** — `4641652` (*E. coli* K-12 MG1655) appears in sliding window and circular arithmetic code
3. **Hardcoded paths** — `writeBED.R` references `/Users/carlstone/Library/CloudStorage/...` and is non-functional for any other user
4. **O(n×m) performance** — `annotateMethylSites()` and `annotateTSS()` use nested R for-loops over every genomic position; this is unacceptably slow
5. **Single organism, single modification** — all bundled data is MG1655-specific; only 6mA (GATC) is supported anywhere
6. **Megalodon-only input** — Megalodon is deprecated; current ONT standard is Dorado + modkit
7. **Zero functional tests** — the single test file contains only `expect_equal(2 * 2, 4)`

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
├── genomeInfo     # chromosome names and sizes
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

### Target file structure

```
R/
├── commaData_class.R         # S4 class def, show(), validity()
├── commaData_constructor.R   # commaData() constructor
├── accessors.R               # methylation(), coverage(), sampleInfo(), etc.
├── parse_modkit.R            # PRIMARY input format
├── parse_dorado.R            # BAM with MM/ML tags
├── parse_megalodon.R         # Backward compatibility
├── load_annotation.R         # GFF3/BED → GRanges
├── find_motif_sites.R        # FASTA + motif regex → GRanges
├── genome_utils.R            # Circular arithmetic, seqinfo helpers
├── annotateSites.R           # Replaces annotateMethylSites/TSS/TTS
├── sliding_window.R          # Generalized methylRollingMedian/Mean
├── methylome_summary.R       # Distribution stats, per-sample QC
├── coverage_analysis.R       # Depth windowing, variance by coverage
├── diffMethyl.R              # Main differential methylation interface
├── beta_binomial.R           # Beta-binomial model
├── methylkit_wrapper.R       # methylKit as alternative method
├── multiple_testing.R        # BH / q-value correction
├── results_methods.R         # results(), filterResults()
├── plot_distribution.R
├── plot_genome_track.R
├── plot_metagene.R
├── plot_volcano.R
├── plot_heatmap.R
├── plot_pca.R
└── plot_coverage.R
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

| Phase | Version | Key deliverable |
|---|---|---|
| 1 — Data Infrastructure | 0.2.0 | `commaData` S4 class + modkit parser + constructor |
| 2 — Genome Generalization | 0.2.0 | Remove all hardcoded MG1655 assumptions |
| 3 — Refactor Functions | 0.3.0 | Vectorized annotation with `GenomicRanges::findOverlaps()` |
| 4 — Differential Methylation | 0.4.0 | `diffMethyl()` with beta-binomial model |
| 5 — Visualization & Release | 0.5.0 | All `plot_*()` functions, real tests, vignettes |
| Bioconductor submission | 1.0.0 | `BiocCheck` passing, full docs |

**Phase 1 is the prerequisite for everything.** Do not implement analysis functions until `commaData` is stable.

---

## 6. Key Implementation Rules

### Always

- Every exported function must accept a `commaData` object as its primary input (after Phase 1)
- Use `GenomicRanges::findOverlaps()` for any genomic interval overlap — never nested for-loops
- Return tidy dataframes (or updated `commaData`) suitable for direct use with ggplot2
- All `plot_*()` functions return a `ggplot` object, not a rendered image
- Use `mod_type` as a parameter or infer it from `commaData@rowData$mod_type`
- Treat genome size as a parameter from `commaData@genomeInfo`, never hardcode
- Document every exported function with full roxygen2: `@param`, `@return`, `@examples`
- Write tests for every exported function using `testthat`

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

### Hard imports (add to `Imports` in DESCRIPTION)

| Package | Purpose |
|---|---|
| `GenomicRanges` | Core genomic interval arithmetic; replaces for-loop annotation |
| `SummarizedExperiment` | Base class infrastructure for `commaData` |
| `IRanges` | Range operations (via GenomicRanges) |
| `Rsamtools` | BAM file parsing for Dorado input |
| `methylKit` | Differential methylation (alternative method) |
| `zoo` | Rolling window operations |
| `ggplot2` | All visualization |
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |
| `BiocGenerics` | Bioconductor generic methods |

### Soft dependencies (add to `Suggests`)

| Package | Purpose |
|---|---|
| `BSgenome` | Genome sequence access for `findMotifSites()` |
| `rtracklayer` | GFF3 import via `import()` |
| `ComplexHeatmap` | Heatmap visualization |
| `ggrepel` | Labeled points in volcano plots |
| `testthat` | Testing |
| `knitr` / `rmarkdown` | Vignettes |

---

## 8. Input Format Reference

### modkit BED (primary target format)

modkit `pileup` output columns:
```
chrom, start, end, mod_code, score, strand, coverage, mod_frequency,
n_mod, n_canonical, n_other_mod, n_delete, n_fail, n_diff, n_no_call
```

`mod_code` values to handle:
- `a` = 6mA (N6-methyladenine)
- `m` = 5mC (5-methylcytosine)
- `21839` = 4mC (N4-methylcytosine)

`mod_frequency` is the beta value (0–1); `coverage` is total read depth.

### Dorado BAM

MM/ML tags in BAM format. Lower priority than modkit (parse after modkit parser is stable). Use `Rsamtools` for reading.

### Megalodon (backward compatibility)

Legacy format from earlier analysis. See `methylKitGATC.R` for the parsing logic that was used historically.

---

## 9. Testing

### Framework

`testthat` (edition 3). Tests live in `tests/testthat/`.

### Test fixture

Build a small synthetic `commaData` object from `comma_example_data` (to be created in Phase 2 — a simulated 100kb genome, 3 samples, 2 conditions, 6mA + 5mC). All tests use this fixture.

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

Write a `?comma` package documentation page explaining the overall workflow.

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

## 12. What to Discard / Move

When implementing Phase 3:

| File | Action |
|---|---|
| `functions.R` (root) | Delete — superseded by `R/` implementations |
| `testscript.R` (root) | Delete — not needed |
| `methylKitGATC.R` (root) | Move to `inst/scripts/methylKitGATC_historical.R` with a header comment explaining it's for reference only |
| `data/*.rda` (all 9 MG1655 files) | Remove — replace with `comma_example_data` synthetic dataset |
| `writeBED.R` | Replace with generalized version or defer to post-1.0 |

---

## 13. Git and CI/CD

### Branches

- Development: `claude/create-claude-documentation-F8lyI` (current working branch)
- Main: `main` — do not push here directly

### Commit style

Use descriptive, imperative commit messages:
```
Add commaData S4 class definition and show() method
Fix circular genome arithmetic to use genomeInfo slot
Replace nested for-loops in annotateSites with findOverlaps
```

### CI/CD

GitHub Actions (`.github/workflows/r.yml`) runs `rcmdcheck` on push/PR against R 3.6.3 and 4.1.1 on macOS. Keep the package passing `R CMD check` throughout development.

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

### Current exported functions

```r
annotateMethylSites(methyl_df, meta_df, location)  # → will be replaced by annotateSites()
annotateTSS(methyl_df, tss_df, window, location)   # → will be replaced by annotateSites()
methylRollingMedian(methyl_df, window, mode)        # → will be generalized in Phase 3
```

---

*Last updated: March 2026*
*See `comma_pm.md` for complete design rationale and implementation details.*
