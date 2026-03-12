# CLAUDE.md — AI Assistant Guide for `comma`

> This file is the primary reference for AI coding agents working on this repository. Read it fully before making any changes. The authoritative project management document is `comma_pm.md` — read that too before implementing anything significant.

---

## 1. Project Overview

**Package name:** `comma` (was `CoMMA`; rename is pending)
**Full name:** Comparative Methylomics for Microbial Analysis
**Author:** Carl Stone, Vanderbilt University (carl.j.stone@vanderbilt.edu)
**Current version:** 0.2.0
**License:** MIT
**Target:** Bioconductor submission at v1.0.0

`comma` is an R package for analyzing bacterial DNA methylation from Oxford Nanopore sequencing data. It characterizes genome-wide methylation patterns, annotates methylation sites relative to genomic features, and identifies differentially methylated sites between conditions. Phases 1 and 2 of the architectural refactor are complete — see `comma_pm.md` for the full design specification and roadmap.

---

## 2. Current Repository State

### What exists now (v0.2.0 — Phases 1 & 2 complete)

```
CoMMA/
├── R/
│   ├── commaData_class.R        # ✅ NEW: S4 class def, show(), validity()
│   ├── commaData_constructor.R  # ✅ NEW: commaData() constructor (~150 lines)
│   ├── accessors.R              # ✅ NEW: methylation(), coverage(), sampleInfo(), etc.
│   ├── parse_modkit.R           # ✅ NEW: PRIMARY input format parser
│   ├── parse_dorado.R           # ⚠️ STUB: errors with helpful message; deferred
│   ├── parse_megalodon.R        # ✅ NEW: backward compatibility parser
│   ├── load_annotation.R        # ✅ NEW: GFF3/BED → GRanges
│   ├── find_motif_sites.R       # ✅ NEW: FASTA + motif regex → GRanges
│   ├── genome_utils.R           # ✅ NEW: genome validation, circular arithmetic
│   ├── comma_example_data.R     # ✅ NEW: roxygen docs for synthetic dataset
│   ├── annotateMethylSites.R    # ⚠️ LEGACY: nested loops (O(n×m)); pending Phase 3
│   ├── annotateTSS.R            # ⚠️ LEGACY: TSS annotation; pending Phase 3
│   ├── annotateTTS.R            # ⚠️ LEGACY: not exported; pending Phase 3
│   ├── methylRollingMedian.R    # ⚠️ LEGACY: hardcoded genome_size default; pending Phase 3
│   ├── methylRollingMean.R      # ⚠️ LEGACY: not exported; pending Phase 3
│   ├── calculateMethylSiteDepth.R # ⚠️ LEGACY: not exported; pending Phase 3
│   ├── varByCoverage.R          # ⚠️ LEGACY: not exported; pending Phase 3
│   └── writeBED.R               # ⚠️ BROKEN: hardcoded developer paths; pending Phase 3
├── data/                        # (empty — MG1655 .rda files removed; example data in data-raw/)
├── data-raw/
│   └── create_example_data.R    # ✅ NEW: generates comma_example_data (set.seed(42))
├── inst/
│   ├── extdata/
│   │   ├── example_modkit.bed   # ✅ NEW: 20-site modkit example (chr_sim)
│   │   └── example.gff3         # ✅ NEW: 5-gene GFF3 example
│   └── scripts/                 # (empty; methylKitGATC.R not yet moved here)
├── tests/testthat/
│   ├── test-annotateMethylSites.R  # ⚠️ Placeholder only (2*2=4) — update in Phase 3
│   ├── test-commaData.R            # ✅ NEW: ~5-6 S4 class tests
│   ├── test-parsers.R              # ✅ NEW: ~10 modkit parser tests
│   └── test-accessors.R            # ✅ NEW: ~8 accessor tests
├── man/                         # Roxygen2-generated docs (commaData, findMotifSites, etc.)
├── .github/workflows/
│   ├── r.yml                    # rcmdcheck on push/PR
│   └── render-rmarkdown.yaml
├── DESCRIPTION                  # v0.2.0; 11 Imports, 6 Suggests
├── NAMESPACE                    # commaData class + all accessors exported
├── NEWS.md                      # v0.2.0 and v0.1.0 entries
├── README.md / README.Rmd       # ⚠️ Still shows v0.1.0 legacy examples; needs update
├── CLAUDE.md                    # ← THIS FILE (AI assistant guide)
├── comma_pm.md                  # ← PROJECT MANAGEMENT DOCUMENT (read this)
├── functions.R                  # Root-level scratchpad — NOT part of package; delete in Phase 3
├── methylKitGATC.R              # 513-line historical script — NOT packaged; move in Phase 3
└── testscript.R                 # Root-level scratch — NOT part of package; delete in Phase 3
```

### Implemented in v0.2.0

- **`commaData` S4 class** — extends `SummarizedExperiment`; slots: `genomeInfo`, `annotation`, `motifSites`; full `validity()` and `show()` methods
- **`commaData()` constructor** — dispatches to parser by `caller` arg; merges multi-sample matrices; applies `min_coverage` thresholding
- **Modkit parser** (`.parseModkit`) — reads 15-column modkit `pileup` BED; maps mod codes (`a`→6mA, `m`→5mC, `21839`→4mC); 0-based→1-based conversion
- **Megalodon parser** (`.parseMegalodon`) — per-read aggregation to per-site beta values; explicit `mod_type` required
- **Dorado parser** (`.parseDorado`) — intentional stub; fails with helpful error recommending `modkit pileup`
- **Accessor S4 methods** — `methylation()`, `coverage()`, `sampleInfo()`, `siteInfo()`, `modTypes()`, `genome()`, `annotation()`, `motifSites()`, `[`, `subset()`
- **`loadAnnotation()`** — GFF3/BED → GRanges with standardized feature_type/name columns
- **`findMotifSites()`** — BSgenome or FASTA + motif regex → GRanges (both strands, IUPAC support)
- **Genome utilities** — `.validateGenomeInfo()`, `.circularIndex()`, `.makeSeqinfo()`
- **`comma_example_data`** — synthetic commaData: 300 sites (200×6mA, 100×5mC), 3 samples, chr_sim (100 kb), differential ground truth
- **Example files** — `inst/extdata/example_modkit.bed` (20 sites) and `inst/extdata/example.gff3` (5 genes)
- **Tests** — ~25–30 real tests across `test-commaData.R`, `test-parsers.R`, `test-accessors.R`

### Remaining legacy issues (to fix in Phase 3+)

1. **O(n×m) annotation** — `annotateMethylSites()` and `annotateTSS()` still use nested R for-loops; will be replaced by `annotateSites()` using `GenomicRanges::findOverlaps()`
2. **Hardcoded genome size** — `methylRollingMedian()` default `genome_size = 4641652` (MG1655); will be replaced with `commaData@genomeInfo`
3. **Hardcoded paths** — `writeBED.R` references developer's local filesystem paths; will be replaced with generalized version
4. **Root-level clutter** — `functions.R`, `testscript.R`, `methylKitGATC.R` remain at root; to be deleted/moved in Phase 3
5. **Placeholder test** — `test-annotateMethylSites.R` still contains only `expect_equal(2 * 2, 4)`; update when `annotateSites()` is implemented
6. **README outdated** — still shows v0.1.0 legacy function examples; needs update for Phase 1 API

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

### Target file structure (✅ = exists, ⏳ = pending)

```
R/
├── commaData_class.R         # ✅ S4 class def, show(), validity()
├── commaData_constructor.R   # ✅ commaData() constructor
├── accessors.R               # ✅ methylation(), coverage(), sampleInfo(), etc.
├── parse_modkit.R            # ✅ PRIMARY input format
├── parse_dorado.R            # ⚠️ Stub only — deferred
├── parse_megalodon.R         # ✅ Backward compatibility
├── load_annotation.R         # ✅ GFF3/BED → GRanges
├── find_motif_sites.R        # ✅ FASTA + motif regex → GRanges
├── genome_utils.R            # ✅ Circular arithmetic, seqinfo helpers
├── annotateSites.R           # ⏳ Phase 3: replaces annotateMethylSites/TSS/TTS
├── sliding_window.R          # ⏳ Phase 3: generalized methylRollingMedian/Mean
├── methylome_summary.R       # ⏳ Phase 3: distribution stats, per-sample QC
├── coverage_analysis.R       # ⏳ Phase 3: depth windowing, variance by coverage
├── diffMethyl.R              # ⏳ Phase 4: main differential methylation interface
├── beta_binomial.R           # ⏳ Phase 4: beta-binomial model
├── methylkit_wrapper.R       # ⏳ Phase 4: methylKit as alternative method
├── multiple_testing.R        # ⏳ Phase 4: BH / q-value correction
├── results_methods.R         # ⏳ Phase 4: results(), filterResults()
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
| 3 — Refactor Functions | 0.3.0 | Vectorized annotation with `GenomicRanges::findOverlaps()` | ⏳ Next |
| 4 — Differential Methylation | 0.4.0 | `diffMethyl()` with beta-binomial model | ⏳ Pending |
| 5 — Visualization & Release | 0.5.0 | All `plot_*()` functions, real tests, vignettes | ⏳ Pending |
| Bioconductor submission | 1.0.0 | `BiocCheck` passing, full docs | ⏳ Pending |

**Phases 1 and 2 are complete.** Phase 3 is the current priority. Do not skip ahead to Phase 4 or 5.

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

### Hard imports (`Imports` in DESCRIPTION — all currently declared)

| Package | Purpose |
|---|---|
| `GenomicRanges` | Core genomic interval arithmetic; replaces for-loop annotation |
| `SummarizedExperiment` | Base class infrastructure for `commaData` |
| `IRanges` | Range operations (via GenomicRanges) |
| `S4Vectors` | DataFrame and other S4 infrastructure used by SummarizedExperiment |
| `Rsamtools` | BAM file parsing for Dorado input |
| `methylKit` | Differential methylation (alternative method) |
| `zoo` | Rolling window operations |
| `ggplot2` | All visualization |
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |
| `BiocGenerics` | Bioconductor generic methods |
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

> **Note:** `ComplexHeatmap` and `ggrepel` are in the future plan but not yet declared in DESCRIPTION. Add them when Phase 5 visualization functions are implemented.

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

Use `comma_example_data` — a synthetic `commaData` object created in Phase 1 (script: `data-raw/create_example_data.R`, `set.seed(42)`):
- **300 sites**: 200 × 6mA, 100 × 5mC
- **3 samples**: ctrl_1, ctrl_2, treat_1
- **2 conditions**: control (n=2), treatment (n=1)
- **Genome**: chr_sim, 100 kb
- **Ground truth**: ~30 of 200 6mA sites are differentially methylated (control ~0.9, treatment ~0.25)
- **Annotation**: 5 simulated genes (GRanges)

All tests should use this fixture. Do not create separate per-test data objects.

### Current test files

| File | Coverage | Tests |
|---|---|---|
| `test-commaData.R` | S4 class validity, bad inputs | ~5-6 |
| `test-parsers.R` | Modkit column mapping, mod codes, coverage filter | ~10 |
| `test-accessors.R` | Matrix shape, value ranges, multi-mod-type | ~8 |
| `test-annotateMethylSites.R` | **PLACEHOLDER ONLY** (`2*2=4`) | 1 |

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

| File | Action | Status |
|---|---|---|
| `functions.R` (root) | Delete — superseded by `R/` implementations | ⏳ Pending |
| `testscript.R` (root) | Delete — not needed | ⏳ Pending |
| `methylKitGATC.R` (root) | Move to `inst/scripts/methylKitGATC_historical.R` with a header comment | ⏳ Pending |
| `data/*.rda` (all 9 MG1655 files) | ✅ Already removed — replaced by `comma_example_data` | ✅ Done |
| `writeBED.R` | Replace with generalized, path-safe version or defer to post-1.0 | ⏳ Pending |
| `parse_dorado.R` stub | Implement full Dorado BAM parser (Phase 3+) | ⏳ Pending |

---

## 13. Git and CI/CD

### Branches

- Development: `claude/add-claude-documentation-phrAo` (current working branch)
- Stable: `master` — do not push here directly; work through PRs
- Tagged: `0.2.0` — snapshot of the Phase 1+2 complete state

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

### Currently exported API (v0.2.0)

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

# Legacy functions — still exported but pending Phase 3 replacement
annotateMethylSites(methyl_df, meta_df, location)  # → replace with annotateSites()
annotateTSS(methyl_df, tss_df, window, location)   # → replace with annotateSites()
methylRollingMedian(methyl_df, window, mode)        # → generalize in Phase 3
```

---

---

## 16. Phase 3 Implementation Guide (Next Phase)

This section provides a focused roadmap for the next developer to pick up.

### Priority order for Phase 3 (v0.3.0)

1. **`annotateSites()`** — the highest-value change
   - Replaces `annotateMethylSites()` and `annotateTSS()` / `annotateTTS()`
   - Accepts a `commaData` object + GRanges annotation
   - Uses `GenomicRanges::findOverlaps()` — no for-loops
   - Returns updated `commaData` with annotation columns added to `rowData`, or a tidy `data.frame`
   - Test file: `tests/testthat/test-annotateSites.R` (create new)

2. **`sliding_window.R`** — generalize `methylRollingMedian()` / `methylRollingMean()`
   - New function signature: `slidingWindow(object, window, stat = c("median", "mean"))`
   - Accept genome size from `object@genomeInfo`, never hardcode
   - Use `zoo::rollapply()` for the window computation
   - Remove the old functions from NAMESPACE after replacement

3. **`methylome_summary.R`** — new function(s)
   - Per-sample distribution stats (mean beta, fraction methylated, etc.)
   - Accepts `commaData`, returns tidy `data.frame` for ggplot2 use

4. **`coverage_analysis.R`** — generalize `calculateMethylSiteDepth()` / `varByCoverage()`
   - Accept `commaData` object; remove hardcoded column name assumptions

5. **Root-level cleanup**
   - Delete `functions.R`, `testscript.R`
   - Move `methylKitGATC.R` → `inst/scripts/methylKitGATC_historical.R`

6. **README update** — add Phase 1 example workflow (commaData construction, accessor usage)

### Key design constraint for Phase 3 functions

All new functions must:
- Accept `commaData` as primary input
- Return either an updated `commaData` or a tidy `data.frame`
- Use `GenomicRanges::findOverlaps()` for any genomic interval query
- Not hardcode any organism-specific values

---

*Last updated: March 2026 (v0.2.0 — Phases 1 & 2 complete)*
*See `comma_pm.md` for complete design rationale and implementation details.*
