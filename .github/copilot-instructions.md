# Copilot Instructions for commaKit

## What This Repository Is

**commaKit** (Comparative Microbial Methylomics Analysis Kit) is an R/Bioconductor package for bacterial DNA methylation analysis from Oxford Nanopore sequencing data. It supports multiple modification types (6mA, 5mC, 4mC), accepts input from modkit/Dorado/Megalodon callers, and provides differential methylation testing, annotation, enrichment analysis, and visualization.

- **Language:** R (>= 4.3.0)
- **CI:** R 4.5 on Ubuntu (pinned — S4Vectors C API breaks in R 4.6.0)
- **Package system:** Bioconductor conventions (not CRAN)
- **R package namespace:** `comma` (rename to `commaKit` planned)
- **Version:** 0.2.0
- **License:** MIT

## Build and Validation Commands

Always run these from the package root directory.

**Regenerate documentation** (always do this after editing roxygen comments in `R/*.R`):
```bash
Rscript -e "pkgload::load_all(); roxygen2::roxygenise()"
```

**Run tests** (two options — `devtools` may not be available in all environments):
```bash
# Preferred
Rscript -e "devtools::test()"

# Fallback if devtools is unavailable
Rscript -e "pkgload::load_all(); testthat::test_dir('tests/testthat')"
```

**Full R CMD check** (CI runs this — must pass with 0 errors, 0 warnings, 0 notes):
```bash
Rscript -e "devtools::check()"
```

**R CMD check without vignettes** (if pandoc unavailable):
```bash
Rscript -e "devtools::check(build_args = c('--no-build-vignettes'))"
```

**Install locally:**
```bash
Rscript -e "devtools::install()"
```

**Important:** `renv` is active. If `devtools` is not installed, use the `pkgload`/`testthat` fallback above. Do not install packages outside `renv`.

## CI Pipeline

Two GitHub Actions workflows run on every push and PR:

1. **R-CMD-check.yaml** — `R CMD check --no-manual --as-cran` on R 4.5 / Ubuntu. This runs with `--run-donttest`, so `\donttest{}` examples ARE executed. Use `\dontrun{}` for examples that need external files.
2. **pkgdown.yaml** — builds the pkgdown site. Deploys to gh-pages on main pushes only.

A PR must have both checks passing before merge.

## Project Layout

```
R/                          # Source code (30 .R files)
  commaData_class.R         # S4 class definition (extends RangedSummarizedExperiment)
  commaData_constructor.R   # Constructor (commaData())
  accessors.R               # Accessor methods (methylation, coverage, etc.)
  diffMethyl.R              # Differential methylation (3 backends)
  annotateSites.R           # Site-to-feature annotation
  enrichment.R              # GO/KEGG ORA + GSEA
  results_methods.R         # results() and filterResults()
  plot_*.R                  # 8 visualization functions
  parse_modkit.R            # modkit BED parser (primary)
  parse_dorado.R            # Dorado BAM parser
  parse_megalodon.R         # Megalodon parser (legacy)
man/                        # Generated Rd docs (do not edit directly)
tests/testthat/             # testthat tests (test-<function>.R naming)
tests/testthat.R            # test runner
data/                       # comma_example_data.rda (bundled dataset)
data-raw/                   # Script to regenerate example data
vignettes/                  # getting-started.Rmd, multiple-modification-types.Rmd
dev/                        # Project management (not in package — .Rbuildignore)
  knowledge/                # Durable findings (known-issues.md, test-quality.md, etc.)
  archive/                  # Historical documents (do not update)
DESCRIPTION                 # Package metadata
NAMESPACE                   # Generated (do not edit directly)
.Rbuildignore               # Excludes dev/, .github/, renv/, etc.
renv/                       # Locked dependency management
AGENTS.md                   # AI coding tool context (see also CLAUDE.md)
```

## Critical Conventions

**Naming:**
- S4 class: camelCase lowercase-first (`commaData`)
- Analysis functions: verbNoun camelCase (`annotateSites`, `diffMethyl`)
- Plot functions: `plot_noun()` snake_case (`plot_volcano`, `plot_metagene`)
- Internal functions: `.` prefix (`.parseBetaValues`)
- Arguments: snake_case (`mod_type`, `min_coverage`)
- Test files: `test-functionName.R`

**Hard rules:**
- Every exported function takes `commaData` as primary input
- Use `GenomicRanges::findOverlaps()` for interval overlap — never nested for-loops
- All `plot_*()` return ggplot/patchwork objects (except `plot_heatmap()` -> ComplexHeatmap)
- Genome size from `seqlengths(object)` (Seqinfo), never hardcoded
- Document every exported function with roxygen2: `@param`, `@return`, `@examples`
- Import `dplyr` and `tidyr` individually — never import `tidyverse` (Bioconductor requirement)
- Do not use `purrr::map_dfr()` (superseded) — use `map(...) |> list_rbind()`
- Do not use `purrr::map_dbl()` (superseded) — use `vapply()`

**Known gotchas:**
- `\donttest{}` examples run in CI — use `\dontrun{}` for examples needing external files
- `diag(x)` with scalar x creates an x*x identity matrix — use `diag(x, nrow=1)`
- `S4Vectors::rename()` masks `dplyr::rename()` — always use `dplyr::rename()` explicitly
- `matrixStats::count()` masks `dplyr::count()` — always use `dplyr::count()` explicitly
- `purrr::map()` / `mclust::map()` collision — use `lapply()` + `purrr::list_rbind()`
- methylKit crashes on zero-variance sites — comma wraps this, assigns p=1
- `org.EcK12.eg.db` requires `::` syntax in examples
- Non-ASCII characters (e.g., x) cause R CMD check notes

**Core design decisions:**
- `diffMethyl()` loops by `mod_context`, not `mod_type` — prevents spurious pooling
- Effect sizes always on beta scale (0-1), not M-value scale
- Multiple testing correction is genome-wide across all mod_contexts
- `annotateSites()` uses list-columns (CharacterList/IntegerList/NumericList) — do not revert to single-match
- `commaData` extends `RangedSummarizedExperiment` — genomic positions in `rowRanges()` (GRanges)
- `mod_context` derived on demand from `mod_type` + `motif`, not stored as a column
- Alignment by genomic position (findOverlaps), not string keys or row names

## Test Data

`comma_example_data` is a synthetic dataset: 588 sites (393 x 6mA GATC, 195 x 5mC CCWGG), 6 samples (3 control + 3 treatment), genome chr_sim (100 kb). 30 of the 393 6mA sites are ground-truth differentially methylated. Created by `data-raw/create_example_data.R` with `set.seed(1312)`.

## Before Submitting Changes

1. Run `roxygen2::roxygenise()` after any roxygen comment changes
2. Run the full test suite and verify 0 failures
3. Run `R CMD check` if possible
4. Check that `\donttest{}` examples don't reference external files
5. Verify no non-ASCII characters in R source or docs
6. Ensure `NAMESPACE` and `man/` are regenerated (they are auto-generated — never edit directly)
