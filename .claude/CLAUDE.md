# commaKit / comma — AI Assistant Guide

**Public name:** commaKit (Comparative Microbial Methylomics Analysis Kit)
**R package namespace:** `commaKit`
**GitHub repo:** `carl-stone/commaKit`
**Author:** Carl Stone, Vanderbilt University (carl.j.stone@vanderbilt.edu)
**Version:** 0.2.0 | **License:** MIT
**Target:** Bioconductor submission is on hold; `0.99.0` is reserved for the eventual submission cycle.

This is an R package (Bioconductor ecosystem). Use R/Bioconductor idioms: S4 classes, tibbles over data.frames, check class compatibility before implementing.

## Scientific Scope

`commaKit` is **modification-type agnostic** — originally conceived for 6mA (GATC/Dam), but 4mC and 5mC have equally important roles. Dorado detects all three simultaneously. Every data structure, function signature, and analysis module must treat methylation type as a first-class parameter, never an assumption. The package must work equally well for 6mA, 4mC, 5mC, and any future modification type.

## Build History (version → key deliverable)

| Version | Deliverable |
|---|---|
| 0.2.0 | Current Schema v2 baseline: `commaData` extends `RangedSummarizedExperiment`; Seqinfo-backed genome metadata; no-rownames alignment; modkit/Dorado/Megalodon import; annotation, differential methylation, enrichment, plots, vignettes, and test data are present. |
| 0.3.0 | Planned next internal feature milestone; exact scope depends on open issues and Carl's decisions. |
| 0.99.0 | Reserved for the eventual Bioconductor submission branch, not the current development line. |

For the most current architecture and issue status, prefer `AGENTS.md` and `dev/ROADMAP.md`.

## Naming Conventions

| Category | Convention | Examples |
|---|---|---|
| S4 class | `camelCase`, lowercase first | `commaData` |
| Constructor | Same as class | `commaData()` |
| Analysis functions | `verbNoun()` camelCase | `annotateSites()`, `diffMethyl()` |
| Plot functions | `plot_noun()` snake_case | `plot_volcano()`, `plot_metagene()` |
| Internal functions | `.` prefix | `.parseBetaValues()`, `.circularIndex()` |
| Arguments | `snake_case` | `mod_type`, `min_coverage` |
| Test files | `test-functionName.R` | `test-annotateSites.R` |

## Core Rules

**Always:**
- Every exported function accepts a `commaData` object as primary input
- Use `GenomicRanges::findOverlaps()` for genomic interval overlap — never nested for-loops
- Return tidy data frames (or updated `commaData`) suitable for ggplot2
- All `plot_*()` functions return a `ggplot`/`patchwork` object, not a rendered image
- Treat genome size as metadata from `Seqinfo`/`seqlengths(object)`, never hardcode
- Document every exported function with full roxygen2: `@param`, `@return`, `@examples`
- Write `testthat` tests for every exported function
- Run `devtools::test()` after changes; run `devtools::document()` after changing roxygen2

**Never:**
- Hardcode genome size, chromosome names, or organism-specific values
- Use nested R for-loops over genomic positions
- Hardcode file paths
- Import `tidyverse` — import `dplyr`, `tidyr` individually (Bioconductor requirement)
- Write stub docs (`"A dataframe."`) for `@param`/`@return`

**Performance:** vectorize with `GenomicRanges`, `zoo::rollapply()`, and matrix operations.

## Keeping Rules Current

When you add, remove, or rename anything below, update the corresponding rule file **in the same commit**:

| Change | Update |
|---|---|
| New/removed exported function | `r-source.md` — API reference |
| New/removed dependency | `r-source.md` — dependency tables |
| New/removed test file | `testing.md` — test file table |
| New/removed vignette | `documentation.md` — vignette list |
| New version / feature set shipped | `CLAUDE.md` — build history table |
| DESCRIPTION / NAMESPACE / biocViews | `bioconductor.md` — requirements checklist |
| CI/CD workflow changes | `ci-cd.md` |

See `.claude/rules/` for file-type-specific guidance loaded on demand.
