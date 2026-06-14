# R/AGENTS.md — commaKit source code rules

This directory contains package implementation files. Read the root `AGENTS.md` first, then follow these rules for any file under `R/`.

## Source architecture

- `commaData_class.R`, `commaData_constructor.R`, and `accessors.R` define the core S4/RangedSummarizedExperiment contract.
- Parser files (`parse_modkit.R`, `parse_dorado.R`, `parse_megalodon.R`) normalize external caller formats into the shared constructor path.
- Analysis files (`diffMethyl.R`, `quasi_f.R`, `limma_wrapper.R`, `methylkit_wrapper.R`, `multiple_testing.R`, `results_methods.R`, `result_layers.R`) must preserve per-`mod_context` analysis and genome-wide multiple-testing semantics.
- Plot files must return ggplot/patchwork objects except `plot_heatmap()`, which returns ComplexHeatmap output.

## Editing rules

1. Preserve S4 validity. If a slot, accessor, or constructor behavior changes, update validity checks, accessors, tests, and docs together.
2. Align genomic sites with `GenomicRanges::findOverlaps()` and explicit `mod_type`/motif checks. Do not align assays by row names or display `site_key` strings.
3. Keep modification-type agnostic behavior. Never assume only 6mA or only one motif/context.
4. Use explicit namespace calls where masks are known: `dplyr::rename()`, `dplyr::count()`, `S4Vectors::DataFrame()`, etc.
5. Prefer base `vapply()`/`lapply()` and `purrr::list_rbind()` over superseded `map_dfr()`/`map_dbl()`.
6. Avoid non-ASCII in R source, roxygen, and examples.
7. Do not hand-edit generated `man/*.Rd`; run `Rscript -e "devtools::document()"` when roxygen changes.

## Roxygen expectations

- Every exported function needs `@param`, `@return`, `@examples`, and `@export` where appropriate.
- Examples that require user files, external annotation databases, internet, or heavy data belong in `\dontrun{}` rather than `\donttest{}` because CI runs donttest examples.
- Keep examples small and based on bundled example data where possible.

## Validation

Use the smallest validation that proves the change:

```bash
Rscript -e "devtools::test(filter = 'diffMethyl')"
Rscript -e "devtools::document()"
Rscript -e "devtools::check(build_args = c('--no-build-vignettes'))"
```

Escalate to full `devtools::check()` when exported API, examples, dependencies, vignettes, or package metadata change.
