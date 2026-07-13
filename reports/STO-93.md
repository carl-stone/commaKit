# STO-93 / GH #258 — Resolve differential-methylation result table once

## Summary

`results()` now uses its already-resolved and filtered `result_data` table for
the data-frame output, matching the GRanges output path. Site metadata is
filtered to the same row index before binding, preserving output columns, row
order, and row names.

## Files changed

- `R/results_methods.R`
- `reports/STO-93.md`

## Validation

| Command | Result | Notes |
|---|---|---|
| `Rscript -e "styler::style_file('R/results_methods.R')"` | passed | Completed without changes or errors. |
| `Rscript dev/precommit.R style R/results_methods.R` | passed | Changed R source meets package style rules. |
| `Rscript dev/precommit.R lint R/results_methods.R` | passed | No lint findings. |
| `Rscript -e "testthat::test_file('tests/testthat/test-results.R')" && Rscript -e "testthat::test_file('tests/testthat/test-resultLayers.R')"` | passed | Both requested targeted test files completed successfully. |
| `git diff --check` | passed | No whitespace errors. |

## Blockers

None. An initial attempt to run both test files concurrently caused a local
`renv` bootstrap lock conflict; the required tests were then run sequentially
and passed.

## Next owner

Factory finalizer should review the unstaged diff, create the PR for
`symphony/STO-93`, and request independent review. No PR was opened by this
worker, per implementer profile constraints.
