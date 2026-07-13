# STO-91 Report

## Summary

Implemented the validated-input contract for internal differential methylation
backend wrappers.

- `diffMethyl()` remains the public validation boundary for backend dependency
  checks and two-level design resolution.
- `.runLimma()`, `.runQuasiF()`, and `.runMethylKit()` now require
  `design_info` from `.resolveDiffMethylDesign()` instead of accepting
  `ref_level` or resolving a fallback design internally.
- Removed duplicate wrapper-level dependency checks; public missing-package
  errors remain in `diffMethyl()`.
- Preserved count-matrix plumbing, reference/treatment direction, backend
  result columns, and existing statistical failure policy.
- Added public API coverage that runs the supported backends through
  `diffMethyl()` and asserts the validated factor-level direction is reflected
  in result metadata and `dm_delta_beta`.

## Intended Branch

`symphony/STO-91`

## Files Changed

- `R/diffMethyl.R`
- `R/limma_wrapper.R`
- `R/methylkit_wrapper.R`
- `R/quasi_f.R`
- `tests/testthat/test-diffMethyl.R`
- `reports/STO-91.md`

## Commands Run

- `git status --short --branch`
  - Result: passed
  - Branch: `symphony/STO-91`.
  - Existing factory artifacts remain untracked: `.symphony/`, `logs/`, and
    `reports/STO-91.codex-final.json`.

- `rg -n "\\.run(MethylKit|Limma|QuasiF)|resolveDiffMethylDesign|requireNamespace\\(\\\"(limma|methylKit)\\\"" R tests/testthat/test-diffMethyl.R`
  - Result: passed
  - Confirmed internal wrapper calls come from `diffMethyl()` and dependency
    checks remain at the public boundary.

- `Rscript -e "styler::style_file(c('R/diffMethyl.R', 'R/limma_wrapper.R', 'R/methylkit_wrapper.R', 'R/quasi_f.R', 'tests/testthat/test-diffMethyl.R'))"`
  - Result: passed
  - Formatted every changed R file without introducing a source diff.

- `Rscript dev/precommit.R style R/diffMethyl.R R/limma_wrapper.R R/methylkit_wrapper.R R/quasi_f.R tests/testthat/test-diffMethyl.R`
  - Result: passed

- `Rscript dev/precommit.R lint R/diffMethyl.R R/limma_wrapper.R R/methylkit_wrapper.R R/quasi_f.R tests/testthat/test-diffMethyl.R`
  - Result: passed

- `git diff --check`
  - Result: passed
  - Output: none

- `Rscript -e "testthat::test_file('tests/testthat/test-diffMethyl.R')"`
  - Result: passed (exit code 0).
  - The initial renv bootstrap completed before this successful rerun.

## Validation

Passed:

- Formatter and pre-commit style/lint checks pass for every changed R file.
- Required focused public API test passes.
- `git diff --check` reports no whitespace errors.

## Blockers

None.

## Next Owner

factory steward / Carl reviewer
