# STO-91 Report

## Summary

Implemented the validated-input contract for internal differential methylation
backend wrappers.

- `diffMethyl()` remains the public validation boundary for backend dependency
  checks and two-level design resolution.
- `.runLimma()`, `.runQuasiF()`, and `.runMethylKit()` now require
  `design_info` from `.resolveDiffMethylDesign()` instead of accepting
  `ref_level` or resolving a fallback design internally.
- Tightened each internal wrapper signature to accept only the validated
  inputs it uses; `diffMethyl()` is the sole internal caller and owns
  dependency validation plus design resolution.
- Removed duplicate wrapper-level dependency checks; public missing-package
  errors remain in `diffMethyl()`.
- Preserved count-matrix plumbing, reference/treatment direction, backend
  result columns, and existing statistical failure policy.
- Added public API coverage that runs the supported backends through
  `diffMethyl()` and asserts the validated factor-level direction is reflected
  in result metadata and `dm_delta_beta`.

## Intended Branch

`symphony/STO-91`

The intended source and test changes are present as unstaged working-tree
changes for trusted finalization. No commits, pushes, or pull requests were
created by this worker.

## Files Changed

- `R/diffMethyl.R`
- `R/limma_wrapper.R`
- `R/methylkit_wrapper.R`
- `R/quasi_f.R`
- `tests/testthat/test-diffMethyl.R`
- `reports/STO-91.md`

## Diff

- Removed wrapper arguments that were only needed for former fallback
  validation paths: `site_df` from limma; `site_df`, `coldata`, and `formula`
  from quasi-F; and `coldata` plus `formula` from methylKit.
- Updated `diffMethyl()` dispatch to pass the minimal validated backend inputs.
- Strengthened the public API test to pass `reference = "WT"` and assert the
  treatment-minus-reference direction and recorded comparison metadata for
  each available backend.

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
- The working tree contains the intended R source and public-test diff; only
  pre-existing factory artifacts remain untracked.

## Blockers

None.

## Next Owner

factory steward / Carl reviewer
