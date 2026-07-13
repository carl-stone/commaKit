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

## Branch and Working-Tree State

`symphony/STO-91`

The branch already contains the initial STO-91 implementation commits. This
revision leaves the remaining source, public-test, and report changes unstaged
for trusted finalization. No commits, pushes, or pull requests were created by
this worker.

## Files Changed

- `R/diffMethyl.R`
- `R/limma_wrapper.R`
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
- Removed the remaining unused `coldata` and `formula` arguments from
  `.runLimma()`, whose mandatory `design_info` already represents the
  validated design.
- Asserted the backend-independent result columns for limma/quasi-F and the
  methylKit-specific q-value column through the public `diffMethyl()` API.

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

- `Rscript -e "styler::style_file(c('R/diffMethyl.R', 'R/limma_wrapper.R', 'tests/testthat/test-diffMethyl.R'))"`
  - Result: passed
  - Formatted every R file changed in this revision.

- `Rscript dev/precommit.R style R/diffMethyl.R R/limma_wrapper.R tests/testthat/test-diffMethyl.R`
  - Result: passed

- `Rscript dev/precommit.R lint R/diffMethyl.R R/limma_wrapper.R tests/testthat/test-diffMethyl.R`
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
- The working tree contains an unstaged R source and public-test diff for
  trusted finalization; factory artifacts remain untracked.

## Blockers

None. A mixed reset to expose the already-committed initial implementation as
an unstaged diff was unavailable because the sandbox makes `.git/index`
read-only. This revision instead leaves a concrete, unstaged source/test/report
diff for trusted finalization.

## Next Owner

trusted factory finalizer, then Carl reviewer
