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
  - Baseline: `## symphony/STO-91...origin/main`
  - Initial status included modified target files and pre-existing untracked
    `logs/`.

- `rg -n "\\.run(MethylKit|Limma|QuasiF)|resolveDiffMethylDesign|requireNamespace\\(\\\"(limma|methylKit)\\\"" R tests/testthat/test-diffMethyl.R`
  - Result: passed
  - Confirmed internal wrapper calls come from `diffMethyl()` and dependency
    checks remain at the public boundary.

- `Rscript --vanilla -e "for (f in c('R/diffMethyl.R','R/limma_wrapper.R','R/methylkit_wrapper.R','R/quasi_f.R','tests/testthat/test-diffMethyl.R')) parse(f); cat('parse ok\n')"`
  - Result: passed
  - Output: `parse ok`

- `git diff --check`
  - Result: passed
  - Output: none

- `Rscript --vanilla -e "cat('R ok\n'); print(requireNamespace('testthat', quietly=TRUE)); print(requireNamespace('limma', quietly=TRUE)); print(requireNamespace('methylKit', quietly=TRUE))"`
  - Result: passed
  - Output showed R runs, but the non-renv/global library lacks
    `testthat`, `limma`, and `methylKit`.

- `Rscript -e "testthat::test_file('tests/testthat/test-diffMethyl.R')"`
  - Result: failed/unavailable
  - Attempt 1 bootstrapped `renv` and then left a stale PTY with no visible R
    process or test result; interrupted after waiting.

- `timeout 120s Rscript -e "testthat::test_file('tests/testthat/test-diffMethyl.R')"`
  - Result: failed/unavailable
  - Output reached `renv` bootstrap/install only:
    `# Bootstrapping renv 1.1.8`, `Downloading renv ... OK`,
    `Installing renv ... OK`.
  - Exit code: `124`
  - The command timed out before `testthat` ran.

- `RENV_CONFIG_CACHE_ENABLED=FALSE RENV_PATHS_LIBRARY=.symphony/renv/library/linux-debian-trixie/R-4.5/x86_64-pc-linux-gnu timeout 600s Rscript --vanilla -e ".libPaths(c(file.path(getwd(), '.symphony/renv/library/linux-debian-trixie/R-4.5/x86_64-pc-linux-gnu'), .libPaths())); library(renv); renv::restore(project = getwd(), prompt = FALSE)"`
  - Result: failed/unavailable
  - Attempted a workspace-local dependency restore with cache disabled.
    The session became stale after loading `renv`; no R process was visible
    and no restore result was produced. Generated `.symphony/` files were
    removed because they were validation noise, not intended source diff.

## Validation

Passed:

- Changed R files and the changed test file parse successfully under
  `Rscript --vanilla`.
- `git diff --check` reports no whitespace errors.

Unavailable:

- Required focused validation
  `Rscript -e "testthat::test_file('tests/testthat/test-diffMethyl.R')"`
  could not complete in this sandbox because repository `renv` activation
  repeatedly timed out or produced a stale PTY during bootstrap before
  `testthat` could run. The global/non-renv R library also lacks `testthat`,
  `limma`, and `methylKit`.

## Blockers

- Focused `test-diffMethyl.R` validation needs a restored commaKit R
  dependency environment. Safest next step is to rerun the required test in
  the factory/finalizer environment where `renv` can complete or dependencies
  are already available.

## Next Owner

factory steward / Carl reviewer

