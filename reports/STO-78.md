# STO-78 - Strengthen Dorado Parser Integration Coverage

## Summary

Added deterministic Dorado BAM integration coverage in
`tests/testthat/test-parse_dorado.R`.

- Added in-test SAM-to-BAM fixture helpers using `Rsamtools::asBam()` and
  `Rsamtools::indexBam()`.
- Added a `.parseDorado()` happy-path fixture that asserts expected site rows,
  modification types, beta values, coverage, and count columns.
- Added parser-level malformed/truncated MM/ML skip coverage.
- Added parser-level CIGAR-unmapped call drop coverage.
- Added a CIGAR helper comparison against `GenomicAlignments` reference/query
  range utilities for supported CIGAR cases.

## Files Changed

- `tests/testthat/test-parse_dorado.R`
- `reports/STO-78.md`

## Branch And PR

- Branch: `symphony/STO-78`
- PR: not opened; `.git` is mounted read-only in this sandbox, so staging,
  commit, push, and PR creation are blocked.
- GitHub issue: `carl-stone/commaKit#279`
- Linear issue: `STO-78`

## Validation

Passed:

- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null Rscript --vanilla -e "parse(file = 'tests/testthat/test-parse_dorado.R'); cat('parse ok\n')"`
  - Result: passed; edited test file parsed successfully.
- `git diff --check`
  - Result: passed; no whitespace errors.

Attempted but blocked by git filesystem permissions:

- `git add tests/testthat/test-parse_dorado.R reports/STO-78.md && git commit -m "Add Dorado parser integration coverage"`
  - Result: failed.
  - Error: `fatal: Unable to create '/home/carl/symphony_workspaces/commakit/STO-78/.git/index.lock': Read-only file system`

Attempted but blocked by local R dependency environment:

- `Rscript -e "testthat::test_file('tests/testthat/test-parse_dorado.R')"`
  - Result: failed before tests ran.
  - Error: renv bootstrap attempted to install into
    `/home/carl/.cache/R/renv/library/STO-78-0aee050e/linux-debian-trixie/R-4.5/x86_64-pc-linux-gnu`,
    which is outside the writable sandbox and did not exist.
- `RENV_PATHS_ROOT=/tmp/renv RENV_PATHS_LIBRARY=/tmp/renv/library timeout 120 Rscript -e "options(testthat.progress.max_fails = Inf); testthat::test_file('tests/testthat/test-parse_dorado.R', reporter = 'summary')"`
  - Result: failed; timed out after 120 seconds with no output after renv
    activation.
- `RENV_PATHS_ROOT=/tmp/renv RENV_PATHS_LIBRARY=/tmp/renv/library RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE timeout 60 Rscript -e "print('started'); print(.libPaths()); print(requireNamespace('testthat', quietly = TRUE))"`
  - Result: failed; timed out before the R expression emitted output.

Not run:

- All parser tests, because the same renv activation/dependency issue blocked
  package test execution in this sandbox.

## Blockers

- Git metadata is read-only in this sandbox, preventing staging, commit, push,
  and PR creation.
- R package validation needs a restored/writable project library for the locked
  renv environment, or a runner where the repository's renv cache is already
  available.

## Risk Notes

The repository diff is test-only and bounded to Dorado parser coverage. Runtime
behavior was not changed.

## Next Owner

Factory steward / Carl reviewer.
