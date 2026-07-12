# STO-78 - Strengthen Dorado Parser Integration Coverage

## Summary

Added deterministic Dorado BAM integration coverage in
`tests/testthat/test-parse_dorado.R` and addressed PR #311 review feedback.

- Added in-test SAM-to-BAM fixture helpers using `Rsamtools::asBam()` and
  `Rsamtools::indexBam()`.
- Dot-prefixed the test helpers as `.make_dorado_test_bam()` and
  `.make_sam_record()` to avoid shared testthat environment name collisions.
- Added a `.parseDorado()` happy-path fixture that asserts expected site rows,
  modification types, beta values, coverage, and count columns.
- Added parser-level malformed/truncated MM/ML skip coverage with an on-sequence
  truncated ML fixture.
- Added parser-level CIGAR-unmapped call drop coverage.
- Added a CIGAR helper comparison against `GenomicAlignments` reference/query
  range utilities for supported CIGAR cases.
- Added `GenomicAlignments` to `Suggests` so the reference-comparison coverage
  is available in normal package test environments.
- Cleaned up the Dorado test section comments and long malformed-CIGAR
  expectation shape to address the failed `lintr` check without changing test
  behavior.

## Files Changed

- `DESCRIPTION`
- `tests/testthat/test-parse_dorado.R`
- `reports/STO-78.md`

## Branch And PR

- Branch: `symphony/STO-78`
- PR: https://github.com/carl-stone/commaKit/pull/311
- GitHub issue: `carl-stone/commaKit#279`
- Linear issue: `STO-78`

## Validation

Passed:

- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null Rscript --vanilla -e "parse(file = 'tests/testthat/test-parse_dorado.R'); read.dcf('DESCRIPTION'); cat('parse/dcf ok\n')"`
  - Result: passed; edited test file parses and `DESCRIPTION` remains valid
    DCF metadata.
- `git diff --check`
  - Result: passed; no whitespace errors.
- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null R_LIBS_USER=$PWD/.symphony/renv-library/linux-debian-trixie/R-4.5/x86_64-pc-linux-gnu:$PWD/.symphony/renv/library/linux-debian-trixie/R-4.5/x86_64-pc-linux-gnu Rscript --vanilla -e "cat('libPaths=', paste(.libPaths(), collapse='|'), '\n'); for (pkg in c('testthat','Rsamtools','GenomicAlignments','lintr')) cat(pkg, requireNamespace(pkg, quietly=TRUE), '\n')"`
  - Result: passed; command completed and confirmed the workspace-local library
    is present but contains none of `testthat`, `Rsamtools`,
    `GenomicAlignments`, or `lintr`.

Attempted but blocked by local R dependency environment:

- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null Rscript --vanilla -e "options(testthat.progress.max_fails = Inf); testthat::test_file('tests/testthat/test-parse_dorado.R', reporter = 'summary')"`
  - Result: failed before tests ran.
  - Error: `there is no package called 'testthat'` in the system library when
    renv activation is disabled.
- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null Rscript --vanilla -e "options(testthat.progress.max_fails = Inf); testthat::test_dir('tests/testthat', filter = 'parse', reporter = 'summary')"`
  - Result: failed before tests ran.
  - Error: `there is no package called 'testthat'` in the system library when
    renv activation is disabled.
- `timeout 180 Rscript -e "options(testthat.progress.max_fails = Inf); testthat::test_file('tests/testthat/test-parse_dorado.R', reporter = 'summary')"`
  - Result: failed; timed out during repository renv startup after bootstrapping
    `renv`.
- `timeout 180 Rscript -e "options(testthat.progress.max_fails = Inf); testthat::test_dir('tests/testthat', filter = 'parse', reporter = 'summary')"`
  - Result: failed during concurrent renv bootstrap.
  - Error: failed to lock the project renv library because the Dorado test
    command was already bootstrapping renv.
- `Rscript -e "cat('libPaths:', paste(.libPaths(), collapse='|'), '\n'); cat('testthat=', requireNamespace('testthat', quietly=TRUE), '\n'); cat('Rsamtools=', requireNamespace('Rsamtools', quietly=TRUE), '\n'); cat('GenomicAlignments=', requireNamespace('GenomicAlignments', quietly=TRUE), '\n')"`
  - Result: stopped after hanging for more than 60 seconds in repository renv
    startup.
- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null Rscript --vanilla -e "cat('libPaths=', paste(.libPaths(), collapse='|'), '\n'); for (pkg in c('testthat','Rsamtools','GenomicAlignments','lintr')) cat(pkg, requireNamespace(pkg, quietly=TRUE), '\n')"`
  - Result: passed; command completed and confirmed the system library lacks
    `testthat`, `Rsamtools`, `GenomicAlignments`, and `lintr` when renv
    activation is disabled.

Not run successfully:

- Dorado parser test file and all parser tests, because `testthat` is not
  available without renv and the repository renv startup did not complete in
  this sandbox.

## Blockers

- R package validation needs a restored/writable project library for the locked
  renv environment, or a runner where the repository's renv cache is already
  available.

## Risk Notes

The repository diff is bounded to test coverage and package test metadata.
Runtime parser behavior was not changed.

## Next Owner

Factory steward / Carl reviewer.
