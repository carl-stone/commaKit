# STO-78 - Strengthen Dorado Parser Integration Coverage

## Summary

The Dorado parser now has deterministic SAM-to-BAM integration coverage. The
fixture reaches `.parseDorado()` and asserts aggregated site rows, modification
types, beta values, coverage, and count columns.

- Added private, dot-prefixed BAM fixture helpers based on `Rsamtools`.
- Compared supported CIGAR mappings against `cigarillo`, the current
  replacement for the deprecated `GenomicAlignments` range helpers; declared
  that test-only dependency in `Suggests`.
- Added explicit parser-level handling tests for truncated ML input, malformed
  MM input, and call-level drops for an insertion mapping.
- Hardened `.parseMmTag()` to reject nonnumeric or negative MM deltas before
  they can reach the read-position map.
- Rewrote the comment at `R/parse_dorado.R:361` that failed
  `commented_code_linter` during trusted finalization.

## Acceptance Evidence

| Criterion | Evidence |
| --- | --- |
| Valid fixture reaches `.parseDorado()` | The synthetic BAM test asserts two known site rows and all output count columns. |
| Error/drop behavior is explicit | Parser-level tests retain valid reads while dropping truncated ML and malformed MM reads; an insertion test retains the mapped call while dropping only the unmappable call. Direct CIGAR tests reject malformed operations. |
| Supported CIGAR cases are compared | The test file compares match, deletion, insertion, soft-clip, and skipped-reference cases with `cigarillo`. |

## Diff

- `DESCRIPTION`: add `cigarillo` to `Suggests` for the CIGAR reference test.
- `R/parse_dorado.R`: reject malformed MM delta vectors with `NULL`.
- `tests/testthat/test-parse_dorado.R`: add deterministic BAM integration and
  malformed-input coverage; resolve all PR #311 automated-review findings.

## Handoff

- Branch: `symphony/STO-78`
- Pull request: https://github.com/carl-stone/commaKit/pull/311
- GitHub issue: https://github.com/carl-stone/commaKit/issues/279
- Linear issue: `STO-78`

## Validation

Passed:

- `git diff --check origin/main`
  - No whitespace errors.
- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null R_LIBS_USER=/home/carl/agent-control-plane/state/software-factory/factory-r-library R_CACHE_ROOTPATH=/tmp/sto-78-r-cache Rscript --vanilla -e "styler::style_file(c('R/parse_dorado.R', 'tests/testthat/test-parse_dorado.R'))"`
  - Both changed R files are already formatted.
- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null R_LIBS_USER=/home/carl/agent-control-plane/state/software-factory/factory-r-library R_CACHE_ROOTPATH=/tmp/sto-78-r-cache Rscript --vanilla dev/precommit.R style R/parse_dorado.R tests/testthat/test-parse_dorado.R`
  - Changed-file style check passed.
- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null R_LIBS_USER=/home/carl/agent-control-plane/state/software-factory/factory-r-library R_CACHE_ROOTPATH=/tmp/sto-78-r-cache Rscript --vanilla dev/precommit.R lint R/parse_dorado.R tests/testthat/test-parse_dorado.R`
  - Changed-file lint passed with no lints, including `commented_code_linter`.
- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "options(testthat.progress.max_fails = Inf); pkgload::load_all('.', quiet = TRUE); testthat::test_file('tests/testthat/test-parse_dorado.R', reporter = 'summary')"`
  - Dorado parser test file passed with no warnings or failures.
- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "options(testthat.progress.max_fails = Inf); pkgload::load_all('.', quiet = TRUE); testthat::test_dir('tests/testthat', filter = 'parse', reporter = 'summary')"`
  - All parser tests passed with no warnings or failures.
- `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "invisible(parse('R/parse_dorado.R')); invisible(parse('tests/testthat/test-parse_dorado.R')); invisible(read.dcf('DESCRIPTION')); cat('R source, test syntax, and DESCRIPTION DCF: OK\\n')"`
  - Both R files parse and `DESCRIPTION` is valid DCF.

Not relied on:

- `_R_CHECK_FORCE_SUGGESTS_=false R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null R_LIBS_USER=/home/carl/R/library R CMD check --no-manual --no-build-vignettes --output=/tmp/sto-78-check .`
  - The check produced pre-existing workspace notes for hidden agent/renv
    directories and stopped while loading its temporary installation without a
    final status. The focused suites above load the current package source and
    passed; run a check from a clean export in CI for a package-level receipt.

## Next Owner

Trusted factory finalizer: publish the unstaged reviewable diff to PR #311. No
files were staged, committed, pushed, or merged.
