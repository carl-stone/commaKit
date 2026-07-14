# STO-91 Report

## Outcome

The assigned `symphony/STO-91` branch implements the validated-input contract
for the internal differential-methylation wrappers. `diffMethyl()` remains the
public preflight boundary; the three backend wrappers accept the resolved
`design_info` and no longer duplicate package checks or fallback design
resolution.

## Implementation Diff

The branch's existing STO-91 commits modify:

- `R/diffMethyl.R`: resolves the two-level design and performs dependency
  checks before dispatch; passes only the validated backend inputs.
- `R/limma_wrapper.R`, `R/quasi_f.R`, and `R/methylkit_wrapper.R`: remove
  duplicate dependency checks, former fallback inputs, and fallback design
  resolution while retaining count-matrix plumbing and result construction.
- `tests/testthat/test-diffMethyl.R`: exercises the wrappers only through
  `diffMethyl()`, checking reference/treatment direction, result columns, and
  metadata for each available backend. Existing absent-package tests target
  the public API.

The limma wrapper constructs its model matrix from the resolved condition
vector with the resolved reference level first, preserving its
treatment-minus-reference contrast after removal of its former raw design
arguments.

## Source Inspection

- GitHub issue #261 was read with
  `gh issue view 261 --repo carl-stone/commaKit --json number,title,body,url,state,labels`.
- `rg -n -C 4 "\\.run(MethylKit|Limma|QuasiF)" R tests` found `diffMethyl()` as
  the only internal caller of all three wrappers.
- The only `requireNamespace()` checks for `limma` and `methylKit` are in
  public `diffMethyl()`; no wrapper has an `is.null(design_info)` fallback.

## Validation Receipts

| Command | Result |
| --- | --- |
| `Rscript -e "styler::style_file(c('R/diffMethyl.R', 'R/limma_wrapper.R', 'R/methylkit_wrapper.R', 'R/quasi_f.R', 'tests/testthat/test-diffMethyl.R'))"` | Passed |
| `Rscript dev/precommit.R style R/diffMethyl.R R/limma_wrapper.R R/methylkit_wrapper.R R/quasi_f.R tests/testthat/test-diffMethyl.R` | Passed |
| `Rscript dev/precommit.R lint R/diffMethyl.R R/limma_wrapper.R R/methylkit_wrapper.R R/quasi_f.R tests/testthat/test-diffMethyl.R` | Passed |
| `Rscript dev/precommit.R roxygen R/diffMethyl.R R/limma_wrapper.R R/methylkit_wrapper.R R/quasi_f.R` | Passed; generated documentation is current |
| `Rscript -e "testthat::test_file('tests/testthat/test-diffMethyl.R')"` | Passed (exit code 0) |
| `git diff --check` | Passed; no whitespace errors |

An attempted `Rscript dev/precommit.R document ...` used a nonexistent
subcommand and exited nonzero without modifying files. It was immediately
replaced with the supported `roxygen` command above, which passed.

## Working Tree and Handoff

- Branch: `symphony/STO-91`.
- No staging, commits, pushes, pull requests, merges, or issue-state changes
  were performed in this worker session.
- The source/test implementation arrived as existing branch history. This
  worker leaves the refreshed `reports/STO-91.md` handoff report unstaged.
- Pre-existing factory artifacts remain untracked: `.symphony/`, `logs/`, and
  `reports/STO-91.codex-final.json`.

## Blockers

None.

## Next Owner

Trusted factory finalizer for PR publication, then independent reviewer.
