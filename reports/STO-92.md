# STO-92 — Inline `p.adjust()` in `diffMethyl()`

## Summary

Removed the no-op `.applyMultipleTesting()` helper and now call
`stats::p.adjust()` directly on `diffMethyl()`'s complete `pvalue_all` vector.
The obsolete helper source, generated Rd file, and helper-specific tests are
removed. Public behavior coverage now compares `dm_padj` with
`stats::p.adjust()` for every supported correction method; existing and
expanded public tests retain the `none` and `NA` contracts.

## Acceptance criteria

- `.applyMultipleTesting()` and wrapper-only tests are removed: satisfied by
  removing `R/multiple_testing.R`, its generated Rd file, its `DESCRIPTION`
  Collate entry, and the direct helper tests. A repository search found no
  remaining source/test/documentation references.
- Public adjusted p-values are unchanged for supported methods: protected by
  `test-diffMethyl.R`, which compares `dm_padj` to `stats::p.adjust()` for each
  value in `stats::p.adjust.methods`. The correction remains applied after all
  per-`mod_context` p-values are collected, preserving genome-wide adjustment.

## Files changed

- `R/diffMethyl.R`
- `R/multiple_testing.R` (removed)
- `DESCRIPTION`
- `R/AGENTS.md`
- `man/dot-applyMultipleTesting.Rd` (removed generated documentation)
- `tests/testthat/test-diffMethyl.R`

## Validation

- `git diff --check` — passed.
- `grep -R -n --exclude-dir=.git --exclude='*.jsonl'
  "\\.applyMultipleTesting" R tests DESCRIPTION NAMESPACE man || true` — passed;
  no remaining references.
- `timeout 30s Rscript -e "testthat::test_file('tests/testthat/test-diffMethyl.R')"`
  — unavailable (exit 124; no test output). The repository `.Rprofile` sources
  `renv/activate.R`, and that activation did not complete in this worker
  environment. `Rscript --vanilla` starts successfully, but its unactivated
  library does not contain `testthat`, so it cannot be used for package test
  validation. No host-level renv cache was configured.

## Handoff

- Intended branch: `symphony/STO-92`.
- No `.symphony/pr_feedback.json` was present.
- No files were staged, committed, pushed, or published.
- Next owner: factory finalizer should run the required targeted test in a
  functional project renv environment, then publish the reviewable PR.

## Blockers

Targeted test execution is unavailable only because project renv activation
hangs in this worker environment. The source diff and static validation are
ready for review.
