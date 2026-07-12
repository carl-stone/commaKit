# STO-85 — Simplify enrichment dispatch while preserving return shape

## Summary

Refactored `enrichMethylation()` enrichment dispatch so ORA/GSEA and GO/KEGG
calls share a single setup path. The public return shape is preserved for
single methods, combined `method = c("ora", "gsea")`, `$go`/`$kegg` slots,
feature-type lists, and `gene_role = "both"`.

## Intended Branch

`symphony/STO-85`

## Files Changed

- `R/enrichment.R`
- `tests/testthat/test-buildKEGGTermGene.R`
- `tests/testthat/test-enrichMethylation.R`
- `reports/STO-85.md`

## Implementation Notes

- Added `.enrichmentDispatchContext()` to collect the shared enrichment
  parameters once per `enrichMethylation()` call.
- Added `.enrichmentCallSpec()` to describe the clusterProfiler function and
  arguments for each GO/KEGG and ORA/GSEA combination before dispatch.
- Added `.runEnrichmentMethod()` and `.runClusterProfilerCall()` so
  `.runEnrichmentForGeneMap()` only coordinates method shape.
- Preserved existing ORA/GSEA warnings for no significant genes and no valid
  gene scores.
- Preserved KEGG `TERM2NAME` behavior where an empty KEGG term-name table is
  treated like `NULL`.
- Updated `enrichMethylation()` call sites to pass the shared dispatch context.
- Follow-up for PR #313 feedback: wrapped long lines in the changed source and
  test files so changed-file lint/style checks no longer trip on line length.

## Acceptance Criteria Mapping

- All supported dispatch combinations retain current return shape:
  - Single method still returns `list(go = ..., kegg = ...)`.
  - Combined methods still return `list(go = list(ora = ..., gsea = ...),
    kegg = list(ora = ..., gsea = ...))`.
  - Feature-type and `gene_role = "both"` wrapping remains outside the shared
    dispatcher and is unchanged.
- Optional-dependency branches have behavior-focused coverage:
  - Tests now assert call specs for custom `TERM2GENE`/`TERM2NAME`.
  - Tests now assert call specs for optional `OrgDb` and `organism` dispatch.
  - Tests now assert empty KEGG `TERM2NAME` tables are omitted.
  - Tests now assert combined short-circuit return shape without needing to call
    clusterProfiler.

## Validation

- `git diff --check`
  - Result: passed.
  - Notes: no whitespace errors.
- `Rscript --vanilla -e "print('r-ok')"`
  - Result: passed.
  - Notes: R starts correctly when repository renv activation is bypassed.
- `timeout 180 Rscript -e "testthat::test_file('tests/testthat/test-enrichMethylation.R')"`
  - Result: failed.
  - Notes: timed out with exit code 124 and no output before test execution
    reached a report.
- `timeout 180 Rscript -e "testthat::test_file('tests/testthat/test-buildKEGGTermGene.R')"`
  - Result: failed.
  - Notes: timed out with exit code 124 and no output before test execution
    reached a report.
- `timeout 120 Rscript dev/precommit.R lint R/enrichment.R tests/testthat/test-enrichMethylation.R tests/testthat/test-buildKEGGTermGene.R`
  - Result: failed.
  - Notes: timed out with exit code 124 and no output before lint execution
    reached a report.
- `timeout 120 Rscript dev/precommit.R style R/enrichment.R tests/testthat/test-enrichMethylation.R tests/testthat/test-buildKEGGTermGene.R`
  - Result: failed.
  - Notes: timed out with exit code 124 and no output before style execution
    reached a report.

## Blockers

Runtime validation through repository-managed R is blocked in this sandbox:
commands that load `.Rprofile` and activate `renv` timed out with no output.
`Rscript --vanilla -e "print('r-ok')"` passes, which narrows the issue to
repository startup/dependency activation rather than the R executable itself.

## Working Tree Notes

Baseline status included untracked `.symphony/`, `logs/`, and
`reports/STO-85.codex-final.json`. This handoff leaves the intended source
diff in place and does not stage, commit, push, or open a PR.

## Next Owner

factory steward / Carl reviewer
