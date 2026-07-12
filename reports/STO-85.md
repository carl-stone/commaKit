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
- `timeout 30 Rscript -e "print('r-ok')"`
  - Result: failed.
  - Notes: timed out after 30 seconds with no output. This indicates Rscript
    does not start in the sandbox before package or test code is reached.
- `timeout 300 Rscript -e "testthat::test_file('tests/testthat/test-enrichMethylation.R', reporter = 'summary')"`
  - Result: failed.
  - Notes: timed out after 300 seconds with no output.
- `timeout 300 Rscript -e "testthat::test_file('tests/testthat/test-buildKEGGTermGene.R', reporter = 'summary')"`
  - Result: failed.
  - Notes: timed out after 300 seconds with no output.

## Blockers

Runtime validation is blocked in this sandbox because `Rscript` does not start:
even a trivial `Rscript -e "print('r-ok')"` command timed out with no output.

## Working Tree Notes

Baseline status already included modified `R/enrichment.R`,
`tests/testthat/test-enrichMethylation.R`, and untracked `.symphony/` plus
`logs/`. This handoff leaves the intended source diff in place and does not
stage, commit, push, or open a PR.

## Next Owner

factory steward / Carl reviewer
