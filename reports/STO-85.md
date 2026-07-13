# STO-85 — Simplify enrichment dispatch while preserving return shape

## Summary

The branch refactors GO/KEGG and ORA/GSEA dispatch through a shared context,
call-spec builder, and method runner. This follow-up resolves the outstanding
PR #313 review finding: duplicate `method` values are de-duplicated before
dispatch, so they cannot make combined-shape assembly index a missing method.

It also repairs each reported style/lint location without changing the public
enrichment API or its universe and role semantics.

## Intended Diff

- `R/enrichment.R`: store unique methods in
  `.enrichmentDispatchContext()` while preserving the first-seen order.
- `tests/testthat/test-enrichMethylation.R`: add a dependency-free regression
  test that passes `c("ora", "ora")`, verifies a single-method return shape,
  and checks the existing no-significant-genes warning.
- `R/enrichment.R`: use snake-case names inside the private dispatch-context
  helper while preserving its existing context fields and documented public
  `OrgDb`, `TERM2GENE`, and `TERM2NAME` arguments.
- `R/enrichment.R`: remove stale commented code and use targeted lint
  suppressions where the repository formatter and indentation linter conflict.
- `tests/testthat/test-enrichMethylation.R`: remove stale/commented-code text
  and reformat the reported assertion continuation.
- `tests/testthat/test-buildKEGGTermGene.R`: remove a stale file-name comment
  so the full changed-file lint check is clean.

The prior branch commits contain the main dispatch refactor and its call-spec
coverage. The current working-tree diff contains only this review follow-up.

## Acceptance Criteria Evidence

- Supported dispatch combinations retain their return shape: the shared
  dispatcher continues to return a single `list(go, kegg)` for one unique
  method and the nested `ora`/`gsea` shape for both unique methods. The new
  regression covers the duplicate single-method path without clusterProfiler.
- Optional-dependency branches have behavior-focused coverage: existing tests
  assert custom `TERM2GENE`/`TERM2NAME`, `OrgDb`, `organism`, KEGG mapping,
  empty KEGG `TERM2NAME`, and combined short-circuit call specs/shapes.

## Validation

- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "pkgload::load_all('.', quiet = TRUE); testthat::test_file('tests/testthat/test-enrichMethylation.R')"`
  - Passed: 141 expectations; 19 expected skips because `clusterProfiler` is
    not installed.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "pkgload::load_all('.', quiet = TRUE); testthat::test_file('tests/testthat/test-buildKEGGTermGene.R')"`
  - Passed: no failures; 20 expected skips because `KEGGREST` is not
    installed.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "R.cache::setCacheRootPath('/tmp/commakit-r-cache'); source('dev/precommit.R')" style R/enrichment.R tests/testthat/test-enrichMethylation.R tests/testthat/test-buildKEGGTermGene.R`
  - Passed: all three branch-touched R files were already formatted.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "source('dev/precommit.R')" lint R/enrichment.R tests/testthat/test-enrichMethylation.R tests/testthat/test-buildKEGGTermGene.R`
  - Passed: no lint diagnostics.
- `Rscript -e "testthat::test_file('tests/testthat/test-enrichMethylation.R')"`
  - Passed: exact required validation command exited with status 0.
- `Rscript -e "testthat::test_file('tests/testthat/test-buildKEGGTermGene.R')"`
  - Passed: exact required validation command exited with status 0.
- `Rscript --vanilla -e "parse('R/enrichment.R'); parse('tests/testthat/test-enrichMethylation.R'); parse('tests/testthat/test-buildKEGGTermGene.R'); cat('parse complete\\n')"`
  - Passed; all branch-touched R files parse successfully.
- `git diff --check`
  - Passed; no whitespace errors.

## Blockers

None. Workspace-source tests used the available shared R library and
`pkgload::load_all()`; style validation used a process-local cache under
`/tmp`. No host-level cache or repository configuration was changed.

## Working Tree Notes

Baseline untracked files were `.symphony/`, `logs/`, and
`reports/STO-85.codex-final.json`; they are untouched. Changes are unstaged.
No commit, push, PR publication, or issue-state change was performed.

## Next Owner

Trusted factory finalizer: review the unstaged follow-up diff and publish the
existing branch/PR when ready.
