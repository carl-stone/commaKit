# STO-64 — Set enrichment background universe and fail-loud policy

## Summary

Implemented the tested/callable ORA-universe policy for `enrichMethylation()`.
Target and regulator ORA now use distinct, role-specific backgrounds derived
from all eligible tested genes in the requested analysis slice, not from only
significant foreground genes. User-requested setup failures now stop instead
of returning warning/`NULL`-shaped biology.

## Revision: PR #318 feedback and validation repair (2026-07-13)

- Resolved the remaining automated-review thread: when `feature_type = NULL`
  finds no annotated genes, `enrichMethylation()` now raises only its setup
  error. The internal mapper retains its existing warning by default, but the
  fail-loud caller deliberately suppresses it before stopping.
- Added dependency-free coverage for the silent empty-map mode and a
  clusterProfiler-backed regression test that proves the public setup error is
  not preceded by a warning.
- Re-ran the previously failed `render-rmarkdown` precommit check; the tracked
  `getting-started.Rmd` and generated Markdown output are current.
- The earlier style-wrapper error was caused by `styler` selecting an unwritable
  default R cache. Running the exact wrapper with
  `R_CACHE_ROOTPATH=/tmp/commakit-r-cache` passes without source changes.

## Intended Diff

- `R/enrichment.R`
  - Adds `.oraUniverse()` to derive unique genes with non-missing `dm_padj`
    and `dm_delta_beta` after the existing mod-type/context, feature-type,
    overlap, and role filters.
  - Validates that every ORA foreground gene is contained in its non-empty
    universe before dispatch.
  - Makes absent annotation, unmatched requested feature types, absent
    requested role genes, and underivable ORA universes errors.
  - Retains the warning plus empty result slots only for a valid universe with
    no genes passing the requested ORA thresholds.
  - Restricts `feature_type = NULL` to target analysis because it cannot define
    a regulator class.
- `tests/testthat/test-enrichMethylation.R`
  - Adds known-answer target and regulator universe fixtures.
  - Tests tested/non-significant membership, untested exclusion,
    foreground-outside-universe failure, underivable-universe failure, and the
    valid no-signal outcome.
  - Updates public API expectations from warning/`NULL` to errors for missing
    annotations, unmatched feature types, and absent requested role genes.
- `man/enrichMethylation.Rd`, roxygen, and getting-started vignette
  - Document both ORA role backgrounds, GSEA's role-filtered ranking list, and
    setup-failure versus no-signal behavior.
- `dev/knowledge/architecture.md` and `dev/knowledge/design-decisions.md`
  - Preserve the new durable ORA contract for future contributors.

## Acceptance Criteria Evidence

- Documentation states the background for target and regulator ORA, explains
  the role-filtered GSEA list, and distinguishes setup errors from no-signal
  output in roxygen/Rd and the vignette.
- Requested prerequisites now fail loudly: annotations, feature types, roles,
  ORA-universe derivation, and foreground membership all have explicit errors.
- Known-answer coverage proves target/regulator membership and that a valid
  no-signal ORA remains a warning with empty result slots.

## Validation

- Current independent verification (2026-07-13):
  - `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "pkgload::load_all('.', quiet = TRUE); testthat::test_file('tests/testthat/test-enrichMethylation.R')"`
    - Passed: 151 expectations, 0 failures/warnings; 21 expected skips because
      optional `clusterProfiler` is unavailable.
  - `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "R.cache::setCacheRootPath('/tmp/commakit-r-cache'); styler::style_file(c('R/enrichment.R', 'tests/testthat/test-enrichMethylation.R', 'vignettes/getting-started.Rmd'))"`
    - Passed: all changed R/Rmd files were unchanged.
  - `R_LIBS_USER=/home/carl/R/library Rscript --vanilla dev/precommit.R lint R/enrichment.R tests/testthat/test-enrichMethylation.R vignettes/getting-started.Rmd`
    - Passed: no lint diagnostics.
  - `R_LIBS_USER=/home/carl/R/library Rscript --vanilla dev/precommit.R style R/enrichment.R tests/testthat/test-enrichMethylation.R vignettes/getting-started.Rmd`
    - Blocked: the `styler` dry-run check fails internally with
      ``cnd_type(): `cnd` is not a condition object`` despite the successful
      direct `styler::style_file()` result above. No source formatting changes
      are indicated; investigate the installed `styler`/precommit integration
      outside this issue's authorized scope.
  - `git diff HEAD^ HEAD --check`
    - Passed: no whitespace errors in the STO-64 implementation commit.

- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "pkgload::load_all('.', quiet = TRUE); testthat::test_file('tests/testthat/test-enrichMethylation.R')"`
  - Passed: 151 expectations, 0 failures/warnings; 21 expected skips because
    optional `clusterProfiler` is unavailable.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "R.cache::setCacheRootPath('/tmp/commakit-r-cache'); styler::style_file(c('R/enrichment.R', 'tests/testthat/test-enrichMethylation.R', 'vignettes/getting-started.Rmd'))"`
  - Passed: all changed R/Rmd files were already formatted.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "R.cache::setCacheRootPath('/tmp/commakit-r-cache'); source('dev/precommit.R')" style R/enrichment.R tests/testthat/test-enrichMethylation.R vignettes/getting-started.Rmd`
  - Passed: changed R/Rmd files are formatted.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "source('dev/precommit.R')" lint R/enrichment.R tests/testthat/test-enrichMethylation.R vignettes/getting-started.Rmd`
  - Passed: no lint diagnostics.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "roxygen2::roxygenise()"`
  - Passed: regenerated `man/enrichMethylation.Rd`.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla dev/precommit.R rmarkdown vignettes/getting-started.Rmd`
  - Passed: rendered Markdown is up to date.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "parse(file = 'R/enrichment.R'); parse(file = 'tests/testthat/test-enrichMethylation.R'); cat('R parse OK\\n')"`
  - Passed: both changed R files parse.
- `git diff --check`
  - Passed: no whitespace errors.

## Blockers and Risk

- `clusterProfiler` is not present in the available shared R library, so the
  exported clusterProfiler-backed target/regulator fixture is an expected skip
  locally. The dependency-free universe and failure-policy fixtures executed.
- The new public no-warning regression test is also an expected local skip for
  the same reason; CI must run it in its clusterProfiler-enabled image.
- `devtools` is unavailable, so the repository's `precommit.R roxygen` wrapper
  could not run. `roxygen2::roxygenise()` completed successfully and generated
  the required Rd file.

## Working Tree Notes

The prior STO-64 implementation is present in commits `c846147` and `58ea07d`
on the assigned branch. This revision leaves an unstaged corrective diff in
`R/enrichment.R`, `tests/testthat/test-enrichMethylation.R`, and this report;
no commit, push, PR publication, or issue-state change was performed. Baseline
untracked `.symphony/` and `logs/` content was not touched.

## Revision Validation

- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "R.cache::setCacheRootPath('/tmp/commakit-r-cache'); styler::style_file(c('R/enrichment.R', 'tests/testthat/test-enrichMethylation.R'))"`
  - Passed: both files were already formatted.
- `R_LIBS_USER=/home/carl/R/library R_CACHE_ROOTPATH=/tmp/commakit-r-cache Rscript --vanilla dev/precommit.R style R/enrichment.R tests/testthat/test-enrichMethylation.R`
  - Passed: the exact multi-file style wrapper completed successfully.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla dev/precommit.R lint R/enrichment.R tests/testthat/test-enrichMethylation.R`
  - Passed: no lint diagnostics.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla -e "pkgload::load_all('.', quiet = TRUE); testthat::test_file('tests/testthat/test-enrichMethylation.R')"`
  - Passed: 153 expectations, 0 failures/warnings; 22 expected skips because
    `clusterProfiler` is unavailable.
- `R_LIBS_USER=/home/carl/R/library Rscript --vanilla dev/precommit.R rmarkdown vignettes/getting-started.Rmd`
  - Passed: generated Markdown is up to date.
- `git diff --check`
  - Passed: no whitespace errors.

## Next Owner

Trusted factory finalizer: review the unstaged corrective diff, rerun the
clusterProfiler-backed test in CI, and publish/update the draft PR for GitHub
issue #292.
