# tests/testthat/AGENTS.md — testthat conventions

## File placement

- Mirror source function names: `R/diffMethyl.R` -> `tests/testthat/test-diffMethyl.R`.
- Shared fixtures belong in `helper-commaData-fixtures.R` when reused across files.
- Do not create hidden dependence on test order; each file should be runnable by itself.

## Expectations

- Check object classes and important values, not just that code runs.
- For S4 objects, assert assays, `rowRanges()`, `sampleInfo()`, `siteInfo()`, `modTypes()`, and `modContexts()` as relevant.
- For differential methylation, verify per-`mod_context` grouping, beta-scale effect sizes, adjusted p-values, and clear behavior on insufficient data.
- For plots, assert returned class, required aesthetics/data, and filtering behavior. Avoid fragile visual snapshots.
- For warnings/errors, match stable message fragments.

## Common traps

- Avoid assumptions that row order is biologically meaningful unless the function guarantees it.
- Do not compare floating point values with exact equality; use tolerance helpers.
- If a test needs generated docs or external files, explain why and keep it out of the fast path when possible.
