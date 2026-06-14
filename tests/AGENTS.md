# tests/AGENTS.md — commaKit test rules

This directory contains package tests. Read root `AGENTS.md` first. For files under `tests/testthat/`, also read `tests/testthat/AGENTS.md`.

## Test philosophy

- Prefer behavior-contract tests over snapshot tests.
- Test both positive paths and biologically plausible edge cases: multiple modification types, missing coverage, empty results, zero variance, circular boundaries, and annotation list-columns.
- Keep fixtures deterministic and small. Reuse helper fixtures before creating new ad hoc data.
- Strengthen weak smoke tests when touching nearby code.

## Commands

```bash
Rscript -e "devtools::test()"
Rscript -e "devtools::test(filter = 'function-or-area')"
Rscript -e "testthat::test_file('tests/testthat/test-example.R')"
```

Record the exact command and result in the worker report or PR body.
