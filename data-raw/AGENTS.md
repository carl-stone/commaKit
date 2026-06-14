# data-raw/AGENTS.md — generated data rules

`data-raw/` contains scripts that create package data.

## Rules

- Preserve deterministic seeds and document them in the script.
- Do not edit generated `.rda`/package data without updating the generating script.
- Keep generated example data small enough for package checks and CI.
- When fixture biology changes, update tests, documentation, and root `AGENTS.md` test-data notes together.

## Validation

After changing data-generation scripts, run the script, inspect the resulting data objects, and run affected tests. Record whether generated files changed intentionally.
