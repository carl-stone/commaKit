# inst/AGENTS.md — installed files rules

Files in `inst/` are installed with the package.

## Rules

- Treat `inst/extdata/` as user-visible example input. Keep files small, documented, and stable.
- Scripts in `inst/scripts/` should be portable and avoid Carl-machine absolute paths.
- If installed data paths change, update examples, vignettes, and tests that use `system.file()`.
