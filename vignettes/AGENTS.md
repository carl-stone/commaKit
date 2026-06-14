# vignettes/AGENTS.md — vignette rules

Vignettes are user-facing teaching material and part of package quality.

## Voice and scope

- Write for a new commaKit user analyzing bacterial ONT methylation data.
- Explain biological meaning before implementation detail.
- Keep examples reproducible from bundled data or clearly mark user-file workflows.
- Prefer concrete code blocks and expected outputs over abstract description.

## Technical rules

- Vignettes must build without internet access unless explicitly documented otherwise.
- Do not rely on local absolute paths.
- Keep long-running or file-dependent chunks disabled or guarded.
- If an API changes, update vignettes in the same PR as roxygen/tests.

## Validation

For vignette-only edits, run at least:

```bash
Rscript -e "devtools::build_vignettes()"
```

If pandoc or dependencies are unavailable, record the exact blocker and run targeted syntax checks where possible.
