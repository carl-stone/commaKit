---
paths:
  - ".github/**/*.yml"
  - ".github/**/*.yaml"
---

# Git and CI/CD

## Branch Strategy

- **Stable:** `main` — do not push here directly; work through PRs
- **Dev branches:** use tool-specific prefixes such as `codex/<description>` or `claude/<description>-<id>` for AI-initiated work

## Commit Style

Use descriptive, imperative messages:

```
Add commaData S4 class definition and show() method
Fix circular genome arithmetic to use genomeInfo slot
Replace nested for-loops in annotateSites with findOverlaps
Implement plot_tss_profile() with loess smooth overlay
```

## CI/CD Workflows

- **`R-CMD-check.yaml`** — runs `R CMD check --no-manual --as-cran` on R 4.5 / Ubuntu. CI runs examples with `--run-donttest`, so use `\dontrun{}` for examples that require user-provided files.
- **`pkgdown.yaml`** — builds the pkgdown site on PRs and deploys to `gh-pages` from `main`.
- **`render-rmarkdown.yaml`** — auto-renders changed `.Rmd` files on push and commits generated `.md` outputs when present.

Keep local checks aligned with the workflows above.
