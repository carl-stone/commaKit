---
paths:
  - ".github/**/*.yml"
  - ".github/**/*.yaml"
---

# Git and CI/CD

## Branch Strategy

- **Stable:** `main` — do not push here directly; work through PRs
- **Dev branches:** `claude/<description>-<id>` naming pattern for AI-initiated work

## Commit Style

Use descriptive, imperative messages:

```
Add commaData S4 class definition and show() method
Fix circular genome arithmetic to use genomeInfo slot
Replace nested for-loops in annotateSites with findOverlaps
Implement plot_tss_profile() with loess smooth overlay
```

## CI/CD Workflows

- **`r.yml`** — runs `rcmdcheck` on push/PR against R 3.6.3 and 4.1.1 on macOS-latest. Keep the package passing `R CMD check` throughout development.
- **`render-rmarkdown.yaml`** — auto-renders `.Rmd` files on push (README.Rmd → README.md).
