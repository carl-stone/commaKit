# AGENTS.md - commaKit Agent Router

This is the first file coding agents should read at the repository root. Keep it short, current, and linked. Detailed knowledge lives in the OKF-style bundle under [dev/knowledge/](dev/knowledge/index.md).

## Start Here

1. Run `git status --short --branch` before editing.
2. Identify whether the work is tied to a Linear issue, GitHub issue, branch, or PR. If none is provided, say so in your final notes.
3. Read the nearest `AGENTS.md` for every directory you touch.
4. Search before assuming: use `rg` or `rg --files` for source, tests, and docs.
5. Make the smallest coherent change that satisfies the request.
6. Run targeted validation and record the exact command/result.
7. Do not stage, revert, or overwrite unrelated user or agent changes.
8. All changes go through PRs. Never push or merge changes directly to `main` without a PR.
9. When a commit fully implements a GitHub issue, include a closing keyword (`Closes #N`, `Fixes #N`, or `Resolves #N`) in the commit message. A bare `(#N)` reference does not close the issue.

## Project Snapshot

- Package: `commaKit` (Comparative Microbial Methylomics Analysis Kit)
- Current version: `0.2.0`
- Language/ecosystem: R >= 4.3.0, Bioconductor style, CI pinned to R 4.5
- License: MIT
- GitHub: `carl-stone/commaKit`
- Primary object: `commaData`, an S4 class extending `RangedSummarizedExperiment`

Read [dev/knowledge/project-status.md](dev/knowledge/project-status.md) for the current source/test inventory and project-state notes.

## Context Map

| Task | Read next |
|---|---|
| R source, S4 class, parser, analysis, plotting, accessors | [R/AGENTS.md](R/AGENTS.md), [dev/knowledge/architecture.md](dev/knowledge/architecture.md) |
| Tests or fixtures | [tests/AGENTS.md](tests/AGENTS.md), [tests/testthat/AGENTS.md](tests/testthat/AGENTS.md), [dev/knowledge/test-quality.md](dev/knowledge/test-quality.md) |
| Vignettes or user tutorials | [vignettes/AGENTS.md](vignettes/AGENTS.md) |
| Generated data | [data-raw/AGENTS.md](data-raw/AGENTS.md) |
| Installed example files | [inst/AGENTS.md](inst/AGENTS.md) |
| Generated Rd docs | [man/AGENTS.md](man/AGENTS.md) |
| GitHub workflows or PR automation | [.github/AGENTS.md](.github/AGENTS.md), [dev/knowledge/git-discipline.md](dev/knowledge/git-discipline.md) |
| Strategy, roadmap, PRD, durable project knowledge | [dev/AGENTS.md](dev/AGENTS.md), [dev/knowledge/index.md](dev/knowledge/index.md) |
| Agent operating model and context hygiene | [dev/knowledge/agent-harness.md](dev/knowledge/agent-harness.md), [dev/knowledge/okf-guidelines.md](dev/knowledge/okf-guidelines.md) |

## Agent skills

Repository-owned Codex procedures live under `.codex/skills/`:

- `linear`: query and mutate the active Linear issue through Symphony's
  `linear_graphql` tool.
- `pull`: merge the latest remote branch and `origin/main` without rebasing.
- `commit`: stage intentional paths and create a scoped, validated commit.
- `push`: publish the branch and create or refresh its pull request.
- `land`: resolve feedback and conflicts, watch checks, and squash-merge only
  when the PR is ready.

Load the matching skill before performing any of these procedures. The skills
supplement, rather than replace, the repository and directory-level
`AGENTS.md` instructions.

### Issue tracker

Internal commaBot work lives in Linear project `commaKit Symphony`. GitHub
Issues remain the public bug and external-contributor tracker, and external PRs
are also a triage request surface. See `docs/agents/issue-tracker.md`.

### Triage labels

Use the repository's current triage labels: `awaiting-factory-triage`,
`needs-info`, `needs-decision`, `ready-for-agent`, `ready-for-human`, and
`wontfix`. See `docs/agents/triage-labels.md`.

### Domain docs

Use the single-context domain-doc layout: root `CONTEXT.md` plus root `docs/adr/` when they exist. See `docs/agents/domain.md`.

## Core Contracts

- `commaData` stores genomic sites in `rowRanges()` as 1-bp `GRanges`; do not move genomic coordinates back into ordinary row-data columns.
- Required raw assays include `methylation` and `coverage`; count-aware paths may use `mod_counts`, `canonical_counts`, and `other_mod_counts`.
- Use `siteCoverage(object)` for the coverage assay. `coverage(commaData)` is a deprecated compatibility wrapper because `IRanges::coverage()` has a different Bioconductor meaning.
- `mod_type` is the chemical modification (`6mA`, `5mC`, `4mC`). `motif` is the sequence context. `mod_context` is computed as `paste(mod_type, motif, sep = "_")`, for example `6mA_GATC`; if motif is unavailable, it falls back to `mod_type`.
- Analyze and test differential methylation by `mod_context`, not just `mod_type`.
- Align sites with `GenomicRanges::findOverlaps()` plus explicit `mod_type`/motif checks. `site_key` is a computed display label such as `chr1:512:+:6mA:GATC`, not an alignment key.
- `annotateSites()` intentionally stores all associations in list-columns (`CharacterList`, `IntegerList`, `NumericList`). Do not collapse to a single closest feature.
- Effect sizes are reported on the beta scale (0-1), and multiple-testing correction is genome-wide across all tested contexts.
- Genome sizes come from `Seqinfo`/`seqlengths(object)`, never from hardcoded organisms.
- Do not add `tidyverse`. Prefer base R, Bioconductor APIs, or explicit narrow dependencies.

## Validation Commands

```bash
Rscript -e "devtools::test()"
Rscript -e "devtools::test(filter = 'annotateSites')"
Rscript -e "testthat::test_file('tests/testthat/test-annotateSites.R')"
Rscript dev/test-timing.R
Rscript -e "devtools::document()"
Rscript -e "styler::style_pkg()"
Rscript -e "devtools::check(build_args = c('--no-build-vignettes'))"
Rscript -e "devtools::check()"
```

Use the narrowest command that proves the change. Run `styler::style_pkg()` before committing to keep code formatted (CI enforces the tidyverse style guide with 2-space indentation and 80-char line width). Run the narrowest relevant tests before committing (for example `Rscript -e "testthat::test_file('tests/testthat/test-annotateSites.R')" or `Rscript -e "devtools::test()"` for broader changes). The repository's `pre-push` hook runs both `devtools::test()` and `dev/test-timing.R`; do not bypass it. CI uploads the timing report as an artifact on every push and PR. Escalate to `devtools::document()` for roxygen edits, and to check-level validation for package metadata, examples, vignettes, dependencies, or broad API changes.

## Knowledge Format

`dev/knowledge/` is organized as an Open Knowledge Format v0.1-style bundle: markdown files with YAML frontmatter, directory `index.md` files for progressive disclosure, ordinary markdown links as graph edges, and optional `log.md` history. See [dev/knowledge/okf-guidelines.md](dev/knowledge/okf-guidelines.md).

When a durable fact changes, update the relevant knowledge card in the same change. If a fact is temporary or issue-specific, keep it out of the durable bundle unless it will help future agents.
