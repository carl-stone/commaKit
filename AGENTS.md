# .github/AGENTS.md — GitHub automation rules

This directory controls GitHub-facing automation and repository
instructions. For branch, PR, and validation expectations, see
`../dev/knowledge/git-discipline.md` and
`../dev/knowledge/branching-releases.md`.

## Rules

- Keep CI aligned with the R version and dependency expectations
  documented in `dev/README.md` and `renv.lock`.
- Do not weaken checks to make a PR pass. If a check is too broad or
  flaky, document the reason and narrow it deliberately.
- Keep pull requests self-contained and link the relevant GitHub issue
  when one exists.
- Avoid committing secrets, tokens, local paths, or machine-specific
  caches.
- Use closing keywords (`Closes #N`, `Fixes #N`, `Resolves #N`) in
  commit messages when a commit fully implements a GitHub issue. A bare
  `(#N)` reference does NOT close the issue. The PR template includes a
  “Closes” section as a reminder.

## Validation

For workflow edits, use `gh workflow list`, `gh run list`, or a PR-based
dry run where possible. YAML-only edits should at least be parsed/linted
locally if tooling is available.
