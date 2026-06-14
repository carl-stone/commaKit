# .github/AGENTS.md — GitHub automation rules

This directory controls GitHub-facing automation and repository instructions.

## Rules

- Keep CI aligned with the R version and dependency expectations documented in `dev/README.md` and `renv.lock`.
- Do not weaken checks to make a PR pass. If a check is too broad or flaky, document the reason and narrow it deliberately.
- PR templates/instructions should point code work back to Linear `commaKit Symphony` when work originates from commaBot.
- Avoid committing secrets, tokens, local paths, or machine-specific caches.

## Validation

For workflow edits, use `gh workflow list`, `gh run list`, or a PR-based dry run where possible. YAML-only edits should at least be parsed/linted locally if tooling is available.
