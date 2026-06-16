---
type: Agent Playbook
title: Agent Harness for commaKit
description: How coding agents should gather context, edit, validate, and report work in this repository.
resource: AGENTS.md
tags: [agents, codex, claude, harness, validation, progressive-disclosure]
timestamp: 2026-06-15T00:00:00Z
status: current
owner: Carl Stone
---

# Summary

Treat agents as maintainers with a context budget. Root `AGENTS.md` routes work; nested `AGENTS.md` files add local rules; this knowledge bundle preserves durable facts that should outlive a single issue.

# Launch Protocol

1. Run `git status --short --branch`.
2. Identify branch, issue, PR, and existing local changes.
3. Read root `AGENTS.md`.
4. Read the nearest `AGENTS.md` for directories touched.
5. Search for the relevant concept in `dev/knowledge/` before opening many files.
6. Inspect source, tests, docs, and generated artifacts before editing.
7. Make scoped changes and run targeted validation.
8. Report exact commands, results, files changed, and residual risk.

For code changes that are meant to ship, follow [Git Discipline](git-discipline.md) and [Branching and Releases](branching-releases.md): intentional commits, PR summary, tests, risks, and any issue identifier.

# Progressive Disclosure

Do not load every markdown file by default. Use this order:

1. Root `AGENTS.md` for routing and non-negotiable contracts.
2. Local `AGENTS.md` for directory rules.
3. `dev/knowledge/index.md` to choose a knowledge card.
4. A small set of source/test files found by `rg`.

When a fact gets repeated in several tool-specific files, move it into one knowledge card and replace the duplicates with links.

# Agent Surfaces

- `AGENTS.md`: durable repo and directory rules that many agents can read.
- `CLAUDE.md`: Claude-facing entry point; keep it as a pointer to shared guidance.
- `.claude/rules/*.md`: Claude on-demand hints; keep them thin and linked.
- `.agents/skills/*/SKILL.md`: Codex repo skills for reusable workflows when the directory is writable.
- `dev/knowledge/`: OKF-style durable concepts.

# Editing Discipline

- Never revert unrelated user changes.
- Do not hand-edit `man/*.Rd`.
- Keep generated files paired with their sources.
- Prefer base R and Bioconductor idioms over new dependencies.
- Add or update tests when behavior changes.
- If validation cannot run, record the blocker exactly.

# Citations

[1] [Codex best practices](https://developers.openai.com/codex/learn/best-practices)

[2] [Custom instructions with AGENTS.md](https://developers.openai.com/codex/guides/agents-md)

[3] [Agent Skills](https://developers.openai.com/codex/skills)

[4] [Introducing the Open Knowledge Format](https://cloud.google.com/blog/products/data-analytics/how-the-open-knowledge-format-can-improve-data-sharing)
