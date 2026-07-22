# dev/AGENTS.md — project knowledge

This directory holds durable project knowledge, not transient worker scratch space. `knowledge/` is an OKF-style knowledge bundle; start with `knowledge/index.md` before editing durable context.

## Sources of truth

- GitHub PRs are the code integration surface.
- GitHub Issues document accepted work, bugs, and external contributor work.
- Durable knowledge belongs in `dev/knowledge/` when it will help future contributors beyond one issue.
- Every non-reserved `knowledge/*.md` concept document needs YAML frontmatter with at least `type`; see `knowledge/okf-guidelines.md`.

## Editing rules

- Update `ROADMAP.md`, `PRD.md`, or `VISION.md` only for strategic changes.
- Update `dev/knowledge/design-decisions.md` when a design choice changes how future agents should act.
- Update `dev/knowledge/known-issues.md` for reproducible bugs, gotchas, or accepted limitations.
- Do not create new task-board files in `dev/`; create GitHub issues instead.
- Archive historical material rather than deleting it when it explains past decisions.
