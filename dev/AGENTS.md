# dev/AGENTS.md — project knowledge and Symphony stewardship

This directory holds durable project knowledge, not transient worker scratch space.

## Source of truth after the Symphony migration

- Linear project `commaKit Symphony` is the steward queue for Carl-facing work, milestones, decomposition, and worker dispatch.
- GitHub PRs are the code integration surface.
- GitHub Issues may still document public bugs or external contributor work, but they are no longer the primary commaBot task queue.
- Durable knowledge belongs in `dev/knowledge/` when it will help future contributors beyond one issue.

## Editing rules

- Update `ROADMAP.md`, `PRD.md`, or `VISION.md` only for strategic changes.
- Update `dev/knowledge/design-decisions.md` when a design choice changes how future agents should act.
- Update `dev/knowledge/known-issues.md` for reproducible bugs, gotchas, or accepted limitations.
- Do not create new task-board files in `dev/`; create Linear issues instead.
- Archive historical material rather than deleting it when it explains past decisions.

## Worker reports

Symphony worker reports are generated in per-issue `reports/` directories in the worker workspace/control plane, not in `dev/`. If a report contains durable package knowledge, distill that knowledge into the appropriate `dev/knowledge/` page.
