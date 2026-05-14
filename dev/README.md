# comma Development Directory

This directory is the project-management and knowledge base for comma.

**Maintained by:** commaBot  
**Audience:** Carl, Claire, future agents, and any human contributor

---

## Start Here

If you want to know...

| Question | Read this |
|----------|-----------|
| What is happening right now? | `STATUS.md` |
| What work needs to be done? | `BACKLOG.md` |
| Where is the package going? | `ROADMAP.md` |
| What is v1.0 supposed to be? | `PRD.md` |
| What is the dream version? | `VISION.md` |
| What do we know about test quality? | `knowledge/test-quality.md` |
| What bugs/gotchas/edge cases are known? | `knowledge/known-issues.md` |
| Why is the package designed this way? | `knowledge/design-decisions.md` |
| How should Claire learn the package? | `ONBOARDING.md` |

---

## File Roles

### `STATUS.md` — Daily View

Short, current, operational.

Contains:
- What is in progress now
- What is next up
- What is blocked
- What was recently completed

Do **not** put every task here. That's what `BACKLOG.md` is for.

---

### `BACKLOG.md` — Single Source of Truth for Work Items

Prioritized list of all tasks, bugs, improvements, documentation work, and submission items.

Each item has:
- ID
- Title
- Priority
- Status
- Size
- Source
- Detailed problem statement

This replaces scattered TODO tables in audit notes.

---

### `ROADMAP.md` — Strategic Direction

Longer-term direction, not tactical task tracking.

Contains:
- v1.0 status
- Bioconductor submission path
- v1.1+ future feature roadmap

If you want to know where comma is going, read this. If you want to know what to do tomorrow, read `BACKLOG.md`.

---

### `knowledge/` — Durable Project Knowledge

Organized by topic, not by date.

- `test-quality.md` — what tests are strong, weak, or missing
- `known-issues.md` — bugs, gotchas, edge cases
- `design-decisions.md` — why important decisions were made

When commaBot learns something durable, it should update one of these files.

---

### `PRD.md` and `VISION.md`

Product context.

- `PRD.md` defines the v1.0 product requirement and scope
- `VISION.md` describes the aspirational long-term package

These are not task boards.

---

### `ONBOARDING.md`

Claire's learning path.

This is written for a second-year PhD student who knows R/dplyr/ggplot2 but is new to S4, roxygen2, and testthat.

---

### `archive/`

Historical documents that are no longer active but should not be deleted.

Examples:
- One-time agent bootstrap instructions
- Resolved implementation specs

---

## Project Management Rules

### 1. One backlog

Every proposed work item goes in `BACKLOG.md`. No scattered TODOs in random notes.

### 2. One status board

`STATUS.md` shows what is happening now. Keep it short.

### 3. Knowledge is topical, not chronological

Do not create endless dated audit files. If an audit produces durable knowledge, absorb it into `knowledge/`.

### 4. Separate strategy from tactics

- Strategy: `ROADMAP.md`, `VISION.md`
- Tactics: `BACKLOG.md`, `STATUS.md`
- Evidence/knowledge: `knowledge/`

### 5. Archive, don't delete, when context may matter later

Move dead-but-useful documents to `archive/`.

---

## For commaBot

When doing project-management work:

1. Document findings in the appropriate `knowledge/` file
2. Add or update work items in `BACKLOG.md`
3. Update `STATUS.md` if work starts, finishes, or gets blocked
4. Update `ROADMAP.md` only for strategic changes
5. Do not silently fix without documenting why

Carl is the product owner. commaBot is the engineering team and PM.
