# comma Status Board

**Last updated:** 2026-05-14  
**Maintained by:** commaBot

This is the daily project view: what is happening now, what is next, what is blocked, and what was recently completed.

For the full prioritized work list, see `dev/BACKLOG.md`.

---

## Current Sprint

**Goal:** Put the project under real project management before doing more code work.

### In Progress

| ID | Task | Owner | Notes |
|----|------|-------|-------|
| T-001 | Add integration test across full pipeline | commaBot | Implemented on branch; awaiting PR review/merge |

---

## Next Up

These are the highest-priority proposed tasks from `BACKLOG.md`. They are not started until Carl accepts them.

| ID | Task | Why it matters | Size |
|----|------|----------------|------|
| P-001 | Add troubleshooting guide for data import | Helps users through the most likely failure point | M |
| P-002 | Unify `mod_type` parameter type across functions | Reduces UX friction and confusion | M |
| P-003 | Add method selection guidance for `diffMethyl()` | Helps users choose the right statistical backend | S |
| T-002 | Add data verification to plot tests | Prevents silent plot mapping bugs | L |

---

## Blocked / Waiting

| ID | Task | Blocked on |
|----|------|------------|
| B-004 | Submit to Bioconductor | Carl registering on support.bioconductor.org; final decision to submit |
| B-001 | Register Zenodo DOI | Carl decision on when package is submission-ready |

---

## Recently Completed

| Date | Task | Notes |
|------|------|-------|
| 2026-05-14 | Reorganized `dev/` project management | Added README, BACKLOG, knowledge base, kanban STATUS, strategic ROADMAP |
| 2026-05-14 | Documentation sync audit | Docs are in sync with code; examples tested and pass |
| 2026-05-14 | Package health audit | UX gaps identified: troubleshooting, method guidance, real-world example |
| 2026-05-14 | Test quality audit | Strong core tests; weak plot tests; no full integration test |
| 2026-05-14 | R CMD check clean | 0 errors, 0 warnings, 0 notes |
| 2026-05-14 | Version reset to 0.1.0.9000 | Fresh-start baseline for serious package development |

---

## Health Snapshot

| Area | Status | Notes |
|------|--------|-------|
| Core features | ✅ Complete | v1.0 feature set is implemented |
| R CMD check | ✅ Clean | 0 errors, 0 warnings, 0 notes |
| Tests | 🟡 Mixed | 938 passing; core stats strong; plots weak; integration missing |
| Documentation sync | ✅ Good | Function docs match code |
| User guidance | 🟡 Needs work | Troubleshooting and method-selection guidance missing |
| Bioconductor readiness | 🟡 Not submitting yet | 0.99.0 reserved for future submission; current baseline is 0.1.0.9000 |

---

## How to Use This File

- **Carl:** Open this file to see what is happening now.
- **commaBot:** Keep this file updated whenever work starts, finishes, or gets blocked.
- **Backlog lives elsewhere:** Do not turn this into a giant task list. Use `dev/BACKLOG.md` for that.
