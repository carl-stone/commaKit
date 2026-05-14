# Git Discipline for comma

**Last updated:** 2026-05-14  
**Maintained by:** commaBot  
**Audience:** Carl, Claire, commaBot, and future contributors

This document defines how we use git on comma from this point forward.

The goal is not bureaucracy. The goal is **trust**: every commit should tell a coherent story, be reviewable, and leave the repository in a working state.

---

## Core Rule

> **Every commit should represent one logical change, and the repository should work after every commit.**

A good commit answers:

1. What changed?
2. Why did it change?
3. What files belong together?
4. Can someone review this without guessing?

If the answer is no, the commit is too big, too mixed, or poorly described.

---

## What Perfect Git Discipline Means Here

### 1. Start from a clean understanding of the repo

Before making changes:

```bash
git status --short
git log --oneline -5
```

Know:
- What branch you are on
- Whether there are existing uncommitted changes
- Whether those changes are yours, Carl's, or another agent's

If the tree is dirty, do **not** blindly edit. Understand what is already there.

---

### 2. One logical change per commit

Good commit boundaries:

- One bug fix
- One documentation update
- One test reorganization
- One feature
- One project-management change

Bad commit boundaries:

- "stuff"
- "updates"
- README changes plus statistical model changes plus unrelated test renames
- Formatting mixed with behavior changes

Examples from this repo:

Good:

```text
Split coverage analysis tests by exported function.
Refresh package metadata and check-safe docs.
Update README and getting-started workflow docs.
Organize development project management docs.
```

Bad:

```text
fix things
updates
more changes
wip
```

---

### 3. Stage intentionally

Prefer staging specific files:

```bash
git add R/diffMethyl.R tests/testthat/test-diffMethyl.R
```

Avoid casual broad staging:

```bash
git add .
git add -A
```

Broad staging is only acceptable after carefully reviewing `git status` and confirming every file belongs in the commit.

Never commit:
- `.env`
- credentials
- API keys
- local scratch directories
- huge generated files unless intentionally part of the package
- agent state directories (`.codex/`, `.letta/`, `.lteams/`)

---

### 4. Review before committing

Before every commit:

```bash
git status --short
git diff --staged
```

Ask:
- Do all staged files belong together?
- Are there accidental generated files?
- Are there secrets?
- Does the diff match the commit message?

For generated files (`man/*.Rd`, `README.md`, README figures), check that they correspond to the source change (`R/*.R`, `README.Rmd`).

---

### 5. Commit messages should explain the why

Format:

```text
Short imperative summary.

One or two sentences explaining why this change exists and what problem it solves.
```

Examples:

```text
Add integration test for core methylation workflow.

Exercise annotateSites, diffMethyl, results, and filterResults together so interface drift between pipeline stages is caught by the test suite.
```

```text
Document statistical backend selection.

Give users practical guidance for choosing methylKit, limma, or quasi_f so the default is not treated as the only valid analysis path.
```

Avoid messages that only describe mechanics:

```text
Changed test file
Edited README
Fixed typo
```

Those may be acceptable for tiny changes, but most commits need the reason.

---

### 6. Test before committing code changes

For R code changes:

```bash
Rscript -e "devtools::test()"
```

For documentation changes that affect roxygen output:

```bash
Rscript -e "devtools::document()"
```

For README changes:

```bash
Rscript -e "devtools::build_readme()"
```

Before PRs or release-like commits:

```bash
Rscript -e "devtools::check()"
```

If tests are not run, say so explicitly in the final summary.

---

### 7. Generated files move with their sources

If editing roxygen comments in `R/*.R`, commit the generated `man/*.Rd` changes in the same commit.

If editing `README.Rmd`, commit the generated `README.md` and `man/figures/*` changes in the same commit.

Do not commit generated files without the source that generated them unless the task is explicitly regeneration-only.

---

### 8. Do not mix formatting with behavior

If reformatting many files, commit that separately from logic changes.

Why:
- Formatting diffs hide real code changes
- Review becomes harder
- Bugs are harder to trace with `git blame`

Good sequence:

1. Commit behavior change
2. Commit formatting cleanup

Bad sequence:

1. Reformat 20 files and change statistical logic in the middle

---

### 9. Keep local/agent state out of the repo

These directories are local tooling state and should remain untracked:

```text
.codex/
.letta/
.lteams/
```

They may be useful locally, but they are not package source. If an agent needs durable knowledge, it should go in:

- commaBot memory, or
- `dev/knowledge/` if it is project knowledge useful to humans

---

### 10. Push only after coherent commits exist

Pushing is publication. Before pushing:

```bash
git log --oneline -5
git status --short
```

Check:
- Commit sequence makes sense
- No accidental files remain staged
- Working tree is either clean or only contains intentionally untracked local state

---

## Special Rules for commaBot

commaBot must follow stricter discipline than a human because it can change many files quickly.

### Before editing

- Inspect current status
- Identify pre-existing changes
- Do not overwrite Carl's work without understanding it

### Before committing

- Review status and staged diff
- Batch files into logical commits
- Never commit secrets or local state
- Never amend unless Carl explicitly asks
- Never force push unless Carl explicitly asks, and never force push `main`

### After committing

Report:
- Commit hashes
- What each commit contains
- Whether tests/docs/checks were run
- What remains untracked or modified

---

## Special Rules for Claire

Claire should use a simple, safe workflow:

```bash
git status --short
# make one focused change
Rscript -e "devtools::test()"
git status --short
git diff
git add specific-file.R tests/testthat/test-specific-file.R
git commit -m "Clear message explaining the change"
```

If unsure whether files belong together, ask commaBot before committing.

---

## How to Decide Commit Boundaries

Ask: "Would I want to revert this independently?"

If yes, separate commit.

Examples:

| Change | Same commit? | Why |
|--------|--------------|-----|
| `R/diffMethyl.R` + `test-diffMethyl.R` | Yes | Code and tests for same behavior |
| `R/diffMethyl.R` + `man/diffMethyl.Rd` | Yes | Roxygen source and generated docs |
| `README.Rmd` + `README.md` + README figures | Yes | Source and generated README outputs |
| `plot_volcano()` bug fix + unrelated README rewrite | No | Different logical changes |
| Test file renames + statistical model change | No | Refactor vs behavior change |
| `.Rbuildignore` update + local agent directories | No | Ignore pattern may be committed; local directories should not |

---

## Current Repo Policy

- `main` is allowed for now because Carl is the only human user and commaBot is the only stateful agent.
- If Claire starts contributing regularly, switch to feature branches and PRs.
- Bioconductor submission work should happen on a branch when the time comes.

---

## The Standard

A clean git history should let future Carl, Claire, or commaBot understand the project without archaeology.

Each commit should feel like a sentence in the story of the package.
