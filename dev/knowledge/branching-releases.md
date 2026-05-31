# Branches, Pull Requests, and Releases for commaKit

**Last updated:** 2026-05-14  
**Maintained by:** commaBot  
**Audience:** Carl, Claire, commaBot, and future contributors

This document defines how commaKit uses branches, pull requests, versions, releases, and tags.

The goal is not office-job ceremony. The goal is to prevent avoidable chaos.

---

## Mental Model

Branches, pull requests, and releases are **risk-management tools**.

They answer different questions:

| Tool | Question it answers |
|------|---------------------|
| Branch | Where can work happen safely before it is ready? |
| Pull request | What changed, why, and how was it reviewed/tested? |
| Release | Which exact commit is stable enough to name and cite? |
| Tag | Which exact commit corresponds to a release version? |

---

## Branch Policy

### `main`

`main` is the current known-good development state.

It does **not** mean perfect. It does mean:

- no half-finished refactors
- no mystery piles of uncommitted work
- no knowingly broken examples or tests unless documented
- usable by Carl or Claire if they clone the repo

### Feature branches

Nontrivial work happens on a branch.

Use branches for:

- exported function changes
- S4 class changes
- statistical behavior changes
- parser changes
- plot semantics changes
- test additions across multiple files
- vignette rewrites
- release preparation
- anything where rollback might matter

Branch name format:

```text
type/short-description
```

Types:

```text
test/
fix/
docs/
api/
feature/
refactor/
release/
pm/
```

Examples:

```text
test/full-pipeline-integration
docs/data-import-troubleshooting
docs/diffmethyl-method-selection
fix/sliding-window-circular-boundary
api/mod-type-vector-consistency
release/0.1.0
```

### Direct-to-main exceptions

Small, low-risk changes may go directly to `main`:

- typo fixes
- small `.gitignore` / `.Rbuildignore` changes
- small planning/documentation updates
- GitHub Issues updates
- focused project-management docs

When in doubt, branch.

---

## Pull Request Policy

A PR is a review artifact, not just a merge button.

Every meaningful code/test/API/docs change should have a PR, even if Carl is the only human reviewer. The PR creates:

- a checkpoint
- CI history
- a readable diff
- a place to discuss decisions
- a clean historical record

### Required PR template

```markdown
## Summary
- What changed
- Why it changed

## Test plan
- [ ] devtools::test()
- [ ] devtools::document() if roxygen changed
- [ ] devtools::check() if release-like or high-risk

## Risk
Low / Medium / High

What could break?

## Notes for Carl
Scientific/domain assumptions, decisions needed, or anything that needs product-owner judgment.
```

### What Carl should review

Carl should not need to inspect every line. commaBot should summarize:

- what changed
- why it matters
- what was tested
- what assumptions need Carl's domain judgment

Carl's review role is product/science judgment, not clerical diff inspection.

---

## Release Policy

A release is a named stability point. It is different from a commit.

A release requires:

- clean git status
- `devtools::test()` passes
- `devtools::check()` passes, or any deviations are explicitly accepted
- `NEWS.md` updated
- `DESCRIPTION` version correct
- README/vignettes consistent with the version
- release tag pushed

### Release branches

Use a release branch when preparing an actual named release:

```bash
git switch -c release/0.1.0
```

Allowed on release branches:

- version bump
- NEWS updates
- documentation fixes
- test/check fixes
- R CMD check / BiocCheck fixes

Not allowed on release branches:

- new features
- API redesigns
- exploratory statistical changes
- unrelated refactors

Release flow:

```text
release/x.y.z → PR → main → tag vx.y.z → GitHub release
```

---

## Tag Policy

Tags are only for real releases.

Example:

```bash
git tag v0.1.0
git push origin v0.1.0
```

Do not tag every merge. Do not tag ordinary development commits.

A tag means:

> This exact commit is version X.

---

## Version Policy

commaKit's canonical current version is `0.2.0`. Earlier version numbers were
informal labels from the script-to-package transition and should not drive
current metadata, badges, URLs, or release planning.

### Version meanings

| Version pattern | Meaning |
|-----------------|---------|
| `0.2.0` | current Schema v2 internal release baseline |
| `0.2.x` | patch releases on the current internal baseline |
| `0.3.0` | next deliberate internal feature milestone |
| `0.x.y` | internal/pre-Bioconductor releases |
| `0.99.0` | Bioconductor submission version |
| `1.0.0` | only after external confidence and real stability |

### Development suffix

The `.9000` suffix means active development. This is common R package practice.

Example:

```text
0.2.0      released checkpoint
0.2.0.9000 development after that checkpoint
0.2.1      patch release
0.2.1.9000 development after patch release
```

---

## Who Owns What

### commaBot owns

- deciding when a branch is needed
- creating branches and PRs
- maintaining clean commit history
- writing PR summaries
- running tests/checks or reporting when they were not run
- updating GitHub Issues and `NEWS.md`
- preventing accidental commits of local state or secrets

### Carl owns

- scientific correctness
- product priorities
- deciding whether a behavior is biologically appropriate
- deciding when something is stable enough to release
- Bioconductor submission timing

### Claire owns, when she starts contributing

- using branches for nontrivial work
- asking when unsure
- running tests before commits
- keeping changes focused

---

## Practical Workflows

### Small direct-to-main docs/config fix

```bash
git status --short
# edit one or two files
git diff
git add specific-file
git commit -m "Clear message"
git push
```

### Normal feature/fix workflow

```bash
git status --short
git switch -c test/full-pipeline-integration
# make changes
Rscript -e "devtools::test()"
git status --short
git diff
git add specific-files
git commit -m "Clear message"
git push -u origin test/full-pipeline-integration
gh pr create
```

### Release workflow

```bash
git switch main
git pull
git switch -c release/0.1.0
# update DESCRIPTION, NEWS, docs
Rscript -e "devtools::test()"
Rscript -e "devtools::check()"
git commit
git push -u origin release/0.1.0
gh pr create
# after merge
git tag v0.1.0
git push origin v0.1.0
```

---

## Default Rule

If commaBot is unsure whether to branch, commaBot should branch.

If Carl is unsure whether a change is release-worthy, do not release.

If Claire is unsure whether files belong in the same commit, ask commaBot.
