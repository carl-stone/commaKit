# commaKit Development Directory

This directory holds project-management context, durable knowledge, and strategic documents for commaKit.

**Audience:** Carl, Claire, future agents, and any human contributor

---

## Tactical Work: GitHub Issues & Pull Requests

**[GitHub Issues](https://github.com/carl-stone/commaKit/issues) track accepted work; [GitHub Pull Requests](https://github.com/carl-stone/commaKit/pulls) are the code-integration surface.**

Use issues for work state and PRs for code review, CI, and merge state.

Labels follow a namespaced scheme:

| Namespace | Labels | Purpose |
|---|---|---|
| `type:` | `bug`, `docs`, `test`, `cleanup`, `api`, `data`, `audit`, `admin` | What kind of work |
| `area:` | `import`, `diffMethyl`, `plots`, `slidingWindow`, `enrichment`, `bioconductor`, `pm` | Which package area |
| `priority:` | `high`, `medium`, `low` | Urgency |
| `status:` | `blocked`, `needs-decision`, `accepted` | GitHub/public workflow state |

## Local R Version

The project records R 4.5.3 in `renv.lock`, `renv/settings.json`, and
`.R-version`. Use `.R-version` with tools such as `rig`, `mise`, or `asdf` to
select the same R minor version as CI before activating `renv`.

`renv.lock` remains the dependency source of truth. `.R-version` only helps
local version managers choose the matching R executable.

On macOS machines where the default framework R has moved ahead of CI,
`dev/run-r-4.5.sh` is a convenience wrapper for running checks against a local
R 4.5 framework.

## VS Code Devcontainer

The preferred low-drama VS Code setup is the repo devcontainer in
`.devcontainer/`. It runs the package inside a Linux container based on the
Bioconductor 3.22 image, which provides the matching R 4.5 runtime independent
of the host machine's R framework. The container runs as the non-root `rstudio`
user under Docker's standard isolation; it does not require privileged mode.

From VS Code, use **Dev Containers: Reopen in Container**. On first create, the
container runs:

```bash
Rscript --vanilla .devcontainer/setup.R
```

That restores `renv.lock` and installs editor helper packages such as
`languageserver`, `httpgd`, and `lintr` into the container environment. The
container image also installs the Air formatter. After the container is ready,
the normal development commands are:

```bash
Rscript -e "devtools::load_all()"
Rscript -e "devtools::test()"
Rscript -e "devtools::check(build_args = c('--no-build-vignettes'))"
```

---

## Strategic & Context Documents

These files are maintained in `dev/`:

| Question | Read this |
|----------|-----------|
| Where is the package going? | `ROADMAP.md` |
| What is v1.0 supposed to be? | `PRD.md` |
| What is the dream version? | `VISION.md` |
| How should Claire learn the package? | `ONBOARDING.md` |

---

## Durable Knowledge (`knowledge/`)

Organized by topic, not by date. Durable project knowledge goes here.

| File | What it covers |
|------|---------------|
| `test-quality.md` | What tests are strong, weak, or missing |
| `known-issues.md` | Bugs, gotchas, edge cases |
| `design-decisions.md` | Why important decisions were made |
| `git-discipline.md` | How commits, staging, generated files, and local state are handled |
| `branching-releases.md` | Branch, PR, release, tag, and version policy |

---

## Archive (`archive/`)

Historical documents that are no longer active but should not be deleted. These are preserved for reference only — do not update them.

| File | Why it was archived |
|------|-------------------|
| `BACKLOG.md` | Pre-GitHub-Issues task tracking. Migrated to GitHub Issues 2026-05-15. |
| `STATUS.md` | Pre-GitHub-Issues sprint board. Migrated to GitHub Issues 2026-05-15. |
| `AGENT_BOOTSTRAP.md` | One-time agent setup instructions, no longer needed. |
| `SPECS.md` | Resolved implementation spec. |

---

## For Contributors

When doing project-management work:

1. Use GitHub Issues to create, label, decompose, and update accepted work.
2. Use GitHub PRs for code integration evidence; do not treat a worker response as Done before PR/CI/review state is resolved or explicitly blocked.
3. Document durable findings in the appropriate `knowledge/` file.
4. Update `ROADMAP.md` only for strategic changes.
5. Do not silently fix without documenting why.
6. Do not create new files in `dev/` for task tracking — use GitHub Issues.

Carl is the product owner.
