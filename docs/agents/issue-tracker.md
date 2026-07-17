# Issue Tracking Surfaces

Internal commaBot work and PRDs live in Linear project `commaKit Symphony`.
Public bugs and external-contributor requests live in GitHub Issues. Code lands
through GitHub pull requests regardless of where the work originated.

## Internal Symphony work

Use the repository `linear` skill and Symphony's `linear_graphql` tool for
Linear reads and mutations. An `STO-...` key identifies a Linear issue. Keep
worker state, decomposition, and Carl-facing task decisions in Linear; attach
the resulting GitHub PR to that issue.

## GitHub issues

Use the `gh` CLI for public bug reports, accepted external requests, and
external-contributor triage.

## Conventions

- **Create an issue**: `gh issue create --title "..." --body "..."`. Use a heredoc for multi-line bodies.
- **Read an issue**: `gh issue view <number> --comments`, filtering comments by `jq` and also fetching labels.
- **List issues**: `gh issue list --state open --json number,title,body,labels,comments --jq '[.[] | {number, title, body, labels: [.labels[].name], comments: [.comments[].body]}]'` with appropriate `--label` and `--state` filters.
- **Comment on an issue**: `gh issue comment <number> --body "..."`
- **Apply or remove labels**: `gh issue edit <number> --add-label "..."` / `--remove-label "..."`
- **Close**: `gh issue close <number> --comment "..."`

Infer the repo from `git remote -v`; `gh` does this automatically when run inside a clone.

## Pull requests as a triage surface

**PRs as a request surface: yes.**

External PRs run through the same labels and states as issues, using the `gh pr` equivalents:

- **Read a PR**: `gh pr view <number> --comments` and `gh pr diff <number>` for the diff.
- **List external PRs for triage**: `gh pr list --state open --json number,title,body,labels,author,authorAssociation,comments` then keep only `authorAssociation` of `CONTRIBUTOR`, `FIRST_TIME_CONTRIBUTOR`, or `NONE` and drop `OWNER`, `MEMBER`, and `COLLABORATOR`.
- **Comment, label, or close**: `gh pr comment`, `gh pr edit --add-label` / `--remove-label`, and `gh pr close`.

GitHub shares one number space across issues and PRs, so a bare `#42` may be either. Resolve with `gh pr view 42` and fall back to `gh issue view 42`.

## Resolving generic tracker instructions

- For internal commaBot/Symphony work, publish to Linear project
  `commaKit Symphony`.
- For public bugs or external-contributor requests, create a GitHub issue in
  `carl-stone/commaKit`.
- An `STO-...` key is Linear; a GitHub `#<number>` is an issue or PR.

## When a skill says "fetch the relevant ticket"

Use the `linear` skill for an `STO-...` key. For a GitHub number, run
`gh issue view <number> --comments`; if it is a PR, use
`gh pr view <number> --comments` instead.
