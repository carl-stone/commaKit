# Domain Docs

How the engineering skills should consume this repo's domain documentation when exploring the codebase.

## Layout

commaKit uses the single-context domain-doc layout:

- `CONTEXT.md` at the repo root for domain vocabulary, if it exists.
- `docs/adr/` at the repo root for architectural decisions, if it exists.

## Before exploring, read these

- `CONTEXT.md` at the repo root, if it exists.
- ADRs under `docs/adr/` that touch the area you are about to work in, if they exist.
- The existing commaKit knowledge bundle under `dev/knowledge/`, especially `dev/knowledge/index.md`, when the task concerns package architecture, project status, tests, durable decisions, or agent operating model.

## File structure

```text
/
|-- CONTEXT.md
|-- docs/
|   `-- adr/
|       |-- 0001-example-decision.md
|       `-- 0002-example-decision.md
|-- dev/
|   `-- knowledge/
`-- R/
```

## Use the glossary's vocabulary

When your output names a domain concept in an issue title, refactor proposal, hypothesis, or test name, use the term as defined in `CONTEXT.md` when that file exists. Do not drift to synonyms the glossary explicitly avoids.

If the concept you need is not in the glossary yet, either reconsider whether you are inventing language the project does not use, or note the gap for `/domain-modeling`.

## Flag ADR conflicts

If your output contradicts an existing ADR, surface it explicitly rather than silently overriding.
