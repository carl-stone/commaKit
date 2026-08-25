---
type: Documentation Standard
title: commaKit OKF Guidelines
description: Local conventions for maintaining the dev/knowledge Open Knowledge Format bundle.
resource: dev/knowledge
tags: [okf, documentation, progressive-disclosure]
timestamp: 2026-06-15T00:00:00Z
status: current
owner: Carl Stone
---

# Summary

`dev/knowledge/` is the project's durable context library. It is intentionally plain: markdown files with YAML frontmatter, ordinary links, and small index files for navigation.

# Local Conformance Rules

Every non-reserved `.md` file in this directory is a concept document and must start with YAML frontmatter containing at least:

```yaml
---
type: Short Concept Type
title: Human-readable title
description: One-sentence summary for search and previews.
tags: [short, useful, tags]
timestamp: 2026-06-15T00:00:00Z
---
```

Use these optional fields when useful:

```yaml
resource: path/or/url/described/by/this/card
status: current
owner: Carl Stone
```

Reserved files:

- `index.md`: directory listing for progressive disclosure.
- `log.md`: chronological update history.

# Writing Rules

- Put one concept per file.
- Use descriptive `type` values; consumers must tolerate unknown types, so clarity beats taxonomy cleverness.
- Prefer absolute repo links or bundle-relative links that remain stable for readers.
- Use headings, tables, and fenced code blocks where structure helps retrieval.
- Include `# Citations` for externally sourced claims.
- Do not store issue scratch notes here unless the lesson is durable.
- Mark historical audits as historical instead of pretending they are live inventories.

# Search Patterns

```bash
rg "^type:" dev/knowledge
rg "mod_context|siteCoverage|result layer" dev/knowledge
rg "^# Citations" dev/knowledge
```

# Citations

[1] [Introducing the Open Knowledge Format](https://cloud.google.com/blog/products/data-analytics/how-the-open-knowledge-format-can-improve-data-sharing)

[2] [Open Knowledge Format v0.1 draft spec](https://github.com/GoogleCloudPlatform/knowledge-catalog/blob/main/okf/SPEC.md)
