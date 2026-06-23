# Update durable knowledge with code changes

When commaKit code or behavior changes make durable project knowledge stale, the same branch should update the relevant docs. This keeps future agents and maintainers from being routed by outdated `dev/knowledge/`, `CONTEXT.md`, or ADR content.

## Consequences

PRs that close issues by changing behavior should update durable knowledge when the issue changes project state, test coverage, known gaps, or parser/user-facing contracts. Splitting the docs into a later PR is acceptable only when the code PR explicitly does not claim to settle the stale knowledge.
