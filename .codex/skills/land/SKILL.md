---
name: land
description:
  Land a PR by monitoring conflicts, resolving them, waiting for checks, and
  squash-merging when green; use when asked to land, merge, or shepherd a PR to
  completion.
---

# Land

## Goals

- Ensure the PR is conflict-free with main.
- Keep CI green and fix failures when they occur.
- Squash-merge the PR once checks pass.
- Do not yield to the user until the PR is merged; keep the watcher loop running
  unless blocked.
- Verify that GitHub removed the same-repository head branch under the enabled
  automatic cleanup setting. Fork branches remain owned by their source repos.

## Preconditions

- `gh` CLI is authenticated.
- The authenticated account has maintainer merge authority.
- You are on the PR branch with a clean working tree.

## Steps

1. Locate the PR for the current branch.
2. Confirm the scope-appropriate validation required by `AGENTS.md` is green
   locally before any push.
3. If the working tree has uncommitted changes, commit with the `commit` skill
   and push with the `push` skill before proceeding.
4. Check mergeability and conflicts against main.
5. If conflicts exist, use the `pull` skill to fetch/merge `origin/main` and
   resolve conflicts, then use the `push` skill to publish the updated branch.
6. Ensure Codex review comments (if present) are acknowledged and any required
   fixes are handled before merging.
7. Watch checks until complete.
8. If checks fail, pull logs, fix the issue, commit with the `commit` skill,
   push with the `push` skill, and re-run checks.
9. When all checks are green and review feedback is addressed, squash-merge
   using the PR title/body for the merge subject/body, then verify the merged
   state and configured same-repository branch cleanup.
10. **Context guard:** Before implementing review feedback, confirm it does not
    conflict with the user’s stated intent or task context. If it conflicts,
    respond inline with a justification and ask the user before changing code.
11. **Pushback template:** When disagreeing, reply inline with: acknowledge +
    rationale + offer alternative.
12. **Ambiguity gate:** When ambiguity blocks progress, use the clarification
    flow (assign PR to current GH user, mention them, wait for response). Do not
    implement until ambiguity is resolved.
    - If you are confident you know better than the reviewer, you may proceed
      without asking the user, but reply inline with your rationale.
13. **Per-comment mode:** For each review comment, choose one of: accept,
    clarify, or push back. Reply inline (or in the issue thread for Codex
    reviews) stating the mode before changing code.
14. **Reply before change:** Always respond with intended action before pushing
    code changes (inline for review comments, issue thread for Codex reviews).

## Commands

```
(
set -euo pipefail

# Ensure branch and PR context
pr_title=$(gh pr view --json title -q .title)
pr_body=$(gh pr view --json body -q .body)
head_ref=$(gh pr view --json headRefName -q .headRefName)
head_sha=$(gh pr view --json headRefOid -q .headRefOid)
is_cross_repo=$(gh pr view --json isCrossRepository -q .isCrossRepository)

# Check mergeability and conflicts
mergeable=$(gh pr view --json mergeable -q .mergeable)

if [ "$mergeable" = "CONFLICTING" ]; then
  echo "Run the pull skill, resolve conflicts, then publish with push." >&2
  exit 5
fi

# Run the single authoritative review/check gate. Stop if it is unavailable.
python3 .codex/skills/land/land_watch.py || exit $?

# Use the ordinary protected-branch path after the watcher proves readiness.
gh pr merge --squash --match-head-commit "$head_sha" \
  --subject "$pr_title" --body "$pr_body" || exit $?
test "$(gh pr view --json state -q .state)" = "MERGED" || exit 6

# GitHub normally removes same-repository heads under this repository's enabled
# automatic cleanup setting. Verify that outcome and use an explicit fallback.
# Fork branches belong to their source repositories and must not be deleted here.
if [ "$is_cross_repo" = "false" ]; then
  if git ls-remote --exit-code --heads origin "refs/heads/$head_ref" >/dev/null 2>&1; then
    git push origin --delete "$head_ref" || {
      echo "PR merged, but remote branch cleanup requires manual follow-up." >&2
      exit 7
    }
  else
    branch_check_status=$?
    if [ "$branch_check_status" -ne 2 ]; then
      echo "PR merged, but remote branch cleanup could not be verified." >&2
      exit "$branch_check_status"
    fi
  fi
fi
)
```

## Landing Watcher

Use the tested asyncio watcher as the only supported review/check gate. If
Python or the helper is unavailable, stop and report that environment blocker;
do not substitute a reduced shell gate.

```
python3 .codex/skills/land/land_watch.py
```

Exit codes:

- 2: Review comments detected (address feedback)
- 3: CI checks failed
- 4: PR head updated (inspect and synchronize the new head)
- 5: PR has merge conflicts (merge `origin/main`, resolve, and push)

## Failure Handling

- If checks fail, pull details with `gh pr checks` and `gh run view --log`, then
  fix locally, commit with the `commit` skill, push with the `push` skill, and
  re-run the watch.
- If the watcher exits 4, fetch and inspect the new PR head. Fast-forward the
  clean local feature branch from its remote counterpart when possible; use
  the `pull` skill for divergent history or conflicts. Re-run local validation
  and restart the watcher against the synchronized head. Stop on an unexpected
  or unauthenticated update.
- If mergeability is `UNKNOWN`, wait and re-check.
- Do not merge while review comments (human or Codex review) are outstanding.
- Codex review jobs retry on failure and are non-blocking. Use a submitted Codex
  review for the current PR head as the completion signal, not job status;
  treat its issue and inline comments as blocking feedback.
- Do not enable auto-merge unless the user requests it and the repository's
  branch-protection behavior is understood.
- Verify that GitHub signed the resulting squash commit. If an ordinary merge
  fails solely because GitHub cannot satisfy the required-signature rule, first
  reconfirm the current head, watcher result, and live protection settings; use
  `gh pr merge --admin` only with explicit repository-administrator authority.
- Verify same-repository head-branch deletion after merge. If automatic cleanup
  fails, remove the merged remote branch explicitly; never delete a fork branch.

## Review Handling

- Codex review results can include a submitted review plus inline comments, or
  `## Codex Review — <persona>` issue comments posted by GitHub Actions. Treat
  every finding as feedback that must be acknowledged before merge.
- Human review comments are blocking and must be addressed (responded to and
  resolved) before requesting a new review or merging.
- If multiple reviewers comment in the same thread, respond to each comment
  (batching is fine) before closing the thread.
- Fetch review comments via `gh api` and reply with a prefixed comment.
- Use review comment endpoints (not issue comments) to find inline feedback:
  - List PR review comments:
    ```
    gh api repos/{owner}/{repo}/pulls/<pr_number>/comments
    ```
  - PR issue comments (top-level discussion):
    ```
    gh api repos/{owner}/{repo}/issues/<pr_number>/comments
    ```
  - Reply to a specific review comment:
    ```
    gh api -X POST /repos/{owner}/{repo}/pulls/<pr_number>/comments \
      -f body='[codex] <response>' -F in_reply_to=<comment_id>
    ```
- `in_reply_to` must be the numeric review comment id (e.g., `2710521800`), not
  the GraphQL node id (e.g., `PRRC_...`), and the endpoint must include the PR
  number (`/pulls/<pr_number>/comments`).
- If GraphQL review reply mutation is forbidden, use REST.
- A 404 on reply typically means the wrong endpoint (missing PR number) or
  insufficient scope; verify by listing comments first.
- All GitHub comments generated by this agent must be prefixed with `[codex]`.
- For a Codex review issue comment, acknowledge that exact result with a root
  issue comment of the form `[codex] Review <comment-id>: <disposition>`. State
  whether you will address the feedback now or defer it (include rationale).
  The result must include a `Reviewed commit` SHA marker for the current head.
  The watcher accepts only results newer than the current head's first
  check run and any later explicit review request. If the result body is edited
  after acknowledgement, inspect the changed content and post a fresh exact-ID
  acknowledgement.
- For a substantive top-level automated review body, acknowledge its findings
  with a root-level `[codex] Review <review-id>: <disposition>` issue comment.
  Generated Copilot “Pull request overview” summaries are informational and do
  not require acknowledgement.
- If feedback requires changes:
  - For inline review comments (human), reply with intended fixes
    (`[codex] ...`) **as an inline reply to the original review comment** using
    the review comment endpoint and `in_reply_to` (do not use issue comments for
    this).
  - Implement fixes, commit, push.
  - Reply with the fix details and commit sha (`[codex] ...`) in the same place
    you acknowledged the feedback (issue comment for Codex reviews, inline reply
    for review comments).
  - The land watcher treats Codex review issue comments as unresolved until a
    newer `[codex]` issue comment is posted acknowledging the findings.
- Only request a new Codex review when you need a rerun (e.g., after new
  commits). Do not request one without changes since the last review.
  - Before requesting a new Codex review, re-run the land watcher and ensure
    there are zero outstanding review comments (all have `[codex]` inline
    replies).
  - After pushing new commits, the Codex review workflow will rerun on PR
    synchronization (or you can re-run the workflow manually). Post a concise
    root-level summary comment so reviewers have the latest delta:
    ```
    [codex] Changes since last review:
    - <short bullets of deltas>
    Commits: <sha>, <sha>
    Tests: <commands run>
    ```
  - Only request a new review if there is at least one new commit since the
    previous request.
  - Wait for the next Codex review result for the current head before merging.

## Scope + PR Metadata

- The PR title and description should reflect the full scope of the change, not
  just the most recent fix.
- If review feedback expands scope, decide whether to include it now or defer
  it. You can accept, defer, or decline feedback. If deferring or declining,
  call it out in the root-level `[codex]` update with a brief reason (e.g.,
  out-of-scope, conflicts with intent, unnecessary).
- Correctness issues raised in review comments should be addressed. If you plan
  to defer or decline a correctness concern, validate first and explain why the
  concern does not apply.
- Classify each review comment as one of: correctness, design, style,
  clarification, scope.
- For correctness feedback, provide concrete validation (test, log, or
  reasoning) before closing it.
- When accepting feedback, include a one-line rationale in the root-level
  update.
- When declining feedback, offer a brief alternative or follow-up trigger.
- Prefer a single consolidated "review addressed" root-level comment after a
  batch of fixes instead of many small updates.
- For doc feedback, confirm the doc change matches behavior (no doc-only edits
  to appease review).
