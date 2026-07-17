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
- No need to delete remote branches after merge; the repo auto-deletes head
  branches.

## Preconditions

- `gh` CLI is authenticated.
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
   using the PR title/body for the merge subject/body. Let the repository's
   configured post-merge behavior delete the remote head branch.
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
# Ensure branch and PR context
branch=$(git branch --show-current)
pr_number=$(gh pr view --json number -q .number)
pr_title=$(gh pr view --json title -q .title)
pr_body=$(gh pr view --json body -q .body)

# Check mergeability and conflicts
mergeable=$(gh pr view --json mergeable -q .mergeable)

if [ "$mergeable" = "CONFLICTING" ]; then
  # Run the `pull` skill to handle fetch + merge + conflict resolution.
  # Then run the `push` skill to publish the updated branch.
fi

# Preferred: use the Async Watch Helper below. The manual loop is a fallback
# when Python cannot run or the helper script is unavailable.
# Wait for either a submitted Codex review or an exactly acknowledged Codex
# issue-review result for the current PR head. Review comments and issue
# comments remain blocking feedback and must be handled before merging.
head_sha=$(gh pr view --json headRefOid -q .headRefOid)
head_marker=${head_sha:0:10}
acknowledger_login=$(gh api user --jq .login)
latest_review_request_at=$(
  gh api repos/{owner}/{repo}/issues/"$pr_number"/comments --paginate \
    --jq '.[] | select(.user.login != "chatgpt-codex-connector[bot]" and ((.body // "") | contains("@codex review"))) | .created_at' \
    | sort | tail -n 1
)
head_check_at=$(gh api repos/{owner}/{repo}/commits/"$head_sha"/check-runs \
  --paginate --jq '.check_runs[].created_at' | sort | head -n 1)
issue_review_not_before=$(
  printf '%s\n%s\n' "$head_check_at" "$latest_review_request_at" \
    | sed '/^$/d' | sort | tail -n 1
)
while true; do
  review_found=$(gh api repos/{owner}/{repo}/pulls/"$pr_number"/reviews \
    --paginate \
    --jq ".[] | select(.user.login == \"chatgpt-codex-connector[bot]\" and .commit_id == \"$head_sha\" and .submitted_at != null and .state != \"DISMISSED\" and .state != \"PENDING\" and (\"$latest_review_request_at\" == \"\" or .submitted_at > \"$latest_review_request_at\")) | .id" \
    | head -n 1)
  [ -n "$review_found" ] && break

  issue_review=$(gh api repos/{owner}/{repo}/issues/"$pr_number"/comments \
    --paginate \
    --jq ".[] | select(((.body // \"\") | (startswith(\"## Codex Review\") or startswith(\"Codex Review:\"))) and ((.body // \"\") | contains(\"**Reviewed commit:** \`$head_marker\`\")) and (\"$issue_review_not_before\" == \"\" or .created_at > \"$issue_review_not_before\")) | [.id, .created_at] | @tsv" \
    | tail -n 1)
  if [ -n "$issue_review" ]; then
    issue_review_id=${issue_review%%$'\t'*}
    issue_review_at=${issue_review#*$'\t'}
    issue_ack=$(gh api repos/{owner}/{repo}/issues/"$pr_number"/comments \
      --paginate \
      --jq ".[] | select(.user.login == \"$acknowledger_login\" and .created_at > \"$issue_review_at\" and ((.body // \"\") | startswith(\"[codex] Review $issue_review_id:\"))) | .id" \
      | head -n 1)
    [ -n "$issue_ack" ] && break
  fi
  sleep 10
done

# Watch checks
if ! gh pr checks --watch; then
  gh pr checks
  # Identify failing run and inspect logs
  # gh run list --branch "$branch"
  # gh run view <run-id> --log
  exit 1
fi

# Squash-merge (remote branches auto-delete on merge in this repo)
gh pr merge --squash --subject "$pr_title" --body "$pr_body"
```

## Async Watch Helper

Preferred: use the asyncio watcher to monitor review comments, CI, and head
updates in parallel:

```
python3 .codex/skills/land/land_watch.py
```

Exit codes:

- 2: Review comments detected (address feedback)
- 3: CI checks failed
- 4: PR head updated (autofix commit detected)
- 5: PR has merge conflicts (merge `origin/main`, resolve, and push)

## Failure Handling

- If checks fail, pull details with `gh pr checks` and `gh run view --log`, then
  fix locally, commit with the `commit` skill, push with the `push` skill, and
  re-run the watch.
- Use judgment to identify flaky failures. If a failure is a flake (e.g., a
  timeout on only one platform), you may proceed without fixing it.
- If CI pushes an auto-fix commit (authored by GitHub Actions), it does not
  trigger a fresh CI run. Detect the updated PR head, pull locally, merge
  `origin/main` if needed, add a real author commit, and force-push to retrigger
  CI, then restart the checks loop.
- If mergeability is `UNKNOWN`, wait and re-check.
- Do not merge while review comments (human or Codex review) are outstanding.
- Codex review jobs retry on failure and are non-blocking. Use a submitted Codex
  review for the current PR head as the completion signal, not job status;
  treat its issue and inline comments as blocking feedback.
- Do not enable auto-merge; this repo has no required checks so auto-merge can
  skip tests.
- If the remote PR branch advanced due to your own prior force-push or merge,
  avoid redundant merges; re-run the formatter locally if needed and
  `git push --force-with-lease`.

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
  The watcher accepts only results newer than the current head's first check
  run and any later explicit review request. If the result body is edited after
  acknowledgement, inspect the changed content and post a fresh exact-ID
  acknowledgement.
- For a substantive top-level automated review body, acknowledge its findings
  with a root-level `[codex]` issue comment. Generated Copilot “Pull request
  overview” summaries are informational and do not require acknowledgement.
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
