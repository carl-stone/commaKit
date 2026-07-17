import asyncio
import importlib.util
import unittest
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch


MODULE_PATH = Path(__file__).with_name("land_watch.py")
SPEC = importlib.util.spec_from_file_location("land_watch", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Unable to load {MODULE_PATH}")
LAND_WATCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LAND_WATCH)

HEAD_SHA = "a" * 40
ISSUE_REVIEW_BODY = (
    "## Codex Review — correctness\n\n"
    f"**Reviewed commit:** `{HEAD_SHA[:10]}`"
)


def utc_time(hour: int, minute: int) -> datetime:
    return datetime(2026, 7, 17, hour, minute, tzinfo=timezone.utc)


def codex_review(**overrides: object) -> dict[str, object]:
    review: dict[str, object] = {
        "user": {
            "login": "chatgpt-codex-connector[bot]",
            "type": "Bot",
        },
        "commit_id": HEAD_SHA,
        "submitted_at": "2026-07-17T04:01:00Z",
        "state": "COMMENTED",
        "body": "",
    }
    review.update(overrides)
    return review


class ReviewGateTests(unittest.TestCase):
    def test_check_boundary_prefers_started_time(self) -> None:
        self.assertEqual(
            LAND_WATCH.check_boundary_timestamp(
                {
                    "started_at": "2026-07-17T04:00:00Z",
                    "completed_at": "2026-07-17T04:02:00Z",
                },
            ),
            utc_time(4, 0),
        )

    def test_accepts_current_head_review(self) -> None:
        self.assertTrue(
            LAND_WATCH.has_current_codex_review(
                [codex_review()],
                HEAD_SHA,
                utc_time(4, 0),
            ),
        )

    def test_rejects_stale_head_review(self) -> None:
        self.assertFalse(
            LAND_WATCH.has_current_codex_review(
                [codex_review()],
                "b" * 40,
                None,
            ),
        )

    def test_rejects_review_not_newer_than_explicit_request(self) -> None:
        self.assertFalse(
            LAND_WATCH.has_current_codex_review(
                [codex_review()],
                HEAD_SHA,
                utc_time(4, 1),
            ),
        )

    def test_rejects_dismissed_review(self) -> None:
        self.assertFalse(
            LAND_WATCH.has_current_codex_review(
                [codex_review(state="DISMISSED")],
                HEAD_SHA,
                None,
            ),
        )

    def test_rejects_pending_unsubmitted_review(self) -> None:
        review = codex_review(
            state="PENDING",
            submitted_at=None,
            created_at="2026-07-17T04:01:00Z",
        )
        self.assertFalse(
            LAND_WATCH.has_current_codex_review(
                [review],
                HEAD_SHA,
                None,
            ),
        )

    def test_generic_github_actions_review_cannot_complete_gate(self) -> None:
        review = codex_review(
            user={"login": "github-actions[bot]", "type": "Bot"},
        )
        self.assertFalse(
            LAND_WATCH.has_current_codex_review(
                [review],
                HEAD_SHA,
                None,
            ),
        )

    def test_review_request_uses_immutable_creation_time(self) -> None:
        request = {
            "user": {"login": "carl-stone", "type": "User"},
            "body": "@codex review",
            "created_at": "2026-07-17T04:00:00Z",
            "updated_at": "2026-07-17T05:00:00Z",
        }
        self.assertEqual(
            LAND_WATCH.latest_review_request_at([request]),
            utc_time(4, 0),
        )

    def test_copilot_overview_is_informational(self) -> None:
        review = {
            "user": {
                "login": "copilot-pull-request-reviewer[bot]",
                "type": "Bot",
            },
            "submitted_at": "2026-07-17T04:01:00Z",
            "state": "COMMENTED",
            "body": "## Pull request overview\nInformational summary",
        }
        self.assertFalse(LAND_WATCH.is_blocking_review(review, None))

    def test_codex_review_wrapper_is_informational(self) -> None:
        review = codex_review(
            body="### 💡 Codex Review\nAutomated review wrapper",
        )
        self.assertFalse(LAND_WATCH.is_blocking_review(review, None))

    def test_substantive_codex_review_requires_acknowledgement(self) -> None:
        review = codex_review(body="Please fix this correctness issue.")
        self.assertTrue(LAND_WATCH.is_blocking_review(review, None))
        self.assertFalse(
            LAND_WATCH.is_blocking_review(
                review,
                None,
                utc_time(4, 2),
            ),
        )

    def test_later_wrapper_does_not_hide_substantive_bot_review(self) -> None:
        finding = codex_review(body="Please fix this correctness issue.")
        wrapper = codex_review(
            body="### 💡 Codex Review\nAutomated review wrapper",
            submitted_at="2026-07-17T04:02:00Z",
        )
        self.assertEqual(
            LAND_WATCH.filter_blocking_reviews(
                [finding, wrapper],
                None,
            ),
            [finding],
        )

    def test_submitted_review_requires_exact_id_acknowledgement(self) -> None:
        finding = codex_review(
            id=55,
            body="Please fix this correctness issue.",
        )
        unrelated_ack = {101: utc_time(4, 2)}
        exact_ack = {55: utc_time(4, 2)}
        self.assertEqual(
            LAND_WATCH.filter_blocking_reviews(
                [finding],
                None,
                unrelated_ack,
            ),
            [finding],
        )
        self.assertEqual(
            LAND_WATCH.filter_blocking_reviews(
                [finding],
                None,
                exact_ack,
            ),
            [],
        )

    def test_later_bot_approval_clears_change_request(self) -> None:
        change_request = codex_review(
            state="CHANGES_REQUESTED",
            body="Please fix this correctness issue.",
        )
        approval = codex_review(
            state="APPROVED",
            body="",
            submitted_at="2026-07-17T04:02:00Z",
        )
        self.assertEqual(
            LAND_WATCH.filter_blocking_reviews(
                [change_request, approval],
                None,
            ),
            [],
        )

    def test_prior_head_bot_change_request_is_nonblocking(self) -> None:
        stale_request = codex_review(
            commit_id="b" * 40,
            state="CHANGES_REQUESTED",
            body="Please fix this old correctness issue.",
        )
        self.assertEqual(
            LAND_WATCH.filter_blocking_reviews(
                [stale_request],
                None,
                head_sha=HEAD_SHA,
            ),
            [],
        )

    def test_substantive_bot_review_requires_acknowledgement(self) -> None:
        review = {
            "user": {
                "login": "copilot-pull-request-reviewer[bot]",
                "type": "Bot",
            },
            "submitted_at": "2026-07-17T04:01:00Z",
            "state": "COMMENTED",
            "body": "Please fix this correctness issue.",
        }
        self.assertTrue(LAND_WATCH.is_blocking_review(review, None))
        self.assertFalse(
            LAND_WATCH.is_blocking_review(
                review,
                None,
                utc_time(4, 2),
            ),
        )

    def test_dismissed_bot_review_is_nonblocking(self) -> None:
        review = {
            "user": {
                "login": "copilot-pull-request-reviewer[bot]",
                "type": "Bot",
            },
            "submitted_at": "2026-07-17T04:01:00Z",
            "state": "DISMISSED",
            "body": "Please fix this correctness issue.",
        }
        self.assertFalse(LAND_WATCH.is_blocking_review(review, None))

    def test_codex_request_does_not_hide_older_copilot_comment(self) -> None:
        comment = {
            "id": 1,
            "user": {
                "login": "copilot-pull-request-reviewer[bot]",
                "type": "Bot",
            },
            "created_at": "2026-07-17T04:00:00Z",
            "body": "Please fix this correctness issue.",
            "pull_request_review_id": 1,
        }
        self.assertEqual(
            LAND_WATCH.filter_codex_comments(
                [comment],
                HEAD_SHA,
                utc_time(4, 1),
                "carl-stone",
            ),
            [comment],
        )

    def test_acknowledgement_uses_immutable_creation_time(self) -> None:
        acknowledgement = {
            "user": {"login": "carl-stone", "type": "User"},
            "body": "[codex] acknowledged",
            "created_at": "2026-07-17T04:00:00Z",
            "updated_at": "2026-07-17T05:00:00Z",
        }
        self.assertEqual(
            LAND_WATCH.latest_codex_issue_reply_time(
                [acknowledgement],
                "carl-stone",
            ),
            utc_time(4, 0),
        )

    def test_acknowledged_issue_review_can_complete_gate(self) -> None:
        comments = [
            {
                "id": 101,
                "user": {
                    "login": "github-actions[bot]",
                    "type": "Bot",
                },
                "body": ISSUE_REVIEW_BODY,
                "created_at": "2026-07-17T04:01:00Z",
            },
            {
                "id": 102,
                "user": {"login": "carl-stone", "type": "User"},
                "body": "[codex] Review 101: acknowledged",
                "created_at": "2026-07-17T04:02:00Z",
            },
        ]
        self.assertTrue(
            LAND_WATCH.has_acknowledged_codex_issue_review(
                comments,
                HEAD_SHA,
                utc_time(4, 0),
                "carl-stone",
            ),
        )
        self.assertFalse(
            LAND_WATCH.has_acknowledged_codex_issue_review(
                comments[:1],
                HEAD_SHA,
                utc_time(4, 0),
                "carl-stone",
            ),
        )

    def test_clean_issue_review_variant_is_recognized(self) -> None:
        self.assertTrue(
            LAND_WATCH.is_codex_review_body(
                "Codex Review: Didn't find any major issues. :tada:",
            ),
        )

    def test_markerless_issue_review_is_not_current_head_feedback(self) -> None:
        comment = {
            "id": 101,
            "user": {
                "login": "github-actions[bot]",
                "type": "Bot",
            },
            "body": "Codex Review: Didn't find any major issues. :tada:",
            "created_at": "2026-07-17T04:01:00Z",
        }
        self.assertEqual(
            LAND_WATCH.filter_codex_comments(
                [comment],
                HEAD_SHA,
                None,
                "carl-stone",
            ),
            [],
        )

    def test_issue_review_must_be_newer_than_head(self) -> None:
        comments = [
            {
                "id": 101,
                "user": {
                    "login": "github-actions[bot]",
                    "type": "Bot",
                },
                "body": ISSUE_REVIEW_BODY,
                "created_at": "2026-07-17T04:01:00Z",
            },
            {
                "id": 102,
                "user": {"login": "carl-stone", "type": "User"},
                "body": "[codex] Review 101: acknowledged",
                "created_at": "2026-07-17T04:02:00Z",
            },
        ]
        self.assertFalse(
            LAND_WATCH.has_acknowledged_codex_issue_review(
                comments,
                HEAD_SHA,
                utc_time(4, 3),
                "carl-stone",
            ),
        )

    def test_bot_identity_can_acknowledge_exact_issue_review(self) -> None:
        comments = [
            {
                "id": 101,
                "user": {
                    "login": "github-actions[bot]",
                    "type": "Bot",
                },
                "body": ISSUE_REVIEW_BODY,
                "created_at": "2026-07-17T04:01:00Z",
            },
            {
                "id": 102,
                "user": {"login": "landing-agent[bot]", "type": "Bot"},
                "body": "[codex] Review 101: acknowledged",
                "created_at": "2026-07-17T04:02:00Z",
            },
        ]
        self.assertTrue(
            LAND_WATCH.has_acknowledged_codex_issue_review(
                comments,
                HEAD_SHA,
                utc_time(4, 0),
                "landing-agent[bot]",
            ),
        )

    def test_unrelated_ack_cannot_complete_issue_review(self) -> None:
        comments = [
            {
                "id": 101,
                "user": {
                    "login": "github-actions[bot]",
                    "type": "Bot",
                },
                "body": ISSUE_REVIEW_BODY,
                "created_at": "2026-07-17T04:01:00Z",
            },
            {
                "id": 102,
                "user": {"login": "carl-stone", "type": "User"},
                "body": "[codex] unrelated status update",
                "created_at": "2026-07-17T04:02:00Z",
            },
        ]
        self.assertFalse(
            LAND_WATCH.has_acknowledged_codex_issue_review(
                comments,
                HEAD_SHA,
                utc_time(4, 0),
                "carl-stone",
            ),
        )

    def test_untrusted_author_cannot_acknowledge_issue_review(self) -> None:
        comments = [
            {
                "id": 101,
                "user": {
                    "login": "github-actions[bot]",
                    "type": "Bot",
                },
                "body": ISSUE_REVIEW_BODY,
                "created_at": "2026-07-17T04:01:00Z",
            },
            {
                "id": 102,
                "user": {"login": "mallory", "type": "User"},
                "body": "[codex] Review 101: acknowledged",
                "created_at": "2026-07-17T04:02:00Z",
            },
        ]
        self.assertFalse(
            LAND_WATCH.has_acknowledged_codex_issue_review(
                comments,
                HEAD_SHA,
                utc_time(4, 0),
                "carl-stone",
            ),
        )

    def test_prior_head_marker_cannot_complete_issue_review(self) -> None:
        comments = [
            {
                "id": 101,
                "user": {
                    "login": "github-actions[bot]",
                    "type": "Bot",
                },
                "body": ISSUE_REVIEW_BODY.replace("a" * 10, "b" * 10),
                "created_at": "2026-07-17T04:01:00Z",
            },
            {
                "id": 102,
                "user": {"login": "carl-stone", "type": "User"},
                "body": "[codex] Review 101: acknowledged",
                "created_at": "2026-07-17T04:02:00Z",
            },
        ]
        self.assertFalse(
            LAND_WATCH.has_acknowledged_codex_issue_review(
                comments,
                HEAD_SHA,
                utc_time(4, 0),
                "carl-stone",
            ),
        )
        self.assertEqual(
            LAND_WATCH.filter_codex_review_issue_comments(
                comments,
                HEAD_SHA,
                None,
                "carl-stone",
            ),
            [],
        )

    def test_newer_issue_review_also_requires_exact_ack(self) -> None:
        comments = [
            {
                "id": 101,
                "user": {
                    "login": "github-actions[bot]",
                    "type": "Bot",
                },
                "body": ISSUE_REVIEW_BODY,
                "created_at": "2026-07-17T04:01:00Z",
            },
            {
                "id": 102,
                "user": {"login": "carl-stone", "type": "User"},
                "body": "[codex] Review 101: acknowledged",
                "created_at": "2026-07-17T04:02:00Z",
            },
            {
                "id": 103,
                "user": {
                    "login": "github-actions[bot]",
                    "type": "Bot",
                },
                "body": ISSUE_REVIEW_BODY.replace("correctness", "security"),
                "created_at": "2026-07-17T04:03:00Z",
            },
        ]
        self.assertFalse(
            LAND_WATCH.has_acknowledged_codex_issue_review(
                comments,
                HEAD_SHA,
                utc_time(4, 0),
                "carl-stone",
            ),
        )

    def test_issue_review_before_latest_request_is_not_blocking(self) -> None:
        comment = {
            "id": 101,
            "user": {
                "login": "github-actions[bot]",
                "type": "Bot",
            },
            "body": ISSUE_REVIEW_BODY,
            "created_at": "2026-07-17T04:01:00Z",
        }
        self.assertEqual(
            LAND_WATCH.filter_codex_review_issue_comments(
                [comment],
                HEAD_SHA,
                utc_time(4, 2),
                "carl-stone",
            ),
            [],
        )

    def test_exact_ack_must_follow_issue_review(self) -> None:
        comments = [
            {
                "id": 100,
                "user": {"login": "carl-stone", "type": "User"},
                "body": "[codex] Review 101: acknowledged",
                "created_at": "2026-07-17T04:01:00Z",
            },
            {
                "id": 101,
                "user": {
                    "login": "github-actions[bot]",
                    "type": "Bot",
                },
                "body": ISSUE_REVIEW_BODY,
                "created_at": "2026-07-17T04:02:00Z",
            },
        ]
        self.assertFalse(
            LAND_WATCH.has_acknowledged_codex_issue_review(
                comments,
                HEAD_SHA,
                utc_time(4, 0),
                "carl-stone",
            ),
        )
        self.assertEqual(
            LAND_WATCH.filter_codex_review_issue_comments(
                comments,
                HEAD_SHA,
                None,
                "carl-stone",
            ),
            [comments[1]],
        )

    def test_issue_review_edit_after_ack_requires_fresh_ack(self) -> None:
        comments = [
            {
                "id": 101,
                "user": {
                    "login": "github-actions[bot]",
                    "type": "Bot",
                },
                "body": ISSUE_REVIEW_BODY,
                "created_at": "2026-07-17T04:01:00Z",
                "updated_at": "2026-07-17T04:03:00Z",
            },
            {
                "id": 102,
                "user": {"login": "carl-stone", "type": "User"},
                "body": "[codex] Review 101: acknowledged",
                "created_at": "2026-07-17T04:02:00Z",
            },
        ]
        self.assertEqual(
            LAND_WATCH.filter_codex_comments(
                comments,
                HEAD_SHA,
                None,
                "carl-stone",
            ),
            [comments[0]],
        )

    def test_green_checks_wait_for_current_head_review(self) -> None:
        async def run() -> None:
            calls = 0

            async def review_context(
                _pr_number: int,
            ) -> tuple[list[object], list[object], list[object], datetime]:
                nonlocal calls
                calls += 1
                reviews = [] if calls == 1 else [codex_review()]
                return [], [], reviews, utc_time(4, 0)

            checks_done = asyncio.Event()
            checks_done.set()
            head_review_boundary = asyncio.get_running_loop().create_future()
            head_review_boundary.set_result(utc_time(4, 0))
            with (
                patch.object(
                    LAND_WATCH,
                    "fetch_review_context",
                    review_context,
                ),
                patch.object(LAND_WATCH, "POLL_SECONDS", 0),
            ):
                await asyncio.wait_for(
                    LAND_WATCH.wait_for_codex(
                        1,
                        HEAD_SHA,
                        head_review_boundary,
                        checks_done,
                        "carl-stone",
                    ),
                    timeout=0.2,
                )
            self.assertEqual(calls, 2)

        asyncio.run(run())

    def test_issue_review_without_request_completes_after_head(self) -> None:
        async def run() -> None:
            issue_comments = [
                {
                    "id": 101,
                    "user": {
                        "login": "github-actions[bot]",
                        "type": "Bot",
                    },
                    "body": ISSUE_REVIEW_BODY,
                    "created_at": "2026-07-17T04:01:00Z",
                },
                {
                    "id": 102,
                    "user": {"login": "carl-stone", "type": "User"},
                    "body": "[codex] Review 101: acknowledged",
                    "created_at": "2026-07-17T04:02:00Z",
                },
            ]

            async def review_context(
                _pr_number: int,
            ) -> tuple[list[object], list[object], list[object], None]:
                return issue_comments, [], [], None

            checks_done = asyncio.Event()
            checks_done.set()
            head_review_boundary = asyncio.get_running_loop().create_future()
            head_review_boundary.set_result(utc_time(4, 0))
            with patch.object(
                LAND_WATCH,
                "fetch_review_context",
                review_context,
            ):
                await asyncio.wait_for(
                    LAND_WATCH.wait_for_codex(
                        1,
                        HEAD_SHA,
                        head_review_boundary,
                        checks_done,
                        "carl-stone",
                    ),
                    timeout=0.2,
                )

        asyncio.run(run())


if __name__ == "__main__":
    unittest.main()
