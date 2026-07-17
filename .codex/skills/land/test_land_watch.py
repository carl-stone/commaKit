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
        "body": "review complete",
    }
    review.update(overrides)
    return review


class ReviewGateTests(unittest.TestCase):
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
            with (
                patch.object(
                    LAND_WATCH,
                    "fetch_review_context",
                    review_context,
                ),
                patch.object(LAND_WATCH, "POLL_SECONDS", 0),
            ):
                await asyncio.wait_for(
                    LAND_WATCH.wait_for_codex(1, HEAD_SHA, checks_done),
                    timeout=0.2,
                )
            self.assertEqual(calls, 2)

        asyncio.run(run())


if __name__ == "__main__":
    unittest.main()
