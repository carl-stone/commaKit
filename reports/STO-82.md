# STO-82 implementation handoff

## Summary

Refactored `.applyKeepFilter()` so `overlap` and `metagene` filter their
retained core association and metadata list-columns through one shared helper.
The helper applies each site's `rel_position == 0L` indices with
`S4Vectors::mendoapply()`, preserving the typed list subclasses and the
parallel association ordering.

Added a `keep = "metagene"` contract test with two overlapping features and
one nearby non-overlapping feature. It verifies metadata passthrough remains a
`CharacterList`, fractional positions remain a `NumericList`, matching
metadata values preserve feature-name order, and the nearby feature is
removed.

## Intended branch

`symphony/STO-82`

## Files changed

- `R/annotateSites.R`
- `tests/testthat/test-annotateSites.R`
- `reports/STO-82.md`

## Validation

- Passed: `Rscript -e "testthat::test_file('tests/testthat/test-annotateSites.R')"`
  - Completed successfully with no testthat failures or warnings emitted.
- Passed: `git diff --check`
  - No whitespace errors.

## Blockers

None.

## Next owner

Factory finalizer: inspect the unstaged diff, create the reviewable PR, and
route it to the factory steward / Carl reviewer. No commit, push, PR, or merge
was performed by this worker.

## Workspace notes

Pre-existing untracked `.symphony/` and `logs/` paths were left untouched.
