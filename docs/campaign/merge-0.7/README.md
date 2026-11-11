# Implementation log: one engine, one API (release 0.7.0)

This log records how the plan in
[`docs/plan-0.7.0-engine-merge.md`](../../plan-0.7.0-engine-merge.md) was
carried out. Each task has its own log in this directory; this page gives the
order, the outcome of each task and every deviation from the plan.

## Execution notes

- All tasks were committed on one integration branch, `merge-0.7`, and merged
  into `main` with a single pull request, followed by a release-preparation
  pull request. This replaces the plan's rule of one pull request per task.
- The stopgap branch `docs/two-engines` (task M0a) was not merged: the merged
  release replaces the explanation of two engines it added.

## Tasks

| Task | Log | Outcome |
| --- | --- | --- |
