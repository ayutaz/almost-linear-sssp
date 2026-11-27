---
name: git-workflow
description: >-
  Commit, integrate and push safely in this repo — local minimum gate, commit granularity, Conventional Commits, offloading the full gate of a task branch to GitHub Actions (scripts/ci/offload.sh + wait.py on wip/<name>), --no-ff integration, push, and fixing CI failures through TDD. Use when asked to commit, push, merge a task branch, open a PR, or check CI.
when_to_use: >-
  "commit this", "push", "merge the branch", "offload CI", "is CI green", "CI failed".
---
# Git workflow

Commit or push only when the user (or the lead acting on the user's request) asked for it.
Forbidden everywhere: `--force` / `--force-with-lease`, rewriting pushed commits,
`git remote add|set-url|remove`, `--no-verify`, `|| true`, skipped tests, disabling a job.

## 1. Local gate (before every commit)
Full gate: `make verify` (+ `make sanitizer` before merging algorithmic changes).
On a task branch whose full gate is offloaded (section 3), the local minimum is:
```bash
<the task's own tests>                         # ctest -R ... / uv run pytest tests/...
make format-check
scripts/ci/changed_tidy.sh main                # clang-tidy on changed src/ tools/ .cpp
uv run python scripts/tdd.py check
```
Never drop clang-tidy from the local check: once it was skipped and main's lint job stayed
red unnoticed. `changed_tidy.sh` uses `make tidy`'s settings (`build/$PRESET`, default debug).

## 2. Commit
- One task ID per commit where practical; tests and implementation of a cycle together.
- Stage explicit paths, never `git add -A` / `git add .`; check `git diff --cached --name-only`.
- Conventional Commits: `feat(paper): ordered merge events (Sec 3.3) [P3-MERGE-EVENTS]`.
  Types: feat, fix, test, refactor, perf, docs, build, ci, chore. Scope = module.
- Then `uv run python scripts/tdd.py done <ID> --commit "$(git rev-parse --short HEAD)"` and
  commit the evidence file as `docs(tdd): close <ID>`.

## 3. Offload the full gate to GitHub Actions (task branches)
```bash
scripts/ci/offload.sh <name> [<sha>]           # prints the wait command below
scripts/ci/wait.py <sha>... --since <ISO>      # exit 0 only if every run is success
```
What `offload.sh` does, and why:
1. `scripts/check_public_safety.py` on the files changed since the merge base with main (as
   of `<sha>`), plus a local-path check of the commit messages. Any finding: nothing pushed.
2. Pushes `<sha>` to `wip/<name>`, never forced. If `wip/<name>` exists and is not an
   ancestor of `<sha>`, it pushes to `wip/<name>-<sha7>` instead.
   Destination: `$SSSP_REPO_URL`, default `origin`. A checkout without a remote gets the URL
   through this variable; never add or change a remote.
3. `gh workflow run ci.yml --ref main -f os=ubuntu -f ref=<sha>` (same for `sanitizer.yml`).
   `--ref main` takes the workflow definition from main; `ref` is the commit tested;
   `os=ubuntu` drops the macOS jobs (that runner queue jams). The run title contains `<sha>`.
   `SSSP_DRY_RUN=1` prints the push and dispatch commands instead.

Waiting (`wait.py`):
- One `gh run list` per workflow every 90 s, shared by all SHAs. Never run several
  `gh run watch` at once: that hit the API limit (5000 requests/hour).
- Runs are matched by the SHA in their title, and only runs created after `--since` count
  (an old cancelled run of the same SHA was once taken for the new one). Use the `since=`
  that `offload.sh` printed.
- Repository: `$SSSP_REPO` (owner/name) or `gh repo view`. Exit 1 = some run not success,
  2 = timeout (`--timeout-min`, default 240) or unknown repository.

Workflow facts the tools rely on (do not break them in `.github/workflows/`):
- `concurrency.group` includes `inputs.ref` (benchmark: runner and suite); without it,
  parallel dispatches cancel each other.
- The Debug build job hit its 45-minute limit; ci.yml splits ctest with
  `matrix.shard: [labels, slow]` over two runners. Add a shard rather than raise the timeout.

## 4. Integrate into main
Only when both `ci` and `sanitizer` are green for the branch head:
```bash
scripts/ci/changed_tidy.sh main                # in the task worktree, before merging
git merge --no-ff -m "merge: <summary> [<ID>]" <task-branch>   # keeps the evidence SHAs valid
git push "${SSSP_REPO_URL:-origin}" --delete wip/<name>        # only the wip branch you pushed
```
Never rebase or squash a task branch: `docs/tdd/tasks/*.json` record its commit SHAs.

## 5. Push main
```bash
URL=${SSSP_REPO_URL:-origin}
old=$(git ls-remote "$URL" refs/heads/main | cut -f1)   # behind? rebase local-only commits only
scripts/ci/changed_tidy.sh "$old"              # changed .cpp since the pushed main
uv run python scripts/check_public_safety.py
git push "${SSSP_REPO_URL:-origin}" main
gh run list --branch main --limit 5            # then one `gh run watch <id> --exit-status`
```
Record the result (run id, conclusion) in `docs/PROGRESS.md` -> "CI".

## 6. CI failed
1. `gh run view <run-id> --log-failed`; reproduce locally with the same preset (`ci` /
   `ci-release` / `sanitizer`): `cmake --preset ci && cmake --build --preset ci && ctest --preset ci`.
2. Platform-only failure (Linux vs macOS): note compiler/version from the log.
3. Fix through `/tdd`: the failing test (or the CI reproduction command) as RED, then GREEN.
4. New commit (`fix(ci): ...` / `fix(<scope>): ...`); never amend a pushed commit, never
   disable the job or test, never `continue-on-error`. Offload or push again and wait.

## 7. Shell pitfalls
- zsh: `"$S:refs/heads/x"` applies the `:r` modifier to `$S`. Write `"${S}:refs/heads/x"`.
- Tools contain no local absolute paths, personal names or e-mail addresses; pass
  machine-specific values through environment variables (`SSSP_REPO_URL`, `SSSP_REPO`).
