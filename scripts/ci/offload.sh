#!/usr/bin/env bash
# usage: scripts/ci/offload.sh <ref-name> [<sha>]      (run inside the repository / worktree)
#
# Offloads the full quality gate of a task branch to GitHub Actions (CLAUDE.md "Quality gate"):
#   1. public-safety audit of the files changed since the merge base with main, and of the
#      commit messages (local paths); stops on any finding;
#   2. pushes <sha> (default HEAD) to wip/<ref-name>, never forced. If wip/<ref-name> already
#      exists and is not an ancestor of <sha>, pushes to wip/<ref-name>-<sha7> instead;
#   3. dispatches the ci and sanitizer workflows with os=ubuntu and ref=<sha>, using main's
#      workflow definition (--ref main). The run titles contain <sha>.
# Then wait with the command it prints (scripts/ci/wait.py <sha> --since <time>).
#
# Environment:
#   SSSP_REPO_URL  push destination (URL); default: the remote `origin`. No remote is ever
#                  added or changed: a checkout without a remote passes the URL here.
#   SSSP_REPO      owner/name for gh; default: derived from a github.com SSSP_REPO_URL, else gh's
#                  own default repository.
#   SSSP_BASE      base branch for the audit (default: main).
#   SSSP_DRY_RUN=1 print the push and dispatch commands instead of running them.
set -euo pipefail

die() { echo "offload: $*" >&2; exit 1; }
run() { if [ "${SSSP_DRY_RUN:-0}" = 1 ]; then echo "[dry-run] $*"; else "$@"; fi; }

[ $# -ge 1 ] && [ $# -le 2 ] || die "usage: scripts/ci/offload.sh <ref-name> [<sha>]"
name=$1
case $name in "" | wip/* | *[[:space:]]* | *..* | -*) die "bad ref-name: '$name'" ;; esac

top=$(git rev-parse --show-toplevel)
cd "$top"
S=$(git rev-parse --verify "${2:-HEAD}^{commit}") || die "not a commit: ${2:-HEAD}"
url=${SSSP_REPO_URL:-origin}
base_branch=${SSSP_BASE:-main}
base=$(git merge-base "$base_branch" "$S") || die "no merge base between $base_branch and $S"

repo=${SSSP_REPO:-}
if [ -z "$repo" ] && [[ $url =~ github\.com[:/]([^/]+/[^/]+)$ ]]; then
  repo=${BASH_REMATCH[1]%.git}
fi
gh_repo=()
[ -n "$repo" ] && gh_repo=(-R "$repo")

# --- 1. public safety (content as of <sha>, not the working tree) ---
if git log --format=%B "${base}..${S}" | grep -qE '/Users/|/home/[a-z]|/private/(tmp|var)/'; then
  die "a commit message in ${base:0:7}..${S:0:7} contains a local path: not pushed"
fi
changed=$(git diff --name-only --diff-filter=d "$base" "$S")
if [ -n "$changed" ]; then
  snap=$(mktemp -d)
  trap 'rm -rf "$snap"' EXIT
  # shellcheck disable=SC2086 # file names are git paths, one per line
  git archive "$S" -- $changed | tar -x -C "$snap"
  # shellcheck disable=SC2086
  uv run python scripts/check_public_safety.py --root "$snap" $changed \
    || die "public-safety findings: not pushed"
fi

# --- 2. push, never forced ---
ref=wip/$name
remote_sha=$(git ls-remote "$url" "refs/heads/${ref}") || die "cannot reach $url"
remote_sha=${remote_sha%%[[:space:]]*}
if [ -n "$remote_sha" ] && ! git merge-base --is-ancestor "$remote_sha" "$S" 2>/dev/null; then
  ref=wip/${name}-${S:0:7}
  echo "offload: wip/${name} is at ${remote_sha:0:7}, not an ancestor of ${S:0:7}; using ${ref}"
fi
# Taken before the dispatch, minus 60 s of clock skew: wait.py ignores runs older than this.
t=$(($(date +%s) - 60))
since=$(date -u -r "$t" +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || date -u -d "@$t" +%Y-%m-%dT%H:%M:%SZ)
run git push "$url" "${S}:refs/heads/${ref}"

# --- 3. dispatch main's workflows on <sha>, Linux only ---
for wf in ci sanitizer; do
  run gh workflow run "${wf}.yml" ${gh_repo[@]+"${gh_repo[@]}"} --ref main -f os=ubuntu -f ref="$S" \
    || die "dispatch of ${wf}.yml failed"
done

echo "ref=${ref} sha=${S} since=${since}"
echo "wait: ${repo:+SSSP_REPO=$repo }scripts/ci/wait.py ${S} --since ${since}"
