"""Tests for .claude/hooks/pre_tool_use.py (PreToolUse guard).

The hook is run exactly like Claude Code runs it (subprocess, event JSON on
stdin). It must exit 0 and either print nothing (no opinion), a deny decision
or an additionalContext warning. It must never answer "allow".
"""

from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path
from typing import Any

import pytest
from harness_util import GREEN_EVENT, RED_EVENT, write_task

TEST_FILE = "tests/unit/test_a.cpp"


# ----------------------------------------------------------------------- helpers
def bash(cmd: str, **extra: Any) -> dict[str, Any]:
    return {
        "hook_event_name": "PreToolUse",
        "tool_name": "Bash",
        "tool_input": {"command": cmd},
        "tool_use_id": "toolu_test",
        "permission_mode": "default",
        **extra,
    }


def file_tool(tool: str, **tool_input: Any) -> dict[str, Any]:
    return {
        "hook_event_name": "PreToolUse",
        "tool_name": tool,
        "tool_input": tool_input,
        "tool_use_id": "toolu_test",
        "permission_mode": "default",
    }


def reason(res: Any) -> str:
    data = res.json or {}
    return (data.get("hookSpecificOutput") or {}).get("permissionDecisionReason") or ""


def assert_deny(res: Any) -> None:
    assert res.code == 0, res.stderr
    assert res.decision == "deny", (res.stdout, res.stderr)
    out = res.json["hookSpecificOutput"]
    assert out["hookEventName"] == "PreToolUse"
    assert reason(res).strip()


def assert_silent(res: Any) -> None:
    assert res.code == 0, res.stderr
    assert res.stdout.strip() == "", res.stdout


def set_enforce(root: Path, enforce: bool) -> None:
    (root / ".tdd" / "config.json").write_text(json.dumps({"enforce": enforce}))


@pytest.fixture
def env_tmp(tmp_path: Path) -> dict[str, str]:
    """A TMPDIR that does not contain the test repository."""
    d = tmp_path / "tmpdir"
    d.mkdir()
    return {"TMPDIR": str(d)}


# ------------------------------------------------------------ (a) Bash: deny
DANGEROUS = [
    # force push
    "git push --force",
    "git push -f origin main",
    "git push origin main -f",
    "git push --force-with-lease",
    "git push --force-with-lease=main:abc123 origin main",
    "git push --force-if-includes origin main",
    "git push -uf origin main",
    "git push origin +main",
    "git push origin +HEAD:main",
    "git push --mirror",
    "git push --forc origin main",
    # reset --hard
    "git reset --hard",
    "git reset --hard HEAD~1",
    "git reset -q --hard origin/main",
    # discarding uncommitted changes
    "git checkout .",
    "git checkout -- src/core/a.cpp",
    "git checkout HEAD~1 -- src/core/a.cpp",
    "git checkout -f main",
    "git checkout --force main",
    "git checkout src/core/a.cpp",
    "git restore .",
    "git restore src/core/a.cpp",
    "git restore --staged --worktree src/core/a.cpp",
    "git restore -SW src/core/a.cpp",
    "git restore",
    # clean
    "git clean -f",
    "git clean --force",
    "git clean -fd",
    "git clean -xdf",
    # branch force delete
    "git branch -D feature",
    "git branch -Dr origin/feature",
    "git branch --delete --force feature",
    "git branch -d -f feature",
    "git branch -df feature",
    # remote changes
    "git remote set-url origin https://example.invalid/x.git",
    "git remote remove origin",
    "git remote rm origin",
    "git remote rename origin upstream",
    "git remote add upstream https://example.invalid/x.git",
    # gh
    "gh repo edit --visibility public",
    "gh repo edit owner/repo --visibility=private --accept-visibility-change-consequences",
    "gh repo delete owner/repo --yes",
    # rm outside the repo / root / .git
    "rm -rf /",
    "rm -rf /*",
    "rm -rf ~",
    "rm -rf ..",
    "rm -rf .",
    "rm -rf ./",
    "rm -rf .git",
    "rm -rf .git/hooks",
    "rm -fr src/../..",
    "rm -f /etc/hosts",
    "rm --recursive --force /usr/local/x",
    "rm -rf $UNDEFINED_PRE_TOOL_USE_VAR/foo",
    "rm -rf $(pwd)",
    "rm -rf /tmp",
    # wrappers / prefixes / paths
    "FOO=1 git push -f",
    "A=1 B=2 git reset --hard",
    "sudo git push --force",
    "sudo -u root git reset --hard",
    "env A=1 git push -f",
    "env -i git reset --hard",
    "command git reset --hard",
    "builtin eval 'git reset --hard'",
    "nohup git reset --hard",
    "time git reset --hard",
    "nice -n 10 git reset --hard",
    "timeout 10 git push -f",
    "timeout -s KILL 10s git push -f",
    "stdbuf -oL git reset --hard",
    "noglob git reset --hard",
    "exec git reset --hard",
    "echo x | xargs git reset --hard",
    "echo / | xargs -I{} rm -rf /",
    "/usr/bin/git push --force",
    "git -C /some/where push --force",
    "git -c core.pager=cat push --force",
    "git --no-pager reset --hard",
    "git -P reset --hard",
    "git --git-dir=.git --work-tree=. reset --hard",
    "git --git-dir .git reset --hard",
    # compound commands
    "cd src && git push --force",
    "echo hi; git reset --hard",
    "echo hi\ngit reset --hard",
    "git status || git clean -fd",
    "git status | git reset --hard",
    "git fetch & git reset --hard",
    "(git reset --hard)",
    "{ git reset --hard; }",
    "if true; then git reset --hard; fi",
    "$(git reset --hard)",
    "echo $(git reset --hard)",
    'echo "$(git reset --hard)"',
    "echo `git reset --hard`",
    "x=$(git push -f origin main)",
    "echo hi # comment\ngit reset --hard",
    "echo hi \\\n; git reset --hard",
    # nested shells / eval
    'bash -c "git push --force"',
    "sh -c 'git reset --hard'",
    "zsh -c 'cd x && git clean -fd'",
    'bash -lc "git push -f"',
    "bash -c \"bash -c 'git reset --hard'\"",
    "eval git reset --hard",
    'eval "git push --force"',
    # unparsable (unterminated quote): conservative regex fallback
    'git push --force "unterminated',
    "echo 'oops ; git reset --hard",
    "rm -rf / 'x",
]


@pytest.mark.parametrize("cmd", DANGEROUS)
def test_dangerous_bash_is_denied(repo: Path, run_hook: Any, env_tmp: dict[str, str], cmd: str):
    res = run_hook("pre_tool_use", bash(cmd), repo, env=env_tmp)
    assert_deny(res)


SAFE = [
    "git status",
    "git push origin main",
    "git push -u origin main",
    "git push --no-force-with-lease origin main",
    "git push --follow-tags origin main",
    "git reset --soft HEAD~1",
    "git reset HEAD src/core/a.cpp",
    "git checkout -b x",
    "git checkout -b src",
    "git checkout main",
    "git switch -c feature",
    "git restore --staged src/core/a.cpp",
    "git restore -S src/core/a.cpp",
    "git clean -n",
    "git clean -nd",
    "git branch -d x",
    "git branch --delete x",
    "git branch -a",
    "git remote -v",
    "git remote",
    "git remote show origin",
    "git log --format=%H -n 1",
    "git commit -m 'reset --hard'",
    'git commit -m "reset --hard"',
    "git commit -m 'git push --force; rm -rf /'",
    'echo "git push --force"',
    "echo 'git reset --hard'",
    "printf '%s\\n' 'rm -rf /'",
    "grep -rn 'rm -rf /' docs",
    "echo '$(git reset --hard)'",
    "gh repo view",
    "gh pr create --title x --body y",
    "gh repo edit --description hello",
    "rm -rf build",
    "rm -rf build/debug src/old",
    "rm -rf ./build",
    "rm -rf /tmp/pre-tool-use-test-dir",
    "rm -rf /private/tmp/pre-tool-use-test-dir",  # public-safety: allow
    "rm -rf $TMPDIR/foo",
    "rm src/core/a.cpp",
    "rm .gitignore",
    "ls -la && echo done",
    "cat <<'EOF' > /dev/null\ngit push --force\nrm -rf /\nEOF",
    "git commit -F - <<'EOF'\nfix: avoid git reset --hard\nEOF",
    "git commit -m \"$(cat <<'EOF'\nfix: never git push --force\n\nrm -rf / is bad\nEOF\n)\"",
    "echo hi 2>&1 | tee /dev/null",
    "echo 'unterminated is fine",
    "cmake --build --preset debug && ctest --preset debug",
    "uv run pytest tests/harness -q",
    "bash scripts/bootstrap.sh",
]


@pytest.mark.parametrize("cmd", SAFE)
def test_safe_bash_is_silent(repo: Path, run_hook: Any, env_tmp: dict[str, str], cmd: str):
    res = run_hook("pre_tool_use", bash(cmd), repo, env=env_tmp)
    assert_silent(res)


def _under_system_tmp(p: Path) -> bool:
    # On Linux CI pytest's tmp_path lives below /tmp, which the hook allows.
    return str(p.resolve()).startswith(("/tmp/", "/private/tmp/"))  # public-safety: allow


def test_rm_inside_scratchpad_is_allowed(
    repo: Path, run_hook: Any, tmp_path: Path, env_tmp: dict[str, str]
):
    scratch = tmp_path / "scratchpad"
    scratch.mkdir()
    payload = bash(f"rm -rf {scratch}/x", scratchpad_dir=str(scratch))
    assert_silent(run_hook("pre_tool_use", payload, repo, env=env_tmp))
    if not _under_system_tmp(tmp_path):
        # without scratchpad_dir the same path is outside the repo -> denied
        assert_deny(run_hook("pre_tool_use", bash(f"rm -rf {scratch}/x"), repo, env=env_tmp))
        # deleting the scratchpad directory itself is not "below" it
        payload = bash(f"rm -rf {scratch}", scratchpad_dir=str(scratch))
        assert_deny(run_hook("pre_tool_use", payload, repo, env=env_tmp))


def test_rm_sibling_outside_repo_is_denied(
    repo: Path, run_hook: Any, tmp_path: Path, env_tmp: dict[str, str]
):
    if _under_system_tmp(tmp_path):
        pytest.skip("tmp_path is below /tmp, where rm -rf is allowed")
    assert_deny(run_hook("pre_tool_use", bash("rm -r ../other"), repo, env=env_tmp))


def test_rm_of_tmpdir_ancestor_of_repo_is_denied(repo: Path, run_hook: Any, tmp_path: Path):
    # The repo lives below TMPDIR: deleting TMPDIR would delete the repo.
    res = run_hook("pre_tool_use", bash(f"rm -rf {tmp_path}"), repo, env={"TMPDIR": str(tmp_path)})
    assert_deny(res)


def test_remote_add_origin_allowed_only_without_remotes(repo: Path, run_hook: Any):
    cmd = "git remote add origin https://example.invalid/x.git"
    assert_silent(run_hook("pre_tool_use", bash(cmd), repo))
    subprocess.run(
        ["git", "-C", str(repo), "remote", "add", "origin", "https://example.invalid/y.git"],
        check=True,
    )
    assert_deny(run_hook("pre_tool_use", bash(cmd), repo))


def test_remote_add_is_denied_when_remote_listing_fails(run_hook: Any, tmp_path: Path):
    nongit = tmp_path / "not-a-repo"
    nongit.mkdir()
    res = run_hook(
        "pre_tool_use", bash("git remote add origin https://example.invalid/x.git"), nongit
    )
    assert_deny(res)


def test_deny_reason_explains_alternative(repo: Path, run_hook: Any):
    res = run_hook("pre_tool_use", bash("git reset --hard"), repo)
    assert_deny(res)
    assert "reset --hard" in reason(res)
    assert "stash" in reason(res)


def test_multiple_reasons_are_combined(repo: Path, run_hook: Any):
    res = run_hook("pre_tool_use", bash(f"git reset --hard && rm {TEST_FILE}"), repo)
    assert_deny(res)
    r = reason(res)
    assert "reset --hard" in r
    assert TEST_FILE in r


# ------------------------------------------------------ (b) test deletion: Bash
TEST_DELETE_CMDS = [
    f"rm {TEST_FILE}",
    f"rm -f {TEST_FILE}",
    "rm -rf tests",
    "rm -rf tests/",
    "rm -r tests/unit",
    "rm tests/unit/*.cpp",
    f"rm ./{TEST_FILE}",
    "cd tests && rm unit/test_a.cpp",
    f"git rm {TEST_FILE}",
    "git rm -r tests/unit",
    f"git rm --cached {TEST_FILE}",
    f"unlink {TEST_FILE}",
    f"truncate -s 0 {TEST_FILE}",
    f"echo > {TEST_FILE}",
    f": > {TEST_FILE}",
    "find tests -name '*.cpp' -delete",
    f"bash -c 'rm {TEST_FILE}'",
]


@pytest.mark.parametrize("cmd", TEST_DELETE_CMDS)
def test_bash_test_deletion_is_denied(repo: Path, run_hook: Any, cmd: str):
    res = run_hook("pre_tool_use", bash(cmd), repo)
    assert_deny(res)
    assert "test" in reason(res).lower()


def test_bash_test_deletion_absolute_path(repo: Path, run_hook: Any):
    assert_deny(run_hook("pre_tool_use", bash(f"rm {repo / TEST_FILE}"), repo))


def test_rm_empty_test_dir_and_append_are_allowed(repo: Path, run_hook: Any):
    (repo / "tests" / "empty").mkdir()
    assert_silent(run_hook("pre_tool_use", bash("rmdir tests/empty"), repo))
    assert_silent(run_hook("pre_tool_use", bash("rm -r tests/empty"), repo))
    assert_silent(run_hook("pre_tool_use", bash(f"echo '// x' >> {TEST_FILE}"), repo))
    assert_silent(run_hook("pre_tool_use", bash("cat tests/unit/test_a.cpp"), repo))


def test_bash_test_deletion_allowed_by_reduction_exception(repo: Path, run_hook: Any):
    write_task(
        repo,
        {
            "id": "T-RED",
            "phase": "GREEN",
            "files": {"impl": [], "tests": ["tests/unit/**"]},
            "events": [RED_EVENT, GREEN_EVENT],
        },
        allow_test_reduction=True,
        allow_reason="duplicate test",
    )
    assert_silent(run_hook("pre_tool_use", bash(f"git rm {TEST_FILE}"), repo))


# ------------------------------------------------- (b) test deletion: file tools
@pytest.mark.parametrize("content", ["", "   \n\t\n"])
def test_write_emptying_test_file_is_denied(repo: Path, run_hook: Any, content: str):
    payload = file_tool("Write", file_path=str(repo / TEST_FILE), content=content)
    assert_deny(run_hook("pre_tool_use", payload, repo))


@pytest.mark.parametrize(
    "content",
    [
        'TEST_CASE("a") {\n  REQUIRE(a() == 1);\n}\n',  # fewer assertions
        "REQUIRE(a() == 1);\nREQUIRE(a() != 2);\n",  # fewer test cases
    ],
)
def test_write_reducing_test_file_is_denied(repo: Path, run_hook: Any, content: str):
    payload = file_tool("Write", file_path=str(repo / TEST_FILE), content=content)
    assert_deny(run_hook("pre_tool_use", payload, repo))


@pytest.mark.parametrize(
    "content",
    [
        'TEST_CASE("a") {\n  CHECK(a() == 1);\n  REQUIRE(a() != 2);\n  REQUIRE(a() != 3);\n}\n',
        'TEST_CASE("a") {\n  CHECK(a() == 1);\n  CHECK(a() != 2);\n}\nTEST_CASE("b") {}\n',
    ],
)
def test_write_keeping_or_adding_tests_is_allowed(repo: Path, run_hook: Any, content: str):
    payload = file_tool("Write", file_path=str(repo / TEST_FILE), content=content)
    assert_silent(run_hook("pre_tool_use", payload, repo))


def test_write_new_test_file_and_non_test_file(repo: Path, run_hook: Any):
    payload = file_tool("Write", file_path=str(repo / "tests/unit/test_new.cpp"), content="")
    assert_silent(run_hook("pre_tool_use", payload, repo))
    (repo / "docs" / "x.md").write_text("REQUIRE(x)\n")
    payload = file_tool("Write", file_path=str(repo / "docs/x.md"), content="")
    assert_silent(run_hook("pre_tool_use", payload, repo))


def test_edit_removing_require_is_denied(repo: Path, run_hook: Any):
    payload = file_tool(
        "Edit",
        file_path=str(repo / TEST_FILE),
        old_string="  REQUIRE(a() != 2);\n",
        new_string="",
        replace_all=False,
    )
    assert_deny(run_hook("pre_tool_use", payload, repo))


def test_edit_require_to_check_is_allowed(repo: Path, run_hook: Any):
    payload = file_tool(
        "Edit",
        file_path=str(repo / TEST_FILE),
        old_string="REQUIRE(",
        new_string="CHECK(",
        replace_all=True,
    )
    assert_silent(run_hook("pre_tool_use", payload, repo))


def test_edit_replace_all_counts_every_occurrence(repo: Path, run_hook: Any):
    # The single occurrence would keep the count; replace_all removes both.
    payload = file_tool(
        "Edit",
        file_path=str(repo / TEST_FILE),
        old_string="REQUIRE(",
        new_string="(void)(",
        replace_all=True,
    )
    assert_deny(run_hook("pre_tool_use", payload, repo))


def test_edit_unreadable_file_compares_strings(repo: Path, run_hook: Any):
    path = repo / "tests/unit/test_missing.cpp"
    deny = file_tool("Edit", file_path=str(path), old_string="REQUIRE(x);", new_string="")
    assert_deny(run_hook("pre_tool_use", deny, repo))
    ok = file_tool("Edit", file_path=str(path), old_string="REQUIRE(x);", new_string="CHECK(x);")
    assert_silent(run_hook("pre_tool_use", ok, repo))


def test_multiedit_applies_edits_in_order(repo: Path, run_hook: Any):
    path = str(repo / TEST_FILE)
    deny = file_tool(
        "MultiEdit",
        file_path=path,
        edits=[
            {"old_string": 'TEST_CASE("a")', "new_string": 'TEST_CASE("b")'},
            {"old_string": 'TEST_CASE("b") {', "new_string": "void f() {"},
        ],
    )
    assert_deny(run_hook("pre_tool_use", deny, repo))
    ok = file_tool(
        "MultiEdit",
        file_path=path,
        edits=[
            {"old_string": 'TEST_CASE("a")', "new_string": 'TEST_CASE("renamed")'},
            {"old_string": "REQUIRE(", "new_string": "CHECK(", "replace_all": True},
        ],
    )
    assert_silent(run_hook("pre_tool_use", ok, repo))


def test_python_test_reduction_is_denied(repo: Path, run_hook: Any):
    path = repo / "tests/unit/test_x.py"
    path.write_text("def test_a():\n    assert 1\n\n\ndef test_b():\n    assert 2\n")
    deny = file_tool(
        "Edit",
        file_path=str(path),
        old_string="\n\ndef test_b():\n    assert 2\n",
        new_string="\n",
    )
    assert_deny(run_hook("pre_tool_use", deny, repo))
    ok = file_tool("Edit", file_path=str(path), old_string="assert 1", new_string="assert 1 == 1")
    assert_silent(run_hook("pre_tool_use", ok, repo))


def test_static_require_counted_once(repo: Path, run_hook: Any):
    path = repo / "tests/unit/test_s.cpp"
    path.write_text('TEST_CASE("s") {\n  STATIC_REQUIRE(true);\n  REQUIRE(1);\n}\n')
    ok = file_tool("Edit", file_path=str(path), old_string="STATIC_REQUIRE(", new_string="CHECK(")
    assert_silent(run_hook("pre_tool_use", ok, repo))


def test_notebook_delete_in_tests_is_denied(repo: Path, run_hook: Any):
    nb = str(repo / "tests/unit/nb.ipynb")
    assert_deny(
        run_hook(
            "pre_tool_use", file_tool("NotebookEdit", notebook_path=nb, edit_mode="delete"), repo
        )
    )
    replace = file_tool("NotebookEdit", notebook_path=nb, edit_mode="replace", new_source="x = 1")
    assert_silent(run_hook("pre_tool_use", replace, repo))
    docs_nb = str(repo / "docs/nb.ipynb")
    assert_silent(
        run_hook(
            "pre_tool_use",
            file_tool("NotebookEdit", notebook_path=docs_nb, edit_mode="delete"),
            repo,
        )
    )


def _reduction_task(repo: Path, **overrides: Any) -> None:
    task: dict[str, Any] = {
        "id": "T-REDUCE",
        "phase": "GREEN",
        "files": {"impl": [], "tests": ["tests/unit/**"]},
        "events": [RED_EVENT, GREEN_EVENT],
    }
    active_fields = overrides.pop("active_fields", {"allow_test_reduction": True})
    task.update(overrides)
    write_task(repo, task, **active_fields, allow_reason="merged duplicate cases")


REDUCE_EDIT = {"old_string": "  REQUIRE(a() != 2);\n", "new_string": ""}


def test_reduction_exception_allows_edit(repo: Path, run_hook: Any):
    _reduction_task(repo)
    payload = file_tool("Edit", file_path=str(repo / TEST_FILE), **REDUCE_EDIT)
    assert_silent(run_hook("pre_tool_use", payload, repo))


@pytest.mark.parametrize(
    "overrides",
    [
        {"active_fields": {}},  # flag not set
        {"active_fields": {"allow_test_reduction": False}},
        {"events": [RED_EVENT]},  # no GREEN evidence
        {"files": {"impl": [], "tests": ["tests/other/**"]}},  # glob does not match
        {"phase": "DONE"},  # finished task is not in progress
    ],
)
def test_reduction_exception_requires_all_conditions(
    repo: Path, run_hook: Any, overrides: dict[str, Any]
):
    _reduction_task(repo, **overrides)
    payload = file_tool("Edit", file_path=str(repo / TEST_FILE), **REDUCE_EDIT)
    assert_deny(run_hook("pre_tool_use", payload, repo))


# ---------------------------------------------------------- (c) TDD enforcement
IMPL_EDIT = {"old_string": "return 1;", "new_string": "return 2;"}


def impl_edit(repo: Path, rel: str = "src/core/a.cpp") -> dict[str, Any]:
    return file_tool("Edit", file_path=str(repo / rel), **IMPL_EDIT)


def test_tdd_warns_when_not_enforced(repo: Path, run_hook: Any):
    res = run_hook("pre_tool_use", impl_edit(repo), repo)
    assert res.code == 0
    assert res.decision is None
    assert "tdd.py start" in res.context
    assert "tdd.py red" in res.context
    assert res.json["hookSpecificOutput"]["hookEventName"] == "PreToolUse"


def test_tdd_denies_when_enforced(repo: Path, run_hook: Any):
    set_enforce(repo, True)
    res = run_hook("pre_tool_use", impl_edit(repo), repo)
    assert_deny(res)
    assert "tdd.py start" in reason(res)


@pytest.mark.parametrize(
    "tool_payload",
    [
        lambda r: file_tool(
            "Write", file_path=str(r / "include/sssp/x.hpp"), content="#pragma once\n"
        ),
        lambda r: file_tool("Write", file_path=str(r / "tools/gen.cpp"), content="int main(){}\n"),
        lambda r: file_tool("MultiEdit", file_path=str(r / "src/core/a.cpp"), edits=[IMPL_EDIT]),
        lambda r: file_tool(
            "NotebookEdit",
            notebook_path=str(r / "src/nb.ipynb"),
            edit_mode="insert",
            new_source="1",
        ),
    ],
)
def test_tdd_covers_all_impl_dirs_and_tools(repo: Path, run_hook: Any, tool_payload: Any):
    set_enforce(repo, True)
    assert_deny(run_hook("pre_tool_use", tool_payload(repo), repo))


def test_tdd_allows_matching_task_with_red(repo: Path, run_hook: Any):
    set_enforce(repo, True)
    write_task(
        repo,
        {
            "id": "T1",
            "files": {"impl": ["src/core/**"], "tests": ["tests/unit/**"]},
            "events": [RED_EVENT],
        },
    )
    assert_silent(run_hook("pre_tool_use", impl_edit(repo), repo))


@pytest.mark.parametrize(
    ("task", "active"),
    [
        ({"files": {"impl": ["src/other/**"], "tests": []}, "events": [RED_EVENT]}, True),
        ({"files": {"impl": ["src/core/**"], "tests": []}, "events": []}, True),
        (
            {
                "files": {"impl": ["src/core/**"], "tests": []},
                "events": [RED_EVENT, GREEN_EVENT],
                "phase": "DONE",
            },
            True,
        ),
        ({"files": {"impl": ["src/core/**"], "tests": []}, "events": [RED_EVENT]}, False),
    ],
)
def test_tdd_rejects_unsuitable_tasks(
    repo: Path, run_hook: Any, task: dict[str, Any], active: bool
):
    set_enforce(repo, True)
    write_task(repo, {"id": "T2", **task}, active=active)
    assert_deny(run_hook("pre_tool_use", impl_edit(repo), repo))


@pytest.mark.parametrize(
    "rel",
    [TEST_FILE, "tests/unit/test_new.cpp", "docs/x.md", "CMakeLists.txt", "scripts/x.py"],
)
def test_tdd_ignores_non_impl_paths(repo: Path, run_hook: Any, rel: str):
    set_enforce(repo, True)
    content = (repo / TEST_FILE).read_text() + "// more\n" if rel == TEST_FILE else "x\n"
    payload = file_tool("Write", file_path=str(repo / rel), content=content)
    assert_silent(run_hook("pre_tool_use", payload, repo))


def test_tdd_ignores_files_outside_repo(repo: Path, run_hook: Any, tmp_path: Path):
    set_enforce(repo, True)
    payload = file_tool(
        "Write", file_path=str(tmp_path / "elsewhere" / "src" / "x.cpp"), content=""
    )
    assert_silent(run_hook("pre_tool_use", payload, repo))


# ------------------------------------------------------------------ robustness
@pytest.mark.parametrize("raw", ["", "   ", "not json", "[1, 2]", "{", "null"])
def test_invalid_stdin_is_silent(repo: Path, run_hook: Any, raw: str):
    assert_silent(run_hook("pre_tool_use", {}, repo, raw_stdin=raw))


@pytest.mark.parametrize(
    "payload",
    [
        {"tool_name": "Read", "tool_input": {"file_path": "/etc/hosts"}},
        {"tool_name": "Bash", "tool_input": {}},
        {"tool_name": "Bash", "tool_input": "git reset --hard"},
        {"tool_name": "Edit"},
        {"tool_name": "Write", "tool_input": {"content": "x"}},
    ],
)
def test_irrelevant_or_malformed_payload_is_silent(repo: Path, run_hook: Any, payload: Any):
    assert_silent(run_hook("pre_tool_use", {"hook_event_name": "PreToolUse", **payload}, repo))


def test_never_answers_allow_and_is_fast(repo: Path, run_hook: Any):
    start = time.monotonic()
    res = run_hook("pre_tool_use", bash("git status"), repo)
    elapsed = time.monotonic() - start
    assert_silent(res)
    assert elapsed < 2.0
    for cmd in ("git push -f", "ls"):
        out = run_hook("pre_tool_use", bash(cmd), repo)
        assert '"allow"' not in out.stdout
