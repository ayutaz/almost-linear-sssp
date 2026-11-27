# Development harness design (Claude Code hooks and TDD state management)

**English** | [日本語](harness.ja.md) | [简体中文](harness.zh-CN.md)

This document describes the design of the five hooks registered in `.claude/settings.json` (`.claude/hooks/*.py`) and of the CLI `scripts/tdd.py`, which manages TDD state and evidence, together with the official specifications they rely on and the verification procedure.

- Target Claude Code version: 2.1.332
- Sources of the specifications (checked on 2026-11-27):
  - Hooks reference: <https://code.claude.com/docs/en/hooks>
  - Hooks guide (troubleshooting): <https://code.claude.com/docs/en/hooks-guide>
  - Permissions (format of deny rules): <https://code.claude.com/docs/en/permissions>

## 1. Structure

| File | Role |
| :- | :- |
| `.claude/settings.json` | Registration of the hooks and `permissions.deny` for dangerous git/gh operations (committed) |
| `.claude/hooks/_common.py` | Shared code (reading stdin, JSON output, root discovery, configuration, task loading, GLOB matching, execution with a timeout) |
| `.claude/hooks/session_start.py` | SessionStart: passes a summary of the situation through `additionalContext` |
| `.claude/hooks/pre_tool_use.py` | PreToolUse: denies or warns about dangerous operations, test deletion and TDD violations |
| `.claude/hooks/post_tool_use.py` | PostToolUse: formats only the target file, and builds and runs only the tests of the relevant category |
| `.claude/hooks/pre_compact.py` | PreCompact: saves a summary of the state into the auto-generated block of `docs/PROGRESS.md` |
| `.claude/hooks/stop_gate.py` | Stop: asks to continue if there are failing tests, a lingering RED, uncommitted changes or a stale PROGRESS |
| `.tdd/config.json` | Harness configuration (committed; starts with `enforce: false`) |
| `.tdd/runtime/` | Runtime state (gitignored) |
| `scripts/tdd.py` | CLI for TDD state and evidence |
| `tests/harness/` | pytest for all of the above (shared helpers in `harness_util.py`, fixtures in `conftest.py`) |

Every hook is a **Python 3.9 compatible script written with the standard library only**, so it runs directly with `python3` even without `uv`.

## 2. Key points of the official specification (what this implementation depends on)

### 2.1 Configuration and launch

- Location: hooks shared by the project go in `.claude/settings.json`. They are **added to, not replaced by,** the hooks in `.claude/settings.local.json` (gitignored) and `~/.claude/settings.json` ("Hook entries merge across settings levels"). The same handler runs only once even if it appears in several settings files.
- Three nested levels: event -> matcher group -> handler.
- matcher: if it consists only of alphanumerics, `_`, `-`, spaces, `,` and `|`, it is an **exact-match list** (`Edit|Write` means Edit or Write). If it contains any other character, it becomes a JavaScript regular expression with partial matching (`Edit.*` also matches `NotebookEdit`). This configuration uses exact-match lists only.
- What the matcher applies to for the main events: the tool name for PreToolUse/PostToolUse, `startup|resume|clear|compact|fork` for SessionStart, and `manual|auto` for PreCompact. Stop does not support a matcher (it is ignored if written).
- `${CLAUDE_PROJECT_DIR}` is the project root at session start. For hooks that contain paths, the **exec form (`command` + `args`) is recommended** (it does not go through a shell, so no quoting is needed). Even in the exec form, `CLAUDE_PROJECT_DIR` is passed as an environment variable.
- Even after entering a worktree, `${CLAUDE_PROJECT_DIR}` stays the original root. The working directory follows the `cwd` of the input JSON. For this reason a hook decides the root in the order "git root of the edited file -> git root of `cwd` -> `$CLAUDE_PROJECT_DIR`" (`_common.project_root`).
- The default timeout is 600 seconds for a command hook. It can be changed individually with `timeout` (seconds).
- Hooks also run inside subagents (`agent_id` and `agent_type` are attached to the input of tool events).

### 2.2 Input and output

- Input: JSON on stdin. The common fields are `session_id`, `transcript_path`, `cwd`, `permission_mode` and `hook_event_name`. `tool_input.file_path` of Write/Edit/Read is always an absolute path.
- exit 0: if stdout is **exactly one JSON object**, it is read as a structured result. If there are other characters before or after it, it is treated as plain text. Plain text from SessionStart goes into Claude's context, and plain text from PreCompact becomes an instruction for the summary.
- exit 2: block (stderr is the reason). Whether blocking is possible depends on the event. PreToolUse stops the tool call, Stop stops the stop (the conversation continues), and **PreCompact stops the compaction**. SessionStart and PostToolUse cannot block.
- Other exit codes: a non-blocking error (processing continues).
- `additionalContext` is written inside `hookSpecificOutput` together with `hookEventName`. Above 10,000 characters it is moved to a file and only a preview is passed, so keep it short. A style that states facts works better than an instruction-like style (it is less likely to be caught by prompt-injection defences).
- How to write decisions:
  - PreToolUse: `hookSpecificOutput.permissionDecision` (`allow`/`deny`/`ask`/`defer`) and `permissionDecisionReason`. `deny` is effective even in bypassPermissions mode. `allow` skips the permission prompt, so **this harness never returns `allow`** (when there is no problem it emits nothing and leaves it to the normal permission flow). The top-level `decision`/`reason` is deprecated for PreToolUse.
  - PostToolUse / Stop / PreCompact: top-level `decision: "block"` and `reason`.
  - Stop: the `reason` of `decision: "block"` is required and becomes the reason Claude continues.
- The input of Stop carries `stop_hook_active` (true while continuing because of a Stop hook), `last_assistant_message`, `background_tasks` and `session_crons`. Claude Code itself also has a cap of "forced end after 8 consecutive continuations" (changeable with `CLAUDE_CODE_STOP_HOOK_BLOCK_CAP`), and this count is reset when a tool call occurs. Stop **does not fire on a user interrupt (Esc)**.
- The JSONL at `transcript_path` is written asynchronously, so the latest turn may lag.

### 2.3 Conditions for activation

- When a settings file is edited directly, it is normally **loaded automatically** by the file watcher ("Direct edits to hooks in settings files are normally picked up automatically by the file watcher"). If it does not show up in `/hooks` after a few seconds, the watcher missed it, so restart the session.
- `/hooks` is a read-only list, and there is no approval operation.
- In an interactive session, no hook from any settings file runs until the **workspace trust dialog** is accepted. Sessions with `-p` or the SDK do not show the dialog and treat the workspace as trusted.
- To stop all hooks temporarily, use `"disableAllHooks": true`, or for a single run only, `claude --settings '{"disableAllHooks": true}'`.
- In this repository, we confirmed in the working session on 2026-11-27 that PreToolUse takes effect right after `.claude/settings.json` was written (a harmless `git checkout -- <nonexistent file>` was denied by the hook).

## 3. Design of each hook

### 3.1 Policy on failure

- Fail-open in principle: an exception inside a hook does not stop the work (exit 0, diagnostics on stderr).
- As an exception, the **dangerous-operation detection** and the **test-deletion detection** of PreToolUse are fail-closed. If an exception occurs during the judgement, it returns `deny` and writes the exception in the reason.
- PreCompact uses neither exit 2 nor `decision: "block"` (so as not to stop the compaction).
- As defence in depth, the main dangerous operations are also in `permissions.deny`. Even when a hook cannot start (for example, python3 is missing), that rule still works.

### 3.2 SessionStart (`startup|resume|clear|compact`, timeout 20 seconds)

It emits the following in `additionalContext` within 40 lines: a summary of `git status --short --branch`, the auto-generated block or the beginning of `docs/PROGRESS.md`, the rows of `docs/TASKS.md` whose Status column is TODO/IN_PROGRESS/BLOCKED (at most 10 plus the total count), the active TDD task and phase, the presence of cmake/ninja/uv/clang-format, and the previous CI. For CI, it calls `gh run list --limit 1` (cut off after 4 seconds) only when there is a `git remote` and `gh` is available, and if it is absent, fails or times out, it is dealt with by a one-line note. The measured time is about 0.2 seconds.

### 3.3 PreToolUse (`Edit|Write|MultiEdit|NotebookEdit|Bash`, timeout 15 seconds)

1. **Dangerous Bash operations** (fail-closed). The command is tokenised with `shlex` (`punctuation_chars`) and split into subcommands at `;` `&&` `||` `|` `&` newline `(` `)` and backticks. A leading `VAR=val` and `sudo env command nohup time nice timeout stdbuf xargs noglob exec` are stripped, and the contents of `bash|sh|zsh -c` and `eval` are inspected recursively. Strings inside quotes (`echo "git push --force"` or `git commit -m "reset --hard"`) are not falsely detected. When tokenisation fails, a conservative regular expression is applied to the raw string.
   - `git push` with `--force`, `-f` (including bundles such as `-uf`), `--force-with-lease`, `--force-if-includes`, `+refspec`, `--mirror`
   - `git reset --hard`
   - Discarding uncommitted changes: `git checkout` with `.`, `-- <path>` or `-f`, and `git restore` (except when it is only `--staged`)
   - `git clean -f` (including bundles such as `-fd` and `-xdf`), `git branch -D` (including `--delete --force`)
   - `git remote set-url/remove/rm/rename/add`. The only exception is `git remote add origin <url>` when there is no remote at all (denied if the list cannot be obtained)
   - `gh repo edit --visibility`, `gh repo delete`
   - A recursive or forced `rm` whose expanded path still contains `$`, is outside the repository, is the root itself, is `/`, or contains `.git`. Paths under the `scratchpad_dir` of the input, `/tmp`, `/private/tmp` and `$TMPDIR` are allowed
2. **Test deletion** (fail-closed). The target is under `tests/`.
   - In Bash: `rm`, `git rm`, `unlink`, `truncate`, and overwriting by redirection (`> tests/...`)
   - Emptying the content with Write; a decrease, after Write/Edit/MultiEdit is applied, in the number of assertions (`REQUIRE*`/`CHECK*`/`STATIC_REQUIRE`/`assert`) or test cases (`TEST_CASE*`/`TEMPLATE_TEST_CASE*`/`SCENARIO*`/`def test_*`); a delete by NotebookEdit
   - It is allowed only when the active task has an explicit permission via `tdd.py allow-test-change <ID> --reason ..`, that task has a GREEN, and the target path matches the test GLOB of that task. The permission is revoked by the next `tdd.py red`.
3. **TDD enforcement** (fail-open). When Edit/Write targets a file under `src/`, `include/` or `tools/`, an active task whose impl GLOB matches must have RED evidence. `tests/` can be edited even before RED (because tests are written first). Documents and the like are out of scope. If `enforce` in `.tdd/config.json` is `false`, it does not deny but warns through `additionalContext`, and if it is `true`, it returns `deny`.

### 3.4 PostToolUse (`Edit|Write|MultiEdit`, timeout 200 seconds)

- Formatting applies only to the target file: C++ uses clang-format (in the order `TDD_CLANG_FORMAT` -> `/opt/homebrew/opt/llvm/bin/clang-format` -> PATH) with `-i --style=file`, and Python uses ruff (in the order `TDD_RUFF` -> `.venv/bin/ruff` -> PATH -> `uv run ruff`) with `format` and `check`.
- The lightweight check runs only when C++ changed and `build/<preset>/build.ninja` already exists (already configured). The category is the cat for `tests/<cat>/...`, and `unit` for `src/`, `include/` and `tools/`. It runs `cmake --build build/<preset> --target sssp_<cat>_tests` and then `ctest --test-dir build/<preset> -L ^<cat>$ -LE bench --output-on-failure --timeout 60`. **It never runs benchmark targets, `all`, or benchmark scripts.**
- The overall timeout is `TDD_POSTCHECK_TIMEOUT` (seconds) > `post_check.timeout_sec` (default 120). Because the build directory is shared with parallel agents, `.tdd/runtime/postcheck.lock` is held exclusively with `flock`, and if it cannot be taken within 10 seconds (`TDD_POSTCHECK_LOCK_WAIT`), the check is marked skipped. `TDD_POSTCHECK=0` turns off the build and tests.
- The result is saved to `.tdd/runtime/last_check.json`, and is reported through `additionalContext` only on a failure, a timeout, a formatting error or a ruff finding. The heavy checks are left to `make verify`.

```json
{"updated": "...", "categories": {"unit": {"status": "pass|fail|timeout|skipped", "timestamp": "...",
  "file": "tests/unit/x.cpp", "command": "...", "exit_code": 0, "output_tail": "...", "reason": "..."}},
 "last": {"category": "unit", "...": "..."}}
```

### 3.5 PreCompact (no matcher, timeout 30 seconds)

It puts the following into Markdown: completed tasks, the active task and phase, unresolved problems (a fail/timeout in last_check, the presence of `.tdd/runtime/stop_gate_unresolved.md`, BLOCKED rows of TASKS.md), the next command per task, the last commit, and the first 10 lines of `docs/PAPER_MAPPING.md` (or `PAPER_MAPPING.ja.md` if it does not exist). It rewrites only the part between `<!-- AUTO:precompact -->` and `<!-- /AUTO -->` in `docs/PROGRESS.md` (the hand-written part outside the block is not changed; if there is no block it appends to the end, and if there is no file it creates one), and saves the same content to `.tdd/runtime/precompact_snapshot.md`. It prints a short instruction for the summary (plain text) to stdout.

### 3.6 Stop (no matcher, timeout 30 seconds)

If any of the following applies, it asks to continue with `{"decision":"block","reason":...}` (each item can be turned off individually through `stop_gate` in `.tdd/config.json`).

1. `last_check.json` has a category with status `fail`
2. The last event of an active task is still RED
3. `git status --porcelain` is not empty (excluding `.tdd/runtime/`)
4. `docs/PROGRESS.md` is missing or older than the uncommitted changed files

For 3 and 4, when the `background_tasks` of the input contains a running subagent / teammate / workflow, it does not block and only warns with `systemMessage` (the changes may be that agent's work in progress, and the lead cannot resolve them).

Prevention of infinite loops:

- When `stop_hook_active` is true, it does not block.
- When blocks for the same reason in a row reach the limit (`max_consecutive_blocks`, default 3), it does not block. The sameness of the reason is **judged by the kind of problem** (for uncommitted changes it is the same reason even if the count or file names change). A stop that was not blocked is not counted.
- When a stop is allowed by the two rules above, it appends the time, the reason, the triggering condition and the reproduction steps to `.tdd/runtime/stop_gate_unresolved.md` (the same reason is not written consecutively; paths are written as relative paths, and absolute paths such as the user's home are not left). This record contains the session ID and is a runtime artifact, so it is placed in `.tdd/runtime/`, which git ignores, and is not written to a tracked path (confirmed by a test). The user is notified briefly with `systemMessage`.
- The state is saved to `.tdd/runtime/stop_gate.json`.

Respecting an explicit stop by the user (the result of considering what can be done within the specification):

- The most reliable way is an **interrupt with Esc**. By the official specification, the Stop hook does not fire at all on an interrupt.
- The input of Stop does not contain the user's utterance, so it reads the JSONL at `transcript_path` from the end and extracts the last "human utterance" (excluding lines that consist only of a tool_result or start with `<system-reminder>` and the like). If it is 200 characters or fewer and contains one of the Japanese stop words (such as 止めて, 止まって, やめて, 中断, ストップ, 終了して) or `stop`, `halt`, `abort`, it does not block. It does not react to a word inside a long instruction. The possibility remains that the transcript lags and the latest utterance is not picked up (even then it blocks at most once, and the next stop passes through `stop_hook_active`).
- To stop it permanently, set `stop_gate.enabled` in `.tdd/config.json` to false; for one time only, use the environment variable `TDD_STOP_GATE=off`.
- An option that was considered and not adopted: recording the utterance into `.tdd/runtime` with a UserPromptSubmit hook. It is reliable, but it adds one more hook to every utterance and goes beyond the scope of this assignment (five hooks), so it was set aside. We will move to it if the transcript lag becomes a real problem.

## 4. `scripts/tdd.py`

```sh
uv run python scripts/tdd.py start <ID> --spec "..." [--paper "Sec 3.2"] --impl 'src/x/*.cpp' --test 'tests/unit/x/*.cpp'
uv run python scripts/tdd.py red <ID> --reason "why it fails" [--tests "case"] -- ctest --preset debug -R x
uv run python scripts/tdd.py green <ID> -- ctest --preset debug -R x
uv run python scripts/tdd.py refactor <ID> -- ctest --preset debug
uv run python scripts/tdd.py allow-test-change <ID> --reason "merge duplicated cases"   # only after GREEN
uv run python scripts/tdd.py done <ID> [--commit SHA]
uv run python scripts/tdd.py status | check [--no-git] | log
```

- `start`: creates `docs/tdd/tasks/<ID>.json` and registers it in `.tdd/runtime/active.json` (several tasks can be active at the same time). For an existing task, it registers it again while keeping the evidence. An ID may contain only alphanumerics and `.` `_` `-`.
- `red`: runs the command and records RED only when the exit code is non-zero (if it is 0, it refuses with an error). `--reason` is required.
- `green` / `refactor`: record only when the exit code is 0. `green` is refused unless it comes after RED, and `refactor` unless it comes after GREEN.
- With one argument (`-- "a && b"`) it runs through a shell, and with several arguments it runs without a shell. The output is passed through to the screen as is, and the last 20 lines are saved in `output_tail`. On saving, the repository root is replaced with `.`, the temporary directory with `$TMPDIR`, the home with `~`, and `pytest-of-<user name>` with `pytest-of-$USER`, so that no user-specific path stays in the evidence. The full text goes to `.tdd/runtime/runs/<ID>-<time>-<PHASE>.log`, and the latest one to `.tdd/runtime/last_run.json`.
- A test-only task (one that only verifies existing behaviour and cannot make an honest RED) is created with `start <ID> --spec .. --test <glob> [--covers <glob>] --test-only "<reason>"` and records from `green`. The format and checks are in "Test-only tasks" of `docs/tdd/README.md`.
- `check`: validates all the evidence (required keys, id matching the file name, the value of phase, a RED exit code that is non-zero with a reason, GREEN/REFACTOR being 0, RED before GREEN, GREEN before REFACTOR, the first event being RED, the last of DONE being GREEN or REFACTOR, phase matching the last event, timestamps in ISO 8601 (both `+0900` and `+09:00` are accepted)). If even one is invalid, the exit code is 1.
  - In addition, for evidence whose `commit` field is non-empty, it verifies that the SHA exists in the repository and is reachable from HEAD (`git merge-base --is-ancestor <sha> HEAD`). A nonexistent SHA and a SHA not reachable from HEAD (rewritten history or an unmerged branch) are invalid. Abbreviated SHAs are also accepted.
  - If git is missing, it is not a repository, it is a shallow clone (`git rev-parse --is-shallow-repository` is `true`), or there is no HEAD, it cannot judge, so it only prints one line `WARN commit fields not checked ...` and does not make it invalid. With `--no-git` this check is skipped (it prints one `NOTE` line).
  - To make this check effective in CI, `actions/checkout` needs `fetch-depth: 0` (with the default depth of 1 it only warns).
- `log`: regenerates `docs/TDD_LOG.md` from all the evidence (the same content gives the same output).
- Reading and writing `active.json` is made exclusive with `.tdd/runtime/tdd.lock`.

## 5. Verification procedure

```sh
# Tests and lint (every hook is verified by passing stdin JSON to a subprocess)
uv run pytest tests/harness
uv run ruff check . && uv run ruff format --check .

# Run one by hand (look at the exit code and the JSON on stdout)
echo '{"hook_event_name":"PreToolUse","tool_name":"Bash","tool_input":{"command":"git push -f"},"cwd":"'"$PWD"'"}' \
  | python3 .claude/hooks/pre_tool_use.py; echo "exit=$?"
echo '{"hook_event_name":"SessionStart","source":"startup","cwd":"'"$PWD"'"}' | python3 .claude/hooks/session_start.py
echo '{"hook_event_name":"Stop","stop_hook_active":false,"cwd":"'"$PWD"'"}' | python3 .claude/hooks/stop_gate.py

# Check that it also works with Python 3.9 (the python3 bundled with macOS)
echo '{}' | /usr/bin/python3 .claude/hooks/pre_tool_use.py; echo "exit=$?"
```

- Checks inside a session: in `/hooks`, hooks derived from the project settings should appear for the five events. If you start with `claude --debug-file <path>`, lines such as `Hooks: PreToolUse:Bash [...] finished with status 0 (31ms)` show the result and duration of each hook.
- PreCompact writes a real file (`docs/PROGRESS.md`), so when trying it by hand, copy the working tree to a temporary directory first and run it there.

## 6. Known limitations

- PreToolUse relies on lexical analysis of Bash, so it cannot detect indirect execution such as `xargs rm -rf < list`, `find ... -delete`, running a command stored in a variable, or running through a script file. `git stash drop` is also out of scope.
- PostToolUse covers only edits made through tools. It does not react to files rewritten by Bash (as in the official specification).
- The result of PostToolUse is the latest one per category, and even when the build fails because of another agent's work in progress, it remains as a fail of that category. Whether it was fixed is confirmed by the next edit in that category or by `make verify`.
- **Uncommitted changes of a running agent (handled, 2026-11-27).** The Stop hook runs in the main session, and `git status --porcelain` looks at the whole working tree of the main session. So previously, the work-in-progress changes of an agent editing directly in the main working tree were also counted as "uncommitted changes" and "PROGRESS.md is stale", and the lead was stopped by changes it could not commit itself (`.claude/worktrees/` is gitignored, so the part of agents in worktrees was not counted to begin with). Now, as long as `background_tasks` in the input of Stop has even one running `subagent` / `teammate` / `workflow` (one whose status is not a terminal state such as completed/failed/killed), these two items do not block and only inform the user with `systemMessage`. The `additionalContext` of Stop would make the conversation continue, so it is not used for the warning. Failing tests and a task left at RED block even while agents are running. A stop that was only warned is not counted in the consecutive count. The remaining limitations and ideas for the next improvement:
  - While an agent is running, the lead's own uncommitted changes also drop down to a warning. To separate them per file, there is a method: in PostToolUse, record the files of edits whose input has an `agent_id` (edits by an agent) in `.tdd/runtime/edit_owners.json`, and in Stop exclude only the files edited by running agents (rewrites by Bash cannot be recorded).
  - `background_tasks` is officially attached only "when the task list can be referenced". In environments where it is not attached, it blocks as before (it can be stopped by setting `stop_gate.check_uncommitted` to false).
  - Operationally too, starting agents with `isolation: "worktree"` in principle makes collisions less likely.
