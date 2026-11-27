# 开发 harness 设计(Claude Code Hooks 与 TDD 状态管理)

[English](harness.md) | [日本語](harness.ja.md) | **简体中文**

本文档汇总了注册在 `.claude/settings.json` 中的 5 个 hook(`.claude/hooks/*.py`)和管理 TDD 状态与证据的 CLI `scripts/tdd.py` 的设计、它们所依据的官方规范,以及验证步骤。

- 目标 Claude Code 版本:2.1.332
- 规范的出处(2026-11-27 确认):
  - Hooks reference:<https://code.claude.com/docs/en/hooks>
  - Hooks guide(故障排查):<https://code.claude.com/docs/en/hooks-guide>
  - Permissions(deny 规则的写法):<https://code.claude.com/docs/en/permissions>

## 1. 结构

| 文件 | 作用 |
| :- | :- |
| `.claude/settings.json` | hook 的注册,以及针对危险 git/gh 操作的 `permissions.deny`(纳入提交) |
| `.claude/hooks/_common.py` | 公共处理(读取 stdin、输出 JSON、查找根目录、配置、读取任务、GLOB 匹配、带超时的执行) |
| `.claude/hooks/session_start.py` | SessionStart:通过 `additionalContext` 传递状况摘要 |
| `.claude/hooks/pre_tool_use.py` | PreToolUse:拒绝或警告危险操作、测试删除和 TDD 违规 |
| `.claude/hooks/post_tool_use.py` | PostToolUse:只格式化目标文件,只构建并运行相应类别的测试 |
| `.claude/hooks/pre_compact.py` | PreCompact:把状态摘要保存到 `docs/PROGRESS.md` 的自动生成块中 |
| `.claude/hooks/stop_gate.py` | Stop:如果存在失败的测试、遗留的 RED、未提交的变更或未更新的 PROGRESS,就要求继续 |
| `.tdd/config.json` | harness 的配置(纳入提交。以 `enforce: false` 开始) |
| `.tdd/runtime/` | 运行时状态(已 gitignore) |
| `scripts/tdd.py` | TDD 状态与证据的 CLI |
| `tests/harness/` | 针对以上所有内容的 pytest(共享辅助函数在 `harness_util.py`,fixture 在 `conftest.py`) |

hook 全部是**只使用标准库、兼容 Python 3.9 的脚本**,即使没有 `uv`,也能直接用 `python3` 运行。

## 2. 官方规范的要点(本实现所依赖的部分)

### 2.1 配置与启动

- 位置:项目共享的 hook 放在 `.claude/settings.json`。它们与 `.claude/settings.local.json`(gitignore 对象)和 `~/.claude/settings.json` 中的 hook **不是替换而是叠加**("Hook entries merge across settings levels")。即使同一个处理程序出现在多个设置文件中,也只会运行一次。
- 三层嵌套:事件 -> matcher 组 -> 处理程序。
- matcher:如果只由字母数字、`_`、`-`、空白、`,`、`|` 组成,则是**完全匹配的列表**(`Edit|Write` 表示 Edit 或 Write)。如果含有其他字符,则成为 JavaScript 的正则表达式并按部分匹配处理(`Edit.*` 也会匹配 `NotebookEdit`)。本配置只使用完全匹配的列表。
- 主要事件的 matcher 对象:PreToolUse/PostToolUse 为工具名,SessionStart 为 `startup|resume|clear|compact|fork`,PreCompact 为 `manual|auto`。Stop 不支持 matcher(写了也会被忽略)。
- `${CLAUDE_PROJECT_DIR}` 是会话开始时的项目根目录。对于含路径的 hook,**推荐使用 exec 形式(`command` + `args`)**(不经过 shell,所以不需要引号)。即使是 exec 形式,`CLAUDE_PROJECT_DIR` 也会作为环境变量传入。
- 即使进入 worktree,`${CLAUDE_PROJECT_DIR}` 仍是原来的根目录。工作目录则跟随输入 JSON 的 `cwd`。因此 hook 按"被编辑文件的 git 根 -> `cwd` 的 git 根 -> `$CLAUDE_PROJECT_DIR`"的顺序决定根目录(`_common.project_root`)。
- command hook 的默认超时为 600 秒。可以用 `timeout`(秒)单独修改。
- hook 在 subagent 内部也会运行(工具类事件的输入中会附带 `agent_id` 和 `agent_type`)。

### 2.2 输入输出

- 输入:stdin 上的 JSON。公共字段为 `session_id`、`transcript_path`、`cwd`、`permission_mode`、`hook_event_name`。Write/Edit/Read 的 `tool_input.file_path` 始终是绝对路径。
- exit 0:如果 stdout **恰好是一个 JSON 对象**,就作为结构化结果读取。如果前后有其他字符,则按纯文本处理。SessionStart 的纯文本会进入 Claude 的上下文,PreCompact 的纯文本会成为对摘要的指示。
- exit 2:阻止(stderr 为理由)。能否阻止取决于事件。PreToolUse 会阻止工具调用,Stop 会阻止停止(对话继续),**PreCompact 会阻止压缩**。SessionStart 和 PostToolUse 无法阻止。
- 其他退出码:不阻止的错误(处理继续)。
- `additionalContext` 与 `hookEventName` 一起写在 `hookSpecificOutput` 中。超过 1 万字会被转存到文件,只传递预览,所以要写得简短。陈述事实的文体比指示式的文体更好(不容易被提示注入防御拦下)。
- 判定的写法:
  - PreToolUse:`hookSpecificOutput.permissionDecision`(`allow`/`deny`/`ask`/`defer`)和 `permissionDecisionReason`。`deny` 即使在 bypassPermissions 模式下也有效。`allow` 会跳过许可提示,所以**本 harness 绝不返回 `allow`**(没有问题时不输出任何内容,交给通常的许可流程)。顶层的 `decision`/`reason` 在 PreToolUse 中已被弃用。
  - PostToolUse / Stop / PreCompact:顶层的 `decision: "block"` 和 `reason`。
  - Stop:`decision: "block"` 的 `reason` 是必需的,并成为 Claude 继续的理由。
- Stop 的输入附带 `stop_hook_active`(由 Stop hook 导致继续时为 true)、`last_assistant_message`、`background_tasks`、`session_crons`。Claude Code 自身也有"连续继续 8 次就强制结束"的上限(可用 `CLAUDE_CODE_STOP_HOOK_BLOCK_CAP` 修改),这个次数在出现工具调用时会被重置。Stop **不会在用户中断(Esc)时触发**。
- `transcript_path` 的 JSONL 是异步写入的,最新的一轮可能会滞后。

### 2.3 启用的条件

- 直接编辑设置文件时,通常会被文件监视器**自动加载**("Direct edits to hooks in settings files are normally picked up automatically by the file watcher")。如果几秒后 `/hooks` 中仍未出现,说明监视器漏掉了,请重启会话。
- `/hooks` 是只读的列表,没有批准操作。
- 在交互式会话中,在**工作区信任对话框**被批准之前,任何设置文件中的 hook 都不会运行。`-p` 和 SDK 的会话不显示对话框,视为已信任。
- 要临时停止所有 hook,用 `"disableAllHooks": true`;只想停一次运行,用 `claude --settings '{"disableAllHooks": true}'`。
- 在本仓库中,我们在 2026-11-27 的工作会话里确认了:写入 `.claude/settings.json` 之后,PreToolUse 立即生效(一个无害的 `git checkout -- <不存在的文件>` 被 hook 拒绝了)。

## 3. 各个 hook 的设计

### 3.1 失败时的方针

- 原则上 fail-open:hook 内部的异常不会中止工作(exit 0,诊断信息输出到 stderr)。
- 作为例外,PreToolUse 的**危险操作检测**和**测试删除检测**是 fail-closed。判定过程中发生异常时返回 `deny`,并把异常写进理由。
- PreCompact 既不使用 exit 2,也不使用 `decision: "block"`(以免阻止压缩)。
- 作为纵深防御,主要的危险操作也写进了 `permissions.deny`。即使 hook 无法启动(例如没有 python3),这一侧仍然有效。

### 3.2 SessionStart(`startup|resume|clear|compact`,timeout 20 秒)

在 `additionalContext` 中以 40 行以内输出以下内容:`git status --short --branch` 的摘要、`docs/PROGRESS.md` 的自动生成块或开头部分、`docs/TASKS.md` 中 Status 列为 TODO/IN_PROGRESS/BLOCKED 的行(最多 10 条和总数)、进行中的 TDD 任务和 phase、cmake/ninja/uv/clang-format 是否存在,以及上一次的 CI。CI 仅在存在 `git remote` 且有 `gh` 时才调用 `gh run list --limit 1`(4 秒后中断),不存在、失败或超时的情况只用一行备注处理。实测约 0.2 秒。

### 3.3 PreToolUse(`Edit|Write|MultiEdit|NotebookEdit|Bash`,timeout 15 秒)

1. **Bash 的危险操作**(fail-closed)。用 `shlex`(`punctuation_chars`)对命令做词法分析,并在 `;` `&&` `||` `|` `&` 换行 `(` `)` 和反引号处拆分为子命令。去掉开头的 `VAR=val` 和 `sudo env command nohup time nice timeout stdbuf xargs noglob exec`,并递归检查 `bash|sh|zsh -c` 和 `eval` 的内容。引号内的字符串(`echo "git push --force"` 或 `git commit -m "reset --hard"`)不会被误检。词法分析失败时,对原始字符串应用保守的正则表达式。
   - `git push` 的 `--force`、`-f`(也包括 `-uf` 这样的组合)、`--force-with-lease`、`--force-if-includes`、`+refspec`、`--mirror`
   - `git reset --hard`
   - 丢弃未提交的变更:`git checkout` 的 `.`、`-- <path>`、`-f`,以及 `git restore`(只有 `--staged` 的情况除外)
   - `git clean -f`(也包括 `-fd`、`-xdf` 等组合)、`git branch -D`(也包括 `--delete --force`)
   - `git remote set-url/remove/rm/rename/add`。但只有在完全没有 remote 时的 `git remote add origin <url>` 是允许的(如果无法获取列表则拒绝)
   - `gh repo edit --visibility`、`gh repo delete`
   - 递归或强制的 `rm`,其展开后的路径仍含有 `$`、位于仓库之外、是根目录本身、是 `/`,或含有 `.git`。但输入中的 `scratchpad_dir`、`/tmp`、`/private/tmp`、`$TMPDIR` 之下是允许的
2. **测试删除**(fail-closed)。对象是 `tests/` 之下。
   - Bash 的 `rm`、`git rm`、`unlink`、`truncate`、用重定向覆盖(`> tests/...`)
   - 用 Write 清空内容;Write/Edit/MultiEdit 应用之后 assertion(`REQUIRE*`/`CHECK*`/`STATIC_REQUIRE`/`assert`)或 test case(`TEST_CASE*`/`TEMPLATE_TEST_CASE*`/`SCENARIO*`/`def test_*`)的数量减少;NotebookEdit 的 delete
   - 只有在进行中的任务通过 `tdd.py allow-test-change <ID> --reason ..` 获得明确许可、该任务已有 GREEN、且目标路径匹配该任务的 test GLOB 时才允许。该许可会在下一次 `tdd.py red` 时被撤销。
3. **TDD 强制**(fail-open)。用 Edit/Write 修改 `src/`、`include/`、`tools/` 下的文件时,impl GLOB 匹配的进行中的任务必须有 RED 证据。`tests/` 即使在 RED 之前也可以编辑(因为要先写测试)。文档等不在对象之内。如果 `.tdd/config.json` 的 `enforce` 为 `false`,则不拒绝,而是通过 `additionalContext` 警告;为 `true` 时则 `deny`。

### 3.4 PostToolUse(`Edit|Write|MultiEdit`,timeout 200 秒)

- 格式化只针对目标文件:C++ 用 clang-format(按 `TDD_CLANG_FORMAT` -> `/opt/homebrew/opt/llvm/bin/clang-format` -> PATH 的顺序)执行 `-i --style=file`,Python 用 ruff(按 `TDD_RUFF` -> `.venv/bin/ruff` -> PATH -> `uv run ruff` 的顺序)执行 `format` 和 `check`。
- 轻量检查只在 C++ 发生变更、且 `build/<preset>/build.ninja` 已经存在(已 configure)时运行。类别对于 `tests/<cat>/...` 是该 cat,对于 `src/`、`include/`、`tools/` 是 `unit`。执行 `cmake --build build/<preset> --target sssp_<cat>_tests`,然后执行 `ctest --test-dir build/<preset> -L ^<cat>$ -LE bench -j <jobs> --output-on-failure --timeout 60`。`<jobs>` 为 `TDD_POSTCHECK_JOBS` > `post_check.jobs`(默认 4)。不是正整数的值会被忽略,改用下一个来源(最后为 4)。**绝不运行基准测试的 target、`all` 和基准测试用的脚本。**
- 整体的超时为 `TDD_POSTCHECK_TIMEOUT`(秒)> `post_check.timeout_sec`(默认 120)。由于 build 目录与并行的 agent 共用,所以用 `flock` 对 `.tdd/runtime/postcheck.lock` 加排他锁,等待 10 秒(`TDD_POSTCHECK_LOCK_WAIT`)仍拿不到就标记为 skipped。`TDD_POSTCHECK=0` 可以关闭构建和测试。
- 结果保存到 `.tdd/runtime/last_check.json`,只有在失败、超时、格式化错误或 ruff 有指摘时才通过 `additionalContext` 通知。繁重的检查交给 `make verify`。

```json
{"updated": "...", "categories": {"unit": {"status": "pass|fail|timeout|skipped", "timestamp": "...",
  "file": "tests/unit/x.cpp", "command": "...", "exit_code": 0, "output_tail": "...", "reason": "..."}},
 "last": {"category": "unit", "...": "..."}}
```

### 3.5 PreCompact(无 matcher,timeout 30 秒)

把以下内容汇总成 Markdown:已完成的任务、进行中的任务和 phase、未解决的问题(last_check 中的 `fail`,以 `Failing check` 列出;是否存在 `.tdd/runtime/stop_gate_unresolved.md`、TASKS.md 中的 BLOCKED 行)、超时的检查(last_check 中的 `timeout`。超时不代表通过或失败,因此不放入 `Failing check`,而是列在专门的标题 `Inconclusive checks (timed out; not a test failure)` 之下;Stop gate 也不会因此阻止)、每个任务的下一条命令、最后一次 commit,以及 `docs/PAPER_MAPPING.md`(如果没有则用 `PAPER_MAPPING.ja.md`)的前 10 行。只改写 `docs/PROGRESS.md` 中 `<!-- AUTO:precompact -->` 到 `<!-- /AUTO -->` 之间的部分(块外手写的部分不会改动。如果没有块则追加到末尾,如果没有文件则创建),并把同样的内容也保存到 `.tdd/runtime/precompact_snapshot.md`。向 stdout 输出对摘要的简短指示(纯文本)。

### 3.6 Stop(无 matcher,timeout 30 秒)

只要符合下面任何一项,就用 `{"decision":"block","reason":...}` 要求继续(各项都可以通过 `.tdd/config.json` 的 `stop_gate` 单独关闭)。

1. `last_check.json` 中有 status 为 `fail` 的类别
2. 进行中任务的最后一个事件仍是 RED
3. `git status --porcelain` 不为空(`.tdd/runtime/` 除外)
4. `docs/PROGRESS.md` 不存在,或比未提交的变更文件旧

对于 3 和 4,当输入的 `background_tasks` 中有运行中的 subagent / teammate / workflow 时,不阻止,只用 `systemMessage` 警告(因为这可能是该 agent 工作到一半的变更,lead 无法解决)。

防止无限循环:

- `stop_hook_active` 为 true 时不阻止。
- 因相同理由连续阻止达到上限(`max_consecutive_blocks`,默认 3)时不阻止。理由是否相同**按问题的种类判定**(对于未提交的变更,即使数量或文件名变化也是相同的理由)。没有阻止的停止不计入次数。
- 在上述两条规则下允许停止时,把时间、理由、发生条件和复现步骤追加到 `.tdd/runtime/stop_gate_unresolved.md`(相同的理由不连续写入;路径写成相对路径,不留下用户主目录等绝对路径)。这份记录包含会话 ID,是运行时的生成物,所以放在被 git 忽略的 `.tdd/runtime/` 中,不写到被追踪的路径(已通过测试确认)。通过 `systemMessage` 简短地通知用户。
- 状态保存在 `.tdd/runtime/stop_gate.json`。

对用户明确要求停止的尊重(在规范范围内能采取的做法的研究结果):

- 最可靠的是**用 Esc 中断**。按官方规范,中断时 Stop hook 根本不会触发。
- Stop 的输入中不包含用户的发言,所以从 `transcript_path` 的 JSONL 末尾往前读,取出最后一条"人类的发言"(只含 tool_result 的行,以及以 `<system-reminder>` 等开头的行除外)。如果在 200 字以内,且包含日语的停止用语(如 止めて、止まって、やめて、中断、ストップ、終了して)或 `stop`、`halt`、`abort`,就不阻止。对长指示文中的单词不会反应。仍然存在 transcript 滞后而没能取到最新发言的可能性(即使那样,最多也只阻止一次,下一次停止会通过 `stop_hook_active` 放行)。
- 要永久停止,把 `.tdd/config.json` 的 `stop_gate.enabled` 设为 false;只停一次,则用环境变量 `TDD_STOP_GATE=off`。
- 研究过但没有采用的方案:用 UserPromptSubmit hook 把发言记录到 `.tdd/runtime`。这样虽然可靠,但每次发言都会多出一个 hook,超出了这次的负责范围(5 个 hook),所以搁置了。如果 transcript 的滞后真的成了问题,再迁移过去。

## 4. `scripts/tdd.py`

```sh
uv run python scripts/tdd.py start <ID> --spec "..." [--paper "Sec 3.2"] --impl 'src/x/*.cpp' --test 'tests/unit/x/*.cpp'
uv run python scripts/tdd.py red <ID> --reason "失败的原因" [--tests "case"] -- ctest --preset debug -R x
uv run python scripts/tdd.py green <ID> -- ctest --preset debug -R x
uv run python scripts/tdd.py refactor <ID> -- ctest --preset debug
uv run python scripts/tdd.py allow-test-change <ID> --reason "合并重复的用例"   # 仅限 GREEN 之后
uv run python scripts/tdd.py done <ID> [--commit SHA]
uv run python scripts/tdd.py status | check [--no-git] | log
```

- `start`:创建 `docs/tdd/tasks/<ID>.json`,并登记到 `.tdd/runtime/active.json`(可以同时进行多个任务)。对于已有的任务,在保留证据的情况下重新登记。ID 只能使用字母数字和 `.` `_` `-`。
- `red`:运行命令,只有退出码非 0 时才记录 RED(为 0 则拒绝并报错)。`--reason` 是必需的。
- `green` / `refactor`:只有退出码为 0 时才记录。`green` 必须在 RED 之后,`refactor` 必须在 GREEN 之后,否则拒绝。
- 只有一个参数(`-- "a && b"`)时通过 shell 运行,有多个参数时不经过 shell 运行。输出原样显示到屏幕,末尾 20 行保存到 `output_tail`。保存时,把仓库根目录替换为 `.`,临时目录替换为 `$TMPDIR`,主目录替换为 `~`,`pytest-of-<用户名>` 替换为 `pytest-of-$USER`,以免在证据中留下用户特有的路径。全文保存在 `.tdd/runtime/runs/<ID>-<时刻>-<PHASE>.log`,最新的一份保存在 `.tdd/runtime/last_run.json`。
- 仅含测试的任务(只验证既有行为,无法制造出诚实的 RED 的任务)用 `start <ID> --spec .. --test <glob> [--covers <glob>] --test-only "<理由>"` 创建,从 `green` 开始记录。格式和检查见 `docs/tdd/README.md` 的 "Test-only tasks"。
- `check`:验证所有证据(必需的键、id 与文件名一致、phase 的取值、RED 的退出码非 0 且有理由、GREEN/REFACTOR 为 0、RED 先于 GREEN、GREEN 先于 REFACTOR、第一个事件为 RED、DONE 的最后一个是 GREEN 或 REFACTOR、phase 与最后一个事件一致、时间戳为 ISO 8601(`+0900` 和 `+09:00` 两种都接受))。只要有一条不合法,退出码就是 1。
  - 此外,对于 `commit` 字段非空的证据,确认该 SHA 存在于仓库中,并且可以从 HEAD 到达(`git merge-base --is-ancestor <sha> HEAD`)。不存在的 SHA,以及无法从 HEAD 到达的 SHA(被改写的历史或未合并的分支)均不合法。缩写的 SHA 也接受。
  - 如果没有 git、不是仓库、是浅克隆(`git rev-parse --is-shallow-repository` 为 `true`)或没有 HEAD,就无法判定,所以只输出一行 `WARN commit fields not checked ...`,不算不合法。加上 `--no-git` 则跳过这项检查(输出一行 `NOTE`)。
  - 要在 CI 中让这项检查生效,`actions/checkout` 需要 `fetch-depth: 0`(默认深度 1 时只会警告)。
- `log`:根据所有证据重新生成 `docs/TDD_LOG.md`(内容相同则输出也相同)。
- `active.json` 的读写用 `.tdd/runtime/tdd.lock` 做互斥。

## 5. 验证步骤

```sh
# 测试和 lint(所有 hook 都通过向 subprocess 传入 stdin JSON 来验证)
uv run pytest tests/harness
uv run ruff check . && uv run ruff format --check .

# 手动单独运行(查看 exit code 和 stdout 的 JSON)
echo '{"hook_event_name":"PreToolUse","tool_name":"Bash","tool_input":{"command":"git push -f"},"cwd":"'"$PWD"'"}' \
  | python3 .claude/hooks/pre_tool_use.py; echo "exit=$?"
echo '{"hook_event_name":"SessionStart","source":"startup","cwd":"'"$PWD"'"}' | python3 .claude/hooks/session_start.py
echo '{"hook_event_name":"Stop","stop_hook_active":false,"cwd":"'"$PWD"'"}' | python3 .claude/hooks/stop_gate.py

# 确认在 Python 3.9 上也能运行(macOS 自带的 python3)
echo '{}' | /usr/bin/python3 .claude/hooks/pre_tool_use.py; echo "exit=$?"
```

- 会话内的确认:在 `/hooks` 中,5 个事件下应出现来自 project settings 的 hook。用 `claude --debug-file <path>` 启动后,可以通过 `Hooks: PreToolUse:Bash [...] finished with status 0 (31ms)` 这样的行看到每个 hook 的结果和耗时。
- PreCompact 会写入真实文件(`docs/PROGRESS.md`),所以手动试验时,请先把工作树复制到临时目录再运行。

## 6. 已知的局限

- PreToolUse 依赖对 Bash 的词法分析,所以无法检测间接的执行,例如 `xargs rm -rf < list`、`find ... -delete`、运行存放在变量中的命令、通过脚本文件运行。`git stash drop` 也不在对象之内。
- PostToolUse 只针对通过工具进行的编辑。对用 Bash 改写的文件不会有反应(与官方规范一致)。
- PostToolUse 的结果是每个类别最新的一条,即使因为别的 agent 工作到一半的变更导致构建失败,也会作为该类别的 fail 保留下来。是否已修复,要通过该类别的下一次编辑或 `make verify` 来确认。
- **运行中的 agent 的未提交变更(已处理,2026-11-27)。** Stop hook 在主会话中运行,`git status --porcelain` 查看的是主会话的整个工作树。因此以前,在主工作树中直接编辑的 agent 工作到一半的变更,也会被算作"未提交的变更"和"PROGRESS.md 过旧",lead 会被自己无法提交的变更拦住(`.claude/worktrees/` 已被 gitignore,所以 worktree 中的 agent 那一部分本来就不计入)。现在,只要 Stop 输入的 `background_tasks` 中哪怕有一个运行中(status 不是 completed/failed/killed 等终止状态)的 `subagent` / `teammate` / `workflow`,这两项就不阻止,只用 `systemMessage` 通知用户。Stop 的 `additionalContext` 会让对话继续,所以不用于警告。失败的测试和停在 RED 的任务,即使 agent 在运行也会阻止。只警告而停止的情况不计入连续次数。剩余的局限和下一步改进的方案:
  - 在 agent 运行期间,lead 自己的未提交变更也会降为警告。要按文件区分,可以采用这样的方法:在 PostToolUse 中,把输入带有 `agent_id` 的编辑(agent 的编辑)所涉及的文件记录到 `.tdd/runtime/edit_owners.json`,在 Stop 中只排除运行中 agent 编辑过的文件(用 Bash 改写的无法记录)。
  - `background_tasks` 按官方说法只在"可以引用任务列表时"才附带。在没有附带的环境中,仍然像以前一样阻止(把 `stop_gate.check_uncommitted` 设为 false 就可以停止它)。
  - 在运维上,原则上用 `isolation: "worktree"` 启动 agent 也更不容易冲突。
