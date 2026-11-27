# 開発ハーネス設計（Claude Code Hooks と TDD 状態管理）

この文書は、`.claude/settings.json` に登録した5つの hook（`.claude/hooks/*.py`）と、TDD の状態・証跡を管理する CLI `scripts/tdd.py` の設計、根拠にした公式仕様、動作確認の手順をまとめる。

- 対象の Claude Code: 2.1.332
- 仕様の出典（2026-11-27 に確認）:
  - Hooks reference: <https://code.claude.com/docs/en/hooks>
  - Hooks guide（トラブルシューティング）: <https://code.claude.com/docs/en/hooks-guide>
  - Permissions（deny ルールの書式）: <https://code.claude.com/docs/en/permissions>

## 1. 構成

| ファイル | 役割 |
| :- | :- |
| `.claude/settings.json` | hooks の登録と、危険な git/gh 操作の `permissions.deny`（コミット対象） |
| `.claude/hooks/_common.py` | 共通処理（stdin 読み込み、JSON 出力、ルート探索、設定、タスク読み込み、GLOB 照合、タイムアウト付き実行） |
| `.claude/hooks/session_start.py` | SessionStart: 状況の要約を `additionalContext` で渡す |
| `.claude/hooks/pre_tool_use.py` | PreToolUse: 危険操作・テスト削除・TDD 違反を拒否または警告 |
| `.claude/hooks/post_tool_use.py` | PostToolUse: 対象ファイルだけ整形し、該当カテゴリのテストだけをビルドして実行 |
| `.claude/hooks/pre_compact.py` | PreCompact: 状態の要約を `docs/PROGRESS.md` の自動生成ブロックへ保存 |
| `.claude/hooks/stop_gate.py` | Stop: 失敗テスト・RED 放置・未コミット変更・PROGRESS 未更新があれば継続を求める |
| `.tdd/config.json` | ハーネスの設定（コミット対象。`enforce: false` で開始） |
| `.tdd/runtime/` | 実行時の状態（gitignore 済み） |
| `scripts/tdd.py` | TDD の状態と証跡の CLI |
| `tests/harness/` | 上のすべてに対する pytest（`harness_util.py` に共有ヘルパー、`conftest.py` に fixture） |

hook はすべて **標準ライブラリだけで書いた Python 3.9 互換のスクリプト**で、`uv` が無くても `python3` で直接動く。

## 2. 公式仕様の要点（この実装が依存している部分）

### 2.1 設定と起動

- 置き場所: プロジェクトで共有する hook は `.claude/settings.json`。`.claude/settings.local.json`（gitignore 対象）や `~/.claude/settings.json` の hook とは**置き換えではなく足し合わされる**（"Hook entries merge across settings levels"）。同じハンドラが複数の設定ファイルにあっても1回しか実行されない。
- 3段の入れ子: イベント → matcher グループ → ハンドラ。
- matcher: 英数字・`_`・`-`・空白・`,`・`|` だけなら**完全一致のリスト**（`Edit|Write` は Edit か Write）。それ以外の文字を含むと JavaScript の正規表現として部分一致になる（`Edit.*` は `NotebookEdit` にも一致する）。本設定は完全一致リストだけを使う。
- 主なイベントの matcher 対象: PreToolUse/PostToolUse はツール名、SessionStart は `startup|resume|clear|compact|fork`、PreCompact は `manual|auto`。Stop は matcher 非対応（書いても無視される）。
- `${CLAUDE_PROJECT_DIR}` はセッション開始時のプロジェクトルート。パスを含む hook は **exec 形式（`command` + `args`）を推奨**（シェルを介さないので引用符が要らない）。exec 形式でも `CLAUDE_PROJECT_DIR` は環境変数として渡される。
- worktree に入っても `${CLAUDE_PROJECT_DIR}` は元のルートのまま。作業ディレクトリは入力 JSON の `cwd` が追従する。このため hook は「編集対象ファイルの git ルート → `cwd` の git ルート → `$CLAUDE_PROJECT_DIR`」の順でルートを決める（`_common.project_root`）。
- 既定のタイムアウトは command hook で 600 秒。`timeout`（秒）で個別に変更できる。
- hook は subagent の中でも動く（ツール系イベントの入力に `agent_id` と `agent_type` が付く）。

### 2.2 入出力

- 入力: stdin に JSON。共通フィールドは `session_id`, `transcript_path`, `cwd`, `permission_mode`, `hook_event_name`。Write/Edit/Read の `tool_input.file_path` は常に絶対パス。
- exit 0: stdout が **JSON オブジェクト1つだけ**なら構造化された結果として読まれる。前後に別の文字があるとプレーンテキスト扱い。SessionStart のプレーンテキストは Claude の文脈に入り、PreCompact のプレーンテキストは要約への指示になる。
- exit 2: ブロック（stderr が理由）。ブロックできるかはイベント次第。PreToolUse はツール呼び出しを止め、Stop は停止を止め（会話が続く）、**PreCompact はコンパクションを止める**。SessionStart と PostToolUse はブロックできない。
- その他の終了コード: ブロックしないエラー（処理は続く）。
- `additionalContext` は `hookSpecificOutput` の中に `hookEventName` と一緒に書く。1万字を超えるとファイルに逃がされプレビューしか渡らないため短くする。指示調より事実を述べる文体のほうがよい（プロンプトインジェクション対策に引っかかりにくい）。
- 判定の書き方:
  - PreToolUse: `hookSpecificOutput.permissionDecision`（`allow`/`deny`/`ask`/`defer`）と `permissionDecisionReason`。`deny` は bypassPermissions モードでも効く。`allow` は許可プロンプトを飛ばすので、**本ハーネスは `allow` を一切返さない**（問題が無ければ無出力で、通常の許可フローに任せる）。トップレベルの `decision`/`reason` は PreToolUse では非推奨。
  - PostToolUse / Stop / PreCompact: トップレベルの `decision: "block"` と `reason`。
  - Stop: `decision: "block"` の `reason` は必須で、Claude が続ける理由になる。
- Stop の入力には `stop_hook_active`（Stop hook による継続中なら true）、`last_assistant_message`、`background_tasks`、`session_crons` が付く。Claude Code 側にも「連続8回の継続で強制終了」の上限があり（`CLAUDE_CODE_STOP_HOOK_BLOCK_CAP` で変更可）、この回数はツール呼び出しが入るとリセットされる。Stop は**ユーザーの割り込み（Esc）では発火しない**。
- `transcript_path` の JSONL は非同期に書かれ、最新ターンが遅れることがある。

### 2.3 有効化の条件

- 設定ファイルを直接編集すると、通常はファイル監視で**自動的に読み込まれる**（"Direct edits to hooks in settings files are normally picked up automatically by the file watcher"）。数秒たっても `/hooks` に出ないときは、監視が取りこぼしたのでセッションを再起動する。
- `/hooks` は読み取り専用の一覧で、承認の操作は無い。
- 対話セッションでは、**ワークスペース信頼ダイアログ**を承認するまで、どの設定ファイルの hook も実行されない。`-p` や SDK のセッションはダイアログを出さず、信頼済みとして扱う。
- 一時的に全 hook を止めるには `"disableAllHooks": true`、1回の実行だけなら `claude --settings '{"disableAllHooks": true}'`。
- 本リポジトリでは、2026-11-27 の作業中セッションで `.claude/settings.json` を書いた直後から PreToolUse が効くことを確認した（無害な `git checkout -- <存在しないファイル>` が hook により拒否された）。

## 3. 各 hook の設計

### 3.1 失敗時の方針

- 原則 fail-open: hook 内部の例外では作業を止めない（exit 0、診断は stderr）。
- 例外として、PreToolUse の**危険操作の検出**と**テスト削除の検出**は fail-closed。判定中に例外が起きたら `deny` し、理由に例外を書く。
- PreCompact は exit 2 も `decision: "block"` も使わない（コンパクションを止めないため）。
- 多層防御として、主要な危険操作は `permissions.deny` にも入れてある。hook が起動できなかった場合（python3 が無いなど）でも、こちらが効く。

### 3.2 SessionStart（`startup|resume|clear|compact`、timeout 20 秒）

`additionalContext` に40行以内で次を出す: `git status --short --branch` の要約、`docs/PROGRESS.md` の自動生成ブロックまたは冒頭、`docs/TASKS.md` の Status 列が TODO/IN_PROGRESS/BLOCKED の行（最大10件と総数）、作業中の TDD タスクと phase、cmake/ninja/uv/clang-format の有無、前回 CI。CI は `git remote` があり `gh` があるときだけ `gh run list --limit 1`（4秒で打ち切り）を呼び、無い・失敗・時間切れでも1行の注記で済ませる。実測は約0.2秒。

### 3.3 PreToolUse（`Edit|Write|MultiEdit|NotebookEdit|Bash`、timeout 15 秒）

1. **Bash の危険操作**（fail-closed）。コマンドを `shlex`（`punctuation_chars`）で字句解析し、`;` `&&` `||` `|` `&` 改行 `(` `)` とバッククォートでサブコマンドに分ける。先頭の `VAR=val` と `sudo env command nohup time nice timeout stdbuf xargs noglob exec` を剥がし、`bash|sh|zsh -c` と `eval` の中身は再帰的に調べる。引用符の中の文字列（`echo "git push --force"` や `git commit -m "reset --hard"`）は誤検出しない。字句解析に失敗したときは、生の文字列に保守的な正規表現をかける。
   - `git push` の `--force`、`-f`（`-uf` のような束も）、`--force-with-lease`、`--force-if-includes`、`+refspec`、`--mirror`
   - `git reset --hard`
   - 未コミット変更の破棄: `git checkout` の `.`・`-- <path>`・`-f`、`git restore`（`--staged` だけの場合を除く）
   - `git clean -f`（`-fd`・`-xdf` などの束も）、`git branch -D`（`--delete --force` も）
   - `git remote set-url/remove/rm/rename/add`。ただし remote が1つも無いときの `git remote add origin <url>` だけは許可（一覧を取れなければ拒否）
   - `gh repo edit --visibility`、`gh repo delete`
   - `rm` の再帰・強制で、展開後に `$` が残るパス、リポジトリ外、ルートそのもの、`/`、`.git` を含むパス。ただし入力の `scratchpad_dir`、`/tmp`、`/private/tmp`、`$TMPDIR` の配下は許可
2. **テスト削除**（fail-closed）。対象は `tests/` 配下。
   - Bash の `rm`、`git rm`、`unlink`、`truncate`、リダイレクトでの上書き（`> tests/...`）
   - Write で中身を空にする、Write/Edit/MultiEdit の適用後に assertion（`REQUIRE*`/`CHECK*`/`STATIC_REQUIRE`/`assert`）か test case（`TEST_CASE*`/`TEMPLATE_TEST_CASE*`/`SCENARIO*`/`def test_*`）の数が減る、NotebookEdit の delete
   - 許可されるのは、作業中タスクに `tdd.py allow-test-change <ID> --reason ..` で明示許可があり、そのタスクに GREEN があり、対象パスがそのタスクの test GLOB に一致するときだけ。許可は次の `tdd.py red` で取り消される。
3. **TDD 強制**（fail-open）。`src/` `include/` `tools/` のファイルを Edit/Write するとき、impl GLOB が一致する作業中タスクに RED 証跡が必要。`tests/` は RED 前でも編集できる（テストを先に書くため）。ドキュメント等は対象外。`.tdd/config.json` の `enforce` が `false` なら拒否せず `additionalContext` で警告し、`true` なら `deny` する。

### 3.4 PostToolUse（`Edit|Write|MultiEdit`、timeout 200 秒）

- 整形は対象ファイルだけ: C++ は clang-format（`TDD_CLANG_FORMAT` → `/opt/homebrew/opt/llvm/bin/clang-format` → PATH の順）で `-i --style=file`、Python は ruff（`TDD_RUFF` → `.venv/bin/ruff` → PATH → `uv run ruff`）で `format` と `check`。
- 軽量チェックは C++ の変更時だけで、`build/<preset>/build.ninja` が既にある（configure 済み）ときに限る。カテゴリは `tests/<cat>/...` ならその cat、`src/` `include/` `tools/` なら `unit`。`cmake --build build/<preset> --target sssp_<cat>_tests` のあと `ctest --test-dir build/<preset> -L ^<cat>$ -LE bench --output-on-failure --timeout 60` を実行する。**ベンチマークのターゲット、`all`、ベンチマーク用スクリプトは一切実行しない。**
- 全体のタイムアウトは `TDD_POSTCHECK_TIMEOUT`（秒）> `post_check.timeout_sec`（既定120）。並列の担当と build ディレクトリを共有するので `.tdd/runtime/postcheck.lock` を `flock` で排他し、10秒（`TDD_POSTCHECK_LOCK_WAIT`）待っても取れなければ skipped にする。`TDD_POSTCHECK=0` でビルドとテストを止められる。
- 結果は `.tdd/runtime/last_check.json` に保存し、失敗・時間切れ・整形エラー・ruff の指摘があるときだけ `additionalContext` で知らせる。重いチェックは `make verify` に任せる。

```json
{"updated": "...", "categories": {"unit": {"status": "pass|fail|timeout|skipped", "timestamp": "...",
  "file": "tests/unit/x.cpp", "command": "...", "exit_code": 0, "output_tail": "...", "reason": "..."}},
 "last": {"category": "unit", "...": "..."}}
```

### 3.5 PreCompact（matcher なし、timeout 30 秒）

完了タスク、作業中タスクと phase、未解決の問題（last_check の fail/timeout、`.tdd/runtime/stop_gate_unresolved.md` の有無、TASKS.md の BLOCKED 行）、タスクごとの次のコマンド、最後の commit、`docs/PAPER_MAPPING.md`（無ければ `PAPER_MAPPING.ja.md`）の先頭10行を Markdown にまとめる。`docs/PROGRESS.md` の `<!-- AUTO:precompact -->`〜`<!-- /AUTO -->` の間だけを書き換え（ブロック外の手書き部分は変えない。ブロックが無ければ末尾に追記、ファイルが無ければ作成）、同じ内容を `.tdd/runtime/precompact_snapshot.md` にも保存する。stdout には要約への短い指示（プレーンテキスト）を出す。

### 3.6 Stop（matcher なし、timeout 30 秒）

次のどれかに当たれば `{"decision":"block","reason":...}` で継続を求める（各項目は `.tdd/config.json` の `stop_gate` で個別に止められる）。

1. `last_check.json` に status `fail` のカテゴリがある
2. 作業中タスクの最終イベントが RED のまま
3. `git status --porcelain` が空でない（`.tdd/runtime/` は除く）
4. `docs/PROGRESS.md` が無い、または未コミットの変更ファイルより古い

3 と 4 は、入力の `background_tasks` に実行中の subagent / teammate / workflow があるときは block せず、`systemMessage` で警告するだけにする（その担当の作業途中の変更かもしれず、リードには解消できないため）。

無限ループの防止:

- `stop_hook_active` が true のときは block しない。
- 同じ理由での block が連続して上限（`max_consecutive_blocks`、既定3）に達したら block しない。理由の同一性は**問題の種類で判定する**（未コミット変更なら件数やファイル名が変わっても同じ理由）。block しなかった停止は回数に数えない。
- 上の2つで停止を許可したときは、時刻・理由・発生条件・再現手順を `.tdd/runtime/stop_gate_unresolved.md` に追記する（同じ理由は連続して書かない。パスは相対で書き、利用者のホームなどの絶対パスは残さない）。この記録はセッション ID を含む実行時の生成物なので、git で無視される `.tdd/runtime/` に置き、追跡対象のパスには書かない（テストで確認している）。利用者には `systemMessage` で短く知らせる。
- 状態は `.tdd/runtime/stop_gate.json` に保存する。

利用者の明示的な停止の尊重（仕様の範囲で取れる方法の検討結果）:

- いちばん確実なのは **Esc での割り込み**。公式仕様で、割り込みでは Stop hook がそもそも発火しない。
- Stop の入力には利用者の発言が含まれないので、`transcript_path` の JSONL を末尾から読み、最後の「人間の発言」（tool_result だけの行や `<system-reminder>` などで始まる行は除く）を取り出す。200字以下で「止めて」「止まって」「やめて」「中断」「ストップ」「終了して」「stop」「halt」「abort」を含むときは block しない。長い指示文の中の単語では反応しない。transcript が遅れて最新の発言を拾えない可能性は残る（そのときも最大で1回 block するだけで、次の停止は `stop_hook_active` で通る）。
- 恒久的に止めるには `.tdd/config.json` の `stop_gate.enabled` を false、1回だけなら環境変数 `TDD_STOP_GATE=off`。
- 検討して採らなかった案: UserPromptSubmit hook で発言を `.tdd/runtime` に記録する方法。確実だが毎回の発言に hook が1本増え、今回の担当範囲（5本の hook）を超えるので見送った。transcript の遅れが実際に問題になったら移行する。

## 4. `scripts/tdd.py`

```sh
uv run python scripts/tdd.py start <ID> --spec "..." [--paper "Sec 3.2"] --impl 'src/x/*.cpp' --test 'tests/unit/x/*.cpp'
uv run python scripts/tdd.py red <ID> --reason "なぜ失敗するか" [--tests "case"] -- ctest --preset debug -R x
uv run python scripts/tdd.py green <ID> -- ctest --preset debug -R x
uv run python scripts/tdd.py refactor <ID> -- ctest --preset debug
uv run python scripts/tdd.py allow-test-change <ID> --reason "重複したケースを統合"   # GREEN 後のみ
uv run python scripts/tdd.py done <ID> [--commit SHA]
uv run python scripts/tdd.py status | check [--no-git] | log
```

- `start`: `docs/tdd/tasks/<ID>.json` を作り、`.tdd/runtime/active.json` に登録する（複数のタスクを同時に作業中にできる）。既にあるタスクなら証跡を残したまま再登録する。ID は英数字と `.` `_` `-` だけ。
- `red`: コマンドを実行し、終了コードが0以外のときだけ RED を記録する（0なら拒否してエラー）。`--reason` は必須。
- `green` / `refactor`: 終了コード0のときだけ記録する。`green` は RED の後、`refactor` は GREEN の後でなければ拒否する。
- 引数が1つ（`-- "a && b"`）ならシェルで、複数ならシェルなしで実行する。出力はそのまま画面に流し、末尾20行を `output_tail` に保存する。保存時にリポジトリのルートを `.`、一時ディレクトリを `$TMPDIR`、ホームを `~`、`pytest-of-<ユーザー名>` を `pytest-of-$USER` に置き換え、利用者固有のパスを証跡に残さない。全文は `.tdd/runtime/runs/<ID>-<時刻>-<PHASE>.log`、最新の1件は `.tdd/runtime/last_run.json`。
- テストだけのタスク（既存の振る舞いを検証するだけで、正直な RED が作れないもの）は `start <ID> --spec .. --test <glob> [--covers <glob>] --test-only "<理由>"` で作り、`green` から記録する。書式と検査は `docs/tdd/README.md` の「Test-only tasks」。
- `check`: すべての証跡を検証する（必須キー、id とファイル名の一致、phase の値、RED の終了コードが0以外で理由がある、GREEN/REFACTOR が0、RED が GREEN より先、GREEN が REFACTOR より先、最初のイベントが RED、DONE の最後が GREEN か REFACTOR、phase と最後のイベントの一致、タイムスタンプが ISO 8601（`+0900` と `+09:00` の両方を受け付ける））。1件でも不正なら終了コード1。
  - さらに、`commit` 欄が空でない証跡について、その SHA がリポジトリに存在し、HEAD から到達できること（`git merge-base --is-ancestor <sha> HEAD`）を確かめる。存在しない SHA と、HEAD から届かない SHA（履歴の書き直しや、マージされていない枝）は不正。短縮 SHA も受け付ける。
  - git が無い、リポジトリでない、浅い clone（`git rev-parse --is-shallow-repository` が `true`）、HEAD が無い、のいずれかなら判定できないので、`WARN commit fields not checked ...` を1行出すだけで不正にはしない。`--no-git` を付けるとこの検査を飛ばす（`NOTE` を1行出す）。
  - CI でこの検査を効かせるには、`actions/checkout` に `fetch-depth: 0` が要る（既定の深さ1では警告になるだけ）。
- `log`: すべての証跡から `docs/TDD_LOG.md` を作り直す（内容が同じなら出力も同じ）。
- `active.json` の読み書きは `.tdd/runtime/tdd.lock` で排他する。

## 5. 動作確認の手順

```sh
# テストと lint（hook はすべて subprocess に stdin JSON を渡して検証している）
uv run pytest tests/harness
uv run ruff check . && uv run ruff format --check .

# 個別に手で動かす（exit code と stdout の JSON を見る）
echo '{"hook_event_name":"PreToolUse","tool_name":"Bash","tool_input":{"command":"git push -f"},"cwd":"'"$PWD"'"}' \
  | python3 .claude/hooks/pre_tool_use.py; echo "exit=$?"
echo '{"hook_event_name":"SessionStart","source":"startup","cwd":"'"$PWD"'"}' | python3 .claude/hooks/session_start.py
echo '{"hook_event_name":"Stop","stop_hook_active":false,"cwd":"'"$PWD"'"}' | python3 .claude/hooks/stop_gate.py

# Python 3.9 でも動くことの確認（macOS 付属の python3）
echo '{}' | /usr/bin/python3 .claude/hooks/pre_tool_use.py; echo "exit=$?"
```

- セッション内での確認: `/hooks` で5つのイベントに project settings 由来の hook が出ていること。`claude --debug-file <path>` で起動すると、`Hooks: PreToolUse:Bash [...] finished with status 0 (31ms)` のような行で各 hook の結果と所要時間が分かる。
- PreCompact は実ファイル（`docs/PROGRESS.md`）を書くので、手で試すときは作業ツリーを一時ディレクトリにコピーしてから実行する。

## 6. 既知の限界

- PreToolUse は Bash の字句解析に頼っているので、`xargs rm -rf < list`、`find ... -delete`、変数に入れたコマンドの実行、スクリプトファイル経由の実行など、間接的な実行は検出できない。`git stash drop` も対象外。
- PostToolUse はツール経由の編集だけが対象。Bash で書き換えたファイルには反応しない（公式仕様どおり）。
- PostToolUse の結果はカテゴリごとの最新1件で、別の担当の作業途中の変更でビルドが失敗した場合もそのカテゴリの fail として残る。直ったかどうかは、そのカテゴリの次の編集か `make verify` で確かめる。
- **作業中の担当の未コミット変更（対処済み、2026-11-27）。** Stop hook はメインのセッションで動き、`git status --porcelain` はメインの作業ツリー全体を見る。そのため以前は、メインの作業ツリーで直接編集している担当の作業途中の変更も「未コミット変更」と「PROGRESS.md が古い」に数えられ、リードが自分ではコミットできない変更で止められていた（`.claude/worktrees/` は gitignore されているので、worktree の担当の分はもともと数えない）。いまは、Stop の入力の `background_tasks` に実行中（status が completed/failed/killed などの終了状態でない）の `subagent` / `teammate` / `workflow` が1つでもある間、この2項目では block せず `systemMessage` で利用者に知らせるだけにしている。Stop の `additionalContext` は会話を継続させてしまうので、警告には使わない。失敗したテストと RED のままのタスクは、担当が動いていても block する。警告だけで止めた停止は連続回数に数えない。残る限界と次の改善案:
  - 担当が動いている間は、リード自身の未コミット変更も警告に下がる。ファイル単位で分けるには、PostToolUse で入力に `agent_id` がある編集（担当による編集）のファイルを `.tdd/runtime/edit_owners.json` に記録し、Stop で実行中の担当が編集したファイルだけを除く方法がある（Bash での書き換えは記録できない）。
  - `background_tasks` は公式に「タスク一覧が参照できるとき」だけ付く。付かない環境では従来どおり block する（`stop_gate.check_uncommitted` を false にすれば止められる）。
  - 運用でも、担当は原則 `isolation: "worktree"` で起動すると衝突しにくい。
