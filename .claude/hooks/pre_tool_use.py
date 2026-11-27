#!/usr/bin/env python3
"""PreToolUse guard: dangerous Bash commands, test deletion, TDD enforcement.

Run by Claude Code as ``python3 .claude/hooks/pre_tool_use.py`` with the event
JSON on stdin (official spec: https://code.claude.com/docs/en/hooks).

Output contract:
- nothing on stdout + exit 0  -> no opinion (normal permission flow)
- {"hookSpecificOutput": {..., "permissionDecision": "deny", ...}} -> blocked
- {"hookSpecificOutput": {..., "additionalContext": ...}} -> warning only
"allow" is never returned (it would skip the user's permission prompt).

Checks:
(a) dangerous Bash (force push, reset --hard, discarding changes, git clean -f,
    branch -D, remote changes, gh repo visibility/delete, rm -rf outside the
    repository)                                         -> fail-closed
(b) deleting tests or reducing assertions / test cases  -> fail-closed
(c) editing src/ include/ tools/ without a RED-evidenced active task
                                                        -> fail-open
Python 3.9 compatible, standard library only.
"""

from __future__ import annotations

import glob as _glob
import os
import re
import sys
from pathlib import Path
from typing import Any

import _common as C

EVENT = "PreToolUse"
FILE_TOOLS = ("Write", "Edit", "MultiEdit", "NotebookEdit")
MAX_DEPTH = 8
SUBST = "$(...)"  # placeholder for a command substitution inside a word


# =========================================================================== lexer
class LexError(ValueError):
    pass


class Command:
    """One simple command: argv words (quotes removed) and redirections."""

    __slots__ = ("argv", "redirs")

    def __init__(self, argv: list[str], redirs: list[tuple[str, str]]) -> None:
        self.argv = argv
        self.redirs = redirs


class _Lexer:
    """Small POSIX-shell tokenizer.

    Splits on ``; && || | & newline ( )``, removes quotes, keeps single-quoted
    text literal, and lexes command substitutions (``$(..)``, backticks, also
    inside double quotes) and process substitutions as separate commands.
    Here-document bodies are skipped. Raises LexError on unterminated quotes.
    """

    def __init__(self, s: str, depth: int = 0) -> None:
        if depth > MAX_DEPTH:
            raise LexError("command substitution nested too deeply")
        self.s = s
        self.n = len(s)
        self.i = 0
        self.depth = depth
        self.cmds: list[Command] = []
        self.heredocs: list[tuple[str, bool]] = []

    def run(self) -> list[Command]:
        self.parse(in_subst=False)
        return self.cmds

    # -- helpers ---------------------------------------------------------------
    def _peek(self, k: int = 0) -> str:
        j = self.i + k
        return self.s[j] if j < self.n else ""

    def _skip_heredocs(self) -> None:
        for delim, strip_tabs in self.heredocs:
            while self.i < self.n:
                end = self.s.find("\n", self.i)
                if end < 0:
                    end = self.n
                line = self.s[self.i : end]
                self.i = end + 1
                if (line.lstrip("\t") if strip_tabs else line) == delim:
                    break
        self.heredocs = []

    def _read_heredoc_delim(self) -> str:
        while self._peek() in (" ", "\t"):
            self.i += 1
        out: list[str] = []
        while self.i < self.n:
            c = self.s[self.i]
            if c in " \t\n;&|<>()":
                break
            if c in "'\"":
                j = self.s.find(c, self.i + 1)
                if j < 0:
                    raise LexError("unterminated here-document delimiter")
                out.append(self.s[self.i + 1 : j])
                self.i = j + 1
                continue
            if c == "\\" and self.i + 1 < self.n:
                out.append(self.s[self.i + 1])
                self.i += 2
                continue
            out.append(c)
            self.i += 1
        return "".join(out)

    def _backtick(self) -> None:
        """self.i is on the opening backtick; lex its content recursively."""
        j = self.i + 1
        buf: list[str] = []
        while j < self.n and self.s[j] != "`":
            if self.s[j] == "\\" and j + 1 < self.n:
                buf.append(self.s[j + 1])
                j += 2
                continue
            buf.append(self.s[j])
            j += 1
        if j >= self.n:
            raise LexError("unterminated backtick")
        self.i = j + 1
        self.cmds.extend(_Lexer("".join(buf), self.depth + 1).run())

    def _subst(self) -> None:
        """self.i is just after ``$(`` / ``<(``: lex until the matching ``)``."""
        if self.depth + 1 > MAX_DEPTH:
            raise LexError("command substitution nested too deeply")
        self.depth += 1
        try:
            self.parse(in_subst=True)
        finally:
            self.depth -= 1

    def _ansi_c(self) -> str:
        """self.i is on the quote of ``$'..'``; returns the (roughly) decoded text."""
        j = self.i + 1
        buf: list[str] = []
        while j < self.n and self.s[j] != "'":
            if self.s[j] == "\\" and j + 1 < self.n:
                buf.append(self.s[j + 1])
                j += 2
                continue
            buf.append(self.s[j])
            j += 1
        if j >= self.n:
            raise LexError("unterminated $'...'")
        self.i = j + 1
        return "".join(buf)

    # -- main loop ---------------------------------------------------------------
    def parse(self, in_subst: bool) -> None:  # noqa: C901 - a tokenizer is a big switch
        s = self.s
        argv: list[str] = []
        redirs: list[tuple[str, str]] = []
        word: list[str] | None = None
        redir_op: str | None = None
        paren = 0

        def finish_word() -> None:
            nonlocal word, redir_op
            if word is None:
                return
            w = "".join(word)
            word = None
            if redir_op is not None:
                redirs.append((redir_op, w))
                redir_op = None
            else:
                argv.append(w)

        def finish_cmd() -> None:
            nonlocal argv, redirs, redir_op
            finish_word()
            if argv or redirs:
                self.cmds.append(Command(argv, redirs))
            argv, redirs, redir_op = [], [], None

        while self.i < self.n:
            c = s[self.i]
            if c in " \t\r":
                finish_word()
                self.i += 1
            elif c == "\n":
                finish_cmd()
                self.i += 1
                self._skip_heredocs()
            elif c == "\\":
                if self._peek(1) == "\n":
                    self.i += 2
                else:
                    if word is None:
                        word = []
                    word.append(self._peek(1))
                    self.i += 2
            elif c == "#" and word is None:
                end = s.find("\n", self.i)
                self.i = self.n if end < 0 else end
            elif c == "'":
                j = s.find("'", self.i + 1)
                if j < 0:
                    raise LexError("unterminated single quote")
                if word is None:
                    word = []
                word.append(s[self.i + 1 : j])
                self.i = j + 1
            elif c == '"':
                if word is None:
                    word = []
                self.i += 1
                while True:
                    if self.i >= self.n:
                        raise LexError("unterminated double quote")
                    d = s[self.i]
                    if d == '"':
                        self.i += 1
                        break
                    if d == "\\" and self._peek(1) in ('"', "\\", "$", "`", "\n"):
                        if self._peek(1) != "\n":
                            word.append(self._peek(1))
                        self.i += 2
                    elif d == "`":
                        self._backtick()
                        word.append(SUBST)
                    elif d == "$" and self._peek(1) == "(":
                        self.i += 2
                        self._subst()
                        word.append(SUBST)
                    else:
                        word.append(d)
                        self.i += 1
            elif c == "`":
                self._backtick()
                if word is None:
                    word = []
                word.append(SUBST)
            elif c == "$":
                if word is None:
                    word = []
                nxt = self._peek(1)
                if nxt == "(":
                    self.i += 2
                    self._subst()
                    word.append(SUBST)
                elif nxt == "'":
                    self.i += 1
                    word.append(self._ansi_c())
                elif nxt == "{":
                    j = s.find("}", self.i)
                    if j < 0:
                        raise LexError("unterminated ${")
                    word.append(s[self.i : j + 1])
                    self.i = j + 1
                else:
                    word.append("$")
                    self.i += 1
            elif c in "<>":
                if self._peek(1) == "(":  # process substitution <(..) >(..)
                    self.i += 2
                    self._subst()
                    if word is None:
                        word = []
                    word.append(SUBST)
                    continue
                if word is not None and word and "".join(word).isdigit():
                    word = None  # file descriptor number (2>&1)
                else:
                    finish_word()
                op = c
                self.i += 1
                if c == "<" and self._peek() == "<":
                    self.i += 1
                    if self._peek() == "<":
                        self.i += 1
                        redir_op = "<<<"
                        continue
                    strip = self._peek() == "-"
                    if strip:
                        self.i += 1
                    self.heredocs.append((self._read_heredoc_delim(), strip))
                    continue
                if self._peek() in (">", "&", "|") or (c == "<" and self._peek() == ">"):
                    op += self._peek()
                    self.i += 1
                redir_op = op
            elif c == "&" and self._peek(1) == ">":
                finish_word()
                self.i += 2
                op = "&>"
                if self._peek() == ">":
                    op = "&>>"
                    self.i += 1
                redir_op = op
            elif c in ";&|":
                finish_cmd()
                while self.i < self.n and s[self.i] in ";&|":
                    self.i += 1
            elif c == "(":
                finish_cmd()
                paren += 1
                self.i += 1
            elif c == ")":
                finish_cmd()
                self.i += 1
                if paren > 0:
                    paren -= 1
                elif in_subst:
                    return
            else:
                if word is None:
                    word = []
                word.append(c)
                self.i += 1
        if in_subst:
            raise LexError("unterminated $( or <(")
        finish_cmd()


def lex(command: str, depth: int = 0) -> list[Command]:
    return _Lexer(command, depth).run()


# ===================================================================== context
class Ctx:
    def __init__(self, data: dict[str, Any], root: Path, cfg: dict[str, Any]) -> None:
        self.data = data
        self.root = root
        self.root_real = os.path.realpath(str(root))
        self.cfg = cfg
        self.test_dirs = [_dir_prefix(d) for d in cfg.get("test_dirs") or [] if d]
        self.impl_dirs = [_dir_prefix(d) for d in cfg.get("impl_dirs") or [] if d]
        cwd = data.get("cwd")
        self.cwd = str(cwd) if isinstance(cwd, str) and os.path.isabs(cwd) else str(root)
        self._tasks: dict[str, dict[str, Any]] | None = None
        self._active: dict[str, dict[str, Any]] | None = None

    # -- tasks -------------------------------------------------------------------
    @property
    def tasks(self) -> dict[str, dict[str, Any]]:
        if self._tasks is None:
            self._tasks = C.load_tasks(self.root)
        return self._tasks

    @property
    def active(self) -> dict[str, dict[str, Any]]:
        if self._active is None:
            self._active = C.load_active(self.root)
        return self._active

    def in_progress(self) -> list[tuple[dict[str, Any], dict[str, Any]]]:
        """(active entry, task json) of registered tasks that are not DONE."""
        out = []
        for tid, entry in self.active.items():
            task = self.tasks.get(str(tid))
            if isinstance(task, dict) and task.get("phase") != "DONE":
                out.append((entry if isinstance(entry, dict) else {}, task))
        return out

    def reduction_allowed(self, rel: str) -> bool:
        """allow_test_reduction exception for one root-relative test path."""
        for entry, task in self.in_progress():
            if entry.get("allow_test_reduction") is not True:
                continue
            globs = (task.get("files") or {}).get("tests") or []
            if C.any_glob(rel, list(globs)) and C.has_event(task, "GREEN"):
                return True
        return False

    def reduction_allowed_path(self, abspath: str, rel: str) -> bool:
        if not os.path.isdir(abspath):
            return self.reduction_allowed(rel)
        files = _files_below(abspath)
        rels = [C.rel_to_root(f, self.root) for f in files]
        return bool(rels) and all(r is not None and self.reduction_allowed(r) for r in rels)

    # -- paths -------------------------------------------------------------------
    def is_test_rel(self, rel: str) -> bool:
        return any(rel == d.rstrip("/") or rel.startswith(d) for d in self.test_dirs)

    def is_impl_rel(self, rel: str) -> bool:
        return any(rel.startswith(d) for d in self.impl_dirs)


def _dir_prefix(d: str) -> str:
    d = str(d).replace("\\", "/")
    if d.startswith("./"):
        d = d[2:]
    return d if d.endswith("/") else d + "/"


def _files_below(path: str) -> list[str]:
    out = []
    for dirpath, _dirs, files in os.walk(path):
        out.extend(os.path.join(dirpath, f) for f in files)
    return out


def _expand(p: str) -> str:
    return os.path.expanduser(os.path.expandvars(p))


def _abspath(p: str, base: str) -> str:
    p = _expand(p)
    if not os.path.isabs(p):
        p = os.path.join(base, p)
    return os.path.normpath(p)


def _has_glob(p: str) -> bool:
    return any(ch in p for ch in "*?[")


def _is_below(path: str, parent: str) -> bool:
    parent = parent.rstrip("/") or "/"
    if parent == "/":
        return path != "/"
    return path.startswith(parent + "/")


# ======================================================================= messages
MSG = {
    "push": (
        "force push 禁止 (`git push --force` / `-f` / `--force-with-lease` / "
        "`--force-if-includes` / `+refspec` / `--mirror`)。リモートの履歴を書き換えない。"
        "代替: 通常の `git push origin <branch>`。"
        "衝突は `git pull --rebase` か追加コミットで解決する。"
    ),
    "reset": (
        "`git reset --hard` 禁止 (未コミット変更を破棄する)。"
        "代替: `git stash`、`git reset --soft` / `--mixed`、取り消しは `git revert`。"
    ),
    "checkout": (
        "`git checkout` による作業ツリー変更の破棄 (`.` / `-- <path>` / `<path>` / `-f`) 禁止。"
        "代替: `git stash` で退避する。ブランチ切替は `git switch <branch>`。"
    ),
    "restore": (
        "`git restore` は `--staged` (`-S`) のみ (index から外すだけ) 許可。"
        "作業ツリーの変更を破棄する `git restore <path>` / `.` / `--worktree` は禁止。"
        "代替: `git stash`、または Edit で必要な箇所だけ戻す。"
    ),
    "clean": (
        "`git clean -f` 禁止 (未追跡ファイルを消す)。"
        "代替: `git clean -n` で一覧を確認し、消すものは利用者に確認する。"
    ),
    "branch": (
        "強制ブランチ削除 (`git branch -D` / `--delete --force`) 禁止。"
        "代替: マージ済みなら `git branch -d <name>`。未マージの削除は利用者に依頼する。"
    ),
    "remote": (
        "`git remote set-url/remove/rm/rename/add` 禁止 (remote が1つも無いときの "
        "`git remote add origin <url>` だけ可)。remote の変更は利用者に依頼する。"
    ),
    "gh": (
        "`gh repo edit --visibility` / `gh repo delete` 禁止 (公開範囲の変更・リポジトリ削除)。"
        "必要なら利用者が自分で行う。"
    ),
    "tests": (
        "テストの削除・弱体化は禁止 (test deletion / assertion reduction)。"
        "代替: テストは残したまま実装を直す。"
        "どうしても減らすなら利用者の承認を得て、作業中タスクの "
        "active.json エントリに `allow_test_reduction: true` と `allow_reason` を設定し "
        "(scripts/tdd.py で登録)、`files.tests` の GLOB が対象に一致し GREEN 証跡がある状態で行う。"
    ),
}


def _rm_msg(path: str, why: str) -> str:
    return (
        f"`rm -r/-f` の対象 `{path}` は{why}ため禁止 (dangerous rm)。"
        "消してよいのはリポジトリ内 (ルート自体と `.git` を除く) と scratchpad / `/tmp` / "
        "`$TMPDIR` 配下だけ。代替: 対象をリポジトリ内の具体的なパスに絞る。"
    )


def _test_msg(rel: str, how: str) -> str:
    return f"`{rel}`: {how}。" + MSG["tests"]


# ============================================================ (a)+(b) Bash checks
ASSIGN_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(?:\[[^]]*\])?\+?=")
KEYWORDS = {"if", "then", "else", "elif", "fi", "do", "done", "while", "until", "!", "{", "}"}
SHELLS = {"bash", "sh", "zsh", "dash", "ksh", "mksh", "ash"}
WRAPPERS: dict[str, set[str]] = {
    "sudo": {
        "-u", "-g", "-C", "-h", "-p", "-U", "-r", "-t", "-D", "-T", "-R",
        "--user", "--group", "--host", "--prompt", "--chdir", "--role", "--type",
        "--other-user", "--close-from", "--command-timeout",
    },
    "env": {"-u", "-C", "--unset", "--chdir"},
    "command": set(),
    "builtin": set(),
    "nohup": set(),
    "time": set(),
    "nice": {"-n", "--adjustment"},
    "timeout": {"-s", "-k", "--signal", "--kill-after"},
    "stdbuf": {"-i", "-o", "-e", "--input", "--output", "--error"},
    "xargs": {
        "-I", "-n", "-P", "-L", "-d", "-E", "-s", "-a",
        "--arg-file", "--delimiter", "--max-args", "--max-procs", "--max-lines",
        "--replace", "--eof", "--max-chars", "--process-slot-var",
    },
    "noglob": set(),
    "exec": {"-a"},
}  # fmt: skip
GIT_OPTS_WITH_ARG = {
    "-C", "-c", "--git-dir", "--work-tree", "--namespace", "--super-prefix",
    "--config-env",
}  # fmt: skip


def strip_wrappers(argv: list[str]) -> tuple[list[str], list[str]]:
    """Drop ``VAR=val`` prefixes, keywords and wrapper commands.

    Returns (remaining argv, command strings to analyse recursively).
    """
    nested: list[str] = []
    i = 0
    while i < len(argv):
        w = argv[i]
        if ASSIGN_RE.match(w) or w in KEYWORDS:
            i += 1
            continue
        name = os.path.basename(w)
        if name not in WRAPPERS:
            break
        with_arg = WRAPPERS[name]
        i += 1
        while i < len(argv):
            a = argv[i]
            if a == "--":
                i += 1
                break
            if name == "env":
                if ASSIGN_RE.match(a) or a == "-":
                    i += 1
                    continue
                if a in ("-S", "--split-string") and i + 1 < len(argv):
                    nested.append(argv[i + 1])
                    i += 2
                    continue
                if a.startswith("--split-string="):
                    nested.append(a.split("=", 1)[1])
                    i += 1
                    continue
                if a.startswith("-S") and len(a) > 2:
                    nested.append(a[2:])
                    i += 1
                    continue
            if not a.startswith("-") or a == "-":
                break
            if name == "nice" and re.match(r"^-\d+$", a):
                i += 1
                continue
            i += 2 if a in with_arg else 1
        if name == "timeout" and i < len(argv):
            i += 1  # the duration
    return argv[i:], nested


def split_opts(
    args: list[str], with_arg: set[str] | None = None
) -> tuple[list[str], list[str], list[str]]:
    """(options, positionals before ``--``, words after ``--``)."""
    with_arg = with_arg or set()
    opts: list[str] = []
    pos: list[str] = []
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--":
            return opts, pos, args[i + 1 :]
        if a.startswith("-") and a != "-":
            opts.append(a)
            if a in with_arg:
                i += 1  # skip the option's value
        else:
            pos.append(a)
        i += 1
    return opts, pos, []


def has_short(opts: list[str], ch: str) -> bool:
    return any(o.startswith("-") and not o.startswith("--") and ch in o[1:] for o in opts)


def has_long(opts: list[str], full: str, minlen: int = 4) -> bool:
    """``--name`` or an unambiguous-looking abbreviation of it (git accepts those)."""
    for o in opts:
        if not o.startswith("--"):
            continue
        name = o.split("=", 1)[0]
        if len(name) >= minlen and full.startswith(name):
            return True
    return False


class BashChecker:
    def __init__(self, ctx: Ctx) -> None:
        self.ctx = ctx
        self.reasons: list[str] = []

    def add(self, msg: str) -> None:
        if msg not in self.reasons:
            self.reasons.append(msg)

    # -- entry -------------------------------------------------------------------
    def check(self, command: str, cwd: str, depth: int = 0) -> None:
        if depth > MAX_DEPTH:
            self.add("シェルの入れ子 (bash -c / eval) が深すぎるため判定できず拒否 (fail-closed)。")
            return
        try:
            cmds = lex(command)
        except LexError as exc:
            self.fallback(command, str(exc))
            return
        for cmd in cmds:
            cwd = self.check_command(cmd, cwd, depth)

    def fallback(self, command: str, err: str) -> None:
        """Conservative regex check of the raw string when lexing failed."""
        hits = [label for label, rx in FALLBACK_PATTERNS if rx.search(command)]
        if hits:
            self.add(
                f"コマンドを字句解析できなかった ({err}) ため生文字列で保守的に判定し、"
                f"危険な操作らしきもの ({', '.join(hits)}) を検出したので拒否 (fail-closed)。"
                "引用符を閉じた正しいコマンドにする。"
            )

    def check_command(self, cmd: Command, cwd: str, depth: int) -> str:
        """Check one simple command; returns the cwd for the following commands."""
        for op, target in cmd.redirs:
            if op in (">", ">|", "&>", "<>"):
                self.check_test_path(
                    target, cwd, "シェルのリダイレクトで上書き", existing_file=True
                )
        argv, nested = strip_wrappers(cmd.argv)
        for s in nested:
            self.check(s, cwd, depth + 1)
        if not argv:
            return cwd
        name = os.path.basename(argv[0])
        args = argv[1:]
        if name in ("cd", "pushd"):
            target = next((a for a in args if not a.startswith("-")), None)
            if target is None:
                return os.path.expanduser("~")
            if "$" in target:
                return cwd
            return _abspath(target, cwd)
        if name in SHELLS:
            script = _shell_c_arg(args)
            if script is not None:
                self.check(script, cwd, depth + 1)
        elif name == "eval":
            self.check(" ".join(args), cwd, depth + 1)
        elif name == "git":
            self.check_git(args, cwd)
        elif name == "gh":
            self.check_gh(args)
        elif name == "rm":
            self.check_rm(args, cwd)
        elif name == "unlink":
            _opts, pos, dd = split_opts(args)
            for p in pos + dd:
                self.check_test_path(p, cwd, "unlink で削除")
        elif name == "truncate":
            _opts, pos, dd = split_opts(args, {"-s", "-r", "--size", "--reference"})
            for p in pos + dd:
                self.check_test_path(p, cwd, "truncate で中身を消去")
        elif name == "find":
            self.check_find(args, cwd)
        return cwd

    # -- git -----------------------------------------------------------------------
    def check_git(self, args: list[str], cwd: str) -> None:
        i = 0
        cdir: str | None = None
        while i < len(args) and args[i].startswith("-"):
            a = args[i]
            if a in GIT_OPTS_WITH_ARG:
                if a == "-C" and i + 1 < len(args):
                    cdir = _abspath(args[i + 1], cdir or cwd)
                i += 2
                continue
            i += 1
        if i >= len(args):
            return
        sub, rest = args[i], args[i + 1 :]
        base = cdir or cwd
        if sub == "push":
            opts, pos, dd = split_opts(rest, {"-o", "--push-option", "--repo", "--receive-pack"})
            if (
                has_short(opts, "f")
                or has_long(opts, "--force")
                or has_long(opts, "--force-with-lease")
                or has_long(opts, "--force-if-includes")
                or has_long(opts, "--mirror")
                or any(p.startswith("+") for p in pos + dd)
            ):
                self.add(MSG["push"])
        elif sub == "reset":
            opts, _pos, _dd = split_opts(rest)
            if has_long(opts, "--hard"):
                self.add(MSG["reset"])
        elif sub == "checkout":
            opts, pos, dd = split_opts(rest, {"-b", "-B", "--orphan"})
            existing = [p for p in pos if os.path.lexists(_abspath(p, base))]
            if has_short(opts, "f") or has_long(opts, "--force") or "." in pos or dd or existing:
                self.add(MSG["checkout"])
        elif sub == "restore":
            opts, _pos, _dd = split_opts(rest, {"-s", "--source"})
            staged = has_short(opts, "S") or has_long(opts, "--staged")
            worktree = has_short(opts, "W") or has_long(opts, "--worktree")
            if not staged or worktree:
                self.add(MSG["restore"])
        elif sub == "clean":
            opts, _pos, _dd = split_opts(rest, {"-e", "--exclude"})
            if has_short(opts, "f") or has_long(opts, "--force"):
                self.add(MSG["clean"])
        elif sub == "branch":
            opts, _pos, _dd = split_opts(rest)
            delete = has_short(opts, "d") or has_long(opts, "--delete")
            force = has_short(opts, "f") or has_long(opts, "--force")
            if has_short(opts, "D") or (delete and force):
                self.add(MSG["branch"])
        elif sub == "remote":
            self.check_git_remote(rest, cdir)
        elif sub == "rm":
            _opts, pos, dd = split_opts(rest, {"--pathspec-from-file"})
            for p in pos + dd:
                self.check_test_path(p, base, "git rm で削除")

    def check_git_remote(self, rest: list[str], cdir: str | None) -> None:
        idx = next((k for k, a in enumerate(rest) if not a.startswith("-")), None)
        if idx is None:
            return
        sub = rest[idx]
        if sub not in ("set-url", "remove", "rm", "rename", "add"):
            return
        if sub == "add":
            opts, pos, dd = split_opts(rest[idx + 1 :], {"-t", "-m"})
            names = pos + dd
            if len(names) == 2 and names[0] == "origin" and not has_long(opts, "--mirror"):
                where = cdir or str(self.ctx.root)
                code, out = C.run(["git", "-C", where, "remote"], timeout=3.0)
                if code != 0:
                    self.add(
                        "remote 一覧 (`git remote`) の取得に失敗したため "
                        f"`git remote add origin` を拒否 (fail-closed, exit {code})。"
                        + MSG["remote"]
                    )
                    return
                if not out.strip():
                    return  # first remote of a fresh repository: allowed
                self.add("既に remote があるため `git remote add` は不可。" + MSG["remote"])
                return
        self.add(MSG["remote"])

    # -- gh ------------------------------------------------------------------------
    def check_gh(self, args: list[str]) -> None:
        pos = [a for a in args if not a.startswith("-")]
        if len(pos) < 2 or pos[0] != "repo":
            return
        visibility = any(a == "--visibility" or a.startswith("--visibility=") for a in args)
        if pos[1] == "delete" or (pos[1] == "edit" and visibility):
            self.add(MSG["gh"])

    # -- rm / find -------------------------------------------------------------------
    def check_rm(self, args: list[str], cwd: str) -> None:
        opts, pos, dd = split_opts(args)
        recursive = has_short(opts, "r") or has_short(opts, "R") or has_long(opts, "--recursive", 5)
        force = has_short(opts, "f") or has_long(opts, "--force")
        for p in pos + dd:
            self.check_test_path(p, cwd, "rm で削除")
            if recursive or force:
                self.check_rm_target(p, cwd)

    def check_rm_target(self, p: str, cwd: str) -> None:
        exp = _expand(p)
        if "$" in exp or "`" in exp:
            self.add(_rm_msg(p, "展開できない変数・コマンド置換を含む"))
            return
        absp = _abspath(exp, cwd)
        cands = [absp]
        if _has_glob(absp):
            cands = sorted(_glob.glob(absp)) or [absp]
        for cand in cands:
            why = self._rm_danger(cand)
            if why:
                self.add(_rm_msg(p, why))
                return

    def _rm_danger(self, absp: str) -> str | None:
        real = os.path.realpath(absp)
        parts = set(Path(absp).parts) | set(Path(real).parts)
        if ".git" in parts:
            return "`.git` を含む"
        if real == "/" or absp == "/":
            return "ルートディレクトリ `/` である"
        root = self.ctx.root_real
        if real == root or _is_below(root, real):
            return "リポジトリのルートそのもの (またはそれを含むディレクトリ) である"
        if _is_below(real, root):
            return None
        for safe in self._safe_dirs():
            if _is_below(real, safe):
                return None
        return "リポジトリの外にある"

    def _safe_dirs(self) -> list[str]:
        dirs = ["/tmp", "/private/tmp"]
        scratch = self.ctx.data.get("scratchpad_dir")
        if isinstance(scratch, str) and scratch:
            dirs.append(scratch)
        tmpdir = os.environ.get("TMPDIR")
        if tmpdir:
            dirs.append(tmpdir)
        out = []
        for d in dirs:
            for v in (os.path.normpath(d), os.path.realpath(d)):
                if v not in out:
                    out.append(v)
        return out

    def check_find(self, args: list[str], cwd: str) -> None:
        deletes = "-delete" in args
        for k, a in enumerate(args):
            if (
                a in ("-exec", "-execdir", "-ok", "-okdir")
                and k + 1 < len(args)
                and os.path.basename(args[k + 1]) in ("rm", "unlink", "truncate")
            ):
                deletes = True
        if not deletes:
            return
        starts = []
        for a in args:
            if a.startswith("-") or a in ("(", "!", ")"):
                break
            starts.append(a)
        for p in starts or ["."]:
            self.check_test_path(p, cwd, "find -delete で削除")

    # -- tests -----------------------------------------------------------------------
    def check_test_path(self, p: str, base: str, how: str, existing_file: bool = False) -> None:
        if not p or p == SUBST:
            return
        absp = _abspath(p, base)
        cands = [absp]
        if _has_glob(absp):
            cands += sorted(_glob.glob(absp))
        for cand in cands:
            rel = C.rel_to_root(cand, self.ctx.root)
            if rel is None or rel in ("", ".") or not self.ctx.is_test_rel(rel):
                continue
            if existing_file and not os.path.isfile(cand):
                continue
            if os.path.isdir(cand) and not _files_below(cand):
                continue  # empty directory
            if self.ctx.reduction_allowed_path(cand, rel):
                continue
            self.add(_test_msg(rel, how))
            return


def _shell_c_arg(args: list[str]) -> str | None:
    """The script of ``bash -c '<script>'`` (also bundles like ``-lc``)."""
    i = 0
    found = False
    while i < len(args):
        a = args[i]
        if a == "--":
            i += 1
            break
        if a in ("-o", "+o", "-O", "+O", "--rcfile", "--init-file"):
            i += 2
            continue
        if (a.startswith("-") or a.startswith("+")) and len(a) > 1:
            if not a.startswith("--") and "c" in a[1:]:
                found = True
            i += 1
            continue
        break
    if found and i < len(args):
        return args[i]
    return None


def _rx(p: str) -> re.Pattern[str]:
    return re.compile(p, re.S)


FALLBACK_PATTERNS = [
    ("git push --force", _rx(r"\bpush\b.*(--force|--mirror|\s-[A-Za-z]*f|\s\+\S)")),
    ("git reset --hard", _rx(r"\breset\b.*--hard")),
    ("git checkout/restore", _rx(r"\bgit\b.*\b(checkout|restore)\b")),
    ("git clean -f", _rx(r"\bclean\b.*(\s-[A-Za-z]*f|--force)")),
    ("git branch -D", _rx(r"\bbranch\b.*(\s-[A-Za-z]*D|--delete)")),
    ("git remote", _rx(r"\bremote\b.*\b(set-url|remove|rm|rename|add)\b")),
    ("gh repo", _rx(r"\bgh\b.*\brepo\b.*\b(delete|edit)\b")),
    ("rm -r/-f", _rx(r"\brm\b.*(\s-[A-Za-z]*[rRf]|--recursive|--force)")),
    ("test deletion", _rx(r"\b(rm|unlink|truncate)\b.*\btests?\b")),
]


# ======================================================== (b) file tool checks
ASSERT_RE = re.compile(r"\b(?:STATIC_REQUIRE\w*|STATIC_CHECK\w*|REQUIRE\w*|CHECK\w*)\s*\(")
PY_ASSERT_RE = re.compile(r"\bassert\b")
STATIC_BARE_RE = re.compile(r"\bSTATIC_REQUIRE\w*\b(?!\s*\()")
CASE_RE = re.compile(r"\b(?:TEMPLATE_\w*TEST_CASE\w*|TEST_CASE\w*|SCENARIO\w*)")
PY_CASE_RE = re.compile(r"\bdef\s+test_\w+")


def count_tests(text: str) -> tuple[int, int]:
    """(assertions, test cases) in a C++ (Catch2) or Python test source."""
    asserts = (
        len(ASSERT_RE.findall(text))
        + len(PY_ASSERT_RE.findall(text))
        + len(STATIC_BARE_RE.findall(text))
    )
    cases = len(CASE_RE.findall(text)) + len(PY_CASE_RE.findall(text))
    return asserts, cases


def _read_text(path: str) -> str | None:
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except OSError:
        return None


def _apply_edits(content: str, edits: list[dict[str, Any]]) -> str | None:
    """Apply Edit/MultiEdit replacements in order; None when one does not apply."""
    for e in edits:
        old = str(e.get("old_string") or "")
        new = str(e.get("new_string") or "")
        if not old or old not in content:
            return None
        content = (
            content.replace(old, new) if e.get("replace_all") else content.replace(old, new, 1)
        )
    return content


def check_file_tool(ctx: Ctx, tool: str, ti: dict[str, Any]) -> list[str]:
    if tool == "NotebookEdit":
        path = ti.get("notebook_path")
        if not isinstance(path, str) or not path:
            return []
        rel = C.rel_to_root(path, ctx.root)
        if rel is None or not ctx.is_test_rel(rel):
            return []
        if ti.get("edit_mode") == "delete" and not ctx.reduction_allowed(rel):
            return [_test_msg(rel, "NotebookEdit でセルを削除 (edit_mode=delete)")]
        return []

    path = ti.get("file_path")
    if not isinstance(path, str) or not path:
        return []
    rel = C.rel_to_root(path, ctx.root)
    if rel is None or not ctx.is_test_rel(rel):
        return []
    old_text = _read_text(path) if os.path.isfile(path) else None

    if tool == "Write":
        if old_text is None:
            return []  # new test file
        before, after = old_text, str(ti.get("content") or "")
        if before.strip() and not after.strip():
            if ctx.reduction_allowed(rel):
                return []
            return [_test_msg(rel, "Write で既存のテストファイルを空にする")]
    else:
        if tool == "Edit":
            edits = [ti]
        else:
            raw = ti.get("edits")
            edits = [e for e in raw if isinstance(e, dict)] if isinstance(raw, list) else []
        applied = _apply_edits(old_text, edits) if old_text is not None else None
        if applied is not None:
            before, after = old_text or "", applied
        else:
            before = "\n".join(str(e.get("old_string") or "") for e in edits)
            after = "\n".join(str(e.get("new_string") or "") for e in edits)

    a0, c0 = count_tests(before)
    a1, c1 = count_tests(after)
    if (a1 < a0 or c1 < c0) and not ctx.reduction_allowed(rel):
        return [
            _test_msg(
                rel,
                f"{tool} で assertion 数 {a0}→{a1} / test case 数 {c0}→{c1} に減る",
            )
        ]
    return []


# ============================================================ (c) TDD enforcement
def check_tdd(ctx: Ctx, tool: str, ti: dict[str, Any]) -> tuple[str, str] | None:
    """('deny'|'warn', message) when an impl file is edited without RED evidence."""
    if tool not in FILE_TOOLS:
        return None
    path = ti.get("notebook_path") if tool == "NotebookEdit" else ti.get("file_path")
    if not isinstance(path, str) or not path:
        return None
    rel = C.rel_to_root(path, ctx.root)
    if rel is None or not ctx.is_impl_rel(rel):
        return None
    for _entry, task in ctx.in_progress():
        impl = (task.get("files") or {}).get("impl") or []
        if C.any_glob(rel, list(impl)) and C.has_event(task, "RED"):
            return None
    msg = (
        f"TDD: 実装ファイル `{rel}` を編集する前に RED 証跡が必要 (no active task with RED "
        "evidence whose files.impl matches)。手順: "
        "`uv run python scripts/tdd.py start <ID> --spec .. --impl <glob> --test <glob>` → "
        "テストを書いて "
        "`uv run python scripts/tdd.py red <ID> --reason .. -- <cmd>` を先に実行する "
        "(tdd skill 参照)。テスト (tests/) や docs は RED 前でも編集できる。"
    )
    return ("deny" if ctx.cfg.get("enforce") is True else "warn", msg)


# ============================================================================ main
def _deny(reasons: list[str]) -> None:
    text = "[pre_tool_use hook] この操作を拒否しました (denied by project hook):\n" + "\n".join(
        f"- {r}" for r in reasons
    )
    C.emit(
        {
            "hookSpecificOutput": {
                "hookEventName": EVENT,
                "permissionDecision": "deny",
                "permissionDecisionReason": text,
            }
        }
    )


def evaluate(data: dict[str, Any]) -> tuple[list[str], list[str]]:
    """(deny reasons, warnings) for one event."""
    tool = data.get("tool_name")
    ti = data.get("tool_input")
    if not isinstance(tool, str) or not isinstance(ti, dict):
        return [], []
    denies: list[str] = []
    warns: list[str] = []

    # (a)+(b): fail-closed
    try:
        if tool == "Bash":
            command = ti.get("command")
            if not isinstance(command, str) or not command.strip():
                return [], []
            root = C.project_root(data)
            ctx = Ctx(data, root, C.load_config(root))
            checker = BashChecker(ctx)
            checker.check(command, ctx.cwd)
            denies.extend(checker.reasons)
        elif tool in FILE_TOOLS:
            hint = ti.get("notebook_path") if tool == "NotebookEdit" else ti.get("file_path")
            hint = hint if isinstance(hint, str) and hint else None
            root = C.project_root(data, path_hint=hint)
            ctx = Ctx(data, root, C.load_config(root))
            denies.extend(check_file_tool(ctx, tool, ti))
        else:
            return [], []
    except Exception as exc:
        denies.append(
            f"安全チェック中に内部エラーが起きたため拒否 (fail-closed): {type(exc).__name__}: {exc}"
        )
        return denies, warns

    # (c): fail-open
    try:
        if tool in FILE_TOOLS:
            res = check_tdd(ctx, tool, ti)
            if res is not None:
                (denies if res[0] == "deny" else warns).append(res[1])
    except Exception as exc:
        C.debug(f"pre_tool_use: TDD check skipped (fail-open): {type(exc).__name__}: {exc}")
    return denies, warns


def main() -> int:
    data = C.read_input()
    if not data:
        return 0
    denies, warns = evaluate(data)
    if denies:
        _deny(denies)
    elif warns:
        C.emit_context(EVENT, "[pre_tool_use hook] 警告 (warning):\n" + "\n".join(warns))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:  # pragma: no cover - last resort, never block by crash
        C.debug(f"pre_tool_use: unexpected error: {type(exc).__name__}: {exc}")
        sys.exit(0)
