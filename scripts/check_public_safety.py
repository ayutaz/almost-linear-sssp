#!/usr/bin/env python3
"""Pre-publication audit: this repository is public OSS.

Rejects, in files that would be committed (git tracked + untracked-not-ignored):
  * local absolute paths (home directories, macOS temp/scratch directories)
  * e-mail addresses (except noreply / example domains)
  * credential-looking tokens and private keys
  * copies of the target paper (arXiv licence forbids redistribution)
A line containing ``public-safety: allow`` is exempt (use sparingly, with a reason).
Exit code: 0 clean, 1 findings.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ALLOW_MARKER = "public-safety: allow"
EMAIL_ALLOW = re.compile(
    r"@(users\.noreply\.github\.com|example\.(com|org|net)|anthropic\.com)$", re.I
)
LINE_RULES: list[tuple[str, re.Pattern[str]]] = [
    ("local home path", re.compile(r"/(Users|home)/[A-Za-z0-9._-]+/")),
    ("local temp path", re.compile(r"/private/(tmp|var)/|/var/folders/")),  # public-safety: allow
    (
        "credential token",
        re.compile(
            r"\b(gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|sk-ant-[A-Za-z0-9_-]{20,}"
            r"|AKIA[0-9A-Z]{16}|xox[baprs]-[A-Za-z0-9-]{10,})"
        ),
    ),
    ("private key", re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----")),
]
EMAIL = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
PAPER_FILE = re.compile(r"(2611\.34911|(^|/)paper\.pdf$)", re.I)
MAX_BYTES = 5 * 1024 * 1024


def candidate_files(root: Path) -> list[str]:
    out = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        cwd=root,
        capture_output=True,
        check=True,
    ).stdout.decode()
    return [f for f in out.split("\0") if f and (root / f).is_file()]


def scan(root: Path, files: list[str]) -> list[str]:
    findings: list[str] = []
    for rel in files:
        path = root / rel
        if PAPER_FILE.search(rel) and rel.lower().endswith((".pdf", ".txt", ".tex")):
            findings.append(f"{rel}: copy of the paper must not be committed")
            continue
        data = path.read_bytes()
        if len(data) > MAX_BYTES:
            findings.append(f"{rel}: file larger than 5 MiB")
            continue
        if b"\0" in data[:4096]:
            continue  # binary
        for no, line in enumerate(data.decode("utf-8", "replace").splitlines(), 1):
            if ALLOW_MARKER in line:
                continue
            for name, pat in LINE_RULES:
                if pat.search(line):
                    findings.append(f"{rel}:{no}: {name}")
            for m in EMAIL.finditer(line):
                if not EMAIL_ALLOW.search(m.group(0)):
                    findings.append(f"{rel}:{no}: e-mail address")
    return findings


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("files", nargs="*")
    args = ap.parse_args(argv)
    root = (args.root or Path.cwd()).resolve()
    files = args.files or candidate_files(root)
    findings = scan(root, files)
    for f in findings:
        print(f)
    print(f"public-safety: scanned {len(files)} files, {len(findings)} finding(s)")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
