"""Tests for scripts/check_public_safety.py (pre-publication audit of tracked files).

Sensitive-looking strings are assembled at runtime so this file itself stays clean.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "check_public_safety.py"
HOME = "/" + "Users" + "/alice/work/x.cpp"
TMP = "/private" + "/tmp/claude-123/scratch"
MAIL = "alice" + "@" + "example-corp.jp"
TOKEN = "gh" + "p_" + "A" * 36


def check(tmp_path: Path, files: dict[str, bytes | str]) -> subprocess.CompletedProcess:
    for name, body in files.items():
        p = tmp_path / name
        p.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(body, str):
            p.write_text(body)
        else:
            p.write_bytes(body)
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--root", str(tmp_path), *files],
        capture_output=True,
        text=True,
    )


def test_clean_files_pass(tmp_path):
    p = check(tmp_path, {"a.md": "relative path src/graph/io.cpp and https://arxiv.org/abs/x\n"})
    assert p.returncode == 0, p.stdout


def test_home_and_tmp_absolute_paths_are_rejected(tmp_path):
    p = check(tmp_path, {"a.md": f"see {HOME}\n", "b.json": f'"{TMP}"\n'})
    assert p.returncode == 1
    assert "a.md:1" in p.stdout and "b.json:1" in p.stdout


def test_email_and_tokens_are_rejected(tmp_path):
    p = check(tmp_path, {"c.txt": f"contact {MAIL}\n", "d.yml": f"token: {TOKEN}\n"})
    assert p.returncode == 1
    assert "c.txt:1" in p.stdout and "d.yml:1" in p.stdout


def test_noreply_and_example_emails_are_allowed(tmp_path):
    body = "x@users.noreply.github.com\nuser" + "@example.com\nnoreply" + "@anthropic.com\n"
    p = check(tmp_path, {"e.md": body})
    assert p.returncode == 0, p.stdout


def test_paper_copies_are_rejected(tmp_path):
    p = check(tmp_path, {"docs/2611.34911v1.pdf": b"%PDF-1.7", "x/paper.pdf": b"%PDF-1.4"})
    assert p.returncode == 1
    assert "2611.34911v1.pdf" in p.stdout and "paper.pdf" in p.stdout


def test_inline_allow_marker_suppresses_a_line(tmp_path):
    p = check(tmp_path, {"f.md": f"{HOME}  <!-- public-safety: allow -->\n"})
    assert p.returncode == 0, p.stdout


def test_default_mode_scans_git_tracked_files_of_this_repo():
    p = subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, capture_output=True, text=True)
    assert p.returncode in (0, 1)
    assert "scanned" in p.stdout
