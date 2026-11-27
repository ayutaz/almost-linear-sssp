"""Tests for scripts/fetch_paper.sh (paper is fetched, never redistributed)."""

import hashlib
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "fetch_paper.sh"


def run(tmp_path: Path, content: bytes, sha: str) -> subprocess.CompletedProcess:
    src = tmp_path / "src.pdf"
    src.write_bytes(content)
    env = dict(
        os.environ, PAPER_URL=src.as_uri(), PAPER_SHA256=sha, PAPER_DIR=str(tmp_path / "out")
    )
    return subprocess.run(["bash", str(SCRIPT)], env=env, capture_output=True, text=True)


def test_fetch_accepts_matching_checksum(tmp_path):
    data = b"%PDF-1.7 fake"
    p = run(tmp_path, data, hashlib.sha256(data).hexdigest())
    assert p.returncode == 0, p.stderr
    assert (tmp_path / "out" / "2611.34911v1.pdf").read_bytes() == data


def test_fetch_rejects_checksum_mismatch(tmp_path):
    p = run(tmp_path, b"%PDF-1.7 tampered", "0" * 64)
    assert p.returncode != 0
    assert not (tmp_path / "out" / "2611.34911v1.pdf").exists()
    assert "checksum" in p.stderr.lower()


def test_fetch_target_is_gitignored():
    p = subprocess.run(["git", "check-ignore", ".cache/paper/2611.34911v1.pdf"], cwd=ROOT)
    assert p.returncode == 0
