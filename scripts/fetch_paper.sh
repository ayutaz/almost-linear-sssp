#!/usr/bin/env bash
# Fetch the target paper (arXiv 2611.34911v1) into .cache/paper/ for local analysis.
# The paper is distributed by arXiv under its non-exclusive distribution licence, so it is
# NOT committed to this repository. The download is verified against a pinned SHA-256.
set -euo pipefail

PAPER_ID="2611.34911v1"
PAPER_URL="${PAPER_URL:-https://arxiv.org/pdf/${PAPER_ID}}"
PAPER_SHA256="${PAPER_SHA256:-0dc54949de1f27014e25609c44146a2718142f2a97205992423e78e692d609bc}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PAPER_DIR="${PAPER_DIR:-${ROOT}/.cache/paper}"
OUT="${PAPER_DIR}/${PAPER_ID}.pdf"

mkdir -p "${PAPER_DIR}"
tmp="$(mktemp "${PAPER_DIR}/.download.XXXXXX")"
trap 'rm -f "${tmp}"' EXIT

curl --fail --silent --show-error --location --output "${tmp}" "${PAPER_URL}"

if command -v sha256sum >/dev/null 2>&1; then
  actual="$(sha256sum "${tmp}" | awk '{print $1}')"
else
  actual="$(shasum -a 256 "${tmp}" | awk '{print $1}')"
fi
if [[ "${actual}" != "${PAPER_SHA256}" ]]; then
  echo "error: checksum mismatch for ${PAPER_URL} (expected ${PAPER_SHA256}, got ${actual})" >&2
  exit 1
fi
mv "${tmp}" "${OUT}"
trap - EXIT

# The verified PDF is the deliverable; the text copy is an optional convenience. A failed
# extraction does not fail the fetch, but it is reported and never leaves an older copy behind.
TXT="${PAPER_DIR}/${PAPER_ID}.txt"
if command -v pdftotext >/dev/null 2>&1; then
  if ! pdftotext -layout "${OUT}" "${TXT}"; then
    rm -f "${TXT}"
    echo "warning: pdftotext failed; ${TXT} was not written (the PDF is fine)" >&2
  fi
fi
echo "paper saved to ${OUT}"
