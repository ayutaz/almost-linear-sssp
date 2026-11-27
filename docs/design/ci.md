# CI design (GitHub Actions)

**English** | [日本語](ci.ja.md) | [简体中文](ci.zh-CN.md)

This document describes the four workflows under `.github/workflows/` and the design of `scripts/bootstrap.sh` and `scripts/coverage.sh`.
The local quality gate is `make verify`, and CI follows the same commands (CMake presets, ctest, ruff, pytest).

## Workflow overview

| File | Trigger | Contents |
|---|---|---|
| `ci.yml` | `push` (main), `pull_request`, manual | lint (clang-format / clang-tidy / ruff), builds of `ci` (Debug, `-Werror`) and `ci-release` (Release, `-Werror`), CTest per label, pytest |
| `sanitizer.yml` | `push` (main), `pull_request`, weekly, manual | `sanitizer` preset (ASan + UBSan). LeakSanitizer is also enabled on Linux |
| `benchmark.yml` | `push` (main, on code changes), manual (with inputs), weekly | push runs smoke, schedule runs standard, manual lets you choose the suite. Saves CSV/JSONL as artifacts and summarises them in the Summary |
| `release.yml` | push of a `v*` tag | tests -> `ci-release` build -> pack the `tools/` binaries into tar.gz -> SHA-256 -> attach to a GitHub Release |

### ci.yml

- **matrix**: `os = [macos-26, ubuntu-26.04]` x `preset = [ci, ci-release]` (`fail-fast: false`).
- **CTest runs as a separate step per label**: `unit`, `integration`, `differential`, `property`, `regression`.
  - Each label is selected with an exact-match regex such as `-L '^unit$'`.
  - `--no-tests=ignore` is passed, so a label with zero tests does not fail.
  - The last step uses `-LE` to pick up tests that carry no label at all (a safety net against forgotten labels).
  - If one label fails, the remaining labels still run (`!cancelled() && steps.build.outcome == 'success'`).
- **pytest** runs only in the job for the `ci` preset (`tests/harness`, `tests/python`). Zero collected tests (exit code 5) counts as success.
- **lint** runs on Ubuntu only.
  - clang-format runs `--dry-run --Werror` on the C++ files from `git ls-files`.
  - clang-tidy runs `make tidy PRESET=ci` (targets `src/` and `tools/`). It fails through `WarningsAsErrors` in `.clang-tidy`.
  - ruff runs `check` and `format --check`.
- **`ci result` job**: it folds the results of all jobs into one. Branch protection should require this job as the required check, not the jobs whose names change with the matrix.
- On failure, `build/<preset>/Testing/Temporary/` is saved as an artifact.

### sanitizer.yml

- `ASAN_OPTIONS` differs per OS.
  - Linux: `detect_leaks=1` (LeakSanitizer enabled), `strict_string_checks`, `detect_stack_use_after_return`, `check_initialization_order`, and so on.
  - macOS: `detect_leaks=0`. ASan in Apple clang does not support LeakSanitizer on macOS. Leak detection is the job of the Linux job.
- `UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1`. The preset also passes `-fno-sanitize-recover=undefined`, so UB always fails the run.
- On Linux, `vm.mmap_rnd_bits` is lowered to 28. Newer kernels use a larger ASLR entropy, and there is a known issue that breaks the ASan shadow mapping (it is unlikely with clang 22, but this is a precaution).

### benchmark.yml

- **Numbers measured on GitHub-hosted runners are not used to judge performance regressions.** The runners are shared VMs, and the CPU model, the load from neighbouring jobs and the turbo frequency change every time. This workflow has only two roles.
  - Check that the benchmark pipeline is not broken (smoke).
  - Look at rough trends.
  - Claims about reproducing the paper and regression judgements are made from results of `make benchmark` on a fixed, dedicated machine (follow the procedure of the `benchmark` skill).
- **How the suite is chosen**: push -> `smoke`, schedule -> `standard`, workflow_dispatch -> the `suite` input (`smoke` / `standard` / `full`).
- **Command that runs**: it calls `./scripts/run_benchmarks.sh <suite>` directly. `--allow-missing-runner` is added only for smoke.
  - The script configures and builds `release` and `release-diag` itself (narrow them down with `BENCH_BUILDS`).
  - If the benchmark runner (`build/<preset>/sssp_bench_runner`) does not exist yet, smoke only checks the build, prints "measurement skipped", and exits 0. CI detects this and writes "not implemented" in the Summary.
  - The `benchmark-smoke` target of the Makefile does not pass this flag, so CI calls the script directly.
  - Exit codes: 0 is normal (timeouts and out-of-memory are only recorded), 1 is an abnormal job exit, a schema violation or a build failure, and 2 is a bad argument or a missing runner (standard / full). Until the runner exists, the weekly standard run fails with exit 2.
- **When measurements ran**, `make report` (`scripts/analyze_results.py`) follows and produces `results/summary/*.csv` and `docs/reports/figures/*.png`. When there is no data, it prints "no data" and exits 0.
- **Inputs of a manual run**: `runner` (`ubuntu-26.04` / `macos-26`), `bench_timeout` (seconds; `BENCH_TIMEOUT`), `algorithms` / `graph_types` / `sizes` / `builds` / `params` / `order_seed` (`BENCH_ALGORITHMS` / `BENCH_GRAPH_TYPES` / `BENCH_SIZES` / `BENCH_BUILDS` / `BENCH_PARAMS` / `BENCH_ORDER_SEED` respectively).
  - Empty inputs are not put into the environment, and the script's defaults apply. The specification is in `docs/design/benchmark_runner_cli.md` section 5.
  - Inputs are always passed through `env:` and never expanded directly into `run:` (protection against shell injection). Before use they are validated to be integers, or to contain only allowed characters (alphanumerics and `_.,=:+-`).
- **When `scripts/run_benchmarks.sh` itself does not exist**, the job does not fail either. It writes "not implemented" in the Summary and finishes successfully.
- **Artifacts** (kept 30 days): `results/raw/` (`<run_id>.jsonl`, `.env.json`, `.env-end.json`, `.plan.tsv`), `results/**/*.csv`, `results/**/*.jsonl`, `docs/reports/figures/*.png`. The Summary lists the CSV / JSONL files generated during the run (row count, bytes) and shows the first 30 lines of the CSVs.
- Only for `full`, `timeout-minutes` is 350 (the limit for GitHub-hosted runners is 360 minutes). Otherwise it is 90 minutes.

### release.yml

1. Validate that the tag name has the form `vX.Y.Z[-suffix]`.
2. On Ubuntu (`linux-x86_64`) and macOS (`macos-arm64`), run `ci-release` in the order configure, build, ctest, pytest.
3. From the names in `tools/*.cpp`, collect the binaries `build/ci-release/<tool>`. Include `README.md`, `LICENSE` and `CITATION.cff` if present, and also `BUILD_INFO.txt` (tag, commit, compiler). Pack these into `almost-linear-sssp-<tag>-<platform>.tar.gz`.
4. Create a `.sha256` for each tar.gz. In the publish job, verify the overall `SHA256SUMS` with `sha256sum --check --strict`.
5. Create the Release with `gh release create --verify-tag --generate-notes`. If a Release for the same tag already exists, replace the assets with `gh release upload --clobber`. Tags containing `-` (such as `v1.0.0-rc1`) become pre-releases.

- A half-finished Release is a problem, so `concurrency` uses `cancel-in-progress: false` (this is the only exception among the four).
- The cache uses `actions/cache/restore` (read-only), and nothing is written to the cache from a tag.
- **Note**: the `ci-release` preset does not set `SSSP_DIAGNOSTICS`, so it is built with the default ON (with diagnostic counters). If you want a distribution without diagnostics, add a `release`-equivalent preset with `-Werror` and switch to it.

## Common policies

### Security

- **Every action is pinned to a commit SHA**, with a `# vX.Y.Z` comment (see the table below). SHAs were obtained with `gh api repos/<owner>/<repo>/git/ref/tags/<tag>`. When the object type was `tag` (an annotated tag), we followed `git/tags/<sha>`, and then confirmed with `gh api repos/<owner>/<repo>/commits/<tag>` that the result matches the commit SHA.
- **`permissions:`** defaults to `contents: read` for each workflow. The only write permission is `contents: write` in the `publish` job of `release.yml`.
- **`pull_request_target` is not used.** An untrusted PR from a fork runs under `pull_request`, so no Secrets are passed and the token is read-only.
- `actions/checkout` uses `persist-credentials: false` (so no token is left for later steps).
- Inputs of `workflow_dispatch` and tag names are passed through `env:` and their format is validated before use.

### Pinned actions (confirmed on 2026-11-27)

| Action | Version | Commit SHA |
|---|---|---|
| `actions/checkout` | v7.0.1 | `3d3c42e5aac5ba805825da76410c181273ba90b1` |
| `actions/cache` (including `/restore`) | v6.1.0 | `55cc8345863c7cc4c66a329aec7e433d2d1c52a9` |
| `actions/upload-artifact` | v7.0.2 | `cf430e030ddbb5b0abf93d22962f4752f3646cd9` |
| `actions/download-artifact` | v8.0.2 | `9000827ccba6bdab643e8b6fd33ac0654aef8333` |
| `astral-sh/setup-uv` | v10.4.0 | `1db35dd8eee11c5894872e0b383ae6c406cf738d` |

All five tags were lightweight tags (the object of `git/ref/tags` was the commit itself). When updating, fetch the SHA again with the same procedure and keep the version in the comment in sync. Adding Dependabot (`package-ecosystem: github-actions`) lets you receive update PRs while keeping the SHA pins.

### Runner labels

| Purpose | Label | Contents |
|---|---|---|
| macOS | `macos-26` | macOS 26, arm64 (Apple Silicon). `macos-latest` currently points to the same image, but we pin the version so that it does not switch on its own |
| Linux | `ubuntu-26.04` | Ubuntu 26.04 LTS, x64. The latest LTS. Where `ubuntu-latest` points to is in transition and differs between sources (see below), so we pin the version |

Sources checked:

- GitHub Docs "GitHub-hosted runners reference": https://docs.github.com/en/actions/reference/runners/github-hosted-runners
  - Source of the table: https://github.com/github/docs/blob/main/data/reusables/actions/supported-github-runners.md
  - As standard runners for public repositories, it lists `ubuntu-latest`, `ubuntu-24.04`, `ubuntu-22.04`, `ubuntu-26.04` for Linux x64 and `macos-latest`, `macos-15`, `macos-26` for macOS arm64.
  - In the Docs, the link target of `ubuntu-latest` is Ubuntu 24.04, while in the README of runner-images it is 26.04, so they disagree (as of 2026-11-27).
- README of actions/runner-images (mapping table between images and labels): https://github.com/actions/runner-images
  - It lists `ubuntu-latest` or `ubuntu-26.04` (x64) and `macos-latest`, `macos-26` or `macos-26-xlarge` (arm64).
  - The Intel macOS labels are `macos-26-intel` / `macos-15-intel`.
- Tool versions in the Ubuntu 26.04 image: https://github.com/actions/runner-images/blob/main/images/ubuntu/Ubuntu2604-Readme.md
  - Clang / clang-format / clang-tidy: 20.1.8, 21.1.8, 22.1.2. CMake 4.5.0. Ninja 1.13.2.
- macOS 26 arm64 image: https://github.com/actions/runner-images/blob/main/images/macos/macos-26-arm64-Readme.md
  - Apple clang 21.0.0 (Xcode). Homebrew llvm@20. CMake 4.5.0. Ninja 1.13.2.

### Compilers

- **Ubuntu**: `clang-22` / `clang-format-22` / `clang-tidy-22` from apt. The major version is set explicitly with `LLVM_MAJOR` in the workflow.
  - If the image already has them, use those, and run `apt-get install` only when they are missing.
  - The output of clang-format changes with the version. We use 22 to match the major version of the local one (Homebrew llvm 22.x). When the major version of Homebrew llvm goes up, raise `LLVM_MAJOR` together with it.
- **macOS**: Apple clang from Xcode (the default `c++`). Run `brew install ninja` only when the image does not have ninja.

### Cache

- **FetchContent** (Catch2, Google Benchmark): `cmake --preset <p> -DFETCHCONTENT_BASE_DIR=$GITHUB_WORKSPACE/.cache/fetchcontent/<p>` puts the download location in a per-preset directory, and it is saved with `actions/cache`.
  - The key is `OS-arch-preset-hash(cmake/Dependencies.cmake)`. When a dependency version or SHA256 changes, the cache switches automatically.
  - `.cache/` is already in `.gitignore` and excluded from ruff.
  - Only `benchmark.yml` cannot pass `-D`, because `run_benchmarks.sh` configures by itself. Instead it caches `build/*/_deps/*-src` and `*-subbuild` (best effort).
- **uv**: `enable-cache: true` and `cache-dependency-glob: uv.lock` of `astral-sh/setup-uv`. `uv sync --locked` fails if the lock is out of sync. Python is pinned with `UV_PYTHON=3.12`.

### concurrency, timeout, paths filter

- `concurrency.group = <workflow>-<ref>` with `cancel-in-progress: true`. When a new push arrives on the same branch or PR, the old run is stopped. Only release is an exception (see above).
- Every job has `timeout-minutes` (lint 25, build 45, sanitizer 60, release 60, benchmark 90 / full 350).
- `ci.yml` and `sanitizer.yml` use `paths-ignore: docs/**, **/*.md, LICENSE, .gitignore`. A documentation-only change does not trigger a build.
  - **Note**: for a PR where the workflow did not start because of the paths filter, the required check stays "pending". If you make `ci result` a required check in branch protection, either configure "do not require the check for docs-only PRs" or move the paths filter into a decision inside the jobs.
- A push to `benchmark.yml` starts it only when code changes (`src/`, `include/`, `tools/`, `benchmarks/`, `scripts/`, CMake, the uv lock).

## Local scripts

### scripts/bootstrap.sh

- It uses `uname -s` / `uname -m` to detect macOS arm64 (Homebrew at `/opt/homebrew`), macOS x86_64 (`/usr/local`; a warning under Rosetta) and Linux.
- macOS: it checks `cmake`, `ninja`, `llvm` and `uv`, and runs `brew install` only for what is missing.
  - Tools installed by other means (for example uv installed to `~/.local/bin` with the official installer) are accepted if they are found on PATH.
  - Homebrew llvm is keg-only and is not on PATH. The Makefile and `coverage.sh` look for it by path.
- Linux: it checks `cmake`, `ninja-build`, `clang-22`, `clang-format-22`, `clang-tidy-22` and `llvm-22` with apt, and installs what is missing. If a versioned package does not exist, it switches to the unversioned `clang` and so on.
  - uv is not in apt. It does not install it by piping a remote script into the shell on its own, and only points to the URL of the procedure.
- With `--check`, it only checks and installs nothing (exit code 1 if something is missing).
- It is idempotent: running it many times leaves the same state. At the end it runs `uv sync`.

### scripts/coverage.sh (`make coverage`)

1. Configure and build the `coverage` preset, then run ctest with `LLVM_PROFILE_FILE=build/coverage/coverage/profraw/%p-%m.profraw`.
2. Merge with `llvm-profdata merge -sparse`, and write the result of `llvm-cov report` to `build/coverage/coverage/summary.txt`. `--html` also produces HTML, and `--lcov` also produces lcov.
3. `_deps/`, `tests/` and `benchmarks/` are excluded from the totals.
4. LLVM tools are searched for in this order.
   - `$LLVM_PROFDATA` / `$LLVM_COV`
   - The same directory as the compiler in CMakeCache
   - For Apple clang, `xcrun --find`
   - On Linux, `llvm-*-<major>` or `/usr/lib/llvm-<major>/bin`
   - Homebrew llvm
   - PATH
   - The `.profraw` format differs between LLVM versions, so you must use tools of the same version as the compiler.

## Verifying with actionlint

```bash
actionlint .github/workflows/*.yml
```

- actionlint 1.7.12 (the latest version as of 2026-11-27) does not yet know the `ubuntu-26.04` label and reports a `runner-label` error. `ubuntu-26.04` is listed as GA in the runner-images table, so this is waiting for an update on the actionlint side.
- Until then, put the following settings in a temporary file and pass it with `-config-file` to verify. Do not put it in the repository.

```yaml
self-hosted-runner:
  labels: [ubuntu-26.04]
```

- If shellcheck is installed, actionlint also checks the contents of `run:`.
