# CI 設計（GitHub Actions）

この文書は `.github/workflows/` の4本の workflow と、`scripts/bootstrap.sh`・`scripts/coverage.sh` の設計をまとめる。
ローカルの品質ゲートは `make verify` で、CI はそれと同じコマンド（CMake preset、ctest、ruff、pytest）をなぞる。

## Workflow 一覧

| ファイル | 起動条件 | 内容 |
|---|---|---|
| `ci.yml` | `push`（main）、`pull_request`、手動 | lint（clang-format / clang-tidy / ruff）、`ci`（Debug, `-Werror`）と `ci-release`（Release, `-Werror`）のビルド、CTest をラベル別、pytest |
| `sanitizer.yml` | `push`（main）、`pull_request`、週1回、手動 | `sanitizer` preset（ASan + UBSan）。Linux では LeakSanitizer も有効 |
| `benchmark.yml` | `push`（main、コード変更時）、手動（入力あり）、週1回 | push は smoke、schedule は standard、手動は suite を選択。CSV/JSONL を artifact に保存し Summary に要約 |
| `release.yml` | `v*` タグの push | テスト → `ci-release` ビルド → `tools/` のバイナリを tar.gz 化 → SHA-256 → GitHub Release に添付 |

### ci.yml

- **matrix**: `os = [macos-26, ubuntu-26.04]` × `preset = [ci, ci-release]`（`fail-fast: false`）。
- **CTest はラベルごとに別ステップ**: `unit`、`integration`、`differential`、`property`、`regression`。
  - `-L '^unit$'` のように完全一致の正規表現で指定する。
  - `--no-tests=ignore` を付けるので、テストが0件のラベルでは失敗しない。
  - 最後に、どのラベルも付いていないテストを拾う `-LE` のステップがある（付け忘れの安全網）。
  - 1つのラベルが落ちても、残りのラベルは実行する（`!cancelled() && steps.build.outcome == 'success'`）。
- **pytest** は `ci` preset のジョブでだけ回す（`tests/harness`、`tests/python`）。テストが0件のとき（終了コード 5）は成功扱いにする。
- **lint** は Ubuntu だけで回す。
  - clang-format は `git ls-files` の C++ ファイルに `--dry-run --Werror`。
  - clang-tidy は `make tidy PRESET=ci`（`src/`、`tools/` が対象）。`.clang-tidy` の `WarningsAsErrors` で失敗させる。
  - ruff は `check` と `format --check`。
- **`ci result` ジョブ**: すべてのジョブの結果を1つにまとめる。ブランチ保護の必須チェックには、matrix で名前が変わるジョブではなく、これを指定する。
- 失敗時は `build/<preset>/Testing/Temporary/` を artifact として保存する。

### sanitizer.yml

- `ASAN_OPTIONS` は OS ごとに分ける。
  - Linux: `detect_leaks=1`（LeakSanitizer 有効）、`strict_string_checks`、`detect_stack_use_after_return`、`check_initialization_order` など。
  - macOS: `detect_leaks=0`。Apple clang の ASan は macOS で LeakSanitizer に対応していない。リークの検出は Linux のジョブが担う。
- `UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1`。preset 側でも `-fno-sanitize-recover=undefined` を指定しているので、UB は必ず失敗になる。
- Linux では `vm.mmap_rnd_bits=28` に下げる。新しいカーネルでは ASLR のエントロピーが大きく、ASan の shadow mapping が壊れる既知の問題がある（clang 22 では起きにくいが、念のため）。

### benchmark.yml

- **GitHub-hosted runner で測った値は、性能回帰の判定に使わない。** runner は共有 VM で、CPU の型番、隣のジョブの負荷、ターボ周波数が毎回変わる。この workflow の役割は次の2つだけである。
  - ベンチマークのパイプラインが壊れていないことを確かめる（smoke）。
  - おおまかな傾向を見る。
  - 論文の再現に関する主張や回帰の判定は、固定した専用マシンで `make benchmark` を回した結果で行う（`benchmark` skill の手順に従う）。
- **suite の決まり方**: push → `smoke`、schedule → `standard`、workflow_dispatch → 入力 `suite`（`smoke` / `standard` / `full`）。
- **実行するコマンド**: `./scripts/run_benchmarks.sh <suite>` を直接呼ぶ。smoke のときだけ `--allow-missing-runner` を付ける。
  - スクリプトが自分で `release` と `release-diag` を configure・ビルドする（`BENCH_BUILDS` で絞れる）。
  - ベンチランナー（`build/<preset>/sssp_bench_runner`）がまだ無いと、smoke はビルドだけ確かめて「計測をスキップしました」と出し、exit 0 で終わる。CI はこれを検出して Summary に「未実装」と書く。
  - Makefile の `benchmark-smoke` はこのフラグを付けないので、CI ではスクリプトを直接呼んでいる。
  - 終了コード: 0 は正常（timeout とメモリ超過は記録するだけ）、1 はジョブの異常終了・スキーマ違反・ビルド失敗、2 は引数の誤りかランナー不在（standard / full）。ランナーができるまでは、週1回の standard は exit 2 で失敗する。
- **計測が走ったとき**は、続けて `make report`（`scripts/analyze_results.py`）で `results/summary/*.csv` と `docs/reports/figures/*.png` を作る。データが無いときは「データなし」と出して exit 0 で終わる。
- **手動実行の入力**: `runner`（`ubuntu-26.04` / `macos-26`）、`bench_timeout`（秒。`BENCH_TIMEOUT`）、`algorithms` / `graph_types` / `sizes` / `builds` / `params` / `order_seed`（それぞれ `BENCH_ALGORITHMS` / `BENCH_GRAPH_TYPES` / `BENCH_SIZES` / `BENCH_BUILDS` / `BENCH_PARAMS` / `BENCH_ORDER_SEED`）。
  - 空の入力は環境変数に入れず、スクリプトの既定値に任せる。仕様は `docs/design/benchmark_runner_cli.md` §5 にある。
  - 入力は必ず `env:` 経由で渡し、`run:` に直接展開しない（シェルインジェクション対策）。整数かどうか、使える文字（英数字と `_.,=:+-`）だけかを検証してから使う。
- **`scripts/run_benchmarks.sh` 自体が無いとき**も失敗させない。Summary に「未実装」と書いて成功で終える。
- **artifact**（30日保存）: `results/raw/`（`<run_id>.jsonl`、`.env.json`、`.env-end.json`、`.plan.tsv`）、`results/**/*.csv`、`results/**/*.jsonl`、`docs/reports/figures/*.png`。Summary には、実行中に生成された CSV / JSONL の一覧（行数、バイト数）と、CSV の先頭30行を載せる。
- `full` のときだけ `timeout-minutes` を 350 にする（GitHub-hosted の上限は 360 分）。それ以外は 90 分。

### release.yml

1. タグ名を `vX.Y.Z[-suffix]` の形か検証する。
2. Ubuntu（`linux-x86_64`）と macOS（`macos-arm64`）で、`ci-release` を configure、ビルド、ctest、pytest の順に回す。
3. `tools/*.cpp` の名前から、`build/ci-release/<tool>` のバイナリを集める。`README.md`、`LICENSE`、`CITATION.cff` があれば同梱し、`BUILD_INFO.txt`（タグ、コミット、コンパイラ）も入れる。これを `almost-linear-sssp-<tag>-<platform>.tar.gz` にまとめる。
4. 各 tar.gz の `.sha256` を作る。publish ジョブで、全体の `SHA256SUMS` を `sha256sum --check --strict` で検証する。
5. `gh release create --verify-tag --generate-notes` で Release を作る。同じタグの Release がすでにあれば `gh release upload --clobber` で差し替える。`-` を含むタグ（`v1.0.0-rc1` など）は pre-release にする。

- 途中で止まった Release は困るので、`concurrency` は `cancel-in-progress: false` にしている（4本のうち、これだけが例外）。
- キャッシュは `actions/cache/restore`（読むだけ）で、タグからはキャッシュに書き込まない。
- **注意**: `ci-release` preset は `SSSP_DIAGNOSTICS` を指定していないので、既定の ON（診断カウンタ入り）でビルドされる。診断を抜いた配布物にしたいなら、`-Werror` を付けた `release` 相当の preset を追加して、それに切り替える。

## 共通の方針

### セキュリティ

- **Action はすべてコミット SHA で固定**し、`# vX.Y.Z` のコメントを付ける（下表）。SHA は `gh api repos/<owner>/<repo>/git/ref/tags/<tag>` で取得した。object の type が `tag`（annotated tag）なら `git/tags/<sha>` を辿り、さらに `gh api repos/<owner>/<repo>/commits/<tag>` で commit SHA と一致することを確かめた。
- **`permissions:`** は workflow の既定を `contents: read` にする。書き込み権限は、`release.yml` の `publish` ジョブの `contents: write` だけ。
- **`pull_request_target` は使わない。** fork からの未信頼な PR は `pull_request` で走るので、Secrets は渡らず、トークンは読み取り専用になる。
- `actions/checkout` は `persist-credentials: false` にする（後のステップにトークンを残さない）。
- `workflow_dispatch` の入力とタグ名は、`env:` 経由で渡し、形式を検証してから使う。

### 固定した Action（2026-11-27 に確認）

| Action | バージョン | コミット SHA |
|---|---|---|
| `actions/checkout` | v7.0.1 | `3d3c42e5aac5ba805825da76410c181273ba90b1` |
| `actions/cache`（`/restore` を含む） | v6.1.0 | `55cc8345863c7cc4c66a329aec7e433d2d1c52a9` |
| `actions/upload-artifact` | v7.0.2 | `cf430e030ddbb5b0abf93d22962f4752f3646cd9` |
| `actions/download-artifact` | v8.0.2 | `9000827ccba6bdab643e8b6fd33ac0654aef8333` |
| `astral-sh/setup-uv` | v10.4.0 | `1db35dd8eee11c5894872e0b383ae6c406cf738d` |

5つともタグは lightweight tag だった（`git/ref/tags` の object がそのまま commit）。更新するときは、同じ手順で SHA を取り直し、コメントのバージョンもそろえる。Dependabot（`package-ecosystem: github-actions`）を入れると、SHA 固定のまま更新 PR を受け取れる。

### Runner ラベル

| 用途 | ラベル | 中身 |
|---|---|---|
| macOS | `macos-26` | macOS 26、arm64（Apple Silicon）。`macos-latest` も現在は同じイメージだが、勝手に切り替わらないようにバージョンで固定する |
| Linux | `ubuntu-26.04` | Ubuntu 26.04 LTS、x64。最新の LTS。`ubuntu-latest` の行き先は移行中で資料によって違う（下記）ため、バージョンで固定する |

確認元:

- GitHub Docs「GitHub-hosted runners reference」: https://docs.github.com/en/actions/reference/runners/github-hosted-runners
  - 表の原稿: https://github.com/github/docs/blob/main/data/reusables/actions/supported-github-runners.md
  - 公開リポジトリ向けの標準 runner として、Linux x64 に `ubuntu-latest`, `ubuntu-24.04`, `ubuntu-22.04`, `ubuntu-26.04`、macOS arm64 に `macos-latest`, `macos-15`, `macos-26` が載っている。
  - Docs では `ubuntu-latest` のリンク先が Ubuntu 24.04、runner-images の README では 26.04 で、食い違っている（2026-11-27 時点）。
- actions/runner-images の README（イメージとラベルの対応表）: https://github.com/actions/runner-images
  - `ubuntu-latest` or `ubuntu-26.04`（x64）、`macos-latest`, `macos-26` or `macos-26-xlarge`（arm64）と記載されている。
  - Intel の macOS は `macos-26-intel` / `macos-15-intel`。
- Ubuntu 26.04 イメージに入っているツールの版: https://github.com/actions/runner-images/blob/main/images/ubuntu/Ubuntu2604-Readme.md
  - Clang / clang-format / clang-tidy: 20.1.8、21.1.8、22.1.2。CMake 4.5.0。Ninja 1.13.2。
- macOS 26 arm64 イメージ: https://github.com/actions/runner-images/blob/main/images/macos/macos-26-arm64-Readme.md
  - Apple clang 21.0.0（Xcode）。Homebrew の llvm@20。CMake 4.5.0。Ninja 1.13.2。

### コンパイラ

- **Ubuntu**: apt の `clang-22` / `clang-format-22` / `clang-tidy-22`。メジャー版は workflow の `LLVM_MAJOR` で明示する。
  - イメージに入っていればそれを使い、無いときだけ `apt-get install` する。
  - clang-format は版によって整形結果が変わる。ローカル（Homebrew llvm 22.x）とメジャー版をそろえるために 22 にしている。Homebrew の llvm のメジャー版が上がったら、`LLVM_MAJOR` も一緒に上げる。
- **macOS**: Xcode の Apple clang（既定の `c++`）。ninja がイメージに無いときだけ `brew install ninja` する。

### キャッシュ

- **FetchContent**（Catch2、Google Benchmark）: `cmake --preset <p> -DFETCHCONTENT_BASE_DIR=$GITHUB_WORKSPACE/.cache/fetchcontent/<p>` で、ダウンロード先を preset ごとのディレクトリにし、それを `actions/cache` で保存する。
  - キーは `OS-arch-preset-hash(cmake/Dependencies.cmake)`。依存のバージョンや SHA256 を変えると、キャッシュも自動で切り替わる。
  - `.cache/` は `.gitignore` 済みで、ruff の対象外。
  - `benchmark.yml` だけは、`run_benchmarks.sh` が自分で configure するので `-D` を渡せない。代わりに `build/*/_deps/*-src` と `*-subbuild` をキャッシュする（ベストエフォート）。
- **uv**: `astral-sh/setup-uv` の `enable-cache: true`、`cache-dependency-glob: uv.lock`。`uv sync --locked` で、lock とずれていたら失敗させる。Python は `UV_PYTHON=3.12` に固定する。

### concurrency、timeout、paths filter

- `concurrency.group = <workflow>-<ref>`、`cancel-in-progress: true`。同じブランチや PR に新しい push が来たら、古い実行を止める。release だけは例外（上記）。
- すべてのジョブに `timeout-minutes` を付ける（lint 25、ビルド 45、sanitizer 60、release 60、benchmark 90 / full 350）。
- `ci.yml` と `sanitizer.yml` は `paths-ignore: docs/**, **/*.md, LICENSE, .gitignore`。ドキュメントだけの変更ではビルドしない。
  - **注意**: paths filter で workflow が起動しなかった PR では、必須チェックが「保留」のままになる。ブランチ保護で `ci result` を必須にするなら、「docs だけの PR はチェックを要求しない」設定にするか、paths filter をジョブ側の判定に移す。
- `benchmark.yml` の push は、コード（`src/`、`include/`、`tools/`、`benchmarks/`、`scripts/`、CMake、uv の lock）が変わったときだけ起動する。

## ローカル用スクリプト

### scripts/bootstrap.sh

- `uname -s` / `uname -m` で、macOS arm64（Homebrew は `/opt/homebrew`）、macOS x86_64（`/usr/local`。Rosetta 上なら警告）、Linux を判定する。
- macOS: `cmake`、`ninja`、`llvm`、`uv` を確認し、無いものだけ `brew install` する。
  - 別の方法で入っているもの（例: 公式インストーラで `~/.local/bin` に入れた uv）は、PATH で見つかればそれで良しとする。
  - Homebrew の llvm は keg-only で PATH に載らない。Makefile と `coverage.sh` はパスで探す。
- Linux: apt で `cmake`、`ninja-build`、`clang-22`、`clang-format-22`、`clang-tidy-22`、`llvm-22` を確認し、無いものを入れる。版つきのパッケージが無いときは、版なしの `clang` などに切り替える。
  - uv は apt に無い。リモートのスクリプトをシェルに流し込む導入は勝手にせず、手順の URL を案内するだけにする。
- `--check` を付けると、確認だけして何も入れない（足りなければ終了コード 1）。
- 冪等で、何度実行しても同じ状態になる。最後に `uv sync` を実行する。

### scripts/coverage.sh（`make coverage`）

1. `coverage` preset を configure してビルドし、`LLVM_PROFILE_FILE=build/coverage/coverage/profraw/%p-%m.profraw` を付けて ctest を回す。
2. `llvm-profdata merge -sparse` でまとめ、`llvm-cov report` の結果を `build/coverage/coverage/summary.txt` に書く。`--html` を付けると HTML、`--lcov` を付けると lcov も出す。
3. `_deps/`、`tests/`、`benchmarks/` は集計から外す。
4. LLVM のツールは次の順で探す。
   - `$LLVM_PROFDATA` / `$LLVM_COV`
   - CMakeCache のコンパイラと同じディレクトリ
   - Apple clang なら `xcrun --find`
   - Linux なら `llvm-*-<major>` または `/usr/lib/llvm-<major>/bin`
   - Homebrew の llvm
   - PATH
   - `.profraw` の形式は LLVM の版ごとに違うので、コンパイラと同じ版のツールを使う必要がある。

## actionlint での検証

```bash
actionlint .github/workflows/*.yml
```

- actionlint 1.7.12（2026-11-27 時点の最新版）は、`ubuntu-26.04` ラベルをまだ知らず、`runner-label` エラーを出す。`ubuntu-26.04` は runner-images の表に GA として載っているので、これは actionlint 側の更新待ちである。
- それまでは、次の設定を一時ファイルに置き、`-config-file` で渡して検証する。リポジトリには置かない。

```yaml
self-hosted-runner:
  labels: [ubuntu-26.04]
```

- shellcheck が入っていれば、actionlint は `run:` の中身も検査する。
