# CI 设计(GitHub Actions)

[English](ci.md) | [日本語](ci.ja.md) | **简体中文**

本文档说明 `.github/workflows/` 下 4 个 workflow 的设计,以及 `scripts/bootstrap.sh` 和 `scripts/coverage.sh` 的设计。
本地的质量门禁是 `make verify`,CI 沿用与它相同的命令(CMake preset、ctest、ruff、pytest)。

## Workflow 一览

| 文件 | 触发条件 | 内容 |
|---|---|---|
| `ci.yml` | `push`(main)、`pull_request`、手动 | lint(clang-format / clang-tidy / ruff)、`ci`(Debug,`-Werror`)与 `ci-release`(Release,`-Werror`)的构建、按标签运行 CTest、pytest |
| `sanitizer.yml` | `push`(main)、`pull_request`、每周一次、手动 | `sanitizer` preset(ASan + UBSan)。在 Linux 上同时启用 LeakSanitizer |
| `benchmark.yml` | `push`(main,代码变更时)、手动(带输入)、每周一次 | push 运行 smoke,schedule 运行 standard,手动可选择 suite。把 CSV/JSONL 保存为 artifact,并在 Summary 中汇总 |
| `release.yml` | 推送 `v*` 标签 | 测试 -> 构建 `ci-release` -> 把 `tools/` 的二进制打成 tar.gz -> SHA-256 -> 附加到 GitHub Release |

### ci.yml

- **matrix**:`os = [macos-26, ubuntu-26.04]` x `preset = [ci, ci-release]`(`fail-fast: false`)。
- **CTest 按标签分成独立的步骤**:`unit`、`integration`、`differential`、`property`、`regression`。
  - 用 `-L '^unit$'` 这样的完全匹配正则来指定。
  - 附带 `--no-tests=ignore`,所以测试数为 0 的标签不会失败。
  - 最后有一个用 `-LE` 的步骤,捞出没有任何标签的测试(防止忘记加标签的保险)。
  - 即使某个标签失败,其余标签仍会运行(`!cancelled() && steps.build.outcome == 'success'`)。
- **pytest** 只在 `ci` preset 的 job 中运行(`tests/harness`、`tests/python`)。测试数为 0 时(退出码 5)视为成功。
- **lint** 只在 Ubuntu 上运行。
  - clang-format 对 `git ls-files` 中的 C++ 文件执行 `--dry-run --Werror`。
  - clang-tidy 执行 `make tidy PRESET=ci`(对象为 `src/`、`tools/`)。通过 `.clang-tidy` 的 `WarningsAsErrors` 使其失败。
  - ruff 执行 `check` 和 `format --check`。
- **`ci result` job**:把所有 job 的结果合并为一个。分支保护的必需检查应指定它,而不是名称会随 matrix 变化的 job。
- 失败时,把 `build/<preset>/Testing/Temporary/` 保存为 artifact。

### sanitizer.yml

- `ASAN_OPTIONS` 按操作系统区分。
  - Linux:`detect_leaks=1`(启用 LeakSanitizer)、`strict_string_checks`、`detect_stack_use_after_return`、`check_initialization_order` 等。
  - macOS:`detect_leaks=0`。Apple clang 的 ASan 在 macOS 上不支持 LeakSanitizer。泄漏检测由 Linux 的 job 负责。
- `UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1`。preset 一侧也指定了 `-fno-sanitize-recover=undefined`,所以出现 UB 一定会失败。
- 在 Linux 上把 `vm.mmap_rnd_bits` 降到 28。较新的内核中 ASLR 的熵更大,存在会破坏 ASan shadow mapping 的已知问题(clang 22 不太会出现,这里只是以防万一)。

### benchmark.yml

- **在 GitHub-hosted runner 上测得的数值,不用于判断性能回归。** runner 是共享的虚拟机,CPU 型号、相邻 job 的负载、睿频每次都不同。这个 workflow 的职责只有以下两项。
  - 确认基准测试的流水线没有损坏(smoke)。
  - 观察大致的趋势。
  - 关于论文复现的主张和回归的判定,依据在固定的专用机器上运行 `make benchmark` 得到的结果(遵循 `benchmark` skill 的步骤)。
- **suite 的决定方式**:push -> `smoke`,schedule -> `standard`,workflow_dispatch -> 输入 `suite`(`smoke` / `standard` / `full`)。
- **运行的命令**:直接调用 `./scripts/run_benchmarks.sh <suite>`。只有 smoke 时才附加 `--allow-missing-runner`。
  - 脚本自己会 configure 并构建 `release` 和 `release-diag`(可用 `BENCH_BUILDS` 缩小范围)。
  - 如果基准测试运行器(`build/<preset>/sssp_bench_runner`)还不存在,smoke 只确认构建,输出"已跳过测量",并以 exit 0 结束。CI 会检测这一点,并在 Summary 中写上"未实现"。
  - Makefile 的 `benchmark-smoke` 不附加这个标志,所以 CI 直接调用脚本。
  - 退出码:0 为正常(超时和内存超限只做记录),1 为 job 异常终止、schema 违规或构建失败,2 为参数错误或运行器不存在(standard / full)。在运行器完成之前,每周一次的 standard 会以 exit 2 失败。
- **运行了测量时**,接着用 `make report`(`scripts/analyze_results.py`)生成 `results/summary/*.csv` 和 `docs/reports/figures/*.png`。没有数据时输出"无数据"并以 exit 0 结束。
- **手动运行的输入**:`runner`(`ubuntu-26.04` / `macos-26`)、`bench_timeout`(秒;`BENCH_TIMEOUT`)、`algorithms` / `graph_types` / `sizes` / `builds` / `params` / `order_seed`(分别对应 `BENCH_ALGORITHMS` / `BENCH_GRAPH_TYPES` / `BENCH_SIZES` / `BENCH_BUILDS` / `BENCH_PARAMS` / `BENCH_ORDER_SEED`)。
  - 空的输入不放入环境变量,交给脚本的默认值处理。规范见 `docs/design/benchmark_runner_cli.md` 第 5 节。
  - 输入一律经由 `env:` 传递,不直接展开到 `run:` 中(防止 shell 注入)。使用前会验证是否为整数,以及是否只含允许的字符(字母数字和 `_.,=:+-`)。
- **`scripts/run_benchmarks.sh` 本身不存在时**也不会失败。在 Summary 中写上"未实现"并以成功结束。
- **artifact**(保存 30 天):`results/raw/`(`<run_id>.jsonl`、`.env.json`、`.env-end.json`、`.plan.tsv`)、`results/**/*.csv`、`results/**/*.jsonl`、`docs/reports/figures/*.png`。Summary 中列出运行期间生成的 CSV / JSONL 文件(行数、字节数),并展示 CSV 的前 30 行。
- 只有 `full` 时 `timeout-minutes` 为 350(GitHub-hosted 的上限是 360 分钟)。其余为 90 分钟。

### release.yml

1. 验证标签名是否为 `vX.Y.Z[-suffix]` 的形式。
2. 在 Ubuntu(`linux-x86_64`)和 macOS(`macos-arm64`)上,按 configure、构建、ctest、pytest 的顺序运行 `ci-release`。
3. 根据 `tools/*.cpp` 的名称,收集 `build/ci-release/<tool>` 的二进制。如果有 `README.md`、`LICENSE`、`CITATION.cff` 则一并打入,并加入 `BUILD_INFO.txt`(标签、提交、编译器)。把这些打包成 `almost-linear-sssp-<tag>-<platform>.tar.gz`。
4. 为每个 tar.gz 生成 `.sha256`。在 publish job 中,用 `sha256sum --check --strict` 校验整体的 `SHA256SUMS`。
5. 用 `gh release create --verify-tag --generate-notes` 创建 Release。如果同一标签的 Release 已存在,则用 `gh release upload --clobber` 替换。含有 `-` 的标签(如 `v1.0.0-rc1`)会成为 pre-release。

- 中途停止的 Release 会带来麻烦,所以 `concurrency` 使用 `cancel-in-progress: false`(4 个 workflow 中只有它是例外)。
- 缓存使用 `actions/cache/restore`(只读),不会从标签向缓存写入。
- **注意**:`ci-release` preset 没有指定 `SSSP_DIAGNOSTICS`,所以按默认的 ON(带诊断计数器)构建。如果想要不含诊断的发布产物,请添加一个附带 `-Werror` 的、相当于 `release` 的 preset 并切换过去。

## 通用方针

### 安全

- **所有 Action 都固定到提交 SHA**,并附上 `# vX.Y.Z` 注释(见下表)。SHA 是用 `gh api repos/<owner>/<repo>/git/ref/tags/<tag>` 获取的。如果 object 的 type 是 `tag`(annotated tag),就沿着 `git/tags/<sha>` 继续追踪,并用 `gh api repos/<owner>/<repo>/commits/<tag>` 确认其与提交 SHA 一致。
- **`permissions:`** 的 workflow 默认值为 `contents: read`。写权限只有 `release.yml` 的 `publish` job 中的 `contents: write`。
- **不使用 `pull_request_target`。** 来自 fork 的不受信任的 PR 通过 `pull_request` 运行,因此不会传入 Secrets,令牌为只读。
- `actions/checkout` 使用 `persist-credentials: false`(不在后续步骤中留下令牌)。
- `workflow_dispatch` 的输入和标签名经由 `env:` 传递,验证格式后才使用。

### 已固定的 Action(2026-11-27 确认)

| Action | 版本 | 提交 SHA |
|---|---|---|
| `actions/checkout` | v7.0.1 | `3d3c42e5aac5ba805825da76410c181273ba90b1` |
| `actions/cache`(含 `/restore`) | v6.1.0 | `55cc8345863c7cc4c66a329aec7e433d2d1c52a9` |
| `actions/upload-artifact` | v7.0.2 | `cf430e030ddbb5b0abf93d22962f4752f3646cd9` |
| `actions/download-artifact` | v8.0.2 | `9000827ccba6bdab643e8b6fd33ac0654aef8333` |
| `astral-sh/setup-uv` | v10.4.0 | `1db35dd8eee11c5894872e0b383ae6c406cf738d` |

这 5 个的标签都是 lightweight tag(`git/ref/tags` 的 object 直接就是提交)。更新时,按同样的步骤重新获取 SHA,并使注释中的版本保持一致。引入 Dependabot(`package-ecosystem: github-actions`)后,可以在保持 SHA 固定的同时接收更新 PR。

### Runner 标签

| 用途 | 标签 | 内容 |
|---|---|---|
| macOS | `macos-26` | macOS 26,arm64(Apple Silicon)。`macos-latest` 目前也指向同一镜像,但为了避免它自行切换,这里按版本固定 |
| Linux | `ubuntu-26.04` | Ubuntu 26.04 LTS,x64。最新的 LTS。`ubuntu-latest` 的指向正在迁移中,各资料的说法不一(见下文),所以按版本固定 |

确认来源:

- GitHub Docs "GitHub-hosted runners reference":https://docs.github.com/en/actions/reference/runners/github-hosted-runners
  - 表格的原稿:https://github.com/github/docs/blob/main/data/reusables/actions/supported-github-runners.md
  - 作为公开仓库的标准 runner,Linux x64 列出了 `ubuntu-latest`、`ubuntu-24.04`、`ubuntu-22.04`、`ubuntu-26.04`,macOS arm64 列出了 `macos-latest`、`macos-15`、`macos-26`。
  - Docs 中 `ubuntu-latest` 的链接目标是 Ubuntu 24.04,而 runner-images 的 README 中是 26.04,两者不一致(截至 2026-11-27)。
- actions/runner-images 的 README(镜像与标签的对应表):https://github.com/actions/runner-images
  - 记载了 `ubuntu-latest` or `ubuntu-26.04`(x64)、`macos-latest`、`macos-26` or `macos-26-xlarge`(arm64)。
  - Intel 的 macOS 为 `macos-26-intel` / `macos-15-intel`。
- Ubuntu 26.04 镜像中安装的工具版本:https://github.com/actions/runner-images/blob/main/images/ubuntu/Ubuntu2604-Readme.md
  - Clang / clang-format / clang-tidy:20.1.8、21.1.8、22.1.2。CMake 4.5.0。Ninja 1.13.2。
- macOS 26 arm64 镜像:https://github.com/actions/runner-images/blob/main/images/macos/macos-26-arm64-Readme.md
  - Apple clang 21.0.0(Xcode)。Homebrew 的 llvm@20。CMake 4.5.0。Ninja 1.13.2。

### 编译器

- **Ubuntu**:apt 的 `clang-22` / `clang-format-22` / `clang-tidy-22`。主版本由 workflow 中的 `LLVM_MAJOR` 明确指定。
  - 镜像里已有就用镜像里的,只有缺失时才 `apt-get install`。
  - clang-format 的格式化结果会随版本变化。为了与本地(Homebrew llvm 22.x)的主版本保持一致,采用 22。Homebrew 的 llvm 主版本升级时,`LLVM_MAJOR` 也要一起升级。
- **macOS**:Xcode 的 Apple clang(默认的 `c++`)。只有镜像里没有 ninja 时才 `brew install ninja`。

### 缓存

- **FetchContent**(Catch2、Google Benchmark):用 `cmake --preset <p> -DFETCHCONTENT_BASE_DIR=$GITHUB_WORKSPACE/.cache/fetchcontent/<p>` 把下载位置设为每个 preset 独立的目录,并用 `actions/cache` 保存。
  - 键为 `OS-arch-preset-hash(cmake/Dependencies.cmake)`。更改依赖的版本或 SHA256 时,缓存也会自动切换。
  - `.cache/` 已在 `.gitignore` 中,并且不在 ruff 的检查范围内。
  - 只有 `benchmark.yml` 因为 `run_benchmarks.sh` 自己会 configure,无法传入 `-D`。作为替代,缓存 `build/*/_deps/*-src` 和 `*-subbuild`(尽力而为)。
- **uv**:`astral-sh/setup-uv` 的 `enable-cache: true`、`cache-dependency-glob: uv.lock`。用 `uv sync --locked`,与 lock 不一致就让它失败。Python 用 `UV_PYTHON=3.12` 固定。

### concurrency、timeout、paths filter

- `concurrency.group = <workflow>-<ref>`,`cancel-in-progress: true`。同一分支或 PR 上有新的 push 到来时,停止旧的运行。只有 release 是例外(见上文)。
- 所有 job 都设置 `timeout-minutes`(lint 25、构建 45、sanitizer 60、release 60、benchmark 90 / full 350)。
- `ci.yml` 和 `sanitizer.yml` 使用 `paths-ignore: docs/**, **/*.md, LICENSE, .gitignore`。只改动文档时不会构建。
  - **注意**:对于因 paths filter 而没有启动 workflow 的 PR,必需检查会一直处于"挂起"状态。如果要在分支保护中把 `ci result` 设为必需检查,就要配置成"只改文档的 PR 不要求该检查",或者把 paths filter 改为在 job 内部判断。
- `benchmark.yml` 的 push 只在代码(`src/`、`include/`、`tools/`、`benchmarks/`、`scripts/`、CMake、uv 的 lock)变更时才启动。

## 本地脚本

### scripts/bootstrap.sh

- 通过 `uname -s` / `uname -m` 判断 macOS arm64(Homebrew 在 `/opt/homebrew`)、macOS x86_64(`/usr/local`;在 Rosetta 上会警告)和 Linux。
- macOS:检查 `cmake`、`ninja`、`llvm`、`uv`,只对缺少的执行 `brew install`。
  - 以其他方式安装的(例如用官方安装程序装到 `~/.local/bin` 的 uv),只要能在 PATH 中找到就视为可以。
  - Homebrew 的 llvm 是 keg-only,不在 PATH 上。Makefile 和 `coverage.sh` 按路径查找。
- Linux:用 apt 检查 `cmake`、`ninja-build`、`clang-22`、`clang-format-22`、`clang-tidy-22`、`llvm-22`,并安装缺少的。如果没有带版本的软件包,就切换到不带版本的 `clang` 等。
  - uv 不在 apt 中。不会擅自采用把远程脚本管道给 shell 的安装方式,只提示步骤的 URL。
- 加上 `--check` 时,只检查而不安装任何东西(缺少则退出码为 1)。
- 具有幂等性,运行多少次都得到相同的状态。最后执行 `uv sync`。

### scripts/coverage.sh(`make coverage`)

1. 配置并构建 `coverage` preset,附带 `LLVM_PROFILE_FILE=build/coverage/coverage/profraw/%p-%m.profraw` 运行 ctest。
2. 用 `llvm-profdata merge -sparse` 合并,并把 `llvm-cov report` 的结果写入 `build/coverage/coverage/summary.txt`。加上 `--html` 还会输出 HTML,加上 `--lcov` 还会输出 lcov。
3. `_deps/`、`tests/`、`benchmarks/` 不计入统计。
4. LLVM 工具按以下顺序查找。
   - `$LLVM_PROFDATA` / `$LLVM_COV`
   - 与 CMakeCache 中的编译器相同的目录
   - 对于 Apple clang,用 `xcrun --find`
   - 对于 Linux,用 `llvm-*-<major>` 或 `/usr/lib/llvm-<major>/bin`
   - Homebrew 的 llvm
   - PATH
   - `.profraw` 的格式随 LLVM 版本而不同,所以必须使用与编译器相同版本的工具。

## 用 actionlint 验证

```bash
actionlint .github/workflows/*.yml
```

- actionlint 1.7.12(截至 2026-11-27 的最新版本)还不认识 `ubuntu-26.04` 标签,会报 `runner-label` 错误。`ubuntu-26.04` 已在 runner-images 的表中作为 GA 列出,所以这是在等待 actionlint 一侧的更新。
- 在此之前,把下面的配置放进临时文件,用 `-config-file` 传入进行验证。不要放进仓库。

```yaml
self-hosted-runner:
  labels: [ubuntu-26.04]
```

- 如果安装了 shellcheck,actionlint 也会检查 `run:` 中的内容。
