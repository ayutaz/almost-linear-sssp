# 先行研究与现有 OSS 调研

[English](prior_work_and_oss.md) | [日本語](prior_work_and_oss.ja.md) | **简体中文**

- 调研日期:2026-11-27
- 目标论文:Iman Gholami, MohammadTaghi Hajiaghayi, "Linear Single-Source Shortest Path, Almost", arXiv:2611.34911v1(2026-11-25 提交) <https://arxiv.org/abs/2611.34911>
  - 问题设定:无向多重图、非负实数边长,求从单一源点到所有顶点的精确距离以及最短路径树。比较-加法模型。期望时间 O(m + n log log log n)。几乎必然终止,且每次终止的运行都是正确的(论文 Theorem 1,正文 Section 1)。
  - 该论文在 arXiv 上的许可为 `arXiv.org perpetual, non-exclusive license 1.0`(abs 页面上的许可链接 <http://arxiv.org/licenses/nonexclusive-distrib/1.0/>)。该许可只授予 arXiv 分发权,因此本仓库不再分发论文的正文、图表和译文。`.cache/` 已在 `.gitignore` 中。
- 方针:事实须附来源。无法确认的内容标记为 "unverified"(未验证)。只写来源中出现的数字,不写推测的数字。
- 克隆下来的代码只阅读,未执行。克隆到本地临时目录(`git clone --depth 1`)。

---

## 0. 对本项目有直接影响的要点

1. **对任何物理上可能的输入,论文的主算法(Algorithm 8 `ShortestPaths`)都一定进入基本情形(base case)。**
   - Algorithm 8 第 4 行写的是 "if `N ≤ L_N^64`,每个顶点对保留一条最轻的边,并用基于扫描的 Dijkstra 搜索求解(Lemma 3.2)"。其中 L_N = ⌈log2(N+2)⌉(定义见论文正文 Section 2 及 Algorithm 8;`.cache/paper/2611.34911v1.txt` 的第 283–290 行和第 2167–2185 行)。
   - Section 2 正文还写道 "The outer and local constructions assume N > L_N^64 and N > 2 L_N^4, respectively. The excluded sizes are absolutely bounded."
   - 求满足 N > L_N^64 的最小 L_N,结果是 589。也就是说,只有当 N ≳ 2^589 时才会进入非平凡路径(此计算由本次调研完成:用 Python 在 L = 2..2000 中扫描,找第一个满足 2^L − 2 > L^64 的 L)。精确值见 [thresholds.md](thresholds.md):最小的 N 是 589^64 + 1 ≈ 1.94×10^177,该判定在 [2^589 − 1, 590^64] 上再次为假,而对所有 N ≥ 590^64 + 1 成立。
   - 因此,采用论文阈值的实现在任何可测量的规模上测到的都始终是 Dijkstra。要真正运行并测量分层构造,需要一个能缩小阈值的参数(例如下文 kingfree 实现中的 `--sep k`)。这超出了论文定理的范围,所以基准测试必须把 "paper thresholds"(论文阈值)和 "shrunk-threshold variant"(缩小阈值的变体)分开报告。
2. **目标论文本身已有一个公开的 OSS 实现**(kingfree/sssp-almost-linear,创建于 2026-11-26)。不过它的 Evaluate 阶段不完整,默认不使用(§1.4)。
3. 所要求的三个仓库(mizchi、lcs147、mjkozicki)实现的都是**有向**的 DMM+25(arXiv:2504.17033),其问题设定(有向对无向)与目标论文不同。

---

## 1. 现有 OSS

### 1.1 概览

| Repository | Exists | License(对照实际文本确认) | Latest commit SHA | Date(committer, UTC) | Main language | Target paper | Target problem |
|---|---|---|---|---|---|---|---|
| [mizchi/check-arxiv-2504-17033](https://github.com/mizchi/check-arxiv-2504-17033) | ○ | **没有 LICENSE 文件**。仅在 README 末尾有 "## License MIT";GitHub API 的 `license` 也是 `null`。`Cargo.toml` 没有 `license` 字段 | `b3f173ca43d4052c0872d824be0bd61456ab7dd3` | 2026-09-23T17:08:40Z | Rust | DMM+25 (arXiv:2504.17033), C-HD (vals.ai) | 有向、非负权重(实验以 f64 保存 [1, 2^20] 内的整数),单源到所有顶点 |
| [lcs147/bmssp](https://github.com/lcs147/bmssp) | ○ | MIT(`LICENSE`:"Copyright (c) 2025 Lucas Castro and Thailsson Clementino") | `0a8eb82805bad7f64dfa1b51e17656c9af4c4681` | 2026-01-24T15:47:14Z | C++20 header-only(按 GitHub 统计,Jupyter Notebook 占比最大) | DMM+25 (STOC 2025)。实验论文:Castro, Clementino, de Freitas, arXiv:2511.03007 | 有向、非负权重(可通过模板选择整数或浮点),单源到所有顶点 |
| [mjkozicki/SSSP](https://github.com/mjkozicki/SSSP) | ○ | MIT(`LICENSE`:"Copyright (c) 2026 Mike Kozicki") | `30231195ef664b7ba63bd9f1190c33ecd2735575` | 2026-10-08T00:35:22Z | Python 等,共 8 种语言(C#, Rust, C++, Go, Java, PHP, Python, TypeScript) | DMM+25 (arXiv:2504.17033v2) | 有向、非负实数(double),单源到所有顶点 |
| [kingfree/sssp-almost-linear](https://github.com/kingfree/sssp-almost-linear)(额外发现) | ○ | MIT(`LICENSE`:"Copyright (c) 2026 kingfree")。README 写道 "C 代码为 MIT;中文翻译的权利归原作者所有" | `151ed68a9771b4211edae7700a07bbe22c5b0c38` | 2026-11-26T10:26:16Z | Typst(译文), C11(实现) | **目标论文 arXiv:2611.34911v1** | 无向、非负实数(IEEE double 被视为 "opaque register"(不透明寄存器)),单源到所有顶点 |

验证方法:`gh api repos/<o>/<r>` 和 `gh api repos/<o>/<r>/commits?per_page=1`,各克隆的 `git rev-parse HEAD`(与 API 的 SHA 一致),以及 `LICENSE` 的文本。

### 1.2 mizchi/check-arxiv-2504-17033

- 来源:README <https://github.com/mizchi/check-arxiv-2504-17033/blob/b3f173ca43d4052c0872d824be0bd61456ab7dd3/README.md>,复现报告 `docs/chd-replication.md`(同一 SHA)
- 实现:Dijkstra(使用二叉堆和惰性删除的 `dijkstra_csr`,使用带索引 4 叉堆的 `dijkstra_dary`)、`src/dmm25.rs`(忠实于论文的 DMM+25 重新实现,含 Algorithms 1–3、Lemma 3.3 的块结构、常数度变换),以及 `src/chd.rs` 中的 C-HD。2025 年的旧实现 V1、V2 和 core 也仍然保留。
- 关于旧实现,该仓库本身在 2026-09 承认了错误:V1 与 Dijkstra 是同一份代码,V2 的 Bellman–Ford 阶段不可达,core 在顶层调用时回退到完整的 Dijkstra。旧 README 中的 "1.1–1.6x faster" 是测量噪声(只测了一次,未比较 ∞ 顶点)。
- 正确性验证:`tests/bmssp_family_test.rs`。在 3000 个 n < 40 的随机有向多重图(含零权重、并列、自环、平行边和不可达顶点)上,以及在中等规模的随机图、网格和 d 正则图上,将所有顶点与 Dijkstra 比较。在比较全部基准用例时同样报告 0 处不一致(`docs/chd-replication.md` §3)。
- 基准条件(同一文档 §4):4 vCPU 云容器,Rust 1.94 `--release`,单线程,3 次运行的中位数。每个顶点有 d 条出边的随机有向图(头部均匀随机),d ∈ {2, ⌈√lg n⌉, ⌊lg^{3/4} n⌋, 16}。网格为 4 邻域且双向。源点为 0。计时只用 `Instant` 包住每次算法调用(`src/bin/compare.rs` 第 88–96 行)。CSR 的构建在计时之外。
- 结果(相对最快 Dijkstra 的慢化倍数,README):DMM+25 慢 14–121 倍(网格上 38–50 倍)。C-HD 含预处理慢 1.47–7.3 倍。n = 2^14–2^22。
- 公平性:
  - (良好)C-HD 同时报告了含预处理和不含预处理的结果。Dijkstra 基线取三种变体中最快的。
  - (注意)Dijkstra 的 CSR 构建在计时之外,但 DMM+25 的常数度变换在计时之内(该函数接收 `Graph`)。这种不对称对算法不利。
  - 比较是在有向图上进行的,因此其数字不能直接与本项目(无向)的数字比较。

### 1.3 lcs147/bmssp

- 来源:README <https://github.com/lcs147/bmssp/blob/0a8eb82805bad7f64dfa1b51e17656c9af4c4681/README.md>,arXiv:2511.03007 <https://arxiv.org/abs/2511.03007>
- 实现:`include/bmssp.hpp`(`spp::bmssp`,最坏情形版本)和 `include/bmssp-expected.hpp`(`spp_expected::bmssp`,内存更低的期望时间版本)。`single_include/bmssp.hpp` 由 CI 生成(quom)。依赖项已打包在仓库中(`unordered_dense.h`、`floyd_rivest_select.h` 等)。
- arXiv:2511.03007 的摘要:最坏情形版本使用 O(n log^{1/3} n) 空间,期望时间版本使用 O(n) 空间,时间复杂度相同。在稀疏随机图、网格和美国道路网(最多 1000 万顶点)上,Dijkstra 在所有设定下都快 3–4 倍。他们估计最坏情形版本需要远多于 10^67 个顶点才能胜过 Dijkstra。
- 正确性验证:doctest(`tests/tests-wc.cpp`、`tests/tests-expected.cpp`)。在网格(5、25、100)和随机 3 度图(32–4096 顶点)上,将距离与打包的 Dijkstra(`helpers/dijkstra.hpp`)比较,并检查由前驱重建的路径的端点(`tests/common.cpp`)。测试了使用和不使用常数度变换两种情况。CI(`.github/workflows/ci.yaml`)运行 `./tests`。
- 基准条件:README 表格为 "random graphs, 3n edges";BMSSP 比 Dijkstra 慢 3.0–4.2 倍(n = 2^12–2^25)。README 没有说明机器和重复次数(unverified;本次未阅读 arXiv:2511.03007 的正文)。在 `experiments/experiment.cpp` 中,计时只包住 `execute(s)` 调用。`prepare_graph(false)`(图的准备)和 Dijkstra 构造函数在计时之外。
- 公平性:
  - 准备(`prepare_graph`)在计时之外。不过实验以 `false`(不做常数度变换)调用它,且图的度为 3,因此影响大概较小(unverified)。是否把常数度变换的成本计入计时,在 vendoring 时必须重新决定。
  - 距离类型是 `long long` 整数权重(`experiment.cpp` 中的 `distT`)。这与实数权重比较-加法模型的条件不同。
  - 面向有向图。若在本项目中使用,无向图需作为双向的弧输入。
  - 它打包了 `tests/graphs/USA-road-t.NY.gr`(DIMACS 道路网)。DIMACS 数据的再分发条款 unverified,因此 vendoring 时要么排除该数据,要么核对条款。

### 1.4 kingfree/sssp-almost-linear(目标论文的实现)

- 来源:README <https://github.com/kingfree/sssp-almost-linear/blob/151ed68a9771b4211edae7700a07bbe22c5b0c38/README.md>,`c/README.md`(同一 SHA)
- 创建于 2026-11-26,2 次提交,0 star(GitHub API)。
- 内容:`typst/` 中有整篇论文的中文译文(Typst 和 PDF),`c/` 中有 C11 实现(约 4000 行)。
- 实现范围(`c/README.md` 自述):
  - 在 `len.h` 中模拟比较-加法模型并统计操作次数。
  - MST(3 趟有界生长)、并行合并、LCA(用 O(log n) 的倍增代替论文的 O(1) 结构)、owner 分解、路由证书(相当于 Algorithms 2 和 3)、owner 坐标 J、外层状态图、阈值分量、层次结构组装。
- **已知偏差(自述)**:
  1. Evaluate 阶段的槽位扫描(`--slot`)在某些输入上会漏掉值或输出过大的值,因此**默认不使用**。默认路径是 "在构造出的状态图上做精确搜索"(`eval.c: state_graph_dist`)。它验证的是 §8.1 的构造能否以原始长度实现测地线;它不是论文评估过程的实现。
  2. 未实现 §7.3 的宽度日历(width calendars)。J 由辅助图上的 Dijkstra 计算。
  3. 未实现 §5.4 的 atom 替换。
  4. O(1) 的 LCA 和区间最小值被替换为 O(log n)。
- 阈值:默认遵循 §0 的阈值 N ≤ L_N^64,因此始终进入基本情形(基于扫描的 Dijkstra)。用 `--sep k` 把分隔符大小缩小为 b = L_N^k,就会运行层次结构路径。
- 正确性验证:`--verify` 将每个顶点与基于堆的 Dijkstra 比较。README 报告 "n = 20–100 的随机图和路径图 12/12 一致,n = 800 也一致"。`hier_validate` 检查层次结构的区间不变式。未发现自动化测试套件(CI 或测试运行器)。
- 基准:`--bench` 为 5 个用例(n = 500–8100)输出 "生成 + 层次结构" 和 "评估" 的时间。README 中没有数字。
- 公平性与待解决点(本次调研通过阅读代码发现;未执行):
  - 在 `c/main.c` 第 323–336 行的 `--bench` 中,当 `name[0]=='g' && name[1]=='r'`(即 "grid-…")时调用 `gen_rand`,而对 "rand-…" 调用 `gen_grid`。**用例名称与图类型似乎被对调了。**
  - 默认路径不是论文的评估过程(上文第 1 项),因此该实现的时间不是论文算法的时间。
  - README 和 c/README 把主算法称为 "算法 9"(Algorithm 9),把评估称为 "算法 7、8"(Algorithms 7 and 8)。但在 arXiv v1 的文本中,主算法是 Algorithm 8(`ShortestPaths`),评估是 Algorithm 3(`Evaluate`)。它们指的是哪个版本的编号 unverified。
  - 中文译文是对 arXiv 非独占许可下论文的翻译,因此是否允许再分发 unverified。它不会被 vendoring 进本项目。

### 1.5 mjkozicki/SSSP

- 来源:README <https://github.com/mjkozicki/SSSP/blob/30231195ef664b7ba63bd9f1190c33ecd2735575/README.md>,`benchmark/README.md`(同一 SHA)
- 实现:8 种语言中每种都有 Dijkstra 和 `duan_mao_shu_yin`。README 自己写道 "the frontier data structure is a simplified version that preserves correctness; using the Lemma 3.3 block structure would improve constants"。
  - Rust 的 `FrontierDS`(`rust/src/sssp.rs` 第 342–411 行)在每次 `insert` 时用 `Vec::retain` 线性删除,在每次 `pull` 时对全部元素排序。它不保证 Lemma 3.3 的复杂度,似乎达不到 O(m log^{2/3} n)(本次调研的阅读结果;未测量)。
- 正确性验证:每种语言都有与 Dijkstra 比较的测试,例如 Rust 的 `matches_dijkstra_small` 和 `matches_dijkstra_parallel_paths`(容差 1e-10)以及 Python 的 `test_matches_dijkstra`。不过它们集中在小图上。
- 基准条件(`benchmark/README.md`):50,000 个顶点,模仿约 20 个区域网络的有向图(区域内每顶点 2–8 条出边,权重 0.5–5;区域间每顶点 0–2 条出边,权重 20–150,seed = 42)。每种语言启动一个 Docker 容器,`wall_sec` 是 **从 `docker run` 到 `docker wait` 的时间**(`benchmark/run_benchmarks.py` 第 209–233 行)。默认只测一次。
- 结果(`benchmark/results.json`):在大多数语言中,DMSY 和 Dijkstra 都约为 1.9–2.0 秒,`peak_mem_mb` 和 `total_cpu_sec` 为 0.0。
- 公平性:
  - 测量包含容器启动和读取图文件。算法之间的差异几乎看不出来,因此不能用于比较算法。
  - 只测一次。许多语言中内存和 CPU 的采集为 0,测量似乎没有生效。
  - `rust/src/bin/compare.rs` 含有另一个取多次重复中位数的比较,但其结果不在仓库中(unverified)。

---

## 2. 先行研究(以 2023–2026 为主)

"Match" 指其设定是否与本项目的问题设定(无向、非负实数、单源到所有顶点、比较-加法模型)一致。

| Abbrev. | Reference | Problem setting | Complexity | Model | Match | Use as a comparison target? |
|---|---|---|---|---|---|---|
| Dij59 / FT87 | Dijkstra 1959, Fredman–Tarjan 1987 (Fibonacci heap)。DOI 10.1007/BF01386390, 10.1145/28869.28874 | 有向和无向,非负实数 | O(m + n log n)。使用二叉堆时为 O((m+n) log n) | 比较-加法 | ○ | **必需(基线)**。自行实现(二叉堆、4 叉堆,可行的话还有配对堆) |
| Tho99 | Thorup, "Undirected Single-Source Shortest Paths with Positive Integer Weights in Linear Time", J. ACM 46(3), 1999。DOI 10.1145/316542.316548 | 无向,**正整数**权重 | O(m)(最坏情形) | word RAM(使用权重的位表示) | △(权重和模型不同) | 可选。仅在整数权重子集上作为参考测量。不纳入实数权重的比较。公开实现 unverified |
| PR05 | Pettie, Ramachandran, "A Shortest Path Algorithm for Real-Weighted Undirected Graphs", SIAM J. Comput. 34(6), 2005。DOI 10.1137/S0097539702419650 | 无向,非负实数 | 预处理 O(MST(m,n) + min{n log n, n log log r}),其中 r 为权重之比 | 比较-加法 | ○ | 理论比较对象。公开实现 unverified。自行实现是一项庞大的工作 |
| DMS+23 | Duan, Mao, Shu, Yin, "A Randomized Algorithm for Single-Source Shortest Path on Undirected Real-Weighted Graphs", FOCS 2023。arXiv:2307.04139 <https://arxiv.org/abs/2307.04139> | 无向,非负实数 | 随机化,O(m √(log n · log log n)) | 比较-加法 | ○ | **有力候选**。问题设定完全一致。摘要称其 "simple in structure and easy to implement"。本次搜索(`gh search repos`)未发现公开实现 |
| Yan26 | Shuyi Yan, "Lossless Derandomization for Undirected Single-Source Shortest Paths and Approximate Distance Oracles", SOSA 2026。arXiv:2510.12598 <https://arxiv.org/abs/2510.12598> | 无向,非负实数 | 确定性,O(m √(log n log log n)) | 比较-加法 | ○ | DMS+23 的确定性版本。若实现了 DMS+23,有余力时再加入 |
| DMM+25 | Duan, Mao, Mao, Shu, Yin, "Breaking the Sorting Barrier for Directed Single-Source Shortest Paths", STOC 2025。arXiv:2504.17033 <https://arxiv.org/abs/2504.17033> | **有向**,非负实数 | 确定性,O(m log^{2/3} n) | 比较-加法 | △(有向) | 次要比较对象。无向图可作为双向的弧来求解,但并未利用无向性。可使用 lcs147/bmssp(§3) |
| DMS+26 | Duan, Mao, Shu, Yin, "A Faster Directed Single-Source Shortest Path Algorithm", ICALP 2026 (LIPIcs 374, 81:1–81:23)。arXiv:2602.07868 <https://arxiv.org/abs/2602.07868> | 有向,非负实数 | 确定性,O(m √(log n) + √(mn log n log log n))。稀疏图上为 O(m √(log n log log n)) | 比较-加法 | △(有向) | 参考。公开实现 unverified |
| GZ26 | Grishin, Zabokritskiy, "Improved Bounds for Directed SSSP in Sparse and Intermediate-Density Graphs", arXiv:2611.26453 <https://arxiv.org/abs/2611.26453>(2026-11-19 提交) | 有向,非负实数 | 确定性,m = O(n) 时为 O(n √(log n) (log log n)^{1/4}) | 比较-加法 | △(有向) | 参考。预印本 |
| KR26 | Kadria, Roditty, "A Faster Undirected Single-Source Shortest Path Algorithm", arXiv:2609.15247 <https://arxiv.org/abs/2609.15247>(2026-09-14 提交) | 无向,非负实数 | O(m log^{1/2} n · (log log n)^{1/4} · (log log log n)^{1/4})。目标论文称其为 "with high probability" | 比较-加法 | ○ | 直接的理论竞争对手(目标论文之前最好的无向结果)。公开实现 unverified。预印本 |
| HHR+24 | Haeupler, Hladík, Rozhoň, Tarjan, Tětek, "Universal Optimality of Dijkstra via Beyond-Worst-Case Heaps", FOCS 2024(pp. 2099–2130,普林斯顿书目记录 <https://collaborate.princeton.edu/en/publications/universal-optimality-of-dijkstra-via-beyond-worst-case-heaps/>)。arXiv:2311.11793 <https://arxiv.org/abs/2311.11793> | 有向和无向,非负实数。主要是 "distance ordering problem"(按距离给顶点排序) | 使用具有 working-set 性质的堆的 Dijkstra 是普适最优的 | 比较-加法 | △(输出的是顺序而不仅是距离值的问题) | 对理解很重要。与目标论文 Section 1 开头紧密相关(计算距离不同于对距离排序)。working-set 堆的公开实现 unverified。FOCS 2024 Best Paper Award 只能通过二手来源(arXiv:2504.17327)确认 |
| Castro+25 | Castro, Clementino, de Freitas, "Implementation and Brief Experimental Analysis of the Duan et al. (2025) Algorithm", arXiv:2511.03007 <https://arxiv.org/abs/2511.03007> | DMM+25 的实现与实验 | — | — | — | 基准方法论的参考(道路网、网格、随机图) |

注:

- 目标论文 Section 1.1(Related Work)所引用的文献,是上表中的 Dij59、FT87、Tho99、PR05、DMS+23、Yan26、DMM+25、DMS+26、GZ26 和 KR26,外加 FRT04(Fakcharoenphol、Rao 和 Talwar 的树嵌入;DOI 10.1016/J.JCSS.2004.04.011)。FRT04 是一个组件(路由估计器),因此不是比较对象。HHR+24 不在目标论文的参考文献中(应要求添加)。
- 目标论文的正文指出,其界是在 "connected undirected inputs or on the reachable subgraph of a directed input"(连通的无向输入,或有向输入的可达子图)上陈述的。识别可达顶点并写出全部 n 个输出标签至多增加 O(m + n)(Section 1.1)。比较时,各方法对不可达顶点(写入 ∞)的处理必须保持一致。
- 目标论文的时间复杂度是**期望**界(随机化)。DMS+23 也是随机化的,而 Yan26、DMM+25、DMS+26 和 GZ26 是确定性的。随机化方法必须用多个种子运行,并以分布形式报告。

---

## 3. vendoring 现有 OSS 与重新实现

### 3.1 建议

**目标论文和基线 Dijkstra 都自行实现。不要 vendoring 现有 OSS 的代码。** 唯一的例外:只有在需要与有向 BMSSP 做次要比较时,才把 lcs147/bmssp 固定版本作为外部依赖。

| Repository | 建议 | 理由 |
|---|---|---|
| kingfree/sssp-almost-linear | **不要 vendoring**。仅作为只读的交叉参考 | 它是 MIT,因此 vendoring 在法律上可行。但是 Evaluate 阶段不完整且默认不使用,基准用例名称与图似乎被对调,没有自动化测试,且创建一天后只有 2 次提交。它不能作为正确复现论文的参考。译文的再分发许可 unverified。本项目的实现完成后,也许可以用它在相同输入上比较距离,以 "交叉检查独立实现"(是否运行另行决定) |
| lcs147/bmssp | **仅在需要次要比较时**,按 SHA 固定为外部依赖 | MIT,header-only C++20,带有与 Dijkstra 比较的测试和 CI。本项目是 C++/CMake,因此易于集成。但其实验假定有向图和整数权重,因此无向图必须以双向的弧给出,且在使用前必须用我们自己的测试检查把权重类型设为 `double` 时的行为 |
| mizchi/check-arxiv-2504-17033 | **不要 vendoring**。仅把其数字和方法论作为文献引用 | 没有 LICENSE 文件(README 里只说了 "MIT")。MIT 要求包含版权声明和许可声明,但该文本缺失。它是 Rust 写的,也与本项目的语言不同。它统计比较次数的方式(`count-cmp`)和基准方法论是有用的参考 |
| mjkozicki/SSSP | **不要 vendoring** | frontier 是简化版,不满足 DMM+25 的复杂度。基准只测一次且包含容器启动,不是算法比较 |

### 3.2 vendoring 时如何固定提交(针对 lcs147/bmssp)

- **写 40 个字符的 SHA,而不是标签或分支**,如 CMake 的 `FetchContent_Declare(bmssp GIT_REPOSITORY https://github.com/lcs147/bmssp.git GIT_TAG 0a8eb82805bad7f64dfa1b51e17656c9af4c4681)`。使用 git submodule 也可以(SHA 记录在 `.gitmodules` 和树中)。
- 只使用 `include/`(包括 `external/`);不要使用 `tests/graphs/` 中的数据(包括源自 DIMACS 的数据)或 `experiments/`。
- 按 MIT 的要求,在发布物中附上其 `LICENSE`,并在 `THIRD_PARTY_NOTICES` 或类似文件中记录名称、SHA 和许可。所打包的 `include/external/*`(`unordered_dense.h`、`floyd_rivest_select.h`、`median_of_ninthers.h` 等)的许可,在 vendoring 之前必须在各文件头部核对(本次调研中 unverified)。
- 本项目自身的仓库中还没有许可(未发现 `LICENSE`)。与 MIT 代码的兼容性,将在确定本项目的许可后最终确认(MIT 通常与几乎所有其他许可兼容)。

### 3.3 比较对象的优先级(提案)

1. Dijkstra(自行实现;二叉堆和 4 叉堆;同时统计比较次数)。如 §0 所述,采用目标论文阈值的实现会进入 Dijkstra 基本情形,所以这实际上成为主要比较。
2. 目标论文的实现(自行实现)。分开报告 "paper thresholds" 和 "shrunk-threshold variant"。
3. DMS+23(若自行实现)。唯一一个既与问题设定完全一致、又易于实现的先行非 Dijkstra 方法。
4. DMM+25(使用固定版本的 lcs147/bmssp)。须明确说明这些是把有向方法应用于无向输入得到的参考值。

---

## 4. 未验证的事项

- arXiv:2511.03007 正文中的测量环境和重复次数(只读了摘要)。
- lcs147/bmssp 的 README 表格是在哪台机器、什么条件下测得的。
- lcs147/bmssp 所打包的 `include/external/*` 的许可。
- DIMACS 道路网数据的再分发条款。
- Tho99、PR05、DMS+23、Yan26、DMS+26、KR26、GZ26 和 HHR+24 是否存在公开实现(用 `gh search repos` 按论文 ID 和标题搜索未发现,但搜索并不详尽)。
- kingfree 实现的 `--bench` 中被对调的用例名称是否出于有意(代码只读未执行)。
- kingfree 实现的算法编号对应论文的哪个版本。
- HHR+24 的 FOCS 2024 Best Paper Award(只能通过二手来源确认)。
