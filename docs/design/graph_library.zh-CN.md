# 图库

[English](graph_library.md) | [日本語](graph_library.ja.md) | **简体中文**

图的表示、文本 IO 和可复现的生成器,由所有 SSSP 实现、基准测试和测试共用。图是**无向图**,权重**非负**;权重类型为 `std::uint64_t`(精确)和 `double`(IEEE-754)。

## 1. 结构

| 组件 | 文件 | 命名空间 | 作用 |
|---|---|---|---|
| CSR 图 | `include/sssp/graph.hpp` | `sssp` | `CsrGraph<W>`、`Edge`、`Arc`、`validate_edges`、`WeightTraits`。被所有算法共用,本库**不修改**它。 |
| 格式 | `include/sssp/graph_formats.hpp`(仅头文件) | `sssp` | `EdgeList<W>`、`AdjacencyList<W>`、转换函数、`timed()` |
| 文本 IO | `include/sssp/graph_io.hpp`、`src/graph/io.cpp` | `sssp::io` | 项目的边列表格式和 DIMACS `.gr`、`ParseError` |
| 生成器 | `include/sssp/generators.hpp`、`src/graph/generators.cpp` | `sssp::gen` | PRNG、确定性数学函数、权重分布、图生成器、目录、`generate()` |
| CLI | `tools/sssp_gen.cpp` | — | `sssp_gen`:生成图并写入文件 |

`io.cpp` 和 `generators.cpp` 中的模板针对 `std::uint64_t` 和 `double` 做了显式实例化。

## 2. 表示与转换

- `EdgeList<W> { std::size_t num_vertices; std::vector<Edge<W>> edges; }`:IO 和生成器之间交换的格式。
- `AdjacencyList<W>::from_edges(n, edges)`:vector 的 vector。它使用与 `CsrGraph<W>` 相同的
  `Arc<W>` 布局和相同的访问接口(`num_vertices`、`num_edges`、`num_arcs`、
  `degree`、`arcs(v)`、`edge(id)`、`edges()`)。它对输入的校验方式与
  `CsrGraph::from_edges` 相同。
- 转换是不带缓存的自由函数:

  | 函数 | 输入 |
  |---|---|
  | `to_csr` | `const EdgeList&`、`EdgeList&&`(移动边的 vector)或 `const AdjacencyList&` |
  | `to_edge_list` | `CsrGraph` 或 `AdjacencyList` |
  | `to_adjacency_list` | `EdgeList` 或 `CsrGraph` |

**顺序约定。** 所有表示都使用相同的边 id:`EdgeId i` 在任何地方都指同一条边。`AdjacencyList::arcs(v)` 列出弧的顺序与 `CsrGraph::arcs(v)` 完全相同:对每条边按 id 顺序,先是 u 侧的弧,再是 v 侧的弧。自环 `{u,u}` 是一条弧。因此,无论采用哪种布局,算法访问弧的顺序都相同。

**单独计时转换。** 基准测试不得把转换的开销计入算法时间:

```cpp
auto [g, convert_seconds] = sssp::timed([&] { return sssp::to_csr(std::move(edge_list)); });
// run and time the algorithm on g; report convert_seconds as a separate column
```

`timed(f)` 返回 `Timed<T>{value, seconds}`,并使用 `std::chrono::steady_clock`。

## 3. 文件格式(`sssp::io`)

### 项目边列表格式(`*.txt`)

```
n m
u v w      # exactly m lines, 0-based vertices, undirected edge {u,v}
```

- 文件中任何位置的空行,以及第一个非空白字符为 `#` 的行,都会被忽略。
- 文件必须恰好包含 `m` 行边。过少或过多都是错误。
- 自环和平行边按原样保留。

### DIMACS 第 9 届挑战赛 `.gr`

```
c comment
p sp n m   # exactly once, before the first arc; m = number of ARC lines
a u v w    # 1-based, directed arc u->v
```

图分两步按无向图读取。

1. `DimacsArcMode` 决定弧如何变成边。
   - `kArcAsEdge`(默认):每条弧变成一条边。因此一对对称的弧会变成两条平行边。距离不变,但 `m` 会加倍。
   - `kMergeReversePairs`:权重逐位相同的弧 `u->v` 和 `v->u` 按文件顺序贪心配对。配成对的两条弧变成**一条**边,放在该对中第一条弧的位置,并采用第一条弧的方向。找不到配对的弧(单向弧,或权重不对称)各自变成一条边。自环弧始终是一条边。DIMACS 道路网络请使用此模式。
2. 然后由 `DuplicateEdgePolicy` 处理剩余的平行边。两条边的无序端点对相同时即为平行;同一顶点上的两个自环也算平行。

   | 策略 | 效果 |
   |---|---|
   | `kKeepAll`(默认) | 保留所有边(多重图) |
   | `kKeepMinWeight` | 保留一条边,位置和方向取其第一次出现时的,权重取最小值 |
   | `kReject` | 在第一条重复边所在的行抛出 `ParseError` |

`write_dimacs` 把每条边写成两条弧(先 `u->v`,再 `v->u`),自环写成一条弧。
用 `kMergeReversePairs` 读回该文件,可以精确还原原来的边列表:顺序相同、方向相同、权重逐位相同。

### API

- 流:`read_edge_list<W>(istream&)`、`write_edge_list(ostream&, el)`、
  `read_dimacs<W>(istream&, DimacsReadOptions{arc_mode, duplicates})`、`write_dimacs(ostream&, el)`。
- 文件:`read_edge_list_file`、`write_edge_list_file`、`read_dimacs_file`、`write_dimacs_file`。
  如果文件无法打开或写入,它们会抛出 `std::runtime_error`。

### 数值

- 解析和输出使用 `std::from_chars` / `std::to_chars`,因此结果不依赖于 locale。
- `uint64_t` 权重必须是 0..2^64-1 范围内的普通无符号十进制数。带符号、小数或指数的写法会被拒绝。
- `double` 权重可以是任何有限的非负十进制或指数形式。它们以 `chars_format::general` 和精度 17 写出(与 `%.17g` 相同),这样每个有限的 double 都能逐位精确地往返还原。`-0` 会被规范化为 `+0`。

### 错误

所有格式错误的输入都会抛出 `io::ParseError`(派生自 `std::runtime_error`)。消息的形式为
`line N: reason`,`line()` 返回从 1 开始的行号(0 表示输入结束或整个输入)。会检测到的错误有:

- 负数、NaN、无穷大或溢出的权重
- 越界的顶点,包括 DIMACS 中的顶点 0
- 缺失、重复或格式错误的头部
- 边行或弧行过少或过多
- 非数字的记号,或一行中记号的个数不对
- `n >= kNoVertex` 或 `m >= kNoEdge`

读取器把预留空间的上限限制在 2^24 个元素,因此头部中很大的 `m` 不会让它预先分配巨大的缓冲区。

## 4. 可复现性约定

每个生成器都是 (参数, 种子) 的纯函数,在所有平台和标准库上产生**逐位相同**的输出。下面的规则保证了这一点。

- **自带 PRNG。** 播种以及派生子流种子使用 `SplitMix64`,随机流使用
  `Xoshiro256ss`(xoshiro256** 1.0)。我们不使用 `<random>`:C++ 标准把 `std::*_distribution`
  各个类的算法留给实现定义,因此同样的引擎输出在 libstdc++、libc++ 和 MSVC 上会得到不同的值。这两个生成器都已对照参考实现的已知答案向量进行了检查。
- **精确的整数。** `uniform_below(b)` 拒绝小于 `2^64 mod b` 的原始值,并返回
  `r % b`。这样没有取模偏差,也不使用浮点数。`uniform_range(lo, hi)` 是闭区间,并能处理完整的 64 位范围。
- **精确的实数。** `uniform01()` 计算 `(next() >> 11) * 2^-53`,得到 [0,1) 中的值,带有 53 位随机性且没有舍入。
- **确定性数学函数。** `det_log` 和 `det_exp` 只使用正确舍入的 IEEE 运算
  (`+ - * /`、`frexp`、`ldexp`、`floor`)。libm 的 `log`、`exp` 和 `pow` 并非正确舍入,且不同平台之间存在差异,所以我们不使用它们。`det_log` 在 [√½, √2) 中的尾数上使用 atanh 级数。`det_exp` 用 ln 2 的 Cody–Waite 拆分
  (fdlibm 的 `ln2_hi` / `ln2_lo`)对参数做归约,然后对泰勒级数求和。两者与 libm 的相对误差都在约 1e-14 以内。`sqrt` 由 IEEE 754 保证正确舍入,因此可以放心使用。
- **禁止 FMA 收缩。** 融合乘加只舍入一次而不是两次,因此对表达式做收缩会改变结果。`generators.cpp` 包含
  `#pragma clang fp contract(off)`。GCC 不需要任何处理,因为在 ISO 模式(`-std=c++20`、
  `CMAKE_CXX_EXTENSIONS OFF`)下它本来就使用 `-ffp-contract=off`。
- **禁止无序遍历。** `unordered_set` 和 `unordered_map` 只用于成员判断和计数。输出顺序从不依赖它们的遍历顺序。
- **独立的子流。** `substream_seed(seed, stream)` 为每种用途派生一条流:拓扑 = 0,权重 = 1,连接边 = 2。因此在种子相同的情况下,无论权重分布如何,也无论 `W` 是 `uint64_t` 还是 `double`,拓扑都是相同的。

### 版本之间的输出变化

- **P2-GRAPH-GEOMETRIC(新类型 `geometric`)。** 新增该生成器没有改变任何已有类型的输出。我们在 2,526 组参数上比较了修改前后的 `sssp_gen`:九种已有类型 × 两种权重类型 × 全部七种权重分布 × 种子 {1, 2, 3, 17, 1000} × n ∈ {50, 1000} × 有无
  `--connected`(regular 在 keep/drop/switch 之间轮换),另加带 `--int-scale` 的 road 以及 DIMACS 输出。所有哈希都相同,`--list` 只是在末尾多了新的条目。

- **P2-DIFF-REGULAR-SWITCH(`random_regular`、`kSwitchToSimple`;`sssp_gen --type regular
  --multi switch`)。** 修复之前,一些合法的参数会抛出 "switching did not converge",
  最常见于小而稠密的情形(d = n − 1,n ≤ 12 时约 3–10% 的种子;例如
  n = 3, d = 2, seed = 20)。这些参数现在会返回一个简单的 d 正则图。
  生成器会从同一条拓扑流中抽取一个新的配置,重新开始(见
  §9)。以前成功的每组参数都给出逐字节相同的输出,因为第一个配置及其交换使用该流的方式与之前完全一致。我们通过在 74,955 组参数上比较新旧库来检查这一点:三种策略、n ≤ 16 × 所有合法的 d ×
  种子 1..200,以及 n ∈ {100, 1000, 5000} × d ∈ {2, 3, 4, 7, 10, 20},有无
  `connected`。以前成功的 74,778 组给出相同的哈希。其余 177 组
  以前会抛出异常,现在会成功。`kKeep` 和 `kDrop` 没有变化。

## 5. 权重分布(`WeightSpec`)

`WeightSpec` 有字段 `dist`、`lo`(默认 1)、`hi`(100)、`mean`(10)、
`zero_probability`(0.25)和 `wide_log2`(50)。如果它们相互矛盾,`validate_weight_spec` 会抛出
`std::invalid_argument`。

所有分布都要求 `lo` 和 `hi` 为有限值、`lo >= 0` 且 `hi <= 2^53`。对于整数
`W`,实数抽样值会四舍五入到最近的整数,恰好居中时向上取整。

| `WeightDist` | CLI 名称 | 抽样 | `uint64_t` | `double` | 额外约束 |
|---|---|---|---|---|---|
| `kUniformInt` | `uniform_int` | 在 [⌈lo⌉, ⌊hi⌋] 的整数上均匀分布 | 精确 | 整数值的 double | 区间内含有整数 |
| `kUniformReal` | `uniform_real` | `lo + (hi-lo)·U`,U ∈ [0,1) | 舍入 | 实数 | lo ≤ hi |
| `kLogUniform` | `log_uniform` | `exp(uniform(log lo, log hi))`,夹在 [lo,hi] 之内 | 舍入 | 实数 | 0 < lo ≤ hi |
| `kExponential` | `exponential` | `-mean·log(1-U)` | 舍入 | 实数 | 0 < mean ≤ 2^47 |
| `kConstant` | `constant` | `lo`(全零权重时用 `lo = 0`) | 舍入 | `lo` | — |
| `kZeroMixed` | `zero_mixed` | 以概率 `zero_probability` 取 0,否则在 [lo,hi] 上均匀分布 | 整数均匀 | 在 [lo,hi) 上的实数均匀 | zero_probability ∈ [0,1];区间内含有整数 |
| `kWideRange` | `wide_range` | 二进制指数 e 在 [0, K] 上均匀分布,然后在 [2^e, 2^(e+1)) 内均匀分布;e = K 时恰好为 2^K;覆盖 [1, 2^K] | 精确 | `2^e·(1+U)` | K = `wide_log2` ∈ [0, 62] |

`to_string(WeightDist)` 和 `parse_weight_dist(name)` 在枚举与 CLI 名称之间转换;遇到未知名称时 `parse_weight_dist` 会抛出异常。`sample_weight<W>(spec, rng)` 做一次抽样。

## 6. 生成器

每个生成器都返回 `EdgeList<W>`,并接受 `GenOptions{weights, seed = 1, connected = false}`。

| CLI `--type` | 函数 | 参数 | 边数 | 简单图 | 连通 |
|---|---|---|---|---|---|
| `er` | `erdos_renyi(n, avg_degree, opt)` / `erdos_renyi_gnm(n, m, opt)` | m = round(n·avg_degree/2) ≤ n(n−1)/2 | m | 是(对点对下标做 Floyd 抽样) | 仅在指定 `connected` 时 |
| `regular` | `random_regular(n, d, policy, opt)` | n·d 为偶数;`kSwitchToSimple` 需要 d < n | n·d/2(`kDrop` 时更少) | `kKeep`:否;`kDrop`:是(度 ≤ d);`kSwitchToSimple`:是(度恰为 d,通过保持度数的双边交换) | 仅在指定 `connected` 时 |
| `ba` | `barabasi_albert(n, m, opt)` | 1 ≤ m < n | m(m+1)/2 + (n−m−1)·m | 是 | 是 |
| `grid` | `grid_2d(rows, cols, opt)`;顶点 id = r·cols + c | rows、cols(默认各为 ⌊√n⌋) | rows(cols−1) + cols(rows−1) | 是 | 是 |
| `road` | `road_like(rows, cols, RoadLikeOptions, opt)` | jitter(默认 0.3,< 0.5)、keep_probability(0.9)、shortcut_probability(0.15)、shortcut_radius(2,≥ 2)、integer_scale(1000) | ≈ grid × keep_probability + 捷径 + 恢复的街道 | 是 | 是 |
| `tree` | `random_tree(n, opt)` | n | n−1 | 是 | 是(由随机 Prüfer 序列得到,在带标号树上均匀分布) |
| `chung_lu` | `chung_lu(n, avg_degree, exponent, opt)` | exponent > 2;m ≤ n(n−1)/4 | m = round(n·avg_degree/2) | 是(自环和重复边会重新抽取) | 仅在指定 `connected` 时 |
| `adversarial_decrease_key` | `adversarial_decrease_key(k)`,CLI 中 k = n−1 | k | k + k(k−1)/2 | 是 | 是 |
| `adversarial_long_path` | `adversarial_long_path(n)` | n | (n−1) + (n−2) | 是 | 是 |
| `geometric` | `random_geometric(n, r, GeometricOptions, opt)` | r ≥ 0 且有限(默认 `geometric_radius(n, avg_degree)` = √(avg_degree/(πn)));权重为 `kFromSpec` 或 `kEuclidean`;integer_scale(1000) | ≈ n·(n−1)·F(r)/2,F 见下文 | 是 | 仅在指定 `connected` 时 |

对不太直观的生成器的更多说明:

- **类道路图(Road-like)。** 网格点的每个坐标轴上都加上 U[−jitter, jitter] 的扰动。每条街道(网格边)以概率 `keep_probability` 保留。以概率
  `shortcut_probability`,一个顶点还会得到一条通向切比雪夫半径
  R 以内、且不是网格邻居的顶点的捷径。然后按 Kruskal 顺序恢复被删除的街道,使图连通,因此每条边都保持局部。权重是扰动后各点之间的欧几里得长度;整数权重为 `max(1, round(len·integer_scale))`。road
  生成器忽略 `opt.weights`。
- **Chung–Lu。** 顶点 i 的权重 w_i ∝ (i+1)^(−1/(exponent−1)),用
  `det_exp`/`det_log` 计算。一条边的两个端点都通过前缀和二分查找,按与 w 成正比的概率抽取。
- **Barabási–Albert。** 从 m+1 个顶点上的团开始。每个新顶点连向 m 个不同的目标,这些目标从边端点的列表中按与度成正比的概率抽取。
- **对抗性生成器。** 它们的权重由构造固定,不使用 `opt.weights`。
- **随机几何图**(`geometric`)。见下一小节。

### 随机几何图

`random_geometric<W>(n, r, GeometricOptions{weights, integer_scale}, opt)` 在单位正方形中构造 G(n, r)。

- **点。** `random_geometric_points(n, seed)` 返回这些点:顶点 i 位于第 i 个点,点取自拓扑子流(0),先 `x = uniform01()`,再 `y = uniform01()`。因此无论半径、无论权重如何选择,点都相同。
- **边。** 对 u < v,`{u, v}` 是边当且仅当在 IEEE double 中 `dx*dx + dy*dy <= r*r`
  (dx = x_u − x_v,dy = y_u − y_v,不使用 FMA)。这个表达式就是约定:测试可以对所有点对重新计算它。边按 (u, v) 的字典序列出,因此输出不依赖于如何找到邻居。r = 0 时没有边(两个点必须重合)。r ≥ √2(正方形的对角线长度)时得到完全图。
- **邻居搜索。** 点被放入 k × k 个单元,其中 k·r ≤ 1 − 10⁻⁹ 且 k² ≤ n。距离不超过 r 的每一对点于是都位于同一个单元或相邻的单元中(10⁻⁹ 的余量吸收了判定式和单元坐标的舍入误差;证明见
  `geometric_cells_per_side`)。计数排序按单元对点分组,每个顶点 v 扫描其 3 × 3 块中 u < v 的点,再按 u 做稳定的计数排序,使点对有序。期望时间为 O(n + m),内存为 O(n + m)。
- **由度数得到半径。** `geometric_radius(n, d)` = √(d/(πn))(n < 2 时为 0)。它忽略边界。期望平均度为 (n−1)·F(r),其中
  F(r) = πr² − (8/3)r³ + r⁴/2(r ≤ 1)是正方形中两个均匀点的距离不超过 r 的概率。实际的平均度比 d 低,相对差约为 8r/(3π)
  (n = 2·10⁴、d = 8 时约为 1%)。
- **权重。** `GeometricWeights::kFromSpec`(CLI `--geo-weights random`,默认)与其他所有生成器一样,按边的顺序从权重子流中根据 `opt.weights` 抽样。
  `kEuclidean`(CLI `--geo-weights euclidean`)使用边长 √(dx² + dy²)。`double`
  原样保留,整数 W 使用 `round(length · integer_scale)`(取最近,恰好居中时向上)。对于非常接近的点,这可能为 0;零权重是合法的。integer_scale 必须在 (0, 2⁵²] 内。
  使用 `kEuclidean` 时,`opt.weights` 被忽略,也不做校验。
- **`connected`。** `connect_components` 与其他生成器一样,用连接子流(2)连接各分量。使用 `kEuclidean` 时,连接用的边也按其长度赋权,因此每个权重都是欧几里得距离。

### 对抗性构造

- **Decrease-key**(k+1 个顶点,源点 0,v_i = i)。边为权重 1 的 `{0,1}`、
  i > 1 时权重 2k 的 `{0,i}`,以及 1 ≤ i < j 时权重 2(j−i)−1 的 `{i,j}`。
  - 权重为 1 的边构成的链给出 dist(v_j) ≤ j,而每条边 `{i,j}` 的代价至少为
    |j−i|,所以 dist(v_j) = j。
  - 因此 Dijkstra 按 1, 2, …, k 的顺序确定顶点。
  - 确定 v_i 会降低**每个** j > i 的 v_j 的键,从 2j−i 降到 2j−i−1。
  - 因此基于堆的 Dijkstra 恰好执行 k(k−1)/2 = Θ(m) 次 decrease-key。这是二叉堆的最坏情形,也是斐波那契堆当初要解决的情形。
- **长路径**(n 个顶点)。路径边 `{i,i+1}` 的权重为 1,诱饵跳边
  `{i,i+2}` 的权重为 3(比它所绕过的两条路径边更大)。从 0 到 n−1 的唯一最短路径有 n−1 跳,每条诱饵边都会被松弛然后被拒绝。这会考验代价随跳数或递归深度增长的方法。

### `connected` 选项与 `connect_components`

`connect_components(el, weights, seed)` 的工作方式如下。

1. 各分量按其最小顶点排序。
2. 每个分量 i > 0 通过一条边连到一个均匀随机选取的更早的分量,这条边连接两者的均匀随机成员。这样就在各分量之上形成一棵随机树。
3. 恰好追加(分量数 − 1)条边。原来的边仍是前缀,保持不变。

生成器把连接子流(2)作为种子传入。按构造就连通的生成器会跳过这一步。`count_components(el)` 统计分量数,每个孤立顶点算一个。

### 目录与通用入口

- `generator_catalog()` 为每种生成器类型返回 `{name, description, parameters}`。
- `generate<W>(GraphSpec)` 按 `GraphSpec::type` 分派。`GraphSpec` 包含 `n`、
  `avg_degree`、`degree`、`attach`、`rows`、`cols`、`exponent`、`multi`、`road`、`radius`
  (`std::optional<double>`;未设置表示 `geometric_radius(n, avg_degree)`)、`geometric` 和
  `options`。未知的类型会抛出 `std::invalid_argument`。
- `parse_multi_edge_policy("keep"|"drop"|"switch")` 转换正则图策略的 CLI 名称。
- `parse_geometric_weights("random"|"euclidean")` 转换几何图权重选择的 CLI 名称。

## 7. 小型测试图

`template <class W> std::vector<NamedGraph<W>> small_graph_suite();` 返回边界情形的图,其中

```cpp
template <class W> struct NamedGraph { std::string name; EdgeList<W> graph; Vertex source = 0; std::string note; };
```

指定大小的构造器位于 `sssp::gen::small`:`single_vertex()`、`no_edges(n)`、
`two_vertices(w)`、`path(n, w)`、`star(n, w)`、`cycle(n, w)`、`complete(n, w)`(边 {i,j}
的权重为 w·(1 + (i + 2j) mod 5))、`binary_tree(n, w)` 和 `many_equal_paths(k)`;
`long_shortest_path` 是 `adversarial_long_path(40)`。其余情形是手工构造的固定
边列表,只存在于该套件内部。该套件仅有头文件,因此适用于任何权重类型。这些图是:

- **结构:** `single_vertex`、`no_edges`、`two_vertices`、`path`、`star`、`cycle`、
  `complete`、`binary_tree`、`disconnected`、`multi_edges`、`self_loops`。
- **权重:** `zero_weights`、`all_equal_weights`、`mixed_magnitudes`、`extreme_ratio`。
- **最短路径形态:** `long_shortest_path`,以及 `many_equal_paths`,即 k 个菱形串成的链,
  有 2^k 条相等的最短路径。

差分测试和性质测试应遍历整个套件,而不是自己列举各个情形。

## 8. CLI(`sssp_gen`)

```sh
sssp_gen --list                                   # generator types and weight distributions
sssp_gen --type er --n 1000 --avg-degree 4 --weights uniform_int --seed 1 --out g.txt
sssp_gen --type road --rows 100 --cols 100 --weight-type double --format dimacs --out road.gr
sssp_gen --type regular --n 1000 --degree 3 --multi switch --weights wide_range --wide-log2 50 --out r.txt
sssp_gen --type er --n 10000 --avg-degree 1 --connected --weights zero_mixed --zero-prob 0.3 --out z.txt
sssp_gen --type geometric --n 100000 --avg-degree 8 --geo-weights euclidean --weight-type double --out geo.txt
sssp_gen --type geometric --n 5000 --radius 0.03 --weights log_uniform --lo 1 --hi 1000 --connected --out geo2.txt
```

| 组 | 选项 |
|---|---|
| 输出 | `--format edgelist`(默认)或 `dimacs`;`--weight-type uint64`(默认)或 `double`;`--out` |
| 权重 | `--weights <dist>`、`--lo`、`--hi`、`--mean`、`--zero-prob`、`--wide-log2` |
| 形状 | `--n`、`--avg-degree`、`--degree`、`--attach`、`--rows`、`--cols`、`--exponent`、`--multi keep\|drop\|switch` |
| Road | `--jitter`、`--keep-prob`、`--shortcut-prob`、`--shortcut-radius`、`--int-scale` |
| Geometric | `--radius`(默认由 `--avg-degree` 得出)、`--geo-weights random\|euclidean`、`--int-scale`(与 road 共用:同时设置两者) |
| 其他 | `--connected`、`--seed`、`--stats`(把 n、m、分量数输出到 stderr)、`--list`、`--help` |

退出码:成功为 0,用法错误(未知选项、数字格式错误)为 2,生成或 IO 错误(参数无效、未知类型、文件不可写)为 1。相同的参数始终产生逐字节相同的输出。

## 9. 已知限制

- 生成器要求 n < 2^32,即 32 位的 `Vertex` 类型。点对键把两个顶点打包进 64 位。
- `chung_lu` 要求 m = round(n·avg_degree/2) ≤ ⌊⌊n(n−1)/2⌋ / 2⌋(整数除法),即至多是所有顶点对的一半,这样重新抽取重复边才能很快终止。否则它会抛出
  `std::invalid_argument`。在 1000 + 100·m 次被拒绝的抽取之后,它会抛出 `std::runtime_error`。
- 使用 `kSwitchToSimple` 的 `random_regular`,当不存在简单的
  d 正则图时(n·d 为奇数,或 d ≥ n 且 d > 0)会抛出 `std::invalid_argument`。对于合法的参数,交换可能
  无法在 1000 + 100·(n·d/2) 次提议之内完成。有些配置永远无法修复,
  例如 n = 3、d = 2 时的三个自环。出现这种情况时,生成器会从同一条带种子的流中抽取一个新的配置并重试。只有 100 个配置都失败时,它才会抛出 `std::runtime_error`。对于每个 (n, d)(n ≤ 12),在 300 个种子中单个配置最多失败了 31 次,因此这种情况尚未出现过。稠密的情形(d 接近 n)很慢,
  因为大多数提议都会被拒绝。
- `random_geometric` 有 m ≈ n²·F(r)/2 条边,因此半径不小时(例如 r ≥ √2 时的完全图)输出会随 n 呈二次增长。稀疏图请使用 `--avg-degree`。
- `DimacsArcMode::kMergeReversePairs` 为未配对的弧保存一个哈希表,因此在大型道路网络上,其内存与弧的数量成正比增长。
- `det_log` 和 `det_exp` 的精度为几个 ulp,并非正确舍入。它们在各平台之间可以精确复现,但与 libm 并不逐位相等。
- 整数权重参数必须不超过 2^53,这样经由 `double` 的转换才是精确的。
