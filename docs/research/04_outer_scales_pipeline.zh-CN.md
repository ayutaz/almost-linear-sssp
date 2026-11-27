# 04. 外层尺度与外层流水线（论文 Sec 8.2–8.3）：实现规格

[English](04_outer_scales_pipeline.md) | [日本語](04_outer_scales_pipeline.ja.md) | **简体中文**

目标：arXiv 2611.34911v1，Sec 8.2（Lemma 8.4，l.2076–2142）与 Sec 8.3（Lemma 8.5 l.2145，Algorithm 8 l.2167–2196，Theorem 1 的证明 l.2158–2217）。行号指 `.cache/paper/2611.34911v1.txt`。
对应任务：**P3-OUTER-SCALES**（`src/paper/outer/scales`）与 **P3-OUTER-PIPELINE**（`src/paper/outer/shortest_paths` 中的 `outer_construction`）。
本笔记承接 `03_pivots_aux_outer.md` §3–§4，并把它整理成实现者可以直接照做的规格。

记号：`[Paper]` 是论文明确陈述的内容。`[Implied]` 是从论文陈述出发、经此处给出的简短论证得到的内容。`[Interpretation]` 是我们对含糊之处的解读（每条在 §5 中有对应的 Q 提案）。`[Supplement]` 是论文未陈述的补充。`[Measured]` 是为本笔记计算或运行得到的结果（大整数脚本，或把 `main` 上的组件链接进临时程序；其中没有任何内容被提交）。

本笔记的范围限制：
- P3-OUTER-STATES（Sec 8.1）的产出尚不在 `main` 上。§2.0 把尺度阶段需要从它得到的内容列为 **假设 A-1**；该任务合入后再对齐名称。
- 这里没有任何东西做过端到端运行：外层构造尚不存在。只运行过到径向状态为止的阶段（§3.5）。

---

## 0. 要点

1. `F = 2^f`，其中 `f = p + 23 + bit_width((2h+1)·h²)`，`P = 2^p`，`h = b + 1`，`K = 16`。只保留指数 `f`；`F²s` 与 `F³a` 通过 `2f` 次和 `3f` 次饱和加倍来形成（D-13）。`h` 必须是**精确的** `b + 1`，不能是饱和后的字（§1.3）。
2. 在 `uint64` 权重与论文的 `P` 下（于是 `f ≥ 45`），对每个 `s ≥ 1`，`F²·s` 都会饱和：**每个已挂接的块都是对齐的，所有支撑合并为一张网格**。在 `uint64` 中只有 `f ≤ 31`（Mode B，`P ≤ 8`，`b ≤ 3`）时才能到达未对齐分支，否则只有在 `double` 且动态范围超过 `2f` 个倍频程时才可能（§1.5）。
3. Sec 8.2 复用的 Sec 5.4 网格 / 合并 / 链接机制位于 `src/paper/local/assemble.cpp` 的匿名命名空间中，并且它为小分量构建精确的 Kruskal 树，而 Sec 8.2 需要的是**划分层次结构的原子叶子**（Def 5.12）。它不能原样调用（§4，M-6，M-7）。
4. 对 `N < 2^64`，Mode A 无法到达 Alg 8 l.8–17（Q-04）。同一份代码用 Mode B 参数来运行。最小的强制设置是：外层阈值 `0·L^0` 与常数分隔符 `b`；其余全部取论文的值（§3.4）。
5. 外层分支变得非平凡的最小收缩后 `N`［Measured，§3.5］：`N = 2`（两个枢轴，`b = 1`），`N = 4`（一个非单点的已挂接块，`b = 2`），**`N = 6`（一个非单点原子，于是 Sec 8.2 的每个阶段都做实际工作，`b = 3`）**。
6. 目前声明的 `outer_construction(c, diag)` 没有生成器也没有参数；`PaperOptions` 没有种子。在能够调用 Alg 8 l.9 之前，两者都是必需的（§4，M-3；Q-70）。

---

## 1. 因子 F（Sec 8.2 l.2080–2081）

### 1.1 定义

[Paper] `F` 是严格大于 `X = 2·P·(2h+1)·(8K²h)²` 的最小的 2 的幂，每个图用整数参数选定一次。`P` 是 NameRoutes 的公因子（Alg 5 l.3），`h = b + 1` 是外层属主分解的跳数上界（Sec 7.1 l.1644），`K = 16`（Sec 8.1）。

[Implied] 各因子的来历：`2P(2h+1)` 界定共享同一名称的两个宏权重之比（Lemma 7.1：`λ/(2h+1) ≤ w < 2Pλ`），`4K²h` 是 Lemma 8.2 中的比值 `θ/B0`。对于带新名称 `a` 的有限连接（§2.1）：`a ≤ B0 ≤ θ < 8P(2h+1)K²h·a = X·a/(8K²h)`，因此 `θ < F·a`，并余下 `8K²h` 的富余因子。

### 1.2 精确的（值，指数）形式（Q-06）

`F` 是一对（值 1，指数 `f`）。`F²` 是 (1, `2f`)，`F³` 是 (1, `3f`)。`F` 从不被实体化为长度或机器整数。

[Implied] 满足 `2^f > X` 的最小 `f` 是 `bit_width(X)`。取 `K = 16` 与 `P = 2^p`：

```
X = 2 · 2^p · (2h+1) · (2^11 · h)²  =  2^(p+23) · (2h+1) · h²
f = p + 23 + bit_width((2h+1) · h²)
```

[Measured] 在 §1.4 的每一行上，闭式都等于 `bit_width(X)`。

### 1.3 无溢出地计算 f

`f` 是对规模（size）做的整数字运算（Sec 2 l.302–304），从不对长度运算。

```cpp
// Inputs: p_log2 (unsigned), h as an exact natural number (multiword).
// Y = (2h + 1) * h * h      -- three multiword operations; Y < 2^(3·bits(h) + 2)
// f = p_log2 + 23 + bit_width(Y)
```

- `h` 必须精确。在 Mode A 中，`b = L_N^64` 对任何 `N` 都放不进一个字（`L_N ≥ 2`），而 `owner::outer_separator_size` / `OwnerDecomposition::hop_bound` 保存的是**饱和**值 `2^64 − 1`。由饱和后的 `h` 推出 `f` 是错的：在 `N = 10^6` 时会得到 `f = 233` 而不是 871［Measured］。要用多字算术来形成 `L^64` 与 `h`（64 次字乘多字；`L ≤ 65`，故 `h < 2^386`）。
- `src/experimental/paper_experimental/params.cpp` 已经这样做了（`OuterDerived::f_log2`，匿名命名空间中的类 `Nat`）。`src/paper/` 不能包含 `src/experimental/`；见 M-2。
- 对每个可表示的输入，`f`、`2f`、`3f` 都放得进 `unsigned`（`f ≤ 29 + 23 + 3·386 + 2`）。

把 `F` 作用于长度（D-13，Q-38；`src/paper/local/saturating.hpp`）：

```cpp
// aligned(s, a)  :=  !(sat::mul_pow2(s, 2 * f) < a)        // F² s ≥ a, Sec 8.2 l.2086
// support_end(a) :=  sat::mul_pow2(a, 3 * f)               // F³ a,    Sec 8.2 l.2081
```

[Implied] 作为判定，两者都是精确的。`sat::mul_pow2` 恰在真实乘积放不下时返回 `WeightTraits<W>::infinity()`，而此时真实乘积至少为 `infinity()`，超过该类型的每个有限长度。所以对每个有限的 `a`，`F²s ≥ a` 的真值与精确算术中的相同；并且饱和后的支撑端点与后面每个支撑相交，当且仅当真实端点如此（后面的支撑从有限值开始）。`sat::mul_pow2` 在到达无穷时停止，因此整数的循环至多执行 `min(2f, 字位数)` 次，`double` 至多 2098 次。
这里不要使用 `LenOps::mul_pow2`：它在溢出时抛出异常而不是饱和。

### 1.4 指数［Measured，Python 大整数］

`u = ⌈log2(t+2)⌉`，`g = 4u+3`，`P = 2^p` = 不小于 `2048·g·u²` 的最小的 2 的幂（Alg 5 l.2–3）。

论文规模与 Mode A 的行（`b = L_N^64`）：

| N | b | t | p | f | 2f | 3f |
|---|---|---|---|---|---|---|
| 589^64 + 1（最小的外层 N） | 589^64 ≈ 2^588.94 | 1（唯一可能的值） | 17 | **1808** | 3616 | 5424 |
| 589^64 + 1（仅公式） | 589^64 | 66 | 22 | 1813 | 3626 | 5439 |
| 589^64 + 1（仅公式） | 589^64 | 2^20 | 27 | 1818 | 3636 | 5454 |
| 590^64 + 1 | 590^64 ≈ 2^589.09 | 1 | 17 | 1809 | 3618 | 5427 |
| 10^3（在 Mode A 中不可达） | 10^64 | 1 | 17 | 679 | 1358 | 2037 |
| 10^6（在 Mode A 中不可达） | 20^64 | 1 | 17 | 871 | 1742 | 2613 |
| 2^63（在 Mode A 中不可达） | 64^64 = 2^384 | 1 | 17 | 1194 | 2388 | 3582 |

Mode B 的行（`t` 列是 `thresholds.md` §5.2 所用的上界 `2⌊N/b⌋ − 1`）：

| N | b | t | p | f | 2f | 3f |
|---|---|---|---|---|---|---|
| 10^6 | 20 | 99999 | 26 | 64 | 128 | 192 |
| 10^6 | 256 | 7811 | 25 | 74 | 148 | 222 |
| 10^6 | 400 | 4999 | 25 | 75 | 150 | 225 |
| 10^6 | 4096 | 487 | 23 | 84 | 168 | 252 |
| 10^6 | 8000 | 249 | 23 | 86 | 172 | 258 |
| 10^6 | 160000 | 11 | 20 | 96 | 192 | 288 |
| 10^8 | 27 | 7407405 | 27 | 66 | 132 | 198 |
| 10^8 | 729 | 274347 | 26 | 79 | 158 | 237 |
| 10^8 | 4096 | 48827 | 26 | 87 | 174 | 261 |
| 10^8 | 19683 | 10159 | 25 | 92 | 184 | 276 |
| 10^8 | 531441 | 375 | 23 | 105 | 210 | 315 |

小分隔符（强制外层的测试）。取论文的 `P` 时，`t = 1..2`（`p = 17`）、`t = 6`（`p = 19`）、`t = 66`（`p = 22`）的 `f`，以及取 D-17 的值 `P = 1`、`P = 2` 时的 `f`：

| b | h | bit_width((2h+1)h²) | f (t ≤ 2) | f (t = 6) | f (t = 66) | f (P = 1) | f (P = 2) |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 5 | 45 | 47 | 50 | 28 | 29 |
| 2 | 3 | 6 | 46 | 48 | 51 | 29 | 30 |
| 3 | 4 | 8 | 48 | 50 | 53 | 31 | 32 |
| 4 | 5 | 9 | 49 | 51 | 54 | 32 | 33 |
| 8 | 9 | 11 | 51 | 53 | 56 | 34 | 35 |
| 16 | 17 | 14 | 54 | 56 | 59 | 37 | 38 |
| 33 | 34 | 17 | 57 | 59 | 62 | 40 | 41 |
| 64 | 65 | 20 | 60 | 62 | 65 | 43 | 44 |

这些与 `tests/unit/paper_experimental/test_experimental_params.cpp` 中已测试的值（871；`b = 160000`、`t = 11`、`P = 2` 时为 77）以及 `thresholds.md` §5.1–5.2 一致。

### 1.5 饱和对各权重类型意味着什么［Implied］

| 类型 | `F²·s ≥ a` | 支撑 `[l, F³a]` |
|---|---|---|
| `uint64`，`2f ≥ 64`（取论文的 `P` 时总是如此） | 对每个 `s ≥ 1` 都为真：每个具有有限挂接的块都是**对齐的** | 对每组参数 `3f ≥ 84 > 64`：每个支撑端点都饱和，所以**所有支撑合并为一张网格** |
| `uint64`，`2f ≤ 63`，即 `f ≤ 31` | 可能为假：需要 `p + bit_width((2h+1)h²) ≤ 8`，即 `b = 1, p ≤ 3`；`b = 2, p ≤ 2`；`b = 3, p = 0` | 同上 |
| `double`，`f ≤ 1022` | 当 `a > 2^(2f)·s` 时为假（跨越超过 `2f` 个倍频程的二进制权重） | 支撑相距超过 `3f` 个倍频程时，网格分量彼此分离 |
| `double`，论文规模（`f ≥ 1808`） | `2^3616` 超出整个 `double` 的范围（最小次正规数到最大有限值之间为 2^2098）：总是对齐 | 总是一张网格 |

对测试的推论：未对齐分支与多分量网格，要么需要在 `uint64` 上使用 `p_log2 = 0` 且 `b ≤ 3` 的 Mode B，要么需要使用二进制权重的 `double`（§6）。

---

## 2. Lemma 8.4，逐阶段说明（Sec 8.2 l.2076–2142）

### 2.0 输入

来自 `main` 上现有组件：

| 输入 | 类型 | 来源 |
|---|---|---|
| 收缩后的图 | `CsrGraph<W>` | `outer::ZeroContraction<W>::graph` |
| 属主分解（外层 `b`） | `OwnerDecomposition<W>`：`radius` (r_v)、`owner`、`pivots`、`corridors`、`hop_bound` | `owner_decomposition(g, b)` |
| 分类 | `owner::OuterEdgeClassification<W>`：`edge_class`、`name`、`dictionary`、`p_log2` | `owner::classify_interowner_edges` |
| 辅助图 | `owner::AuxGraph<W>`：`aux_vertex`、`radius`、`owner_distance` | `owner::build_aux_graph` |
| 日历（用于分离的网格） | `owner::Calendars<W>`：`graph_piece`、`radius_index`、`graph_pieces[].origin`、`.grid` | `owner::OwnerDistancesResult<W>::calendars` |

**假设 A-1**（P3-OUTER-STATES，不在 `main` 上；字段名为暂定）。尺度阶段需要：

| 需要 | 含义 |
|---|---|
| `StateGraph<W> graph`、`projection` | 每条状态弧（普通与特殊），以及每个状态的物理像 |
| 每个状态的 `block`、`num_blocks`、块种类 | 种类 ∈ {正的普通块, 枢轴副本, 特殊}；特殊状态与枢轴副本是单点块（Sec 8.1 l.1996） |
| `arc_kind`、`arc_edge` | 用于识别窗口内的弧及其输入边 |
| 每个顶点的 `native_inward` | Alg 8 l.14 的种子；对枢轴而言是它的向内枢轴副本 |
| 每个状态的 `theta`、`b0`、`b0_edge` | Lemma 8.2 的值与达到它的宏边（“有限的 B0 保存达到它的宏边的指针”，l.2078–2079） |
| 连接记录 | 块 MST 边（取其权重），状态到属主的连接（取 θ），走廊连接（取走廊权重，Lemma 8.3），形如 `ThresholdJoins<W>::joins` |

如果 P3-OUTER-STATES 复用 `build_radial_states`，它可以传入 `EdgeClassification::macro_key[e] = cls.name[e]`；这样 `ThresholdJoins::b0_key` 就是达到该值的宏边的名称在字典中的索引，这恰好是 §2.1 所需要的。
单位：Alg 7 还会返回 `scaled_distance = 2^7·J`（D-21）。本节中的一切都用**输入单位**（r、θ、MST 权重、名称）。

### 2.1 阶段 S0 —— 新名称（l.2077–2079）

[Paper] 把每条宏边的名称替换为携带该名称的所有宏边中最小的物理权重。

```cpp
// in : g, cls (edge_class, name, dictionary)
// out: std::vector<W> name_weight   // per dictionary index; infinity() if no macro edge has it
for (EdgeId e : edges with cls.edge_class[e] == kMacro)
  name_weight[cls.name[e]] = min(name_weight[cls.name[e]], g.edge(e).weight);
```

- [Interpretation, Q-65] “物理权重”指输入权重 `w_e`（Lemma 7.1 就是针对 `w` 陈述的），而不是投影长度 `A_e`；并且最小值只在宏边上取。
- 现有函数：无；一次线性扫描。长度只做比较。
- 不变量（可测试）：对每条宏边 `e`，设 `a = name_weight[name[e]]`：`a ≤ w_e`；`w_e < 2P(2h+1)·a`（饱和乘积）；`a` 是带该名称的某条宏边的权重。对每个 θ 有限的状态 `z`，设 `a_z = name_weight[name[b0_edge[z]]]`：`a_z ≤ b0[z] ≤ theta[z]` 且 `theta[z] < F·a_z`。

### 2.2 阶段 S1 —— F

`f = outer_scale_factor_log2(p_log2, exact h)`（§1.3），每个图一次。[Supplement] 当没有宏边（`t = 1`）时，`f` 会被计算但从不被使用。

### 2.3 阶段 S2 —— 块尺度与挂接（Sec 8.1 l.1999；Sec 8.2 l.2083–2086）

```cpp
// in : states (A-1), od.radius, name_weight
// out: std::vector<OuterBlockScale<W>> blocks      // per block
struct OuterBlockScale {
  Kind kind;            // kSingleton, kAligned, kUnaligned, kUnattached
  std::size_t size;
  W s;                  // s_U = min over states x of U of od.radius[projection[x]]   (Q-35)
  Vertex s_vertex;      // a vertex attaining s_U (least id) -- needed for the separate grid
  StateId z;            // a state minimising theta over U (least id); kNoState if all infinite
  std::uint32_t name;   // dictionary index of the name of b0_edge[z]; a_U = name_weight[name]
  W t;                  // t_U, set in S5
};
```

- [Paper] 单点块（每个枢轴副本、每个特殊状态，以及任何只有一个状态的普通块）只保留它自己的有限 θ 连接，并跳过块尺度测试、网格、MST 分类与原子替换（l.2083–2085）。
- [Paper] 对具有有限挂接的非单点正普通块，选取一个取最小值的状态及其名称 `a_U`（l.2085–2086）。
- [Interpretation, Q-35] `s_U` 使用原始半径 `r_v = b_T(v, p(v))`，而不是 `J`。
- [Supplement] 并列情形（θ 有多个最小者、半径为 `s_U` 的顶点有多个）取编号最小者。任何选择都有效：每个最小者 `z` 都满足 `a_z ≤ θ(z) < F·a_z`。
- 不变量：对 `U` 的每个状态 `x`，`0 < s_U ≤ r_x`；对每个 `z ∈ U`，`s_U ≤ θ(z)`（Lemma 8.2 `r ≤ θ`）；`θ(z_U) = min_U θ`；块内每个物理顶点只有一个副本（块的状态单射地投影）。

### 2.4 阶段 S3 —— 对齐或未对齐（l.2086–2087；Lemma 8.4 证明 (1)、(2)）

```cpp
// in : blocks (kind not yet split), name_weight, f
// out: blocks[U].kind in {kAligned, kUnaligned} for attached nonsingleton blocks
// group the attached nonsingleton blocks by blocks[U].name; in a group with name weight a:
//   s_star = pivots::median_prune(keys = {s_U}, feasible(s) = !(sat::mul_pow2(s, 2f) < a))
//   U is aligned  iff  s_star exists and !(s_U < *s_star)
```

- [Paper] 对齐当且仅当 `F²·s_U ≥ a_U`。论文按名称把测试分组，并使用中位数剪枝，每个被测中位数做 `O(log F)` 次加倍（l.2130–2131）。
- 现有函数：`pivots::median_prune`（`src/paper/pivots/name_routes.hpp`；该谓词关于 `s` 单调）、`sat::mul_pow2`。按名称分组是在字典索引上的计数排序。
- [Supplement] 逐块测试给出相同的分类，每块至多做 `2f` 次加倍；它只改变工作量界。若使用，请把它记为仅涉及时间的偏离。
- 不变量：
  - 未对齐块：每条块 MST 边（实际上 `U` 的每条窗口内的弧）的权重都 `< a_U/2`，以 `sat::add(w, w) < a_U` 测试。[Paper, l.2086–2087] 理由：`w ≤ 4K²·J_v ≤ 4K²h·s_U < a_U/2`。
  - 未对齐块：`U` 的所有状态的 θ 相同（Lemma 8.4 证明 (1)，“equal attachment threshold”）。
  - 对齐块：`s_U ≤ θ(z_U) < F·a_U` 且 `a_U ≤ F²·s_U`。
  - 差分：分组得到的结果等于逐块测试 `!(sat::mul_pow2(s_U, 2f) < a_U)`。
  - 比较次数：对大小为 `n_a` 的每个名称组，谓词测试至多 `⌊log2 n_a⌋ + 1` 次。

### 2.5 阶段 S4 —— 支撑与全局网格（l.2081–2083、l.2087–2091；Lemma 8.4 证明 (2)）

键：每个为某个有限连接命名的字典索引一个，`T_A` 的每条走廊一个（`od.corridors`）。

```cpp
// in : name_weight, od.corridors, aligned blocks, theta joins, f
// out: std::vector<W> grid;                       // increasing grid points of all components
//      std::vector<std::uint32_t> component_of;   // per grid point
//      std::vector<GridSupport<W>> support;       // per key: lower, upper, lower_index, upper_index
// 1. lower[key] = a (name weight, or the corridor weight); upper[key] = sat::mul_pow2(a, 3f)
// 2. for every aligned U: lower[name_U] = min(lower[name_U], s_U)
// 3. sort the supports by lower; merge two neighbours when they overlap or when
//    next.lower < 2 * current_component.upper (gap of ratio below two)
// 4. per merged component: grid points lo, 2lo, 4lo, ... until the first point strictly above the
//    component's upper end (one final doubling as sentinel); stop at infinity()
// 5. per key: lower_index = largest grid point <= lower; upper_index = first grid point > upper
```

- [Paper] 支撑 `[a, F³a]`；对对齐块把起点降到 `s_U`；合并相交的支撑以及间隙小于二倍的支撑；索引从 `l` 的网格下取整直到严格大于 `u` 的第一个点；对 `k = O(t log²(N+2))` 个名称，共 `O(k log F)` 个倍频程。
- [Interpretation, Q-65] 走廊连接使用它自身权重的支撑，而不是它的宏名称的支撑。
- [Interpretation, Q-68] “间隙小于二倍”是严格测试 `next.lower < 2·upper`。经此合并之后，没有任何类能跨越分量的起点，所以每个类恰为 `[g, 2g)`（不同于 Sec 5.4，那里一个类可以止于下一个分量，l.1279）。
- [Supplement, 同 `assemble.cpp`] 高于最大有限连接权重的网格点不必生成：那里没有任何东西被分类。
- 现有函数：无可调用者。相同的步骤存在于 `src/paper/local/assemble.cpp`（匿名命名空间）中，以 `N^16σ` 代替 `F³a`；见 M-6。
- 不变量：
  - 每个有限的 `θ(z)` 都落在它自己的键的 `[lower, upper]` 之内；每个走廊权重都落在它自己的支撑之内。
  - 对键为 `κ` 的对齐块 `U`：`lower[κ] ≤ s_U`，并且 `U` 的每条满足 `w ≥ t_U` 的块 MST 边都满足 `grid[lower_index[κ]] ≤ w < grid[upper_index[κ]]`（Lemma 8.4 证明 (1)：`w ≤ 4K²h·s_U < F²·a_U`）。
  - 对对齐块 `U` 与每个 θ 有限的 `z ∈ U`：`name(z)` 的支撑与 `name_U` 的支撑位于同一个合并后的分量中。[Implied] `s_U ≤ θ(z) ≤ max(θ(z_U), 4K²h·s_U) < F²·a_U`，所以 `θ(z)` 同时落在两个支撑内。
  - 网格：严格递增；分量内部 `grid[i+1] = 2·grid[i]`；分量之间 `grid[下一个分量的第一个点] ≥ 2·（前一个分量的上端）`。
  - 规模：合并前每个输入支撑至多 `3f + 2 + log2(降低量)` 个点，其中降低量至多为 `F²`（对齐时 `s_U ≥ a_U/F²`），所以至多 `5f + 2`。

### 2.6 阶段 S5 —— 截断值 t_U（l.2091–2093、l.2100–2104）

```cpp
// aligned U (key k):   t_U = largest grid point <= s_U, searched in [lower_index[k], upper_index[k]]
// every other nonsingleton positive ordinary block (unaligned, or without finite attachment):
//   v  = blocks[U].s_vertex;  x = aux.aux_vertex[v];  C = cal.graph_piece[x]
//   i  = cal.radius_index[x]              // cal grid is 2^7-scaled: grid[i] <= 2^7 r_v < 2 grid[i]
//   t_U = sat::mul_pow2(cal.graph_pieces[C].origin, i - 7)     // origin * 2^(i-7), input units
```

- [Paper] 两种选择都给出 `s_U/2 < t_U ≤ s_U`（Q-36：左侧为严格不等号，不同于 Sec 5.4）。
- [Paper] 对分离的网格：选取达到 `s_U` 的 `v`；块中每个 `J` 至多为 `K²·J_v`；最大的局部边至多为 `4K²·R_v`，最不显著的边超过 `r_v/2`；使用 `v` 的辅助网格指针，每个片段的表格一次性多延伸固定个数的倍频程。
- [Interpretation, Q-67] “辅助网格指针”指 `v` 的 Sec 7.3 图片段网格（`Calendars`）。该网格以 `2^7` 的倍数缩放存放，`grid[7] = 2^7·origin`，所以输入单位下 `r_v` 之下的网格点为 `origin·2^(i−7)`，`i ≥ 7`（因为 `r_v ≥ origin`）。它通过加倍到达；不需要减半。
- [Interpretation, Q-66] “其他每个非单点正普通块”包括没有有限挂接的块（如 Sec 5.4 l.1213）。
- [Supplement] `U` 的分离网格的类是 `[t_U·2^j, t_U·2^(j+1))`，`j ≥ 0`，直到最大的块 MST 边所在的类（`≤ 4K²·R_v`，所以 `j ≤ 10 + ⌈log2 h⌉`）。这些就是片段网格的点。现有的片段网格只保证到达 `4·max R`；按论文所述共享同一张表，需要在顶端多延伸 9 次加倍（`K²` 占 8 次，哨兵占 1 次）。按块从 `t_U` 开始加倍来生成这些点，得到相同的类，并且每块花费 `O(log(K²h))` 次加倍（仅涉及时间；M-8）。
- 不变量：`t_U ≤ s_U < 2·t_U`（以 `!(s_U < t_U) && s_U < sat::add(t_U, t_U)` 测试）；对齐：`t_U` 是 `name_U` 所在分量的网格点，且 `grid[lower_index] ≤ t_U`；分离：`t_U = origin·2^j`，其中 origin 属于 `s_vertex` 所在的片段。

### 2.7 阶段 S6 —— 原子（l.2000–2001、l.2133–2142；Def 5.11）

```cpp
// in : block MST join records, blocks[U].t
// out: per state: part id; std::vector<ReplaceableAtom<W>> atoms
// DisjointSets over the states; unite the endpoints of every block MST edge with weight < t_U.
// A component with one state is a singleton part. A component Y with >= 2 states is an atom:
//   states        = members of Y
//   edges (H)     = every within-window arc pair with both endpoints in Y, once per input edge,
//                   with its physical weight (heavier chords included)
//   spanning_tree = the block MST edges < t_U inside Y
//   cutoff        = t_U
```

- [Paper] 由权重小于 `t_U` 的边构成的分量被局部表示所取代；每个原子位于 `T − A` 的一个分量内，且顶点数少于 `b`；在原子内部，属主与原始半径为常数（Lemma 8.4 证明 (3)）。
- [Interpretation, Q-53 扩展；Q-74] 内部图包含两端都在该原子内的每条窗口内的边。`substitute_atoms` 会拒绝 `D` 中两端都在同一原子内的任何弧，所以这些边必须全部放入 `H`。
- [Implied] 其他种类的弧不可能两端都在同一原子内：转移弧与切换弧离开窗口或阶段，枢轴弧与特殊弧涉及单点块，而宏边连接两个属主单元格，而一个原子只有一个属主。检查这一点，违反时抛出 `std::logic_error`。
- [Interpretation, Q-74] 只有一个状态的分量是单点部分（Def 5.12，Lemma 5.15 “every nonsingleton part”）；不为它们调用 Alg 4。
- 现有函数：`DisjointSets`（`src/paper/primitives/dsu.hpp`）；类型 `ReplaceableAtom<W>`、`kNoAtom`（`src/paper/local/atom_replacement.hpp`）。
- 不变量：`|Y| < b`；`Y` 上有一个属主、一个 `r`，且 `r ≥ s_U ≥ t_U`；`Y` 位于一个块与 `T − A` 的一个分量之内；`Y` 单射地投影；每条生成树边都 `< t_U`（所以 `< 2·cutoff`）；`H` 连通；每个 `z ∈ Y` 的有限 `θ(z)` 都 `≥ t_U`。

### 2.8 阶段 S7 —— 事件分类（l.2091、l.2094；Lemma 8.4 证明 (2)）

| 事件 | 类 |
|---|---|
| 块 MST 边 `w < t_U` | 无（在原子内部） |
| 块 MST 边 `w ≥ t_U`，`U` 对齐 | 全局类 `i`：`grid[i] ≤ w < grid[i+1]`，在 `name_U` 的索引范围内做二分查找 |
| 块 MST 边 `w ≥ t_U`，`U` 未对齐或未挂接 | `U` 的分离类 `j`：`t_U·2^j ≤ w < t_U·2^(j+1)` |
| 在有限 θ(z) 处的状态到属主的连接 | 在 `name(z)` 的索引范围内的全局类 |
| 走廊连接 | 在其走廊键的索引范围内的全局类 |

- [Paper] `O(N)` 个事件，每个 `O(1 + log log F) = O(τ(N))` 次比较。
- 不变量：所找到的类括住该权重；每个事件的二分查找步数 `≤ ⌈log2(upper_index − lower_index + 1)⌉`；对未对齐的 `U`，设 `p_top` 为它的最高分离类的下端，`g` 为 `θ(z_U)` 的全局类：`2·p_top ≤ g`［Implied：`2·p_top ≤ 2·w_max < 8K²h·a_U/F² < a_U/2 < g`］，所以整个块在它的第一个全局事件之前是同一个节点（l.2093–2094）。

### 2.9 阶段 S8 —— 划分层次结构（l.2094–2099；Def 5.12；Lemma 5.10 后处理）

```cpp
// in : events with classes, parts (singletons and atoms), the state graph
// out: PartitionHierarchy<W> ph   // graph = outside arcs only; nodes; leaf_atom; atoms; c0 = 2
// 1. leaves: one per singleton part and one per atom
// 2. join vertices: states and one auxiliary vertex per pivot; DisjointSets (Lemma 3.6)
// 3. per separate-grid block, classes in increasing j; then global classes in increasing index.
//    All events of one class are merged simultaneously; the node's height is the class lower
//    boundary. Auxiliary nodes without states are dropped and unary nodes compressed.
// 4. links: every outside arc is assigned to the pair of children it joins at the LCA of its
//    endpoints' parts; at a node of class [g, 2g) the assigned arcs of weight < 2g form the child
//    adjacency; a DFS takes one arc per discovery edge (endpoints are the actual states).
// 5. order: at each node, topological order of the children by the assigned arcs of weight < g;
//    roots in condensation order of the arcs between roots.
```

- [Paper] 按整数网格类处理并集、抑制辅助节点、并对比每个节点更轻的弧排序，就得到一个划分层次结构；连接树来自 Sec 5.4 的子节点商搜索；对类 `[a, 2a)`，Lemma 8.3 的严格形式给出由严格小于 `2a` 的弧所形成的连通性。
- [Interpretation, Q-68] 高度 = 类的下边界，链接 = 低于上边界的弧，因此 `c0 = 2`，与 Lemma 5.10 对 Sec 5.4 的处理完全一致。“先分离网格、后全局类”的顺序是我们的；任何与每条根路径上各高度一致的顺序都等价。
- 现有函数：`DisjointSets`、`ForestQueries`（LCA 与层级祖先）、`add_internal_node`、`collect_roots`；类型 `PartitionHierarchy<W>`。`assemble.cpp` 的合并构建器与链接 / 排序步骤不能调用（M-6，M-7）。
- 不变量：
  - `substitute_atoms` 接受 `ph`（各部分划分状态集，`D` 中没有弧在原子内部，每个原子的截断值 `≤` 其父节点的高度）。
  - 对高度为 `g`、类为 `[g, g')` 的每个内部节点 `X`：`X` 之下的状态是长度 `< g'` 的状态弧的一个 SCC，且 `X` 的子节点是其中长度 `< g` 的弧的 SCC（Lemma 8.3 的严格形式；在小输入上用 Tarjan 检查）。
  - 高度朝根方向不减；每条链接的长度 `< 2·height`。
  - 每个原子的父节点高度 `≥ t_U`。

### 2.10 阶段 S9 —— 原子替换（l.2099；Lemma 5.15）

```cpp
AtomSubstitution<W> sub = substitute_atoms(ph, &diag, params.local);   // one Alg 4 call per atom
```

- 现有函数：`substitute_atoms`（`src/paper/local/atom_replacement.hpp`）。它对每个原子调用 `local_hierarchy`，当 Alg 4 进入其基本情形时，会推入 `FallbackRecord{"local: k <= 2*L_k^4 (paper base case)", k, m}`。
- 不变量：`validate_ordered_hierarchy(sub.hierarchy, contracted graph)` 在 `c0 ≤ 2` 下不报告任何违规（`docs/TASKS.md` 中的验收标准）；`sub.old_state` 与 `projection` 的复合给出每个新状态的物理像；每个原子恰有一次 `local_hierarchy` 调用（`atom_replace_atoms`）。

### 2.11 声称的复杂度

[Paper, Lemma 8.4] 对显著的局部事件与有限连接（包括所有全局网格）的分类，工作量为 `O(m + N·τ(N))`。Lemma 5.15 为替换另加 `O(m + N log log(b+4))`。

---

## 3. Algorithm 8 l.8–17（Sec 8.3）

### 3.1 阶段、输入与输出

`c` 是 Alg 8 l.2 得到的 `outer::ZeroContraction<W>`（权重为正、无环、连通、保留平行边），`g = c.graph`，`N = g.num_vertices()`。

| Alg 8 | 步骤 | 调用（除非另有标注，均在 `main` 上） | 输入 → 输出 |
|---|---|---|---|
| l.8 | MST（Lemma 3.5）、属主单元格、删除被树严格支配的非树边 | `owner_decomposition(g, b, &diag)` | `g`、`b` → `OwnerDecomposition<W> od`（`tree_edges`、`reduced`、`pivots`、`owner`、`radius`、`owner_distance`、`corridors`） |
| l.9 | 对投影后的属主间边运行 NameRoutes | `owner::classify_interowner_edges(g, od, rng, name_routes_params, &diag)` | → `OuterEdgeClassification<W> cls`（`name`、`dictionary`、`p_log2`、`detour` d_e） |
| l.10 | 宏 / 径向；`G'`、`T'` | （分类是上面那次调用的一部分）`owner::build_aux_graph(g, od, cls, &diag)` | → `AuxGraph<W> aux` |
| l.11 | OwnerDistances（Alg 7） | `owner::owner_distances(aux, &diag, {.local = local_params})` | → `J`（`distance`、`scaled_distance`）、`calendars` |
| l.12 | 删除满足 `w_e > J_x + J_y + d_e` 的径向边 `e = xy` | **P3-OUTER-STATES** | → 幸存的边集 |
| l.13a | 外层状态、θ、连接（Sec 8.1） | **P3-OUTER-STATES** | → A-1 |
| l.13b | 名称、F、尺度、网格、划分层次结构（Sec 8.2） | **P3-OUTER-SCALES**（§2，S0–S8） | → `PartitionHierarchy<W> ph` |
| l.13c | 带截断的局部层次结构 | `substitute_atoms(ph, &diag, local_params)` | → `AtomSubstitution<W> sub` |
| l.14 | 种子 | `redirect_offers(sub, {{native_inward[c.source], 0}})` | → 在 `s` 的原生向内副本处的一个 `StateOffer<W>`；若 `s` 位于某个原子内，则在它的原生局部子副本处 |
| l.15 | Evaluate，逐顶点取最小值 | `evaluate(SeededOrderedHierarchy<W>{sub.hierarchy, {{c.source, 0}}, offers})` | → 每个新状态的 `EvalResult<W>::value`；`dist[v] = min`，取遍满足 `projection[sub.old_state[x]] = v` 的新状态 `x` |
| l.16–17 | 还原、紧弧、首次发现树 | `outer::restore_predecessor_tree(g_original, s, c, dist)`（`paper_sssp` 已经在调用） | → `SsspResult<W>` |

说明：
- l.8：论文中的“对每个顶点对保留一条最轻的边”只属于基本情形（l.5）。外层分支保留平行边；NameRoutes 自己对它的实例去重（Alg 5 l.4）。
- l.9 是唯一的随机化步骤（Lemma 6.8）。生成器是显式传入的 `model::WordRng&`。
- l.15［Interpretation, Q-73］：最小值取遍每个投影到 `v` 的状态，包括特殊状态与所有局部子副本。Lemma 4.6 保证不会低估，且至少有一个精确的副本。若某个顶点没有有限的发出值，就意味着不变量被违反：抛出 `std::logic_error`，绝不为源分量中的顶点返回无穷。
- `outer_construction` 只返回收缩后的距离；树不是由 Evaluate 产生的。

### 3.2 Lemma 8.5（l.2145–2157）

[Paper] 收缩零权边分量之后，收缩后的源分量的精确距离就确定所有原始距离与一棵最短路径树，耗时 `O(m + n)`：当 `d(u) + w(uv) = d(v) < ∞` 时，把 `uv` 从 `u` 定向到 `v`，从 `s` 出发遍历，保留首次发现的弧，并沿树向下输出原始权重之和。

这已在基础分支上实现并为 VERIFIED（`restore_predecessor_tree`，PAPER_MAPPING 中的 “Sec 8.3, Lemma 8.5” 行）。它只接受距离，所以 Alg 8 l.16 不需要新代码。外层分支必须遵守的是：
- 输入标签必须精确。当某个有限标签不是由一条紧的原始弧达到时，`restore_predecessor_tree` 抛出 `std::invalid_argument`。对 `double`，外层的距离是按另一种顺序求得的和，与树上的和不同，所以非二进制的输入在这里可能失败（Q-06，Q-55 选项 A，Q-63）。在 Mode B 中，这是**失败**，绝不是替换成另一种搜索的理由。
- 不可达的原始顶点位于 `c.vertex_map` 之外，在还原中得到无穷；外层构造从不会看到它们。

### 3.3 回退在哪里记录

| 位置 | 记录 | 由谁推入 |
|---|---|---|
| Alg 8 l.4–5（`N ≤ L_N^64`） | `kOuterBaseCaseScanReason` / `kOuterBaseCaseD1Reason`，`(N, m')` | `outer_base_case*`（已有） |
| Alg 7 l.6：对每个 Sec 7.3 原子运行 Alg 4 | `"local: k <= <c>*L_k^<e> (paper base case)"`，`(k, m_Y)` | 经 `cap_local_hierarchy` 的 `local_hierarchy`（已有） |
| Alg 8 l.13 / Lemma 5.15：对每个 Sec 8.2 原子运行 Alg 4 | 相同的原因，`(k, m_Y)` | 经 `substitute_atoms` 的 `local_hierarchy`（已有） |
| Alg 5 重试上限（D-10，仅 Mode B） | NameRoutes 的记录 | `name_routes`（已有） |
| `outer_construction` 本身 | **无** | 它从不用直接搜索取代任何阶段；不支持的输入会抛出异常 |

- [Measured] Alg 7 对单顶点原子也推入记录（`k = 1`，`m = 0`）：在 §3.5 的 6 顶点例子上有四条这样的记录。
- 一个顶点可能被多条记录覆盖：在 Alg 7 中一次，以及对包含它的某个副本的每个 Sec 8.2 原子各一次（至多四个副本：`in_i`、`in_{i−1}`、`out_i`、`out_{i+1}`）。在 6 顶点例子中，点对 `{x, y}` 在三个块里构成一个原子。因此 `FallbackRecord::vertices` 之和可能超过 `N`；回退占比需要一个定义（Q-69）。
- [Supplement] 提议的计数器（不是回退）：`outer_constructions`、`outer.pivots`、`outer.macro_edges`、`outer.radial_edges_deleted`、`outer.states`、`outer.blocks_{singleton,aligned,unaligned,unattached}`、`outer.atoms`、`outer.atom_states`、`outer.grid_points`、`outer.grid_components`、`outer.global_classes`、`outer.separate_classes`、`outer.scale_factor_log2`。
- 退化情形 `t = 1`（无宏边，每个块都未挂接，`J` 是从枢轴出发的完整距离）是论文的构造，不是回退；通过 `outer.pivots = 1` 来报告它（Q-72）。

### 3.4 Mode A 与 Mode B

Mode A：对每个 `N < 2^64`，`takes_outer_construction(N)` 都为假（`L_N ≥ 2`，所以 `L_N^64 ≥ 2^64 > N`）。Alg 8 l.8–17 无法通过 `paper_sssp` 到达；能到达它的最小 `N` 是 `589^64 + 1`（Q-04，`thresholds.md` §2）。因此这段代码只由 (i) Mode B 参数与 (ii) 直接调用 `outer_construction` 来运行。

[Implied] 在小图上以论文的常数直接调用 `outer_construction` 是有定义的：`b` 饱和到 `N` 之上，分解只有一个枢轴，没有宏边，每个非单点块都未挂接，Alg 7 把 `J` 算成从枢轴出发的完整距离。这是对 Mode A 实例化的一个有用测试（§6，PL-2），但它除未挂接路径之外并不触及 Sec 8.2。

强制外层分支做实际工作的最小 Mode B 参数集（`paper_experimental::PaperParams`）：

| 参数 | 值 | 论文 | 偏离 |
|---|---|---|---|
| `outer.threshold` | `SizeThreshold{0, 0}`（对每个 `N ≥ 1` 都适用） | `1·L^64` | D-7 |
| `outer.separator` | `SeparatorRule::constant(b)`，`b ∈ {1, 2, 3, …}` | `L^64` | D-7 |
| 其余全部 | 论文的值：样本 `16t`，余量 `32t`，论文的 `P`，论文的 `H`，无重试上限，局部阈值 `2·L^4`，局部 `b = L^4` | — | 无 |

推出：`h = b + 1`，`F` 来自 §1（例如 `b = 3`、`t = 2` 时 `f = 48`）。在这些取值下，每个外层原子的顶点数都少于 `b` 并进入 Alg 4 的基本情形（已记录），而且 NameRoutes 从不采样（当 `t ≤ 65` 时 `|U| ≤ 32t`），所以每条属主间边都证明它自己的权重，`d_e = P·A_e > w_e`，并且**每条属主间边都是宏边**；径向边恰是属主内部的边。

可选的、用来扩大覆盖面的开关（每个都是一项已列出的偏离）：

| 目标 | 开关 | 偏离 |
|---|---|---|
| NameRoutes 采样轮次，从而得到径向的属主间边与 Alg 8 l.12 的删除 | `sample_factor = 1`、`residual_factor = 2`、小的 `p_log2`、`retry_cap` | D-17、D-10 |
| `uint64` 上的未对齐块 | `p_log2 = 0` 且 `b ≤ 3`（`f ≤ 31`） | D-17 |
| 外层原子内部 Alg 4 的径向分支 | 外层 `b ≥ 33`，以及某个原子的局部阈值与局部 `b` 满足 `k ≥ max(32, b_local + 1)`（Q-03） | D-7 |

### 3.5 外层构造开始运行的最小 N

“运行”有几种含义。`N` 是收缩后的规模；外层阈值被禁用；NameRoutes 与局部参数取论文的值。

| 层级 | 含义 | 最小 N | b | 例子 |
|---|---|---|---|---|
| L0 | 进入该分支 | 1 | 任意 | 一个顶点；一个枢轴；仅此而已 |
| L1 | 两个枢轴、一条走廊、一条宏边、一个带边的 NameRoutes 实例 | **2** | 1 | 一条边；两端都是枢轴 |
| L2 | L1 外加一个非枢轴顶点（正的普通状态、非平凡的 `J`） | **4** | 2 | 路径 0–1–2–3 |
| L3 | L2 外加一个具有有限挂接的非单点块（对齐测试得以运行） | **4** | 2 | 路径 0–1–2–3（权重 10, 3, 3）加边 1–3（权重 5） |
| L4 | L3 外加一个非单点原子（Alg 4、截断与替换做实际工作）：**Sec 8.2 的每个阶段都是非平凡的**。[Correction, Q-78] 在此例中唯一的宏边连接两个枢轴，所以每个非单点块都未挂接，对齐测试不会运行（L3 未满足）；SC-8 的输入（加上边 2–5）才提供这一点 | **6** | 3 | 边 0–1 (4)、1–2 (4)、0–3 (100)、3–4 (8)、4–5 (1) |

下界［Implied］：非单点原子位于 `T − A` 的一个分量内，该分量的顶点数少于 `b`，所以 `b ≥ 3`。两个枢轴需要两个被标记的顶点（只有一个标记时没有分支顶点），而每个标记至少消耗它自己的 `b` 个顶点（Sec 5.1 标记规则），所以 `N ≥ 2b`：L4 需要 `N ≥ 6`。当 `b = 1` 时每个顶点都被标记，所以非枢轴顶点需要 `b ≥ 2`，因此 L2 与 L3 需要 `N ≥ 4`。

［Measured］`main` 上的各组件（`owner_decomposition`、`classify_interowner_edges`、`build_aux_graph`、`owner_distances`，以及取坐标 `J`、`K = 16`、经过 l.12 删除之后的 `build_radial_states`）在这些例子上、以及对每个 `(N, b)`（`N ≤ 7`，`b ≤ 4`）各 60,000 个随机连通多重图上运行过：
- L4 例子：枢轴 {0, 3}；一条宏边（0–3）；顶点 4 与 5 处 `r = 8`，`J = 8, 9`；块 `{4, 5}` 位于向内窗口 0 与向外窗口 0 和 1，各自 `s_U = 8`，且权重为 1 的边低于 `s_U/2`；块 `{1, 2}`，`s_U = 4` 且有一条权重为 4 的边（显著，无原子）。[Correction, Q-78] 这六个块全部未挂接：唯一的宏边连接两个枢轴，所以没有任何块具有有限挂接。
- L3 例子：枢轴 {0, 2}；宏边 0–1；顶点 1 与 3 处 `r = J = 3`；三个块 `{1, 3}`，`s_U = 3`，以及权重为 5 的边。
- 随机搜索：`N ≤ 3` 时没有任何实例同时具有两个枢轴与一个非单点块；`N ≤ 5` 时没有任何实例同时具有两个枢轴与一条低于 `s_U` 的窗口内的边；两者都出现在 `N = 4, b = 2` 与 `N = 6, b = 3`。
Sec 8.1 的特殊状态、θ 及其后的一切都没有运行过（不在 `main` 上）。

L4 之外（未测量；仅给出界）：
- 以论文的常数做 NameRoutes 采样，需要多于 `32t` 个不同的枢轴对，因此 `t ≥ 66`（`C(66, 2) = 2145 > 2112`）。取 `b = 1` 时这是 `N = 66`（完全图；没有非枢轴顶点）。与 L4（`b = 3`）合起来：`102 ≤ N ≤ 198`（`t ≤ 2|A0| − 1` 与 `|A0| ≤ N/b` 给出下界；66 个各含 3 个顶点的标记给出上界）。精确的最小值没有计算。
- 取 D-17（`sample 1t`、`residual 2t`）时，采样需要 `t ≥ 6`：与 L4 合起来 `12 ≤ N ≤ 18`。
- 外层原子内部 Alg 4 的径向分支需要至少 32 个顶点的原子，所以外层 `b ≥ 33`，两个枢轴时 `N ≥ 66`。

---

## 4. 现有组件缺少什么

| ID | 缺少 | 需要它的部分 | 细节 |
|---|---|---|---|
| M-1 | Sec 8.1 的输出类型 | 尺度、流水线 | P3-OUTER-STATES 正在进行；假设 A-1（§2.0）列出了字段。特别是每个状态的 `b0_edge`（或 `b0_key` = 字典索引）与块种类 |
| M-2 | `src/paper/` 中的精确 `f` | 尺度 | 一个多字自然数（`experimental/paper_experimental/params.cpp` 里的 `Nat` 是文件内私有的）与 `outer_scale_factor_log2(p_log2, exact h)`；以及精确的 `b = L_N^64`（今天 `outer_separator_size` 会饱和）。之后 Mode B 的 `OuterDerived::f_log2` 应调用同一个函数（Q-04：只有一份实现） |
| M-3 | 外层分支的参数与生成器 | 流水线 | `paper::outer_construction(const ZeroContraction<W>&, Diagnostics&)` 两者都没有。需要：`src/paper/outer/` 中的一个参数结构体（分隔符规则、`pivots::NameRoutesParams`、`LocalHierarchyParams`）与一个 `model::WordRng&`；`PaperOptions` 中的一个种子；把 `PaperParams` 映射到它的 Mode B 驱动程序（P5-EXP-THRESHOLDS） |
| M-4 | 新名称 | 尺度 | 每个字典索引的 `name_weight`（§2.1） |
| M-5 | 以 `EdgeClassification` 表示的外层分类 | 外层状态 | `OuterEdgeClass` → `EdgeClass`，省略 l.12 的删除，且 `macro_key = name`（大概放在 P3-OUTER-STATES 里） |
| M-6 | 可复用的网格 / 支撑 / 同时合并 / 链接 / 子节点排序引擎 | 尺度 | Sec 5.4 的代码位于 `local/assemble.cpp` 的匿名命名空间中（`MergeBuilder`、`floor_index`、`bucket`），并与 `OwnerDecomposition` 的走廊以及 `N^c`、`N^16` 绑定。Mode B 已经有第二份拷贝（`local_assembly.cpp`，Q-59，P5-LOCAL-UNIFY）。要么提取共享引擎，要么明确接受第三份拷贝（由负责人决定） |
| M-7 | 划分层次结构的输出 | 尺度 | `assemble_local_hierarchy` 返回的 `OrderedHierarchy` 对小分量使用精确的 Kruskal 树；Sec 8.2 需要带原子叶子的 `PartitionHierarchy<W>`（`leaf_atom`、`ReplaceableAtom`） |
| M-8 | 输入单位下的分离网格 | 尺度 | `Calendars::graph_pieces[].grid` 以 `2^7` 缩放并到达 `4·max R`；Sec 8.2 需要输入单位下直到 `4K²·R_v` 的点（§2.6）。还需要用于查找的顶点 → 辅助顶点映射（`aux.aux_vertex`） |
| M-9 | 检查器 | 测试 | `check_outer_scales(...)`，对 §2 中每条被违反的不变量返回一条消息，风格同 `check_aux_metric` / `check_calendars` |
| M-10 | 外层原子的局部变体 | Mode B | `substitute_atoms` 与 `owner_distances` 只接受 `LocalHierarchyParams`（规模测试与分隔符）；原子无法选用 D-8 与 D-11 的引擎。只要原子走基本情形，这就无关紧要 |
| M-11 | 诊断信息的接线 | 流水线 | `absorb` 是 `sssp.cpp` 内文件私有的；`evaluate` 返回它自己的 `EvalResult::diag`，必须被合并 |
| M-12 | 对各副本取最小值 | 流水线 | 新状态 → `sub.old_state` → `projection` → 顶点；一次线性遍历，外加 Q-73 的“无发出值”检查 |

原样即可使用的现有部分：`owner_decomposition`、`classify_interowner_edges`、`build_aux_graph`、`owner_distances`、`sat::add` / `sat::mul_pow2`、`pivots::median_prune`、`DisjointSets`、`ForestQueries`、`PartitionHierarchy` / `ReplaceableAtom` / `substitute_atoms` / `redirect_offers`、`validate_ordered_hierarchy`、`evaluate`、`contract_zero_edges`、`restore_predecessor_tree`。

---

## 5. 含糊之处：拟议的开放问题（自 Q-65 起）

这些是给 `open-questions.md` 的提案（该文件由负责人编辑）。均未定案。

| ID | 位置 | 问题 | 候选 | 提议 |
|---|---|---|---|---|
| Q-65 | Sec 8.2 l.2077–2083 | (a) “携带该名称的任一宏边的最小物理权重”是输入权重 `w_e` 还是投影长度 `A_e`，在哪些边上取？(b) 走廊连接使用它自身权重的支撑，还是它的宏名称的支撑？ | (a1) 在宏边上取 `w_e`；(a2) `A_e`。(b1) 自身权重，每条走廊一个键；(b2) 名称 | (a1)：Lemma 7.1 界定的是 `w`。(b1)：正文把“every pivot corridor weight a”与名称并列列出 |
| Q-66 | Sec 8.2 l.2085–2093 | “其他每个非单点正普通块”涵盖哪些块？ | (1) 未对齐块与没有有限挂接的块；(2) 仅未对齐块，未挂接块未作规定 | (1)，如 Sec 5.4 l.1213 对局部情形所述 |
| Q-67 | Sec 8.2 l.2100–2104 | “v 的辅助网格指针”是什么，用什么单位，增加多少个倍频程？ | (1) Sec 7.3 的图片段网格，从 origin 出发加倍换算为输入单位，在顶端延伸 9 次加倍；(2) 同一张网格在 `2^7` 缩放单位下使用（每个块权重都缩放，D-20 的范围限制）；(3) 从 `s_U` 出发的按块网格，`t_U = s_U`，如 Sec 5.4（不是 Sec 8.2 所说的；属于一项偏离） | (1)。倍频程的个数是我们的；论文说的是“a fixed number” |
| Q-68 | Sec 8.2 l.2089–2099；Lemma 8.4 证明 (2) | (a) 间隙测试是严格的吗（`next.lower < 2·upper`）？(b) 类 `[a, 2a)` 的节点高度与链接界。(c) 处理分离类与全局类的顺序 | (a) 严格 / 非严格。(b) 高度 `a`，链接 `< 2a`，`c0 = 2`（如 Lemma 5.10）/ 高度 `2a` 且链接 `< 2a`。(c) 先分离网格 / 按值交错 | (a) 严格（“smaller than a factor two”）。(b) 高度 `a`。(c) 先分离；与任何和各高度一致的顺序等价 |
| Q-69 | CLAUDE.md 的回退策略；Alg 7 l.6；Lemma 5.15 | 当一个顶点被多个局部基本情形覆盖时（Alg 7 的原子，以及通过它的副本被至多四个 Sec 8.2 原子覆盖），回退占比如何定义？ | (1) 按原因的原始总和，注明“可能超过 N”；(2) 不同的物理顶点（需要额外记账）；(3) 基本情形原子中的状态数除以全部状态数，按阶段报告 | (3) 按阶段报告并附上 (1) 的原始值；需要负责人同意 |
| Q-70 | Alg 8 l.9；`include/sssp/paper.hpp` | 外层分支的生成器与参数如何进入 Mode A 的 API？ | (1) `PaperOptions::seed`，带固定默认值与内部的 `WordRng`；(2) 一个接受 `WordRng&` 的重载；(3) 没有公开途径，测试直接调用 `outer::outer_construction` | (1) 与 (3) 并用；这会改动公开头文件，所以由负责人决定 |
| Q-71 | Sec 8.2 l.2080–2081；Sec 7.1 l.1644 | 在 Mode A 中 `b = L_N^64` 是以饱和形式保存的。`F` 由哪个 `h` 推出？ | (1) 精确的多字 `h`（M-2）；(2) 当饱和后的 `b ≥ N` 时跳过 `F`（此时 `t = 1`，`F` 不被使用） | (1)；(2) 增加了论文里没有的分支 |
| Q-72 | Alg 8 l.8–13，`t = 1` | 单枢轴的运行（无宏边，`J` 即 Alg 7 给出的完整距离）算回退吗？ | (1) 否，它是论文的构造；只做计数器；(2) 记录一条 FallbackRecord | (1) |
| Q-73 | Alg 8 l.15 | “每个顶点的最小发出值”是在哪些状态上取的？若某个顶点没有怎么办？ | (1) 所有投影到该顶点的状态（普通、特殊、局部子副本）；没有则 → `std::logic_error`；(2) 仅普通副本 | (1)：特殊状态投影为物理游走，所以它们不会低估 |
| Q-74 | Sec 8.1 l.2000–2001；Def 5.11–5.12 | (a) 哪些边构成 Sec 8.2 原子的内部图？(b) 只有一个状态的分量是原子吗？ | (a1) 两端都在该原子内的每条窗口内的边（包括较重的弦，如 Q-53）；(a2) 仅 `< t_U` 的边。(b1) 单点部分，不调用 Alg 4；(b2) 单状态原子（如 Alg 7 目前所做） | (a1)，`substitute_atoms` 要求；(b1)，依 Def 5.12 |

已检查并发现一致的地方（无问题）：Lemma 8.4 证明 (1) 的对齐界（`a_U/F² ≤ s_U ≤ F·a_U`）由 Lemma 8.2 与 Lemma 7.1 得出，并余下 `8K²h` 的富余因子；两种网格下 `s_U/2 < t_U ≤ s_U` 都成立；原子截断值 `t_U` 满足 Def 5.11（生成树的边 `< t_U < 2t_U`，父节点高度 `≥ t_U`）。

---

## 6. 拟议的验收测试

权重：除另有说明外为 `uint64`；`double` 只用二进制权重（Q-55 选项 A）。种子固定，并在 NameRoutes 会采样的地方变化。预言机相互独立：对原始图用 Bellman–Ford，对状态弧用 Tarjan SCC，对 `f` 用 Python 大整数。

### 6.1 P3-OUTER-SCALES

拟议的文件：`tests/unit/paper/outer/test_outer_scales.cpp`、`tests/property/paper/test_outer_scales_property.cpp`。

| 测试 | 输入 | 检查 |
|---|---|---|
| SC-1 `scale_factor_log2 matches big-integer values` | 来自 §1.4 的 `(p, b)`，包括 `b = 589^64`（1808）、`20^64`（871）、`b = 3, p = 17`（48）、`b = 160000, p = 1`（77） | `f` 等于表中的值；在测试中用第二个多字例程重新计算 `2^f > X ≥ 2^(f−1)` |
| SC-2 `scale factor is not derived from a saturated h` | `N = 10^6`，论文的 `b` | 871，而不是 233 |
| SC-3 `F^2 s >= a is decided exactly` | 随机的 `s`、`a`、`f`，在 `uint64` 与二进制 `double` 中，围绕饱和边界 | 等于精确答案（`unsigned __int128` 移位，或对 2 的幂用指数运算） |
| SC-4 `new names are the least macro weight` | 手工分类，两个名称、四条宏边 | `name_weight`；`a ≤ w < 2P(2h+1)a`；每个有限连接上 `a ≤ B0 ≤ θ < F·a` |
| SC-5 `hand example L3` | §3.5 L3，`b = 2` | 三个非单点块 `{1, 3}`，`s_U = 3`，对齐，`t_U ∈ (1.5, 3]`，无原子，权重为 5 的边被分到括住它的全局类。[Correction, Q-78] 三个块中只有两个是对齐的（向内窗口 0 与向外窗口 0，`t_U = 3`）；向外窗口 1 没有进入弧，所以 I0 为无穷，该块未挂接（分离网格，`t_U = 3`） |
| SC-6 `hand example L4` | §3.5 L4，`b = 3` | 原子 `{4, 5}` 位于三个块中，截断值 `t_U ∈ (4, 8]`，`H` = 一条权重为 1 的边；块 `{1, 2}` 无原子；`ph` 被 `substitute_atoms` 接受；`validate_ordered_hierarchy` 在 `c0 ≤ 2` 下无问题。[Correction, Q-78] 唯一的宏边连接两个枢轴，所以六个非单点块全部未挂接（分离网格；`t_U = 8` 的原子 `{4, 5}` 仍然存在） |
| SC-7 `unaligned block on uint64` | L3 图，走廊权重 `2^60`，`p_log2 = 0`，`b = 2`（`f = 29`） | 块 `{1, 3}` 未对齐；它所有窗口内的边 `< a_U/2`；块上 θ 相同；`t_U` 在片段网格上；`2·p_top ≤ g` |
| SC-8 `unaligned block with an atom` | L4 图，权重 0–3：`2^63 + 1`，3–4：2，4–5：1，`p_log2 = 0`，`b = 3`（`f = 31`）。[Correction, Q-78] 用这个输入时块 `{4, 5}` 是未挂接的，而不是未对齐的；测试加入非树边 2–5（权重 `2^63 + 2`），使这些块有一个宏端点，名称手工给定 | 块 `{4, 5}` 未对齐（`F²·2 = 2^63 < a_U`）；每个是一个截断值为 2 的原子；块 `{1, 2}` 对齐。[Correction, Q-78] 加入边 2–5 之后：`{4, 5}` 的向内块与向外窗口 0 未对齐，向外窗口 1 未挂接，各自是一个截断值为 2 的原子；块 `{1, 2}` 中两个已挂接的是对齐的，向外窗口 1 未挂接 |
| SC-9 `several grid components (double)` | 二进制权重，两个走廊权重相距超过 `3f` 个倍频程 | 两个分量；各分量内部 `grid[i+1] = 2·grid[i]`；没有类跨越分量的起点；间隙小于二倍的支撑被合并 |
| SC-10 `unattached blocks` | 一个枢轴（`b > N`） | 每个非单点块都未挂接，只有分离网格，若干个根按凝聚（condensation）顺序排列 |
| SC-11 性质：Lemma 8.4 | 随机连通图，`N ≤ 200`，`b ∈ {1, 2, 3, 4, 8, ⌈log2 N⌉}`，测试框架的生成器族，多个种子 | 通过 `check_outer_scales` 检查 §2 的每条不变量：未对齐的 MST 边 `< a_U/2`；对齐的显著权重位于其支撑之内；`s_U` = 最小的原始半径；`s_U/2 < t_U ≤ s_U`；原子 `< b` 个顶点，一个属主，一个 `r`；层次结构在 `c0 ≤ 2` 下通过 Def 4.1 的验证 |
| SC-12 性质：SCC 相等 | 小 `N`（≤ 40） | 对每个内部节点：其下的状态 = 长度 `<` 类上界的弧的一个 SCC；子节点 = 长度 `<` 高度的弧的 SCC（Tarjan） |
| SC-13 `alignment by median pruning equals the per-block test` | 随机分组 | 相同的分类；每组的谓词测试次数 `≤ ⌊log2 n⌋ + 1` |
| SC-14 `distances on the hierarchy` | SC-11 的输入 | 在替换后的层次结构上以 Alg 8 l.14 的种子运行 `evaluate`：逐顶点的最小值等于 Bellman–Ford |

### 6.2 P3-OUTER-PIPELINE

拟议的文件：`tests/unit/paper/outer/test_shortest_paths.cpp`、`tests/differential/paper_experimental/test_outer_forced_vs_dijkstra.cpp`（该名称已在 `PAPER_MAPPING.md` 中）。

| 测试 | 输入 | 检查 |
|---|---|---|
| PL-1 `Mode A never takes the outer branch` | 差分测试的各图族 | `takes_outer_construction` 为假；恰有一条基本情形 FallbackRecord（保持现有行为） |
| PL-2 `outer_construction with the paper's constants` | 在小图上直接调用（`N ≤ 300`） | 距离等于 Bellman–Ford；`outer.pivots = 1`；无宏边；没有带 `outer:` 原因的 FallbackRecord |
| PL-3 `forced outer, levels L0–L4` | §3.5 的五个例子 | 距离与树等于预言机；`verify_sssp`；该层级的特征出现在计数器中（枢轴、宏边、非单点块、原子） |
| PL-4 `forced outer vs oracle` | 每个生成器族，`N ≤ 2000`，`b ∈ {1, 2, 3, 4, 8, 16, 64}`，种子固定并变化 | 距离等于 Bellman–Ford 与基线 Dijkstra；`verify_sssp`（根为 0，每条弧可行，树弧紧，父节点到达 `s`，不可达 = 无穷） |
| PL-5 `stage checkers pass on the intermediate states` | PL-4 的输入 | `check_route_names`、`check_aux_metric`、`check_calendars`、`check_atoms`、P3-OUTER-STATES 的检查器、`check_outer_scales`、`validate_seeded_ordered_hierarchy` 全部为空 |
| PL-6 `fallback records` | PL-4 的输入 | 没有外层基本情形记录；每条记录都是局部基本情形的原因（设置了上限时则为 D-10 的原因）；每个进入基本情形的 Alg 4 调用一条记录；参数集被记录（`PaperParams::record`） |
| PL-7 `zero edges, unreachable vertices, loops, parallel edges` | 含零权环与多个分量的输入 | 还原正确；收缩后的 `N` 就是外层分支所看到的 |
| PL-8 `seed independence` | 16 个种子，输入使 NameRoutes 采样（D-17 参数，`t ≥ 6`） | 每个种子的距离都相同；至少有一个实例满足 `aux_metric.radial_interowner_edges > 0` 且 `outer.radial_edges_deleted > 0` |
| PL-9 `source position` | 源在枢轴处、在不处于任何原子内的非枢轴处、在原子内 | 全都等于预言机（Alg 8 l.14 的重定向） |
| PL-10 `double weights` | 二进制 `double`；非二进制 `double` | 二进制：与精确参考完全相等；非二进制：要么在记录的容差内相等，要么抛出异常，绝不静默地给出错误结果（Mode B 的失败规则，Q-63） |
| PL-11 `unaligned end to end` | SC-7 与 SC-8 的输入 | 距离等于预言机 |
| PL-12 `sanitizer` | 缩小规模的 PL-3、PL-4 | 在 `make sanitizer` 下干净（P4-SANITIZER） |

对 `docs/TASKS.md` 中验收标准的建议补充（供负责人参考）：
- P3-OUTER-SCALES：“`f` 等于笔记 04 §1.4 的大整数表，且从不由饱和的 `h` 推出；未对齐分支由一个 `p_log2 = 0` 的 `uint64` 测试与一个二进制 `double` 测试覆盖；`check_outer_scales` 在每个性质测试输入上为空。”
- P3-OUTER-PIPELINE：“覆盖笔记 04 §3.5 的 L0–L4 层级；以论文的常数直接调用（单枢轴）等于预言机；除局部基本情形外没有其他 FallbackRecord；没有发出值的顶点会抛出异常。”
