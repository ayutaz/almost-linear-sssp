# Dijkstra 基线

[English](baselines.md) | [日本語](baselines.ja.md) | **简体中文**

基线是针对论文算法的一切性能主张的参照点,因此编写时追求的是*公平*而不是巧妙:只有一个共享的扫描循环,堆在本仓库内实现(不使用 `std::priority_queue`),所有堆的存储空间都在第一次 push 之前分配好,并且各处使用相同的操作计数器。

公开 API:`include/sssp/baseline.hpp`。源码:`src/baseline/`。

| 函数 | 堆(`src/baseline/heaps/`) | Decrease-key | 权重类型 | 时间 | 堆的空间 |
|---|---|---|---|---|---|
| `dijkstra_binary_heap` | `BinaryHeap<K, V>`(延迟删除) | 无(重新 push) | `uint64_t`、`double` | O((n + m) log m) | O(m) 个条目,预留 `num_arcs + 1` 个 |
| `dijkstra_4ary_heap` | `IndexedDaryHeap<K, 4>` | 有 | `uint64_t`、`double` | O((n + m) log₄ n) | O(n) |
| `dijkstra_fibonacci_heap` | `FibonacciHeap<K>` | 有(cut + cascading cut) | `uint64_t`、`double` | O(m + n log n) 摊还 | O(n) 的节点池 |
| `dijkstra_radix_heap` | `RadixHeap`(65 个桶) | 有(桶移动) | 仅 `uint64_t` | O(m + 64 n) | O(n) 的节点池 + 65 个桶头 |

`n` = 顶点数,`m` = 无向边数;CSR 的 `num_arcs` 为 2m 减去自环对应的弧数。
模板在各个 `.cpp` 文件中针对 `std::uint64_t` 和 `double` 显式实例化。

## 共享的扫描循环

`src/baseline/dijkstra_common.hpp`:

1. 如果 `source >= n`(空图也包含在内),`init_result` 抛出 `std::out_of_range`;否则设置 `dist = WeightTraits<W>::infinity()`、`parent = kNoVertex`、`parent_edge = kNoEdge`、`dist[source] = 0`。
2. 弹出最小的 `(d_u, u)`,扫描 CSR 中 `u -> v` 的**每一条**弧,包括指向已确定顶点的弧和自环。不查询"已确定"位图:在权重非负时,这样的弧不可能改善 `dist[v]`,而跳过它们会使各实现之间的计数器不一致。
3. 计算 `nd = WeightTraits<W>::checked_add(d_u, w)`;如果 `nd < dist[v]`(严格不等),就设置 `dist[v] = nd`、`parent[v] = u`、`parent_edge[v] = arc.edge`,并执行 push 或 decrease-key。

`dijkstra_binary_heap` 对延迟堆使用同一个循环:被改善的顶点会再次 push,弹出的条目如果 `dist[u] < d`,就是过期条目并被跳过。因为只有在严格改善时才会 push,所以顶点的最终距离恰好被 push 一次,从而每个顶点最多被扫描一次。

## 并列时的处理

堆只比较键(`a < b`),没有第二键。键相等时的弹出顺序由堆的结构决定,对于给定的图和源点是确定的(CSR 的弧序就是边列表的顺序),但不同的堆之间会不同。因此:

- 对于整数权重,`dist` 在所有基线中完全相同;
- 当存在多条最短路时,`parent` / `parent_edge` 在不同基线之间可能不同。每一个都是合法的最短路树;测试从结构上检查树(`parent_edge` 存在,连接 `parent[v]` 与 `v`,`dist[v] == dist[parent] + w`,并且沿父节点链能到达源点),而不是比较父节点。
- 先到先得:之后才发现的同样长度的路径不会替换父节点。

## 零权重

允许零权重的边和零权重的环。严格改善(`nd < dist[v]`)保证了零权重的环不会产生父节点的环:顶点的父节点总是在该顶点最后一次被改善之前就已弹出(已确定)的顶点。radix 堆接受等于上一次弹出键的键(桶 0),而这正是零权重的弧所产生的键。

## 整数权重与溢出

对于 `std::uint64_t`,当 `d_u + w >= 2^64 - 1`(和将达到 `infinity()`)时,`WeightTraits::checked_add` 抛出 `std::overflow_error`。这个检查在扫描的**每一条**弧上都会执行,所以只要有任何候选和溢出,运行就会抛出异常,即使该候选并不会改善 `dist[v]`(例如指回源点的一条巨大的边)。这是有意为之的:结果绝不会悄悄出错。因此,权重为 `2^64 - 1` 的边,只要它的某个端点被扫描,就一定会抛出异常。

哪些输入会抛出异常(`double` 也适用同样的规则,见下文):

- **每一个 Dijkstra 基线**都在完全相同的输入上抛出异常:当且仅当存在某个可达的 `u`,其距离 `d(u)` 可表示,且它有一条权重为 `w` 的弧使得 `d(u) + w` 不可表示。该条件不依赖于堆,因为顶点只会以其最终距离被扫描,而每个距离可表示的可达顶点都会被扫描,除非运行已经抛出了异常。
- **`bellman_ford`**(独立的预言机,不调用 `checked_add`)会跳过溢出的候选,当且仅当某个可达顶点的距离不可表示时才抛出异常。所以只要 `bellman_ford` 抛出异常,每个 Dijkstra 也都会抛出,但反过来不成立:对于单条边 `{0, 1, 2^64 - 2}`,`bellman_ford` 返回 `{0, 2^64 - 2}`,而每个 Dijkstra 在扫描顶点 1(指回 0 的弧)时都会抛出异常。

这两条规则由 `tests/differential/test_baselines_differential.cpp`(uint64,带独立的饱和运算预言机)和 `tests/differential/test_double_overflow_baselines.cpp`(double,带独立的舍入预言机)固定下来。

## 浮点权重

对于 `double`,`checked_add` 是两个有限距离的 IEEE-754 加法(舍入到最近;从不启用 `-ffast-math`),如果舍入后的和不是有限值,就抛出 `std::overflow_error`。`+inf` 是"不可达"的标记,所以舍入为 `+inf` 的报价(例如 `1e308 + 1e308`,或 `DBL_MAX + 2^970`,即 `DBL_MAX` 的半个 ulp)若不处理,就会被严格的 `nd < dist[v]` 丢弃,可达顶点会被报告为不可达(回归测试 `tests/regression/cases/20261127-double-sum-overflow.json`)。而被吸收的和(`DBL_MAX + 2^969 = DBL_MAX`)是有限的,不会抛出异常。抛出异常的规则与上一节相同,只是把"可表示"理解为"舍入后有限";`bellman_ford` 用它自己的 `!= +inf` 检查来发现这种情况。所依赖的性质如下:

- 当 `w >= 0` 时 `a + w >= a`,所以指向已确定顶点的弧不会改善它们,堆的键按非递减顺序弹出;
- `dist[v]` 存储的恰好是被比较过的值 `dist[parent[v]] + w`,所以树的检查 `dist[v] == dist[parent] + w` *精确地*成立;
- 距离是沿路径按路径顺序求的和。不同的算法(例如把子路径的和相加的 Floyd–Warshall)舍入可能不同,所以除非权重是二进有理数,否则跨算法比较 `double` 的结果需要容差(测试中精确比较使用 1/16 的倍数,一般实数使用 1e-12 的相对容差)。

radix 堆仅支持整数,不会针对 `double` 实例化。

## 诊断(`SsspResult::diag.ops`)

所有计数器都通过 `SSSP_COUNT` 递增,当 `SSSP_DIAGNOSTICS=0`(release / 计时构建)时编译为空,此时它们保持为 0。

| 计数器 | 定义 |
|---|---|
| `vertex_scans` | 扫描了邻接表的顶点(= 可达顶点) |
| `edge_scans` | 检查过的弧(被扫描顶点的 CSR 度数之和;无向边对每个被扫描的端点各计一次,自环计一次) |
| `relaxations` | 对 `dist[v]` 成功的严格改善 |
| `additions` | `d_u + w` 的求值(= `edge_scans`) |
| `comparisons` | `nd < dist[v]`(每条弧一次)+ 每次 pop 一次过期检查(仅二叉堆)+ 堆内部进行的键比较(`heap.comparisons()`) |
| `queue_pushes` | 堆的插入(二叉:1 + `relaxations`;带索引:1 + 首次到达的顶点数) |
| `queue_pops` | extract-min 的调用,包括过期的(二叉:= `queue_pushes`;带索引:= `vertex_scans`) |
| `decrease_keys` | decrease-key 的调用(二叉:始终为 0;带索引:`relaxations + 1 - queue_pushes`) |

堆内部的键比较:

- `BinaryHeap`、`IndexedDaryHeap`:sift-up / sift-down 中的每一次 `key < key`(d 叉堆的 sift-down:找出最小子节点最多 D-1 次,再与被移动的元素比较 1 次)。
- `FibonacciHeap`:插入时(第一个根除外)和 decrease-key 时与最小值的比较,decrease-key 时与父节点的比较,consolidation 期间每次 link 一次,以及在 consolidation 之后的根中查找最小值。
- `RadixHeap`:只有在被重新分配的桶中查找最小值(k 个元素需要 k − 1 次比较)。桶的下标用 `bit_width(key XOR last)` 计算,这不算作比较。因此 radix 堆的比较次数不能与基于比较的堆直接比较;请比较 `edge_scans` / 时间。

堆内部的参数检查(id 的范围、成员关系、键增大、radix 的单调性)不计数。

手算示例(`tests/unit/baseline/dijkstra_test_support.hpp` 中的经典 6 顶点图,源点 0):顶点扫描 6 次,弧扫描和加法 18 次,relaxation 7 次;二叉堆 push 8 次 / pop 8 次(其中过期 2 次);带索引的堆 push 6 次,decrease-key 2 次,pop 6 次。

## 分配策略

运行期间,任何堆操作都不分配内存:

- `BinaryHeap` 在第一次 push 之前预留 `num_arcs + 1` 个条目(push 次数的上界)。
- `IndexedDaryHeap` 在构造函数中分配其堆数组、位置映射和键数组(大小为 n)。
- `FibonacciHeap` 在构造函数中为每个顶点分配一个节点,并分配一个 n 项的根暂存缓冲区;度数表是固定的 `std::array<id, 64>`(当 n < 2³² 时,最大度数 < 1.4405 log₂ n + 2 < 64)。
- `RadixHeap` 在构造函数中为每个顶点分配一个节点(侵入式双向链接的桶链表);桶头是固定的 `std::array<id, 65>`。

结果向量(`dist`、`parent`、`parent_edge`)由 `init_result` 一次性分配。

## 堆的误用

对空堆执行 pop 或使用 `[0, n)` 之外的 id 时,所有堆都抛出 `std::out_of_range`;push 已存在的 id、对不存在的 id 做 decrease-key 或把键改为更大的值,以及(radix 中)键小于上一次弹出的键时,抛出 `std::invalid_argument`。Dijkstra 不会触发这些情况。

## 测试

- 堆:`tests/unit/baseline/heaps/`(手算的弹出顺序、cascading cut、consolidation、桶的行为、计数器)和 `tests/property/baseline/`(随机操作序列与排好序的 `std::multiset` / `std::map` 参照对比,种子 {1, 7, 42, 2026, 0xC0FFEE},操作后运行不变量校验器)。
- Dijkstra:`tests/unit/baseline/test_dijkstra_*.cpp`,共享的检查在 `dijkstra_test_support.hpp` 中:手算的图(经典图、零权重 / 重边 / 自环 / 不连通 / 孤立顶点、零权重的环、单顶点、二进有理数的 double)、无效的源点、溢出、计数器,以及在随机多重图上使用仅供测试的 Floyd–Warshall 进行的检查(种子 {1, 2, 3, 42, 2026},从每个源点出发的 125 个 n ≤ 30 的图,外加每个种子一个 n = 200 的图)。
