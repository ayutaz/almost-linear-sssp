# 数值设计:精确模型与 IEEE-754 double

[English](numerics.md) | [日本語](numerics.ja.md) | **简体中文**

论文(arXiv 2611.34911v1,Sec 2)在**比较-加法模型**中讨论:边长是非负实数,算法只能对它们做精确的比较和加法。CPU 上的 `double` 每次加法都会舍入。本文档规定本项目如何区分两者、对 `double` 结果使用什么容差进行检查及其理由,以及为什么对整数权重做精确检查。

代码:`include/sssp/verify.hpp`(`Tolerance`、`verify_sssp`、`write_result_json`)、
`include/sssp/counted.hpp`(`Counted<W>`)、`tests/python/reference/`(精确有理数 oracle)。

## 1. 两种权重类型,两种正确性标准

| 权重类型        | 算术                            | 正确性标准                                    |
|-----------------|---------------------------------|-----------------------------------------------|
| `std::uint64_t` | 精确(与模型一致)              | **精确**相等,溢出视为错误                    |
| `double`        | IEEE-754 binary64,就近舍入      | 在 Sec 3 的容差之内,或在舍入代数中精确相等(Sec 2) |

整数是正确性的主要目标:对整数输入,只要没有任何和达到 `WeightTraits<uint64_t>::infinity() = 2^64 - 1`,模型的每个运算都是精确的,因此论文的实现必须复现精确的距离。`double` 的结果则用舍入误差界来判定。

## 2. 松弛类算法在 `double` 中计算的是什么

用 `fl(x)` 表示舍入到最近的 double,记 `a ⊕ w = fl(a + w)`。对 `w ≥ 0`:

1. **单调:** `a ≤ b ⇒ a ⊕ w ≤ b ⊕ w`(舍入是单调的);
2. **不减:** `a ⊕ w ≥ a`,但**不是严格的**——存在吸收:`2^60 ⊕ 1 = 2^60`。

(1)与(2)使 `(double, min, ⊕)` 成为单调、不减的路径代数,Dijkstra 和 Bellman–Ford 在其上都是正确的。因此两者计算出完全相同的值

    d̂(v) = 所有 s–v 路径 P 上,P 从左到右舍入求和的最小值,

并且该最小值可由简单路径取到(由(1)与(2),去掉一个环不会使左折叠变大)。测试所用的推论如下:

* 任意两个基于松弛的算法(Bellman–Ford、各种 Dijkstra 变体)在 `dist` 上必须**逐位**一致;只有父节点树可以不同(平局)。
* 对这类结果,`verify_sssp` 的局部检查**精确**成立:每条弧上 `dist[v] == fl(dist[parent] + w)` 且 `dist[v] ≤ fl(dist[u] + w)`。因此即使对 `double`,默认也是 `Tolerance::exact()`。
* 使用严格的 `<` 松弛时,性质(2)保证即使存在吸收,父节点图也保持无环。

以其他结合顺序做加法的算法(例如预先计算的捷径长度 `d(u) + (w1 + w2)`、Floyd–Warshall 的 `d[i][k] + d[k][j]`)得到的舍入值不同,因此可能与松弛的结果相差一个舍入量。容差正是为此而设。

## 3. 容差

`Tolerance{rel, abs}` 在下式成立时接受两个距离值 `a`、`b`:

    |a − b| ≤ rel · max(|a|, |b|) + abs                  (单侧:a ≤ b + 同样的余量)

`Tolerance::for_path_hops(k)` 取 `rel = γ_k`、`abs = 0`,其中

    u = 2^-53(binary64 的单位舍入),     γ_k = k·u / (1 − k·u)     (要求 k·u < 1/2)。

`Tolerance::for_graph(g)` 使用 `k = n − 1`。

### 3.1 为什么是 γ_k

对 `m` 个非负 double `x_1..x_m` 以**任意**加括号方式(任意由 `m − 1` 次加法构成的二叉树)求和,有(Higham,*Accuracy and Stability of Numerical Algorithms*,第 2 版,Lemma 3.1 与式 (4.4)):

    |ŝ − s| ≤ γ_{m−1} · Σ|x_i| = γ_{m−1} · s.

最短路径是简单路径,因此至多有 `n − 1` 条边;由于 `dist[source] = 0`,第一次加法 `0 + w_1` 是精确的,所以至多有 `n − 2` 次加法发生舍入。因此 `γ_{n−1}` 对任意计算顺序都是安全的(保守一步的)界。设 `d*` 为精确距离,`d̂ = ŝ(P)` 为计算得到的距离,由某条路径 `P`(树上的路径)实现:

* 下界:`d̂ = ŝ(P) ≥ (1 − γ) L(P) ≥ (1 − γ) d*`;
* 上界:`d̂ ≤ ŝ(P*) ≤ (1 + γ) L(P*) = (1 + γ) d*`(算法对路径取最小);

于是 **`|d̂ − d*| ≤ γ_{n−1} · d*`**——这就是 `check_double.py` 针对精确有理数 `d*` 所检查的界。

### 3.2 比较两个舍入结果

若 `a` 与 `b` 都是舍入结果(例如论文的实现与 double 版 Floyd–Warshall),则各自都在 `γ_k d*` 之内,所以 `|a − b| ≤ 2γ_k d* ≤ 2γ_k/(1 − γ_k) · max(a, b) = γ_{2k} · max(a, b)`(恒等式 `2γ_k / (1 − γ_k) = γ_{2k}` 是精确的)。使用 `Tolerance::for_path_hops(2(n − 1))`。与最接近 `d*` 的 double 比较时,再加上半个 ulp:`γ_k + u ≤ γ_{k+1}`。

### 3.3 为什么默认没有绝对项

IEEE-754 加法不会因下溢而损失精度:两个 double 之和若为非规格化数,则该和是精确的。因此对非负的和,误差纯粹是相对误差,`abs = 0`。`abs` 是为那些有意离开比较-加法模型的调用方(例如用减法 `d − base` 算出的桶键),或需要显式余量的调用方而设;使用处必须说明理由。

### 3.4 在容差下,局部检查能得出和不能得出的结论

树检查 `dist[v] ≈ fl(dist[p] + w)` 与证书 `dist[v] ≲ fl(dist[u] + w)` 对每条弧使用相同的余量。有了余量,证书的论证(`dist[v] ≤ d*(v)`)每一跳累积一份余量,因此仅靠带容差的局部检查,只能证明 `dist ≤ d*(1 + γ)^{hops}` 这一类的界。非松弛类算法的全局精度要通过与参考结果比较来确立:参考是精确有理数(`check_double.py`),或是使用 `for_path_hops(2(n − 1))` 的松弛结果。

`check_double.py` 还会限制所报告的树路径的精确长度 `L(v)`:由于 `dist(v)` 是该路径的舍入和,有 `(1 − rel) L ≤ dist(v) ≤ (1 + rel) d* + abs`,因此 `d* ≤ L ≤ ((1 + rel) d* + abs) / (1 − rel)`。

## 4. 模型与 double 真正不同的地方

对实数有效、对 double 却不成立的推理,不得在 `double` 代码中依赖:

* **严格增长:** `w > 0 ⇒ a + w > a`(存在吸收时为假;见 `tests/support/graphs.hpp` 中的 `extreme_ratio`,权重为 `2^60` 和 `1`);
* **结合律:** `(a + b) + c = a + (b + c)`(为假;使结果依赖于顺序);
* **抵消:** `(a + w) − w = a`(模型中没有减法;`Counted<W>` 禁止减法)。

终止性或桶的论证用到严格增长的算法,必须用 `double` 在 `extreme_ratio` 类型的输入上测试,并主要在 `uint64_t` 上验证。

## 5. 为什么禁止 fast-math

`cmake/Warnings.cmake` 向每个目标传入 `-fno-fast-math`。`-ffast-math` 会启用:

* `-ffinite-math-only`:编译器可以假定不存在无穷大,于是 `dist[v] == inf`(我们的“不可达”哨兵)可能被折叠为 `false`——结果悄无声息地出错;
* `-fassociative-math` / 重新结合:破坏 Sec 2 的舍入路径代数,使松弛类算法不再逐位一致,结果依赖优化级别;
* flush-to-zero / denormals-are-zero:非规格化数的和不再精确,Sec 3.3 不再成立;
* `-fno-signed-zeros`、倒数运算:此处无害,但也没有好处。

浮点收缩(FMA)无关紧要,因为距离只使用加法。假定 `FLT_EVAL_METHOD == 0`(SSE2 / AArch64,没有 x87 扩展精度)。

## 6. 整数权重:精确,并检测溢出

* `infinity() = 2^64 − 1` 为保留值;每个有限距离都 `< 2^64 − 1`。
* 若 `a + b ≥ infinity()`,`WeightTraits<uint64_t>::checked_add` 抛出 `std::overflow_error`。
* `bellman_ford` 跳过会溢出的候选(它可能位于非最短路径上),仅当收敛之后,与某个有限顶点相邻的顶点仍为无穷大时才抛出 `std::overflow_error`——即当且仅当某个真实距离无法表示。
* `double` 遵循同样的规则,只需把“达到无穷大”读作“舍入后的和为 `+inf`”:若两个有限距离的舍入和不是有限值(例如 `1e308 + 1e308`),`WeightTraits<double>::checked_add` 抛出 `std::overflow_error`,而不是对一个可达顶点返回表示“不可达”的 `+inf` 哨兵;当且仅当 Sec 2 的某个舍入距离 `d̂(v)` 为 `+inf` 时,`bellman_ford` 抛出异常(详见:`docs/design/baselines.md`)。
* `verify_sssp` 对整数做精确比较,忽略 `Tolerance`;溢出的 `dist[parent] + w` 报告为 `arithmetic_overflow`;在证书上,溢出的界是空洞的(真正的约束位于可表示范围之外)。

## 7. 对 C++ `double` 结果做精确验证

1. 把图写成边列表 `n m` / `u v w`,权重以 `%.17g` 输出(使 `Fraction(float(w))` 能恢复 C++ 代码所用的精确 double)。
2. 用 `sssp::write_result_json(os, source, result)` 导出结果:
   `{"n":..,"source":..,"dist":[%.17g or null],"parent":[-1..],"parent_edge":[-1..]}`。
3. `uv run python tests/python/reference/check_double.py graph.txt result.json`
   (`--hops K`、`--rel R`、`--abs A`、`--weights exact` 表示把十进制权重当作实数)。
   它用 `fractions.Fraction` 重新计算 `d*`,把导出的每个数读作其 double 的精确二进制值,并检查 Sec 3.1 与树(退出码 0 / 1 / 2 = 接受 / 拒绝 / 无法读取)。

导出格式在两侧都被固定:`tests/unit/verify/test_verify.cpp` 断言精确的字符串,`tests/python/test_exact_reference.py::test_cpp_dump_format_is_readable_and_exact` 解析同一个字符串。

## 8. 统计模型运算

`sssp::Counted<W>`(W = `double` 或 `uint64_t`)用线程局部计数器(`counted_stats()`、`reset_counted_stats()`、`CountedScope`)统计 `+`/`+=` 和六种比较。它没有 `-`、`*`、`/`,也不能转换回 `W`,因此能用 `Counted<W>` 编译通过的代码,可证明对权重只使用比较-加法运算。`checked_add` 中的溢出防护不计入(属于模型之外的簿记)。局限:`validate_edges` 只校验有符号或浮点的 `W`,因此负数或 NaN 的 `Counted<double>` 权重在建图时不会被拒绝;应先用经过校验的原始权重建图。
