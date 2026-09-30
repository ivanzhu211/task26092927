# freq_bounds — 跨媒体曝光频次区间精确求界

纯 Python 标准库实现。输入只有**汇总数据**（总体人数、四个渠道各自触达人数、
六个渠道两两重叠人数），不含任何逐人记录。输出去重总触达人数，以及曝光频次
（一个人覆盖的渠道数；同一渠道不重复计次）为 1/2/3/4 的人数各自的真实最小值
与最大值。

> 两两重叠不能唯一确定频次直方图。本库对每个指标分别求其在所有与输入完全一致
> 的**整数人群**上的真实极值；不要求同一个人群同时达到所有指标的端点。

## 安装 / 导入

无需第三方依赖，直接把仓库目录放到 `PYTHONPATH` 即可：

```python
from freq_bounds import Summary, bounds, InfeasibleError
```

## 输入约定

`Summary.from_mapping({...})` 接受一个字典：

- `n`：总体人数（非负整数），**包含完全未曝光者**。
- `reaches`：四个渠道 `TV / Digital / Social / Print` 各自触达人数。
- `overlaps`：六个渠道对的重叠人数（同时覆盖两个渠道的人，**也包含**覆盖了
  第三、第四个渠道的人）。键可以是有序元组 `("TV","Digital")`，也可以是
  `"TV|Digital"` 形式的字符串；六个对缺一不可，顺序无关。

所有数字必须是非负整数。基础局部一致性（触达不超过 `n`、重叠不超过任一渠道
触达、键合法）在构造时立即校验。

## 输出约定

`bounds(summary)` 返回 `BoundsResult`：

- `.reach`：`MetricBounds(minimum, maximum)`，去重总触达人数的精确区间。
- `.frequencies[k]`：恰好覆盖 `k` 个渠道（频次 k，`k ∈ {1,2,3,4}`）人数的
  `MetricBounds(minimum, maximum)`。

每个端点都是某个与全部输入完全一致的整数人群真实达到的值（精确极值，不是
保守界）。若局部看似合理的汇总实际上无法由任何整数人群同时满足，抛出
`InfeasibleError`。

## 使用示例

```python
from freq_bounds import Summary, bounds, InfeasibleError

summary = Summary.from_mapping({
    "n": 4,
    "reaches": {"TV": 2, "Digital": 2, "Social": 2, "Print": 0},
    "overlaps": {
        ("TV", "Digital"): 1, ("TV", "Social"): 1, ("TV", "Print"): 0,
        ("Digital", "Social"): 1, ("Digital", "Print"): 0, ("Social", "Print"): 0,
    },
})
result = bounds(summary)
print(result.reach)                       # MetricBounds(minimum=3, maximum=4)
print(result.frequencies[3])              # MetricBounds(minimum=0, maximum=1)

fractional_only = Summary.from_mapping({
    "n": 4,
    "reaches": {c: 2 for c in ("TV", "Digital", "Social", "Print")},
    "overlaps": {
        ("TV", "Digital"): 1, ("TV", "Social"): 1, ("TV", "Print"): 1,
        ("Digital", "Social"): 1, ("Digital", "Print"): 1, ("Social", "Print"): 1,
    },
})
try:
    bounds(fractional_only)
except InfeasibleError:
    print("fractionally feasible but no integer population")  # 实际会走到这里
```

若需要验证某个端点对应的人群，可用 `freq_bounds.bounds.realize(summary, point)`
把最优点展开为完整的 16 格 Venn 计数（含未曝光格）。

## 精确性依据

四渠道的 Venn 图有 15 个非空格。用五个非负整数变量即可参数化全部格子：

- `q`：覆盖全部四个渠道的人数（频次 4）；
- `y_i`：覆盖除渠道 `i` 外其余三个渠道、但不覆盖 `i` 的人数（三元独占格）。

其余格子由输入唯一确定：

- 渠道对 `(i,j)` 的二元独占格

  `e_ij = overlap(i,j) - q - y_a - y_b`，其中 `{a,b}` 是另外两个渠道；

- 渠道 `i` 的一元格

  `s_i = reach_i − Σ(与 i 相关的三个重叠) + (y_j + y_k + y_l) + 2q`。

可行性恰好等价于三组线性不等式：

1. `q + y_a + y_b ≤ overlap(i,j)`（六个，二元独占格非负）；
2. `−(y_j + y_k + y_l) − 2q ≤ reach_i − Σ相关重叠`（四个，一元格非负）；
3. `(y_0+…+y_3) + 3q ≤ N − Σ触达 + Σ重叠`（总并集不超过 N）。

这些条件**充分且必要**：任何满足它们的非负整数 `(q, y_0…y_3)` 都使 15 个
Venn 格全部为非负整数，从而定义一个真实人群；反之任何真实人群都给出这样一个
点。所有目标指标（去重触达及各频次人数）都是这五个变量的线性函数，因此对上述
整数可行域分别做最小化/最大化，得到的就是真实极值。

整数优化用**两阶段单纯形（`fractions.Fraction` 精确有理运算）+ 分支定界**
完成，全程无浮点误差；整数不可行（仅分数模型可行）会被精确识别并报告。
测试中对 N ≤ 3 的全部 Venn 分布做了穷举交叉校验，并覆盖题目指定的两个场景。

## 运行测试

```bash
python3 -m unittest discover -s tests -v
```

仅使用标准库，兼容 Python 3.10+（开发验证于 Python 3.14）。
