# 跨媒体曝光频次区间计算库（frequency_bounds）

仅使用 Python 标准库实现。输入只接受**汇总数据**（总体人数、4 个渠道各自触达、
6 个渠道两两重叠），不接受逐人记录。对去重总触达人数与频次 1/2/3/4 的人数，
给出**真实可达的整数最小值与最大值**；不存在任何与汇总一致的整数人群时，
明确报告不可行。

## 输入约定

`Summary(n, tv, digital, social, print, tv_digital, tv_social, tv_print,
digital_social, digital_print, social_print)`

- 全部为非负整数（拒绝 `bool`、`float`、负数）。
- `n`：总体人数，**包含完全未曝光者**。
- `tv / digital / social / print`：各渠道触达人数。
- 六个重叠字段：同时覆盖这两个渠道的人数；**也包括**同时覆盖第三、第四个渠道的人。
- 频次 = 一个人覆盖的**不同渠道数**（同一渠道不重复计次）。
- 基本边际检查：单渠道触达 `<= n`；渠道对重叠 `<= min(两渠道触达)`。
  通过这些检查**不代表**一定存在真实人群，整体一致性由精确求解器判定。

## 输出约定

`bounds(summary)` 返回 `BoundsResult`：

- 不可行：`feasible=False`，其余字段为 `None`。
- 可行：`feasible=True`，并给出
  `reach`（去重总触达）与 `frequency_1 .. frequency_4` 的
  `MetricBounds(minimum, maximum)`。
- 每个端点都是某个与**全部**输入汇总一致的整数人群实际能达到的值；
  各指标**分别**取极值，不要求同一个人群同时达到所有端点。

另有 `is_feasible(summary) -> bool`。

## 使用方法

```python
import sys
sys.path.insert(0, "src")
from frequency_bounds import Summary, bounds

s = Summary(
    n=4,
    tv=2, digital=2, social=2, print=0,
    tv_digital=1, tv_social=1, tv_print=0,
    digital_social=1, digital_print=0, social_print=0,
)
r = bounds(s)
print(r.feasible)               # True
print(r.reach)                  # MetricBounds(minimum=3, maximum=4)
print(r.frequency_3)            # MetricBounds(minimum=0, maximum=1)

s2 = Summary(4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1)
print(bounds(s2).feasible)      # False（分数模型可行、整数人群不存在）
```

也可运行：`PYTHONPATH=src python3 examples.py`

## 精确求界的正确性依据

记 4 个渠道下 15 个非空渠道子集的人数为“原子计数”。仅给出单渠道与两两边际时，
这些原子计数（因而频次直方图）**不唯一**。本库做如下精确约简：

1. 定义 `D` = 覆盖全部 4 个渠道的人数；`c_i` = **缺少渠道 `i`** 的那个
   三人覆盖原子的人数（4 个三人原子，下标 `i` 表示被缺的渠道）。
   令 `C = c_0+c_1+c_2+c_3`，即频次恰为 3 的总人数。
2. 其余原子由输入边际与 `(D, c_0..c_3)` 唯一线性确定：
   - 仅覆盖渠道 `i` 的人数
     `a_i = s_i - Σ_{j≠i} p_ij + C - c_i + 2D`；
   - 恰好覆盖 `{i,j}` 的人数 `b_ij = p_ij - c_k - c_l - D`
     （`k,l` 为另外两个渠道）。
3. 因而“存在与汇总一致的整数人群”当且仅当 5 个非负整数变量
   `(D, c_0..c_3)` 满足全部 `a_i ≥ 0`、`b_ij ≥ 0` 以及总人数约束
   `R = S1 - P2 + C + 3D ≤ N`。这是一个只有 5 个变量的**整数线性规划**。
4. 各指标都是这些变量的线性函数：
   - `R   = S1 - P2 + C + 3D`
   - `F1  = S1 - 2·P2 + 3C + 8D`
   - `F2  = P2 - 3C - 6D`
   - `F3  = C`
   - `F4  = D`
   （其中 `S1 = Σ s_i`，`P2 = Σ_{i<j} p_ij`。
   恒等校验：`S1 = F1+2F2+3F3+4F4`，`P2 = F2+3F3+6F4`，`R = F1+F2+F3+F4`。）
5. 对 5 个目标分别做“最大化 / 最小化”，用 **Fraction 精确算术的两阶段单纯形**
   解线性松弛（Bland 规则，不循环；无浮点误差），再用**分支定界**求整数最优解。
   ILP 的整数最优值即“真实可达的整数极值”，而不是保守估计；不可行由
   第一阶段精确判定。分支定界对整数目标的 LP 上界做精确下取整剪枝，保证终止。
6. 整数性不可省：`tests/test_bounds.py` 中场景 2 表明，4 个变量各对之和 `≤ 1`
   而总和 `≥ 2` 时，分数解 `x_i = 1/2` 可行，但不存在整数解——这正是该汇总
   “允许分数人数模型、却不存在整数人群”的原因。

自动化测试另外用**独立的穷举 oracle**（枚举 N 较小时 15 个原子计数的所有分配，
推导对应汇总并统计每个指标的真实最小/最大值）与本库结果逐一比对：

- 默认运行 N=3 全量穷举；
- `RUN_SLOW_TESTS=1` 时额外运行 N=4 全量穷举。

## 运行测试

```bash
python3 -m unittest discover -s tests          # 快速套件（含 N=3 穷举对照）
RUN_SLOW_TESTS=1 python3 -m unittest discover -s tests   # 再加 N=4 全量穷举
```

## 代码结构

- `src/frequency_bounds/types.py`：`Summary` / `MetricBounds` / `BoundsResult`
  与 `bounds`、`is_feasible` 入口、输入校验。
- `src/frequency_bounds/ilp.py`：纯标准库精确 LP（两阶段单纯形）与 ILP（分支定界）。
- `src/frequency_bounds/core.py`：汇总 → 5 变量线性模型，6 组精确区间。
- `tests/`：规定场景、输入校验、ILP 单元测试、穷举对照。
