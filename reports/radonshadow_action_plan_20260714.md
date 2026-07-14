# RadonShadow 行动蓝图

**制定日期**: 2026-07-14 | **版本**: v1.0 | **预计周期**: 4-6 周

---

## 总览

| 阶段 | 内容 | 预估工时 | 产出 | 里程碑 |
|------|------|---------|------|--------|
| **W1** | 补齐实验缺口 | 2晚 | 完整数据表 | 实验矩阵填满 |
| **W2** | 理论整理 | 3晚 | 主要定理形式化 | 定理+证明草稿 |
| **W3** | 论文初稿 | 3-4晚 | 可投初稿 | 发给导师审阅 |
| **W4+** | 迭代+探索 | 灵活 | 论文终稿 | 准备投 arXiv |

---

## Week 1：补齐实验（2-3 个实验夜）

### 缺口清单

1. **Qutrit × 混合态**（态类型效应在 d=3 的验证）
   - 已有：qutrit 纯态（ratio 0.91-0.97）
   - 需要：qutrit 混合态（预测 ratio 更接近 1）
   - 作用：验证 3-qubit 的态类型效应在 qudit 上是否有普遍性

2. **大 d 扫描（d=13, 17, 19）的单项可观测量**
   - 已有：d=5,7,11（方向3）
   - 需要：扩展到 d=13,17,19 检验 d↗ 时确定性优势的收敛点
   - 作用：找到 d_max 使得确定性优势消失

3. **2-qubit 全态重建的 CS vs DPRT head-to-head**（含修正公式）
   - 已有：qudit (d=2,3,5,7) 全态 head-to-head ✅
   - 需要：2-qubit (d=4 非素数！) 的 MUB 线性反演（d=4 有 5 组 MUB）
   - 难点：d=4 不是素数，MUB 构造是否仍然完备？→ 需要查 Sibasish-Halder 2019 关于 d=4 的 MUB 完备性的文献

4. **Post-processing 复杂度对比**（DPRT 加法 vs CS 矩阵乘法）
   - 已有：复杂度理论分析（DPRT Ops 表格）
   - 需要：Python/NumPy 实际的 wall-clock 对比（不是大 O，是实际 ms）
   - 需要：内存占用对比（DPRT 只需 O(d²)、CS 需要存 N 个 d×d 矩阵）

### 每天实验夜产出
- 代码（追加到 `radon_shadow_phase2.py`）
- 结果记录（追加到 `qdprt_experiment_phase2.json`）

---

## Week 2：理论整理（3 个理论夜）

### 核心定理（需要形式化）

**定理 1：MUB 完备性 = DPRT 投影完备性**
- 陈述：对素数 d，(a) d+1 组 MUB 的投影算子张成算子空间 (b) d+1 个 DPRT 投影方向完全覆盖 $\mathbb{Z}_d^2$ (c) DPRT 反演和 MUB 线性反演在代数上等价
- 证明骨架：张成性来自 MUB 互不偏性质的 rank 论证；完全覆盖来自 DPRT 的离散 Fourier 切片定理；等价性来自 $R[m,k] \leftrightarrow M_m(k)$ 的映射

**定理 2：MUB 反演方差上界**
- 陈述：用 d+1 组 MUB 做确定性测量，有限样本下的密度矩阵估计 $\hat{\rho}$ 满足
  $$\mathbb{E}\|\hat{\rho} - \rho\|_F^2 \leq \frac{2d(d+1)}{N}$$
  其中 N 是总测量次数
- 证明骨架：霍夫丁不等式 + MUB 投影算子的正交性 + 矩不等式

**定理 3：确定性 vs 随机的相变点**
- 陈述：存在阈值 $d^* \approx 7$ 和 $k^* \approx 2$ 使得
  - 当 $d \leq d^*$ 或任意 d 下观测量的局域性 $k \leq k^*$ 时，确定性方案优
  - 当 $d > d^*$ 且 $k > k^*$ 时，随机方案优
- 证明需求：建立 DPRT 投影的方差分解公式（固有噪声 vs 基抽样方差）

### 文献阅读（每个理论夜配合）

1. **Sibasish-Halder 2019**：d=4,8 的 MUB 构造（非素数维 → 与 DPRT 的对应需要调整）
2. **Kingston 2007**：模 N 投影几何的精确重建条件（可以转写成 DPRT×MUB 的合数维推广定理）
3. **Huang-Kueng-Preskill 2020 SI**：与 DPRT 方差公式的数学对应

---

## Week 3：论文初稿（3-4 个写作夜）

### 论文结构

```
Title: Deterministic Classical Shadows via Complete
       Mutually Unbiased Basis Tomography

1. Introduction (1页)
   - Classical Shadow 的随机性：为什么随机？
   - 随机性的代价：基抽样方差
   - 我们提出：确定性方案（DPRT/完全MUB）

2. Preliminaries (2页)
   - Classical Shadow 回顾 (HKP 2020)
   - MUB 回顾 (Wootters-Fields 1989)
   - DPRT 回顾 (Hsung-Lun-Siu 1996)

3. The MUB-DPRT Correspondence (2页)  ← 核心贡献 1
   - 定理 1：MUB 完备性 = DPRT 语义
   - ρ = Σ pΠ - I：修正的反演公式
   - MUB 线性反演 vs Classical Shadow 单快照

4. Deterministic Shadow Protocol (2页)
   - 协议描述：输入 ρ、选维 d、测所有 d+1 MUB、用公式反演
   - 方差分析：定理 2（上界）

5. Numerical Experiments (3页)  ← 核心贡献 2
   - 实验矩阵：8 组对比实验
   - 3-qubit 的局域性分析（图 1）
   - d=5,7,11 的维度趋势（图 2）
   - 全态重建 DPRT vs CS head-to-head（表 1）
   - 复杂度对比（表 2）

6. Discussion (1页)
   - 限制：d 必须为素数（当前协议）
   - 两个推广方向：CRT 合数维 / 量子实现
   - 开放问题：d>11 的全面实验

7. References
```

### 每写作夜产出
- 论文骨架的 2 节
- 图/表（从实验数据导出）

---

## Week 4+：迭代与反馈

### 反馈循环

1. **发给导师审阅**（Week 3 末或 Week 4 初）
   - 附上 1 页中英文摘要
   - 标注需要导师研判的三个点：
     a. Kingston 2007 后为什么没有后续？
     b. 合数维的 DPRT 推广在信号处理领域还有人在做吗？
     c. DPRT 的 ODT（正交离散变换）性质能否用到量子层析的 proof？*

2. **根据导师反馈修改**（Week 4）
   - 可能需要调整实验设置（导师可能要求更大的 N 或更多的重复）
   - 可能需要重新形式化定理（导师可能指出 MUB-DPRT 代数对应上的遗漏）

3. **arXiv 发布准备**（Week 5-6）
   - LaTeX 排版（ieeeconf 或 revtex 格式）
   - 补充 SI（Supplementary Information）
   - 代码开源（GitHub repo 整理）

### 并行推进

在等待导师反馈时：
- 继续 d=13,17,19 的扫描（占坑实验）
- 阅读 Kingston 2007 的精确重建条件 → 尝试形式化推广定理
- 如果 R(5,5) 或 πe 有紧急进展，量子项目可以弹性暂停

---

## 风险与缓解

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| d=4 的 MUB 不完备 | 中 | 2-qubit 全态实验费劲 | 用 tensor product of d=2 MUBs（已知 3³=9>5 组约束） |
| 导师不看好 | 低 | 论文延期 | 作为学术探索记录，不投稿当练习 |
| 实验差距太小无法发表 | 中 | 需要更精巧的统计检验 | 用 bootstrap 或 Bayesian 重取样增强统计 power |
| Wigner 映射纠错时间太长 | 低 | 拖延理论进展 | 不走 Wigner 路线，直接用 MUB 完备性 |

---

## 每日检查点（自学版的节奏）

```
今晚  → 补齐 qutrit 混合态 + d=13 扫描
明晚  → 2-qubit 全态重建 + post-processing 复杂度
周三  → 读文献（Sibasish-Halder 2019 + Kingston 2007 关键页）
周四  → 定理 1 形式化（MUB = DPRT 代数等价）
周五  → 定理 2 形式化（方差上界）
周末  → 论文 Introduction + Preliminaries 初稿
```

---

## 附录 A：实验代码管理

| 文件 | 内容 | 状态 |
|------|------|------|
| `radon_shadow_experiments.py` | Phase 1: 2Q/3Q/qutrit 确定性 vs 随机 | ✅ 完成 |
| `radon_shadow_phase2.py` | Phase 2: DPRT 反演 + 大 d 扫描 | ✅ 完成 |
| `qdprt_experiment_phase2.json` | Phase 2 结果 | ✅ 保存 |
| `qdprt_experiment_phase3.json` | Phase 3 结果 | ✅ 保存 |
| `qdprt_*_analysis_*.md` | 分析文档 | ✅ 3 份 |

## 附录 B：关键公式速查

```
MUB 完备性反演:    ρ = Σ_{b,k} p(b,k)·|v_{b,k}⟩⟨v_{b,k}| - I_d
DPRT 反演:         f[i,j] = (Σ_m R[m,(j-mi)%d] + R[d,i] - Σ_k R[0,k]) / d
CS 单快照:          ρ̂_snap = (d+1)·|ψ⟩⟨ψ| - I_d
影子范数:           ||P||²_shadow = 3^k for k-local Pauli
MUB-DPRT 对应:     d+1 组 MUB ↔ d+1 个 DPRT 方向 (素数 d)
```

## 附录 C：给导师的讨论清单

约谈时准备的三个核心问题：

1. **"DPRT 的 ODT 性质和 MUB 的互不偏性质，在代数上是相同的数学结构（都是仿射群在有限域上的不变性）。这个对应关系之前有没有人注意过？"**

2. **"Kingston 2005/2007 把 DPRT 推到了 Galois 环和模 N 几何。这个方向为什么停在 2007 没有后续？是理论障碍还是应用场景不够？"**

3. **"如果把 DPRT 的确定性投影方向用在量子态测量上，本质上就是用一组完备的 MUB 做线性反演。数值实验在所有维度上优势都成立（3-18%）。这个方向值得深入吗？"**
