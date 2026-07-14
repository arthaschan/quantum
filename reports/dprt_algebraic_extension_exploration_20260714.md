# 离散Radon变换的代数扩张：从 ℤ_p 到 Galois 环 GR(p^n, m)

**研究方向探索报告 — 2026-07-14**
**与导师熊体操教授（Dr. Richard T.C. Hsung）DPRT 工作一脉相承**

---

## 0. 问题回顾

DPRT 的代数核心是：**$\mathbb{Z}_N$ 必须是域**，否则零因子使部分投影方向退化。

- 素数 $N=p$ → $\mathbb{Z}_p$ 是域 → 每条"直线"完整覆盖 → 逆变换可行
- 合数 $N$ → $\mathbb{Z}_N$ 是环 → 零因子方向 → 直线退化 → 信息丢失

**你的追问**：能否沿 Kingston 等人开辟的方向，用更深刻的代数结构（Galois环、有限域扩张、有限几何）突破素数限制？

**答案**：**可以，而且这条路线有多个维度的潜力尚未挖掘。** 以下系统展开。

---

## 1. 已知路线图：四条代数扩张路径

### 1.1 路径A：Galois环 GR(p^n, m) —— Kingston 2005 的推广

**核心思想**：将 $\mathbb{Z}_p$（域）提升为 Galois 环 $\operatorname{GR}(p^n, m)$，后者是环论中最接近域的结构。

**Galois 环的定义**：
$$\operatorname{GR}(p^n, m) = \mathbb{Z}_{p^n}[x]/(f(x))$$

其中 $f(x)$ 是 $\mathbb{Z}_{p^n}$ 上的 $m$ 次首一基本不可约多项式。$\operatorname{GR}(p^n, m)$ 有 $p^{nm}$ 个元素，**所有零因子都在极大理想 $(p)$ 中**。

**关键性质**（使 Kington 的推广可行）：
1. 全体单位构成乘法群 $\operatorname{GR}(p^n, m)^\times$，大小为 $p^{nm} - p^{n(m-1)}$
2. 每个非零元要么是单位，要么是 $p$ 的倍数
3. 存在 Teichmüller 代表元集，同构于 $\mathbb{F}_{p^m}$
4. 有一个 Galois 自同构群，由 Frobenius 映射生成

**Kingston 2005 的做法**：在 $p^n \times p^n$ 图像上定义 ODPRT，投影方向取自 $\operatorname{GR}(p^n, 1) = \mathbb{Z}_{p^n}$，但**投射空间的几何来自 $\mathbb{F}_p$ 上的仿射平面**。关键技巧是对每个投影切片运用 **Gray 映射**（从 $\mathbb{Z}_{p^n}$ 到 $\mathbb{F}_p^n$ 的保距映射）来展开为二进制平面，然后对每个平面重建。

**现状与局限**：
- Kingston 2005 解决了 $p^n \times p^n$，但对一般合数 $N = p_1^{e_1} p_2^{e_2} \cdots p_k^{e_k}$ 未解决
- 多个素因子需要 CRT 分解，然后处理混合 Galois 环结构

### 1.2 路径B：模 N 投影几何 —— Kingston & Svalbe 2006-2007

**2006年论文**："A Discrete Modulo N Projective Radon Transform for N×N Images"（DGCI 会议）
**2007年论文**："Generalised finite radon transform for N×N images"（Image and Vision Computing, Vol. 25, pp. 1620–1630）

**核心思想**：直接在 **$\mathbb{Z}_N$（任意 $N$，不一定素数）** 上定义 Radon 变换，但通过**巧妙选择投影方向**避开零因子问题。

**关键创新**：不取所有 $N(N+1)$ 个投影方向，改为选取一组**平移方向集合** $\Theta$，其中每个方向 $(a,b)$ 满足：
$$\gcd(a, b, N) = 1$$

即 $a$ 和 $b$ 在 $\mathbb{Z}_N$ 中不能有公共零因子。这确保了每条投影直线上的采样点数是 $N$ 而不是更少。

**投影定义**：对于方向 $(a,b) \in \Theta$，投影切片为：
$$R_{(a,b)}(t) = \sum_{x=0}^{N-1} \sum_{y=0}^{N-1} I(x,y) \cdot \delta(\langle ax + by - t \rangle_N)$$

**精确重建条件**：需要选取足够的投影方向使得 $\mathbb{Z}_N$ 上的离散 Fourier 切片定理覆盖所有频点。具体来说，对每个 $(u,v) \in \mathbb{Z}_N^2$，需要存在至少一个方向 $(a,b) \in \Theta$ 和某个整数 $\lambda$ 满足：
$$(u,v) = \lambda \cdot (a,b) \pmod N$$

**优点**：
- 对**任意 $N$** 成立（不限于素数或素数幂）
- 保持了 DPRT 的离散 Fourier 切片定理
- 投影数 $|\Theta|$ 为 $O(N \cdot \psi(N))$，其中 $\psi(N) = N \prod_{p|N}(1+1/p)$ 是 Dedekind ψ 函数

### 1.3 路径C：Mojette 变换 —— Guédon 1995 至今

**核心思想**：放弃模 $N$ 的"循环"几何，改用**有理斜率**的投影方向。

**投影方向**：$(p_i, q_i)$ 满足 $\gcd(p_i, q_i) = 1$（互素）。投影定义在 $\mathbb{Z}^2$ 上：
$$M_{p,q}(b) = \sum_{k} I(k, \lfloor (b - pk)/q \rfloor)$$

**精确重建条件**（Katz 准则）：
$$\sum_i |p_i| \ge N-1 \quad \text{或} \quad \sum_i |q_i| \ge N-1$$

**与 DPRT 的关系**：
- Mojette 不要求 $N$ 为素数 → 天然支持任意尺寸
- 但它没有 DPRT 的**无乘法逆变换**——Mojette 的逆变换需要求解稀疏线性方程组
- Mojette 的投影数是 $O(N)$（而非 $N$ 或 $N+1$），需要更多投影来重建

### 1.4 路径D：有限域扩张 𝔽_{p^k} —— Wootters 1987 方向

**核心思想**：将信号定义为 $\mathbb{F}_{p^k}$ 上的函数，而非 $\mathbb{Z}_N$ 上的函数。$\mathbb{F}_{p^k}$ 天然是域，所有非零元有逆元。

信号：$x: \mathbb{F}_q \to \mathbb{C}$（或更一般地 $x: \mathbb{F}_q \times \mathbb{F}_q \to \mathbb{C}$），其中 $q = p^k$。

投影方向：$\mathbb{F}_q$ 中的元素 $m$。
投影定义为 $\mathbb{F}_q$ 上的仿射直线：
$$X_m(d) = \sum_{n \in \mathbb{F}_q} x(n) \cdot \delta_{\mathbb{F}_q}(d - m \cdot n)$$

其中 $\delta_{\mathbb{F}_q}$ 是 $\mathbb{F}_q$ 上的 Kronecker delta。

**离散 Fourier 切片定理**：利用 $\mathbb{F}_q$ 上的加法特征 $\chi(a) = e^{2\pi i \operatorname{Tr}(a)/p}$，可以得到完整的 Fourier 切片定理。

**优点**：$q$ 可取任意素数幂，$\mathbb{F}_q$ 是完美的域。
**缺点**：信号长度只能是素数幂，不是任意合数。加法特征是复指数，引入了**浮点运算**（丢失了 DPRT 的无乘法性质）。

---

## 2. 尚未探索的前沿方向

以上四条路径均有公开文献，以下是我认为**值得进一步探索的未充分开辟方向**：

### 2.1 方向①：混合模数 —— 用 CRT 分解合数 N，但保持方向耦合

**现状**：Kingston & Svalbe 2007 直接在 $\mathbb{Z}_N$ 上工作，通过限制投影方向规避零因子。但这种方式丢失了 $\mathbb{Z}_N$ 上的某些"自然"投影方向。

**新想法**：
$$\mathbb{Z}_N \cong \prod_{i=1}^k \mathbb{Z}_{p_i^{e_i}}$$

不是简单地在每个分量上独立做 DPRT 再组合，而是利用 CRT 的**显式同构映射**构造一个"提升"——在分量的 DPRT 域之间引入**交叉约束**。

具体来说，CRT 给出了双射：
$$\Phi: \mathbb{Z}_N \to \mathbb{Z}_{p_1^{e_1}} \times \cdots \times \mathbb{Z}_{p_k^{e_k}}$$

在混合域中定义投影方向 $(m_1, \ldots, m_k)$，每个分量 $m_i \in \mathbb{Z}_{p_i^{e_i}}$。但关键的创新点是：**投影位移 $d$ 应全局映射**而非分量独立。

$$X_{(m_1,\ldots,m_k)}(d) = \sum_{n=0}^{N-1} x(n) \cdot \delta_N(d - M \cdot n)$$

其中 $M = \Phi^{-1}(m_1,\ldots,m_k)$ 是 CRT 逆映射给出的全局乘法器。

**核心问题与机会**：
- $M$ 在 $\mathbb{Z}_N$ 中是否是零因子？→ 只要有一个 $m_i$ 在其分量中是零因子，$M$ 在全局就是零因子
- 但如果**不同时选择全部零因子方向**，可以通过约束 $m_i$ 的选择来构造"几乎完整"的投影集
- 这比 Kingston 2007 的穷举方向更精细——用分量结构做筛选

### 2.2 方向②：Galois 环上的一般 Radon 变换（超越 2D）

**现状**：Kingston 2005 仅处理了 $p^n \times p^n$ 的 2D 图像。**对 $d$ 维推广到 $\operatorname{GR}(p^n, m)^d$ 尚无人做。**

**新想法**：
在 $d$ 维 Galois 环网格上，投影超平面定义为：
$$\{(x_1,\ldots,x_d) \in \operatorname{GR}^d : a_1 x_1 + \cdots + a_d x_d = t\}$$

其中 $(a_1,\ldots,a_d)$ 是投影方向。对于 $\operatorname{GR}(p^n, m)$，**只要至少有一个 $a_i$ 是单位（非零因子）**，该超平面的大小就是 $|\operatorname{GR}|^{d-1}$，且每个 $(a_1,\ldots,a_d)$ 生成一个正则纤维化。

**潜在应用**：
- 4D 时空信号（3D + 时间）的精确层析重建
- 多维量子态层析中的离散 Wigner 函数
- 高维张量数据的压缩感知

### 2.3 方向③：DPRT 的量子推广 —— Quantum DPRT

**这个方向与你正在编写的量子层析教材直接衔接！**

**灵感来源**：Wootters 1987 在 $\mathbb{F}_q$ 上定义了离散 Wigner 函数。DPRT 可以自然地推广到量子设定：

**qDPRT**：在 $d$ 维量子系统（$d$ 不一定素数）上，定义离散相空间 $\mathbb{Z}_d \times \mathbb{Z}_d$（位置×动量）。广义 Pauli 算子：
$$X = \sum_{k=0}^{d-1} |k+1\rangle\langle k|, \quad Z = \sum_{k=0}^{d-1} \omega^k |k\rangle\langle k|$$

其中 $\omega = e^{2\pi i/d}$。相空间点算子：
$$A_{(q,p)} = X^q Z^p \omega^{-qp/2}$$

DPRT 方向 $(a,b)$ 对应的投影相当于对相空间沿方向 $(a,b)$ 做"切片"：
$$\Pi_{(a,b)}(t) = \sum_{q,p: aq+bp=t} A_{(q,p)}$$

**这构建了一个"量子版本的投影算子"**，其测量结果的概率分布就是量子态在该方向上的离散投影。

**与 Classical Shadows 的联系**：
- Classical Shadows (Huang, Kueng, Preskill 2020) 的核心是随机 Pauli 测量
- qDPRT 视角下，Classical Shadows 可以理解为在**随机 Radon 投影方向**上做测量
- DPRT 的正交性保证了这些投影方向的信息完备性
- **一个新的影子层析协议**：用 DPRT 投影方向集合替代随机 Pauli → 确定性投影 → 可能达到更优的样本复杂度

**这是一个具有原创性的交叉方向，值得深入。**

### 2.4 方向④：用 Tropical 几何统一零因子分析

当 $N$ 为合数时，$\mathbb{Z}_N$ 上的"直线"退化现象可以统一用 Tropical 几何描述。

**基本观察**：方程 $d \equiv a x + b y \pmod{N}$ 的解集大小取决于 $\gcd(a,b,N)$。
令 $g = \gcd(a,b,N)$，则解集大小为 $N/g$，一共有 $g$ 个可到达的 $d$ 值。

**Tropical 视角**：可以将 $\mathbb{Z}_N$ 上的仿射几何视为一个 **Tropical 簇**上的截面图。当 $g>1$ 时，解集从 $N$ 个点退化到 $N/g$ 个点，对应于 Tropical 退化的维数下降。

潜在数学工具：
- Tropical 矩阵的秩和零空间
- Baker-Norine 除子理论的离散类比
- 这可能在 $N$ 的素因子之间建立"信息流"的代数描述

### 2.5 方向⑤：卷积神经网络（CNN）视角 —— 学习自适应 DPRT 方向

**这是一个 AI+DPRT 的交叉方向**：

传统 DPRT 的投影方向是人工选定的（模 $N$ 的仿射直线）。可以训练一个轻量 CNN 来**学习最优投影方向集合**：

1. 初始化一个可微分的"投影算子"（参数化投影方向）
2. 针对特定任务（去噪/重建/压缩），端到端训练投影方向
3. 学到的方向可能落在传统 DPRT 方向之间（非整数斜率）

**与导师 DPRT 工作的自然衔接**：
- DPRT 的无乘法逆变换可以作为网络中的**固定层**（不可训练，保证重建精度）
- 学习的投影方向可以补足固定 DPRT 方向 → 混合 DPRT + Learnable Projections
- 这对量子态层析中的自适应测量基选择有直接应用

---

## 3. 最优先推荐：两个高潜力方向

### 🥇 推荐一：qDPRT × Classical Shadows —— 量子+Radon 交叉

**为什么**：
- 与你正在编写的量子层析教材**直接衔接**
- 与导师的 DPRT 工作**一脉相承**（DPRT → 量子推广 → 层析应用）
- Classical Shadows 是当前量子计算最热门的方向之一（2020年以来 > 2000 引用）
- 这是一个**真正有原创性**的交叉点，目前文献中无人做过

**可发表的潜在贡献**：
1. 将 DPRT 投影方向集合引入 Classical Shadows → 新的确定性影子协议
2. 比较随机 Pauli vs DPRT 投影的样本复杂度下界
3. DPRT 的正交性 → 保证 shadow norm 的最优性（或否证）

### 🥈 推荐二：混合模数 CRT 提升 DPRT

**为什么**：
- 这是一个干净的数学问题，适合硕士阶段
- 比 qDPRT 更"纯数学"，不需要量子力学前置
- 结果是 DPRT 理论本身的直接推广
- 可以与导师讨论后快速形成论文框架

### 🥉 其他推荐

| 方向 | 难度 | 原创性 | 与导师关系 | 适合阶段 |
|------|------|--------|-----------|---------|
| qDPRT × Shadows | 中高 | ⭐⭐⭐⭐⭐ | 间接衔接 | 论文核心 |
| 混合模数 CRT | 中 | ⭐⭐⭐⭐ | 直接延续 | 硕士论文 |
| Galois 环 d-D | 高 | ⭐⭐⭐ | 直接延续 | 博士方向 |
| Tropical 统一 | 高 | ⭐⭐⭐⭐ | 间接 | 理论探索 |
| 学习 DPRT 方向 | 低中 | ⭐⭐⭐ | 间接 | 应用论文 |

---

## 4. 写给导师的讨论要点

如果要和 Prof. Hsung 讨论这个方向，可以这样切入：

1. **"Kingston 2005/2007 把 DPRT 推到了 Galois 环和模 N 几何，但似乎停在 2007 就没有后续了。"**
   - 导师可能知道 Kingston（澳大利亚 Monash 大学，物理学院）的后续工作
   - 可以问导师：这个方向为什么没有再推进？是遇到了理论障碍，还是应用场景不够？

2. **"我在想，DPRT 的正交性和无乘法性质，在量子层析中有没有用？"**
   - 这引出 qDPRT 的想法
   - 导师 1996 的 DPRT 在量子域中的自然类比是一个自然的学术传承

3. **"合数的 CRT 分解能不能不丢失方向耦合？"**
   - 这是混合模数方向的核心问题
   - 导师对 CRT 和 DPRT 几何的理解可能直接给出答案或反例

---

## 参考文献

1. **Hsung, Lun, Siu** (1996). The discrete periodic Radon transform. *IEEE TSP*, 44(10), 2651–2657.
2. **Lun, Hsung, Shen** (2003). Orthogonal DPRT Part I: Theory and realization. *Signal Processing*, 83(5), 939–955.
3. **Lun, Hsung, Shen** (2003). Orthogonal DPRT Part II: Applications. *Signal Processing*, 83(5), 957–971.
4. **Kingston** (2005). Orthogonal discrete Radon transform over $p^n \times p^n$ images. *Signal Processing*, 86(8), 2040–2050.
5. **Kingston, Svalbe** (2006). A discrete modulo N projective Radon transform for N×N images. *DGCI 2006*, LNCS 4245, pp. 136–147.
6. **Kingston, Svalbe** (2007). Generalised finite Radon transform for N×N images. *Image and Vision Computing*, 25(10), 1620–1630.
7. **Guédon et al.** (1995–). The Mojette transform: The first ten years. *DGCI 2005*, LNCS 3429, pp. 79–91.
8. **Wootters** (1987). A Wigner-function formulation of finite-state quantum mechanics. *Annals of Physics*, 176(1), 1–21.
9. **Huang, Kueng, Preskill** (2020). Predicting many properties of a quantum system from very few measurements. *Nature Physics*, 16, 1050–1057.
10. **Matúš, Flusser** (1993). Image representation via a finite Radon transform. *IEEE TPAMI*, 15(10), 996–1006.
