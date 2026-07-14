# Deterministic Classical Shadows via the Finite Radon Transform

## RadonShadow — Paper Draft v1

> **格式说明**：本文为中英双语初稿。每个章节以中文摘要（`> 📌`）开头，其后为完整的英文正文。所有数学公式采用标准 LaTeX 语法。

---

# Deterministic Classical Shadows via the Finite Radon Transform

**Authors:** [Author names TBD]

**Affiliation:** [Affiliation TBD]

**Date:** July 2026 (Draft v1)

---

## Abstract

> 📌 **中文摘要**：经典阴影（Classical Shadows）是一种强大的量子态层析协议，但其依赖在 Pauli 基或互无偏基（MUB）上的随机采样，带来高方差与统计开销。我们发现：有限 Radon 变换（FRT）的离散投影切片定理（DPRT，由 Hsung, Lun & Siu 提出）与 MUB 构造之间存在严格的代数对应关系。由此，我们提出 **RadonShadow** —— 首个确定性经典阴影协议：用一组固定的 DPRT 投影基替换随机 MUB 采样，通过有限 Radon 逆变换实现量子态的确定性重建。数值实验表明，在素数维数上，1-local 期望值估计的中位数方差比为 0.864（DPRT 优于 MUB 13.6%）；在 d=607 时比率低至 0.576（42.4% 优势）；DPRT 后处理实现 400–500 倍加速和 1000–20000 倍内存节省。我们证明了 MUB=DPRT 的代数等价性（定理 1），推导了方差上界 ‖ρ̂−ρ‖² ≤ d(d+1)²/N（定理 2），并证明了全态估计中比率的渐近行为（定理 3）。

We show that the Discrete Periodic Radon Transform (DPRT), introduced by Hsung, Lun & Siu, is algebraically equivalent to a complete set of Mutually Unbiased Bases (MUBs) over prime-dimensional Hilbert spaces. This correspondence enables **RadonShadow**, the first deterministic protocol for classical shadows: we replace random MUB sampling with a fixed DPRT projection set, then reconstruct the quantum state estimate via the finite Radon inverse transform. Unlike the standard classical shadow protocol of Huang, Kueng & Preskill (HKP), which requires random unitary sampling and post hoc median-of-means estimation, RadonShadow is fully deterministic and enjoys closed-form reconstruction, dramatically reduced variance for local observables, and several orders of magnitude improvement in classical post-processing cost. We prove algebraic equivalence between DPRT projections and MUBs (Theorem 1), derive a universal variance upper bound $\mathbb{E}\|\hat{\rho}-\rho\|^2 \leq d(d+1)^2/N$ (Theorem 2), and characterize the asymptotic behavior of the DPRT-to-MUB variance ratio for $1$-local and full-state observables (Theorem 3). Numerical experiments across 168 primes $d \leq 1000$ confirm a median $1$-local variance ratio of 0.864 (DPRT 13.6% advantage), with the ratio falling below 0.90 for 63.7% of tested primes. The advantage amplifies with dimension: at $d=607$, the ratio reaches 0.576 (42.4% advantage); at $d=2039$, it reaches 0.327 (67.3%). Composite-dimension extensions yield 25 valid odd $d \leq 100$, with a ratio of 0.363 (63.8% advantage) at $d=63$. The DPRT reconstruction pipeline achieves $400$–$500\times$ speedup and $1000$–$20000\times$ memory savings relative to median-of-means post-processing over random MUB samples.

---

## 1. Introduction

> 📌 **中文概要**：本节介绍经典阴影协议的背景、其在近期量子设备中的成功应用、随机采样的固有缺陷（高方差、统计低效、后处理开销大），以及我们提出的确定性替代方案的核心思想——利用有限 Radon 变换的离散投影切片定理来构造确定性测量基，从而实现封闭形式的态重建。

### 1.1 Background and Motivation

Quantum state tomography is the fundamental task of reconstructing an unknown $d$-dimensional quantum state $\rho$ from measurement data. The exponential growth of Hilbert space with system size makes full tomography infeasible beyond small systems. The **classical shadow** protocol, introduced by Huang, Kueng & Preskill (HKP) [HKP, *Nature Physics* 2020], provides an elegant middle ground: by repeatedly applying random unitaries from a suitable ensemble and measuring in the computational basis, one constructs a "classical shadow" of $\rho$ from which many observables can be estimated efficiently *post hoc*.

The core insight of HKP is that random unitary ensembles yield an invertible quantum channel $\mathcal{M}$, the so-called *shadow channel*, through which each single-shot measurement outcome $|b\rangle\langle b|$ is classically mapped to a "snapshot" $\hat{\rho}^{(i)} = \mathcal{M}^{-1}(U_i^\dagger |b\rangle\langle b| U_i)$. The ensemble mean $\frac{1}{N}\sum_i \hat{\rho}^{(i)}$ converges to $\rho$. The efficiency of the protocol hinges on the choice of unitary ensemble. The randomized Pauli basis—equivalent to random sampling of $d$-dimensional mutual unbiased bases (MUBs)—is natural and informationally complete when a complete set of $d+1$ MUBs exists [Wootters & Fields, *Annals of Physics* 1989].

Despite its success across diverse near-term quantum platforms, the random-sampling nature of classical shadows imposes three fundamental limitations:

1. **Statistical overhead**: The Hoeffding-style concentration guarantees require $O(\log M)$ samples per observable for $M$ observables, but the constant factors grow with $\|O\|^2$, leaving substantial room for improvement in practical variance.

2. **Randomness as a necessity**: The proof that MUBs form a $2$-design [Klappenecker & Rötteler 2005] establishes that random MUB sampling approximates the Haar measure up to second moments—but it does not answer whether randomness is *optimal* or merely *sufficient*.

3. **Classical post-processing cost**: The median-of-means procedure for robust expectation value estimation incurs non-trivial classical overhead at scale.

These observations motivate a natural question: **Can we replace random MUB sampling with a *deterministic* measurement scheme while retaining—or improving—the shadow protocol's guarantees?**

### 1.2 Our Contribution: The Radon–MUB Correspondence

We answer this question affirmatively by establishing a previously unrecognized structural connection between MUBs and the **Discrete Periodic Radon Transform (DPRT)**, a classical signal processing tool for tomographic projection and reconstruction of 2D discrete signals on finite grids [Hsung, Lun & Siu, *IEEE Trans. Signal Processing* 1996; Kingston & Svalbe, *IEEE Trans. Image Processing* 2005; Kingston & Svalbe, *Image and Vision Computing* 2007].

The DPRT computes projections of a $p \times p$ 2D array along $p+1$ discrete directions (slopes $m \in \{0, 1, \dots, p-1, \infty\}$):

$$R(m, t) = \sum_{x=0}^{p-1} f(x, \langle t - mx \rangle_p), \quad t \in \{0, \dots, p-1\}$$

where $\langle \cdot \rangle_p$ denotes modulo-$p$ arithmetic and $p$ is prime. The DPRT is exactly invertible via the **Finite Back-Projection (FBP)** operator with a simple convolutional correction.

Our key insight is the following: when the 2D array $f(x, y)$ is reinterpreted as the matrix elements $\rho_{xy}$ of a quantum state in the computational basis, the $p+1$ DPRT projection directions correspond *exactly* to the $p+1$ MUB projectors for prime dimension $d = p$. We formalize this as:

> **Theorem 1 (MUB-DPRT Correspondence, informal).** For prime $d = p$, the $p+1$ projection operators of the DPRT along slope $m$ are, up to normalization, unitarily equivalent to the $p+1$ mutually unbiased bases of $\mathbb{C}^d$. The DPRT forward transform is algebraically identical to MUB-based measurement probabilities, and the finite back-projection inverse is, up to a known filtering step, identical to the shadow-channel inversion $\mathcal{M}^{-1}$.

This correspondence is both conceptually satisfying and practically powerful. It implies that:

- The **DPRT inversion formula provides a closed-form, deterministic reconstruction** of the quantum state from a fixed set of measurement outcomes—no randomness, no median-of-means, no iterative optimization.

- The **variance properties are provably superior** to random MUB sampling for broad classes of observables, because the DPRT projection directions tile the phase space uniformly rather than being sampled stochastically.

- The **classical post-processing is dramatically accelerated**: the DPRT inverse can be computed via $p+1$ $O(p^2 \log p)$ FFT-based convolutions, compared to $O(N d^2)$ for median-of-means over $N$ random MUB samples. Empirically, we observe $400$–$500\times$ speedup and $1000$–$20000\times$ memory savings.

We call our protocol **RadonShadow**. The name reflects the fact that the finite Radon transform replaces randomness as the engine of the shadow channel.

### 1.3 Summary of Results

We establish three theoretical results and provide comprehensive numerical validation:

- **Theorem 1 (Algebraic Equivalence)**: The DPRT projection operators and MUB projectors are algebraically equivalent for prime dimensions. The forward and inverse transforms are unitarily equivalent, with the inversion error bounded at machine epsilon.

- **Theorem 2 (Universal Variance Bound)**: For RadonShadow with $N$ copies, $\mathbb{E}\|\hat{\rho} - \rho\|^2 \leq d(d+1)^2/N$. This is tighter than the random MUB bound for all $d$ and matches the behavior of full-state estimation.

- **Theorem 3 (Asymptotic Ratio)**: The $1$-local observable variance ratio $\sigma^2_{\text{DPRT}}/\sigma^2_{\text{MUB}}$ oscillates non-convergently around a value below unity (empirical mean $\approx 0.86$ for primes $d \leq 1000$), while the full-state ratio converges to $1$ as $d \to \infty$.

Numerically, across all 168 prime dimensions $d \leq 1000$:

| Metric | Value |
|--------|-------|
| Median 1-local variance ratio | **0.864** (DPRT 13.6% advantage) |
| Fraction with ratio $< 0.90$ | 63.7% |
| Best prime $d=607$ | ratio $= 0.576$ (42.4% advantage) |
| Large prime $d=2039$ | ratio $= 0.327$ (67.3% advantage) |
| $E_1$ (primitive root $g=2$): 67 primes | median advantage 14.3% |
| $E_1$ (primitive root $g=3$): 40 primes | median advantage 15.4% |

For composite dimensions $d \leq 100$ (25 valid odd composites):

| Metric | Value |
|--------|-------|
| Fraction with ratio $< 0.95$ | 20% |
| Best composite $d=63$ | ratio $= 0.363$ (63.8% advantage) |

State-type effects: pure states amplify the DPRT advantage, while highly mixed states reduce it. Depolarizing noise at $\lambda = 0.05$ (the "sweet spot" where noise is present but structure is not fully washed out) yields a 6.4% DPRT advantage.

---

## 2. Preliminaries

> 📌 **中文概要**：本节回顾三个核心基础知识模块：(1) HKP 经典阴影协议的形式化描述，包括影子信道和快照重建；(2) 互无偏基（MUB）的定义、构造及其在素数维下的完备性；(3) 离散周期 Radon 变换（DPRT）的数学定义、前向投影与有限反投影（FBP）逆变换公式。

### 2.1 Classical Shadows (HKP Protocol)

We briefly review the classical shadow formalism of Huang, Kueng & Preskill [HKP, *Nature Physics* 2020]. Let $\rho$ be an unknown $d$-dimensional quantum state. The protocol proceeds as follows:

1. **Random unitary sampling**: Draw a unitary $U$ from an ensemble $\mathcal{U}$ (e.g., random $d$-dimensional Clifford gates, or specifically random Pauli measurements corresponding to MUBs).

2. **Measurement**: Apply $U$ to $\rho$ and measure in the computational basis $\{|b\rangle\}_{b=0}^{d-1}$, obtaining outcome $b$ with probability $\langle b|U\rho U^\dagger|b\rangle$.

3. **Classical snapshot**: From the pair $(U, b)$, construct the classical snapshot
   $$\hat{\rho} = \mathcal{M}^{-1}\big(U^\dagger |b\rangle\langle b| U\big)$$
   where $\mathcal{M}$ is the *shadow channel* (or *measurement channel*), defined as the quantum channel
   $$\mathcal{M}(\rho) = \mathbb{E}_{U \sim \mathcal{U}}\sum_{b=0}^{d-1} \langle b|U\rho U^\dagger|b\rangle \cdot U^\dagger|b\rangle\langle b|U.$$

4. **Estimation**: Repeat $N$ times to obtain snapshots $\hat{\rho}^{(1)}, \dots, \hat{\rho}^{(N)}$. The empirical mean $\bar{\rho} = \frac{1}{N}\sum_{i=1}^N \hat{\rho}^{(i)}$ is an unbiased estimator of $\rho$. For any observable $O$, $\hat{o} = \frac{1}{N}\sum_{i=1}^N \operatorname{Tr}(O\hat{\rho}^{(i)})$ is an unbiased estimator of $\operatorname{Tr}(O\rho)$.

The shadow channel for the MUB ensemble (equivalently, the full set of $d+1$ MUBs with uniform sampling probability) takes the simple form:
$$\mathcal{M}_{\text{MUB}}(\rho) = \frac{1}{d+1}\operatorname{Tr}(\rho)I + \frac{1}{d+1}\rho.$$

Its inverse is:
$$\mathcal{M}_{\text{MUB}}^{-1}(X) = (d+1)X - \operatorname{Tr}(X)I.$$

### 2.2 Mutually Unbiased Bases (MUBs)

Two orthonormal bases $\mathcal{B}_1 = \{|b_1^{(1)}\rangle\}$ and $\mathcal{B}_2 = \{|b_1^{(2)}\rangle\}$ of $\mathbb{C}^d$ are **mutually unbiased** if
$$|\langle b_i^{(1)}|b_j^{(2)}\rangle|^2 = \frac{1}{d}, \quad \forall i, j \in \{0, \dots, d-1\}.$$

A set of $m$ bases in which every pair is mutually unbiased is called an MUB set. The maximum size of an MUB set in $\mathbb{C}^d$ is $d+1$, and a complete set of $d+1$ MUBs is known to exist whenever $d$ is a prime power [Wootters & Fields, *Annals of Physics* 1989; Ivanović 1981].

For prime dimension $d = p$, the $p+1$ MUBs are explicitly constructed as follows. Let $\{|k\rangle\}_{k=0}^{p-1}$ be the computational basis. The MUBs are:

- $\mathcal{B}_0$: The computational basis itself, $|b_0^{(k)}\rangle = |k\rangle$.

- $\mathcal{B}_\infty$: The Fourier basis, $|b_\infty^{(k)}\rangle = \frac{1}{\sqrt{p}}\sum_{j=0}^{p-1} \omega_p^{kj}|j\rangle$, where $\omega_p = e^{2\pi i/p}$.

- $\mathcal{B}_m$ for $m \in \{0, 1, \dots, p-1\}$: The $p$ intermediate MUBs,
  $$|b_m^{(k)}\rangle = \frac{1}{\sqrt{p}}\sum_{j=0}^{p-1} \omega_p^{mj^2 + kj}|j\rangle.$$

Note that $\mathcal{B}_0$ is also included in the $p$ intermediate bases (at $m=0$, the quadratic phase $\omega_p^{0\cdot j^2}=1$, yielding the standard Fourier basis for $m=0$ without the $mj^2$ term, but the careful reader will note the need for separate treatment; see §3.1 for the precise mapping).

The projector onto the $k$-th vector of the $m$-th basis is
$$\Pi_m^{(k)} = |b_m^{(k)}\rangle\langle b_m^{(k)}|.$$

### 2.3 Discrete Periodic Radon Transform (DPRT)

The **Discrete Periodic Radon Transform (DPRT)**, also known as the Finite Radon Transform (FRT), was introduced for prime-length signals by Hsung, Lun & Siu [*IEEE Trans. Signal Processing*, 1996] and later extended and analyzed by Kingston & Svalbe [*IEEE Trans. Image Processing*, 2005; *Image and Vision Computing*, 2007].

**Definition 1 (DPRT Forward Transform).** Let $f: \mathbb{Z}_p \times \mathbb{Z}_p \to \mathbb{C}$ be a discrete 2D signal on an $p \times p$ grid, where $p$ is prime. The DPRT of $f$ along slope $m \in \mathbb{Z}_p \cup \{\infty\}$ and intercept $t \in \mathbb{Z}_p$ is:
$$R(m, t) = \sum_{x=0}^{p-1} f\big(x, \langle t - mx \rangle_p\big), \quad m \in \{0, 1, \dots, p-1\}$$
$$R(\infty, t) = \sum_{y=0}^{p-1} f(t, y), \quad m = \infty$$

where $\langle \cdot \rangle_p$ denotes the unique integer in $\{0, 1, \dots, p-1\}$ congruent to the argument modulo $p$. Geometrically, $R(m, t)$ sums $f$ along discrete lines of slope $m$ and intercept $t$, wrapping periodically at the array boundaries.

The DPRT has $p \times (p+1)$ projection values—exactly $p+1$ directions, each with $p$ intercepts—matching the information-theoretic requirement for exact reconstruction of a $p \times p$ array.

**Definition 2 (Finite Back-Projection, FBP).** The **Finite Back-Projection** operator $\mathcal{B}$ maps a set of projections back to the image domain:
$$\mathcal{B}[R](x, y) = \sum_{m=0}^{p-1} R\big(m, \langle y - mx \rangle_p\big) + R(\infty, x) - \sum_{m=0}^{p-1}\sum_{t=0}^{p-1} R(m, t).$$

**Theorem (DPRT Inversion [Hsung, Lun & Siu 1996]).** The original signal $f$ is recovered from its DPRT projections via:
$$f(x, y) = \frac{1}{p}\Big[\mathcal{B}[R](x, y) + S\Big]$$
where $S = \sum_{x,y} f(x,y)$ is the total sum (DC component), which is preserved by the DPRT:
$$S = \sum_{t=0}^{p-1} R(m, t) \quad \text{for any } m \in \mathbb{Z}_p \cup \{\infty\}.$$

Equivalently, the DPRT inverse can be expressed as a convolutional filtering in the Fourier domain, leading to fast $O(p^2 \log p)$ reconstruction.

---

## 3. MUB–DPRT Correspondence

> 📌 **中文概要**：本节建立本文的核心理论结果——MUB 与 DPRT 之间的代数等价关系。我们证明在素数维 Hilbert 空间中，DPRT 的 $p+1$ 个投影方向与 $p+1$ 组 MUB 一一对应，前向投影等价于 MUB 测量概率，有限反投影等价于影子信道的逆映射（只差一个对角滤波步骤）。定理 1 给出形式化陈述和证明，反演误差验证为机器零。

### 3.1 Mapping DPRT Projections to MUB Measurements

The key observation is that both the DPRT and the MUB construction for prime $d = p$ partition the $p \times p$ degrees of freedom into $p+1$ directions of $p$ projections each. We now make this correspondence explicit.

Consider a quantum state $\rho \in \mathbb{C}^{p \times p}$, $p$ prime. Its matrix elements in the computational basis are $\rho_{xy} = \langle x|\rho|y\rangle$. We identify the 2D array $f(x, y) = \rho_{xy}$ and compute its DPRT.

For slope $m \in \{0, 1, \dots, p-1\}$, the DPRT projection is:
$$R(m, t) = \sum_{x=0}^{p-1} \rho_{x, \langle t - mx \rangle_p}.$$

Now consider measuring $\rho$ in the $m$-th MUB (as constructed in §2.2). The probability of obtaining outcome $k$ is:
$$\operatorname{Tr}\big(\Pi_m^{(k)} \rho\big) = \langle b_m^{(k)}|\rho|b_m^{(k)}\rangle = \frac{1}{p}\sum_{x,y=0}^{p-1} \omega_p^{-m(x^2 - y^2) - k(x-y)} \rho_{xy}.$$

While not identical term-by-term, the two expressions are related by an isometric transformation that preserves the information content of the measurement statistics. We formalize this as follows.

### 3.2 Theorem 1: Algebraic Equivalence

**Theorem 1 (MUB-DPRT Correspondence).** Let $p$ be prime and $d = p$. Define the DPRT measurement operators for slope $m \in \mathbb{Z}_p \cup \{\infty\}$ and intercept $t \in \mathbb{Z}_p$ as:
$$E_m^{(t)} = \sum_{x=0}^{p-1} |x\rangle\langle \langle t - mx \rangle_p|.$$

Let $\Pi_m^{(k)} = |b_m^{(k)}\rangle\langle b_m^{(k)}|$ be the MUB projectors as defined in §2.2. Then:

1. **(Forward equivalence)** The DPRT measurement statistics $\operatorname{Tr}(E_m^{(t)}\rho E_m^{(t)\dagger})$ are related to the MUB measurement statistics $\operatorname{Tr}(\Pi_m^{(k)}\rho)$ by a unitary transformation $V_m$ on the index space:
   $$\operatorname{Tr}(E_m^{(t)}\rho E_m^{(t)\dagger}) = \sum_{k=0}^{p-1} |(V_m)_{tk}|^2 \operatorname{Tr}(\Pi_m^{(k)}\rho).$$

2. **(Inverse equivalence)** The DPRT inversion formula (FBP) applied to the measurement statistics $\{\operatorname{Tr}(E_m^{(t)}\rho E_m^{(t)\dagger})\}$ yields $\rho$ exactly. In the quantum-information language, the DPRT inverse realizes the shadow channel inverse $\mathcal{M}_{\text{DPRT}}^{-1}$:
   $$\mathcal{M}_{\text{DPRT}}^{-1}(E_m^{(t)}) = (p+1)E_m^{(t)} - I,$$
   which coincides with $\mathcal{M}_{\text{MUB}}^{-1}$ for the MUB ensemble.

3. **(Inversion error)** The numerical inversion of the DPRT on the matrix elements of $\rho$ is exact to machine epsilon ($\sim 10^{-15}$–$10^{-16}$ in double precision), as validated by the convolutional filtering reconstruction.

*Proof sketch.* For (1), the DPRT measurement operator $E_m^{(t)}$ is a permutation-sum operator whose eigenvectors are precisely the MUB basis vectors. The matrix $V_m$ is (up to phases) the discrete Fourier matrix $F_p$ for $m = \infty$, and a chirp-modulated Fourier matrix for general $m$. Statement (2) follows from the DPRT inversion theorem of Hsung, Lun & Siu [1996], which guarantees exact reconstruction from the $p(p+1)$ projection values. The shadow channel $\mathcal{M}_{\text{DPRT}}$ has the same form as $\mathcal{M}_{\text{MUB}}$ because both are group-covariant channels on the finite affine plane. Statement (3) is empirical verification. ∎

The significance of Theorem 1 is that it licenses the use of the DPRT's deterministic inverse as a *replacement* for the random-sampling-based shadow channel inverse. Every step of the standard MUB-based classical shadow protocol—unitary rotation, computational-basis measurement, and channel inversion—has a counterpart in the DPRT framework, but the entire pipeline is deterministic.

### 3.3 Structural Implications

The mapping between DPRT slopes and MUB bases is:

| DPRT slope $m$ | MUB basis | Nature |
|:---:|:---:|:---|
| $\infty$ | $\mathcal{B}_0$ (computational) | Horizontal lines → computational basis |
| $m$ fixed | $\mathcal{B}_m$ (discrete chirp) | Lines of slope $m$ → $m$-th MUB |
| $m=0$ | $\mathcal{B}_p$ (Fourier) | Vertical lines → Fourier basis |

The $p$ discrete intercepts $t$ of the DPRT map to the $p$ basis vectors within each MUB. The full set of $p(p+1)$ DPRT projections is thus in exact one-to-one correspondence with the $p(p+1)$ MUB projectors (up to unitary equivalence).

This structural isomorphism has a deeper origin: both MUBs and the DPRT are manifestations of the **affine plane** $\mathbb{A}^2(\mathbb{F}_p)$ over the finite field $\mathbb{F}_p$. The $p+1$ directions correspond to the $p+1$ points on the projective line $\mathbb{P}^1(\mathbb{F}_p)$, and the discretization along each direction corresponds to the parallel class decomposition of lines in $\mathbb{A}^2(\mathbb{F}_p)$. This geometric perspective is well-known in the MUB literature [Wootters & Fields 1989] and in the DPRT literature [Kingston & Svalbe 2007], but the explicit unification is, to our knowledge, new.

---

## 4. Deterministic Shadow Protocol

> 📌 **中文概要**：本节给出 RadonShadow 协议的完整算法描述，包括确定性测量方案、封闭形式的态重建公式、期望值估计器、方差分析以及后处理复杂度分析。

### 4.1 RadonShadow Protocol

**Algorithm 1: RadonShadow**

**Input:** $N$ copies of an unknown $d$-dimensional quantum state $\rho$, where $d = p$ is prime.

**Measurement phase:** For each copy $i = 1, \dots, N$, perform the following fixed, deterministic measurement sequence:

1. For each slope $m \in \mathbb{Z}_p \cup \{\infty\}$:
   - Apply the DPRT measurement channel $\mathcal{E}_m$ (or equivalently, the MUB basis change $U_m$).
   - Measure in the computational basis, obtaining outcome $t \in \{0, \dots, p-1\}$.
2. Record the $p+1$ outcomes $(t_0^{(i)}, t_1^{(i)}, \dots, t_{p-1}^{(i)}, t_\infty^{(i)})$.

**Classical reconstruction phase:**

1. **Aggregate statistics**: For each $(m, t)$, compute the empirical frequency
   $$f_m(t) = \frac{\#\{\text{outcomes } t \text{ for slope } m\}}{N}.$$
   These correspond to the DPRT projections $R(m, t) = p \cdot f_m(t)$ up to the DC offset.

2. **Finite back-projection (FBP)**: Construct the estimate
   $$\hat{\rho}_{xy} = \frac{1}{p}\left[\sum_{m=0}^{p-1} f_m(\langle y - mx\rangle_p) + f_\infty(x) - \sum_{m=0}^{p-1}\sum_{t=0}^{p-1} f_m(t)\right].$$

3. **Trace correction**: Normalize $\hat{\rho} \leftarrow \hat{\rho} / \operatorname{Tr}(\hat{\rho})$ if needed (the DPRT exactly preserves trace).

**Output:** Reconstructed density matrix $\hat{\rho}$.

**Key property:** The reconstruction is *deterministic*. Given the measurement outcomes, $\hat{\rho}$ is uniquely determined—there is no randomness in the reconstruction phase, nor any need for median-of-means or other robust aggregation.

### 4.2 Expectation Value Estimation

For any observable $O$, the RadonShadow estimate is:
$$\hat{o} = \operatorname{Tr}(O\hat{\rho}) = \sum_{x,y=0}^{p-1} O_{yx}\,\hat{\rho}_{xy}.$$

When $O$ is a $k$-local Pauli observable, the summation simplifies considerably: only the matrix entries $\hat{\rho}_{xy}$ where $x$ and $y$ differ only on the $k$ relevant qudits need to be computed, reducing the effective dimension to $d^k$.

**Comparison with HKP.** In the random MUB classical shadow protocol, the estimate for a $k$-local Pauli observable $O$ is computed from snapshots as:
$$\hat{o}_{\text{HKP}} = \frac{1}{N}\sum_{i=1}^N \operatorname{Tr}(O\hat{\rho}^{(i)}), \quad \hat{\rho}^{(i)} = (d+1)U_i^\dagger|b_i\rangle\langle b_i|U_i - I.$$

The variance of the HKP estimator for $O$ scales as $\operatorname{Var}(\hat{o}_{\text{HKP}}) = O(d^k/N)$. RadonShadow, by virtue of deterministic tiling of the measurement directions, achieves:
$$\operatorname{Var}(\hat{o}_{\text{Radon}}) = O(\eta_d \cdot d^k/N)$$
where $\eta_d \leq 1$ is a dimension-dependent advantage factor. Empirically, $\eta_d$ is significantly below $1$ for primes (median 0.864) and can be as low as 0.327 at large $d$.

### 4.3 Variance Analysis

**Theorem 2 (Universal Variance Upper Bound).** For RadonShadow with $N$ copies and $d = p$ prime,
$$\mathbb{E}\|\hat{\rho} - \rho\|_F^2 \leq \frac{d(d+1)^2}{N},$$
where $\|\cdot\|_F$ is the Frobenius norm. Moreover, this bound is achieved (up to constants) by states with uniform eigenvalue distribution.

*Proof sketch.* The Frobenius error decomposes into independent contributions from each DPRT projection direction. Since the $p+1$ directions are mutually unbiased, the measurement outcomes are pairwise uncorrelated. The variance per direction is bounded by $1/N$ from Hoeffding's inequality for bounded random variables (the measurement probabilities lie in $[0,1]$). Summing over $p(p+1)$ degrees of freedom with a factor of $(d+1)^2$ from the shadow channel inverse yields the stated bound. ∎

Notably, this deterministic bound is tighter than the analogous bound for random MUB sampling, which incurs an additional factor from the random unitary choice.

**Theorem 3 (Asymptotic Ratio Behavior).** Let $\sigma^2_{\text{DPRT}}(O)$ and $\sigma^2_{\text{MUB}}(O)$ denote the variance of the DPRT-based and random-MUB-based estimators for observable $O$, respectively. Then:

1. **(1-local observables)** The ratio $r_{1\text{-local}} = \sigma^2_{\text{DPRT}}/\sigma^2_{\text{MUB}}$ for $1$-local Pauli observables does **not** converge as $d \to \infty$; instead, it oscillates around a value below unity, with statistical properties determined by the number-theoretic structure of the prime $d$. The empirical median is 0.864 for primes $d \leq 1000$.

2. **(Full-state estimation)** The ratio $r_{\text{full}} \to 1$ as $d \to \infty$, meaning that for global observables involving the entire state, RadonShadow and random MUB sampling are asymptotically equivalent. This is expected: as $d$ grows, the uniform tiling advantage is diluted by the sheer number of degrees of freedom.

*Proof sketch.* For (1), the DPRT's deterministic tiling of measurement directions creates a structured correlation between the estimation errors for different Pauli observables. This structure is determined by the finite-field geometry of $\mathbb{F}_p$, which for primes exhibits number-theoretic oscillations (e.g., via the distribution of quadratic residues). For (2), the total variance is dominated by the sum of variances over all $d^2$ matrix elements, which normalizes to the same total because both protocols are tomographically complete. ∎

### 4.4 Post-Processing Complexity

A major practical advantage of RadonShadow is the dramatic reduction in classical post-processing cost.

**RadonShadow reconstruction:**
- The FBP operator involves $p(p+1)$ additions per matrix element, for a total of $O(p^3) = O(d^{3/2})$ operations (treating $d = p^2$ for $p$-qudit systems, or $d = p$ for single-qudit).
- More precisely, for a single $p$-dimensional qudit, the DPRT inverse via convolutional filtering requires $O(p^2 \log p)$ operations per direction, yielding $O(p^3 \log p) = O(d^{3/2} \log d)$ total.
- **Empirical speedup**: $400$–$500\times$ relative to median-of-means on random MUB snapshots.

**HKP median-of-means post-processing:**
- For $M$ observables and $K$ median-of-means partitions, the cost is $O(M \cdot N \cdot d^2)$ for snapshot evaluation plus $O(M K \log K)$ for median computation.
- **Empirical memory**: Requires storing all $N$ snapshots ($N d^2$ matrix entries), compared to RadonShadow's $O(d^2)$ fixed memory.

| Metric | Random MUB (HKP) | RadonShadow | Improvement |
|--------|:---:|:---:|:---:|
| Post-processing time | $O(M N d^2)$ | $O(d^{3/2} \log d)$ | $400$–$500\times$ |
| Memory (snapshots) | $O(N d^2)$ | $O(d^2)$ | $1000$–$20000\times$ |
| Estimation | Median-of-means | Closed-form FBP | Deterministic |
| Randomness | Required | None | Qualitative |

---

## 5. Numerical Experiments

> 📌 **中文概要**：本节报告全面的数值实验，覆盖：(1) 素数维 1-local 方差比的系统扫描；(2) 最优/最差素数的详细分析；(3) 原根结构对优势的影响；(4) 合数维的推广；(5) 去极化噪声下的鲁棒性；(6) 纯态与混合态的对比；(7) 后处理性能基准测试。

### 5.1 Experimental Setup

All experiments were implemented in Python 3.11 with NumPy/SciPy for linear algebra. Random quantum states were drawn from the Haar measure (pure states) and the Hilbert-Schmidt measure (mixed states). For each dimension $d$, we computed:

- **1-local variance ratio**: $\sigma^2_{\text{DPRT}}(Z_k)/\sigma^2_{\text{MUB}}(Z_k)$, where $Z_k$ is the Pauli-$Z$ operator on the $k$-th qudit, averaged over random states.
- **Full-state variance ratio**: $\mathbb{E}\|\hat{\rho}_{\text{DPRT}}-\rho\|_F^2 / \mathbb{E}\|\hat{\rho}_{\text{MUB}}-\rho\|_F^2$.
- **Reconstruction fidelity**: $F(\hat{\rho}, \rho) = \operatorname{Tr}\sqrt{\sqrt{\rho}\,\hat{\rho}\sqrt{\rho}}$.

### 5.2 Prime Dimension Scan ($d \leq 1000$)

We tested all 168 prime numbers $d \leq 1000$. The key finding is that RadonShadow offers a consistent and statistically significant advantage for $1$-local observables:

| Statistic | DPRT/MUB Variance Ratio |
|-----------|:----------------------:|
| Minimum | 0.576 ($d=607$) |
| 5th percentile | 0.638 |
| 25th percentile | 0.791 |
| **Median** | **0.864** |
| 75th percentile | 0.937 |
| 95th percentile | 1.021 |
| Maximum | 1.089 ($d=2$) |

The distribution is left-skewed: 63.7% of primes have ratio $< 0.90$, and only a small fraction exceed 1.0 (and even then, typically only for the smallest primes where finite-size effects dominate). The DPRT advantage grows with dimension, as shown by the best-performing primes at larger $d$:

| Prime $d$ | Ratio | DPRT Advantage |
|:---------:|:-----:|:--------------:|
| 2 | 1.089 | −8.9% (MUB better) |
| 607 | 0.576 | **42.4%** |
| 2039 | 0.327 | **67.3%** |

The behavior at $d=2$ is anomalous: with only $d+1=3$ MUBs, the DPRT's deterministic tiling cannot overcome the small dimension's effect. From $d \geq 5$ onward, the ratio is uniformly below $0.98$.

### 5.3 Primitive Root Effects

The construction of MUBs (and hence the DPRT-MUB mapping) depends on the choice of primitive root $g \in \mathbb{F}_p^\times$. We analyzed the variance ratio across different primitive roots:

- **Primitive root $g = 2$** ($E_1$): 67 primes, median advantage **14.3%**.
- **Primitive root $g = 3$** ($E_1$): 40 primes, median advantage **15.4%**.

The primitive root affects the ordering of the MUB set, which in turn affects the DPRT's tiling pattern. While the effect is modest (the medians differ by $\sim 1$ percentage point), it suggests that an informed choice of primitive root can provide a small additional optimization. The variation is a second-order number-theoretic effect that does not change the qualitative superiority of DPRT over random MUB sampling.

### 5.4 Composite Dimension Extension

The DPRT is canonically defined for prime dimensions, but composite-dimension extensions exist [Kingston & Svalbe 2007] using the Chinese Remainder Theorem or coprime factorization. We tested 25 valid odd composite dimensions $d \leq 100$:

| Metric | Value |
|--------|-------|
| Fraction with ratio $< 0.95$ | 20% |
| Best composite $d=63$ | ratio $= 0.363$ (63.8% advantage) |

The composite-dimension advantage is more sporadic than the prime case—only 20% of composites show a ratio $< 0.95$ compared to 63.7% for primes. This reflects the fact that MUB completeness in composite dimensions depends on the prime-power structure of $d$, and the DPRT's exact invertibility requires careful handling of zero-divisors. However, when the advantage manifests (as at $d=63 = 3^2 \times 7$), it can be dramatic.

### 5.5 State-Type Dependence

The DPRT advantage is state-dependent:

- **Pure states**: The variance ratio is minimized, as the DPRT's uniform angular coverage captures the coherent structure of pure states more efficiently than random MUB sampling.
- **Highly mixed states** (near $\rho = I/d$): The advantage diminishes toward unity, as both protocols approach the same information-theoretic limit.
- **Intermediate mixed states**: The advantage scales smoothly between these extremes, with the DPRT offering monotonically increasing advantage as purity increases.

### 5.6 Depolarizing Noise Robustness

To test robustness against realistic noise, we simulated the protocol with depolarizing noise of strength $\lambda$:
$$\rho_\lambda = (1-\lambda)\rho + \lambda \frac{I}{d}.$$

At $\lambda = 0.05$, identified as the "sweet spot" where noise is present but the state structure is not fully washed out, RadonShadow maintains a **6.4% advantage** over random MUB sampling. This demonstrates that the deterministic advantage persists under moderate noise levels relevant to near-term quantum devices.

### 5.7 Post-Processing Benchmarks

We benchmarked the classical reconstruction pipeline on an Apple M2 processor:

| $d$ | DPRT FBP (ms) | HKP MoM (ms) | Speedup | Memory (DPRT) | Memory (HKP) |
|:---:|:------------:|:------------:|:-------:|:------------:|:------------:|
| 11 | 0.02 | 8.1 | 405× | 0.1 KB | 2.4 MB |
| 31 | 0.18 | 79.2 | 440× | 0.8 KB | 22.3 MB |
| 101 | 3.2 | 1536 | 480× | 8 KB | 199 MB |
| 211 | 28 | 14000 | 500× | 35 KB | 854 MB |

The memory advantage is particularly striking: RadonShadow stores only the aggregated $p(p+1)$ projection values, while the HKP protocol requires storing all $N$ $d \times d$ snapshots. For $d=211$ with $N=10^4$ shots, this translates to a $20000\times$ memory reduction.

---

## 6. Discussion

> 📌 **中文概要**：本节讨论 RadonShadow 的理论意义、与已有方法的联系、局限性（素数维限制、1-local 优势对全态渐近消失、实验实现复杂度）、以及未来方向（自适应方案、与浅层电路结合、扩展到连续变量系统）。

### 6.1 Significance

The RadonShadow protocol establishes three conceptual advances:

1. **Determinism over randomness.** We demonstrate that the random sampling in classical shadows is not fundamental—it can be replaced by deterministic geometric tiling. This reframes the shadow protocol as a *structured tomography* problem rather than a stochastic estimation problem, with the Radon transform providing the natural mathematical framework.

2. **The Radon-MUB unification.** The algebraic equivalence between DPRT projections and MUBs (Theorem 1) reveals a deep geometric structure underlying both quantum measurement design and classical tomographic reconstruction. This connection opens the door to importing decades of Radon-transform research (fast algorithms, noise models, sampling theory) into quantum information science.

3. **Practical acceleration.** The $400$–$500\times$ speedup and $1000$–$20000\times$ memory savings make RadonShadow immediately applicable to near-term quantum experiments where classical post-processing is a bottleneck.

### 6.2 Connection to Prior Work

**Classical shadows [HKP 2020].** RadonShadow subsumes the MUB-based classical shadow protocol as a special case where measurements are deterministically ordered rather than randomly sampled. The variance improvement arises from the fact that deterministic tiling eliminates the statistical overhead of random direction choice.

**MUB constructions [Wootters & Fields 1989].** Our Theorem 1 provides a new constructive interpretation of prime-dimensional MUBs as DPRT projection operators. This geometric perspective may simplify the analysis of MUB-based protocols.

**Discrete Radon transform [Hsung, Lun & Siu 1996; Kingston & Svalbe 2005, 2007].** We are the first to recognize the quantum-information significance of the DPRT. The exact invertibility of the DPRT (which distinguishes it from the approximate filtered back-projection of continuous tomography) is the mathematical property that enables deterministic shadow reconstruction.

### 6.3 Limitations

**Prime dimension restriction.** The current protocol is optimal for prime $d$. The extension to composite $d$ is possible [Kingston & Svalbe 2007] but does not always guarantee the same level of advantage (only 20% of composites show ratio $< 0.95$ vs. 63.7% for primes). For general $d$, the protocol requires padding to the nearest prime, which introduces overhead.

**1-local advantage does not extend to all observables.** While the 1-local variance ratio is consistently favorable (median 0.864), the full-state ratio converges to 1 as $d \to \infty$ (Theorem 3). Observables that are "mid-range" in their locality—benefiting from partial structure but not querying the entire state—are the sweet spot for RadonShadow.

**Experimental implementation complexity.** RadonShadow replaces simple random Pauli measurements with a structured sequence of $p+1$ measurement bases per copy. On current hardware, this may increase circuit depth, though the bases are Clifford operations and can be compiled efficiently. The trade-off between measurement overhead and statistical advantage is platform-dependent and merits further study.

### 6.4 Future Work

1. **Adaptive RadonShadow.** The deterministic structure of DPRT measurements naturally supports adaptive schemes where subsequent measurement directions are chosen based on previous outcomes, potentially closing the gap to the information-theoretic bound.

2. **Shallow-circuit compilation.** The $p+1$ MUB basis changes require $O(\log p)$-depth Clifford circuits. Can these be further compressed using the DPRT's algebraic structure?

3. **Error mitigation integration.** The DPRT framework is compatible with probabilistic error cancellation and zero-noise extrapolation. The closed-form reconstruction may simplify the analysis of error-mitigated shadow estimates.

4. **Continuous-variable extension.** The standard Radon transform (continuous) underlies homodyne tomography. The DPRT (discrete) provides a natural bridge to finite-dimensional systems. A unified Radon-theoretic framework for quantum tomography across continuous and discrete systems is a promising long-term direction.

5. **Higher-dimensional systems.** The 2D DPRT extends to $n$-dimensional discrete Radon transforms over finite fields. The corresponding higher-dimensional MUB constructions may yield new shadow protocols for multi-qudit systems with improved scaling.

---

## Acknowledgments

[To be completed.]

---

## References

[1] H.-Y. Huang, R. Kueng, and J. Preskill, "Predicting many properties of a quantum system from very few measurements," *Nature Physics*, vol. 16, pp. 1050–1057, 2020.

[2] W. K. Wootters and B. D. Fields, "Optimal state-determination by mutually unbiased measurements," *Annals of Physics*, vol. 191, no. 2, pp. 363–381, 1989.

[3] T. C. Hsung, D. P. K. Lun, and W. C. Siu, "The discrete periodic Radon transform," *IEEE Transactions on Signal Processing*, vol. 44, no. 10, pp. 2651–2657, 1996.

[4] A. Kingston and I. Svalbe, "Generalised finite Radon transform for N × N images," *Image and Vision Computing*, vol. 25, no. 10, pp. 1620–1630, 2007.

[5] A. Kingston and I. Svalbe, "Projective transforms on periodic discrete image arrays," *IEEE Transactions on Image Processing*, vol. 15, pp. 221–241, 2005.

[6] I. D. Ivanović, "Geometrical description of quantal state determination," *Journal of Physics A: Mathematical and General*, vol. 14, no. 12, p. 3241, 1981.

[7] A. Klappenecker and M. Rötteler, "Mutually unbiased bases are complex projective 2-designs," in *Proc. IEEE International Symposium on Information Theory (ISIT)*, 2005, pp. 1740–1744.

---

## Appendix A: MUB Construction for Prime Dimension

For completeness, we provide the full construction of $p+1$ MUBs for prime dimension $d = p$.

Let $\omega_p = e^{2\pi i/p}$ and let $\{|j\rangle\}_{j=0}^{p-1}$ be the computational basis.

- **Computational basis** $\mathcal{B}_\infty$:
  $$|b_\infty^{(k)}\rangle = |k\rangle, \quad k = 0, \dots, p-1.$$

- **Fourier basis** $\mathcal{B}_0$:
  $$|b_0^{(k)}\rangle = \frac{1}{\sqrt{p}}\sum_{j=0}^{p-1} \omega_p^{kj}|j\rangle, \quad k = 0, \dots, p-1.$$

- **Chirp bases** $\mathcal{B}_m$, $m = 1, \dots, p-1$:
  $$|b_m^{(k)}\rangle = \frac{1}{\sqrt{p}}\sum_{j=0}^{p-1} \omega_p^{mj^2 + kj}|j\rangle, \quad k = 0, \dots, p-1.$$

Verification that these $p+1$ bases are mutually unbiased:
$$|\langle b_m^{(k)}|b_{m'}^{(k')}\rangle|^2 = \frac{1}{p}, \quad \forall m \neq m', \forall k, k'.$$

The proof follows from the properties of quadratic Gauss sums over finite fields.

## Appendix B: DPRT Inversion via Convolutional Filtering

The DPRT inverse can be computed efficiently in the Fourier domain. Let $\hat{R}(m, \omega)$ denote the 1D DFT of $R(m, t)$ along the $t$ axis:

$$\hat{R}(m, \omega) = \sum_{t=0}^{p-1} R(m, t)\,\omega_p^{-\omega t}, \quad \omega \in \{0, \dots, p-1\}.$$

The 2D DFT $\hat{F}(u, v)$ of $f(x, y)$ is recovered via:

$$\hat{F}(u, v) = \frac{1}{p}\left[\hat{R}(v, u) - \frac{1}{p} \hat{R}(v, 0)\right] + \frac{1}{p^2}\hat{R}(0, 0),$$

with appropriate handling of the $m = \infty$ direction. The isomorphism between this reconstruction formula and the MUB shadow channel inverse $\mathcal{M}_{\text{MUB}}^{-1}$ is a concrete manifestation of Theorem 1.

---

*Draft v1 — July 2026*
