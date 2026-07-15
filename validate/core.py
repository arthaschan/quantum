#!/usr/bin/env python3
"""
RadonShadow 独立验证 — 核心模块
================================
包含 MUB 构造、DPRT、Wigner 函数、态生成等基础功能。
完全独立于原始 scripts/ 目录，避免循环验证。

遵循论文中描述的公式和协议设计。
"""

import numpy as np
from numpy.linalg import norm
import math

# ============================================================
# 1. 数论工具
# ============================================================

def primes_upto(n: int) -> list:
    """Eratosthenes 筛法，返回所有 ≤ n 的素数"""
    sieve = bytearray(b'\x01') * (n + 1)
    sieve[:2] = b'\x00\x00'
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i*i:n+1:i] = b'\x00' * ((n - i*i) // i + 1)
    return [i for i in range(2, n + 1) if sieve[i]]

def prime_factors(n: int) -> list:
    """返回 n 的所有质因子（去重）"""
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            if d not in factors:
                factors.append(d)
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factors.append(n)
    return factors

def primitive_roots(p: int) -> list:
    """返回素数 p 的所有原根"""
    if p <= 1:
        return []
    phi = p - 1
    factors = prime_factors(phi)
    roots = []
    for g in range(2, p):
        is_root = True
        for q in factors:
            if pow(g, phi // q, p) == 1:
                is_root = False
                break
        if is_root:
            roots.append(g)
    return roots

def smallest_primitive_root(p: int) -> int:
    """返回最小的原根"""
    roots = primitive_roots(p)
    return roots[0] if roots else None

# ============================================================
# 2. 量子态生成
# ============================================================

def random_pure_state(d: int, seed: int = None) -> np.ndarray:
    """生成随机纯态 |ψ⟩ ∈ C^d，从复正态分布采样后归一化"""
    rng = np.random.RandomState(seed) if seed is not None else np.random
    psi = rng.randn(d) + 1j * rng.randn(d)
    return psi / norm(psi)

def random_mixed_state(d: int, seed: int = None, rank: int = None) -> np.ndarray:
    """生成随机混合态 ρ，使用 Ginibre 系综方法"""
    rng = np.random.RandomState(seed) if seed is not None else np.random
    if rank is None:
        rank = d
    G = rng.randn(d, rank) + 1j * rng.randn(d, rank)
    rho = G @ G.conj().T
    return rho / np.trace(rho)

def depolarized_state(psi: np.ndarray, lam: float) -> np.ndarray:
    """对纯态 |ψ⟩ 施加去极化噪声: ρ_λ = (1-λ)|ψ⟩⟨ψ| + λ I/d"""
    d = len(psi)
    rho_pure = np.outer(psi, psi.conj())
    return (1 - lam) * rho_pure + lam * np.eye(d) / d

def maximally_mixed_state(d: int) -> np.ndarray:
    """最大混合态 I_d / d"""
    return np.eye(d) / d

# ============================================================
# 3. MUB 构造（Wootters-Fields 1989 构造）
# ============================================================

def build_mubs_wf(d: int) -> list:
    """
    为素数 d=p 构造完备的 p+1 组 MUB。
    
    使用 Wootters-Fields 1989 构造:
    - 第 0 组: 计算基 {|0⟩, |1⟩, ..., |d-1⟩}
    - 第 a 组 (a=1..d-1): |e^(a)_b⟩ = (1/√d) Σ_j ω^{a·j(j-1)/2 + b·j} |j⟩
    - 第 d 组 (∞ 组): |e^(∞)_b⟩ = (1/√d) Σ_j ω^{b·j} |j⟩ (DFT 基)
    
    返回: list of (list of np.ndarray) — 共 d+1 组，每组 d 个 d 维向量。
    """
    omega = np.exp(2j * np.pi / d)
    mubs = []
    
    # 第 0 组: 计算基
    mubs.append([np.eye(d, dtype=complex)[k] for k in range(d)])
    
    # 第 a 组 (a = 1, 2, ..., d-1)
    for a in range(1, d):
        basis = []
        for b in range(d):
            vec = np.array([omega**(a * j * (j-1) // 2 + b * j) for j in range(d)], dtype=complex)
            basis.append(vec / np.sqrt(d))
        mubs.append(basis)
    
    # 第 d 组 (∞ 组): DFT 基
    basis_inf = []
    for b in range(d):
        vec = np.array([omega**(b * j) for j in range(d)], dtype=complex)
        basis_inf.append(vec / np.sqrt(d))
    mubs.append(basis_inf)
    
    return mubs


def verify_mub_property(mubs: list, d: int, tol: float = 1e-12) -> dict:
    """
    验证 MUB 性质:
    1. 每组基是正交归一的
    2. 不同组之间满足互无偏: |⟨e_i^(a)|e_j^(b)⟩| = 1/√d
    
    返回验证结果字典。
    """
    n_bases = len(mubs)
    expected = 1.0 / np.sqrt(d)
    results = {
        'n_bases': n_bases,
        'd': d,
        'orthonormal_ok': True,
        'mutually_unbiased_ok': True,
        'max_inner_product_error': 0.0,
        'inner_products': []
    }
    
    # 检查每组基的正交归一性
    for i in range(n_bases):
        for a in range(d):
            for b in range(d):
                ip = np.abs(np.vdot(mubs[i][a], mubs[i][b]))
                if a == b:
                    if abs(ip - 1.0) > tol:
                        results['orthonormal_ok'] = False
                else:
                    if abs(ip) > tol:
                        results['orthonormal_ok'] = False
    
    # 检查互无偏性
    for i in range(n_bases):
        for j in range(i + 1, n_bases):
            for a in range(d):
                for b in range(d):
                    ip = np.abs(np.vdot(mubs[i][a], mubs[j][b]))
                    error = abs(ip - expected)
                    results['max_inner_product_error'] = max(results['max_inner_product_error'], error)
                    if error > tol:
                        results['inner_products'].append({
                            'basis_i': i, 'basis_j': j,
                            'vec_a': a, 'vec_b': b,
                            'inner_product': float(ip),
                            'error': float(error)
                        })
                        results['mutually_unbiased_ok'] = False
    
    return results


# ============================================================
# 4. MUB 投影概率
# ============================================================

def mub_projection_probabilities(rho: np.ndarray, mubs: list) -> np.ndarray:
    """
    计算密度矩阵 ρ 在完备 MUB 上的投影概率。
    
    返回: (d+1) × d 矩阵 P[a][b] = ⟨e^(a)_b|ρ|e^(a)_b⟩
    """
    d = rho.shape[0]
    n_bases = len(mubs)
    probs = np.zeros((n_bases, d))
    for a in range(n_bases):
        for b in range(d):
            vec = mubs[a][b]
            probs[a, b] = np.real(np.conj(vec) @ rho @ vec)
    # 确保概率归一化
    probs = np.clip(probs, 0, None)
    row_sums = probs.sum(axis=1)
    for a in range(n_bases):
        if row_sums[a] > 0:
            probs[a] /= row_sums[a]
    return probs


# ============================================================
# 5. 离散 Wigner 函数
# ============================================================

def discrete_wigner(psi: np.ndarray) -> np.ndarray:
    """
    计算纯态 |ψ⟩ 在 F_p×F_p 上的离散 Wigner 函数。
    
    W_ψ(x, y) = (1/p) Σ_j c_j c_{j+x}^* ω^{-y(j + x/2)}
    
    返回: p×p 实矩阵
    """
    d = len(psi)
    omega = np.exp(2j * np.pi / d)
    W = np.zeros((d, d), dtype=complex)
    
    for x in range(d):
        for y in range(d):
            s = 0.0
            for j in range(d):
                jx = (j + x) % d
                # 使用整数运算避免浮点模运算
                phase = -y * (j + x/2)
                s += psi[j] * np.conj(psi[jx]) * omega**phase
            W[x, y] = s / d
    
    return np.real(W)  # Wigner 函数是实值的


# ============================================================
# 6. DPRT（离散投影 Radon 变换）
# ============================================================

def dprt_forward(f: np.ndarray) -> np.ndarray:
    """
    对 p×p 网格函数 f 执行 DPRT。
    
    方向:
    - m = 0, 1, ..., p-1: L(m, t) = {(x, y): y ≡ mx + t (mod p)}
    - m = p: L(p, t) = {(t, y): y = 0, ..., p-1} (垂直方向)
    
    返回: (p+1) × p 矩阵 R[m][t]
    """
    p = f.shape[0]
    R = np.zeros((p + 1, p))
    
    # 方向 m = 0, 1, ..., p-1
    for m in range(p):
        for t in range(p):
            for x in range(p):
                y = (m * x + t) % p
                R[m, t] += f[x, y]
    
    # 方向 m = p (垂直)
    for t in range(p):
        for y in range(p):
            R[p, t] += f[t, y]
    
    return R


def dprt_inverse(R: np.ndarray) -> np.ndarray:
    """
    DPRT 反演: 从投影数据 R 恢复原始 p×p 函数 f。
    
    反演公式 (layperson_guide §1.6):
    f(x,y) = (1/p) [ Σ_{m=0}^{p-1} R(m, (y-mx) mod p) + R(p, x) - Σ_{k=0}^{p-1} R(0, k) ]
    """
    p = R.shape[1]  # R is (p+1) × p
    f = np.zeros((p, p))
    
    # 总投影和（用于归一化）
    total = np.sum(R[0, :])
    
    for x in range(p):
        for y in range(p):
            s = 0.0
            for m in range(p):
                t = (y - m * x) % p
                s += R[m, t]
            s += R[p, x] - total
            f[x, y] = s / p
    
    return f


# ============================================================
# 7. Classical Shadow 估计
# ============================================================

def classical_shadow_snapshot(rho: np.ndarray, mubs: list, m_idx: int = None,
                               rng: np.random.RandomState = None) -> np.ndarray:
    """
    单次 Classical Shadow 快照（MUB 变体）。
    
    若 m_idx 给定，使用确定性方向；否则随机选择方向。
    
    ρ̂_snap = (d+1) |e^(a)_b⟩⟨e^(a)_b| - I_d
    
    其中 b 按投影概率 ∝ ⟨e^(a)_b|ρ|e^(a)_b⟩ 采样。
    """
    d = rho.shape[0]
    n_bases = len(mubs)
    
    if rng is None:
        rng = np.random
    
    if m_idx is None:
        m_idx = rng.randint(n_bases)
    
    # 在当前 MUB 上的投影概率
    probs = np.array([np.real(np.conj(mubs[m_idx][b]) @ rho @ mubs[m_idx][b]) for b in range(d)])
    probs = np.clip(probs, 0, None)
    probs /= probs.sum()
    
    # 按概率采样结果
    b = rng.choice(d, p=probs)
    
    # 构造影子快照
    vec = mubs[m_idx][b]
    snap = n_bases * np.outer(vec, vec.conj()) - np.eye(d)
    return snap


def estimate_1local_expectations(rho: np.ndarray, mubs: list,
                                   shots_per_dir: int = 1,
                                   delta: int = 15) -> dict:
    """
    确定性遍历方案（DPRT-like）的 1-local 可观测量估计。
    
    遍历全部 d+1 个 MUB 方向，每个方向测量 shots_per_dir 次。
    
    返回: {
        'truth': 真实期望值列表,
        'estimates': 估计值列表,
        'mae': MAE,
        'n_total': 总测量次数
    }
    """
    d = rho.shape[0]
    n_bases = len(mubs)
    n_total = n_bases * shots_per_dir
    
    # 构造 1-local 可观测量列表
    obs_list = []
    count = 0
    for a in range(d):
        if count >= delta: break
        for b in range(a + 1, d):
            if count >= delta: break
            obs_list.append((a, b))
            count += 1
    
    # 真实值
    truth = []
    for a, b in obs_list:
        truth.append(2 * np.real(rho[a, b]))  # tr((|a⟩⟨b|+|b⟩⟨a|) ρ) = 2Re(ρ_ab)
    
    # 确定性收集影子快照
    snaps = []
    rng = np.random.RandomState(42)
    for m in range(n_bases):
        for _ in range(shots_per_dir):
            snap = classical_shadow_snapshot(rho, mubs, m, rng)
            snaps.append(snap)
    
    # 估计
    estimates = []
    for a, b in obs_list:
        O = np.zeros((d, d), dtype=complex)
        O[a, b] = 1.0
        O[b, a] = 1.0  # |a⟩⟨b| + |b⟩⟨a|
        est = np.mean([np.real(np.trace(O @ s)) for s in snaps])
        estimates.append(est)
    
    mae = np.mean(np.abs(np.array(estimates) - np.array(truth)))
    
    return {
        'truth': truth,
        'estimates': estimates,
        'mae': mae,
        'n_total': n_total
    }


def estimate_1local_random(rho: np.ndarray, mubs: list,
                             n_total: int,
                             delta: int = 15,
                             seed: int = 42) -> dict:
    """
    随机 Classical Shadow 方案的 1-local 可观测量估计。
    
    随机选择 MUB 方向 n_total 次。
    """
    d = rho.shape[0]
    n_bases = len(mubs)
    
    # 构造 1-local 可观测量列表
    obs_list = []
    count = 0
    for a in range(d):
        if count >= delta: break
        for b in range(a + 1, d):
            if count >= delta: break
            obs_list.append((a, b))
            count += 1
    
    # 真实值
    truth = []
    for a, b in obs_list:
        truth.append(2 * np.real(rho[a, b]))
    
    # 随机收集影子快照
    snaps = []
    rng = np.random.RandomState(seed)
    for _ in range(n_total):
        snap = classical_shadow_snapshot(rho, mubs, None, rng)
        snaps.append(snap)
    
    # 估计
    estimates = []
    for a, b in obs_list:
        O = np.zeros((d, d), dtype=complex)
        O[a, b] = 1.0
        O[b, a] = 1.0
        est = np.mean([np.real(np.trace(O @ s)) for s in snaps])
        estimates.append(est)
    
    mae = np.mean(np.abs(np.array(estimates) - np.array(truth)))
    
    return {
        'truth': truth,
        'estimates': estimates,
        'mae': mae,
        'n_total': n_total
    }


# ============================================================
# 8. 噪声通道
# ============================================================

def apply_noise(rho: np.ndarray, channel: str, lam: float) -> np.ndarray:
    """
    对密度矩阵 ρ 施加噪声通道。
    
    通道:
    - 'depol': 去极化 ε(ρ) = (1-λ)ρ + λ I/d
    - 'amp': 振幅阻尼 (简化版，d=2)
    - 'phase': 相位阻尼 / 纯退相干
    - 'bitflip': Pauli-X 噪声 (广义)
    """
    d = rho.shape[0]
    
    if channel == 'depol':
        return (1 - lam) * rho + lam * np.eye(d) / d
    
    elif channel == 'phase':
        # 相位阻尼: 非对角元衰减
        rho_out = rho.copy()
        for i in range(d):
            for j in range(d):
                if i != j:
                    rho_out[i, j] *= (1 - lam)
        return rho_out
    
    elif channel == 'bitflip':
        # 广义比特翻转: 以概率 λ 应用 X 算符
        # X|j⟩ = |j+1 mod d⟩
        rho_out = (1 - lam) * rho.copy()
        for k in range(1, d):
            Xk = np.roll(np.eye(d), k, axis=1)
            rho_out += (lam / (d - 1)) * Xk @ rho @ Xk.conj().T
        return rho_out
    
    elif channel == 'amp':
        # 振幅阻尼 (仅支持 d=2 的精确实现)
        if d == 2:
            E0 = np.array([[1, 0], [0, np.sqrt(1 - lam)]], dtype=complex)
            E1 = np.array([[0, np.sqrt(lam)], [0, 0]], dtype=complex)
            return E0 @ rho @ E0.conj().T + E1 @ rho @ E1.conj().T
        else:
            # 广义振幅阻尼近似
            rho_out = (1 - lam) * rho.copy()
            rho_out[0, 0] += lam * np.trace(rho)
            return rho_out / np.trace(rho_out)
    
    else:
        raise ValueError(f"Unknown channel: {channel}")


# ============================================================
# 9. Spearman 相关系数
# ============================================================

def spearman_rho(x: list, y: list) -> float:
    """计算 Spearman 秩相关系数"""
    from scipy.stats import spearmanr
    return spearmanr(x, y)[0]

def pearson_r(x: list, y: list) -> float:
    """计算 Pearson 相关系数"""
    from scipy.stats import pearsonr
    return pearsonr(x, y)[0]


# ============================================================
# 10. 统计工具
# ============================================================

def bootstrap_median_ci(data: list, n_bootstrap: int = 10000, ci: float = 0.95) -> tuple:
    """Bootstrap 中位数置信区间"""
    medians = []
    n = len(data)
    rng = np.random.RandomState(42)
    for _ in range(n_bootstrap):
        sample = rng.choice(data, size=n, replace=True)
        medians.append(np.median(sample))
    alpha = (1 - ci) / 2
    lower = np.percentile(medians, 100 * alpha)
    upper = np.percentile(medians, 100 * (1 - alpha))
    return np.median(medians), lower, upper
