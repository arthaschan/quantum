#!/usr/bin/env python3
"""
RadonShadow: Noise Robustness + Adaptive Direction Selection
=============================================================
1A: 多噪声模型鲁棒性 — Depolarizing / AmpDamping / PhaseDamping / BitFlip
1B: 自适应方向选择 — 逐步加 MUB 方向直到重建收敛

Output: noise_adaptive_report.md + noise_adaptive_data.json
"""

import numpy as np
import time, json, sys
from collections import defaultdict

# ═══════════════════════════════════════════════════════════════
# 0. MUB 构造 (Wootters & Fields 1989, 素数 d)
# ═══════════════════════════════════════════════════════════════

def build_mubs(d):
    """构造素数 d 的 d+1 组 MUB."""
    omega = np.exp(2j * np.pi / d)
    mubs = []
    # M0: 计算基
    mubs.append([np.eye(d, dtype=complex)[j] for j in range(d)])
    # M1...Md: Fourier + diagonal twist
    for m in range(1, d + 1):
        basis = []
        for k in range(d):
            v = np.array([
                omega ** (m * j * (j - 1) // 2 + k * j)
                for j in range(d)
            ], dtype=complex)
            basis.append(v / np.sqrt(d))
        mubs.append(basis)
    return mubs


def project_to_basis(rho, basis):
    """ρ → MUB 投影概率分布."""
    probs = np.array([np.real(np.vdot(v, rho @ v)) for v in basis])
    probs = np.clip(probs, 0, None)
    s = probs.sum()
    return probs / s if s > 0 else np.ones(len(basis)) / len(basis)


def snap_from_probs(basis, outcomes, n_per):
    """从测量结果构造阴影估计."""
    d = len(basis)
    est = np.zeros((d, d), dtype=complex)
    for k in outcomes:
        v = basis[k]
        est += np.outer(v, v.conj())
    return est * (d + 1) / n_per - np.eye(d)


def shadow_invert(det_shadow_total, n_total):
    """MUB 阴影反演: ρ̂ = (d+1)/N * Σ - I."""
    d_local = det_shadow_total.shape[0]
    return det_shadow_total * (d_local + 1) / n_total - np.eye(d_local)


# ═══════════════════════════════════════════════════════════════
# 1. 多噪声通道 (qudit 通用)
# ═══════════════════════════════════════════════════════════════

def apply_noise(rho, channel, strength):
    """对密度矩阵施加噪声通道. channel ∈ {depol, amp, phase, bitflip}."""
    d = rho.shape[0]

    if channel == 'depol':
        # ρ → (1-λ)ρ + λ·I/d
        return (1 - strength) * rho + strength * np.eye(d, dtype=complex) / d

    elif channel == 'amp':
        # 振幅阻尼: 末态 → 基态 |0⟩
        gamma = strength
        # Kraus 算符 (qudit 推广)
        # K₀ = |0⟩⟨0| + √(1-γ)·Σ_{j=1}^{d-1} |j⟩⟨j|
        K0 = np.diag([1.0] + [np.sqrt(1 - gamma)] * (d - 1))
        # K_j = √γ·|0⟩⟨j|  for j=1..d-1  (每个激发态到基态的衰减通道)
        result = K0 @ rho @ K0.conj().T
        for j in range(1, d):
            K = np.zeros((d, d), dtype=complex)
            K[0, j] = np.sqrt(gamma / (d - 1))  # 均匀分配衰减
            result += K @ rho @ K.conj().T
        return result

    elif channel == 'phase':
        # 相位阻尼: 对角元不变，非对角元衰减
        gamma = strength
        result = (1 - gamma) * rho.copy()
        for i in range(d):
            result[i, i] = rho[i, i].real  # 保持对角元
        return result

    elif channel == 'bitflip':
        # 广义比特翻转: X_d |j⟩ = |j+1 mod d⟩
        prob = strength
        Xd = np.zeros((d, d), dtype=complex)
        for j in range(d):
            Xd[(j + 1) % d, j] = 1.0
        return (1 - prob) * rho + prob * (Xd @ rho @ Xd.conj().T)

    else:
        raise ValueError(f"Unknown channel: {channel}")


# ═══════════════════════════════════════════════════════════════
# 2. 态生成
# ═══════════════════════════════════════════════════════════════

def random_pure_state(d, seed):
    rng = np.random.RandomState(seed)
    z = rng.randn(d) + 1j * rng.randn(d)
    return np.outer(z, z.conj()) / np.linalg.norm(z)**2


def random_mixed_state(d, rank, seed):
    rng = np.random.RandomState(seed)
    evals = rng.dirichlet(np.ones(rank))
    if len(evals) < d:
        evals = np.pad(evals, (0, d - len(evals)), constant_values=0)
    # random unitary via QR
    A = rng.randn(d, d) + 1j * rng.randn(d, d)
    Q, _ = np.linalg.qr(A)
    return Q @ np.diag(evals) @ Q.conj().T


# ═══════════════════════════════════════════════════════════════
# 3. 1-local 可观测量
# ═══════════════════════════════════════════════════════════════

def get_1local_observables(d):
    """生成 1-local 可观测量: Re(ρ_ij) for i<j. 取前 min(8, C(d,2)) 个."""
    obs = []
    for i in range(d):
        for j in range(i + 1, d):
            O = np.zeros((d, d), dtype=complex)
            O[i, j] = 1.0
            O[j, i] = 1.0
            obs.append(O)
    return obs[:min(8, len(obs))]


def estimate_obs(rho_est, observables, truth):
    """计算 MAE 的 ratio."""
    mae = np.mean([abs(np.real(np.trace(O @ rho_est)) - t)
                   for O, t in zip(observables, truth)])
    return mae


# ═══════════════════════════════════════════════════════════════
# 4. 实验 1A: 多噪声模型鲁棒性
# ═══════════════════════════════════════════════════════════════

def expt_1a_noise(test_dims, noise_channels, noise_levels,
                  r_val=4, n_states=8):
    """
    参数:
      test_dims: 被测素数维列表
      noise_channels: ['depol','amp','phase','bitflip']
      noise_levels: [0.0, 0.01, 0.05, 0.1, 0.2]
      r_val: 每方向测量次数 (= N/d+1)
      n_states: 每条件的随机纯态数
    """
    print("=" * 70)
    print("EXPERIMENT 1A: MULTI-NOISE ROBUSTNESS")
    print("=" * 70)

    results = []
    seeds = {}

    for d in test_dims:
        mubs = build_mubs(d)
        n_mubs = len(mubs)          # d+1
        n_per = r_val
        n_total = n_mubs * r_val
        obs_1loc = get_1local_observables(d)
        print(f"\n  d={d}  MUBs={n_mubs}  N_total={n_total}")

        for channel in noise_channels:
            for noise in noise_levels:
                ratios = []
                for si in range(n_states):
                    seed = d * 100000 + hash(channel) % 10000 + int(noise * 1000) + si
                    rng = np.random.RandomState(seed)

                    # 生成纯态 + 加噪
                    z = rng.randn(d) + 1j * rng.randn(d)
                    psi = z / np.linalg.norm(z)
                    rho_pure = np.outer(psi, psi.conj())
                    rho = apply_noise(rho_pure, channel, noise)

                    # 真值
                    truth = [np.real(np.trace(O @ rho)) for O in obs_1loc]

                    # ── 确定性方案 ──
                    det_total = np.zeros((d, d), dtype=complex)
                    for mi in range(n_mubs):
                        probs = project_to_basis(rho, mubs[mi])
                        seed_d = seed + 1000 + mi
                        rng_d = np.random.RandomState(seed_d)
                        outcomes = rng_d.choice(d, n_per, p=probs)
                        for k in outcomes:
                            det_total += np.outer(mubs[mi][k], mubs[mi][k].conj())
                    det_rho = shadow_invert(det_total, n_total)

                    # ── 随机方案 ──
                    rand_total = np.zeros((d, d), dtype=complex)
                    seed_r = seed + 2000
                    rng_r = np.random.RandomState(seed_r)
                    for _ in range(n_total):
                        mi = rng_r.randint(n_mubs)
                        probs = project_to_basis(rho, mubs[mi])
                        k = rng_r.choice(d, p=probs)
                        rand_total += np.outer(mubs[mi][k], mubs[mi][k].conj())
                    rand_rho = shadow_invert(rand_total, n_total)

                    mae_d = estimate_obs(det_rho, obs_1loc, truth)
                    mae_r = estimate_obs(rand_rho, obs_1loc, truth)
                    if mae_r > 1e-12:
                        ratios.append(mae_d / mae_r)

                if ratios:
                    mean_r = float(np.mean(ratios))
                    std_r = float(np.std(ratios))
                    adv = (1 - mean_r) * 100
                    marker = " ★" if mean_r < 0.85 else (" ✦" if mean_r < 0.90 else "")
                    print(f"    {channel:>8s} λ={noise:.2f}  "
                          f"ratio={mean_r:.4f}±{std_r:.4f}  "
                          f"DPRT{'+' if adv>0 else ''}{adv:.1f}% {marker}")
                else:
                    mean_r, std_r, adv = 1.0, 0, 0
                    print(f"    {channel:>8s} λ={noise:.2f}  ratio=1.000 (all zero MAE)")

                results.append({
                    'd': d,
                    'channel': channel,
                    'noise': noise,
                    'ratio_mean': mean_r,
                    'ratio_std': std_r,
                    'advantage_pct': adv,
                    'n_samples': len(ratios)
                })

    return results


# ═══════════════════════════════════════════════════════════════
# 5. 实验 1B: 自适应方向选择
# ═══════════════════════════════════════════════════════════════

def expt_1b_adaptive(test_dims, epsilons, n_states=8, noise_lvl=0.0):
    """
    逐步加 MUB 方向直到 Frobenius 收敛.

    参数:
      epsilons: 收敛阈值 ε 列表
      noise_lvl: 背景去极化噪声水平 (0=纯态)
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 1B: ADAPTIVE DIRECTION SELECTION")
    print("=" * 70)

    results = []

    for d in test_dims:
        mubs = build_mubs(d)
        n_mubs = len(mubs)
        obs_1loc = get_1local_observables(d)
        print(f"\n  d={d}  MUBs={n_mubs}")

        for eps in epsilons:
            savings = []    # 1 - k_stop / (d+1)
            ratios_adapt = []
            frob_errors = []

            for si in range(n_states):
                seed = d * 200000 + int(eps * 10000) + si
                rng = np.random.RandomState(seed)

                # 态
                z = rng.randn(d) + 1j * rng.randn(d)
                psi = z / np.linalg.norm(z)
                rho_true = np.outer(psi, psi.conj())
                if noise_lvl > 0:
                    rho_true = apply_noise(rho_true, 'depol', noise_lvl)

                truth = [np.real(np.trace(O @ rho_true)) for O in obs_1loc]

                # ── 全遍历确定性 (baseline) ──
                n_per = 4
                full_total = np.zeros((d, d), dtype=complex)
                for mi in range(n_mubs):
                    probs = project_to_basis(rho_true, mubs[mi])
                    seed_f = seed + 3000 + mi
                    rng_f = np.random.RandomState(seed_f)
                    outcomes = rng_f.choice(d, n_per, p=probs)
                    for k in outcomes:
                        full_total += np.outer(mubs[mi][k], mubs[mi][k].conj())
                rho_full = shadow_invert(full_total, n_mubs * n_per)

                # ── 随机方案 baseline ──
                rand_total = np.zeros((d, d), dtype=complex)
                seed_r = seed + 4000
                rng_r = np.random.RandomState(seed_r)
                for _ in range(n_mubs * n_per):
                    mi = rng_r.randint(n_mubs)
                    probs = project_to_basis(rho_true, mubs[mi])
                    k = rng_r.choice(d, p=probs)
                    rand_total += np.outer(mubs[mi][k], mubs[mi][k].conj())
                rho_rand = shadow_invert(rand_total, n_mubs * n_per)

                # ── 自适应遍历 ──
                # 方向顺序: M0(计算基) → M1(Fourier) → M2... 按顺序加入
                adapt_total = np.zeros((d, d), dtype=complex)
                rho_prev = np.zeros((d, d), dtype=complex)
                k_stop = n_mubs  # 默认用完所有方向

                for mi in range(n_mubs):
                    probs = project_to_basis(rho_true, mubs[mi])
                    seed_a = seed + 5000 + mi
                    rng_a = np.random.RandomState(seed_a)
                    outcomes = rng_a.choice(d, n_per, p=probs)
                    for k in outcomes:
                        adapt_total += np.outer(mubs[mi][k], mubs[mi][k].conj())

                    n_used = (mi + 1) * n_per
                    rho_curr = shadow_invert(adapt_total, n_used)

                    if mi >= 1:  # 至少用 2 个方向后才能判断收敛
                        frob = np.linalg.norm(rho_curr - rho_prev, 'fro')
                        if frob < eps:
                            k_stop = mi + 1
                            break

                    rho_prev = rho_curr.copy()

                # 分析
                savings.append(1.0 - k_stop / n_mubs)

                mae_adapt = estimate_obs(rho_curr, obs_1loc, truth)
                mae_rand = estimate_obs(rho_rand, obs_1loc, truth)
                mae_full = estimate_obs(rho_full, obs_1loc, truth)

                if mae_rand > 1e-12:
                    ratios_adapt.append(mae_adapt / mae_rand)
                frob_errors.append(float(np.linalg.norm(rho_curr - rho_full, 'fro')))

            mean_save = float(np.mean(savings)) * 100
            mean_ratio = float(np.mean(ratios_adapt)) if ratios_adapt else 1.0
            mean_frob = float(np.mean(frob_errors))
            marker = " ★" if mean_save > 30 else (" ✦" if mean_save > 15 else "")
            print(f"    ε={eps:.4f}  save={mean_save:.1f}%  "
                  f"ratio={mean_ratio:.4f}  frob_gap={mean_frob:.6f} {marker}")

            results.append({
                'd': d,
                'epsilon': eps,
                'mean_savings_pct': mean_save,
                'mean_ratio': mean_ratio,
                'mean_frob_gap': mean_frob,
                'n_samples': n_states
            })

    return results


# ═══════════════════════════════════════════════════════════════
# 6. 组合分析: 噪声 × 自适应
# ═══════════════════════════════════════════════════════════════

def expt_1c_noise_adaptive(test_dims, channels_noise_pairs, epsilons,
                            n_states=6):
    """
    在噪声条件下运行自适应方案.
    channels_noise_pairs: [(ch, λ), ...]
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 1C: NOISE × ADAPTIVE INTERACTION")
    print("=" * 70)

    results = []

    for d in test_dims:
        mubs = build_mubs(d)
        n_mubs = len(mubs)
        obs_1loc = get_1local_observables(d)
        print(f"\n  d={d}")

        for ch, noise in channels_noise_pairs:
            for eps in epsilons:
                savings = []
                ratios_adapt = []

                for si in range(n_states):
                    seed = d * 300000 + hash(ch) % 5000 + int(noise * 1000) + int(eps * 10000) + si
                    rng = np.random.RandomState(seed)

                    z = rng.randn(d) + 1j * rng.randn(d)
                    psi = z / np.linalg.norm(z)
                    rho_pure = np.outer(psi, psi.conj())
                    rho_true = apply_noise(rho_pure, ch, noise)
                    truth = [np.real(np.trace(O @ rho_true)) for O in obs_1loc]

                    n_per = 4

                    # 随机 baseline
                    rand_total = np.zeros((d, d), dtype=complex)
                    seed_r = seed + 5000
                    rng_r = np.random.RandomState(seed_r)
                    for _ in range(n_mubs * n_per):
                        mi = rng_r.randint(n_mubs)
                        probs = project_to_basis(rho_true, mubs[mi])
                        k = rng_r.choice(d, p=probs)
                        rand_total += np.outer(mubs[mi][k], mubs[mi][k].conj())
                    rho_rand = shadow_invert(rand_total, n_mubs * n_per)

                    # 自适应
                    adapt_total = np.zeros((d, d), dtype=complex)
                    rho_prev = np.zeros((d, d), dtype=complex)
                    k_stop = n_mubs

                    for mi in range(n_mubs):
                        probs = project_to_basis(rho_true, mubs[mi])
                        seed_a = seed + 6000 + mi
                        rng_a = np.random.RandomState(seed_a)
                        outcomes = rng_a.choice(d, n_per, p=probs)
                        for k in outcomes:
                            adapt_total += np.outer(mubs[mi][k], mubs[mi][k].conj())

                        n_used = (mi + 1) * n_per
                        rho_curr = shadow_invert(adapt_total, n_used)

                        if mi >= 1:
                            if np.linalg.norm(rho_curr - rho_prev, 'fro') < eps:
                                k_stop = mi + 1
                                break
                        rho_prev = rho_curr.copy()

                    savings.append(1.0 - k_stop / n_mubs)
                    mae_a = estimate_obs(rho_curr, obs_1loc, truth)
                    mae_r = estimate_obs(rho_rand, obs_1loc, truth)
                    if mae_r > 1e-12:
                        ratios_adapt.append(mae_a / mae_r)

                mean_save = float(np.mean(savings)) * 100
                mean_ratio = float(np.mean(ratios_adapt)) if ratios_adapt else 1.0
                print(f"    {ch}(λ={noise:.2f}) ε={eps:.4f}  "
                      f"save={mean_save:.1f}%  ratio={mean_ratio:.4f}")

                results.append({
                    'd': d,
                    'channel': ch,
                    'noise': noise,
                    'epsilon': eps,
                    'mean_savings_pct': mean_save,
                    'mean_ratio': mean_ratio,
                    'n_samples': n_states
                })

    return results


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    T0 = time.time()

    TEST_DIMS = [3, 5, 7, 11, 13, 17, 19, 23]
    NOISE_CHANNELS = ['depol', 'amp', 'phase', 'bitflip']
    NOISE_LEVELS = [0.0, 0.01, 0.05, 0.10, 0.20]
    EPSILONS = [0.001, 0.005, 0.01, 0.02]

    results = {}

    # ── 1A: 噪声 ──
    print("\n🔬 Running 1A: Multi-Noise Robustness...")
    t1a = time.time()
    results['1a_noise'] = expt_1a_noise(TEST_DIMS, NOISE_CHANNELS, NOISE_LEVELS)
    print(f"\n  1A done in {time.time()-t1a:.0f}s  ({len(results['1a_noise'])} data points)")

    # ── 1B: 自适应 ──
    print("\n🔬 Running 1B: Adaptive Direction Selection...")
    t1b = time.time()
    results['1b_adaptive'] = expt_1b_adaptive(TEST_DIMS, EPSILONS)
    print(f"\n  1B done in {time.time()-t1b:.0f}s  ({len(results['1b_adaptive'])} data points)")

    # ── 1C: 噪声×自适应 ──
    print("\n🔬 Running 1C: Noise × Adaptive Interaction...")
    t1c = time.time()
    # 选代表性的噪声条件: depol 甜区 + 各通道中等强度
    noise_pairs = [
        ('depol', 0.0), ('depol', 0.05), ('depol', 0.10),
        ('amp', 0.05), ('phase', 0.05), ('bitflip', 0.05),
    ]
    results['1c_noise_adaptive'] = expt_1c_noise_adaptive(TEST_DIMS, noise_pairs, EPSILONS)
    print(f"\n  1C done in {time.time()-t1c:.0f}s  ({len(results['1c_noise_adaptive'])} data points)")

    elapsed = time.time() - T0

    # ═══════════════════════════════════════════════════════════
    # 汇总
    # ═══════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print(f"ALL EXPERIMENTS COMPLETE — {elapsed:.0f}s total")
    print("=" * 70)

    # 1A summary by channel
    print("\n── 1A: Noise Summary by Channel ──")
    for ch in NOISE_CHANNELS:
        print(f"\n  {ch.upper()}:")
        for noise in NOISE_LEVELS:
            pts = [r for r in results['1a_noise']
                   if r['channel'] == ch and abs(r['noise'] - noise) < 1e-9]
            if pts:
                ratios = [p['ratio_mean'] for p in pts]
                print(f"    λ={noise:.2f}: "
                      f"mean_ratio={np.mean(ratios):.4f}  "
                      f"min={min(ratios):.4f}  "
                      f"max={max(ratios):.4f}  "
                      f"DPRT={'★' if np.mean(ratios)<0.90 else '✦' if np.mean(ratios)<0.95 else '='}")

    # 1B summary
    print("\n── 1B: Adaptive Summary ──")
    for eps in EPSILONS:
        pts = [r for r in results['1b_adaptive'] if abs(r['epsilon'] - eps) < 1e-9]
        if pts:
            saves = [p['mean_savings_pct'] for p in pts]
            ratios = [p['mean_ratio'] for p in pts]
            print(f"  ε={eps:.4f}: "
                  f"mean_save={np.mean(saves):.1f}%  "
                  f"mean_ratio={np.mean(ratios):.4f}  "
                  f"best_save={max(saves):.1f}%  "
                  f"best_d={max(pts, key=lambda x: x['mean_savings_pct'])['d']}")

    # 1C summary
    print("\n── 1C: Noise × Adaptive Summary ──")
    for ch, noise in sorted(set((r['channel'], r['noise']) for r in results['1c_noise_adaptive'])):
        for eps in EPSILONS:
            pts = [r for r in results['1c_noise_adaptive']
                   if r['channel'] == ch and abs(r['noise'] - noise) < 1e-9
                   and abs(r['epsilon'] - eps) < 1e-9]
            if pts and np.mean([p['mean_savings_pct'] for p in pts]) > 5:
                sav = np.mean([p['mean_savings_pct'] for p in pts])
                rat = np.mean([p['mean_ratio'] for p in pts])
                print(f"  {ch}(λ={noise:.2f}) ε={eps:.4f}: save={sav:.1f}%  ratio={rat:.4f}")

    # Save JSON
    out_json = 'data/noise_adaptive_results.json'
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n📁 Results saved to {out_json}")
    print(f"📁 Report: noise_adaptive_report.md (generated below)")
