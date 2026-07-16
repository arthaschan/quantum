#!/usr/bin/env python3
"""
完整 Classical Shadow 对比实验
================================
对比三种方案:
1. DPRT 确定性遍历（论文原方案）
2. 简化 CS（论文当前基线，纯取平均）
3. 完整 CS（HKP 2020 标准方案: MoM + Physical Projection）

验证论文的核心 question: 完整 CS 是否会缩小或消除 DPRT 的优势？
"""

import numpy as np
import json, os, time, sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'validate'))
from core import build_mubs_wf, random_pure_state, primes_upto

# ============================================================
# 完整 CS 后处理
# ============================================================

def median_of_means(snapshots, K=10):
    """
    Median-of-Means (MoM) 估计。
    将 N 个快照随机分成 K 组，每组取均值，返回 K 个均值的 element-wise 中位数。
    """
    N = len(snapshots)
    # 至少需要 K 个快照才能分 K 组
    K = min(K, N)
    group_size = N // K
    
    # 随机打乱
    idx = np.random.permutation(N)
    
    group_means = []
    for k in range(K):
        start = k * group_size
        end = start + group_size if k < K - 1 else N
        if end <= start:  # 空组，跳过
            continue
        group_snaps = [snapshots[i] for i in idx[start:end]]
        group_mean = np.mean(group_snaps, axis=0)
        group_means.append(group_mean)
    
    if not group_means:
        return np.mean(snapshots, axis=0)  # fallback: simple mean
    
    # Element-wise 中位数
    stacked = np.stack(group_means)
    return np.median(np.real(stacked), axis=0) + 1j * np.median(np.imag(stacked), axis=0)


def physical_projection(rho, tol=1e-12):
    """
    Physical Projection: 将任意矩阵投影到合法密度矩阵空间。
    1. 对角化 → 负特征值截断为 0
    2. 重新归一化 trace=1
    """
    d = rho.shape[0]
    # 确保 Hermitian
    rho = (rho + rho.conj().T) / 2
    
    eigvals, eigvecs = np.linalg.eigh(rho)
    # 截断负特征值
    eigvals = np.maximum(eigvals, 0)
    # 重新归一化
    s = np.sum(eigvals)
    if s > tol:
        eigvals /= s
    
    return eigvecs @ np.diag(eigvals) @ eigvecs.conj().T


def classical_shadow_full(rho, mubs, n_total, K=None, seed=42):
    """
    完整 Classical Shadow (HKP 2020 标准方案)。
    K=None 时自适应：K = max(2, int(sqrt(n_total)))
    """
    d = rho.shape[0]
    n_bases = len(mubs)
    rng = np.random.RandomState(seed)
    
    if K is None:
        K = max(2, int(np.sqrt(n_total)))
    
    # 收集影子快照
    snaps = []
    for _ in range(n_total):
        m = rng.randint(n_bases)
        basis = mubs[m]
        probs = np.array([np.real(np.conj(v) @ rho @ v) for v in basis])
        probs = np.clip(probs, 0, None); probs /= probs.sum()
        b = rng.choice(d, p=probs)
        vec = basis[b]
        snaps.append(np.outer(vec, vec.conj()))
    
    # 构建影子矩阵（使用 MUB shadow 公式）
    shadow_snaps = [n_bases * s - np.eye(d) for s in snaps]
    
    # MoM
    rho_mom = median_of_means(shadow_snaps, K=K)
    
    # Physical Projection
    rho_proj = physical_projection(rho_mom)
    
    return rho_proj


def classical_shadow_simple(rho, mubs, n_total, seed=42):
    """简化 CS（论文原基线）：简单取平均，无 MoM/Projection。"""
    d = rho.shape[0]
    n_bases = len(mubs)
    rng = np.random.RandomState(seed)
    
    snaps = []
    for _ in range(n_total):
        m = rng.randint(n_bases)
        basis = mubs[m]
        probs = np.array([np.real(np.conj(v) @ rho @ v) for v in basis])
        probs = np.clip(probs, 0, None); probs /= probs.sum()
        b = rng.choice(d, p=probs)
        vec = basis[b]
        snaps.append(np.outer(vec, vec.conj()))
    
    return np.mean(snaps, axis=0) * n_bases - np.eye(d)


def deterministic_dprt(rho, mubs, r_per_dir=4):
    """DPRT 确定性遍历（论文方案）。"""
    d = rho.shape[0]
    n_bases = len(mubs)
    
    snaps = []
    rng = np.random.RandomState(42)
    for m in range(n_bases):
        basis = mubs[m]
        probs = np.array([np.real(np.conj(v) @ rho @ v) for v in basis])
        probs = np.clip(probs, 0, None); probs /= probs.sum()
        for _ in range(r_per_dir):
            b = rng.choice(d, p=probs)
            vec = basis[b]
            snaps.append(np.outer(vec, vec.conj()))
    
    return np.mean(snaps, axis=0) * n_bases - np.eye(d)


# ============================================================
# 实验运行
# ============================================================

def run_comparison(test_dims, n_states=5, r_per_dir=2, K_mom=10, n_obs=8):
    """对每个维度，比较三种方案。"""
    results = []
    
    for d in test_dims:
        print(f"  d={d:>4d}...", end=" ", flush=True)
        t0 = time.time()
        
        mubs = build_mubs_wf(d)
        n_total = len(mubs) * r_per_dir
        
        ratios_simple = []
        ratios_full = []
        ratios_dprt_vs_full = []
        
        for si in range(n_states):
            psi = random_pure_state(d, seed=d * 1000 + si)
            rho = np.outer(psi, psi.conj())
            
            # 构造 1-local 可观测量
            obs_list = []
            count = 0
            for a in range(d):
                if count >= n_obs: break
                for b in range(a + 1, d):
                    if count >= n_obs: break
                    O = np.zeros((d, d), dtype=complex)
                    O[a, b] = 1.0; O[b, a] = 1.0
                    obs_list.append(O)
                    count += 1
            
            truth = np.array([np.real(np.trace(O @ rho)) for O in obs_list])
            
            # DPRT 确定性
            rho_d = deterministic_dprt(rho, mubs, r_per_dir=r_per_dir)
            est_d = np.array([np.real(np.trace(O @ rho_d)) for O in obs_list])
            mae_d = np.mean(np.abs(est_d - truth))
            
            # 简化 CS
            rho_cs = classical_shadow_simple(rho, mubs, n_total, seed=d*1000+si+1)
            est_cs = np.array([np.real(np.trace(O @ rho_cs)) for O in obs_list])
            mae_cs = np.mean(np.abs(est_cs - truth))
            
            # 完整 CS (MoM + Projection)
            rho_csf = classical_shadow_full(rho, mubs, n_total, K=K_mom, seed=d*1000+si+2)
            est_csf = np.array([np.real(np.trace(O @ rho_csf)) for O in obs_list])
            mae_csf = np.mean(np.abs(est_csf - truth))
            
            if mae_cs > 0:
                ratios_simple.append(mae_d / mae_cs)
            if mae_csf > 0:
                ratios_full.append(mae_d / mae_csf)
            ratios_dprt_vs_full.append(mae_csf / mae_cs if mae_cs > 0 else 1.0)
        
        dt = time.time() - t0
        r = {
            'd': d,
            'ratio_dprt_vs_simple_cs': float(np.median(ratios_simple)),
            'ratio_dprt_vs_full_cs': float(np.median(ratios_full)),
            'ratio_full_cs_vs_simple_cs': float(np.median(ratios_dprt_vs_full)),
            'dprt_adv_over_simple_pct': float((1 - np.median(ratios_simple)) * 100),
            'dprt_adv_over_full_pct': float((1 - np.median(ratios_full)) * 100),
            'n_states': n_states,
            'time_s': float(dt)
        }
        
        marker = "★" if r['ratio_dprt_vs_full_cs'] < 0.95 else ""
        print(f"DPRT/simple={r['ratio_dprt_vs_simple_cs']:.3f}  "
              f"DPRT/full={r['ratio_dprt_vs_full_cs']:.3f}  "
              f"full/simple={r['ratio_full_cs_vs_simple_cs']:.3f}  "
              f"[{dt:.1f}s] {marker}", flush=True)
        
        results.append(r)
    
    return results


if __name__ == '__main__':
    # 代表性素数
    test_dims = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--fast', action='store_true', help='Only test small dims')
    parser.add_argument('--states', type=int, default=5)
    parser.add_argument('--K', type=int, default=10, help='MoM groups')
    args = parser.parse_args()
    
    if args.fast:
        test_dims = [3, 5, 7, 11, 13, 17, 19]
    
    print("=" * 75)
    print(f"FULL CS COMPARISON: {len(test_dims)} dims × {args.states} states")
    print(f"  MoM groups K={args.K}")
    print("=" * 75)
    print(f"  {'d':>4s}  {'DPRT/simple':>11s}  {'DPRT/full':>11s}  {'full/simple':>12s}")
    print(f"  {'-'*4}  {'-'*11}  {'-'*11}  {'-'*12}")
    
    T0 = time.time()
    all_results = run_comparison(test_dims, n_states=args.states, K_mom=args.K)
    total_time = time.time() - T0
    
    # 统计分析
    ratios_s = [r['ratio_dprt_vs_simple_cs'] for r in all_results]
    ratios_f = [r['ratio_dprt_vs_full_cs'] for r in all_results]
    
    print(f"\n{'='*75}")
    print(f"SUMMARY ({len(test_dims)} dims, {total_time:.0f}s)")
    print(f"{'='*75}")
    print(f"  DPRT vs Simple CS:")
    print(f"    median ratio: {np.median(ratios_s):.4f}  (advantage: {(1-np.median(ratios_s))*100:.1f}%)")
    print(f"    mean ratio:   {np.mean(ratios_s):.4f}")
    print(f"    ratio < 1.00: {sum(1 for r in ratios_s if r<1.00)}/{len(ratios_s)}")
    print(f"    ratio < 0.90: {sum(1 for r in ratios_s if r<0.90)}/{len(ratios_s)}")
    print()
    print(f"  DPRT vs Full CS (MoM + Projection):")
    print(f"    median ratio: {np.median(ratios_f):.4f}  (advantage: {(1-np.median(ratios_f))*100:.1f}%)")
    print(f"    mean ratio:   {np.mean(ratios_f):.4f}")
    print(f"    ratio < 1.00: {sum(1 for r in ratios_f if r<1.00)}/{len(ratios_f)}")
    print(f"    ratio < 0.90: {sum(1 for r in ratios_f if r<0.90)}/{len(ratios_f)}")
    
    # 关键对比
    delta = np.median(ratios_f) - np.median(ratios_s)
    print(f"\n  Impact of using Full CS:")
    print(f"    Simple CS median ratio: {np.median(ratios_s):.4f}")
    print(f"    Full CS median ratio:   {np.median(ratios_f):.4f}")
    print(f"    Change:                 {delta:+.4f} (DPRT advantage {'increased' if delta<0 else 'decreased'} by {abs(delta)*100:.1f} pp)")
    
    # 保存
    out = {
        'config': {'test_dims': test_dims, 'n_states': args.states, 'K_mom': args.K, 'r_per_dir': 2},
        'results': all_results,
        'summary': {
            'dprt_vs_simple_median': float(np.median(ratios_s)),
            'dprt_vs_full_median': float(np.median(ratios_f)),
            'advantage_delta_pp': float(delta * 100),
            'dprt_wins_over_simple': sum(1 for r in ratios_s if r<1.0),
            'dprt_wins_over_full': sum(1 for r in ratios_f if r<1.0),
        }
    }
    
    out_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'full_cs_comparison.json')
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    
    print(f"\n→ {out_path}")
