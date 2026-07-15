#!/usr/bin/env python3
"""
验证 E3: 大素数扫描 + 乘法子群光滑度分析
===========================================
独立验证:
1. d=2039, 2851, 4289 的 ratio
2. Spearman ρ(p-1 质因子数, ratio) ≈ -0.8913

注意: 大素数计算非常耗时，使用采样模式。
"""

import numpy as np
import json
import sys
import os
import time

sys.path.insert(0, os.path.dirname(__file__))
from core import *

def run_e3_verification():
    results = {}
    
    print("=" * 70)
    print("E3 VERIFICATION: Large Primes + p-1 Smoothness")
    print("=" * 70)
    
    # ────────────────────────────────────────────────
    # Part 1: 关键大素数验证
    # ────────────────────────────────────────────────
    print("\n── Part 1: Key Large Primes ──")
    test_primes = [2039, 2851, 4289]
    # 同时加入中等素数作为对照
    medium_primes = [607, 941, 857, 829]
    all_test = test_primes + medium_primes
    
    key_results = {}
    
    for d in all_test:
        print(f"  d={d}...", end=" ", flush=True)
        t0 = time.time()
        mubs = build_mubs_wf(d)
        n_mubs = len(mubs)
        
        ratios = []
        for r in [1, 2]:
            n_total = n_mubs * r
            psi = random_pure_state(d, seed=d * 100)
            rho = np.outer(psi, psi.conj())
            
            det = estimate_1local_expectations(rho, mubs, shots_per_dir=r, delta=10)
            rand = estimate_1local_random(rho, mubs, n_total=n_total, delta=10, seed=d * 100 + 1)
            
            ratio = det['mae'] / rand['mae'] if rand['mae'] > 0 else 1.0
            ratios.append(ratio)
        
        best_ratio = min(ratios)
        dt = time.time() - t0
        key_results[d] = {
            'best_ratio': float(best_ratio),
            'advantage_pct': float((1 - best_ratio) * 100),
            'ratio_r1': float(ratios[0]),
            'ratio_r2': float(ratios[1]),
            'computation_time_s': float(dt)
        }
        print(f"ratio={best_ratio:.4f} ({(1-best_ratio)*100:.1f}% adv) [{dt:.0f}s]")
    
    results['key_primes'] = key_results
    
    # ────────────────────────────────────────────────
    # Part 2: Spearman 相关性分析
    # ────────────────────────────────────────────────
    print("\n── Part 2: p-1 Smoothness vs Ratio Correlation ──")
    
    # 我们需要 ratio 数据。如果验证模式是 sample，用 E1 数据
    e1_data_path = os.path.join(os.path.dirname(__file__), 'e1_verification.json')
    
    if os.path.exists(e1_data_path):
        with open(e1_data_path) as f:
            e1_data = json.load(f)
        
        ratios_dict = {item['d']: item['ratio'] for item in e1_data['top10']}
        for d_str, val in e1_data['all_details'].items():
            d = int(d_str)
            ratios_dict[d] = val['best_ratio']
    else:
        # 如果没有 E1 数据，快速计算一些素数
        print("  No E1 data found, computing ratios for correlation analysis...")
        ratios_dict = {}
        test_primes_corr = primes_upto(200)  # d ≤ 200 for speed
        for d in test_primes_corr:
            if d < 3: continue
            mubs = build_mubs_wf(d)
            n_mubs = len(mubs)
            psi = random_pure_state(d, seed=d * 100)
            rho = np.outer(psi, psi.conj())
            
            det = estimate_1local_expectations(rho, mubs, shots_per_dir=2, delta=8)
            rand = estimate_1local_random(rho, mubs, n_total=n_mubs * 2, delta=8, seed=d * 100 + 1)
            ratio = det['mae'] / rand['mae'] if rand['mae'] > 0 else 1.0
            ratios_dict[d] = ratio
    
    # 计算 p-1 质因子数与 ratio 的相关性
    pf_counts = []
    ratio_vals = []
    details = []
    
    for d in sorted(ratios_dict.keys()):
        if d <= 2: continue
        pf = prime_factors(d - 1)
        pf_counts.append(len(pf))
        ratio_vals.append(ratios_dict[d])
        details.append({
            'd': d,
            'p_minus_1': d - 1,
            'prime_factors': pf,
            'n_factors': len(pf),
            'ratio': float(ratios_dict[d])
        })
    
    # 计算相关系数
    try:
        from scipy.stats import spearmanr, pearsonr
        sp_rho, sp_p = spearmanr(pf_counts, ratio_vals)
        pr_r, pr_p = pearsonr(pf_counts, ratio_vals)
    except ImportError:
        # Fallback: 手动计算近似
        sp_rho = np.corrcoef(np.argsort(np.argsort(pf_counts)).astype(float),
                             np.argsort(np.argsort(ratio_vals)).astype(float))[0, 1]
        pr_r = np.corrcoef(pf_counts, ratio_vals)[0, 1]
        sp_p, pr_p = None, None
    
    correlation = {
        'n_primes_analyzed': len(details),
        'spearman_rho': float(sp_rho),
        'spearman_p_value': float(sp_p) if sp_p is not None else None,
        'pearson_r': float(pr_r),
        'pearson_p_value': float(pr_p) if pr_p is not None else None,
        'details': details[:20]  # 只保留前 20 个
    }
    
    results['correlation'] = correlation
    
    # 按质因子数分组
    group_stats = {}
    for d in details:
        nf = d['n_factors']
        if nf not in group_stats:
            group_stats[nf] = []
        group_stats[nf].append(d['ratio'])
    
    print(f"\n  Spearman ρ(unique_pf, ratio) = {sp_rho:.4f}  (p={sp_p if sp_p else 'N/A'})")
    print(f"  Pearson r(unique_pf, ratio) = {pr_r:.4f}  (p={pr_p if pr_p else 'N/A'})")
    print(f"\n  By number of unique prime factors of p-1:")
    
    for nf in sorted(group_stats.keys()):
        vals = group_stats[nf]
        print(f"    nf={nf}: n={len(vals)}, mean_ratio={np.mean(vals):.4f}, "
              f"median_ratio={np.median(vals):.4f}")
    
    group_summary = {}
    for nf in sorted(group_stats.keys()):
        vals = group_stats[nf]
        group_summary[f'nf={nf}'] = {
            'n': len(vals),
            'mean': float(np.mean(vals)),
            'median': float(np.median(vals)),
            'min': float(np.min(vals)),
            'max': float(np.max(vals))
        }
    results['group_by_factors'] = group_summary
    
    # ────────────────────────────────────────────────
    # Part 3: Safe Prime 分析
    # ────────────────────────────────────────────────
    print("\n── Part 3: Safe Prime Candidates ──")
    safe_primes = []
    for d in details:
        if d['n_factors'] == 2 and d['d'] > 5:
            # Check if p-1 = 2 × q where q is prime
            factors = d['prime_factors']
            if 2 in factors:
                other = [f for f in factors if f != 2]
                if len(other) == 1:
                    safe_primes.append(d)
    
    print(f"  Found {len(safe_primes)} Safe Prime candidates:")
    for sp in safe_primes[:10]:
        print(f"    p={sp['d']}, p-1=2×{sp['prime_factors'][-1]}, ratio={sp['ratio']:.4f}")
    
    results['safe_primes'] = safe_primes[:20]
    
    # ────────────────────────────────────────────────
    # 保存
    # ────────────────────────────────────────────────
    out_path = os.path.join(os.path.dirname(__file__), 'e3_verification.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n  Results saved to e3_verification.json")
    
    return results

if __name__ == '__main__':
    run_e3_verification()
