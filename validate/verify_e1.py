#!/usr/bin/env python3
"""
验证 E1: 素数扫描 d ≤ 1000 — 1-local 可观测量
=================================================
独立验证 168 个素数的 ratio 数据：
- 中位数 ratio ≈ 0.864
- ratio < 0.90 占比 ~63.7%
- Top-10 表现素数

与原始实验结果交叉对比。
"""

import numpy as np
import json
import sys
import os
import time

sys.path.insert(0, os.path.dirname(__file__))
from core import *

def run_e1_verification(max_d: int = 1000, n_states_per_d: int = 3,
                         r_list: list = [1, 2, 4], delta: int = 15,
                         verification_mode: str = 'sample'):
    """
    验证 E1 实验。
    
    参数:
    - max_d: 最大维度
    - n_states_per_d: 每个维度测试的随机态数量
    - r_list: 每个方向的测量次数列表
    - delta: 1-local 可观测量数量上限
    - verification_mode: 'sample' 快速抽样 或 'full' 全量 168 素数
    """
    
    if verification_mode == 'sample':
        # 快速模式：代表性素数抽样
        test_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                       53, 59, 61, 67, 71, 73, 79, 83, 89, 97,  # d ≤ 100
                       107, 409, 547, 607, 733, 829, 857, 887, 941, 997]  # Top primes
    else:
        test_primes = primes_upto(max_d)
    
    print("=" * 70)
    print(f"E1 VERIFICATION: Prime Scan d ≤ {max_d} — 1-local Observables")
    print(f"Mode: {verification_mode}, {len(test_primes)} primes")
    print("=" * 70)
    
    all_ratios = {}  # d -> best ratio
    all_details = {}  # d -> full details
    
    T0 = time.time()
    
    for i, d in enumerate(test_primes):
        # 报告进度
        if i % 20 == 0 and verification_mode == 'full':
            elapsed = time.time() - T0
            print(f"  Progress: {i}/{len(test_primes)} primes, {elapsed:.0f}s elapsed", flush=True)
        
        mubs = build_mubs_wf(d)
        n_mubs = len(mubs)
        
        ratios_1loc = []
        
        for r in r_list:
            n_total = n_mubs * r
            ratio_vals = []
            
            for si in range(n_states_per_d):
                seed = d * 10000 + r * 100 + si * 10
                psi = random_pure_state(d, seed=seed)
                rho = np.outer(psi, psi.conj())
                
                # 确定性估计
                det = estimate_1local_expectations(rho, mubs, shots_per_dir=r, delta=delta)
                
                # 随机估计
                rand = estimate_1local_random(rho, mubs, n_total=n_total, delta=delta,
                                              seed=seed + 1)
                
                ratio = det['mae'] / rand['mae'] if rand['mae'] > 0 else 1.0
                ratio_vals.append(ratio)
            
            ratios_1loc.append(np.mean(ratio_vals))
        
        best_ratio = min(ratios_1loc)
        all_ratios[d] = best_ratio
        all_details[d] = {
            'ratios_per_r': {f'r={r}': float(v) for r, v in zip(r_list, ratios_1loc)},
            'best_ratio': float(best_ratio),
            'advantage_pct': float((1 - best_ratio) * 100)
        }
    
    elapsed = time.time() - T0
    
    # ────────────────────────────────────────────────
    # 统计分析
    # ────────────────────────────────────────────────
    ratio_list = list(all_ratios.values())
    ratio_list.sort()
    
    n = len(ratio_list)
    median_ratio = ratio_list[n // 2]
    mean_ratio = np.mean(ratio_list)
    
    count_under_09 = sum(1 for r in ratio_list if r < 0.90)
    count_under_085 = sum(1 for r in ratio_list if r < 0.85)
    count_under_08 = sum(1 for r in ratio_list if r < 0.80)
    count_under_10 = sum(1 for r in ratio_list if r < 1.00)
    
    # Bootstrap 中位数 CI
    med_boot, med_lower, med_upper = bootstrap_median_ci(ratio_list)
    
    stats = {
        'n_primes': n,
        'median_ratio': float(median_ratio),
        'mean_ratio': float(mean_ratio),
        'median_advantage_pct': float((1 - median_ratio) * 100),
        'ratio_lt_090': {'count': count_under_09, 'pct': float(count_under_09 / n * 100)},
        'ratio_lt_085': {'count': count_under_085, 'pct': float(count_under_085 / n * 100)},
        'ratio_lt_080': {'count': count_under_08, 'pct': float(count_under_08 / n * 100)},
        'ratio_lt_100': {'count': count_under_10, 'pct': float(count_under_10 / n * 100)},
        'dprt_wins': count_under_10,
        'dprt_loses': n - count_under_10,
        'bootstrap_median_95ci': [float(med_lower), float(med_upper)],
        'computation_time_s': float(elapsed)
    }
    
    # Top-10
    ranked = sorted(all_ratios.items(), key=lambda x: x[1])
    top10 = [{'d': d, 'ratio': float(r), 'advantage_pct': float((1-r)*100)} 
             for d, r in ranked[:10]]
    
    # 按维度范围分组的统计
    range_stats = []
    ranges = [(0, 50), (50, 100), (100, 200), (200, 500), (500, 1000)]
    for lo, hi in ranges:
        in_range = [r for d, r in all_ratios.items() if lo <= d < hi]
        if in_range:
            range_stats.append({
                'range': f'[{lo},{hi})',
                'n': len(in_range),
                'median': float(np.median(in_range)),
                'mean': float(np.mean(in_range)),
                'min': float(np.min(in_range)),
                'max': float(np.max(in_range))
            })
    
    # ────────────────────────────────────────────────
    # 与原始数据对比（如果可用）
    # ────────────────────────────────────────────────
    comparison = None
    orig_data_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'primes_1000_results.json')
    if os.path.exists(orig_data_path):
        with open(orig_data_path) as f:
            orig_data = json.load(f)
        
        orig_ratios = {}
        for key, val in orig_data.items():
            d = val['d']
            best = min(r['ratio_1loc'] for r in val['runs'])
            orig_ratios[d] = best
        
        # 找交集
        common_dims = set(all_ratios.keys()) & set(orig_ratios.keys())
        diffs = []
        for d in sorted(common_dims):
            diff = all_ratios[d] - orig_ratios[d]
            diffs.append(abs(diff))
        
        comparison = {
            'n_compared': len(common_dims),
            'mean_abs_diff': float(np.mean(diffs)) if diffs else None,
            'max_abs_diff': float(np.max(diffs)) if diffs else None,
            'median_abs_diff': float(np.median(diffs)) if diffs else None,
            'correlation': float(np.corrcoef(
                [all_ratios[d] for d in sorted(common_dims)],
                [orig_ratios[d] for d in sorted(common_dims)]
            )[0, 1]) if len(common_dims) > 1 else None,
            'orig_median': float(np.median([orig_ratios[d] for d in sorted(common_dims)])),
            'our_median': float(np.median([all_ratios[d] for d in sorted(common_dims)])),
        }
    
    # ────────────────────────────────────────────────
    # 输出
    # ────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"E1 SUMMARY ({n} primes, {elapsed:.1f}s)")
    print(f"{'='*70}")
    print(f"  Median ratio:  {median_ratio:.4f}  (DPRT advantage: {(1-median_ratio)*100:.1f}%)")
    print(f"  Mean ratio:    {mean_ratio:.4f}")
    print(f"  Bootstrap 95% CI for median: [{med_lower:.4f}, {med_upper:.4f}]")
    print(f"  ratio < 0.90:  {count_under_09}/{n} ({count_under_09/n*100:.1f}%)")
    print(f"  ratio < 0.85:  {count_under_085}/{n} ({count_under_085/n*100:.1f}%)")
    print(f"  ratio < 0.80:  {count_under_08}/{n} ({count_under_08/n*100:.1f}%)")
    print(f"  DPRT wins:     {count_under_10}/{n} ({count_under_10/n*100:.1f}%)")
    
    print(f"\n  Top-10 Best Dimensions:")
    for item in top10:
        print(f"    d={item['d']:>4d}  ratio={item['ratio']:.4f}  advantage={item['advantage_pct']:.1f}%")
    
    print(f"\n  Ratio by d-Range:")
    for rs in range_stats:
        print(f"    d∈{rs['range']}: n={rs['n']:>3d}  median={rs['median']:.4f}  "
              f"mean={rs['mean']:.4f}  min={rs['min']:.4f}  max={rs['max']:.4f}")
    
    if comparison:
        print(f"\n  Comparison with Original Data ({comparison['n_compared']} common dims):")
        print(f"    Original median:  {comparison['orig_median']:.4f}")
        print(f"    Our median:       {comparison['our_median']:.4f}")
        print(f"    Mean absolute difference: {comparison['mean_abs_diff']:.5f}")
        print(f"    Max absolute difference:  {comparison['max_abs_diff']:.5f}")
        print(f"    Correlation:      {comparison['correlation']:.4f}")
    
    # ────────────────────────────────────────────────
    # 保存
    # ────────────────────────────────────────────────
    output = {
        'metadata': {
            'mode': verification_mode,
            'n_states_per_d': n_states_per_d,
            'r_list': r_list,
            'max_d': max_d,
            'delta': delta,
        },
        'statistics': stats,
        'top10': top10,
        'range_breakdown': range_stats,
        'comparison_with_original': comparison,
        'all_details': all_details
    }
    
    out_path = os.path.join(os.path.dirname(__file__), 'e1_verification.json')
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n  Results saved to e1_verification.json")
    
    return output

if __name__ == '__main__':
    run_e1_verification(verification_mode='sample')
