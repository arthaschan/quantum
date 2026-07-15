#!/usr/bin/env python3
"""
验证 E7: 多通道噪声鲁棒性分析
===============================
独立验证 4 种噪声通道下 DPRT 优势的存续情况。
"""

import numpy as np
import json
import sys
import os
import time

sys.path.insert(0, os.path.dirname(__file__))
from core import *

def run_e7_verification():
    results = {}
    
    print("=" * 70)
    print("E7 VERIFICATION: Noise Robustness Analysis")
    print("=" * 70)
    
    # 实验参数
    test_dims = [3, 5, 7, 11, 13, 17, 19, 23]
    noise_levels = [0.0, 0.01, 0.05, 0.10, 0.20]
    channels = ['depol', 'phase', 'bitflip', 'amp']
    
    n_states_per = 5  # 每个条件测试的随机态数
    delta = 10  # 1-local 可观测量数
    
    all_data = []  # 所有数据点
    
    print(f"\n  Config: {len(test_dims)} dims × {len(channels)} channels "
          f"× {len(noise_levels)} noise levels × {n_states_per} states")
    print(f"  Total data points: {len(test_dims) * len(channels) * len(noise_levels) * n_states_per}")
    print()
    
    T0 = time.time()
    
    for d in test_dims:
        mubs = build_mubs_wf(d)
        n_mubs = len(mubs)
        n_total_det = n_mubs * 2  # r=2 per direction
        
        print(f"  d={d}: ", end="", flush=True)
        d_start = time.time()
        
        for ch in channels:
            best_advantages = []
            
            for lam in noise_levels:
                ratio_vals = []
                
                for si in range(n_states_per):
                    seed = d * 10000 + hash(ch) % 10000 + int(lam * 1000) + si * 100
                    psi = random_pure_state(d, seed=seed)
                    rho_pure = np.outer(psi, psi.conj())
                    
                    # 施加噪声
                    if lam > 0:
                        rho = apply_noise(rho_pure, ch, lam)
                    else:
                        rho = rho_pure
                    
                    # 确定性估计
                    det = estimate_1local_expectations(rho, mubs, shots_per_dir=2, delta=delta)
                    
                    # 随机估计
                    rand = estimate_1local_random(rho, mubs, n_total=n_total_det, 
                                                   delta=delta, seed=seed + 1)
                    
                    ratio = det['mae'] / rand['mae'] if rand['mae'] > 0 else 1.0
                    ratio_vals.append(ratio)
                
                mean_ratio = np.mean(ratio_vals)
                advantage = (1 - mean_ratio) * 100
                
                all_data.append({
                    'd': d,
                    'channel': ch,
                    'lambda': lam,
                    'mean_ratio': float(mean_ratio),
                    'advantage_pct': float(advantage),
                    'ratio_std': float(np.std(ratio_vals)),
                    'winner': 'DPRT' if advantage > 0 else 'MUB/CS'
                })
                
                best_advantages.append(advantage)
            
            # 该通道的最佳表现
            best_lam_idx = np.argmax([a['advantage_pct'] for a in all_data[-5:]])
        
        dt = time.time() - d_start
        print(f"{dt:.0f}s", flush=True)
    
    total_time = time.time() - T0
    print(f"\n  Total: {total_time:.1f}s")
    
    # ────────────────────────────────────────────────
    # 分析
    # ────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("NOISE ANALYSIS SUMMARY")
    print("=" * 70)
    
    # 每个通道的平均表现
    print("\n── Per-Channel Summary (excluding λ=0) ──")
    channel_summary = {}
    for ch in channels:
        ch_data = [d for d in all_data if d['channel'] == ch and d['lambda'] > 0]
        mean_adv = np.mean([d['advantage_pct'] for d in ch_data])
        pct_positive = sum(1 for d in ch_data if d['advantage_pct'] > 0) / len(ch_data) * 100
        
        channel_summary[ch] = {
            'mean_advantage_pct': float(mean_adv),
            'pct_positive': float(pct_positive),
            'n_points': len(ch_data)
        }
        
        print(f"  {ch:12s}: mean_adv={mean_adv:+.1f}%, positive={pct_positive:.0f}%, "
              f"n={len(ch_data)}")
    
    results['channel_summary'] = channel_summary
    
    # Top-10 最强数据点
    print("\n── Top-10 Strongest DPRT Advantages ──")
    sorted_data = sorted(all_data, key=lambda x: x['advantage_pct'], reverse=True)
    top10 = sorted_data[:10]
    
    for i, d in enumerate(top10):
        print(f"  {i+1:2d}. ch={d['channel']:8s}  λ={d['lambda']:.2f}  d={d['d']:3d}  "
              f"ratio={d['mean_ratio']:.3f}  adv={d['advantage_pct']:+.1f}%")
    
    results['top10_advantages'] = top10
    
    # 噪声作为正则化器：检查是否 λ>0 的 ratio < λ=0 的 ratio
    print("\n── Noise as Regularizer: λ>0 vs λ=0 ──")
    regularizer_cases = []
    
    for ch in channels:
        for d in test_dims:
            lam0 = [x for x in all_data if x['channel'] == ch and x['d'] == d and x['lambda'] == 0.0]
            if not lam0: continue
            base_ratio = lam0[0]['mean_ratio']
            
            for lam in noise_levels[1:]:
                lam_data = [x for x in all_data if x['channel'] == ch and x['d'] == d 
                           and abs(x['lambda'] - lam) < 0.001]
                if lam_data and lam_data[0]['mean_ratio'] < base_ratio:
                    regularizer_cases.append({
                        'd': d,
                        'channel': ch,
                        'lambda': lam,
                        'base_ratio': float(base_ratio),
                        'noisy_ratio': float(lam_data[0]['mean_ratio']),
                        'improvement_pct': float((base_ratio - lam_data[0]['mean_ratio']) * 100)
                    })
    
    print(f"  Found {len(regularizer_cases)} cases where noise improved DPRT advantage")
    for case in regularizer_cases[:10]:
        print(f"    d={case['d']}, {case['channel']}, λ={case['lambda']:.2f}: "
              f"ratio {case['base_ratio']:.3f} → {case['noisy_ratio']:.3f} "
              f"({case['improvement_pct']:.1f}% improvement)")
    
    results['regularizer_cases'] = regularizer_cases
    
    # ────────────────────────────────────────────────
    # 保存
    # ────────────────────────────────────────────────
    out_path = os.path.join(os.path.dirname(__file__), 'e7_verification.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n  Results saved to e7_verification.json")
    
    return results

if __name__ == '__main__':
    run_e7_verification()
