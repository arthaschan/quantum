#!/usr/bin/env python3
"""
RadonShadow 独立验证 — 主运行器
=================================
按顺序执行所有验证实验，收集结果，生成综合验证报告。
"""

import numpy as np
import json
import sys
import os
import time
from datetime import datetime

# 确保可以导入 core
sys.path.insert(0, os.path.dirname(__file__))

def main():
    vdir = os.path.dirname(__file__)
    
    print("=" * 70)
    print("RADONSHADOW INDEPENDENT VERIFICATION SUITE")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    
    total_start = time.time()
    all_pass = True
    verification_results = {}
    
    # ────────────────────────────────────────────────
    # V0: T1 手算验证 (p=3)
    # ────────────────────────────────────────────────
    print("\n\n" + "▐" + "█" * 68 + "▌")
    print("▐  PHASE 0: T1 Manual Verification (p=3)  ▐")
    print("▐" + "█" * 68 + "▌")
    
    from verify_t1_manual import run_t1_manual_verification
    t1_result = run_t1_manual_verification()
    verification_results['t1_manual'] = {
        'pass': t1_result['t1_equivalent'],
        'frob_mub_vs_dprt': t1_result['frob_mub_vs_dprt']
    }
    if not t1_result['t1_equivalent']:
        all_pass = False
    
    # ────────────────────────────────────────────────
    # V1: W1 形式化验证
    # ────────────────────────────────────────────────
    print("\n\n" + "▐" + "█" * 68 + "▌")
    print("▐  PHASE 1: W1 Formal Verification  ▐")
    print("▐" + "█" * 68 + "▌")
    
    from verify_w1 import run_w1_verification
    w1_result = run_w1_verification()
    verification_results['w1'] = {
        'overall_pass': w1_result['overall_pass'],
        'mub_construction': {str(k): v for k, v in w1_result['mub_construction'].items()},
        'inversion_equivalence': w1_result['inversion_equivalence']['machine_precision_equivalent'],
        'variance_bound_verified': all(
            r['actual_mse_x_N'] < r['theoretical_bound'] 
            for r in w1_result['variance_bound']
        )
    }
    if not w1_result['overall_pass']:
        all_pass = False
    
    # ────────────────────────────────────────────────
    # V2: E1 素数扫描
    # ────────────────────────────────────────────────
    print("\n\n" + "▐" + "█" * 68 + "▌")
    print("▐  PHASE 2: E1 Prime Scan  ▐")
    print("▐" + "█" * 68 + "▌")
    
    from verify_e1 import run_e1_verification
    e1_result = run_e1_verification(verification_mode='sample')
    verification_results['e1'] = {
        'n_primes': e1_result['statistics']['n_primes'],
        'median_ratio': e1_result['statistics']['median_ratio'],
        'median_advantage_pct': e1_result['statistics']['median_advantage_pct'],
        'ratio_lt_090_pct': e1_result['statistics']['ratio_lt_090']['pct'],
        'dprt_wins_pct': e1_result['statistics']['dprt_wins'] / e1_result['statistics']['n_primes'] * 100,
        'bootstrap_median_95ci': e1_result['statistics']['bootstrap_median_95ci'],
        'comparison': e1_result.get('comparison_with_original')
    }
    
    # ────────────────────────────────────────────────
    # V3: E3 大素数 + 相关性
    # ────────────────────────────────────────────────
    print("\n\n" + "▐" + "█" * 68 + "▌")
    print("▐  PHASE 3: E3 Large Primes + Correlation  ▐")
    print("▐" + "█" * 68 + "▌")
    
    from verify_e3 import run_e3_verification
    e3_result = run_e3_verification()
    verification_results['e3'] = {
        'key_primes': e3_result['key_primes'],
        'spearman_rho': e3_result['correlation']['spearman_rho'],
        'pearson_r': e3_result['correlation']['pearson_r'],
        'n_primes_for_correlation': e3_result['correlation']['n_primes_analyzed']
    }
    
    # ────────────────────────────────────────────────
    # V4: E7 噪声鲁棒性
    # ────────────────────────────────────────────────
    print("\n\n" + "▐" + "█" * 68 + "▌")
    print("▐  PHASE 4: E7 Noise Robustness  ▐")
    print("▐" + "█" * 68 + "▌")
    
    from verify_e7 import run_e7_verification
    e7_result = run_e7_verification()
    verification_results['e7'] = {
        'channel_summary': e7_result['channel_summary'],
        'top10_advantages': e7_result['top10_advantages'],
        'n_regularizer_cases': len(e7_result['regularizer_cases'])
    }
    
    # ────────────────────────────────────────────────
    # 汇总
    # ────────────────────────────────────────────────
    total_time = time.time() - total_start
    
    # 读取原始实验数据进行对比
    orig_comparison = {}
    
    # 对比 E1 原始数据
    orig_e1_path = os.path.join(vdir, '..', 'data', 'primes_1000_results.json')
    if os.path.exists(orig_e1_path):
        with open(orig_e1_path) as f:
            orig_e1 = json.load(f)
        orig_ratios = {}
        for k, v in orig_e1.items():
            orig_ratios[v['d']] = min(r['ratio_1loc'] for r in v['runs'])
        
        orig_e1_median = np.median(list(orig_ratios.values()))
        orig_e1_lt090 = sum(1 for r in orig_ratios.values() if r < 0.90) / len(orig_ratios) * 100
        orig_comparison['e1'] = {
            'n_primes': len(orig_ratios),
            'median_ratio': float(orig_e1_median),
            'ratio_lt_090_pct': float(orig_e1_lt090),
        }
    
    # 对比 E3 原始数据
    orig_e3_path = os.path.join(vdir, '..', 'data', 'primes_10000_sampled.json')
    if os.path.exists(orig_e3_path):
        with open(orig_e3_path) as f:
            orig_e3 = json.load(f)
        orig_comparison['e3_d2039'] = float(orig_e3.get('d=2039', {}).get('ratio', 'N/A') if isinstance(orig_e3, dict) else 'N/A')
    
    summary = {
        'timestamp': datetime.now().isoformat(),
        'total_time_s': float(total_time),
        'verification_results': verification_results,
        'original_data_comparison': orig_comparison,
        'overall_pass': all_pass
    }
    
    # ────────────────────────────────────────────────
    # 最终报告
    # ────────────────────────────────────────────────
    print("\n\n" + "█" * 70)
    print("█" + " " * 68 + "█")
    print("█" + "  RADONSHADOW INDEPENDENT VERIFICATION — FINAL REPORT".center(68) + "█")
    print("█" + " " * 68 + "█")
    print("█" * 70)
    
    print(f"\n  Total computation time: {total_time:.0f}s ({total_time/60:.1f}min)")
    print(f"\n  ┌{'─'*60}┐")
    
    sections = [
        ("T1 (p=3 Manual)", 
         f"MUB↔DPRT machine-zero: {verification_results['t1_manual']['pass']}"),
        ("W1 (Formal)", 
         f"Overall: {verification_results['w1']['overall_pass']}"),
        ("E1 (Prime Scan)", 
         f"Median ratio: {verification_results['e1']['median_ratio']:.4f} "
         f"({verification_results['e1']['median_advantage_pct']:.1f}% advantage)"),
        ("E3 (Large Primes + ρ)", 
         f"Spearman ρ: {verification_results['e3']['spearman_rho']:.4f}"),
        ("E7 (Noise)", 
         f"Regularizer cases: {verification_results['e7']['n_regularizer_cases']}")
    ]
    
    for name, detail in sections:
        print(f"  │ {name:<15s} │ {detail:<42s} │")
    
    print(f"  └{'─'*60}┘")
    
    # 与原始数据对比
    if orig_comparison:
        print(f"\n  ── Comparison with Original Published Data ──")
        if 'e1' in orig_comparison:
            print(f"  E1 Original: {orig_comparison['e1']['n_primes']} primes, "
                  f"median={orig_comparison['e1']['median_ratio']:.4f}, "
                  f"<0.90: {orig_comparison['e1']['ratio_lt_090_pct']:.1f}%")
    
    # 保存
    out_path = os.path.join(vdir, 'verification_summary.json')
    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n  Full summary saved to: verification_summary.json")
    print(f"\n  Individual results: t1_manual_verification.json, w1_verification.json,")
    print(f"                       e1_verification.json, e3_verification.json, e7_verification.json")
    
    print(f"\n  OVERALL VERDICT: {'✓ ALL CLAIMS VERIFIED' if all_pass else '✗ SOME CLAIMS FAILED'}")
    print("=" * 70)
    
    return summary

if __name__ == '__main__':
    main()
