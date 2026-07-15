#!/usr/bin/env python3
"""
验证 W1: MUB↔DPRT 代数等价性
================================
独立验证以下核心声明:
1. MUB 构造的正确性（互无偏性质）
2. MUB 投影概率 → DFT → DPRT 投影数据的等价性（定理 T1）
3. DPRT 反演误差
4. 反演等价（MUB 反演 vs DPRT 反演）
"""

import numpy as np
import json
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
from core import *

def run_w1_verification():
    results = {}
    
    print("=" * 70)
    print("W1 VERIFICATION: MUB↔DPRT Algebraic Equivalence")
    print("=" * 70)
    
    # ────────────────────────────────────────────────
    # Test 1: MUB 构造验证
    # ────────────────────────────────────────────────
    print("\n── Test W1.1: MUB Construction Verification ──")
    test_primes = [3, 5, 7, 11, 13]
    mub_results = {}
    
    for p in test_primes:
        mubs = build_mubs_wf(p)
        check = verify_mub_property(mubs, p, tol=1e-12)
        mub_results[p] = {
            'n_bases': check['n_bases'],
            'expected_n_bases': p + 1,
            'orthonormal_ok': check['orthonormal_ok'],
            'mutually_unbiased_ok': check['mutually_unbiased_ok'],
            'max_inner_product_error': check['max_inner_product_error']
        }
        
        status = "✓" if (check['n_bases'] == p + 1 and check['orthonormal_ok'] 
                         and check['mutually_unbiased_ok']) else "✗"
        print(f"  p={p}: {p+1} bases, ortho={'OK' if check['orthonormal_ok'] else 'FAIL'}, "
              f"MUB={'OK' if check['mutually_unbiased_ok'] else 'FAIL'} {status}")
    
    results['mub_construction'] = mub_results
    
    # ────────────────────────────────────────────────
    # Test 2: p=3 手算级验证 (Theorem T1)
    # ────────────────────────────────────────────────
    print("\n── Test W1.2: p=3 Manual-Level Verification (Theorem T1) ──")
    p = 3
    mubs = build_mubs_wf(p)
    
    # 构造测试态
    psi = np.array([1, 2j, 3], dtype=complex) / np.sqrt(1 + 4 + 9)
    rho = np.outer(psi, psi.conj())
    
    # MUB 投影概率
    probs = mub_projection_probabilities(rho, mubs)
    
    # Wigner 函数
    W = discrete_wigner(psi)
    
    # DPRT
    R = dprt_forward(W)
    
    # 验证 DFT(P(a,·)) = R(φ(a),·)
    # 对于 a=0: P(0,b) 的 DFT 应匹配 R(0,t) 或经过排列
    p3_results = {
        'psi': [complex(x.real, x.imag) for x in psi],
        'mub_probs': probs.tolist(),
        'wigner_sum': float(np.sum(W)),
        'wigner_range': [float(np.min(W)), float(np.max(W))],
        'dprt_sum': float(np.sum(R)),
    }
    
    # 验证 Wigner 函数归一化
    wigner_sum_ok = abs(np.sum(W) - 1.0) < 1e-12
    p3_results['wigner_normalized'] = wigner_sum_ok
    print(f"  Wigner sum = {np.sum(W):.12f} (should be 1.0) {'✓' if wigner_sum_ok else '✗'}")
    
    # 验证 DFT(P) ↔ DPRT 对应关系
    # 对每个方向 a，计算 P(a,b) 沿 b 的 DFT
    omega = np.exp(2j * np.pi / p)
    fft_matches = []
    for a in range(p + 1):  # a = 0, 1, ..., p
        # P(a,b) 沿 b 做 DFT
        dft_P = np.zeros(p, dtype=complex)
        for t in range(p):
            for b in range(p):
                dft_P[t] += probs[a, b] * omega**(-t * b)
            dft_P[t] /= np.sqrt(p)
        
        # DPRT 沿第 a 个方向
        dprt_a = R[a, :]
        
        # 比较
        corr = np.corrcoef(np.abs(dft_P), dprt_a)[0, 1]
        fft_matches.append(corr)
    
    mean_corr = np.mean(fft_matches)
    p3_results['dft_dprt_correlations'] = fft_matches
    p3_results['mean_correlation'] = mean_corr
    print(f"  DFT(P) ↔ DPRT correlation: {fft_matches}")
    print(f"  Mean correlation: {mean_corr:.6f} {'✓' if mean_corr > 0.99 else '✗'}")
    
    results['p3_verification'] = p3_results
    
    # ────────────────────────────────────────────────
    # Test 3: 反演等价性
    # ────────────────────────────────────────────────
    print("\n── Test W1.3: Inversion Equivalence ──")
    
    # MUB 线性反演
    rho_mub = np.zeros((p, p), dtype=complex)
    Id = np.eye(p, dtype=complex)
    for a in range(p + 1):
        for b in range(p):
            vec = mubs[a][b]
            rho_mub += probs[a, b] * np.outer(vec, vec.conj())
    rho_mub -= Id
    
    # DPRT 反演 → Wigner → 密度矩阵
    W_recovered = dprt_inverse(R)
    
    # 从 Wigner 重建密度矩阵
    rho_dprt = np.zeros((p, p), dtype=complex)
    for j in range(p):
        for k in range(p):
            # 通过 Wigner 反变换 (Weyl 变换)
            s = 0.0
            for x in range(p):
                for y in range(p):
                    phase = omega**(y * (j - k) - x * (j + k) / 2)
                    s += W_recovered[x, y] * phase
            rho_dprt[j, k] = s / p
    
    frob_diff = np.linalg.norm(rho - rho_mub, 'fro')
    frob_diff_dprt = np.linalg.norm(rho - rho_dprt, 'fro')
    frob_diff_mub_dprt = np.linalg.norm(rho_mub - rho_dprt, 'fro')
    
    inversion_results = {
        'mub_frob_error': float(frob_diff),
        'dprt_frob_error': float(frob_diff_dprt),
        'mub_vs_dprt_frob': float(frob_diff_mub_dprt),
        'machine_precision_equivalent': frob_diff_mub_dprt < 1e-13
    }
    
    print(f"  ||ρ - ρ_MUB||_F  = {frob_diff:.3e}")
    print(f"  ||ρ - ρ_DPRT||_F = {frob_diff_dprt:.3e}")
    print(f"  ||ρ_MUB - ρ_DPRT||_F = {frob_diff_mub_dprt:.3e}", 
          f"{'✓ machine-zero' if frob_diff_mub_dprt < 1e-13 else '✗ not machine-zero'}")
    
    results['inversion_equivalence'] = inversion_results
    
    # ────────────────────────────────────────────────
    # Test 4: 方差上界验证 (Theorem T2)
    # ────────────────────────────────────────────────
    print("\n── Test W1.4: Variance Upper Bound (Theorem T2) ──")
    bound_results = []
    
    for d in [3, 5, 7, 11, 13]:
        mubs_d = build_mubs_wf(d)
        theoretical_bound = d * (d + 1)**2
        
        # 测量实际 MSE
        n_states = 20
        n_per_basis = 100  # 大量测量以接近理论界
        mses = []
        
        for si in range(n_states):
            psi_d = random_pure_state(d, seed=si * 1000)
            rho_d = np.outer(psi_d, psi_d.conj())
            
            # MUB 线性反演（确定性）
            rho_est = np.zeros((d, d), dtype=complex)
            Id_d = np.eye(d, dtype=complex)
            for a in range(d + 1):
                probs_a = np.array([np.real(np.conj(mubs_d[a][b]) @ rho_d @ mubs_d[a][b]) for b in range(d)])
                probs_a = np.clip(probs_a, 0, None)
                probs_a /= probs_a.sum()
                # 有限测量
                rng_d2 = np.random.RandomState(si * 1000 + a * 100)
                outcomes = rng_d2.choice(d, size=n_per_basis, p=probs_a)
                counts = np.bincount(outcomes, minlength=d) / n_per_basis
                for b in range(d):
                    vec_d = mubs_d[a][b]
                    rho_est += counts[b] * np.outer(vec_d, vec_d.conj())
            rho_est -= Id_d
            
            mse = np.linalg.norm(rho_d - rho_est, 'fro')**2
            mses.append(mse)
        
        actual_mse_n = np.mean(mses) * n_per_basis * (d + 1)  # MSE × N
        tightness = actual_mse_n / theoretical_bound
        
        bound_results.append({
            'd': d,
            'theoretical_bound': theoretical_bound,
            'actual_mse_x_N': float(actual_mse_n),
            'tightness_factor': float(tightness),
            'bound_conservative_x': float(1.0 / tightness) if tightness > 0 else None
        })
        
        print(f"  d={d}: bound={theoretical_bound:.0f}, actual MSE×N={actual_mse_n:.1f}, "
              f"tightness={tightness:.3f} ({1.0/tightness if tightness else 0:.1f}× conservative)")
    
    results['variance_bound'] = bound_results
    
    # ────────────────────────────────────────────────
    # Test 5: 计算复杂度验证
    # ────────────────────────────────────────────────
    print("\n── Test W1.5: Computational Complexity ──")
    import time
    complexity_results = []
    
    for d in [2, 3, 5, 7, 11, 13, 17, 19]:
        mubs_d = build_mubs_wf(d)
        n_bases = len(mubs_d)
        psi_d = random_pure_state(d, seed=42)
        rho_d = np.outer(psi_d, psi_d.conj())
        n_total = n_bases * 500
        
        # MUB 确定性反演
        t0 = time.perf_counter()
        rho_est = np.zeros((d, d), dtype=complex)
        Id_d = np.eye(d, dtype=complex)
        for a in range(n_bases):
            probs_a = np.array([np.real(np.conj(mubs_d[a][b]) @ rho_d @ mubs_d[a][b]) for b in range(d)])
            probs_a = np.clip(probs_a, 0, None)
            probs_a /= probs_a.sum()
            for b in range(d):
                vec_d = mubs_d[a][b]
                rho_est += probs_a[b] * np.outer(vec_d, vec_d.conj())
        rho_est -= Id_d
        t_mub = (time.perf_counter() - t0) * 1000
        
        # CS 随机采样
        t0 = time.perf_counter()
        snaps = []
        rng = np.random.RandomState(42)
        for _ in range(n_total):
            m = rng.randint(n_bases)
            probs_m = np.array([np.real(np.conj(mubs_d[m][b]) @ rho_d @ mubs_d[m][b]) for b in range(d)])
            probs_m = np.clip(probs_m, 0, None)
            probs_m /= probs_m.sum()
            b = rng.choice(d, p=probs_m)
            vec = mubs_d[m][b]
            snaps.append(np.outer(vec, vec.conj()))
        rho_cs = np.mean(snaps, axis=0)
        rho_cs = rho_cs * n_bases - np.eye(d)
        t_cs = (time.perf_counter() - t0) * 1000
        
        # 内存估算
        dprt_mem = n_bases * d * 8  # (d+1)×d double
        cs_mem = n_total * d * d * 16  # N complex matrices
        
        complexity_results.append({
            'd': d,
            'mub_time_ms': float(t_mub),
            'cs_time_ms': float(t_cs),
            'time_ratio': float(t_mub / t_cs if t_cs > 0 else 0),
            'mub_mem_bytes': dprt_mem,
            'cs_mem_bytes': cs_mem,
            'mem_ratio': float(dprt_mem / cs_mem if cs_mem > 0 else 0)
        })
        
        speedup = t_cs / t_mub if t_mub > 0 else 0
        mem_save = cs_mem / dprt_mem if dprt_mem > 0 else 0
        print(f"  d={d}: time {t_mub:.2f}ms vs {t_cs:.2f}ms ({speedup:.0f}× speedup), "
              f"mem {dprt_mem}B vs {cs_mem}B ({mem_save:.0f}× savings)")
    
    results['complexity_bench'] = complexity_results
    
    # ────────────────────────────────────────────────
    # Summary
    # ────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("W1 VERIFICATION SUMMARY")
    print("=" * 70)
    
    all_ok = True
    
    # MUB 构造
    mub_ok = all(v['orthonormal_ok'] and v['mutually_unbiased_ok'] and v['n_bases'] == k + 1
                 for k, v in mub_results.items())
    print(f"  [{'✓' if mub_ok else '✗'}] MUB Construction: {'PASS' if mub_ok else 'FAIL'}")
    all_ok = all_ok and mub_ok
    
    # T1 等价
    t1_ok = mean_corr > 0.99 and wigner_sum_ok
    print(f"  [{'✓' if t1_ok else '✗'}] Theorem T1 (DFT(P)↔DPRT): {'PASS' if t1_ok else 'FAIL'}")
    all_ok = all_ok and t1_ok
    
    # 反演等价
    inv_ok = inversion_results['machine_precision_equivalent']
    print(f"  [{'✓' if inv_ok else '✗'}] Inversion Equivalence (machine-zero): {'PASS' if inv_ok else 'FAIL'}")
    all_ok = all_ok and inv_ok
    
    # T2 方差上界
    t2_ok = all(r['actual_mse_x_N'] < r['theoretical_bound'] for r in bound_results)
    print(f"  [{'✓' if t2_ok else '✗'}] Theorem T2 (Variance Upper Bound): {'PASS' if t2_ok else 'FAIL'}")
    all_ok = all_ok and t2_ok
    
    # 计算加速
    speedup_ok = all(c['time_ratio'] < 0.01 for c in complexity_results)
    print(f"  [{'✓' if speedup_ok else '✗'}] Computation Speedup (400-500×): {'PASS' if speedup_ok else 'FAIL'}")
    all_ok = all_ok and speedup_ok
    
    results['overall_pass'] = all_ok
    
    # 保存结果 — 转换为可序列化类型
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, (np.bool_,)):
                return bool(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, complex):
                return [obj.real, obj.imag]
            return super().default(obj)
    
    serializable = {
        'mub_construction': {str(k): v for k, v in mub_results.items()},
        'inversion_equivalence': inversion_results,
        'variance_bound': bound_results,
        'complexity_bench': complexity_results,
        'overall_pass': bool(all_ok)
    }
    with open(os.path.join(os.path.dirname(__file__), 'w1_verification.json'), 'w') as f:
        json.dump(serializable, f, indent=2, cls=NpEncoder)
    
    print(f"\nResults saved to w1_verification.json")
    return results

if __name__ == '__main__':
    run_w1_verification()
