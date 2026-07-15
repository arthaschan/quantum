#!/usr/bin/env python3
"""
E3 多态验证: 对关键大素数进行多随机态 ratio 分析
======================================================
修正原始实验的单态采样缺陷，提供统计可靠的 ratio 估计。
"""

import numpy as np, json, time, os, sys

def run_multistate(d, n_states=10, n_per=4, n_obs=8, seed_base=42):
    """
    对单个维度 d 运行 n_states 个随机态的独立实验。
    返回 ratio 列表。
    """
    omega = np.exp(2j * np.pi / d)
    n_mubs = d + 1
    n_total = n_mubs * n_per
    sqrt_d = np.sqrt(d)
    
    ratios = []
    
    for si in range(n_states):
        # 生成随机态
        rng = np.random.RandomState(d * 100000 + seed_base * 100 + si)
        psi = rng.randn(d) + 1j * rng.randn(d)
        psi /= np.linalg.norm(psi)
        
        # 稀疏 1-local 可观测量
        obs_indices = []
        step = max(1, d // n_obs)
        for i in range(n_obs):
            a = i * step
            b = (i * step + step // 2) % d
            if a == b: b = (b + 1) % d
            obs_indices.append((a, b))
        
        truth = np.array([2 * np.real(np.conj(psi[a]) * psi[b]) for a, b in obs_indices])
        
        # FFT 计算所有 MUB 概率
        all_probs = [np.abs(psi)**2]
        for m in range(1, n_mubs):
            phase = np.array([omega**(-m * j * (j-1) // 2) for j in range(d)], dtype=complex)
            v = psi * phase
            f = np.fft.fft(v)
            probs = np.abs(f)**2 / d
            probs = np.clip(probs, 0, None); probs /= probs.sum()
            all_probs.append(probs)
        
        # 确定性测量
        det_m_k = []
        rng_d = np.random.RandomState(d * 200000 + seed_base * 100 + si)
        for m in range(n_mubs):
            cum = np.cumsum(all_probs[m])
            for _ in range(n_per):
                k = int(np.searchsorted(cum, rng_d.rand()))
                det_m_k.append((m, k))
        
        # 随机测量
        rand_m_k = []
        rng_r = np.random.RandomState(d * 300000 + seed_base * 100 + si)
        for _ in range(n_total):
            m = rng_r.randint(n_mubs)
            cum = np.cumsum(all_probs[m])
            k = int(np.searchsorted(cum, rng_r.rand()))
            rand_m_k.append((m, k))
        
        # 估计
        est_det = np.zeros(n_obs)
        for m_s, k_s in det_m_k:
            for oi, (a_val, b_val) in enumerate(obs_indices):
                if m_s == 0:
                    ca = 1.0 if a_val == k_s else 0.0
                    cb = 1.0 if b_val == k_s else 0.0
                else:
                    phase_a = m_s * a_val * (a_val - 1) // 2 + k_s * a_val
                    phase_b = m_s * b_val * (b_val - 1) // 2 + k_s * b_val
                    ca = omega**phase_a / sqrt_d
                    cb = omega**phase_b / sqrt_d
                est_det[oi] += 2 * np.real(np.conj(ca) * cb)
        est_det *= (d + 1) / n_total
        
        est_rand = np.zeros(n_obs)
        for m_s, k_s in rand_m_k:
            for oi, (a_val, b_val) in enumerate(obs_indices):
                if m_s == 0:
                    ca = 1.0 if a_val == k_s else 0.0
                    cb = 1.0 if b_val == k_s else 0.0
                else:
                    phase_a = m_s * a_val * (a_val - 1) // 2 + k_s * a_val
                    phase_b = m_s * b_val * (b_val - 1) // 2 + k_s * b_val
                    ca = omega**phase_a / sqrt_d
                    cb = omega**phase_b / sqrt_d
                est_rand[oi] += 2 * np.real(np.conj(ca) * cb)
        est_rand *= (d + 1) / n_total
        
        mae_d = np.mean(np.abs(est_det - truth))
        mae_r = np.mean(np.abs(est_rand - truth))
        r = mae_d / mae_r if mae_r > 0 else 1.0
        ratios.append(r)
    
    return ratios


if __name__ == '__main__':
    # 关键大素数
    key_primes = [2039, 2851, 4289, 1759, 3169, 6737, 7681]
    n_states = 10
    
    print("=" * 70)
    print(f"E3 MULTI-STATE VERIFICATION: {n_states} states per prime")
    print("=" * 70)
    
    all_results = {}
    
    for d in key_primes:
        print(f"\n  d={d} (p-1={d-1})...", flush=True)
        t0 = time.time()
        ratios = run_multistate(d, n_states=n_states, n_per=4, n_obs=8)
        dt = time.time() - t0
        
        median_r = np.median(ratios)
        mean_r = np.mean(ratios)
        min_r = np.min(ratios)
        max_r = np.max(ratios)
        dprt_win_rate = sum(1 for r in ratios if r < 1.0) / len(ratios) * 100
        
        print(f"    ratios: {[f'{r:.3f}' for r in ratios]}")
        print(f"    median={median_r:.4f}  mean={mean_r:.4f}  min={min_r:.4f}  max={max_r:.4f}")
        print(f"    DPRT wins: {dprt_win_rate:.0f}%  [{dt:.0f}s]", flush=True)
        
        all_results[d] = {
            'ratios': [float(r) for r in ratios],
            'median': float(median_r),
            'mean': float(mean_r),
            'min': float(min_r),
            'max': float(max_r),
            'dprt_win_rate_pct': float(dprt_win_rate),
            'time_s': float(dt)
        }
    
    # 汇总
    print(f"\n{'='*70}")
    print(f"MULTI-STATE SUMMARY ({n_states} states each)")
    print(f"{'='*70}")
    print(f"  {'d':>6s}  {'median':>8s}  {'mean':>8s}  {'min':>8s}  {'max':>8s}  {'DPRT%':>6s}")
    print(f"  {'-'*6}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*6}")
    
    for d in key_primes:
        r = all_results[d]
        print(f"  {d:>6d}  {r['median']:>8.4f}  {r['mean']:>8.4f}  "
              f"{r['min']:>8.4f}  {r['max']:>8.4f}  {r['dprt_win_rate_pct']:>5.0f}%")
    
    # 与原始单态数据的对比
    orig_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'primes_10000_sampled.json')
    if os.path.exists(orig_path):
        with open(orig_path) as f:
            orig_data = json.load(f)
        print(f"\n  Comparison with original single-state data:")
        print(f"  {'d':>6s}  {'orig':>8s}  {'our_median':>10s}  {'our_mean':>10s}  {'our_min':>10s}")
        print(f"  {'-'*6}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*10}")
        for r in orig_data['results']:
            if r['d'] in key_primes:
                d = r['d']
                print(f"  {d:>6d}  {r['ratio']:>8.4f}  {all_results[d]['median']:>10.4f}  "
                      f"{all_results[d]['mean']:>10.4f}  {all_results[d]['min']:>10.4f}")
    
    out_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'primes_multistate_results.json')
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\n→ {out_path}")
