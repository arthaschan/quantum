#!/usr/bin/env python3
"""
M6 + M8 实验数据生成
M6: 用真实实验数据生成噪声热力图
M8: 对 d≤7 加入 MLE 对比
"""
import numpy as np, json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'validate'))
from core import build_mubs_wf, random_pure_state

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

def m6_real_noise_heatmap():
    """用真实 MUB 构造计算噪声通道下的 DPRT ratio"""
    print("M6: REAL NOISE HEATMAP")
    
    test_dims = [3, 5, 7, 11, 13]
    channels = {
        'Depolarizing': lambda rho, lam: (1-lam)*rho + lam*np.eye(len(rho))/len(rho),
        'Phase Damping': lambda rho, lam: (1-lam)*rho + lam*np.diag(np.diag(rho)),
    }
    noise_levels = [0.0, 0.01, 0.05, 0.1, 0.2]
    n_states = 10
    
    data = np.zeros((len(channels), len(test_dims), len(noise_levels)))
    
    for ci, (cname, channel) in enumerate(channels.items()):
        for di, d in enumerate(test_dims):
            mubs = build_mubs_wf(d)
            n_bases = len(mubs)
            r_per_dir = 2
            n_total = n_bases * r_per_dir
            
            for ni, lam in enumerate(noise_levels):
                ratios = []
                for si in range(n_states):
                    psi = random_pure_state(d, seed=d*10000+si)
                    rho = np.outer(psi, psi.conj())
                    rho_n = channel(rho, lam)
                    
                    # 1-local observable (0,1)
                    truth = 2 * np.real(rho_n[0,1])
                    
                    # DPRT
                    rng = np.random.RandomState(d*2000+si)
                    est_d = 0
                    for m in range(n_bases):
                        basis = mubs[m]
                        probs = np.abs([np.conj(v)@rho_n@v for v in basis])
                        probs = np.clip(probs, 0, None); probs /= probs.sum()
                        for _ in range(r_per_dir):
                            b = rng.choice(d, p=probs)
                            val = 2*np.real(np.conj(basis[b][0])*basis[b][1])
                            est_d += (n_bases*val) / n_total
                    
                    # CS
                    rng2 = np.random.RandomState(d*3000+si)
                    est_c = 0
                    for _ in range(n_total):
                        m = rng2.randint(n_bases)
                        basis = mubs[m]
                        probs = np.abs([np.conj(v)@rho_n@v for v in basis])
                        probs = np.clip(probs, 0, None); probs /= probs.sum()
                        b = rng2.choice(d, p=probs)
                        val = 2*np.real(np.conj(basis[b][0])*basis[b][1])
                        est_c += (n_bases*val) / n_total
                    
                    mae_d = np.mean(np.abs(est_d - truth))
                    mae_c = np.mean(np.abs(est_c - truth))
                    if mae_c > 0:
                        ratios.append(mae_d / mae_c)
                
                data[ci, di, ni] = np.median(ratios) if ratios else 1.0
                print(f"  {cname:>15s} d={d:>3d} λ={lam:.2f}: ratio={data[ci,di,ni]:.4f}", flush=True)
    
    # Save
    result = {
        'channels': list(channels.keys()),
        'dims': test_dims,
        'noise_levels': noise_levels,
        'data': data.tolist(),
    }
    with open(os.path.join(DATA_DIR, 'real_noise_heatmap.json'), 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"  → {os.path.join(DATA_DIR, 'real_noise_heatmap.json')}")
    return result


def m8_mle_comparison():
    """MLE 对比 — d≤7，使用凸优化的最大似然估计"""
    print("\nM8: MLE COMPARISON")
    
    from scipy.optimize import minimize
    
    def mle_reconstruct(measurements, d):
        """简化的 MLE 重建：最小化负对数似然"""
        # measurements: list of (basis, outcome_count) tuples
        def neg_log_likelihood(params):
            # 参数化为 Cholesky 因子（保证半正定）
            L = np.zeros((d,d), dtype=complex)
            idx = 0
            for i in range(d):
                for j in range(i+1):
                    if i == j:
                        L[i,j] = params[idx]  # 对角元为实
                        idx += 1
                    else:
                        L[i,j] = params[idx] + 1j*params[idx+1]
                        idx += 2
            
            rho = L @ L.conj().T
            rho = rho / np.trace(rho)  # 归一化
            
            # Compute log-likelihood
            nll = 0
            for basis, counts in measurements:
                probs = np.abs([np.conj(v)@rho@v for v in basis])
                probs = np.clip(probs, 1e-12, None)
                probs /= probs.sum()
                nll -= np.sum(counts * np.log(probs))
            return nll
        
        n_params = d*d  # d² real parameters
        init = np.eye(d).flatten().real[:n_params] * 0.1 + 0.5/np.sqrt(d)
        
        result = minimize(neg_log_likelihood, init, method='L-BFGS-B', 
                         options={'maxiter': 1000})
        
        # Reconstruct rho from optimized params
        L = np.zeros((d,d), dtype=complex)
        idx = 0
        for i in range(d):
            for j in range(i+1):
                if i == j:
                    L[i,j] = result.x[idx]; idx += 1
                else:
                    L[i,j] = result.x[idx] + 1j*result.x[idx+1]; idx += 2
        rho_mle = L @ L.conj().T
        return rho_mle / np.trace(rho_mle)
    
    test_dims = [2, 3, 5, 7]
    results_mle = []
    
    for d in test_dims:
        mubs = build_mubs_wf(d)
        n_bases = len(mubs)
        r_per_dir = 20  # more shots for MLE
        
        for si in range(5):  # fewer states (MLE is slow)
            psi = random_pure_state(d, seed=d*1000+si)
            rho = np.outer(psi, psi.conj())
            
            # Generate measurements for MLE
            measurements = []
            for m in range(n_bases):
                basis = mubs[m]
                probs = np.abs([np.conj(v)@rho@v for v in basis])
                probs = np.clip(probs, 0, None); probs /= probs.sum()
                counts = np.random.multinomial(r_per_dir, probs)
                measurements.append((basis, counts))
            
            # MLE
            rho_mle = mle_reconstruct(measurements, d)
            
            # DPRT estimate (same data)
            rho_d = np.zeros((d,d), dtype=complex)
            for m, (basis, counts) in enumerate(measurements):
                for b, cnt in enumerate(counts):
                    vec = basis[b]
                    rho_d += cnt * np.outer(vec, vec.conj())
            rho_d = rho_d * (d+1) / (n_bases * r_per_dir) - np.eye(d)
            
            # CS estimate
            n_total = n_bases * r_per_dir
            rho_c = np.zeros((d,d), dtype=complex)
            rng = np.random.RandomState(d*3000+si)
            for _ in range(n_total):
                m = rng.randint(n_bases)
                basis = mubs[m]
                probs = np.abs([np.conj(v)@rho@v for v in basis])
                probs = np.clip(probs, 0, None); probs /= probs.sum()
                b = rng.choice(d, p=probs)
                vec = basis[b]
                rho_c += np.outer(vec, vec.conj())
            rho_c = rho_c * (d+1) / n_total - np.eye(d)
            
            # 1-local observable
            a, b = 0, d//2
            truth = 2 * np.real(rho[a,b])
            
            mae_d = abs(2*np.real(rho_d[a,b]) - truth)
            mae_c = abs(2*np.real(rho_c[a,b]) - truth)
            mae_mle = abs(2*np.real(rho_mle[a,b]) - truth)
            
            entry = {
                'd': d, 'state': si,
                'mae_dprt': float(mae_d),
                'mae_cs': float(mae_c),
                'mae_mle': float(mae_mle),
                'ratio_dprt_cs': float(mae_d/mae_c if mae_c>0 else 1),
                'ratio_dprt_mle': float(mae_d/mae_mle if mae_mle>0 else 1),
            }
            results_mle.append(entry)
            print(f"  d={d} s={si}: DPRT={mae_d:.4f} CS={mae_c:.4f} MLE={mae_mle:.4f} "
                  f"DPRT/CS={entry['ratio_dprt_cs']:.3f} DPRT/MLE={entry['ratio_dprt_mle']:.3f}", flush=True)
    
    # Summary
    for d in test_dims:
        entries = [e for e in results_mle if e['d']==d]
        avg_dc = np.mean([e['ratio_dprt_cs'] for e in entries])
        avg_dm = np.mean([e['ratio_dprt_mle'] for e in entries])
        print(f"\n  d={d} summary: DPRT/CS={avg_dc:.3f} DPRT/MLE={avg_dm:.3f}")
    
    with open(os.path.join(DATA_DIR, 'mle_comparison.json'), 'w') as f:
        json.dump(results_mle, f, indent=2)
    
    return results_mle


if __name__ == '__main__':
    m6_real_noise_heatmap()
    m8_mle_comparison()
