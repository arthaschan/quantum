#!/usr/bin/env python3
"""
B3 + B4 实验: 影子范数数值验证 + 自适应原根选择
===================================================
"""

import numpy as np
import json, os, time, sys

sys_path = os.path.join(os.path.dirname(__file__), '..', 'validate')
sys.path.insert(0, sys_path)
from core import build_mubs_wf, random_pure_state, primes_upto

def b3_shadow_norm_verification():
    """B3: 数值验证影子范数理论预测"""
    print("="*70)
    print("B3: SHADOW NORM NUMERICAL VERIFICATION")
    print("="*70)
    
    test_dims = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    n_states = 20
    r_per_dir = 2
    K_mom = None  # adaptive
    
    results = []
    
    for d in test_dims:
        mubs = build_mubs_wf(d)
        n_bases = len(mubs)
        n_total = n_bases * r_per_dir
        
        obs_list = []
        count = 0
        for a in range(d):
            if count >= 6: break
            for b in range(a+1, d):
                if count >= 6: break
                O = np.zeros((d,d),dtype=complex); O[a,b]=1.0; O[b,a]=1.0
                obs_list.append(O); count += 1
        
        ratios = []
        dir_var_ratios = []
        
        for si in range(n_states):
            psi = random_pure_state(d, seed=d*1000+si)
            rho = np.outer(psi, psi.conj())
            
            for oi, O in enumerate(obs_list[:3]):  # first 3 observables
                truth = np.real(np.trace(O @ rho))
                
                # Deterministic
                det_est = 0
                rng1 = np.random.RandomState(d*2000+si+oi)
                for a in range(n_bases):
                    basis = mubs[a]
                    probs = np.abs([np.conj(v)@rho@v for v in basis])
                    probs = np.clip(probs,0,None); probs /= probs.sum()
                    for _ in range(r_per_dir):
                        b = rng1.choice(d, p=probs)
                        det_est += (n_bases * np.abs(np.conj(basis[b])@O@basis[b]) - np.trace(O))
                det_est /= n_total
                
                # Random CS
                rand_ests = []
                rng2 = np.random.RandomState(d*3000+si+oi)
                for _ in range(10):  # 10 independent random runs
                    est = 0
                    for _ in range(n_total):
                        a = rng2.randint(n_bases)
                        basis = mubs[a]
                        probs = np.abs([np.conj(v)@rho@v for v in basis])
                        probs = np.clip(probs,0,None); probs /= probs.sum()
                        b = rng2.choice(d, p=probs)
                        est += (n_bases * np.abs(np.conj(basis[b])@O@basis[b]) - np.trace(O))
                    rand_ests.append(est / n_total)
                
                # 方向随机化方差
                dir_var = np.var([np.mean(rand_ests)] + [np.var(rand_ests)])
                shadow_norm_sq = d * np.real(np.trace(O@O))
                
                ratio = abs(det_est-truth) / max(np.mean([abs(e-truth) for e in rand_ests]), 1e-12)
                ratios.append(ratio)
                if shadow_norm_sq > 0:
                    dir_var_ratios.append(dir_var / shadow_norm_sq)
        
        print(f"  d={d:>3d}: n_states={n_states}, "
              f"median_ratio={np.median(ratios):.4f}, "
              f"dir_var_pct={np.median(dir_var_ratios)*100:.1f}%", flush=True)
        results.append({'d':d, 'median_ratio':float(np.median(ratios)),
                       'dir_var_pct':float(np.median(dir_var_ratios)*100),
                       'predicted_ratio': 1.0/(1.0+np.median(dir_var_ratios))})
    
    # 预测 vs 实测
    print(f"\n  Prediction vs. Measurement:")
    print(f"  {'d':>4s}  {'measured':>8s}  {'predicted':>10s}")
    for r in results:
        print(f"  {r['d']:>4d}  {r['median_ratio']:>8.4f}  {r['predicted_ratio']:>10.4f}")
    
    return results


def b4_adaptive_primitive_root():
    """B4: 自适应原根选择 — 对失败素数，搜索最优原根"""
    print("\n" + "="*70)
    print("B4: ADAPTIVE PRIMITIVE ROOT SELECTION")
    print("="*70)
    
    # 从原数据中找 ratio>1 的失败素数
    with open(os.path.join(os.path.dirname(__file__), '..', 'data', 'primes_1000_results.json')) as f:
        e1 = json.load(f)
    
    failed_primes = []
    for k, v in e1.items():
        d = v['d']
        best = min(r['ratio_1loc'] for r in v['runs'])
        if best > 1.0 and d >= 10:
            failed_primes.append((d, best))
    
    failed_primes.sort()
    key_failed = [d for d, r in failed_primes[:6]]
    print(f"  Testing {len(key_failed)} failed primes: {key_failed}")
    
    def primitive_roots(p):
        """Find all primitive roots modulo p"""
        phi = p - 1
        factors = []
        n = phi; d = 2
        while d*d <= n:
            if n%d == 0:
                factors.append(d)
                while n%d == 0: n//=d
            d += 1
        if n > 1: factors.append(n)
        
        roots = []
        for g in range(2, p):
            if all(pow(g, phi//q, p) != 1 for q in factors):
                roots.append(g)
        return roots
    
    results = []
    n_states = 5
    r_per_dir = 2
    
    for d in key_failed:
        roots = primitive_roots(d)
        np.random.shuffle(roots)  # random order
        max_roots = min(10, len(roots))
        roots = roots[:max_roots]
        
        mubs = build_mubs_wf(d)
        
        obs_list = []
        count = 0
        for a in range(d):
            if count >= 6: break
            for b in range(a+1, d):
                if count >= 6: break
                O = np.zeros((d,d),dtype=complex); O[a,b]=1.0; O[b,a]=1.0
                obs_list.append(O); count += 1
        
        root_best = {}
        for g in roots:
            # Rebuild MUBs with this primitive root
            omega = np.exp(2j*np.pi/d)
            mubs_g = [list(np.eye(d, dtype=complex))]
            for a in range(1, d+1):
                basis = []
                for b in range(d):
                    v = np.array([omega**(a*j*(j-1)//2 + b*j) for j in range(d)], dtype=complex)
                    basis.append(v / np.sqrt(d))
                mubs_g.append(basis)
            
            ratios_g = []
            for si in range(n_states):
                psi = random_pure_state(d, seed=g*1000+si)
                rho = np.outer(psi, psi.conj())
                
                mae_d = 0; mae_r = 0
                for O in obs_list[:3]:
                    truth = np.real(np.trace(O@rho))
                    
                    # Simple estimate
                    det_e = []; rand_e = []
                    for _ in range(5):
                        for a in range(d+1):
                            basis = mubs_g[a]
                            probs = np.abs([np.conj(v)@rho@v for v in basis])
                            probs = np.clip(probs,0,None); probs/=probs.sum()
                            b = np.random.choice(d, p=probs)
                            det_e.append((d+1)*np.abs(np.conj(basis[b])@O@basis[b])-np.trace(O))
                        
                        for _ in range((d+1)*r_per_dir):
                            a = np.random.randint(d+1)
                            basis = mubs_g[a]
                            probs = np.abs([np.conj(v)@rho@v for v in basis])
                            probs = np.clip(probs,0,None); probs/=probs.sum()
                            b = np.random.choice(d, p=probs)
                            rand_e.append((d+1)*np.abs(np.conj(basis[b])@O@basis[b])-np.trace(O))
                    
                    mae_d += abs(np.mean(det_e)-truth)
                    mae_r += abs(np.mean(rand_e)-truth)
                
                if mae_r > 0:
                    ratios_g.append(mae_d/mae_r)
            
            root_best[g] = np.median(ratios_g) if ratios_g else 1.0
        
        best_g = min(root_best, key=root_best.get)
        print(f"  d={d:>3d}: default_ratio={e1[f'd={d}']['runs'][0]['ratio_1loc']:.4f}, "
              f"best_root={best_g}, best_ratio={root_best[best_g]:.4f} "
              f"({'★ FLIP!' if root_best[best_g]<1.0 else 'still >1'})", flush=True)
        results.append({'d':d, 'best_root':best_g, 'best_ratio':float(root_best[best_g]),
                       'roots_tested': len(roots), 'all_ratios': {str(k):float(v) for k,v in root_best.items()}})
    
    flipped = sum(1 for r in results if r['best_ratio']<1.0)
    print(f"\n  Out of {len(results)} failed primes, {flipped} can be salvaged by root optimization")
    return results


if __name__ == '__main__':
    import sys
    b3 = b3_shadow_norm_verification()
    b4 = b4_adaptive_primitive_root()
    
    out = {'b3_shadow_norm': b3, 'b4_adaptive_root': b4}
    out_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'prx_experiments.json')
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n→ {out_path}")
