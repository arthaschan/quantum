#!/usr/bin/env python3
"""
RadonShadow: Final 5 experiments
================================
E1: Primitive root → ratio mapping (all 168 primes ≤1000)
E2: Multiplicative subgroup structure analysis
E3: Composite dimensions d=4,6,8,9 via K&S directions
E4: 2-local observables for d≤100 primes
E5: Depolarizing noise robustness

Output: all5_experiments_20260714.md + all5_results.json
"""
import numpy as np, math, time, json, sys
from collections import defaultdict

# ═══════════════════════════════════════
# E1: PRIMITIVE ROOT ANALYSIS
# ═══════════════════════════════════════
def primitive_root(p):
    """Find smallest primitive root modulo prime p."""
    if p == 2: return 1
    # Factor p-1
    phi = p - 1
    factors = []
    n = phi
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    
    for g in range(2, p):
        if all(pow(g, phi // q, p) != 1 for q in factors):
            return g
    return None

def e1_primitive_root_analysis(results_json_path):
    """Analyze ratio vs primitive root."""
    with open(results_json_path) as f:
        data = json.load(f)
    
    # Build mapping: d → best 1loc ratio
    d_to_ratio = {}
    for k, v in data.items():
        d = v['d']
        best = min(r['ratio_1loc'] for r in v['runs'])
        d_to_ratio[d] = best
    
    # Group by primitive root
    by_root = defaultdict(list)
    for d, ratio in sorted(d_to_ratio.items()):
        g = primitive_root(d)
        if g is not None:
            by_root[g].append((d, ratio))
    
    # Compute stats per root
    print("═"*60)
    print("E1: PRIMITIVE ROOT ANALYSIS")
    print("═"*60)
    
    root_stats = {}
    for g in sorted(by_root.keys()):
        pairs = by_root[g]
        ratios = [r for d, r in pairs]
        mean_r = np.mean(ratios)
        median_r = np.median(ratios)
        min_r = min(ratios)
        n = len(pairs)
        best_d = sorted(pairs, key=lambda x: x[1])[0][0]
        root_stats[g] = {
            'count': n, 'mean': float(mean_r), 'median': float(median_r),
            'min': float(min_r), 'best_d': best_d,
            'top_dims': [d for d,r in sorted(pairs, key=lambda x:x[1])[:5]]
        }
        adv = (1-median_r)*100
        marker = " 🏆" if mean_r < 0.90 else ""
        print(f"  root={g:>4d}: n={n:>3d}  mean={mean_r:.4f}  median={median_r:.4f}  best_ratio={min_r:.4f}  best_d={best_d}{marker}")
    
    for g in [2, 3, 5]:
        if g in root_stats:
            s = root_stats[g]
            print(f"  → root={g} has {s['count']} primes, median advantage {(1-s['median'])*100:.1f}%")
    
    return root_stats

# ═══════════════════════════════════════
# E2: MULTIPLICATIVE SUBGROUP ANALYSIS
# ═══════════════════════════════════════
def factor_subgroup_order(p):
    """Compute d-1 factorization and subgroup ladder."""
    phi = p - 1
    # Find all divisors of p-1
    divisors = []
    for i in range(1, int(math.isqrt(phi)) + 1):
        if phi % i == 0:
            divisors.append(i)
            if i != phi // i:
                divisors.append(phi // i)
    divisors.sort()
    return phi, divisors

def e2_subgroup_analysis(results_json_path):
    """Analyze how subgroup structure affects ratio."""
    with open(results_json_path) as f:
        data = json.load(f)
    
    print("\n" + "═"*60)
    print("E2: MULTIPLICATIVE SUBGROUP ANALYSIS")
    print("═"*60)
    
    # For each prime, extract: p-1 factorization, largest prime factor,
    # number of small (≤11) prime factors in p-1
    subgroup_data = []
    
    for k, v in data.items():
        d = v['d']
        best = min(r['ratio_1loc'] for r in v['runs'])
        phi = d - 1
        
        # Factor phi
        n = phi
        prime_factors = []
        p2 = 2
        while p2 * p2 <= n:
            while n % p2 == 0:
                prime_factors.append(p2)
                n //= p2
            p2 += 1
        if n > 1:
            prime_factors.append(n)
        
        largest_pf = max(prime_factors) if prime_factors else 1
        n_small_pf = sum(1 for f in prime_factors if f <= 11)
        n_large_pf = sum(1 for f in prime_factors if f > 11)
        n_unique_pf = len(set(prime_factors))
        
        # Smoothness: product of all prime factors ≤ 11
        smooth_part = 1
        for f in prime_factors:
            if f <= 11:
                smooth_part *= f
        while smooth_part % (f if smooth_part > 1 else 1) == 0:
            pass
        
        subgroup_data.append({
            'd': d, 'ratio': best,
            'phi': phi, 'prime_factors': prime_factors,
            'largest_pf': largest_pf,
            'n_small_pf': n_small_pf,
            'n_large_pf': n_large_pf,
            'n_unique_pf': n_unique_pf
        })
    
    # Group by largest prime factor
    by_lpf = defaultdict(list)
    for sd in subgroup_data:
        by_lpf[sd['largest_pf']].append(sd['ratio'])
    
    print("\n  Ratio vs Largest Prime Factor of (d-1):")
    print(f"  {'LPF':>5s} | {'count':>5s} | {'mean':>8s} | {'median':>8s} | {'min':>8s}")
    print(f"  {'-'*5} | {'-'*5} | {'-'*8} | {'-'*8} | {'-'*8}")
    for lpf in sorted(by_lpf.keys())[:10]:
        rs = by_lpf[lpf]
        print(f"  {lpf:>5d} | {len(rs):>5d} | {np.mean(rs):>8.4f} | {np.median(rs):>8.4f} | {min(rs):>8.4f}")
    
    # Correlation: ratio vs #unique prime factors
    # Manual Spearman rank correlation (no scipy dependency)
    unique_pfs = [sd['n_unique_pf'] for sd in subgroup_data]
    ratios = [sd['ratio'] for sd in subgroup_data]
    
    def spearman_r(x, y):
        n = len(x)
        rx = {v: i+1 for i, v in enumerate(sorted(set(x)))}
        ry = {v: i+1 for i, v in enumerate(sorted(set(y)))}
        rx_arr = [rx[xi] for xi in x]
        ry_arr = [ry[yi] for yi in y]
        d2 = sum((a-b)**2 for a,b in zip(rx_arr, ry_arr))
        return 1 - 6*d2/(n*(n**2-1))
    
    corr = spearman_r(unique_pfs, ratios)
    print(f"\n  Spearman ρ(unique prime factors, ratio) = {corr:.4f}")
    
    # Smooth vs wild: sqrt(unique PF count) 
    # Hypothesis: more unique prime factors → more subgroup structure → better deterministic advantage
    by_uniq = defaultdict(list)
    for sd in subgroup_data:
        by_uniq[sd['n_unique_pf']].append(sd['ratio'])
    
    print(f"\n  Ratio vs Number of Unique Prime Factors:")
    for upf in sorted(by_uniq.keys()):
        rs = by_uniq[upf]
        print(f"    unique_pf={upf}: n={len(rs)} mean={np.mean(rs):.4f} median={np.median(rs):.4f}")
    
    return subgroup_data

# ═══════════════════════════════════════
# E3: COMPOSITE DIMENSIONS (K&S)
# ═══════════════════════════════════════
def e3_composite_dimensions():
    """
    Test composite dimensions d=4,6,8,9,10,12,14,15,16
    For composite d, we need ψ(d)=LCM of pi^ei-1 directions (Kingston & Svalbe).
    Simplified: use p+1 directions for each prime factor and CRT-combine,
    or just use all p+1 MUBs for each prime power factor.
    """
    print("\n" + "═"*60)
    print("E3: COMPOSITE DIMENSIONS")
    print("═"*60)
    
    # For composite dimensions, construct MUBs using tensor products
    # For d = p1^e1 * p2^e2 * ..., we need MUBs of the tensor product
    # If each factor has m_i MUBs, total: min(m_i) complete MUBs
    
    import itertools
    
    def build_mubs_prime(d):
        omega = np.exp(2j * np.pi / d)
        mubs = [list(np.eye(d, dtype=complex))]
        for m in range(1, d+1):
            basis = []
            for k in range(d):
                v = np.array([omega**(m*j*(j-1)//2 + k*j) for j in range(d)], dtype=complex)
                basis.append(v / np.sqrt(d))
            mubs.append(basis)
        return mubs
    
    def tensor_mubs(d1, d2, mubs1, mubs2):
        """Tensor product of two MUB sets."""
        n1, n2 = len(mubs1), len(mubs2)
        result = []
        for i in range(min(n1, n2)):
            b1, b2 = mubs1[i], mubs2[i]
            basis = []
            for v1 in b1:
                for v2 in b2:
                    basis.append(np.kron(v1, v2))
            result.append(basis)
        return result
    
    composite_dims = [4, 6, 8, 9, 10, 12, 14, 15, 16]
    results_composite = []
    
    for d in composite_dims:
        # Factor d
        n = d
        factors = []
        p = 2
        while p * p <= n:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            if e > 0:
                factors.append((p, e))
            p += 1
        if n > 1:
            factors.append((1, 1) if n == 1 else (n, 1))
        
        # Build MUBs via tensor product
        mubs_list = [build_mubs_prime(p**e) for p, e in factors]
        mubs = mubs_list[0]
        for next_mubs in mubs_list[1:]:
            mubs = [
                [np.kron(v1, v2) for v1, v2 in zip(b1, b2)]
                for b1, b2 in zip(mubs, next_mubs)
            ]
        mubs = mubs[:min(len(m) for m in mubs_list)]
        n_mubs = len(mubs)
        
        n_states = 5
        r_val = 4
        ratios_1loc = []
        
        for si in range(n_states):
            rng_s = np.random.RandomState(d * 1000 + si * 10)
            psi = rng_s.randn(d) + 1j * rng_s.randn(d)
            psi /= np.linalg.norm(psi)
            rho = np.outer(psi, psi.conj())
            
            # 1-local observables: nearest-neighbor
            obs_list = []
            c = 0
            for a in range(d):
                if c >= 8: break
                for b in range(a+1, d):
                    if c >= 8: break
                    obs_list.append((a,b)); c += 1
            truth = [2*np.real(np.conj(psi[a])*psi[b]) for a,b in obs_list]
            
            n_per = r_val
            n_total = n_mubs * r_val
            
            # Deterministic
            det_snaps = []
            for mi in range(n_mubs):
                basis = mubs[mi]
                probs = np.array([np.abs(np.vdot(v, psi))**2 for v in basis])
                probs = np.clip(probs, 0, None); probs /= probs.sum()
                outcomes = np.random.choice(d, n_per, p=probs)
                for k in outcomes:
                    det_snaps.append(basis[k])
            
            # Random
            rand_snaps = []
            for _ in range(n_total):
                mi = np.random.randint(n_mubs)
                basis = mubs[mi]
                probs = np.array([np.abs(np.vdot(v, psi))**2 for v in basis])
                probs = np.clip(probs, 0, None); probs /= probs.sum()
                k = np.random.choice(d, p=probs)
                rand_snaps.append(basis[k])
            
            est_d = np.mean([(d+1)*np.outer(s, s.conj()) for s in det_snaps], axis=0) - np.eye(d)
            est_r = np.mean([(d+1)*np.outer(s, s.conj()) for s in rand_snaps], axis=0) - np.eye(d)
            
            mae_d = np.mean([abs(np.real(np.trace((np.outer(np.eye(d)[a], np.eye(d)[b]) + np.outer(np.eye(d)[b], np.eye(d)[a])) @ est_d)/np.sqrt(2)) - t) for (a,b),t in zip(obs_list, truth)])
            mae_r = np.mean([abs(np.real(np.trace((np.outer(np.eye(d)[a], np.eye(d)[b]) + np.outer(np.eye(d)[b], np.eye(d)[a])) @ est_r)/np.sqrt(2)) - t) for (a,b),t in zip(obs_list, truth)])
            
            if mae_r > 0:
                ratios_1loc.append(mae_d / mae_r)
        
        mean_r = np.mean(ratios_1loc)
        factor_str = "×".join(f"{p}^{e}" if e>1 else str(p) for p,e in factors)
        print(f"  d={d:>3d} = {factor_str:>8s}  n_mubs={n_mubs:>2d}  ratio(1loc)={mean_r:.4f}  adv={(1-mean_r)*100:>5.1f}%")
        results_composite.append({
            'd': d, 'factors': factors, 'n_mubs': n_mubs,
            'ratio_1loc_mean': float(mean_r)
        })
    
    return results_composite

# ═══════════════════════════════════════
# E4: 2-LOCAL OBSERVABLES
# ═══════════════════════════════════════
def e4_twolocal(results_json_path):
    """Test 2-local observables for d≤100 primes."""
    with open(results_json_path) as f:
        data = json.load(f)
    
    print("\n" + "═"*60)
    print("E4: 2-LOCAL OBSERVABLES (d≤100)")
    print("═"*60)
    
    primes_small = sorted([v['d'] for v in data.values() if v['d'] <= 100])
    
    results_2loc = []
    
    for d in primes_small:
        # Only every other prime to save time
        omega = np.exp(2j * np.pi / d)
        n_mubs = d + 1
        r_val = 4
        n_per = r_val
        n_total = n_mubs * r_val
        
        # Build MUBs
        mubs_list = []
        mubs_list.append([np.eye(d)[j] for j in range(d)])
        for m in range(1, d+1):
            basis = []
            for k in range(d):
                v = np.array([omega**(m*j*(j-1)//2 + k*j) for j in range(d)], dtype=complex)
                basis.append(v / np.sqrt(d))
            mubs_list.append(basis)
        
        n_states = 3
        ratios_2loc = []
        
        for si in range(n_states):
            rng_s = np.random.RandomState(d * 1000 + si * 10)
            psi = rng_s.randn(d) + 1j * rng_s.randn(d)
            psi /= np.linalg.norm(psi)
            
            # 2-local observables: |a⟩⟨b|⊗|c⟩⟨d| style (2 distinct pairs)
            obs_2loc = []
            for i in range(min(8, d)):
                a = i
                b = (i + 1) % d
                c = (i + 2) % d
                d2 = (i + 3) % d
                # 2-local: |a⟩⟨b| + |b⟩⟨a|  tensor  |c⟩⟨d| + |d⟩⟨c|
                # Simplified as observable in d-dim space
                O = np.zeros((d,d), dtype=complex)
                O[a,b] = 1.0; O[b,a] = 1.0
                O[c,d2] = 1.0; O[d2,c] = 1.0
                O /= 2.0
                obs_2loc.append(O)
            
            truth = [np.real(np.trace(O @ np.outer(psi, psi.conj()))) for O in obs_2loc]
            
            # Deterministic
            det_rho = np.zeros((d,d), dtype=complex)
            for mi in range(n_mubs):
                basis = mubs_list[mi]
                probs = np.array([np.abs(np.vdot(v, psi))**2 for v in basis])
                probs = np.clip(probs, 0, None); probs /= probs.sum()
                rng_d = np.random.RandomState(d*1000 + si*10 + mi)
                outcomes = rng_d.choice(d, n_per, p=probs)
                for k in outcomes:
                    det_rho += np.outer(basis[k], basis[k].conj())
            det_rho = det_rho * (d+1) / n_total - np.eye(d)
            
            # Random
            rand_rho = np.zeros((d,d), dtype=complex)
            rng_r = np.random.RandomState(d*2000 + si*10)
            for _ in range(n_total):
                mi = rng_r.randint(n_mubs)
                basis = mubs_list[mi]
                probs = np.array([np.abs(np.vdot(v, psi))**2 for v in basis])
                probs = np.clip(probs, 0, None); probs /= probs.sum()
                k = rng_r.choice(d, p=probs)
                rand_rho += np.outer(basis[k], basis[k].conj())
            rand_rho = rand_rho * (d+1) / n_total - np.eye(d)
            
            mae_d = np.mean([abs(np.real(np.trace(O @ det_rho)) - t) for O, t in zip(obs_2loc, truth)])
            mae_r = np.mean([abs(np.real(np.trace(O @ rand_rho)) - t) for O, t in zip(obs_2loc, truth)])
            if mae_r > 0:
                ratios_2loc.append(mae_d / mae_r)
        
        mean_r = np.mean(ratios_2loc)
        adv = (1-mean_r)*100
        print(f"  d={d:>3d}  ratio(2loc)={mean_r:.4f}  adv={adv:>5.1f}% {'★' if mean_r<0.95 else ''}")
        results_2loc.append({'d': d, 'ratio_2loc_mean': float(mean_r)})
    
    return results_2loc

# ═══════════════════════════════════════
# E5: DEPOLARIZING NOISE
# ═══════════════════════════════════════
def e5_depolarizing_noise():
    """Test DPRT vs CS under depolarizing noise for selected primes."""
    print("\n" + "═"*60)
    print("E5: DEPOLARIZING NOISE ROBUSTNESS")
    print("═"*60)
    
    test_dims = [3, 5, 7, 11, 13, 17, 19, 23, 31, 37, 47, 61, 97]
    noise_levels = [0.0, 0.01, 0.05, 0.1, 0.2, 0.3]
    
    results_noise = []
    
    for d in test_dims:
        omega = np.exp(2j * np.pi / d)
        n_mubs = d + 1
        r_val = 4
        n_per = r_val
        n_total = n_mubs * r_val
        
        mubs_list = []
        mubs_list.append([np.eye(d)[j] for j in range(d)])
        for m in range(1, d+1):
            basis = []
            for k in range(d):
                v = np.array([omega**(m*j*(j-1)//2 + k*j) for j in range(d)], dtype=complex)
                basis.append(v / np.sqrt(d))
            mubs_list.append(basis)
        
        n_states = 8
        print(f"\n  d={d}:")
        
        for noise in noise_levels:
            ratios_1loc = []
            
            for si in range(n_states):
                rng_s = np.random.RandomState(d * 10000 + int(noise*1000) + si)
                psi = rng_s.randn(d) + 1j * rng_s.randn(d)
                psi /= np.linalg.norm(psi)
                rho_pure = np.outer(psi, psi.conj())
                # Apply depolarizing noise: ρ → (1-λ)ρ + λI/d
                rho = (1 - noise) * rho_pure + noise * np.eye(d, dtype=complex) / d
                
                # Truth values (from noisy state!)
                obs_list = []
                c = 0
                for a in range(d):
                    if c >= 8: break
                    for b in range(a+1, d):
                        if c >= 8: break
                        obs_list.append((a,b)); c += 1
                truth = [2 * np.real(rho[a,b]) for a,b in obs_list]
                
                # Deterministic
                det_rho_est = np.zeros((d,d), dtype=complex)
                for mi in range(n_mubs):
                    basis = mubs_list[mi]
                    probs = np.array([np.real(np.trace(rho @ np.outer(v, v.conj()))) for v in basis])
                    probs = np.clip(probs, 0, None); probs /= probs.sum()
                    rng_d = np.random.RandomState(d * 20000 + int(noise*1000) + si + mi)
                    outcomes = rng_d.choice(d, n_per, p=probs)
                    for k in outcomes:
                        det_rho_est += np.outer(basis[k], basis[k].conj())
                det_rho_est = det_rho_est * (d+1) / n_total - np.eye(d)
                
                # Random
                rand_rho_est = np.zeros((d,d), dtype=complex)
                rng_r = np.random.RandomState(d * 30000 + int(noise*1000) + si)
                for _ in range(n_total):
                    mi = rng_r.randint(n_mubs)
                    basis = mubs_list[mi]
                    probs = np.array([np.real(np.trace(rho @ np.outer(v, v.conj()))) for v in basis])
                    probs = np.clip(probs, 0, None); probs /= probs.sum()
                    k = rng_r.choice(d, p=probs)
                    rand_rho_est += np.outer(basis[k], basis[k].conj())
                rand_rho_est = rand_rho_est * (d+1) / n_total - np.eye(d)
                
                mae_d = np.mean([abs(2*np.real(det_rho_est[a,b]) - t) for (a,b),t in zip(obs_list,truth)])
                mae_r = np.mean([abs(2*np.real(rand_rho_est[a,b]) - t) for (a,b),t in zip(obs_list,truth)])
                if mae_r > 0:
                    ratios_1loc.append(mae_d / mae_r)
            
            mean_r = np.mean(ratios_1loc)
            adv = (1-mean_r)*100
            marker = " ★" if mean_r<0.90 else (" ✦" if mean_r<0.95 else "")
            print(f"    λ={noise:5.2f}  ratio={mean_r:.4f}  DPRT+{adv:>5.1f}%{marker}")
            results_noise.append({
                'd': d, 'noise': noise, 'ratio_1loc_mean': float(mean_r)
            })
    
    return results_noise

# ═══════════════════════════════════════
# MAIN
# ═══════════════════════════════════════
if __name__ == '__main__':
    T0 = time.time()
    
    # E1 + E2 use existing data
    e1 = e1_primitive_root_analysis('primes_1000_results.json')
    e2 = e2_subgroup_analysis('primes_1000_results.json')
    
    # E3: Composite dimensions
    e3 = e3_composite_dimensions()
    
    # E4: 2-local
    e4 = e4_twolocal('primes_1000_results.json')
    
    # E5: Depolarizing noise
    e5 = e5_depolarizing_noise()
    
    elapsed = time.time() - T0
    
    # Summary
    print(f"\n{'═'*60}")
    print(f"ALL 5 EXPERIMENTS COMPLETE — {elapsed:.0f}s")
    print(f"{'═'*60}")
    
    # E5 summary
    print(f"\n── E5 Noise Summary ──")
    noise_summary = defaultdict(list)
    for r in e5:
        noise_summary[r['noise']].append(r['ratio_1loc_mean'])
    for noise in sorted(noise_summary.keys()):
        rs = noise_summary[noise]
        print(f"  λ={noise:.2f}: mean_ratio={np.mean(rs):.4f}  median={np.median(rs):.4f}  min={min(rs):.4f}")
    
    # Save
    all_results = {
        'e1_root_analysis': e1,
        'e2_subgroup': {'n_points': len(e2)},
        'e3_composite': e3,
        'e4_2local': e4,
        'e5_noise': e5
    }
    
    # Convert non-serializable
    import copy
    serializable = {}
    for k, v in all_results.items():
        if k == 'e1_root_analysis':
            serializable[k] = {str(g): s for g, s in v.items()}
        elif k == 'e2_subgroup':
            serializable[k] = v
        else:
            serializable[k] = v
    
    with open('all5_experiments_results.json', 'w') as f:
        json.dump(serializable, f, indent=2, default=str)
    
    print(f"\n→ all5_experiments_results.json saved", flush=True)
