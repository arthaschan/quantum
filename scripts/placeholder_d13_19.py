#!/usr/bin/env python3
"""
RadonShadow: d=13,17,19 占坑实验
=================================
- Full-state reconstruction DPRT vs CS head-to-head
- k-locality observable breakdown
- Deterministic vs random at multiple shot counts
"""

import numpy as np
from collections import defaultdict
import time, json

def build_mubs_prime(d):
    """d+1 MUBs for prime d."""
    omega = np.exp(2j * np.pi / d)
    mubs = [list(np.eye(d, dtype=complex))]
    for m in range(1, d+1):
        basis = []
        for k in range(d):
            v = np.array([omega**(m*j*(j-1)//2 + k*j) for j in range(d)], dtype=complex)
            basis.append(v / np.sqrt(d))
        mubs.append(basis)
    return mubs

def random_pure_state(d, seed):
    np.random.seed(seed)
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    return psi / np.linalg.norm(psi)

def snap_mub(rho, mi, mubs, d):
    basis = mubs[mi]
    projs = [np.outer(v, v.conj()) for v in basis]
    probs = np.array([np.real(np.trace(rho @ P)) for P in projs])
    probs = np.clip(probs, 0, None); probs /= probs.sum()
    outcome = np.random.choice(d, p=probs)
    return (d+1) * projs[outcome] - np.eye(d, dtype=complex)

def full_state_reconstruct(rho, d, n_per_basis, mubs):
    """DPRT/MUB linear inversion for full state."""
    Id = np.eye(d, dtype=complex)
    projs_all = [[np.outer(v, v.conj()) for v in basis] for basis in mubs]
    n_mubs = len(mubs)
    
    rho_est = np.zeros((d,d), dtype=complex)
    for i in range(n_mubs):
        probs = np.array([np.real(np.trace(rho @ P)) for P in projs_all[i]])
        probs = np.clip(probs, 0, None); probs /= probs.sum()
        outcomes = np.random.multinomial(n_per_basis, probs)
        stats = outcomes / n_per_basis
        for k in range(d):
            rho_est += stats[k] * projs_all[i][k]
    rho_est -= Id
    return np.linalg.norm(rho - rho_est, 'fro')

def cs_full_state_reconstruct(rho, d, n_total, mubs, seed):
    np.random.seed(seed)
    n_mubs = len(mubs)
    Id = np.eye(d, dtype=complex)
    snaps = []
    for _ in range(n_total):
        mi = np.random.randint(n_mubs)
        basis = mubs[mi]
        projs = [np.outer(v, v.conj()) for v in basis]
        probs = np.array([np.real(np.trace(rho @ P)) for P in projs])
        probs = np.clip(probs, 0, None); probs /= probs.sum()
        o = np.random.choice(d, p=probs)
        snaps.append((d+1) * projs[o] - Id)
    rho_est = np.mean(snaps, axis=0)
    return np.linalg.norm(rho - rho_est, 'fro')

def k_local_observables(d, k):
    """Generate k-local observables for a d-dim system.
    For qudit d, 'k-local' is simulated as k-norm of Gell-Mann matrices."""
    obs = []
    # 1-local (nearest-neighbor Gell-Mann pairs)
    count = 0
    max_obs = min(50, d*(d-1))
    for j in range(d):
        for l in range(j+1, d):
            if count >= max_obs: break
            # Symmetric
            M = np.zeros((d,d), dtype=complex)
            M[j,l] = M[l,j] = 1.0
            obs.append(('1-local', M / np.sqrt(2)))
            count += 1
    # 2-local: pair correlation
    count = 0
    for j in range(d):
        for l in range(j+1, d):
            for p in range(d):
                for q in range(p+1, d):
                    if count >= 30: break
                    if set([j,l,p,q]) != set([j,l,p,q]): pass
                    M = np.eye(d, dtype=complex)
                    M[j,j] = M[l,l] = 0
                    obs.append(('2-local', M))
                    count += 1
            if count >= 30: break
        if count >= 30: break
    # If too few 2-local, add random observables
    if len([o for t,o in obs if t=='2-local']) < 10:
        for _ in range(10):
            H = np.random.randn(d,d) + 1j * np.random.randn(d,d)
            H = (H + H.conj().T) / 2
            obs.append(('2-local', H))
    return obs

def run_placeholder(d, n_states=20, r_list=[1,2,4,8]):
    """Run full experiments for a given d."""
    print(f"\n{'='*60}")
    print(f"d={d}: Full-state + k-local obs + timing")
    print(f"{'='*60}")
    
    is_prime = d in [13,17,19]
    mubs = build_mubs_prime(d) if is_prime else None
    n_mubs = len(mubs)
    
    results = {'d': d, 'n_mubs': n_mubs, 'experiments': []}
    
    for r in r_list:
        n_per_basis = r
        n_total = n_mubs * r
        
        full_dprt = []
        full_cs = []
        k1_dprt = []
        k1_cs = []
        k2_dprt = []
        k2_cs = []
        times_dprt = []
        times_cs = []
        
        observables = k_local_observables(d, k=1)
        
        for si in range(n_states):
            s_seed = d * 10000 + r * 1000 + si * 100
            psi = random_pure_state(d, s_seed)
            rho = np.outer(psi, psi.conj())
            
            # Full-state DPRT
            t0 = time.perf_counter()
            e1 = full_state_reconstruct(rho, d, n_per_basis, mubs)
            times_dprt.append((time.perf_counter()-t0)*1000)
            full_dprt.append(e1)
            
            # Full-state CS
            t0 = time.perf_counter()
            e2 = cs_full_state_reconstruct(rho, d, n_total, mubs, s_seed+1)
            times_cs.append((time.perf_counter()-t0)*1000)
            full_cs.append(e2)
            
            # k-local observable estimation
            truth = {i: np.real(np.trace(O @ rho)) for i, (_, O) in enumerate(observables)}
            
            # Deterministic: same mubs for all observable estimation
            snaps_d = []
            for _ in range(r):
                for mi in range(n_mubs):
                    snaps_d.append(snap_mub(rho, mi, mubs, d))
            
            snaps_r = []
            for _ in range(n_total):
                mi = np.random.randint(n_mubs)
                snaps_r.append(snap_mub(rho, mi, mubs, d))
            
            for i, (ktype, O) in enumerate(observables):
                pred_d = np.mean([np.real(np.trace(O @ s)) for s in snaps_d])
                pred_r = np.mean([np.real(np.trace(O @ s)) for s in snaps_r])
                err_d = abs(pred_d - truth[i])
                err_r = abs(pred_r - truth[i])
                if ktype == '1-local':
                    k1_dprt.append(err_d); k1_cs.append(err_r)
                else:
                    k2_dprt.append(err_d); k2_cs.append(err_r)
        
        fd = np.mean(full_dprt); fc = np.mean(full_cs)
        td = np.mean(times_dprt); tc = np.mean(times_cs)
        k1d = np.mean(k1_dprt) if k1_dprt else 0; k1c = np.mean(k1_cs) if k1_cs else 0
        k2d = np.mean(k2_dprt) if k2_dprt else 0; k2c = np.mean(k2_cs) if k2_cs else 0
        
        print(f"  r={r:>2d} N={n_total:>4d}  |  Full: DPRT={fd:.5f} CS={fc:.5f} ratio={fd/fc:.4f}")
        print(f"                |  1-local: DET={k1d:.5f} RAND={k1c:.5f} ratio={k1d/k1c:.4f}")
        print(f"                |  2-local: DET={k2d:.5f} RAND={k2c:.5f} ratio={k2d/k2c:.4f}")
        print(f"                |  Time:  DPRT={td:.2f}ms CS={tc:.1f}ms  speedup={tc/td:.0f}x")
        
        results['experiments'].append({
            'r': r, 'n_per_basis': n_per_basis, 'n_total': n_total,
            'full_dprt': float(fd), 'full_cs': float(fc), 'full_ratio': float(fd/fc),
            'k1_dprt': float(k1d), 'k1_cs': float(k1c), 'k1_ratio': float(k1d/k1c) if k1c>0 else 0,
            'k2_dprt': float(k2d), 'k2_cs': float(k2c), 'k2_ratio': float(k2d/k2c) if k2c>0 else 0,
            'dprt_time_ms': float(td), 'cs_time_ms': float(tc), 'speedup': float(tc/td)
        })
    return results

if __name__ == '__main__':
    T0 = time.time()
    all_results = {}
    
    for d in [13, 17, 19]:
        r = run_placeholder(d, n_states=15, r_list=[1,2,4,8])
        all_results[f'd={d}'] = r
    
    elapsed = time.time() - T0
    print(f"\n{'='*60}")
    print(f"PLACEHOLDER EXPERIMENTS COMPLETE — {elapsed:.1f}s")
    print(f"{'='*60}")
    
    # Key summary table
    print(f"\n── Full-State Ratio Summary ──")
    print(f"{'d':>4s} {'r=1':>8s} {'r=2':>8s} {'r=4':>8s} {'r=8':>8s} {'Best':>8s}")
    for d in [13,17,19]:
        exps = all_results[f'd={d}']['experiments']
        ratios = [e['full_ratio'] for e in exps]
        best = min(ratios)
        print(f"{d:>4d} {ratios[0]:>8.4f} {ratios[1]:>8.4f} {ratios[2]:>8.4f} {ratios[3]:>8.4f} {best:>8.4f}")
    
    with open('d13_19_placeholder_results.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\n→ d13_19_placeholder_results.json saved")
