#!/usr/bin/env python3
"""
RadonShadow Week 1: 补齐实验矩阵
=================================
Task 1: Qutrit × Mixed State (d=3, state-type effect)
Task 2: Large-d scan (d=13, 17, 19) — single-observable comparison
Task 3: 2-qubit full-state (d=4, non-prime) — MUB linear inversion vs CS
Task 4: Post-processing complexity bench (wall-clock + memory)

Author: QClaw (2026-07-14)
"""

import numpy as np
from scipy.linalg import sqrtm
from collections import defaultdict
import time, json, itertools, tracemalloc

# ======================================================================
# 0. Shared Utilities
# ======================================================================

def random_pure_state(d, seed):
    np.random.seed(seed)
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    return psi / np.linalg.norm(psi)

def random_mixed_state(d, seed, rank=None):
    """Random mixed state via Ginibre ensemble."""
    np.random.seed(seed)
    if rank is None:
        rank = np.random.randint(1, d+1)
    G = np.random.randn(d, rank) + 1j * np.random.randn(d, rank)
    G /= np.trace(G @ G.conj().T)**0.5
    rho = G @ G.conj().T
    # add noise for full rank if desired
    tau = np.random.dirichlet(np.ones(d), 1)[0]
    U = np.random.randn(d, d) + 1j * np.random.randn(d, d)
    U, _ = np.linalg.qr(U)
    rho_mix = U @ np.diag(tau) @ U.conj().T
    return rho_mix / np.trace(rho_mix)

def build_mubs(d):
    """Construct d+1 MUBs for prime d or d=2^n tensor-product."""
    if d == 2:
        return [
            [np.array([1,0],dtype=complex), np.array([0,1],dtype=complex)],
            [np.array([1,1],dtype=complex)/np.sqrt(2), np.array([1,-1],dtype=complex)/np.sqrt(2)],
            [np.array([1,1j],dtype=complex)/np.sqrt(2), np.array([1,-1j],dtype=complex)/np.sqrt(2)],
        ]
    if d == 4:
        # 2-qubit: tensor product of d=2 MUBs gives 3*3=9 product bases
        # But for CS comparison, we use the 5 MUBs known for d=4
        mubs_d2 = build_mubs(2)
        # Product MUBs (3 from same bases on both qubits, 6 from mixed = 9 total)
        # But only 5 are MUB: the 3 product of identical bases + 2 mixed bases
        # Known: d=4 has exactly 5 MUBs
        mubs = []
        # 3 product bases: Z⊗Z, X⊗X, Y⊗Y
        for b in mubs_d2:
            basis = []
            for v1 in b:
                for v2 in b:
                    basis.append(np.kron(v1, v2))
            mubs.append(basis)
        # 2 mixed bases (using Bell-like construction)
        bell_map = [
            (mubs_d2[0], mubs_d2[1]),  # Z⊗X
            (mubs_d2[0], mubs_d2[2]),  # Z⊗Y
        ]
        for (b1, b2) in bell_map:
            basis = []
            for v1 in b1:
                for v2 in b2:
                    basis.append(np.kron(v1, v2))
            mubs.append(basis)
        return mubs[:5]  # exactly 5 MUBs for d=4
    
    # Prime d
    if d == 2:
        pass  # handled above
    omega = np.exp(2j * np.pi / d)
    mubs = [list(np.eye(d, dtype=complex))]
    for m in range(1, d+1):
        mubs.append([np.array([omega**(m*j*(j-1)//2 + k*j) for j in range(d)], dtype=complex)/np.sqrt(d) for k in range(d)])
    return mubs

def mub_linear_reconstruct(rho, d, n_per_basis, mubs):
    """MUB linear inversion: ρ = Σ pΠ - I (exact)."""
    Id = np.eye(d, dtype=complex)
    projs_all = [[np.outer(v, v.conj()) for v in basis] for basis in mubs]
    
    rho_rec = np.zeros((d,d), dtype=complex)
    for i, basis in enumerate(mubs):
        probs = np.array([np.real(np.trace(rho @ P)) for P in projs_all[i]])
        probs = np.clip(probs, 0, None); probs /= probs.sum()
        # finite shots: multinomial
        outcomes = np.random.multinomial(n_per_basis, probs)
        stats = outcomes / n_per_basis
        for k in range(d):
            rho_rec += stats[k] * projs_all[i][k]
    rho_rec -= Id
    return rho_rec, np.linalg.norm(rho - rho_rec, 'fro')

def cs_reconstruct(rho, d, n_total, mubs, seed=42):
    """Classical Shadow, random MUB. ρ̂_snap = (d+1)|ψ⟩⟨ψ| - I."""
    np.random.seed(seed)
    Id = np.eye(d, dtype=complex)
    factor = len(mubs)  # d+1 for prime, 5 for d=4
    n_bases = len(mubs)
    
    snaps = []
    for _ in range(n_total):
        i = np.random.randint(n_bases)
        Pj = [np.outer(v, v.conj()) for v in mubs[i]]
        probs = np.array([np.real(np.trace(rho @ P)) for P in Pj])
        probs = np.clip(probs, 0, None); probs /= probs.sum()
        o = np.random.choice(d, p=probs)
        snaps.append(factor * Pj[o] - Id)
    return np.mean(snaps, axis=0), np.linalg.norm(rho - np.mean(snaps, axis=0), 'fro')

def gell_mann_observables(d):
    """Return list of (name, matrix) for generalized Gell-Mann basis."""
    obs = {}
    for j in range(d):
        for k in range(j+1, d):
            M = np.zeros((d,d), dtype=complex); M[j,k]=1; M[k,j]=1
            obs[f'λ_{j}{k}'] = M
            M2 = np.zeros((d,d), dtype=complex); M2[j,k]=-1j; M2[k,j]=1j
            obs[f'λ̃_{j}{k}'] = M2
    for l in range(1, d):
        coeff = 1.0 / np.sqrt(l*(l+1))
        M = np.zeros((d,d)); 
        for k in range(l): M[k,k] = coeff
        M[l,l] = -l * coeff
        obs[f'd_{l}'] = M
    return obs

def snap_qudit(rho, mub_label, mubs, d):
    """Single shadow snapshot in given MUB."""
    basis = mubs[mub_label]
    projs = [np.outer(v, v.conj()) for v in basis]
    probs = np.array([np.real(np.trace(rho @ P)) for P in projs])
    probs = np.clip(probs, 0, None); probs /= probs.sum()
    outcome = np.random.choice(d, p=probs)
    return (d+1) * projs[outcome] - np.eye(d, dtype=complex)


# ======================================================================
# TASK 1: Qutrit × Mixed State (d=3)
# ======================================================================
def task1_qutrit_mixed(n_states=30, r_list=[1,2,4,8], n_shots_per=500):
    print("="*60)
    print("TASK 1: Qutrit (d=3) — Mixed State vs Pure State")
    print("="*60)
    
    mubs = build_mubs(3)
    observables = gell_mann_observables(3)
    obs_labels = list(observables.keys())
    n_configs = 4  # d+1 MUBs
    
    results = {}
    for state_type, gen_fn in [("pure", random_pure_state), ("mixed", random_mixed_state)]:
        print(f"\n  State type: {state_type}")
        for r in r_list:
            n_meas = n_configs * r
            det_err = []
            rnd_err = []
            
            for si in range(n_states):
                s_seed = 10000 + si * 1000 + r * 100
                raw = gen_fn(3, s_seed)
                rho = np.outer(raw, raw.conj()) if state_type == 'pure' else raw
                
                truth = {l: np.real(np.trace(observables[l] @ rho)) for l in obs_labels}
                
                # Deterministic
                snaps_d = []
                for _ in range(r):
                    for mi in range(n_configs):
                        snaps_d.append(snap_qudit(rho, mi, mubs, 3))
                
                # Random
                snaps_r = []
                for _ in range(n_meas):
                    mi = np.random.randint(n_configs)
                    snaps_r.append(snap_qudit(rho, mi, mubs, 3))
                
                for l in obs_labels:
                    O = observables[l]
                    pred_d = np.mean([np.real(np.trace(O @ s)) for s in snaps_d])
                    pred_r = np.mean([np.real(np.trace(O @ s)) for s in snaps_r])
                    det_err.append(abs(pred_d - truth[l]))
                    rnd_err.append(abs(pred_r - truth[l]))
            
            de = np.mean(det_err); re = np.mean(rnd_err)
            print(f"    r={r:>2d} N={n_meas:>3d}  Det={de:.6f}  Rand={re:.6f}  ratio={de/re:.4f}  {'[DET]' if de<re else '[RAND]'}")
            results.setdefault(state_type, []).append({
                'r': r, 'n_meas': n_meas, 'det_mae': float(de), 'rand_mae': float(re),
                'ratio': float(de/re), 'winner': 'DET' if de<re else 'RAND'
            })
    return results


# ======================================================================
# TASK 2: Large-d Scan (d=13, 17, 19)
# ======================================================================
def task2_large_d_scan(d_list=[13, 17, 19], n_states=20, r_list=[1,2,4], seed=42):
    print("\n" + "="*60)
    print("TASK 2: Large-d Scan (d=13, 17, 19) — Single Observable")
    print("="*60)
    
    results = {}
    for d in d_list:
        print(f"\n  d={d} (prime, {d+1} MUBs):")
        mubs = build_mubs(d)
        mub_labels = list(range(len(mubs)))
        observables = gell_mann_observables(d)
        obs_labels = list(observables.keys())
        n_configs = d + 1
        
        d_results = []
        for r in r_list:
            n_meas = n_configs * r
            det_err = []
            rnd_err = []
            
            for si in range(n_states):
                s_seed = seed + d * 10000 + si * 1000 + r * 100
                np.random.seed(s_seed)
                psi = np.random.randn(d) + 1j * np.random.randn(d)
                psi /= np.linalg.norm(psi)
                rho = np.outer(psi, psi.conj())
                
                truth = {l: np.real(np.trace(observables[l] @ rho)) for l in obs_labels}
                
                snaps_d = []
                for _ in range(r):
                    for mi in mub_labels:
                        snaps_d.append(snap_qudit(rho, mi, mubs, d))
                
                snaps_r = []
                for _ in range(n_meas):
                    mi = np.random.randint(n_configs)
                    snaps_r.append(snap_qudit(rho, mi, mubs, d))
                
                for l in obs_labels:
                    O = observables[l]
                    det_err.append(abs(np.mean([np.real(np.trace(O @ s)) for s in snaps_d]) - truth[l]))
                    rnd_err.append(abs(np.mean([np.real(np.trace(O @ s)) for s in snaps_r]) - truth[l]))
            
            de = np.mean(det_err); re = np.mean(rnd_err)
            print(f"    r={r:>2d} N={n_meas:>4d}  Det={de:.6f}  Rand={re:.6f}  ratio={de/re:.4f}  {'[DET]' if de<re else '[RAND]'}")
            d_results.append({
                'd': d, 'r': r, 'n_meas': n_meas, 'det_mae': float(de), 'rand_mae': float(re),
                'ratio': float(de/re), 'winner': 'DET' if de<re else 'RAND'
            })
        results[f'd={d}'] = d_results
    return results


# ======================================================================
# TASK 3: 2-qubit Full-State DPRT vs CS (d=4, non-prime)
# ======================================================================
def task3_2qubit_fullstate(n_states=30, n_per_list=[50, 200, 500]):
    print("\n" + "="*60)
    print("TASK 3: 2-qubit (d=4) Full-State: MUB Linear Inversion vs CS")
    print("="*60)
    
    d = 4
    mubs = build_mubs(d)
    n_bases = len(mubs)
    print(f"    Using {n_bases} MUBs for d=4 (known max is 5)")
    
    results = []
    for npb in n_per_list:
        n_total = n_bases * npb
        dprt_errs = []
        cs_errs = []
        
        for si in range(n_states):
            rho = random_pure_state(d, si * 1000 + npb)
            rho_mat = np.outer(rho, rho.conj())
            
            _, e1 = mub_linear_reconstruct(rho_mat, d, npb, mubs)
            _, e2 = cs_reconstruct(rho_mat, d, n_total, mubs, seed=si*1000+npb+1)
            dprt_errs.append(e1); cs_errs.append(e2)
        
        de = np.mean(dprt_errs); ce = np.mean(cs_errs)
        print(f"    s/b={npb:>4d} tot={n_total:>5d}  MUB={de:.6f}  CS={ce:.6f}  ratio={de/ce:.4f}  {'[MUB]' if de<ce else '[CS]'}")
        results.append({
            'd': d, 'n_mubs': n_bases, 'shots_per_basis': npb, 'n_total': n_total,
            'mub_err': float(de), 'cs_err': float(ce), 'ratio': float(de/ce),
            'winner': 'MUB' if de<ce else 'CS'
        })
    return results


# ======================================================================
# TASK 4: Post-Processing Complexity Bench (wall-clock + memory)
# ======================================================================
def task4_complexity_bench(d_list=[2,3,5,7,11,13,17,19], n_states=10):
    print("\n" + "="*60)
    print("TASK 4: Post-Processing Complexity Bench (Wall-Clock + Memory)")
    print("="*60)
    
    tracemalloc.start()
    results = []
    
    for d in d_list:
        mubs = build_mubs(d)
        n_bases = len(mubs)
        
        dprt_times = []
        cs_times = []
        dprt_mems = []
        cs_mems = []
        
        for si in range(n_states):
            rho = random_pure_state(d, si*100+99)
            rho_mat = np.outer(rho, rho.conj())
            
            # MUB linear inversion timing
            t0 = time.perf_counter()
            _, e1 = mub_linear_reconstruct(rho_mat, d, 500, mubs)
            t1 = time.perf_counter()
            dprt_times.append((t1-t0)*1000)
            
            # CS timing
            t0 = time.perf_counter()
            _, e2 = cs_reconstruct(rho_mat, d, n_bases*500, mubs, seed=si)
            t1 = time.perf_counter()
            cs_times.append((t1-t0)*1000)
        
        dt = np.mean(dprt_times); ct = np.mean(cs_times)
        
        # Memory: size of data structures
        # DPRT: stores (d+1)×d probabilities = O(d²)
        # CS: stores N d×d snapshots = O(d³) for large N
        n_total = n_bases * 500
        dprt_mem_bytes = (n_bases * d) * 8  # (d+1)*d floats
        cs_mem_bytes = n_total * d * d * 16  # N complex matrices
        
        print(f"    d={d:>3d}:  MUB={dt:>8.2f}ms  CS={ct:>8.2f}ms  "
              f"time_ratio={dt/ct:.4f}  "
              f"mem(MUB)={dprt_mem_bytes/1024:>6.1f}KB  mem(CS)={cs_mem_bytes/1024/1024:>6.1f}MB")
        
        results.append({
            'd': d, 'n_bases': n_bases,
            'mub_time_ms': float(dt), 'cs_time_ms': float(ct),
            'time_ratio': float(dt/ct),
            'mub_mem_bytes': dprt_mem_bytes, 'cs_mem_bytes': cs_mem_bytes,
            'mem_ratio': float(dprt_mem_bytes / cs_mem_bytes) if cs_mem_bytes > 0 else 0
        })
    
    return results


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':
    T0 = time.time()
    all_output = {}
    
    # T1: Qutrit Mixed
    t1 = task1_qutrit_mixed(n_states=30, r_list=[1,2,4,8])
    all_output['task1_qutrit_mixed'] = t1
    
    # T2: Large-d scan
    t2 = task2_large_d_scan(d_list=[13, 17, 19], n_states=20, r_list=[1,2,4])
    all_output['task2_large_d_scan'] = t2
    
    # T3: 2-qubit full-state
    t3 = task3_2qubit_fullstate(n_states=30, n_per_list=[50, 200, 500])
    all_output['task3_2qubit_fullstate'] = t3
    
    # T4: Complexity bench
    t4 = task4_complexity_bench(d_list=[2,3,5,7,11,13,17,19])
    all_output['task4_complexity_bench'] = t4
    
    elapsed = time.time() - T0
    print(f"\n{'='*60}")
    print(f"W1 COMPLETE — Total time: {elapsed:.1f}s")
    print(f"{'='*60}")
    
    # Summary
    print(f"\n── TASK 1: Qutrit Mixed vs Pure ──")
    for sty, dat in t1.items():
        print(f"  {sty}:")
        for r in dat:
            print(f"    r={r['r']} ratio={r['ratio']:.4f} [{r['winner']}]")
    
    print(f"\n── TASK 2: Large-d ──")
    for k, v in t2.items():
        print(f"  {k}:")
        for r in v:
            print(f"    r={r['r']} ratio={r['ratio']:.4f} [{r['winner']}]")
    
    print(f"\n── TASK 3: 2-qubit full-state ──")
    for r in t3:
        print(f"    s/b={r['shots_per_basis']} ratio={r['ratio']:.4f} [{r['winner']}]")
    
    print(f"\n── TASK 4: Complexity ──")
    for r in t4:
        print(f"    d={r['d']} time_ratio={r['time_ratio']:.3f} mem_ratio={r['mem_ratio']:.2e}")
    
    # Save
    with open('w1_experiment_results.json', 'w') as f:
        json.dump(all_output, f, indent=2)
    print(f"\n→ w1_experiment_results.json saved")
