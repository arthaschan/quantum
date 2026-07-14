#!/usr/bin/env python3
"""
RadonShadow: All primes 23..97 full-state + 1-local + timing
Flush every print, use fewer states for huge d.
"""
import numpy as np, math, time, json, sys

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

def random_pure_state(d, seed):
    rng = np.random.RandomState(seed)
    psi = rng.randn(d) + 1j * rng.randn(d)
    return psi / np.linalg.norm(psi)

def full_state_reconstruct(rho, d, n_per_basis, mubs, rng):
    projs_all = [[np.outer(v, v.conj()) for v in basis] for basis in mubs]
    n_mubs = len(mubs)
    rho_est = np.zeros((d,d), dtype=complex)
    for i in range(n_mubs):
        probs = np.array([np.real(np.trace(rho @ P)) for P in projs_all[i]])
        probs = np.clip(probs, 0, None); probs /= probs.sum()
        outcomes = rng.multinomial(n_per_basis, probs)
        stats = outcomes / n_per_basis
        for k in range(d):
            rho_est += stats[k] * projs_all[i][k]
    rho_est -= np.eye(d, dtype=complex)
    return np.linalg.norm(rho - rho_est, 'fro')

def cs_full_state_reconstruct(rho, d, n_total, mubs, rng):
    n_mubs = len(mubs)
    snaps = []
    for _ in range(n_total):
        mi = rng.randint(n_mubs)
        basis = mubs[mi]
        projs = [np.outer(v, v.conj()) for v in basis]
        probs = np.array([np.real(np.trace(rho @ P)) for P in projs])
        probs = np.clip(probs, 0, None); probs /= probs.sum()
        o = rng.choice(d, p=probs)
        snaps.append((d+1) * projs[o] - np.eye(d, dtype=complex))
    return np.linalg.norm(rho - np.mean(snaps, axis=0), 'fro')

def snap_mub(rho, mi, mubs, d, rng):
    basis = mubs[mi]
    projs = [np.outer(v, v.conj()) for v in basis]
    probs = np.array([np.real(np.trace(rho @ P)) for P in projs])
    probs = np.clip(probs, 0, None); probs /= probs.sum()
    outcome = rng.choice(d, p=probs)
    return (d+1) * projs[outcome] - np.eye(d, dtype=complex)

def one_local_obs(d, max_obs=15):
    obs = []
    count = 0
    for j in range(d):
        if count >= max_obs: break
        for k in range(j+1, d):
            if count >= max_obs: break
            M = np.zeros((d,d), dtype=complex)
            M[j,k] = M[k,j] = 1.0
            obs.append(M / np.sqrt(2))
            count += 1
    return obs

# Only primes 23..97
primes = [p for p in range(23, 98) if all(p%i for i in range(2, int(math.isqrt(p))+1))]
print(f"Testing {len(primes)} primes: {primes}", flush=True)

# Adaptive state count: fewer for large d
def n_states_for(d):
    if d <= 31: return 8
    if d <= 59: return 5
    return 3

all_results = {}
T0 = time.time()

for d in primes:
    dt0 = time.perf_counter()
    n_st = n_states_for(d)
    mubs = build_mubs_prime(d)
    n_mubs = len(mubs)
    obs_list = one_local_obs(d, max_obs=15)
    
    d_result = {'d': d, 'n_mubs': n_mubs, 'n_states': n_st, 'runs': []}
    
    for r in [1, 2, 4]:
        n_per = r
        n_total = n_mubs * r
        
        full_d = []; full_c = []
        k1_d = []; k1_c = []
        time_d = []; time_c = []
        
        for si in range(n_st):
            rng_d = np.random.RandomState(d*10000 + r*100 + si*10 + 1)
            rng_c = np.random.RandomState(d*10000 + r*100 + si*10 + 2)
            rng_snap = np.random.RandomState(d*10000 + r*100 + si*10 + 3)
            
            psi = random_pure_state(d, d*1000 + r*100 + si)
            rho = np.outer(psi, psi.conj())
            
            t0 = time.perf_counter()
            e1 = full_state_reconstruct(rho, d, n_per, mubs, rng_d)
            time_d.append((time.perf_counter()-t0)*1000)
            full_d.append(e1)
            
            t0 = time.perf_counter()
            e2 = cs_full_state_reconstruct(rho, d, n_total, mubs, rng_c)
            time_c.append((time.perf_counter()-t0)*1000)
            full_c.append(e2)
            
            # 1-local
            truth = [np.real(np.trace(O @ rho)) for O in obs_list]
            snaps_d = [snap_mub(rho, mi, mubs, d, rng_snap) for _ in range(r) for mi in range(n_mubs)]
            snaps_r = [snap_mub(rho, rng_snap.randint(n_mubs), mubs, d, rng_snap) for _ in range(n_total)]
            
            for i, O in enumerate(obs_list):
                k1_d.append(abs(np.mean([np.real(np.trace(O@s)) for s in snaps_d]) - truth[i]))
                k1_c.append(abs(np.mean([np.real(np.trace(O@s)) for s in snaps_r]) - truth[i]))
        
        fd, fc = np.mean(full_d), np.mean(full_c)
        td, tc = np.mean(time_d), np.mean(time_c)
        k1d, k1c = np.mean(k1_d), np.mean(k1_c)
        
        ratio_f = fd/fc if fc>0 else 0
        ratio_1 = k1d/k1c if k1c>0 else 0
        
        d_result['runs'].append({
            'r': r, 'n_total': n_total,
            'full_dprt': float(fd), 'full_cs': float(fc), 'full_ratio': float(ratio_f),
            'k1_dprt': float(k1d), 'k1_cs': float(k1c), 'k1_ratio': float(ratio_1),
            'time_dprt_ms': float(td), 'time_cs_ms': float(tc)
        })
    
    all_results[f'd={d}'] = d_result
    dt = time.perf_counter()-dt0
    best_r = min(r['full_ratio'] for r in d_result['runs'])
    spd = d_result['runs'][-1]['time_cs_ms']/max(d_result['runs'][-1]['time_dprt_ms'],0.001)
    print(f"  d={d:>3d} done {dt:>5.0f}s | full best={best_r:.4f} | 1loc best={min(r['k1_ratio'] for r in d_result['runs']):.4f} | speedup={spd:.0f}x", flush=True)

elapsed = time.time() - T0
print(f"\n═══ ALL {len(primes)} PRIMES 23–97 COMPLETE — {elapsed:.0f}s ═══", flush=True)

# Grand summary
print(f"\n{'d':>4s} | {'r=1 full':>9s} {'r=2 full':>9s} {'r=4 full':>9s} | {'r=1 1loc':>9s} {'r=2 1loc':>9s} {'r=4 1loc':>9s} | {'spd':>6s} | best")
print('-'*90)
for d in primes:
    r = all_results[f'd={d}']['runs']
    f1,f2,f4 = r[0]['full_ratio'], r[1]['full_ratio'], r[2]['full_ratio']
    k1,k2,k4 = r[0]['k1_ratio'], r[1]['k1_ratio'], r[2]['k1_ratio']
    spd = r[-1]['time_cs_ms']/max(r[-1]['time_dprt_ms'],0.001)
    best = min(f1,f2,f4)
    m = " ★" if best < 0.95 else ""
    print(f"{d:>4d} | {f1:>9.4f} {f2:>9.4f} {f4:>9.4f} | {k1:>9.4f} {k2:>9.4f} {k4:>9.4f} | {spd:>5.0f}x | {best:.4f}{m}")

with open('primes_23_97_results.json', 'w') as f:
    json.dump(all_results, f, indent=2)
print(f"\n→ primes_23_97_results.json saved", flush=True)
