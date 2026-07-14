#!/usr/bin/env python3 -u
"""E3: Composite dimensions d=4,6,8,9,10,12,14,15,16"""
import numpy as np, json, sys

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

def factorize(d):
    n = d; factors = []; p = 2
    while p*p <= n:
        e = 0
        while n % p == 0: n //= p; e += 1
        if e > 0: factors.append((p, e))
        p += 1
    if n > 1: factors.append((n, 1))
    return factors

print("E3: COMPOSITE DIMENSIONS", flush=True)
composite_dims = [4, 6, 8, 9, 10, 12, 14, 15, 16]
for d in composite_dims:
    factors = factorize(d)
    factor_str = "×".join(f"{p}^{e}" if e>1 else str(p) for p,e in factors)
    
    # Build MUBs via tensor product of prime power factors
    mubs_parts = [build_mubs_prime(p**e) for p,e in factors]
    mubs = mubs_parts[0]
    for nm in mubs_parts[1:]:
        mubs = [[np.kron(v1, v2) for v1, v2 in zip(b1, b2)] for b1, b2 in zip(mubs, nm)]
    n_mubs = min(len(m) for m in mubs_parts)  # min complete set
    mubs = mubs[:n_mubs]
    
    n_states, r_val = 5, 4
    n_per = r_val; n_total = n_mubs * r_val
    ratios_1loc = []
    
    for si in range(n_states):
        rng_s = np.random.RandomState(d * 1000 + si * 10)
        psi = rng_s.randn(d) + 1j * rng_s.randn(d)
        psi /= np.linalg.norm(psi)
        
        obs_list = []
        c = 0
        for a in range(d):
            if c >= 8: break
            for b in range(a+1, d):
                if c >= 8: break
                obs_list.append((a,b)); c += 1
        truth = [2*np.real(np.conj(psi[a])*psi[b]) for a,b in obs_list]
        
        # Deterministic
        det_rho = np.zeros((d,d), dtype=complex)
        for mi in range(n_mubs):
            basis = mubs[mi]
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
            basis = mubs[mi]
            probs = np.array([np.abs(np.vdot(v, psi))**2 for v in basis])
            probs = np.clip(probs, 0, None); probs /= probs.sum()
            k = rng_r.choice(d, p=probs)
            rand_rho += np.outer(basis[k], basis[k].conj())
        rand_rho = rand_rho * (d+1) / n_total - np.eye(d)
        
        mae_d = np.mean([abs(2*np.real(det_rho[a,b]) - t) for (a,b),t in zip(obs_list,truth)])
        mae_r = np.mean([abs(2*np.real(rand_rho[a,b]) - t) for (a,b),t in zip(obs_list,truth)])
        if mae_r > 0:
            ratios_1loc.append(mae_d / mae_r)
    
    mean_r = np.mean(ratios_1loc)
    adv = (1-mean_r)*100
    m = " ★DPRT" if mean_r<0.95 else (" ==" if abs(mean_r-1)<0.05 else " CS")
    print(f"  d={d:>3d} = {factor_str:>10s}  MUBs={n_mubs:>2d}  ratio={mean_r:.4f}  {adv:>+5.1f}% {m}", flush=True)

print("\nE3 COMPLETE", flush=True)
