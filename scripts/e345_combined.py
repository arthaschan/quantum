#!/usr/bin/env python3 -u
"""E3+E4+E5: Composite, 2-local, Noise — all-in-one, robust"""
import numpy as np, json, sys, math

# ═══ HELPERS ═══
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

def safe_probs(basis, psi):
    """Safe probability computation for a basis."""
    probs = np.array([np.abs(np.vdot(v, psi))**2 for v in basis])
    probs = np.maximum(probs.real, 0)
    s = probs.sum()
    if s < 1e-15:
        probs = np.ones(len(probs)) / len(probs)
    else:
        probs /= s
    return probs

# ═══════════════════════════════════════
# E3: COMPOSITE DIMENSIONS
# ═══════════════════════════════════════
print("="*50)
print("E3: COMPOSITE DIMENSIONS")
print("="*50, flush=True)

composite_dims = [4, 6, 8, 9, 10, 12, 14, 15, 16]
e3_results = []

for d in composite_dims:
    factors = factorize(d)
    factor_str = "×".join(f"{p}^{e}" if e>1 else str(p) for p,e in factors)
    
    # Build MUBs: for prime power, use direct Wootters; for composite, tensor product
    mubs_parts = [build_mubs_prime(p**e) for p,e in factors]
    mubs = mubs_parts[0]
    for nm in mubs_parts[1:]:
        mubs = [[np.kron(v1, v2) for v1, v2 in zip(b1, b2)] for b1, b2 in zip(mubs, nm)]
    
    n_mubs = min(len(m) for m in mubs_parts)
    mubs = mubs[:n_mubs]
    
    n_states, r_val = 5, 4
    n_per = r_val
    n_total = n_mubs * r_val
    ratios = []
    
    for si in range(n_states):
        rng = np.random.RandomState(d * 1000 + si * 10)
        psi = rng.randn(d) + 1j * rng.randn(d)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        
        obs = []
        for a in range(min(d, 4)):
            b = (a+1) % d
            if a != b:
                obs.append((a, b))
        truth = [2*np.real(rho[a,b]) for a,b in obs]
        n_obs = len(obs)
        
        # Deterministic
        det_est = np.zeros((d,d), dtype=complex)
        rng_d = np.random.RandomState(d*1000 + si*10 + 777)
        for mi in range(n_mubs):
            basis = mubs[mi]
            probs = safe_probs(basis, psi)
            outcomes = rng_d.choice(len(basis), n_per, p=probs)
            for k in outcomes:
                vk = basis[k]
                det_est += np.outer(vk, vk.conj())
        det_est = det_est * (d+1) / n_total - np.eye(d)
        
        # Random
        rand_est = np.zeros((d,d), dtype=complex)
        rng_r = np.random.RandomState(d*1000 + si*10 + 999)
        for _ in range(n_total):
            mi = rng_r.randint(n_mubs)
            basis = mubs[mi]
            probs = safe_probs(basis, psi)
            k = rng_r.choice(len(basis), p=probs)
            vk = basis[k]
            rand_est += np.outer(vk, vk.conj())
        rand_est = rand_est * (d+1) / n_total - np.eye(d)
        
        mae_d = np.mean([abs(2*np.real(det_est[a,b]) - t) for (a,b),t in zip(obs,truth)])
        mae_r = np.mean([abs(2*np.real(rand_est[a,b]) - t) for (a,b),t in zip(obs,truth)])
        if mae_r > 0:
            ratios.append(mae_d / mae_r)
    
    mr = np.mean(ratios)
    adv = (1-mr)*100
    m = "★" if mr<0.95 else ("≈" if abs(mr-1)<0.05 else "  ")
    print(f"  d={d:>3d} = {factor_str:>10s}  MUBs={n_mubs:>2d}  ratio={mr:.4f}  {adv:>+5.1f}% {m}", flush=True)
    e3_results.append({'d': d, 'ratio': float(mr), 'n_mubs': n_mubs, 'factors': factor_str})

# ═══════════════════════════════════════
# E4: 2-LOCAL (d≤100 primes)
# ═══════════════════════════════════════
print("\n" + "="*50)
print("E4: 2-LOCAL OBSERVABLES")
print("="*50, flush=True)

e4_results = []
primes_small = [p for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] if p not in [2]]

for d in primes_small:
    mubs = build_mubs_prime(d)
    n_mubs = len(mubs)
    n_states, r_val = 3, 4
    n_per = r_val
    n_total = n_mubs * r_val
    
    ratios_2loc = []
    
    for si in range(n_states):
        rng = np.random.RandomState(d * 20000 + si)
        psi = rng.randn(d) + 1j * rng.randn(d)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        
        # 2-local: use 8 observables that grab non-adjacent matrix elements
        obs_2loc = []
        n_o = min(8, d//2)
        step = max(1, d//n_o)
        for idx in range(n_o):
            a, b = idx*2, (idx*2+1) % d
            c, d2 = (idx*2+2)%d, (idx*2+3)%d
            obs_2loc.append((a, b, c, d2))
        
        truth = []
        for a,b,c,d2 in obs_2loc:
            # 2-local ≈ use 4-element tensor structure in d-dim space
            t = np.real(rho[a,b] * rho[c,d2])  # cross-correlation
            truth.append(t)
        
        # Deterministic
        det_snaps = []
        for mi in range(n_mubs):
            basis = mubs[mi]
            probs = safe_probs(basis, psi)
            rng_d = np.random.RandomState(d*20000 + si + mi + 1000)
            outcomes = rng_d.choice(d, n_per, p=probs)
            for k in outcomes:
                det_snaps.append(basis[k])
        
        # Random
        rand_snaps = []
        rng_r = np.random.RandomState(d*20000 + si + 2000)
        for _ in range(n_total):
            mi = rng_r.randint(n_mubs)
            basis = mubs[mi]
            probs = safe_probs(basis, psi)
            k = rng_r.choice(d, p=probs)
            rand_snaps.append(basis[k])
        
        # 2-local estimation from shadows
        # Shadow predicts: tr(O·3(d+1)|v⟩⟨v|-I)/3(d+1) ... this is wrong for 2-local
        # For shadow: E[tr(O·(d+1)|v⟩⟨v|-I)] = tr(O·ρ), but for 2-local with k=2 basis measurement:
        # shadow norm blows up, so use direct state estimation
        
        # Use naive: estimate ρ from deterministic shadows, then compute 2-local
        det_rho = np.zeros((d,d), dtype=complex)
        for vk in det_snaps:
            det_rho += np.outer(vk, vk.conj())
        det_rho = det_rho * (d+1) / len(det_snaps) - np.eye(d)
        
        rand_rho = np.zeros((d,d), dtype=complex)
        for vk in rand_snaps:
            rand_rho += np.outer(vk, vk.conj())
        rand_rho = rand_rho * (d+1) / len(rand_snaps) - np.eye(d)
        
        mae_d = np.mean([abs(np.real(det_rho[a,b]*det_rho[c,d2]) - t) for (a,b,c,d2),t in zip(obs_2loc,truth)])
        mae_r = np.mean([abs(np.real(rand_rho[a,b]*rand_rho[c,d2]) - t) for (a,b,c,d2),t in zip(obs_2loc,truth)])
        if mae_r > 1e-10:
            ratios_2loc.append(mae_d / mae_r)
    
    if ratios_2loc:
        mr = np.mean(ratios_2loc)
        adv = (1-mr)*100
        m = "★" if mr<0.95 else ("≈" if abs(mr-1)<0.05 else "  ")
        print(f"  d={d:>3d}  ratio(2loc)={mr:.4f}  {adv:>+5.1f}% {m}", flush=True)
        e4_results.append({'d': d, 'ratio_2loc': float(mr)})
    else:
        print(f"  d={d:>3d}  SKIPPED (no valid data)", flush=True)

# ═══════════════════════════════════════
# E5: DEPOLARIZING NOISE
# ═══════════════════════════════════════
print("\n" + "="*50)
print("E5: DEPOLARIZING NOISE")
print("="*50, flush=True)

noise_levels = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3]
test_dims = [3, 5, 7, 11, 13, 17, 19, 23, 31, 37, 47, 61, 97]
e5_results = []

for noise in noise_levels:
    ratios_noise = []
    for d in test_dims:
        mubs = build_mubs_prime(d)
        n_mubs = len(mubs)
        n_states, r_val = 5, 4
        n_per = r_val
        n_total = n_mubs * r_val
        ratios_d = []
        
        for si in range(n_states):
            rng = np.random.RandomState(d * 10000 + int(noise*1000) + si)
            psi = rng.randn(d) + 1j * rng.randn(d)
            psi /= np.linalg.norm(psi)
            rho_pure = np.outer(psi, psi.conj())
            rho = (1-noise) * rho_pure + noise * np.eye(d, dtype=complex)/d
            
            obs = []
            for a in range(d):
                if len(obs) >= 8: break
                b = (a+1)%d
                if a != b: obs.append((a,b))
            truth = [2*np.real(rho[a,b]) for a,b in obs]
            
            # Deterministic
            det_est = np.zeros((d,d), dtype=complex)
            for mi in range(n_mubs):
                basis = mubs[mi]
                probs = safe_probs(basis, psi)
                # Measurement from noisy state:
                noisy_probs = np.array([np.real(np.trace(rho @ np.outer(v, v.conj()))) for v in basis])
                noisy_probs = np.maximum(noisy_probs, 0)
                noisy_probs /= noisy_probs.sum()
                rng_d = np.random.RandomState(d*10000 + int(noise*1000) + si + mi)
                outcomes = rng_d.choice(d, n_per, p=noisy_probs)
                for k in outcomes:
                    det_est += np.outer(basis[k], basis[k].conj())
            det_est = det_est * (d+1) / n_total - np.eye(d)
            
            # Random
            rand_est = np.zeros((d,d), dtype=complex)
            rng_r = np.random.RandomState(d*20000 + int(noise*1000) + si)
            for _ in range(n_total):
                mi = rng_r.randint(n_mubs)
                basis = mubs[mi]
                noisy_probs = np.array([np.real(np.trace(rho @ np.outer(v, v.conj()))) for v in basis])
                noisy_probs = np.maximum(noisy_probs, 0)
                noisy_probs /= noisy_probs.sum()
                k = rng_r.choice(d, p=noisy_probs)
                rand_est += np.outer(basis[k], basis[k].conj())
            rand_est = rand_est * (d+1) / n_total - np.eye(d)
            
            mae_d = np.mean([abs(2*np.real(det_est[a,b]) - t) for (a,b),t in zip(obs,truth)])
            mae_r = np.mean([abs(2*np.real(rand_est[a,b]) - t) for (a,b),t in zip(obs,truth)])
            if mae_r > 1e-10:
                ratios_d.append(mae_d / mae_r)
        
        if ratios_d:
            ratios_noise.extend(ratios_d)
    
    if ratios_noise:
        mr = np.mean(ratios_noise)
        med = np.median(ratios_noise)
        adv = (1-mr)*100
        m = " ★DPRT" if mr<0.95 else (" ≈" if abs(mr-1)<0.05 else "  CS")
        print(f"  λ={noise:5.2f}  mean_ratio={mr:.4f}  median={med:.4f}  {adv:>+5.1f}%{m}", flush=True)
        e5_results.append({'noise': noise, 'mean_ratio': float(mr), 'median_ratio': float(med)})
    else:
        print(f"  λ={noise:5.2f}  SKIPPED", flush=True)

# ═══ SAVE ═══
all_out = {'e3_composite': e3_results, 'e4_2local': e4_results, 'e5_noise': e5_results}
with open('all5_part2_results.json', 'w') as f:
    json.dump(all_out, f, indent=2)

# ═══ FINAL SUMMARY ═══
print("\n" + "="*50)
print("ALL E3-E5 COMPLETE")
print("="*50)

# E3 summary
print("\n── E3 Composite Summary ──")
better = [r for r in e3_results if r['ratio'] < 0.95]
print(f"  DPRT better (<0.95): {len(better)}/{len(e3_results)}")
for r in better:
    print(f"    d={r['d']} = {r['factors']:>10s} ratio={r['ratio']:.4f}")

# E4 summary
print("\n── E4 2-local Summary ──")
r2l = [r['ratio_2loc'] for r in e4_results]
print(f"  mean={np.mean(r2l):.4f}  median={np.median(r2l):.4f}")
print(f"  <0.95: {sum(1 for r in r2l if r<0.95)}/{len(r2l)}")
print(f"  <0.90: {sum(1 for r in r2l if r<0.90)}/{len(r2l)}")

# E5 summary
print("\n── E5 Noise Summary ──")
for r in e5_results:
    adv = (1-r['mean_ratio'])*100
    print(f"  λ={r['noise']:.2f}: mean={r['mean_ratio']:.4f}  DPRT+{adv:+5.1f}%")

print("\n→ all5_part2_results.json saved", flush=True)
