#!/usr/bin/env python3 -u
"""E2: Multiplicative subgroup analysis — standalone with flush"""
import numpy as np, json, sys, math

with open('primes_1000_results.json') as f:
    data = json.load(f)

subgroup_data = []
for k, v in data.items():
    d = v['d']
    best = min(r['ratio_1loc'] for r in v['runs'])
    phi = d - 1
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
    n_unique_pf = len(set(prime_factors))
    
    subgroup_data.append({
        'd': d, 'ratio': best, 'phi': phi,
        'largest_pf': largest_pf, 'n_unique_pf': n_unique_pf
    })

print("E2: SUBGROUP ANALYSIS", flush=True)

# By largest PF
by_lpf = {}
for sd in subgroup_data:
    lpf = sd['largest_pf']
    if lpf not in by_lpf: by_lpf[lpf] = []
    by_lpf[lpf].append(sd['ratio'])

print("\nLargest PF groups:", flush=True)
for lpf in sorted(by_lpf.keys())[:12]:
    rs = by_lpf[lpf]
    print(f"  LPF={lpf:>4d}: n={len(rs):>3d} mean={np.mean(rs):.4f} median={np.median(rs):.4f} min={min(rs):.4f}", flush=True)

# By unique PF count
by_uniq = {}
for sd in subgroup_data:
    upf = sd['n_unique_pf']
    if upf not in by_uniq: by_uniq[upf] = []
    by_uniq[upf].append(sd['ratio'])

print("\nUnique prime factors vs ratio:", flush=True)
for upf in sorted(by_uniq.keys()):
    rs = by_uniq[upf]
    print(f"  unique_pf={upf}: n={len(rs)} mean={np.mean(rs):.4f} median={np.median(rs):.4f}", flush=True)

# Spearman
x = [sd['n_unique_pf'] for sd in subgroup_data]
y = [sd['ratio'] for sd in subgroup_data]
n = len(x)
rx = {v: i+1 for i, v in enumerate(sorted(set(x)))}
ry = {v: i+1 for i, v in enumerate(sorted(set(y)))}
d2 = sum((rx[xi]-ry[yi])**2 for xi,yi in zip(x,y))
rho = 1 - 6*d2/(n*(n**2-1))
print(f"\nSpearman ρ(unique_pf, ratio) = {rho:.4f}", flush=True)

# Smoothness: ratio vs #small primes (≤11) in factorization
by_smooth = {}
for sd in subgroup_data:
    n = sd['phi']
    small_count = 0
    for p in [2,3,5,7,11]:
        while n % p == 0:
            small_count += 1
            n //= p
    if small_count not in by_smooth: by_smooth[small_count] = []
    by_smooth[small_count].append(sd['ratio'])

print("\nSmoothness (#small prime factors ≤11):", flush=True)
for sc in sorted(by_smooth.keys()):
    rs = by_smooth[sc]
    print(f"  smoothness={sc}: n={len(rs)} mean={np.mean(rs):.4f} median={np.median(rs):.4f}", flush=True)

print("\nE2 COMPLETE", flush=True)
