#!/usr/bin/env python3 -u
"""
E3 Extended: ALL composite dimensions ≤ 100
============================================
- Factor each composite d into prime factors
- Build MUBs via tensor product of prime MUBs
- Validate MUB property numerically
- Test 1-local DPRT deterministic vs random CS
- Categorize by factorization type
"""
import numpy as np, json, sys, math, time, itertools

# ═══ Helpers ═══
def primes_upto(n):
    s = bytearray(b'\x01')*(n+1); s[:2] = b'\x00\x00'
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i:n+1:i] = b'\x00'*((n-i*i)//i+1)
    return [i for i in range(2, n+1) if s[i]]

PRIMES = set(primes_upto(100))

# Cache for prime MUBs
MUB_CACHE = {}
def get_prime_mubs(p):
    if p not in MUB_CACHE:
        omega = np.exp(2j * np.pi / p)
        mubs = [list(np.eye(p, dtype=complex))]
        for m in range(1, p+1):
            basis = []
            for k in range(p):
                v = np.array([omega**(m*j*(j-1)//2 + k*j) for j in range(p)], dtype=complex)
                basis.append(v / np.sqrt(p))
            mubs.append(basis)
        MUB_CACHE[p] = mubs
    return MUB_CACHE[p]

def factor_into_primes(d):
    """Factor d into list of prime factors (with multiplicity)."""
    factors = []
    n = d
    for p in [2,3,5,7,11,13]:
        while n % p == 0:
            factors.append(p)
            n //= p
    if n > 1:
        # Remaining factor must be prime (d ≤ 100)
        factors.append(n)
    return factors

def build_composite_mubs(d):
    """Build MUBs for composite d via prime factor tensor product."""
    p_factors = factor_into_primes(d)
    if len(p_factors) == 1:
        return get_prime_mubs(p_factors[0])
    
    mubs_parts = [get_prime_mubs(p) for p in p_factors]
    # Start with first prime's MUBs
    result_mubs = mubs_parts[0]
    for next_mubs in mubs_parts[1:]:
        n_common = min(len(result_mubs), len(next_mubs))
        new_mubs = []
        for i in range(n_common):
            basis = []
            for v1 in result_mubs[i]:
                for v2 in next_mubs[i]:
                    basis.append(np.kron(v1, v2))
            new_mubs.append(basis)
        result_mubs = new_mubs
    return result_mubs

def validate_mubs(mubs, d, tol=1e-10):
    """Check orthonormality and MUB property."""
    n = len(mubs)
    for i in range(n):
        basis = mubs[i]
        if len(basis) != d:
            return False, f"basis {i} has {len(basis)} vectors, expected {d}"
        # Orthonormality
        for a in range(d):
            for b in range(d):
                ip = np.vdot(basis[a], basis[b])
                expected = 1.0 if a==b else 0.0
                if abs(ip - expected) > tol:
                    return False, f"basis {i}: ⟨{a}|{b}⟩ = {ip:.2e} ≠ {expected}"
        # MUB with other bases
        for j in range(i+1, n):
            for a in range(d):
                for b in range(d):
                    ip = abs(np.vdot(mubs[i][a], mubs[j][b]))
                    expected = 1.0 / np.sqrt(d)
                    if abs(ip - expected) > tol:
                        return False, f"bases {i},{j}: |⟨{a}|{b}⟩| = {ip:.6f} ≠ {expected:.6f}"
    return True, "OK"

def safe_probs(basis, psi):
    probs = np.array([np.abs(np.vdot(v, psi))**2 for v in basis])
    probs = np.maximum(probs.real, 0)
    s = probs.sum()
    return np.ones(len(basis))/len(basis) if s < 1e-15 else probs/s

def classify_factorization(p_factors):
    """Classify composite by structure."""
    unique = set(p_factors)
    counts = {p: p_factors.count(p) for p in unique}
    
    if len(unique) == 1:
        p = list(unique)[0]
        e = counts[p]
        return f"ppow"  # prime power p^e
    else:
        return f"mixed"  # mixed primes

# ═══ Main scan ═══
composites = [d for d in range(4, 101) if d not in PRIMES]
print(f"Testing {len(composites)} composite dimensions: {composites[:5]}...{composites[-3:]}", flush=True)

results = []
T0 = time.time()
skipped = []
valid_count = 0

for d in composites:
    p_factors = factor_into_primes(d)
    mubs = build_composite_mubs(d)
    n_mubs = len(mubs)
    
    # Validate MUBs
    valid, msg = validate_mubs(mubs, d)
    if not valid:
        skipped.append({'d': d, 'reason': msg, 'factors': p_factors, 'n_mubs': n_mubs})
        print(f"  d={d:>3d} = {'×'.join(str(p) for p in p_factors):>20s}  MUBs={n_mubs:>2d}  SKIP: {msg}", flush=True)
        continue
    
    valid_count += 1
    n_states, r_val = 4, 4
    n_per = r_val
    n_total = n_mubs * r_val
    ratios = []
    
    dt0 = time.perf_counter()
    for si in range(n_states):
        rng = np.random.RandomState(d * 10000 + si)
        psi = rng.randn(d) + 1j * rng.randn(d)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        
        # 8 1-local observables: nearest-neighbor
        obs = []
        step = max(1, d//8)
        for idx in range(8):
            a = idx * step % d
            b = (idx * step + step//2) % d
            if a != b: obs.append((a, b))
        truth = [2*np.real(rho[a,b]) for a,b in obs]
        
        # Deterministic
        det_est = np.zeros((d,d), dtype=complex)
        rng_d = np.random.RandomState(d*10000 + si + 777)
        for mi in range(n_mubs):
            basis = mubs[mi]
            probs = safe_probs(basis, psi)
            outcomes = rng_d.choice(d, n_per, p=probs)
            for k in outcomes:
                det_est += np.outer(basis[k], basis[k].conj())
        det_est = det_est * (d+1) / n_total - np.eye(d)
        
        # Random
        rand_est = np.zeros((d,d), dtype=complex)
        rng_r = np.random.RandomState(d*10000 + si + 999)
        for _ in range(n_total):
            mi = rng_r.randint(n_mubs)
            basis = mubs[mi]
            probs = safe_probs(basis, psi)
            k = rng_r.choice(d, p=probs)
            rand_est += np.outer(basis[k], basis[k].conj())
        rand_est = rand_est * (d+1) / n_total - np.eye(d)
        
        mae_d = np.mean([abs(2*np.real(det_est[a,b]) - t) for (a,b),t in zip(obs,truth)])
        mae_r = np.mean([abs(2*np.real(rand_est[a,b]) - t) for (a,b),t in zip(obs,truth)])
        if mae_r > 1e-10:
            ratios.append(mae_d / mae_r)
    
    dt = time.perf_counter()-dt0
    if ratios:
        mr = np.mean(ratios)
        mr_std = np.std(ratios)
        adv = (1-mr)*100
        cls = classify_factorization(p_factors)
        m = " ★DPRT" if mr < 0.95 else (" ≈" if abs(mr-1) < 0.05 else "  CS")
        fstr = '×'.join(str(p) for p in p_factors)
        print(f"  d={d:>3d} = {fstr:>25s}  MUBs={n_mubs:>2d}  ratio={mr:.4f}±{mr_std:.3f}  {adv:+5.1f}%  [{cls}]{m}", flush=True)
        results.append({
            'd': d, 'factors': p_factors, 'factor_str': fstr,
            'n_mubs': n_mubs, 'class': cls,
            'ratio_mean': float(mr), 'ratio_std': float(mr_std),
            'n_states': n_states
        })
    else:
        print(f"  d={d:>3d} SKIP (no valid data)", flush=True)

elapsed = time.time() - T0
print(f"\n═══ {valid_count} valid + {len(skipped)} skipped composites — {elapsed:.0f}s ═══", flush=True)

# ═══ Analysis ═══
print(f"\n{'='*60}")
print(f"SUMMARY BY CLASS")
print(f"{'='*60}")

# By class
by_class = {}
for r in results:
    cls = r['class']
    if cls not in by_class: by_class[cls] = []
    by_class[cls].append(r)

for cls in sorted(by_class.keys()):
    rs = [r['ratio_mean'] for r in by_class[cls]]
    n_m = sum(1 for r in by_class[cls] if r['ratio_mean'] < 0.95)
    print(f"\n  [{cls}] n={len(rs)}:", flush=True)
    print(f"    mean={np.mean(rs):.4f}  median={np.median(rs):.4f}", flush=True)
    print(f"    DPRT_better(<0.95): {n_m}/{len(rs)} ({100*n_m/len(rs):.0f}%)", flush=True)
    print(f"    min={min(rs):.4f}  max={max(rs):.4f}", flush=True)
    
    # Best 5
    ranked = sorted(by_class[cls], key=lambda x: x['ratio_mean'])[:5]
    print(f"    Best:", " | ".join(f"d={r['d']}={r['factor_str']}:{r['ratio_mean']:.3f}" for r in ranked))

# By n_mubs
print(f"\n{'='*60}")
print(f"BY NUMBER OF MUBS")
print(f"{'='*60}")
by_mubs = {}
for r in results:
    nm = r['n_mubs']
    if nm not in by_mubs: by_mubs[nm] = []
    by_mubs[nm].append(r)

for nm in sorted(by_mubs.keys()):
    rs = [r['ratio_mean'] for r in by_mubs[nm]]
    n_m = sum(1 for r in by_mubs[nm] if r['ratio_mean'] < 0.95)
    print(f"  MUBs={nm:>2d}: n={len(rs):>3d}  mean={np.mean(rs):.4f}  DPRT<0.95={n_m}/{len(rs)}")

# Overall
all_ratios = [r['ratio_mean'] for r in results]
print(f"\n{'='*60}")
print(f"OVERALL ({len(results)} composites)")
print(f"{'='*60}")
print(f"  mean={np.mean(all_ratios):.4f}  median={np.median(all_ratios):.4f}")
print(f"  <0.95: {sum(1 for r in all_ratios if r<0.95)}/{len(all_ratios)} ({100*sum(1 for r in all_ratios if r<0.95)/len(all_ratios):.0f}%)")
print(f"  <0.90: {sum(1 for r in all_ratios if r<0.90)}/{len(all_ratios)} ({100*sum(1 for r in all_ratios if r<0.90)/len(all_ratios):.0f}%)")
print(f"  min={min(all_ratios):.4f}  max={max(all_ratios):.4f}")

# Top 10
ranked = sorted(results, key=lambda x: x['ratio_mean'])[:10]
print(f"\n── Top 10 Composites ──")
for r in ranked:
    print(f"  d={r['d']:>3d} = {r['factor_str']:>25s}  MUBs={r['n_mubs']:>2d}  ratio={r['ratio_mean']:.4f}  DPRT+{(1-r['ratio_mean'])*100:.1f}%")

# Save
with open('composite_100_results.json', 'w') as f:
    json.dump({'results': results, 'skipped': skipped, 'summary': {
        'total_testable': len(composites),
        'valid': valid_count,
        'skipped': len(skipped),
        'mean': float(np.mean(all_ratios)),
        'median': float(np.median(all_ratios)),
        'min': float(min(all_ratios)),
        'dprt_better_pct': float(100*sum(1 for r in all_ratios if r<0.95)/len(all_ratios))
    }}, f, indent=2)
print(f"\n→ composite_100_results.json saved", flush=True)
