#!/usr/bin/env python3
"""
RadonShadow: Sampled primes 1000–10000 (1-local only, sparse)
- ~30 primes uniformly sampled from 168 total in 1000-10000
- 1 state, r=4, 1-local only (8 observables)
- Single-pass FFT, no full-state
"""
import numpy as np, math, time, json, sys

def primes_upto(n):
    sieve = bytearray(b'\x01')*(n+1)
    sieve[:2] = b'\x00\x00'
    for i in range(2, int(n**0.5)+1):
        if sieve[i]:
            sieve[i*i:n+1:i] = b'\x00'*((n-i*i)//i+1)
    return [i for i in range(2, n+1) if sieve[i]]

def run_d_fast(d, seed_offset=0):
    """
    Ultra-fast single-pass for 1 d:
    - 1 random pure state
    - r=4 per MUB
    - 8 1-local observables
    """
    omega = np.exp(2j * np.pi / d)
    n_mubs = d + 1
    n_per = 4
    n_total = n_mubs * n_per
    
    # Generate pure state
    rng_state = np.random.RandomState(d * 100000 + seed_offset)
    psi = rng_state.randn(d) + 1j * rng_state.randn(d)
    psi /= np.linalg.norm(psi)
    
    # Select 8 sparse 1-local observables: |a⟩⟨b| + |b⟩⟨a|
    obs_indices = []
    step = max(1, d // 8)
    for i in range(8):
        a = i * step
        b = (i * step + step // 2) % d
        if a == b: b = (b + 1) % d
        obs_indices.append((a, b))
    
    # Truth values
    truth = np.array([2 * np.real(np.conj(psi[a]) * psi[b]) for a,b in obs_indices])
    
    # --- Precompute all MUB probabilities (single pass) ---
    all_probs = []
    for m in range(n_mubs):
        if m == 0:
            probs = np.abs(psi)**2
        else:
            # Precompute ω^{-m·j(j-1)/2} · ψ_j
            phase = np.array([omega**(-m * j * (j-1) // 2) for j in range(d)], dtype=complex)
            v = psi * phase
            f = np.fft.fft(v)
            probs = np.abs(f)**2 / d
        probs = np.clip(probs, 0, None)
        probs /= probs.sum()
        all_probs.append(probs)
    
    # --- DETERMINISTIC: n_per shots per MUB ---
    det_m_k = []
    rng_d = np.random.RandomState(d * 200000 + seed_offset)
    for m in range(n_mubs):
        cum = np.cumsum(all_probs[m])
        for _ in range(n_per):
            k = int(np.searchsorted(cum, rng_d.rand()))
            det_m_k.append((m, k))
    
    # --- RANDOM: n_total random MUB selections ---
    rand_m_k = []
    rng_r = np.random.RandomState(d * 300000 + seed_offset)
    for _ in range(n_total):
        m = rng_r.randint(n_mubs)
        cum = np.cumsum(all_probs[m])
        k = int(np.searchsorted(cum, rng_r.rand()))
        rand_m_k.append((m, k))
    
    # --- Estimate 1-local observables ---
    sqrt_d = np.sqrt(d)
    
    def estimate(snap_list, n_snaps):
        est = np.zeros(8)
        for m_s, k_s in snap_list:
            # Compute only needed basis vector components
            for oi, (a_val, b_val) in enumerate(obs_indices):
                if m_s == 0:
                    ca = 1.0 if a_val == k_s else 0.0
                    cb = 1.0 if b_val == k_s else 0.0
                else:
                    phase_a = m_s * a_val * (a_val - 1) // 2 + k_s * a_val
                    phase_b = m_s * b_val * (b_val - 1) // 2 + k_s * b_val
                    ca = omega**phase_a / sqrt_d
                    cb = omega**phase_b / sqrt_d
                est[oi] += 2 * np.real(np.conj(ca) * cb)
        est *= (d + 1) / n_snaps
        return est
    
    est_det = estimate(det_m_k, n_total)
    est_rand = estimate(rand_m_k, n_total)
    
    mae_det = np.mean(np.abs(est_det - truth))
    mae_rand = np.mean(np.abs(est_rand - truth))
    ratio = mae_det / mae_rand if mae_rand > 0 else 1.0
    
    return ratio

if __name__ == '__main__':
    all_primes = primes_upto(10000)
    # Sample primes from 1000-10000 range
    # Total: ~1060 primes in 1000-10000, take every 32nd → ~33 samples
    primes_above_1000 = [p for p in all_primes if p >= 1000]
    step = max(1, len(primes_above_1000) // 30)
    sampled = primes_above_1000[::step]
    # Always include max
    if primes_above_1000[-1] not in sampled:
        sampled.append(primes_above_1000[-1])
    
    print(f"Total primes ≤10000: {len(all_primes)}")
    print(f"Primes ≥1000: {len(primes_above_1000)}")
    print(f"Sampling every {step}th → {len(sampled)} primes: {sampled[:5]}...{sampled[-3:]}\n", flush=True)
    
    T0 = time.time()
    results = []
    
    for i, d in enumerate(sampled):
        dt0 = time.perf_counter()
        ratio = run_d_fast(d, seed_offset=i)
        dt = time.perf_counter()-dt0
        adv = (1-ratio)*100
        marker = " ★" if ratio < 0.85 else (" ✦" if ratio < 0.90 else "")
        print(f"  [{i+1:>2d}/{len(sampled)}] d={d:>5d} {dt:>5.1f}s | ratio={ratio:.4f} | DPRT+{adv:>5.1f}%{marker}", flush=True)
        results.append({'d': d, 'ratio': float(ratio), 'time_s': float(dt)})
    
    elapsed = time.time()-T0
    ratios = [r['ratio'] for r in results]
    
    print(f"\n═══ SAMPLED PRIMES 1000–10000 COMPLETE — {elapsed:.0f}s ═══", flush=True)
    
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"  Sampled:    {len(results)} primes")
    print(f"  Mean:       {np.mean(ratios):.4f}")
    print(f"  Median:     {sorted(ratios)[len(ratios)//2]:.4f}")
    print(f"  Min:        {min(ratios):.4f}  (DPRT advantage {(1-min(ratios))*100:.1f}%)")
    print(f"  Max:        {max(ratios):.4f}")
    print(f"  < 0.90:     {sum(1 for r in ratios if r<0.90)}/{len(ratios)} ({100*sum(1 for r in ratios if r<0.90)/len(ratios):.1f}%)")
    print(f"  < 0.85:     {sum(1 for r in ratios if r<0.85)}/{len(ratios)} ({100*sum(1 for r in ratios if r<0.85)/len(ratios):.1f}%)")
    
    # Trend by sub-range
    print(f"\n── Trend by d-Range ──")
    ranges = [(1000,2000),(2000,4000),(4000,6000),(6000,8000),(8000,10000)]
    for lo, hi in ranges:
        in_r = [r['ratio'] for r in results if lo <= r['d'] < hi]
        if in_r:
            print(f"  [{lo},{hi}): n={len(in_r):>2d}  mean={np.mean(in_r):.4f}  median={np.median(in_r):.4f}  min={min(in_r):.4f}")
    
    # Top 5
    ranked = sorted(results, key=lambda x: x['ratio'])
    print(f"\n── Top 5 ──")
    for r in ranked[:5]:
        print(f"  d={r['d']:>5d}  ratio={r['ratio']:.4f}  advantage={(1-r['ratio'])*100:.1f}%")
    
    # Combine with ≤1000 results for grand summary
    # Load previous results
    with open('primes_1000_results.json') as f:
        old = json.load(f)
    
    old_ratios_1loc = [min(r['ratio_1loc'] for r in v['runs']) for v in old.values()]
    all_ratios = old_ratios_1loc + ratios
    all_med = sorted(all_ratios)[len(all_ratios)//2]
    all_under_090 = sum(1 for r in all_ratios if r < 0.90)
    all_under_085 = sum(1 for r in all_ratios if r < 0.85)
    
    print(f"\n═══ GRAND TOTAL (168 + {len(results)} samples) ═══")
    print(f"  Total points:  {len(all_ratios)}")
    print(f"  Overall median: {all_med:.4f}")
    print(f"  < 0.90:         {all_under_090}/{len(all_ratios)} ({100*all_under_090/len(all_ratios):.1f}%)")
    print(f"  < 0.85:         {all_under_085}/{len(all_ratios)} ({100*all_under_085/len(all_ratios):.1f}%)")
    
    with open('primes_10000_sampled.json', 'w') as f:
        json.dump({'sampled_count': len(results), 'results': results, 'summary': {
            'mean': float(np.mean(ratios)), 'median': float(np.median(ratios)),
            'min': float(min(ratios)), 'max': float(max(ratios)),
            'under_090': sum(1 for r in ratios if r<0.90),
            'under_085': sum(1 for r in ratios if r<0.85)
        }}, f, indent=2)
    print(f"\n→ primes_10000_sampled.json saved", flush=True)
