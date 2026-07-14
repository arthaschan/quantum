#!/usr/bin/env python3
"""
RadonShadow: ALL primes ≤ 1000 — 1-local observable scan (FFT-optimized)
======================================================================
Strategy:
- FFT for probabilities p_m(k) = |fft(ψ·ω^{-m·j(j-1)/2})|²/d  → O(d² log d)
- 1-local observables: only need 2 components per basis vector
- d ≤ 199: full-state + 1-local (5 states)
- d > 199: 1-local only (3 states → 2 states for d > 500)
"""
import numpy as np, math, time, json, sys

def primes_upto(n):
    """Sieve for primes ≤ n."""
    sieve = bytearray(b'\x01')*(n+1)
    sieve[:2] = b'\x00\x00'
    for i in range(2, int(n**0.5)+1):
        if sieve[i]:
            sieve[i*i:n+1:i] = b'\x00'*((n-i*i)//i+1)
    return [i for i in range(2, n+1) if sieve[i]]

def random_pure_state(d, seed):
    rng = np.random.RandomState(seed)
    psi = rng.randn(d) + 1j * rng.randn(d)
    return psi / np.linalg.norm(psi)

def probs_mub_fft(psi, m, d, omega, precomp):
    """Fast p_m(k) via FFT: p_m(k) = |fft(ψ_j·ω^{-m·j(j-1)/2})[k]|²/d"""
    if m == 0:
        return np.abs(psi)**2
    v = psi * precomp[m]  # precomp[m][j] = ω^{-m·j(j-1)/2}
    f = np.fft.fft(v)
    return np.abs(f)**2 / d

def sample_outcomes(probs, n_shots, rng):
    """Sample n_shots outcomes from categorical distribution."""
    cum = np.cumsum(probs)
    return np.searchsorted(cum, rng.rand(n_shots))

def estimate_1local(d, n_total, mubs, psi, shots_data, rng, delta=15):
    """
    Deterministic or random estimation of 1-local observables.
    shots_data: list of (m, k) pairs defining which MUB and which outcome.
    """
    obs_list = []
    count = 0
    for a in range(d):
        if count >= delta: break
        for b in range(a+1, d):
            if count >= delta: break
            obs_list.append((a, b))  # 1-local: |a⟩⟨b| + |b⟩⟨a|
            count += 1
    
    # Pre-compute basis vector components for all needed (m,k,a) and (m,k,b)
    omega = np.exp(2j * np.pi / d)
    estimates = np.zeros(len(obs_list))
    
    for m_sampled, k_sampled in shots_data:
        # Compute |b_m^(k)⟩ components for all needed indices
        comps = {}
        needed = set()
        for a,b in obs_list:
            needed.add(a); needed.add(b)
        
        for idx in needed:
            if m_sampled == 0:
                comps[idx] = 1.0 if idx == k_sampled else 0.0
            else:
                phase = m_sampled * idx * (idx-1) // 2 + k_sampled * idx
                comps[idx] = omega**phase / np.sqrt(d)
        
        for oi, (a,b) in enumerate(obs_list):
            val = 2 * np.real(np.conj(comps[a]) * comps[b])  # ⟨ψ|O|ψ⟩
            estimates[oi] += (d+1) * val  # HKP: (d+1)⟨ψ|O|ψ⟩ - tr(O), tr(O)=0
    
    return estimates / len(shots_data)

def run_d(d, obs_delta=15):
    """Run experiment for single prime d."""
    omega = np.exp(2j * np.pi / d)
    
    # Precompute ω^{-m·j(j-1)/2} for m=1..d, j=0..d-1
    precomp = {m: np.array([omega**(-m*j*(j-1)//2) for j in range(d)], dtype=complex) 
               for m in range(1, d+1)}
    
    n_states = 5 if d <= 199 else (3 if d <= 500 else 2)
    r_list = [1, 2, 4]
    n_mubs = d + 1
    
    result = {'d': d, 'n_mubs': n_mubs, 'n_states': n_states, 'runs': []}
    
    dt0 = time.perf_counter()
    
    for r in r_list:
        n_per = r
        n_total = n_mubs * r
        
        ratios_1loc = []
        times_det = []
        times_rand = []
        
        # Full-state (only for d ≤ 199)
        ratios_full = [] if d <= 199 else None
        
        for si in range(n_states):
            # Generate random pure state
            psi = random_pure_state(d, d*10000 + r*100 + si*10)
            
            # --- DETERMINISTIC ---
            t0 = time.perf_counter()
            det_snaps = []
            for m in range(n_mubs):
                probs = probs_mub_fft(psi, m, d, omega, precomp)
                probs = np.clip(probs, 0, None); probs /= probs.sum()
                rng_d = np.random.RandomState(d*20000 + r*100 + si*10 + m)
                outcomes = sample_outcomes(probs, n_per, rng_d)
                for k in outcomes:
                    det_snaps.append((m, int(k)))
            time_det = (time.perf_counter()-t0)*1000
            
            # --- RANDOM (Classical Shadow) ---
            t0 = time.perf_counter()
            rand_snaps = []
            # Precompute all MUB probs for random sampling
            all_probs = []
            for m in range(n_mubs):
                probs = probs_mub_fft(psi, m, d, omega, precomp)
                probs = np.clip(probs, 0, None); probs /= probs.sum()
                all_probs.append(probs)
            
            rng_r = np.random.RandomState(d*30000 + r*100 + si*10)
            for _ in range(n_total):
                m = rng_r.randint(n_mubs)
                k = int(sample_outcomes(all_probs[m], 1, rng_r)[0])
                rand_snaps.append((m, k))
            time_rand = (time.perf_counter()-t0)*1000
            
            # --- 1-local estimation ---
            obs_list = []
            count = 0
            for a in range(d):
                if count >= obs_delta: break
                for b in range(a+1, d):
                    if count >= obs_delta: break
                    obs_list.append((a, b))
                    count += 1
            
            # True values
            truth = []
            for a, b in obs_list:
                truth.append(2 * np.real(np.conj(psi[a]) * psi[b]))
            
            # Deterministic estimates
            omega_arr = omega
            n_det = len(det_snaps)
            est_det = np.zeros(len(obs_list))
            for m_s, k_s in det_snaps:
                comps = {}
                needed = set()
                for a,b in obs_list:
                    needed.add(a); needed.add(b)
                for idx in needed:
                    if m_s == 0:
                        comps[idx] = 1.0 if idx == k_s else 0.0
                    else:
                        phase = m_s * idx * (idx-1) // 2 + k_s * idx
                        comps[idx] = omega_arr**phase / np.sqrt(d)
                for oi, (a,b) in enumerate(obs_list):
                    est_det[oi] += 2 * np.real(np.conj(comps[a]) * comps[b])
            est_det *= (d+1) / n_det
            
            # Random estimates
            est_rand = np.zeros(len(obs_list))
            for m_s, k_s in rand_snaps:
                comps = {}
                needed = set()
                for a,b in obs_list:
                    needed.add(a); needed.add(b)
                for idx in needed:
                    if m_s == 0:
                        comps[idx] = 1.0 if idx == k_s else 0.0
                    else:
                        phase = m_s * idx * (idx-1) // 2 + k_s * idx
                        comps[idx] = omega_arr**phase / np.sqrt(d)
                for oi, (a,b) in enumerate(obs_list):
                    est_rand[oi] += 2 * np.real(np.conj(comps[a]) * comps[b])
            est_rand *= (d+1) / n_total
            
            mae_det = np.mean(np.abs(est_det - np.array(truth)))
            mae_rand = np.mean(np.abs(est_rand - np.array(truth)))
            ratios_1loc.append(mae_det / mae_rand if mae_rand > 0 else 1.0)
            
            times_det.append(time_det)
            times_rand.append(time_rand)
            
            # Full-state reconstruction (d ≤ 199 only)
            if ratios_full is not None:
                # Deterministic full-state via MUB linear inversion
                rho_est_d = np.zeros((d,d), dtype=complex)
                for m in range(n_mubs):
                    probs = probs_mub_fft(psi, m, d, omega, precomp)
                    probs = np.clip(probs, 0, None); probs /= probs.sum()
                    rng_d2 = np.random.RandomState(d*40000 + r*100 + si*10 + m)
                    outcomes = sample_outcomes(probs, n_per, rng_d2)
                    counts = np.bincount(outcomes, minlength=d) / n_per
                    # Reconstruct using basis vectors
                    for k in range(d):
                        if counts[k] > 0:
                            if m == 0:
                                v = np.zeros(d, dtype=complex); v[k] = 1.0
                            else:
                                v = np.array([omega**(m*j*(j-1)//2 + k*j) for j in range(d)], dtype=complex) / np.sqrt(d)
                            rho_est_d += counts[k] * np.outer(v, v.conj())
                rho_est_d -= np.eye(d, dtype=complex)
                
                # CS full-state
                rho_est_c = np.zeros((d,d), dtype=complex)
                rng_c2 = np.random.RandomState(d*50000 + r*100 + si*10)
                for _ in range(n_total):
                    m = rng_c2.randint(n_mubs)
                    probs = all_probs[m]
                    k = int(sample_outcomes(probs, 1, rng_c2)[0])
                    if m == 0:
                        v = np.zeros(d, dtype=complex); v[k] = 1.0
                    else:
                        v = np.array([omega**(m*j*(j-1)//2 + k*j) for j in range(d)], dtype=complex) / np.sqrt(d)
                    rho_est_c += np.outer(v, v.conj())
                rho_est_c = rho_est_c * (d+1) / n_total - np.eye(d, dtype=complex)
                
                rho_true = np.outer(psi, psi.conj())
                err_d = np.linalg.norm(rho_true - rho_est_d, 'fro')
                err_c = np.linalg.norm(rho_true - rho_est_c, 'fro')
                ratios_full.append(err_d / err_c if err_c > 0 else 1.0)
        
        r1l = float(np.mean(ratios_1loc))
        td = float(np.mean(times_det))
        tr = float(np.mean(times_rand))
        
        run_data = {
            'r': r, 'n_total': n_total,
            'ratio_1loc': r1l,
            'time_det_ms': td, 'time_rand_ms': tr
        }
        if ratios_full is not None:
            rf = float(np.mean(ratios_full))
            run_data['ratio_full'] = rf
            run_data['full_best_ratio'] = rf
        else:
            rf = None
        
        result['runs'].append(run_data)
    
    dt = time.perf_counter()-dt0
    r1l_best = min(r['ratio_1loc'] for r in result['runs'])
    if ratios_full is not None:
        rf_best = min(r.get('ratio_full', 1) for r in result['runs'])
        print(f"  d={d:>4d} {dt:>5.1f}s | 1loc={r1l_best:.4f} | full={rf_best:.4f} | states={n_states}", flush=True)
    else:
        print(f"  d={d:>4d} {dt:>5.1f}s | 1loc={r1l_best:.4f} | states={n_states}", flush=True)
    
    return result

if __name__ == '__main__':
    all_primes = primes_upto(1000)
    print(f"Testing {len(all_primes)} primes: {all_primes[:5]}...{all_primes[-3:]}\n", flush=True)
    
    T0 = time.time()
    all_results = {}
    
    for d in all_primes:
        all_results[f'd={d}'] = run_d(d, obs_delta=15)
    
    elapsed = time.time() - T0
    print(f"\n═══ ALL {len(all_primes)} PRIMES ≤ 1000 COMPLETE — {elapsed:.0f}s ═══", flush=True)
    
    # Export results
    with open('primes_1000_results.json', 'w') as f:
        # Simplify for JSON size
        slim = {}
        for k,v in all_results.items():
            slim[k] = {
                'd': v['d'],
                'n_mubs': v['n_mubs'],
                'runs': [{'r': r['r'], 'ratio_1loc': r['ratio_1loc']} for r in v['runs']]
            }
            if v['runs'][0].get('ratio_full') is not None:
                from_runs = v['runs']
                slim[k]['ratio_full_best'] = min(r.get('ratio_full', 1) for r in from_runs)
        json.dump(slim, f, indent=2)
    
    # Summary statistics
    r1l_vals = [min(r['ratio_1loc'] for r in v['runs']) for v in all_results.values()]
    r1l_best = min(r1l_vals)
    r1l_worst = max(r1l_vals)
    r1l_median = sorted(r1l_vals)[len(r1l_vals)//2]
    count_under_09 = sum(1 for x in r1l_vals if x < 0.90)
    count_under_085 = sum(1 for x in r1l_vals if x < 0.85)
    count_under_08 = sum(1 for x in r1l_vals if x < 0.80)
    
    print(f"\n{'='*60}")
    print(f"SUMMARY STATISTICS")
    print(f"{'='*60}")
    print(f"  Total primes:     {len(all_primes)}")
    print(f"  1loc best:        {r1l_best:.4f}  (DPRT advantage {(1-r1l_best)*100:.1f}%)")
    print(f"  1loc worst:       {r1l_worst:.4f}")
    print(f"  1loc median:      {r1l_median:.4f}")
    print(f"  ratio < 0.90:     {count_under_09}/{len(all_primes)} ({100*count_under_09/len(all_primes):.1f}%)")
    print(f"  ratio < 0.85:     {count_under_085}/{len(all_primes)} ({100*count_under_085/len(all_primes):.1f}%)")
    print(f"  ratio < 0.80:     {count_under_08}/{len(all_primes)} ({100*count_under_08/len(all_primes):.1f}%)")
    
    # Top-10 best dimensions
    print(f"\n── Top 10 Best Dimensions (1-local) ──")
    ranked = sorted([(min(r['ratio_1loc'] for r in v['runs']), v['d']) for v in all_results.values()])
    for ratio, d in ranked[:10]:
        adv = (1-ratio)*100
        print(f"  d={d:>4d}  ratio={ratio:.4f}  advantage={adv:.1f}%")
    
    # Pattern: ratio vs d range
    print(f"\n── Ratio by d-Range ──")
    ranges = [(0,50), (50,100), (100,200), (200,500), (500,1000)]
    for lo, hi in ranges:
        in_range = [min(r['ratio_1loc'] for r in v['runs']) 
                    for v in all_results.values() if lo <= v['d'] < hi]
        if in_range:
            m = np.mean(in_range)
            med = sorted(in_range)[len(in_range)//2]
            mn = min(in_range)
            mx = max(in_range)
            print(f"  d∈[{lo},{hi}): n={len(in_range):>3d}  mean={m:.4f}  median={med:.4f}  min={mn:.4f}  max={mx:.4f}")
    
    print(f"\n→ primes_1000_results.json saved ({elapsed:.0f}s)", flush=True)
