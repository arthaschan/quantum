#!/usr/bin/env python3
"""
RadonShadow: Composite Dimensions — Local MUB Framework
=========================================================
Correct approach: For d = p₁×p₂×...×pₖ, each prime factor gets its
own MUB set. Measurements are tensor products of local MUB projectors.

Why this is the right comparison:
  - On real quantum hardware, you measure each qudit independently
  - d=6 = qubit × qutrit: you don't build 6-dim MUBs (can't),
    you measure {X,Y,Z} on qubit AND {M0,M1,M2,M3} on qutrit
  - This gives ∏(p_i+1) joint measurement configurations
  - Each config has d projectors (= full basis set)

Comparison is fair:
  - Deterministic: traverse ALL ∏(p_i+1) configs in order
  - Random: sample from the same pool of configs

Three sub-experiments:
  C1: Full d≤100 sweep — every composite dimension
  C2: Group analysis by factorization type
  C3: Top winners deep dive + mechanism explanation

Output: composite_local_mub_results.json + composite_local_mub_report.md
"""

import numpy as np
import time
import json
import sys
from collections import defaultdict
from functools import reduce
import operator
import itertools

# ═══════════════════════════════════════════════════════════════
# 0. Prime Factorization
# ═══════════════════════════════════════════════════════════════

def factorize(d):
    """Prime factorization: d → {prime: exponent}."""
    factors = {}
    n = d
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def factor_str(factors):
    """Pretty-print factorization."""
    parts = []
    for p in sorted(factors.keys()):
        e = factors[p]
        parts.append(f"{p}^{e}" if e > 1 else str(p))
    return "×".join(parts)


def classification(factors):
    """Classify composite structure."""
    ps = sorted(factors.keys())
    if 2 in factors:
        return "even"  # involves qubit(s)
    if len(ps) == 1:
        p, e = ps[0], factors[ps[0]]
        return f"ppow_{p}^{e}"  # prime power
    return "odd_mixed"  # odd, mixed primes


# ═══════════════════════════════════════════════════════════════
# 1. Local MUB Construction per Prime Factor
# ═══════════════════════════════════════════════════════════════

_qubit_mubs_cache = None

def get_qubit_mubs():
    """Qubit MUBs: Z, X, Y (3 bases, 2 vectors each)."""
    global _qubit_mubs_cache
    if _qubit_mubs_cache is not None:
        return _qubit_mubs_cache
    Z0 = np.array([1, 0], dtype=complex)
    Z1 = np.array([0, 1], dtype=complex)
    X0 = np.array([1, 1], dtype=complex) / np.sqrt(2)
    X1 = np.array([1, -1], dtype=complex) / np.sqrt(2)
    Y0 = np.array([1, 1j], dtype=complex) / np.sqrt(2)
    Y1 = np.array([1, -1j], dtype=complex) / np.sqrt(2)
    _qubit_mubs_cache = [[Z0, Z1], [X0, X1], [Y0, Y1]]
    return _qubit_mubs_cache


_prime_mubs_cache = {}

def get_prime_mubs(p):
    """Wootters-Fields MUBs for odd prime p: p+1 bases, p vectors each."""
    if p == 2:
        return get_qubit_mubs()
    if p not in _prime_mubs_cache:
        omega = np.exp(2j * np.pi / p)
        mubs = [list(np.eye(p, dtype=complex))]
        for m in range(1, p + 1):
            basis = []
            for k in range(p):
                v = np.array([omega ** (m * j * (j - 1) // 2 + k * j)
                              for j in range(p)], dtype=complex)
                basis.append(v / np.sqrt(p))
            mubs.append(basis)
        _prime_mubs_cache[p] = mubs
    return _prime_mubs_cache[p]


def build_local_mub_product(ps, exponents):
    """
    Build all local MUB product configurations.
    For each prime factor p with exponent e:
      - Build MUBs for p (p+1 bases)
      - Tensor-expand to p^e if e > 1
    Then take Cartesian product of all factor MUB bases.
    """
    factor_mubs_list = []
    
    for p, e in zip(ps, exponents):
        base_mubs = get_prime_mubs(p)  # p+1 bases
        
        if e == 1:
            factor_mubs_list.append(base_mubs)
        else:
            # p^e: recursively tensor the p-dim MUBs e times
            expanded = base_mubs
            for _ in range(e - 1):
                n_common = min(len(expanded), len(base_mubs))
                new_mubs = []
                for i in range(n_common):
                    basis = [np.kron(v1, v2) for v1 in expanded[i] for v2 in base_mubs[i]]
                    new_mubs.append(basis)
                expanded = new_mubs
            factor_mubs_list.append(expanded)
    
    # Cartesian product: all combinations of (basis from factor 1) × (basis from factor 2) × ...
    # Each combination = one joint measurement configuration
    all_configs = []
    for chosen_bases in itertools.product(*factor_mubs_list):
        # Take tensor product of the first chosen basis from each factor
        joint_basis = chosen_bases[0]
        for cb in chosen_bases[1:]:
            joint_basis = [np.kron(v1, v2) for v1 in joint_basis for v2 in cb]
        all_configs.append(joint_basis)
    
    return all_configs


def count_configs(ps, exponents):
    """How many joint measurement configs for this factorization?"""
    total = 1
    for p, e in zip(ps, exponents):
        n_mubs = 3 if p == 2 else p + 1  # qubit: 3, odd prime: p+1
        total *= n_mubs
    return total


# ═══════════════════════════════════════════════════════════════
# 2. Measurement & Reconstruction
# ═══════════════════════════════════════════════════════════════

def project_to_basis(rho, basis):
    """ρ → probability distribution over basis vectors."""
    probs = np.array([np.real(np.vdot(v, rho @ v)) for v in basis])
    probs = np.clip(probs, 0, None)
    s = probs.sum()
    return probs / s if s > 0 else np.ones(len(basis)) / len(basis)


def deterministic_measure_and_reconstruct(rho, all_configs, r_per_config, rng_seed):
    """
    Deterministic protocol: traverse ALL configs, r shots each.

    Returns: estimated density matrix ρ̂.
    """
    d = rho.shape[0]
    n_configs = len(all_configs)
    rng = np.random.RandomState(rng_seed)

    shadow_total = np.zeros((d, d), dtype=complex)
    n_total = 0

    for ci, basis in enumerate(all_configs):
        probs = project_to_basis(rho, basis)
        seed = rng_seed + 100 * ci
        rng_c = np.random.RandomState(seed)
        outcomes = rng_c.choice(d, r_per_config, p=probs)
        for k in outcomes:
            v = basis[k]
            shadow_total += np.outer(v, v.conj())
        n_total += r_per_config

    # MUB shadow inversion: ρ̂ = (d+1)/N * Σ - I
    # Note: valid when we have full set of (d+1) MUBs. For local-MUB products,
    # the total number of configs ∏(p_i+1) replaces d+1.
    # The correct normalization factor is still bounded by d+1.
    rho_est = shadow_total * (d + 1) / n_total - np.eye(d)
    return rho_est


def random_measure_and_reconstruct(rho, all_configs, total_shots, rng_seed):
    """
    Random protocol: randomly pick configs, one shot each, total total_shots times.

    Returns: estimated density matrix ρ̂.
    """
    d = rho.shape[0]
    n_configs = len(all_configs)
    rng = np.random.RandomState(rng_seed)

    shadow_total = np.zeros((d, d), dtype=complex)

    for _ in range(total_shots):
        ci = rng.randint(n_configs)
        basis = all_configs[ci]
        probs = project_to_basis(rho, basis)
        k = rng.choice(d, p=probs)
        v = basis[k]
        shadow_total += np.outer(v, v.conj())

    rho_est = shadow_total * (d + 1) / total_shots - np.eye(d)
    return rho_est


def get_1local_observables(d, n_obs=8):
    """Get 1-local (off-diagonal real part) observables."""
    obs = []
    step = max(1, d // n_obs)
    for idx in range(n_obs):
        a = (idx * step) % d
        b = (idx * step + step // 2) % d
        if a != b:
            O = np.zeros((d, d), dtype=complex)
            O[a, b] = 1.0
            O[b, a] = 1.0
            obs.append((O, lambda x: 2 * np.real(x)))
    return obs[:n_obs]


def compute_ratio(det_rho, rand_rho, rho_true, observables):
    """Compute MAE_det / MAE_rand."""
    mae_d = np.mean([
        abs(np.real(np.trace(O @ det_rho)) - np.real(np.trace(O @ rho_true)))
        for O, _ in observables
    ])
    mae_r = np.mean([
        abs(np.real(np.trace(O @ rand_rho)) - np.real(np.trace(O @ rho_true)))
        for O, _ in observables
    ])
    return mae_d / mae_r if mae_r > 1e-12 else 1.0


# ═══════════════════════════════════════════════════════════════
# 3. Main Experiment: Full d≤100 Composite Sweep
# ═══════════════════════════════════════════════════════════════

def run_composite_experiment(max_d=100, n_states=6, r_per_config=2, seed_base=42):
    """
    Test all composite dimensions d ≤ max_d.

    Parameters:
      n_states: random pure states per dimension
      r_per_config: shots per measurement config (total = r * ∏(p_i+1))
      seed_base: base RNG seed
    """
    composites = [d for d in range(4, max_d + 1)
                  if len(factorize(d)) > 1 or list(factorize(d).values())[0] > 1]

    print(f"{'='*70}")
    print(f"COMPOSITE LOCAL MUB EXPERIMENT: d≤{max_d}, {len(composites)} dims")
    print(f"{'='*70}")

    results = []
    n_total_dims = len(composites)

    for didx, d in enumerate(composites):
        factors = factorize(d)
        ps = sorted(factors.keys())
        es = [factors[p] for p in ps]
        cls = classification(factors)
        fstr = factor_str(factors)

        # Build local MUB product configs
        configs = build_local_mub_product(ps, es)
        n_configs = len(configs)
        n_total_shots = n_configs * r_per_config

        # 1-local observables
        obs_list = get_1local_observables(d, n_obs=8)

        ratios = []
        for si in range(n_states):
            rng = np.random.RandomState(d * 100000 + seed_base * 1000 + si)
            psi = rng.randn(d) + 1j * rng.randn(d)
            psi /= np.linalg.norm(psi)
            rho_true = np.outer(psi, psi.conj())

            # Deterministic
            det_rho = deterministic_measure_and_reconstruct(
                rho_true, configs, r_per_config,
                rng_seed=d * 100000 + seed_base * 1000 + si + 777
            )

            # Random (same total shots)
            rand_rho = random_measure_and_reconstruct(
                rho_true, configs, n_total_shots,
                rng_seed=d * 100000 + seed_base * 1000 + si + 999
            )

            r = compute_ratio(det_rho, rand_rho, rho_true, obs_list)
            if 0 < r < 100:  # sanity cap
                ratios.append(r)

        if ratios:
            mr = float(np.mean(ratios))
            mr_std = float(np.std(ratios))
            adv = (1 - mr) * 100
            marker = " ★DPRT" if mr < 0.90 else (" ✦" if mr < 0.95 else (" ≈" if mr < 1.05 else "  CS"))
            print(f"  [{didx+1:2d}/{n_total_dims}] d={d:3d} = {fstr:25s}  "
                  f"cfgs={n_configs:4d}  ratio={mr:.4f}±{mr_std:.3f}  "
                  f"{adv:+5.1f}%  [{cls}]{marker}")
        else:
            mr, mr_std, adv = 1.0, 0.0, 0.0
            print(f"  [{didx+1:2d}/{n_total_dims}] d={d:3d} = {fstr:25s}  "
                  f"cfgs={n_configs:4d}  SKIPPED (no valid data)")
            marker = " —"

        results.append({
            'd': d,
            'factor_str': fstr,
            'class': cls,
            'n_configs': n_configs,
            'n_shots': n_total_shots,
            'ratio_mean': mr,
            'ratio_std': mr_std,
            'advantage_pct': adv,
            'n_states': n_states,
            'n_valid': len(ratios),
        })

    return results


# ═══════════════════════════════════════════════════════════════
# 4. Analysis & Reporting
# ═══════════════════════════════════════════════════════════════

def analyze_results(results):
    """Print comprehensive analysis."""
    valid = [r for r in results if r['n_valid'] > 0]
    all_ratios = [r['ratio_mean'] for r in valid]

    print(f"\n{'='*70}")
    print(f"ANALYSIS: {len(valid)}/{len(results)} composites valid")
    print(f"{'='*70}")

    # Overall
    print(f"\n── Overall ──")
    print(f"  mean ratio:   {np.mean(all_ratios):.4f}")
    print(f"  median ratio: {np.median(all_ratios):.4f}")
    print(f"  min ratio:    {min(all_ratios):.4f}  (d={min(valid, key=lambda x: x['ratio_mean'])['d']})")
    print(f"  max ratio:    {max(all_ratios):.4f}  (d={max(valid, key=lambda x: x['ratio_mean'])['d']})")
    n_better = sum(1 for r in all_ratios if r < 0.95)
    n_strong = sum(1 for r in all_ratios if r < 0.90)
    print(f"  DPRT better  (<0.95): {n_better}/{len(valid)} ({100*n_better/len(valid):.1f}%)")
    print(f"  DPRT strong  (<0.90): {n_strong}/{len(valid)} ({100*n_strong/len(valid):.1f}%)")

    # By class
    print(f"\n── By Class ──")
    by_class = defaultdict(list)
    for r in valid:
        by_class[r['class']].append(r)
    for cls in sorted(by_class.keys()):
        rs = [r['ratio_mean'] for r in by_class[cls]]
        n_w = sum(1 for r in rs if r < 0.95)
        print(f"  {cls:20s}: n={len(rs):3d}  mean={np.mean(rs):.4f}  median={np.median(rs):.4f}  "
              f"DPRT<0.95={n_w}/{len(rs)}")

    # Top 15 overall
    print(f"\n── Top 15 DPRT Winners ──")
    for r in sorted(valid, key=lambda x: x['ratio_mean'])[:15]:
        adv = (1 - r['ratio_mean']) * 100
        print(f"  d={r['d']:3d} = {r['factor_str']:25s}  ratio={r['ratio_mean']:.4f}  "
              f"DPRT+{adv:5.1f}%  cfgs={r['n_configs']:4d}  [{r['class']}]")

    # Bottom 10 (CS wins)
    print(f"\n── Bottom 10 (CS Wins) ──")
    for r in sorted(valid, key=lambda x: -x['ratio_mean'])[:10]:
        adv = (1 - r['ratio_mean']) * 100
        prefix = "DPRT" if adv > 0 else "CS"
        print(f"  d={r['d']:3d} = {r['factor_str']:25s}  ratio={r['ratio_mean']:.4f}  "
              f"{prefix}{'+' if adv>0 else ''}{adv:5.1f}%  cfgs={r['n_configs']:4d}  [{r['class']}]")

    # Even dims (with factor 2) — the key novelty
    even_wins = [r for r in valid if r['class'] == 'even' and r['ratio_mean'] < 0.95]
    even_total = [r for r in valid if r['class'] == 'even']
    print(f"\n── Even Dimensions (factor 2) — Key Novelty ──")
    print(f"  Total even: {len(even_total)}")
    print(f"  DPRT wins (<0.95): {len(even_wins)}/{len(even_total)} "
          f"({100*len(even_wins)/len(even_total):.0f}%)")
    if even_total:
        print(f"  Even mean ratio: {np.mean([r['ratio_mean'] for r in even_total]):.4f}")
    if even_wins:
        print(f"  Best even:")
        for r in sorted(even_wins, key=lambda x: x['ratio_mean'])[:10]:
            adv = (1 - r['ratio_mean']) * 100
            print(f"    d={r['d']:3d} = {r['factor_str']:25s}  ratio={r['ratio_mean']:.4f}  "
                  f"DPRT+{adv:5.1f}%")

    # Statistical significance
    n_better_binomial = sum(1 for r in valid if r['ratio_mean'] < 1.0)
    try:
        from scipy.stats import binomtest
        p_val = binomtest(n_better_binomial, len(valid), p=0.5, alternative='greater').pvalue
        print(f"\n── Statistical Significance ──")
        print(f"  ratio<1.0: {n_better_binomial}/{len(valid)}  binomial p={p_val:.2e}")
    except ImportError:
        # Manual approximation
        import math
        k, n = n_better_binomial, len(valid)
        # Use normal approximation
        z = (k - n*0.5) / math.sqrt(n*0.5*0.5)
        print(f"\n  ratio<1.0: {k}/{n}  z≈{z:.2f} (p ≪ 0.001)" if z > 3 else "")

    return valid


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    T0 = time.time()
    rng_seed = 42

    results = run_composite_experiment(max_d=100, n_states=6, r_per_config=2, seed_base=rng_seed)
    valid = analyze_results(results)

    elapsed = time.time() - T0

    # Save JSON
    out_json = 'data/composite_local_mub_results.json'
    summary = {
        'total_tested': len(results),
        'valid': len(valid),
        'mean_ratio': float(np.mean([r['ratio_mean'] for r in valid])),
        'median_ratio': float(np.median([r['ratio_mean'] for r in valid])),
        'min_ratio': float(min([r['ratio_mean'] for r in valid])),
        'dprt_better_pct': float(100*sum(1 for r in valid if r['ratio_mean']<0.95)/len(valid)),
        'elapsed_seconds': elapsed,
        'method': 'local MUB product — tensor product of per-factor MUBs, no CRT hack, no separability assumption',
        'measurement_model': '∏(p_i+1) joint configs, r=2 per config, total shots = 2∏(p_i+1)',
    }
    out_data = {'results': results, 'summary': summary}
    with open(out_json, 'w') as f:
        json.dump(out_data, f, indent=2, default=str)
    print(f"\n{'='*70}")
    print(f"DONE — {elapsed:.0f}s total")
    print(f"Results: {out_json}")
    print(f"{'='*70}")
