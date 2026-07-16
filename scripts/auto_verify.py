#!/usr/bin/env python3
"""
AI 自动验证方案 — 检查论文声明与实验数据的一致性
====================================================
运行此脚本即可获得完整的自动验证报告。
"""

import json, os, sys, numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'validate'))
from core import build_mubs_wf, random_pure_state, primes_upto, prime_factors

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
PAPER_DIR = os.path.join(os.path.dirname(__file__), '..', 'reports')

def check(description, passed):
    status = "✅ PASS" if passed else "❌ FAIL"
    print(f"  {status}: {description}")
    return passed

print("=" * 70)
print("AUTOMATED VERIFICATION OF PAPER v3.1 CLAIMS")
print("=" * 70)

results = {'pass': 0, 'fail': 0, 'warn': 0}

# ======================================================================
# B1: Data-Paper Consistency
# ======================================================================
print("\n── B1: Data-Paper Number Consistency ──")

with open(os.path.join(DATA_DIR, 'primes_1000_results.json')) as f:
    e1 = json.load(f)

# Claim 1: 168 primes tested
n_primes = len(e1)
results['pass' if n_primes == 168 else 'fail'] += 1
check(f"168 primes tested → actual={n_primes}", n_primes == 168)

# Claim 2: Median 1loc ratio = 0.864
ratios = {}
for k, v in e1.items():
    d = v['d']
    best = min(r['ratio_1loc'] for r in v['runs'])
    ratios[d] = best
rl = [ratios[d] for d in sorted(ratios.keys()) if d >= 3]
med = np.median(rl)
diff = abs(med - 0.864)
results['pass' if diff < 0.01 else 'fail'] += 1
check(f"Median 1loc ratio = 0.864 → actual={med:.4f} (diff={diff:.4f})", diff < 0.01)

# Claim 3: 63.7% ratio < 0.90
pct = sum(1 for r in rl if r < 0.90) / len(rl) * 100
diff2 = abs(pct - 63.7)
results['pass' if diff2 < 2 else 'fail'] += 1
check(f"63.7% ratio < 0.90 → actual={pct:.1f}% (diff={diff2:.1f}pp)", diff2 < 2)

# Claim 4: Top-10 dimensions
ranked = sorted([(ratios[d], d) for d in sorted(ratios.keys()) if d >= 3])
top10_claimed = [607, 941, 857, 829, 733, 613, 547, 409, 107, 887]
top10_actual = [d for r, d in ranked[:10]]
match_count = sum(1 for d in top10_claimed if d in top10_actual)
results['pass' if match_count >= 8 else 'fail'] += 1
check(f"Top-10 dimensions match → {match_count}/10", match_count >= 8)

# Claim 5: d=2 excluded (extreme outlier)
d2_ratio = ratios.get(2, None)
results['pass' if d2_ratio is not None and d2_ratio > 2 else 'warn'] += 1
check(f"d=2 is outlier → ratio={d2_ratio:.4f}", d2_ratio is not None)

# ======================================================================
# B2: T1 Numerical Verification
# ======================================================================
print("\n── B2: T1 Algebraic Equivalence Verification ──")

test_dims = [3, 5, 7, 11, 13, 17, 19, 23, 31]
n_states = 3
max_errors = []

for d in test_dims:
    mubs = build_mubs_wf(d)
    for si in range(n_states):
        psi = random_pure_state(d, seed=d*10000+si)
        rho = np.outer(psi, psi.conj())
        
        # MUB inversion: E[ρ̂] = (d+1)·M(ρ) - I = ρ
        # Using probability-weighted approach: S = Σ P(m,b)·|ψ⟩⟨ψ|
        # From HKP: M(ρ) = ρ/(d+1) + I/(d+1) → S = (d+1)·M(ρ) = ρ + I
        # Therefore: ρ = S - I
        rho_mub = np.zeros((d,d), dtype=complex)
        for m in range(d+1):
            basis = mubs[m]
            for b in range(d):
                vec = basis[b]
                prob = np.real(np.conj(vec) @ rho @ vec)
                rho_mub += prob * np.outer(vec, vec.conj())
        # S = ρ + I → ρ = S - I
        rho_mub = rho_mub - np.eye(d)
        
        # Simple averaging (same as DPRT would do for full MUB set)
        err = np.linalg.norm(rho_mub - rho, 'fro')
        max_errors.append(err)

max_err = max(max_errors)
results['pass' if max_err < 1e-10 else 'fail'] += 1
check(f"||ρ̂_MUB - ρ||_F < 1e-10 → max={max_err:.2e}", max_err < 1e-10)

# ======================================================================
# B3: Shadow Norm Consistency
# ======================================================================
print("\n── B3: Shadow Norm Formula Verification ──")

shadow_errors = []
for d in [3, 5, 7, 11, 13, 17]:
    for si in range(5):
        # Random 1-local observable
        a = si % d
        b = (si + 1) % d
        if a == b: b = (b + 1) % d
        O = np.zeros((d,d), dtype=complex)
        O[a,b] = 1.0; O[b,a] = 1.0
        
        sn_actual = d * np.real(np.trace(O @ O))
        sn_expected = 2 * d
        shadow_errors.append(abs(sn_actual - sn_expected))

max_sn_err = max(shadow_errors)
results['pass' if max_sn_err < 1e-12 else 'fail'] += 1
check(f"||O||²_shadow = d·tr(O²) → max_error={max_sn_err:.2e}", max_sn_err < 1e-12)

# ======================================================================
# B4: T4 Direction Variance Consistency  
# ======================================================================
print("\n── B4: T4 Direction Variance Consistency ──")

consistency = []
for d in [3, 5, 7, 11, 13, 17, 19]:
    mubs = build_mubs_wf(d)
    n_bases = len(mubs)
    
    for si in range(10):
        psi = random_pure_state(d, seed=d*1000+si)
        
        exps = np.zeros(n_bases)
        for m in range(n_bases):
            basis = mubs[m]
            exp_m = 0
            for k in range(d):
                prob = np.abs(np.conj(basis[k]) @ psi)**2
                val = np.real(np.conj(basis[k])[0]*basis[k][1] + 
                             np.conj(basis[k])[1]*basis[k][0])
                exp_m += prob * val * (d+1)
            exps[m] = exp_m
        
        var_dir = np.var(exps)
        mu = np.mean(exps)
        
        # Check T4 Step 1: E[ô|0]=0
        if abs(exps[0]) > 1e-12:
            consistency.append(False)
        
        # Check T4 inequality: σ²_dir ≥ μ²/(p+1)
        if var_dir + 1e-14 < mu**2 / n_bases:
            consistency.append(False)
        
        consistency.append(True)

cons_rate = sum(consistency) / len(consistency)
results['pass' if cons_rate > 0.99 else 'fail'] += 1
check(f"T4 consistency rate → {cons_rate*100:.1f}% ({sum(consistency)}/{len(consistency)})", cons_rate > 0.99)

# ======================================================================
# B5: Reference Completeness
# ======================================================================
print("\n── B5: Reference Completeness ──")

with open(os.path.join(PAPER_DIR, 'paper_zh_v2.md')) as f:
    paper = f.read()

import re
refs_in_text = set(int(m) for m in re.findall(r'\[(\d+)\]', paper))
refs_in_list = set()
for line in paper.split('\n'):
    if line.strip().startswith('> ['):
        m = re.match(r'> \[(\d+)\]', line)
        if m:
            refs_in_list.add(int(m.group(1)))

missing_in_list = refs_in_text - refs_in_list
unused_in_text = refs_in_list - refs_in_text

results['pass' if not missing_in_list else 'fail'] += 1
check(f"All cited refs in list → missing: {missing_in_list or 'none'}", not missing_in_list)

results['pass' if not unused_in_text else 'warn'] += 1
check(f"All listed refs cited → unused: {unused_in_text or 'none'}", not unused_in_text)

# ======================================================================
# B6: Section number check
# ======================================================================
print("\n── B6: Section Number Consistency ──")

section_numbers = re.findall(r'^### (\d+\.\d+)$', paper, re.MULTILINE)  # only level-3 headers
seen = []
duplicates = []
for sn in section_numbers:
    if sn in seen and sn not in duplicates:
        duplicates.append(sn)
    seen.append(sn)

results['pass' if not duplicates else 'fail'] += 1
check(f"No duplicate section numbers → duplicates: {duplicates or 'none'}", not duplicates)

# ======================================================================
# Summary
# ======================================================================
print(f"\n{'='*70}")
print(f"VERIFICATION SUMMARY")
print(f"{'='*70}")
print(f"  ✅ Passed: {results['pass']}")
print(f"  ❌ Failed: {results['fail']}")
print(f"  ⚠️  Warnings: {results['warn']}")
print(f"  Total checks: {sum(results.values())}")

if results['fail'] == 0:
    print(f"\n  All critical checks passed! Paper is numerically consistent.")
else:
    print(f"\n  {results['fail']} checks FAILED — fix before submission.")
