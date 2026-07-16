#!/usr/bin/env python3
"""
Path B: 解析证明 σ²_dir > 0 → ratio < 1
============================================
核心洞察: 对1-local可观测量O_ab(a≠b), 计算基方向E[ô|m=p]=0,
         Fourier基方向E[ô|m=0]≠0 → σ²_dir > 0 严格成立。

三步证明:
1. 证明 E[ô|m=p] = 0（平凡，计算基不访问非对角元）
2. 证明 E[ô|m=0] ≠ 0 对几乎所有纯态成立（Fourier基包含全信息）
3. 下界: σ²_dir ≥ (E[ô|0] - μ)²/(p+1) > 0
"""

import numpy as np
from collections import defaultdict

print("="*70)
print("THEOREM: σ²_dir > 0 FOR ALL PRIME p, GENERIC PURE STATES")
print("="*70)

print("""
Proof structure:

Step 1: E[ô|p] = 0
  Computational basis |e^{(p)}_k⟩ = |k⟩
  For O_ab = |a⟩⟨b|+|b⟩⟨a| with a≠b:
  ⟨k|O_ab|k⟩ = δ_{ka}δ_{kb} + δ_{kb}δ_{ka} = 0  (a≠b, can't have k=a AND k=b)
  ∴ E[ô|p] = 0

Step 2: E[ô|0] ≠ 0 generically
  Fourier basis |e^{(0)}_k⟩ = (1/√p) Σ_j ω^{-k·j} |j⟩
  ⟨e^{(0)}_k|O_ab|e^{(0)}_k⟩ = (2/p) cos(2πk(b-a)/p)
  This is non-zero for k ≠ 0, p/2 when (b-a) ≠ 0 (mod p)
  Since b≠a and 0<b-a<p, there exists k such that value ≠ 0.
  
  For Haar-random |ψ⟩, the probability distribution P(0,k) = |⟨e^{(0)}_k|ψ⟩|²
  is uniform on average. Therefore:
  E[ô|0] = (p+1)·(2/p)·Σ_k P(0,k)·cos(2πk(b-a)/p) ≠ 0 almost surely.
  
  The set of |ψ⟩ for which E[ô|0] = 0 has measure zero
  (it requires exact cancellation of the weighted cosines).

Step 3: Lower bound on σ²_dir
  μ = (1/(p+1)) Σ_m E[ô|m]  (mean of conditional expectations)
  Since E[ô|p] = 0 and E[ô|0] = v₀ ≠ 0:
  μ = (v₀ + Σ_{m≠0,p} E[ô|m])/(p+1)
  
  σ²_dir = (1/(p+1)) Σ_m (E[ô|m] - μ)²
         ≥ (E[ô|p] - μ)²/(p+1)       [by dropping all but m=p term]
         = μ²/(p+1)                    [since E[ô|p]=0]
  
  If μ = 0, then v₀ + Σ_{m≠0,p} E[ô|m] = 0 ⇒ v₀ = -Σ_{m≠0,p} E[ô|m]
  But v₀ depends on |ψ⟩ through P(0,k), while the RHS depends on |ψ⟩ through 
  different MUB directions. For Haar-random |ψ⟩, the probability that this exact
  cancellation occurs is measure zero.
  
  Therefore μ ≠ 0 almost surely ⇒ σ²_dir ≥ μ²/(p+1) > 0. ∎

Corollary (T4): For any prime p and generic pure state |ψ⟩,
  ratio = 1/(1 + σ²_dir/||O||²_shadow) < 1
  → DPRT strictly outperforms random MUB sampling.
  
  The ratio is approximately 1/(1 + 0.7/p²) based on numerical fitting [§3.6].
""")

# 数值验证上述证明的每一步
print("\n" + "="*70)
print("NUMERICAL VERIFICATION OF THE PROOF")
print("="*70)

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'validate'))
from core import build_mubs_wf, random_pure_state

def verify_step1(d=7):
    """验证 E[ô|p] = 0"""
    mubs = build_mubs_wf(d)
    basis_p = mubs[d]  # computational basis
    
    violations = 0
    for si in range(100):
        psi = random_pure_state(d, seed=si)
        for a in range(d):
            for b in range(a+1, d):
                exp_p = 0
                for k in range(d):
                    prob = np.abs(np.conj(basis_p[k]) @ psi)**2
                    val = np.real(np.conj(basis_p[k])[a] * basis_p[k][b] + 
                                 np.conj(basis_p[k])[b] * basis_p[k][a])
                    exp_p += prob * val * (d+1)
                if abs(exp_p) > 1e-12:
                    violations += 1
    print(f"  Step 1 (E[ô|p]=0): {violations} violations in {100*d*(d-1)//2} tests")

def verify_step2(d=7):
    """验证 E[ô|0] ≠ 0 generically"""
    mubs = build_mubs_wf(d)
    basis_0 = mubs[0]  # Fourier basis
    
    n_nonzero = 0
    n_test = 100
    for si in range(n_test):
        psi = random_pure_state(d, seed=si)
        exp_0 = 0
        for k in range(d):
            prob = np.abs(np.conj(basis_0[k]) @ psi)**2
            val = (2/d) * np.cos(2*np.pi*k/d)  # for a=0,b=1
            exp_0 += prob * val * (d+1)
        if abs(exp_0) > 1e-6:
            n_nonzero += 1
    
    print(f"  Step 2 (E[ô|0]≠0): {n_nonzero}/{n_test} non-zero (should be ~100)")

def verify_step3(d=7):
    """验证 σ²_dir ≥ μ²/(p+1)"""
    mubs = build_mubs_wf(d)
    n_bases = len(mubs)
    
    n_violations = 0
    n_test = 100
    for si in range(n_test):
        psi = random_pure_state(d, seed=si)
        
        exps = np.zeros(n_bases)
        for m in range(n_bases):
            basis = mubs[m]
            exp_m = 0
            for k in range(d):
                prob = np.abs(np.conj(basis[k]) @ psi)**2
                val = np.real(np.conj(basis[k])[0] * basis[k][1] + 
                             np.conj(basis[k])[1] * basis[k][0])
                exp_m += prob * val * (d+1)
            exps[m] = exp_m
        
        var_dir = np.var(exps)
        mu = np.mean(exps)
        
        if mu**2 / n_bases > var_dir + 1e-12:
            n_violations += 1
    
    print(f"  Step 3 (σ²_dir ≥ μ²/(p+1)): {n_violations}/{n_test} violations (should be 0)")

    # 更重要的: μ²/(p+1) 与 σ²_dir 的实际比例
    ratios = []
    for si in range(n_test):
        psi = random_pure_state(d, seed=si)
        exps = np.zeros(n_bases)
        for m in range(n_bases):
            basis = mubs[m]
            exp_m = 0
            for k in range(d):
                prob = np.abs(np.conj(basis[k]) @ psi)**2
                val = np.real(np.conj(basis[k])[0] * basis[k][1] + 
                             np.conj(basis[k])[1] * basis[k][0])
                exp_m += prob * val * (d+1)
            exps[m] = exp_m
        var_dir = np.var(exps)
        mu = np.mean(exps)
        if var_dir > 0:
            ratios.append(mu**2 / (n_bases * var_dir))
    
    print(f"  Tightness: μ²/(p+1) / σ²_dir = {np.median(ratios):.4f} (median)")
    print(f"    (This measures how tight the bound is — 1.0 means the bound is exact)")

verify_step1(7)
verify_step2(7)
verify_step3(7)

# 全量验证
print(f"\n  Full spectrum verification (d=3..31):")
for d in [3,5,7,11,13,17,19,23,29,31]:
    verify_step1(d)
    verify_step2(d)
    verify_step3(d)
