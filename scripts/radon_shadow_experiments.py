#!/usr/bin/env python3
"""
RadonShadow Experiment: Phase 1.5 & Phase 2
============================================
- 2-qubit fair comparison (clean rewrite)
- qutrit (d=3) comparison: 4 MUBs vs random
- DPRT post-processing acceleration demo

Author: QClaw (2026-07-14)
"""

import numpy as np
from scipy.stats import unitary_group
import itertools
from collections import defaultdict
import json
import time

# ======================================================================
# 1. Constants
# ======================================================================
I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)
PAULI = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

Z0 = np.array([1,0], dtype=complex)
Z1 = np.array([0,1], dtype=complex)
X0 = np.array([1,1], dtype=complex)/np.sqrt(2)
X1 = np.array([1,-1], dtype=complex)/np.sqrt(2)
Y0 = np.array([1,1j], dtype=complex)/np.sqrt(2)
Y1 = np.array([1,-1j], dtype=complex)/np.sqrt(2)

BASIS_PROJECTORS = {
    'Z': [np.outer(Z0, Z0.conj()), np.outer(Z1, Z1.conj())],
    'X': [np.outer(X0, X0.conj()), np.outer(X1, X1.conj())],
    'Y': [np.outer(Y0, Y0.conj()), np.outer(Y1, Y1.conj())],
}
BASES = ['Z', 'X', 'Y']

# ======================================================================
# 2. Qutrit MUBs (Wootters & Fields 1989)
# ======================================================================
def make_qutrit_mubs():
    omega = np.exp(2j * np.pi / 3)
    # M0: standard
    m0 = [np.array(v, dtype=complex) for v in [[1,0,0],[0,1,0],[0,0,1]]]
    # M1: Fourier
    F = np.array([[1,1,1],[1,omega,omega**2],[1,omega**2,omega]], dtype=complex)/np.sqrt(3)
    m1 = [F[:,j] for j in range(3)]
    # M2: diagonal-twisted Fourier
    D = np.diag([1, omega, omega**2])
    F2 = F @ D
    m2 = [F2[:,j] for j in range(3)]
    # M3: squared-diagonal
    D2 = np.diag([1, omega**2, omega])
    F3 = F @ D2
    m3 = [F3[:,j] for j in range(3)]
    return {'M0': m0, 'M1': m1, 'M2': m2, 'M3': m3}

QUTRIT_MUBS = make_qutrit_mubs()
QUTRIT_PROJ = {k: [np.outer(v, v.conj()) for v in vecs] for k, vecs in QUTRIT_MUBS.items()}
QUTRIT_LABELS = ['M0', 'M1', 'M2', 'M3']

# Gell-Mann matrices (traceless, tr(Li Lj)=2δij)
GELL_MANN = {
    'λ1': np.array([[0,1,0],[1,0,0],[0,0,0]], dtype=complex),
    'λ2': np.array([[0,-1j,0],[1j,0,0],[0,0,0]], dtype=complex),
    'λ3': np.array([[1,0,0],[0,-1,0],[0,0,0]], dtype=complex),
    'λ4': np.array([[0,0,1],[0,0,0],[1,0,0]], dtype=complex),
    'λ5': np.array([[0,0,-1j],[0,0,0],[1j,0,0]], dtype=complex),
    'λ6': np.array([[0,0,0],[0,0,1],[0,1,0]], dtype=complex),
    'λ7': np.array([[0,0,0],[0,0,-1j],[0,1j,0]], dtype=complex),
    'λ8': np.array([[1,0,0],[0,1,0],[0,0,-2]], dtype=complex)/np.sqrt(3),
}

# ======================================================================
# 3. State generation
# ======================================================================
def random_pure_state(dim):
    z = np.random.randn(dim) + 1j * np.random.randn(dim)
    return z / np.linalg.norm(z)

def random_mixed_state(dim, rank=None):
    if rank is None: rank = dim
    d = np.random.dirichlet(np.ones(rank))
    d = np.pad(d, (0, dim - len(d)), constant_values=0)
    U = unitary_group.rvs(dim)
    return U @ np.diag(d) @ U.conj().T

# ======================================================================
# 4. 2-qubit shadow snapshots
# ======================================================================
def snap_2qubit(rho, b0, b1):
    """Classical shadow snapshot, 2-qubit local Pauli."""
    # Partial traces
    rho_A = rho[0::2, 0::2] + rho[1::2, 1::2]
    rho_B = rho[0:2, 0:2] + rho[2:4, 2:4]
    # qubit 0
    p0 = np.real(np.trace(rho_A @ BASIS_PROJECTORS[b0][0]))
    o0 = np.random.choice([0,1], p=[np.clip(p0,0,1), 1-np.clip(p0,0,1)])
    s0 = 3 * BASIS_PROJECTORS[b0][o0] - I2
    # qubit 1
    p1 = np.real(np.trace(rho_B @ BASIS_PROJECTORS[b1][0]))
    o1 = np.random.choice([0,1], p=[np.clip(p1,0,1), 1-np.clip(p1,0,1)])
    s1 = 3 * BASIS_PROJECTORS[b1][o1] - I2
    return np.kron(s0, s1)

# ======================================================================
# 5. Qutrit shadow snapshots
# ======================================================================
def snap_qutrit(rho, basis_label):
    projs = QUTRIT_PROJ[basis_label]
    probs = np.clip([np.real(np.trace(rho @ P)) for P in projs], 0, None)
    probs = probs / probs.sum()
    o = np.random.choice([0,1,2], p=probs)
    return 4 * projs[o] - np.eye(3, dtype=complex)

# ======================================================================
# 6. Observable prediction
# ======================================================================
def predict(snapshots, observable):
    return np.mean([np.real(np.trace(observable @ s)) for s in snapshots])

# ======================================================================
# 7. 2-Qubit Experiment
# ======================================================================
def run_2qubit(n_states=100, r_list=[1,2,4,8], seed=42):
    dim = 4
    obs = {}
    for labels in itertools.product(['I','X','Y','Z'], repeat=2):
        o = np.eye(1, dtype=complex)
        for l in labels: o = np.kron(o, PAULI[l])
        obs[''.join(labels)] = o
    obs_labels = [l for l in obs if l != 'II']
    
    results = []
    print("="*60)
    print("2-QUBIT: DETERMINISTIC (9 configs) vs RANDOM")
    print("="*60)
    
    for r in r_list:
        n_meas = 9 * r
        det_err = defaultdict(list)
        rnd_err = defaultdict(list)
        diffs = []
        
        for si in range(n_states):
            s_seed = seed + si * 1000 + r * 100
            np.random.seed(s_seed)
            psi = random_pure_state(dim)
            rho = np.outer(psi, psi.conj())
            truth = {l: np.real(np.trace(obs[l] @ rho)) for l in obs_labels}
            
            # Deterministic
            snaps = [snap_2qubit(rho, b0, b1) for _ in range(r) 
                     for b0 in ['Z','X','Y'] for b1 in ['Z','X','Y']]
            for l in obs_labels:
                det_err[l].append(abs(predict(snaps, obs[l]) - truth[l]))
            
            # Random
            snaps = [snap_2qubit(rho, np.random.choice(BASES), np.random.choice(BASES))
                     for _ in range(n_meas)]
            for l in obs_labels:
                rnd_err[l].append(abs(predict(snaps, obs[l]) - truth[l]))
            
            # Paired diff for t-test
            for l in obs_labels:
                diffs.append(det_err[l][-1] - rnd_err[l][-1])
        
        det_mae = np.mean([np.mean(det_err[l]) for l in obs_labels])
        rnd_mae = np.mean([np.mean(rnd_err[l]) for l in obs_labels])
        t_stat = np.mean(diffs) / (np.std(diffs)/np.sqrt(len(diffs)))
        
        # 1-local vs 2-local split
        local1 = [l for l in obs_labels if l.count('I') == 1]
        local2 = [l for l in obs_labels if l.count('I') == 0]
        d1 = np.mean([np.mean(det_err[l]) for l in local1])
        r1 = np.mean([np.mean(rnd_err[l]) for l in local1])
        d2 = np.mean([np.mean(det_err[l]) for l in local2])
        r2 = np.mean([np.mean(rnd_err[l]) for l in local2])
        
        winner = "DET" if det_mae < rnd_mae else "RAND"
        print(f"\nr={r:>2d}  N={n_meas:>3d}  "
              f"Det={det_mae:.5f}  Rand={rnd_mae:.5f}  "
              f"Ratio={det_mae/rnd_mae:.4f}  t={t_stat:+.3f}  [{winner}]")
        print(f"  1-local: det={d1:.5f} rand={r1:.5f}  r={d1/r1:.3f}")
        print(f"  2-local: det={d2:.5f} rand={r2:.5f}  r={d2/r2:.3f}")
        
        results.append({
            'r': r, 'n_meas': n_meas,
            'det_mae': float(det_mae), 'rand_mae': float(rnd_mae),
            'ratio': float(det_mae/rnd_mae), 't_stat': float(t_stat),
            'winner': winner,
            'local1_ratio': float(d1/r1), 'local2_ratio': float(d2/r2),
        })
    return results

# ======================================================================
# 8. Qutrit Experiment
# ======================================================================
def run_qutrit(n_states=50, r_list=[1,2,4,8], seed=42):
    dim = 3
    obs = GELL_MANN
    obs_labels = list(obs.keys())
    
    results = []
    print("\n" + "="*60)
    print("QUTRIT (d=3): 4-MUB DETERMINISTIC vs RANDOM MUB")
    print("="*60)
    
    for r in r_list:
        n_meas = 4 * r
        det_err = defaultdict(list)
        rnd_err = defaultdict(list)
        
        for si in range(n_states):
            s_seed = seed + 5000 + si * 1000 + r * 100
            np.random.seed(s_seed)
            # Use pure states for fair comparison with 2-qubit
            psi = random_pure_state(dim)
            rho = np.outer(psi, psi.conj())
            truth = {l: np.real(np.trace(O @ rho)) for l, O in obs.items()}
            
            # Deterministic: all 4 MUBs × r rounds
            snaps = [snap_qutrit(rho, ml) for _ in range(r) for ml in QUTRIT_LABELS]
            for l in obs_labels:
                det_err[l].append(abs(predict(snaps, obs[l]) - truth[l]))
            
            # Random
            snaps = [snap_qutrit(rho, np.random.choice(QUTRIT_LABELS)) for _ in range(n_meas)]
            for l in obs_labels:
                rnd_err[l].append(abs(predict(snaps, obs[l]) - truth[l]))
        
        det_mae = np.mean([np.mean(det_err[l]) for l in obs_labels])
        rnd_mae = np.mean([np.mean(rnd_err[l]) for l in obs_labels])
        winner = "DET" if det_mae < rnd_mae else "RAND"
        
        print(f"\nr={r:>2d}  N={n_meas:>3d}  "
              f"Det={det_mae:.5f}  Rand={rnd_mae:.5f}  "
              f"Ratio={det_mae/rnd_mae:.4f}  [{winner}]")
        for l in obs_labels:
            print(f"  {l}: det={np.mean(det_err[l]):.5f}  rand={np.mean(rnd_err[l]):.5f}  "
                  f"r={np.mean(det_err[l])/np.mean(rnd_err[l]):.3f}")
        
        results.append({
            'r': r, 'n_meas': n_meas,
            'det_mae': float(det_mae), 'rand_mae': float(rnd_mae),
            'ratio': float(det_mae/rnd_mae), 'winner': winner,
        })
    return results

# ======================================================================
# 9. DPRT Post-Processing (MUB → Density Matrix)
# ======================================================================
def dprt_inverse_3x3(p0, p1, p2, p3):
    """
    DPRT inverse for d=3.
    Forward: R[m,k] = Σ_i f[i, (k + m*i) mod 3], m=0,1,2
             R[3,i] = Σ_j f[i,j]  (row sums)
    Inverse: f[i,j] = (R[0,j] + R[1,(j-i)%3] + R[2,(j-2i)%3] + R[3,i] - total) / 3
    """
    total = np.sum(p0)
    f = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            s = p3[i]  # row sum of row i
            s += p0[j] + p1[(j - i) % 3] + p2[(j - 2*i) % 3]
            s -= total
            f[i, j] = s / 3.0
    return f

def test_dprt_correctness():
    """Verify DPRT forward/inverse correctness."""
    np.random.seed(42)
    f = np.random.randint(0, 100, (3, 3)).astype(float)
    d = 3
    
    # Forward DPRT
    R = np.zeros((d+1, d))
    for m in range(d):
        for k in range(d):
            s = 0
            for i in range(d):
                j = (k + m*i) % d
                s += f[i, j]
            R[m, k] = s
    for i in range(d):
        R[d, i] = np.sum(f[i, :])
    
    # Inverse
    frec = dprt_inverse_3x3(R[0], R[1], R[2], R[3])
    err = np.max(np.abs(f - frec))
    print(f"\nDPRT inversion correctness test: max error = {err:.2e}")
    return err < 1e-10

# ======================================================================
# 10. 2-Qubit Mixed State Experiment
# ======================================================================
def run_2qubit_mixed(n_states=50, r_list=[1,2,4,8], seed=999):
    """2-qubit + MIXED states (vs pure-state experiment to isolate state-type effect)."""
    dim = 4
    obs = {}
    for labels in itertools.product(['I','X','Y','Z'], repeat=2):
        o = np.eye(1, dtype=complex)
        for l in labels: o = np.kron(o, PAULI[l])
        obs[''.join(labels)] = o
    obs_labels = [l for l in obs if l != 'II']
    
    results = []
    print("\n" + "="*60)
    print("2-QUBIT MIXED: Deterministic vs Random (mixed states)")
    print("="*60)
    
    for r in r_list:
        n_meas = 9 * r
        det_err = defaultdict(list)
        rnd_err = defaultdict(list)
        
        for si in range(n_states):
            s_seed = seed + si * 1000 + r * 100
            np.random.seed(s_seed)
            rho = random_mixed_state(dim, rank=np.random.randint(1, 5))
            truth = {l: np.real(np.trace(obs[l] @ rho)) for l in obs_labels}
            
            # Deterministic
            snaps = [snap_2qubit(rho, b0, b1) for _ in range(r) 
                     for b0 in ['Z','X','Y'] for b1 in ['Z','X','Y']]
            for l in obs_labels:
                det_err[l].append(abs(predict(snaps, obs[l]) - truth[l]))
            
            # Random
            snaps = [snap_2qubit(rho, np.random.choice(BASES), np.random.choice(BASES))
                     for _ in range(n_meas)]
            for l in obs_labels:
                rnd_err[l].append(abs(predict(snaps, obs[l]) - truth[l]))
        
        det_mae = np.mean([np.mean(det_err[l]) for l in obs_labels])
        rnd_mae = np.mean([np.mean(rnd_err[l]) for l in obs_labels])
        winner = "DET" if det_mae < rnd_mae else "RAND"
        
        print(f"r={r:>2d}  N={n_meas:>3d}  "
              f"Det={det_mae:.5f}  Rand={rnd_mae:.5f}  "
              f"Ratio={det_mae/rnd_mae:.4f}  [{winner}]")
        
        results.append({
            'r': r, 'n_meas': n_meas,
            'det_mae': float(det_mae), 'rand_mae': float(rnd_mae),
            'ratio': float(det_mae/rnd_mae), 'winner': winner,
        })
    return results

# ======================================================================
# 11. 3-Qubit Experiment
# ======================================================================
def snap_3qubit(rho, b0, b1, b2):
    """3-qubit local Pauli shadow snapshot."""
    bases = [b0, b1, b2]
    snaps = []
    n = 3
    for q in range(3):
        rho_q = np.zeros((2,2), dtype=complex)
        shift = 2 - q
        for a in [0,1]:
            for b in [0,1]:
                for others in range(2**(n-1)):
                    bits = [(others >> k) & 1 for k in range(n-1)]
                    bits.insert(shift, a)
                    idx_a = sum(bits[k] << (n-1-k) for k in range(n))
                    bits[shift] = b
                    idx_b = sum(bits[k] << (n-1-k) for k in range(n))
                    rho_q[a, b] += rho[idx_a, idx_b] if a <= b else np.conj(rho[idx_b, idx_a])
        basis = bases[q]
        p0 = np.real(np.trace(rho_q @ BASIS_PROJECTORS[basis][0]))
        p0 = np.clip(p0, 0, 1)
        o = np.random.choice([0,1], p=[p0, 1-p0])
        snaps.append(3 * BASIS_PROJECTORS[basis][o] - I2)
    snap = np.eye(1, dtype=complex)
    for s in snaps: snap = np.kron(snap, s)
    return snap

def run_3qubit(n_states=30, r_list=[1,2,4], seed=777):
    """3-qubit pure states: 27 configs (3^3) × r rounds vs random."""
    dim = 8
    obs = {}
    for labels in itertools.product(['I','X','Y','Z'], repeat=3):
        o = np.eye(1, dtype=complex)
        for l in labels: o = np.kron(o, PAULI[l])
        obs[''.join(labels)] = o
    obs_labels = [l for l in obs if l not in ['III']]
    
    all_configs = list(itertools.product(BASES, repeat=3))
    
    results = []
    print("\n" + "="*60)
    print("3-QUBIT: DETERMINISTIC (27 configs) vs RANDOM")
    print("="*60)
    
    for r in r_list:
        n_meas = 27 * r
        det_err = defaultdict(list)
        rnd_err = defaultdict(list)
        
        for si in range(n_states):
            s_seed = seed + si * 1000 + r * 100
            np.random.seed(s_seed)
            psi = random_pure_state(dim)
            rho = np.outer(psi, psi.conj())
            truth = {l: np.real(np.trace(obs[l] @ rho)) for l in obs_labels}
            
            # Deterministic
            snaps = [snap_3qubit(rho, b0, b1, b2) for _ in range(r)
                     for b0, b1, b2 in all_configs]
            for l in obs_labels:
                det_err[l].append(abs(predict(snaps, obs[l]) - truth[l]))
            
            # Random (same N)
            snaps_r = [snap_3qubit(rho, np.random.choice(BASES),
                                    np.random.choice(BASES), np.random.choice(BASES))
                       for _ in range(n_meas)]
            for l in obs_labels:
                rnd_err[l].append(abs(predict(snaps_r, obs[l]) - truth[l]))
        
        det_mae = np.mean([np.mean(det_err[l]) for l in obs_labels])
        rnd_mae = np.mean([np.mean(rnd_err[l]) for l in obs_labels])
        winner = "DET" if det_mae < rnd_mae else "RAND"
        
        print(f"r={r:>2d}  N={n_meas:>3d}  "
              f"Det={det_mae:.5f}  Rand={rnd_mae:.5f}  "
              f"Ratio={det_mae/rnd_mae:.4f}  [{winner}]")
        
        # Break down by locality for insight
        for ll, gname in [(2, '1-local'), (1, '2-local'), (0, '3-local')]:
            grp = [l for l in obs_labels if l.count('I') == ll]
            if not grp: continue
            d = np.mean([np.mean(det_err[l]) for l in grp])
            rd = np.mean([np.mean(rnd_err[l]) for l in grp])
            print(f"  {gname}: det={d:.5f} rand={rd:.5f} r={d/rd:.3f}")
        
        results.append({
            'r': r, 'n_meas': n_meas,
            'det_mae': float(det_mae), 'rand_mae': float(rnd_mae),
            'ratio': float(det_mae/rnd_mae), 'winner': winner,
        })
    return results

# ======================================================================
# 12. Shadow Norm Analysis
# ======================================================================
def shadow_norm_analytic(O, d=2):
    """
    ||O||²_shadow for local Pauli ensemble.
    For k-local Pauli observable P: ||P||²_shadow = 3^k
    For general O with Pauli decomposition O = Σ α_P P:
    ||O||²_shadow = Σ 3^k_P |α_P|²
    where k_P is the locality (weight) of Pauli string P.
    """
    from itertools import product
    
    n = int(np.log2(O.shape[0]))
    pauli_basis = {}
    for labels in product(['I','X','Y','Z'], repeat=n):
        o = np.eye(1, dtype=complex)
        for l in labels: o = np.kron(o, PAULI[l])
        pauli_basis[''.join(labels)] = o
    
    shadow_sq = 0
    for label, P in pauli_basis.items():
        coeff = np.real(np.trace(P @ O)) / O.shape[0]
        k = sum(1 for c in label if c != 'I')
        shadow_sq += 3**k * coeff**2
    
    return shadow_sq

def analyze_shadow_norms():
    """Compare shadow norms of deterministic vs random: theoretical analysis."""
    print("\n" + "="*60)
    print("SHADOW NORM ANALYSIS")
    print("="*60)
    
    # For the inverse channel map M_pauli:
    # M⁻¹ : A → ⊗ⱼ (3 Tr_Bⱼ(A) ⊗ I/2 - I/2 Tr(A)/2)
    # Wait, this is wrong. Let's be very careful.
    
    # HKP 2020 local Pauli:
    # ρ̂ = ⊗ⱼ (3 |bⱼ⟩⟨bⱼ| - I)
    # E[ρ̂] = ρ  (unbiased)
    # Shadow norm: ||O||²_shadow = ||E_{U,b} [⟨b|U M⁻¹(O) U†|b⟩²]||
    
    # For Pauli observables P (k-local):
    # ||P||²_shadow = 3^k
    # This is the key result from HKP Eq (S16).
    
    # Now compare:
    # Deterministic: N = 3^n · r measurements, each config used r times.
    # Random: N = 3^n · r measurements, random configs.
    
    # The variance for observable O with N snapshots:
    # Var(ô) = (1/N) · (||O||²_shadow - tr(ρ O)²)
    #            + (1/N²) · Σ_{i≠j} Cov(ô_i, ô_j)
    #
    # For RANDOM: Cov term = 0 (independent samples)
    # For DETERMINISTIC: Cov(ô_i, ô_j) ≠ 0 because configs repeat!
    #   But the covariance depends on the observable and state.
    
    # Let's compute the exact variance for a 2-qubit example.
    
    dim = 4  # 2-qubit
    
    # Pauli observables
    obs = {}
    for labels in itertools.product(['I','X','Y','Z'], repeat=2):
        o = np.eye(1, dtype=complex)
        for l in labels: o = np.kron(o, PAULI[l])
        obs[''.join(labels)] = o
    
    print(f"\n  Observable  │ k  │ ||O||²_shadow │ tr(O²) │ ratio")
    print(f"  {'─'*12}│{'─'*4}│{'─'*15}│{'─'*9}│{'─'*7}")
    for label, O in obs.items():
        k = sum(1 for c in label if c != 'I')
        if k == 0: continue
        ssq = shadow_norm_analytic(O, d=2)
        trO2 = np.real(np.trace(O @ O))
        print(f"  {label:>12s} │ {k:>2d} │ {ssq:>15.3f} │ {trO2:>9.3f} │ {ssq/trO2:>.3f}")
    
    # Deterministic trace sampling: repeat same config r times
    # For r repeats of config c, the r shadow snapshots are NOT independent
    # because they share the same basis choice.
    # Cov(ρ̂_i(c), ρ̂_j(c)) ~ O(1 - δ_{ij}) when same config.
    #
    # Key: In deterministic protocol, N_det = 3^n * r, but we have
    # only 3^n independent configs, each used r times.
    # The effective sample size is closer to 3^n, not 3^n * r.
    #
    # For random protocol with same N: all N are i.i.d.
    # So for the same N, random has higher effective sample size!
    #
    # However: deterministic covers all bases equally, which
    # guarantees good coverage of the measurement space.
    # For observables with support on specific bases, this helps.
    
    print(f"\n  Theoretical insight:")
    print(f"  - Random: Var ~ (1/N) · ||O||²_shadow (i.i.d., effective N = actual N)")
    print(f"  - Deterministic: Var ~ (1/N_det) · ||O||²_shadow · (1 + r·f)")
    print(f"    where f = basis-covariance factor (0 for Pauli-commuting obs, >0 for others)")
    print(f"  - Trade-off: det eliminates basis-sampling variance but reduces effective N")

# ======================================================================
# 13. All Experiments
# ======================================================================
if __name__ == '__main__':
    t0 = time.time()

    results_2q = run_2qubit(n_states=200, r_list=[1, 2, 4, 8])
    results_q = run_qutrit(n_states=60, r_list=[1, 2, 4, 8])
    results_2q_mixed = run_2qubit_mixed(n_states=60, r_list=[1, 2, 4, 8])
    results_3q = run_3qubit(n_states=20, r_list=[1, 2, 4])

    test_dprt_correctness()
    analyze_shadow_norms()

    t1 = time.time()
    print(f"\nTotal: {t1-t0:.1f}s")

    # Master comparison table
    print("\n" + "="*60)
    print("MASTER COMPARISON: Det/Rand Ratio (ratio<1 = DET better)")
    print("="*60)
    print(f"\n{'N_meas':>7s} | {'2Q Pure':>9s} | {'2Q Mixed':>9s} | {'QuT Pure':>9s} | {'3Q Pure':>9s}")
    print("-" * 47)
    for i in range(max(len(d) for d in [results_2q, results_2q_mixed, results_q, results_3q])):
        row = []
        for ds in [results_2q, results_2q_mixed, results_q, results_3q]:
            if i < len(ds):
                row.append(f"{ds[i]['ratio']:>9.4f}")
            else:
                row.append(f"{'─':>9s}")
        if i < len(results_2q):
            print(f"{results_2q[i]['n_meas']:>7d} | {row[0]} | {row[1]} | {row[2]} | {row[3]}")

    all_results = {
        'two_qubit_pure': results_2q,
        'two_qubit_mixed': results_2q_mixed,
        'qutrit_pure': results_q,
        'three_qubit': results_3q,
    }
    with open('qdprt_experiment_phase2.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print("\n→ qdprt_experiment_phase2.json")
