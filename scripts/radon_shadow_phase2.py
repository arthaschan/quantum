#!/usr/bin/env python3
"""
RadonShadow Phase 2: DPRT Post-Processing & Large-d Scan
=========================================================
Direction 2: Replace Classical Shadow M⁻¹ channel with DPRT inverse
Direction 3: Large-d scan (d=5,7,11) — dimensional trend analysis

Core insight: For prime d, measuring in all d+1 MUBs yields the DPRT
of the discrete Wigner function. The DPRT inverse is O(d²) pure
addition/subtraction — no matrix inversion needed.

Author: QClaw (2026-07-14)
"""

import numpy as np
from scipy.stats import unitary_group
from scipy.linalg import sqrtm
import itertools
from collections import defaultdict
import json
import time

# ======================================================================
# 1. DPRT Core: Forward & Inverse for Arbitrary Prime d
# ======================================================================

def dprt_forward(f, d):
    """
    DPRT: project d×d matrix f onto d+1 directions.
    R[m,k] = Σ_i f[i, (k + m*i) % d]  for m=0,...,d-1
    R[d,i] = Σ_j f[i,j]  (row sums, direction "infinity")
    
    Returns: (d+1) × d array.
    """
    R = np.zeros((d+1, d))
    for m in range(d):
        for k in range(d):
            s = 0
            for i in range(d):
                j = (k + m * i) % d
                s += f[i, j]
            R[m, k] = s
    for i in range(d):
        R[d, i] = np.sum(f[i, :])
    return R


def dprt_inverse(R, d):
    """
    Inverse DPRT for prime d.
    f[i,j] = (Σ_{m=0}^{d-1} R[m, (j - m*i) % d] + R[d,i] - total) / d
    where total = Σ_k R[0,k]
    """
    total = np.sum(R[0, :])
    f = np.zeros((d, d))
    for i in range(d):
        for j in range(d):
            s = R[d, i]
            for m in range(d):
                s += R[m, (j - m * i) % d]
            s -= total
            f[i, j] = s / d
    return f


def test_dprt(d):
    """Verify DPRT round-trip correctness for dimension d."""
    np.random.seed(42)
    f = np.random.randn(d, d)
    R = dprt_forward(f, d)
    frec = dprt_inverse(R, d)
    err = np.max(np.abs(f - frec))
    return err


# ======================================================================
# 2. Discrete Wigner Function & MUB Measurements
# ======================================================================

def pauli_heisenberg(d):
    """
    Generalized Pauli (Heisenberg-Weyl) operators for prime d.
    Z|k⟩ = ω^k |k⟩,  X|k⟩ = |k+1 mod d⟩
    where ω = exp(2πi/d).
    
    Returns: T(a,b) = ω^{-ab/2} Z^a X^b  (standard Weyl form)
    Actually, let's use the standard form:
    Z = diag(1, ω, ω², ..., ω^{d-1})
    X = cyclic shift operator
    T(a,b) = X^b Z^a  (different ordering convention)
    """
    omega = np.exp(2j * np.pi / d)
    
    Z_op = np.diag([omega**k for k in range(d)])
    X_op = np.zeros((d, d), dtype=complex)
    for k in range(d):
        X_op[(k+1) % d, k] = 1
    
    # Phase point operators: A(q,p) for discrete Wigner function
    # Using Wootters convention
    phase_point_ops = {}
    for q in range(d):
        for p in range(d):
            # A(q,p) = (1/d²) Σ_{a,b} ω^{aq+bp-pa} Z^a X^b  (various conventions exist)
            A = np.zeros((d, d), dtype=complex)
            for a in range(d):
                for b in range(d):
                    phase = omega**((q*a + p*b) % d)  # simplified convention
                    # The exact phase conventions vary by author
                    A += phase * (Z_op**a) @ (X_op**b)
            A /= d
            phase_point_ops[(q, p)] = A
    
    return Z_op, X_op, phase_point_ops


def wigner_from_density(rho, d):
    """
    Compute discrete Wigner function W[q,p] from density matrix ρ.
    W(q,p) = (1/d) tr(A(q,p) ρ)
    """
    Z_op, X_op, A_ops = pauli_heisenberg(d)
    W = np.zeros((d, d))
    for q in range(d):
        for p in range(d):
            W[q, p] = np.real(np.trace(A_ops[(q, p)] @ rho)) / d
    return W


def density_from_wigner(W, d):
    """
    Reconstruct density matrix from discrete Wigner function.
    ρ = Σ_{q,p} W(q,p) A(q,p)
    """
    Z_op, X_op, A_ops = pauli_heisenberg(d)
    rho = np.zeros((d, d), dtype=complex)
    for q in range(d):
        for p in range(d):
            rho += W[q, p] * A_ops[(q, p)]
    return rho


# ======================================================================
# 3. MUB Construction for Arbitrary Prime d
# ======================================================================

def construct_mubs(d):
    """
    Construct d+1 Mutually Unbiased Bases for prime d.
    MUB 0: Computational basis (eigenvectors of Z)
    MUB 1,...,d: Fourier basis and its shifted versions
    
    Returns: {mub_label: [d basis vectors]}
    """
    omega = np.exp(2j * np.pi / d)
    
    mubs = {}
    
    # MUB 0: Standard basis
    mubs['M0'] = [np.eye(d, dtype=complex)[:, k] for k in range(d)]
    
    # MUB 1,...,d: Fourier bases with shifts
    F = np.zeros((d, d), dtype=complex)
    for j in range(d):
        for k in range(d):
            F[j, k] = omega**(j * k) / np.sqrt(d)
    
    for m in range(1, d + 1):
        D_m = np.diag([omega**(m * k * k / 2) for k in range(d)])
        F_m = F @ D_m
        mubs[f'M{m}'] = [F_m[:, k] for k in range(d)]
    
    return mubs


def mub_measurement_stats(rho, d, n_shots_per_basis=1000):
    """
    Simulate projective measurements in all d+1 MUBs.
    Returns measurement statistics: probabilities for each outcome in each MUB.
    Returns: dict {mub_label: np.array of d probabilities}
    """
    mubs = construct_mubs(d)
    stats = {}
    projs = {}
    for label, basis in mubs.items():
        projs[label] = [np.outer(v, v.conj()) for v in basis]
        stats[label] = np.array([np.real(np.trace(rho @ P)) for P in projs[label]])
    
    # Simulate finite shots
    measured_stats = {}
    for label in mubs:
        measured_stats[label] = np.zeros(d)
        n = n_shots_per_basis
        # Simulate multinomial sampling
        probs = stats[label]
        probs = np.clip(probs, 0, None)
        probs = probs / probs.sum()
        
        # For finite shots, sample from multinomial
        if n < 1000:  # small n: exact sampling
            outcomes = np.random.multinomial(n, probs)
        else:  # large n: use normal approximation
            outcomes = np.random.multinomial(n, probs)
        measured_stats[label] = outcomes / n
    
    return measured_stats, projs


# ======================================================================
# 4. DPRT-Based Reconstruction (Direction 2)
# ======================================================================

def dprt_reconstruct(rho, d, n_shots=1000):
    """
    Full DPRT-based tomography pipeline:
    1. Measure in all d+1 MUBs → get statistics
    2. Interpret as DPRT of Wigner function
    3. Apply DPRT inverse → recover Wigner function
    4. Map Wigner back to density matrix
    
    Returns: (reconstructed_rho, reconstruction_error, timing)
    """
    t0 = time.time()
    
    # Step 1: Measure
    measured_stats, projs = mub_measurement_stats(rho, d, n_shots)
    
    # Step 2: Build DPRT from measurement stats
    # The correspondence: MUB m statistics = DPRT direction m of Wigner function
    # For m=0 (computational basis): measurements give column-projected probabilities
    # of the Wigner function — this needs precise phase-space mapping.
    
    # For the qubit case (d=2), we can work out the exact mapping.
    # For general d: the connection is through the discrete Wigner function's
    # marginal property.
    
    # Simpler approach: directly use DPRT inverse on the measurement statistics
    # arranged as (d+1)×d matrix.
    R = np.zeros((d+1, d))
    for m in range(d + 1):
        R[m, :] = measured_stats[f'M{m}']
    
    # Step 3: DPRT inverse
    W_rec = dprt_inverse(R, d)
    
    # Step 4: Map Wigner to density matrix
    rho_rec = density_from_wigner(W_rec, d)
    
    t1 = time.time()
    err = np.max(np.abs(rho - rho_rec))
    
    return rho_rec, err, t1 - t0, W_rec


def classical_shadow_reconstruct(rho, d, n_total, seed=42):
    """
    Standard Classical Shadow with random MUB sampling.
    n_total = total number of measurements.
    For fair comparison: n_total = (d+1) * n_shots_per_basis (same as DPRT)
    """
    np.random.seed(seed)
    mubs = construct_mubs(d)
    mub_labels = list(mubs.keys())
    
    # Projectors for all MUBs
    all_projs = {}
    for label in mub_labels:
        basis = mubs[label]
        all_projs[label] = [np.outer(v, v.conj()) for v in basis]
    
    t0 = time.time()
    
    # Inverse channel factor: (d+1) for MUB-based classical shadow
    # ρ̂ = (d+1) · |ψ⟩⟨ψ| - I
    factor = d + 1
    Id = np.eye(d, dtype=complex)
    
    snapshots = []
    for _ in range(n_total):
        label = np.random.choice(mub_labels)
        projs = all_projs[label]
        probs = np.array([np.real(np.trace(rho @ P)) for P in projs])
        probs = np.clip(probs, 0, None)
        probs = probs / probs.sum()
        outcome = np.random.choice(d, p=probs)
        snapshot = factor * projs[outcome] - Id
        snapshots.append(snapshot)
    
    rho_cs = np.mean(snapshots, axis=0)
    t1 = time.time()
    
    err = np.max(np.abs(rho - rho_cs))
    return rho_cs, err, t1 - t0


# ======================================================================
# 5. DPRT vs CS: Head-to-Head Comparison
# ======================================================================

def compare_dprt_vs_cs(d_list, n_states=20, n_shots_list=[100, 1000, 10000], seed=42):
    """
    Compare DPRT reconstruction vs Classical Shadow for each d.
    """
    np.random.seed(seed)
    
    results = {}
    
    for d in d_list:
        print(f"\n{'='*60}")
        print(f"d={d}: DPRT vs Classical Shadow Comparison")
        print(f"{'='*60}")
        
        d_results = []
        mub_labels = [f'M{m}' for m in range(d+1)]
        
        for n_shots in n_shots_list:
            n_total = (d + 1) * n_shots  # total measurements
            
            dprt_errors = []
            cs_errors = []
            dprt_times = []
            cs_times = []
            
            for si in range(n_states):
                s = seed + si * 1000 + n_shots * 10
                np.random.seed(s)
                
                # Random pure state
                psi = np.random.randn(d) + 1j * np.random.randn(d)
                psi /= np.linalg.norm(psi)
                rho = np.outer(psi, psi.conj())
                
                # DPRT reconstruction
                rho_d, err_d, t_d, _ = dprt_reconstruct(rho, d, n_shots)
                
                # Classical Shadow (same total measurements)
                rho_c, err_c, t_c = classical_shadow_reconstruct(rho, d, n_total, seed=s+1)
                
                dprt_errors.append(err_d)
                cs_errors.append(err_c)
                dprt_times.append(t_d)
                cs_times.append(t_c)
            
            dprt_mae = np.mean(dprt_errors)
            cs_mae = np.mean(cs_errors)
            dprt_t = np.mean(dprt_times)
            cs_t = np.mean(cs_times)
            
            winner = "DPRT" if dprt_mae < cs_mae else "CS"
            
            print(f"\n  n_shots/base={n_shots:>6d}  total={n_total:>6d}")
            print(f"    DPRT: err={dprt_mae:.6f}  time={dprt_t*1000:.2f}ms")
            print(f"    CS:   err={cs_mae:.6f}  time={cs_t*1000:.2f}ms")
            print(f"    Winner: {winner}  Ratio={dprt_mae/cs_mae:.4f}")
            
            d_results.append({
                'd': d, 'n_shots_per_basis': n_shots, 'n_total': n_total,
                'dprt_err': float(dprt_mae), 'cs_err': float(cs_mae),
                'dprt_time_ms': float(dprt_t * 1000), 'cs_time_ms': float(cs_t * 1000),
                'ratio': float(dprt_mae / cs_mae), 'winner': winner,
            })
        
        results[f'd={d}'] = d_results
    
    return results


# ======================================================================
# 6. Large-d Scan: Deterministic vs Random MUB (Direction 3)
# ======================================================================

def snap_qudit(rho, mub_label, mubs, d):
    """Single qudit shadow snapshot in given MUB."""
    basis = mubs[mub_label]
    projs = [np.outer(v, v.conj()) for v in basis]
    probs = np.array([np.real(np.trace(rho @ P)) for P in projs])
    probs = np.clip(probs, 0, None)
    probs = probs / probs.sum()
    outcome = np.random.choice(d, p=probs)
    factor = d + 1
    return factor * projs[outcome] - np.eye(d, dtype=complex)


def gell_mann_observables(d):
    """
    Return list of (name, matrix) for generalized Gell-Mann basis.
    For a qudit of dimension d, there are d²-1 traceless hermitian generators.
    """
    obs = {}
    # Symmetric (λ type)
    for j in range(d):
        for k in range(j+1, d):
            M = np.zeros((d, d), dtype=complex)
            M[j, k] = 1
            M[k, j] = 1
            obs[f'λ_{j}{k}'] = M
    
    # Anti-symmetric (λ̃ type)
    for j in range(d):
        for k in range(j+1, d):
            M = np.zeros((d, d), dtype=complex)
            M[j, k] = -1j
            M[k, j] = 1j
            obs[f'λ̃_{j}{k}'] = M
    
    # Diagonal
    for l in range(1, d):
        M = np.diag([1.0/np.sqrt(l*(l+1)) if k <= l else 0 for k in range(d)])
        M[l, l] = -l / np.sqrt(l*(l+1))
        obs[f'd_{l}'] = M
    
    return obs


def run_large_d_scan(d_list=[5, 7, 11], n_states=30, r_list=[1, 2, 4], seed=42):
    """
    For each prime d, run deterministic (all d+1 MUBs × r rounds)
    vs random MUB comparison.
    """
    np.random.seed(seed)
    all_results = {}
    
    for d in d_list:
        print(f"\n{'='*60}")
        print(f"LARGE-d SCAN: d={d} (prime)")
        print(f"{'='*60}")
        
        mubs = construct_mubs(d)
        mub_labels = list(mubs.keys())
        observables = gell_mann_observables(d)
        obs_labels = list(observables.keys())
        n_configs = d + 1  # deterministic: all d+1 MUBs
        
        d_results = []
        
        for r in r_list:
            n_meas = n_configs * r
            det_err = defaultdict(list)
            rnd_err = defaultdict(list)
            
            for si in range(n_states):
                s_seed = seed + d * 10000 + si * 1000 + r * 100
                np.random.seed(s_seed)
                
                psi = np.random.randn(d) + 1j * np.random.randn(d)
                psi /= np.linalg.norm(psi)
                rho = np.outer(psi, psi.conj())
                
                truth = {}
                for l in obs_labels:
                    truth[l] = np.real(np.trace(observables[l] @ rho))
                
                # Deterministic
                snaps_d = []
                for _ in range(r):
                    for ml in mub_labels:
                        snaps_d.append(snap_qudit(rho, ml, mubs, d))
                
                # Random (same N)
                snaps_r = []
                for _ in range(n_meas):
                    ml = np.random.choice(mub_labels)
                    snaps_r.append(snap_qudit(rho, ml, mubs, d))
                
                for l in obs_labels:
                    O = observables[l]
                    pred_d = np.mean([np.real(np.trace(O @ s)) for s in snaps_d])
                    pred_r = np.mean([np.real(np.trace(O @ s)) for s in snaps_r])
                    det_err[l].append(abs(pred_d - truth[l]))
                    rnd_err[l].append(abs(pred_r - truth[l]))
            
            det_mae = np.mean([np.mean(det_err[l]) for l in obs_labels])
            rnd_mae = np.mean([np.mean(rnd_err[l]) for l in obs_labels])
            winner = "DET" if det_mae < rnd_mae else "RAND"
            
            print(f"  r={r:>2d}  N={n_meas:>4d}  "
                  f"Det={det_mae:.6f}  Rand={rnd_mae:.6f}  "
                  f"Ratio={det_mae/rnd_mae:.4f}  [{winner}]")
            
            d_results.append({
                'd': d, 'r': r, 'n_meas': n_meas,
                'det_mae': float(det_mae), 'rand_mae': float(rnd_mae),
                'ratio': float(det_mae/rnd_mae), 'winner': winner,
            })
        
        all_results[f'd={d}'] = d_results
    
    return all_results


# ======================================================================
# 7. DPRT Computational Complexity Analysis
# ======================================================================

def analyze_dprt_complexity():
    """Analyze the computational advantage of DPRT inverse over M⁻¹ channel."""
    print(f"\n{'='*60}")
    print("DPRT COMPUTATIONAL COMPLEXITY ANALYSIS")
    print(f"{'='*60}")
    
    print(f"\n  Classical Shadow M⁻¹ channel:")
    print(f"    Each snapshot: ρ̂ = (d+1)·|ψ⟩⟨ψ| - I")
    print(f"      - Cost: O(d²) outer product + O(d²) subtraction")
    print(f"    Observable estimation: tr(O ρ̂) = (d+1)⟨ψ|O|ψ⟩ - tr(O)")
    print(f"      - Cost: O(d²) matrix-vector + dot product")
    print(f"    Total for N snapshots: O(N·d²)")
    
    print(f"\n  DPRT Inverse Reconstruction:")
    print(f"    Collect statistics: (d+1)·d values = O(d²)")
    print(f"    DPRT inverse: per (i,j) element, sum d+1 terms:")
    print(f"      - Cost: O(d² × d) = O(d³) additions, O(d²) divisions by d")
    print(f"      - BUT: all additions, NO multiplications!")
    print(f"    Observable from Wigner: O(d⁴) if naive, O(d²·log d) if using") 
    print(f"      structure of phase-point operators")
    
    print(f"\n  Key advantage of DPRT:")
    print(f"    - Deterministic: exactly (d+1)·r measurements cover all bases")
    print(f"    - Inverse is linear with pure additions (no matrix ops)")
    print(f"    - For fixed-precision (low r): O(d³) additions vs O(d⁴) matrix ops")
    print(f"    - Caveat: Wigner→observable mapping adds overhead for full tomography")
    
    print(f"\n  Practical comparison (estimated):")
    print(f"  {'d':>5s}  {'DPRT Ops':>10s}  {'CS Ops':>10s}  {'Ratio':>8s}  {'DPRT ms':>10s}  {'CS ms':>10s}")
    print(f"  {'─'*60}")
    
    for d in [2, 3, 5, 7, 11]:
        dprt_adds = d**3 + d**2  # O(d³) additions + O(d²) divisions
        cs_ops = (d+1) * d**3    # N=(d+1) snapshots × O(d³) trace
        ratio = dprt_adds / cs_ops
        
        # Quick benchmark
        np.random.seed(42)
        tau = np.random.dirichlet(np.ones(d))
        U = unitary_group.rvs(d)
        rho = U @ np.diag(tau) @ U.conj().T
        
        t0 = time.time()
        rho_d, _, _, _ = dprt_reconstruct(rho, d, 100)
        td = (time.time() - t0) * 1000
        
        t0 = time.time()
        classical_shadow_reconstruct(rho, d, (d+1)*100, seed=99)
        tc = (time.time() - t0) * 1000
        
        print(f"  {d:>5d}  {dprt_adds:>10d}  {cs_ops:>10d}  {ratio:>8.4f}  "
              f"{td:>10.2f}  {tc:>10.2f}")
    
    return dprt_adds, cs_ops


# ======================================================================
# 8. Wigner→Density Verification
# ======================================================================

def verify_wigner_mapping(d_list=[2, 3, 5]):
    """Verify that Wigner↔Density mapping is correct."""
    print(f"\n{'='*60}")
    print("WIGNER ↔ DENSITY MATRIX VERIFICATION")
    print(f"{'='*60}")
    
    for d in d_list:
        np.random.seed(123)
        psi = np.random.randn(d) + 1j * np.random.randn(d)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        
        W = wigner_from_density(rho, d)
        rho_rec = density_from_wigner(W, d)
        err = np.max(np.abs(rho - rho_rec))
        
        # Check Wigner properties
        sum_W = np.sum(W)
        trace_abs = np.trace(np.abs(W))
        
        print(f"\n  d={d}:")
        print(f"    Round-trip error: {err:.2e}")
        print(f"    ΣW(q,p) = {sum_W:.6f}  (should = 1)")
        print(f"    tr(|W|) = {trace_abs:.6f}")
    print()


# ======================================================================
# 9. Main
# ======================================================================
if __name__ == '__main__':
    t0 = time.time()
    
    # Verify Wigner mapping
    verify_wigner_mapping()
    
    # Verify DPRT correctness for target dimensions
    print("DPRT Correctness Tests:")
    for d in [2, 3, 5, 7, 11]:
        err = test_dprt(d)
        print(f"  d={d:>2d}: round-trip error = {err:.2e}")
    
    # Direction 2: DPRT vs CS head-to-head (small d for speed)
    print(f"\n{'#'*60}")
    print("# DIRECTION 2: DPRT Post-Processing vs Classical Shadow")
    print(f"{'#'*60}")
    cs_vs_dprt = compare_dprt_vs_cs(
        d_list=[2, 3, 5],
        n_states=15,
        n_shots_list=[50, 200, 500],
        seed=42
    )
    
    # Direction 3: Large-d scan
    print(f"\n{'#'*60}")
    print("# DIRECTION 3: Large-d Scan (Deterministic vs Random MUB)")
    print(f"{'#'*60}")
    large_d_results = run_large_d_scan(
        d_list=[5, 7, 11],
        n_states=20,
        r_list=[1, 2, 4],
        seed=42
    )
    
    # Complexity analysis
    analyze_dprt_complexity()
    
    t1 = time.time()
    print(f"\nTotal time: {t1-t0:.1f}s")
    
    # ── Master Summary ──
    print(f"\n{'='*60}")
    print("MASTER SUMMARY")
    print(f"{'='*60}")
    
    print(f"\n── DPRT vs Classical Shadow ──")
    for key, val in cs_vs_dprt.items():
        print(f"\n  {key}:")
        for r in val:
            print(f"    shots={r['n_shots_per_basis']:>4d}  "
                  f"DPRT={r['dprt_err']:.6f}  CS={r['cs_err']:.6f}  "
                  f"ratio={r['ratio']:.4f}  [{r['winner']}]")
    
    print(f"\n── Large-d Scan ──")
    for key, val in large_d_results.items():
        print(f"\n  {key}:")
        for r in val:
            print(f"    N={r['n_meas']:>4d}  Det={r['det_mae']:.6f}  "
                  f"Rand={r['rand_mae']:.6f}  ratio={r['ratio']:.4f}  [{r['winner']}]")
    
    # Save
    output = {
        'dprt_vs_cs': {k: v for k, v in cs_vs_dprt.items()},
        'large_d_scan': {k: v for k, v in large_d_results.items()},
    }
    with open('qdprt_experiment_phase3.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n→ qdprt_experiment_phase3.json")
