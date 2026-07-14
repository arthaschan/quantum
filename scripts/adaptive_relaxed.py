#!/usr/bin/env python3
"""
RadonShadow: Adaptive Relaxed — 放宽 ε 阈值扫描
=================================================
确认自适应提前停止实际能省多少测量次数。

ε ∈ {0.05, 0.10, 0.20, 0.50}
d ∈ {3, 5, 7, 11, 13, 17, 19, 23}
纯态 + depol(λ=0, 0.05)
"""

import numpy as np, time, json

def build_mubs(d):
    omega = np.exp(2j * np.pi / d)
    mubs = []
    mubs.append([np.eye(d, dtype=complex)[j] for j in range(d)])
    for m in range(1, d + 1):
        basis = []
        for k in range(d):
            v = np.array([omega**(m*j*(j-1)//2 + k*j) for j in range(d)], dtype=complex)
            basis.append(v / np.sqrt(d))
        mubs.append(basis)
    return mubs

def apply_noise(rho, channel, lam):
    d = rho.shape[0]
    if channel == 'depol':
        return (1 - lam) * rho + lam * np.eye(d) / d
    elif channel == 'amp':
        gamma = lam
        K0 = np.diag([1.0] + [np.sqrt(1 - gamma)] * (d - 1))
        result = K0 @ rho @ K0.conj().T
        for j in range(1, d):
            K = np.zeros((d, d), dtype=complex)
            K[0, j] = np.sqrt(gamma / (d - 1))
            result += K @ rho @ K.conj().T
        return result
    elif channel == 'phase':
        result = (1 - lam) * rho.copy()
        for i in range(d):
            result[i,i] = rho[i,i].real
        return result
    elif channel == 'bitflip':
        Xd = np.zeros((d,d), dtype=complex)
        for j in range(d):
            Xd[(j+1)%d, j] = 1.0
        return (1-lam)*rho + lam*(Xd@rho@Xd.conj().T)

def get_1local(d):
    obs = []
    for i in range(d):
        for j in range(i+1, d):
            O = np.zeros((d,d), dtype=complex)
            O[i,j] = O[j,i] = 1.0
            obs.append(O)
    return obs[:min(8, len(obs))]

def shadow_invert(total, N):
    d = total.shape[0]
    return total * (d+1) / N - np.eye(d)

TEST_DIMS = [3, 5, 7, 11, 13, 17, 19, 23]
EPSILONS_RELAXED = [0.05, 0.10, 0.20, 0.50]
N_STATES = 8
R_VAL = 4

results = []

for d in TEST_DIMS:
    mubs = build_mubs(d)
    nm = len(mubs)
    obs = get_1local(d)
    print(f"\n  d={d}  MUBs={nm}")

    for noise_ch in [('depol', 0.0), ('depol', 0.05), ('amp', 0.05), ('phase', 0.05)]:
        ch, lam = noise_ch
        for eps in EPSILONS_RELAXED:
            savings = []
            ratios = []
            frobs = []

            for si in range(N_STATES):
                seed = d*500000 + hash(ch)%7000 + int(lam*1000) + int(eps*10000) + si
                rng = np.random.RandomState(seed)

                z = rng.randn(d) + 1j*rng.randn(d)
                psi = z / np.linalg.norm(z)
                rho_pure = np.outer(psi, psi.conj())
                rho_true = apply_noise(rho_pure, ch, lam)

                truth = [np.real(np.trace(O@rho_true)) for O in obs]

                # 全遍历
                full_total = np.zeros((d,d), dtype=complex)
                for mi in range(nm):
                    probs = np.array([np.real(np.vdot(v, rho_true@v)) for v in mubs[mi]])
                    probs = np.clip(probs,0,None); probs /= probs.sum()
                    rng_f = np.random.RandomState(seed+3000+mi)
                    outcomes = rng_f.choice(d, R_VAL, p=probs)
                    for k in outcomes:
                        full_total += np.outer(mubs[mi][k], mubs[mi][k].conj())
                rho_full = shadow_invert(full_total, nm*R_VAL)

                # 随机 baseline
                rand_total = np.zeros((d,d), dtype=complex)
                rng_r = np.random.RandomState(seed+4000)
                for _ in range(nm*R_VAL):
                    mi = rng_r.randint(nm)
                    probs = np.array([np.real(np.vdot(v, rho_true@v)) for v in mubs[mi]])
                    probs = np.clip(probs,0,None); probs /= probs.sum()
                    k = rng_r.choice(d, p=probs)
                    rand_total += np.outer(mubs[mi][k], mubs[mi][k].conj())
                rho_rand = shadow_invert(rand_total, nm*R_VAL)

                # 自适应
                adapt_total = np.zeros((d,d), dtype=complex)
                rho_prev = np.zeros((d,d), dtype=complex)
                k_stop = nm

                for mi in range(nm):
                    probs = np.array([np.real(np.vdot(v, rho_true@v)) for v in mubs[mi]])
                    probs = np.clip(probs,0,None); probs /= probs.sum()
                    rng_a = np.random.RandomState(seed+5000+mi)
                    outcomes = rng_a.choice(d, R_VAL, p=probs)
                    for k in outcomes:
                        adapt_total += np.outer(mubs[mi][k], mubs[mi][k].conj())
                    n_used = (mi+1)*R_VAL
                    rho_curr = shadow_invert(adapt_total, n_used)
                    if mi >= 1:
                        if np.linalg.norm(rho_curr - rho_prev, 'fro') < eps:
                            k_stop = mi + 1
                            break
                    rho_prev = rho_curr.copy()

                savings.append(1.0 - k_stop/nm)
                mae_a = np.mean([abs(np.real(np.trace(O@rho_curr))-t) for O,t in zip(obs, truth)])
                mae_r = np.mean([abs(np.real(np.trace(O@rho_rand))-t) for O,t in zip(obs, truth)])
                if mae_r > 1e-12:
                    ratios.append(mae_a/mae_r)
                frobs.append(float(np.linalg.norm(rho_curr - rho_full, 'fro')))

            mean_save = float(np.mean(savings))*100
            mean_ratio = float(np.mean(ratios)) if ratios else 1.0
            mean_frob = float(np.mean(frobs))
            triggered = sum(1 for s in savings if s > 0.001)
            marker = " ★★★" if mean_save > 30 else (" ★★" if mean_save > 10 else (" ★" if mean_save > 0 else ""))
            print(f"    {ch}(λ={lam:.2f}) ε={eps:.2f}  "
                  f"save={mean_save:.1f}%  "
                  f"trig={triggered}/{N_STATES}  "
                  f"ratio={mean_ratio:.4f}  "
                  f"frob={mean_frob:.4f}{marker}")

            results.append({
                'd': d, 'channel': ch, 'noise': lam, 'epsilon': eps,
                'mean_savings_pct': mean_save, 'trigger_count': triggered,
                'mean_ratio': mean_ratio, 'mean_frob_gap': mean_frob,
                'n_states': N_STATES
            })

# Save
with open('data/adaptive_relaxed_results.json', 'w') as f:
    json.dump(results, f, indent=2)

# Summary
print("\n" + "="*70)
print("SUMMARY: RELAXED ADAPTIVE")
print("="*70)

for eps in EPSILONS_RELAXED:
    pts = [r for r in results if abs(r['epsilon']-eps) < 1e-9]
    if pts:
        saves = [p['mean_savings_pct'] for p in pts]
        trigs = [p['trigger_count'] for p in pts]
        ratios = [p['mean_ratio'] for p in pts]
        triggered_pct = sum(trigs)/len(trigs)/N_STATES*100
        print(f"  ε={eps:.2f}: save={np.mean(saves):.1f}%  "
              f"triggered={triggered_pct:.0f}%  "
              f"ratio={np.mean(ratios):.4f}  "
              f"best_save={max(saves):.1f}%  "
              f"best_ratio={min(ratios):.4f}")

# Best combos
print("\n── Top 6 Combinations ──")
top = sorted(results, key=lambda x: (x['mean_savings_pct'], -x['mean_ratio']), reverse=True)[:10]
for r in top:
    print(f"  d={r['d']:<3d}  {r['channel']}(λ={r['noise']:.2f})  "
          f"ε={r['epsilon']:.2f}  save={r['mean_savings_pct']:.1f}%  "
          f"ratio={r['mean_ratio']:.4f}  trig={r['trigger_count']}/{N_STATES}")

print("\n✅ Done. See data/adaptive_relaxed_results.json")
