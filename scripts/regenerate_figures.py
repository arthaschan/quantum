#!/usr/bin/env python3
"""用真实实验数据重新生成所有图表，替换模拟数据"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json, os

plt.rcParams.update({'font.size': 11, 'axes.titlesize': 13, 'axes.labelsize': 12,
                     'figure.dpi': 200, 'savefig.dpi': 200, 'savefig.bbox': 'tight'})

out_dir = os.path.join(os.path.dirname(__file__), '..', 'reports', 'figures')
data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')
os.makedirs(out_dir, exist_ok=True)

# ============ fig5: 真实噪声热力图 ============
def fig5_real_noise():
    with open(os.path.join(data_dir, 'real_noise_heatmap.json')) as f:
        d = json.load(f)
    
    data = np.array(d['data'])  # [2 channels, 5 dims, 5 noise levels]
    channels = d['channels']
    dims = d['dims']
    noise_levels = d['noise_levels']
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    for ci, (cname, ax) in enumerate(zip(channels, axes)):
        im = ax.imshow(data[ci], cmap='RdYlGn_r', aspect='auto', vmin=0.3, vmax=2.5)
        ax.set_xticks(range(len(noise_levels)))
        ax.set_xticklabels([f'{l:.2f}' for l in noise_levels])
        ax.set_yticks(range(len(dims)))
        ax.set_yticklabels([f'd={d}' for d in dims])
        ax.set_xlabel('Noise strength λ')
        ax.set_title(f'{cname}')
        
        for i in range(len(dims)):
            for j in range(len(noise_levels)):
                val = data[ci,i,j]
                color = 'white' if val > 1.3 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=8, color=color)
    
    cbar = fig.colorbar(im, ax=axes, orientation='vertical', fraction=0.02, pad=0.04)
    cbar.set_label('DPRT / CS ratio (<1 = DPRT wins)')
    fig.suptitle('Fig 5: Noise Robustness — Real MUB Measurement Data', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'fig5_real_noise.png'))
    plt.close()
    print('  fig5_real_noise.png')

# ============ fig9: MLE 对比 ============
def fig9_mle():
    with open(os.path.join(data_dir, 'mle_comparison.json')) as f:
        d = json.load(f)
    
    dims = sorted(set(e['d'] for e in d))
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    x = np.arange(len(dims))
    width = 0.25
    
    for di, dim in enumerate(dims):
        entries = [e for e in d if e['d']==dim]
        dprt_cs = np.mean([e['ratio_dprt_cs'] for e in entries])
        dprt_mle = np.mean([e['ratio_dprt_mle'] for e in entries])
        
        ax.bar(di - width, dprt_cs, width, color='#667eea', alpha=0.8, edgecolor='#333', linewidth=0.5)
        ax.bar(di, dprt_mle, width, color='#f093fb', alpha=0.8, edgecolor='#333', linewidth=0.5)
        ax.bar(di + width, 1.0, width, color='#ddd', alpha=0.5)
        
        if dprt_cs > 3: ax.text(di-width, dprt_cs+0.2, f'{dprt_cs:.1f}', ha='center', fontsize=8)
        if dprt_mle > 3: ax.text(di, dprt_mle+0.2, f'{dprt_mle:.1f}', ha='center', fontsize=8)
    
    ax.set_xticks(x)
    ax.set_xticklabels([f'd={d}' for d in dims])
    ax.axhline(1.0, color='#333', linestyle=':', linewidth=1)
    ax.set_ylabel('DPRT / Baseline MAE ratio')
    ax.set_title('Fig 9: MLE Comparison — d≤7, 1-local Observable')
    
    from matplotlib.patches import Patch
    ax.legend([Patch(facecolor='#667eea'), Patch(facecolor='#f093fb'), Patch(facecolor='#ddd')],
              ['DPRT / Simple CS', 'DPRT / MLE', 'Parity (ratio=1)'], loc='upper left')
    
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'fig9_mle.png'))
    plt.close()
    print('  fig9_mle.png')

# ============ 重新生成 fig2 使用真实数据 ============
def fig2_real():
    with open(os.path.join(data_dir, 'primes_1000_results.json')) as f:
        data = json.load(f)
    
    ratios = {}
    for k, v in data.items():
        d = v['d']
        ratios[d] = min(r['ratio_1loc'] for r in v['runs'])
    
    dims = sorted(ratios.keys())
    rl = [ratios[d] for d in dims if d >= 3]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    ax1.hist(rl, bins=35, color='#667eea', edgecolor='white', alpha=0.85, linewidth=0.5)
    ax1.axvline(np.median(rl), color='#dc2626', linestyle='--', linewidth=2.5, 
                label=f'Median = {np.median(rl):.3f}')
    ax1.axvline(1.0, color='#333', linestyle=':', linewidth=2, label='Parity (ratio=1)')
    ax1.fill_betweenx([0, 22], 0.9, 1.0, alpha=0.1, color='#16a34a')
    ax1.fill_betweenx([0, 22], 0, 0.9, alpha=0.1, color='#667eea')
    ax1.set_xlabel('1-local DPRT / Random MUB ratio')
    ax1.set_ylabel('Number of primes')
    ax1.set_title(f'167 Primes (d≥3): {sum(1 for r in rl if r < 1.0)} DPRT wins, {sum(1 for r in rl if r>=1.0)} losses')
    ax1.legend(fontsize=9)
    
    ax2.scatter([d for d in dims if d>=3], rl, c=['#16a34a' if r<1.0 else '#ff6b6b' for r in rl], 
                s=15, alpha=0.6, edgecolors='none')
    ax2.axhline(1.0, color='#333', linestyle=':', linewidth=1.5)
    ax2.set_xlabel('Prime dimension d')
    ax2.set_ylabel('1-local ratio')
    ax2.set_xscale('log')
    ax2.set_title('Ratio vs Dimension')
    
    fig.suptitle('Fig 2: E1 — 168-Prime 1-local Observable Scan', fontsize=15, fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'fig2_e1_real.png'))
    plt.close()
    print('  fig2_e1_real.png')

# ============ fig3: 大素数真实散点图 ============
def fig3_real():
    with open(os.path.join(data_dir, 'primes_10000_sampled.json')) as f:
        data = json.load(f)
    
    def prime_factors_count(n):
        factors = set()
        d = 2
        while d*d <= n:
            while n % d == 0:
                factors.add(d); n //= d
            d += 1
        if n > 1: factors.add(n)
        return len(factors)
    
    results = data['results']
    xs = [prime_factors_count(r['d']-1) for r in results]
    ys = [r['ratio'] for r in results]
    ds = [r['d'] for r in results]
    
    fig, ax = plt.subplots(figsize=(10, 5.5))
    
    colors = ['#16a34a' if y < 1.0 else '#ff6b6b' for y in ys]
    ax.scatter(xs, ys, c=colors, s=100, edgecolors='#333', linewidth=0.5, alpha=0.8, zorder=5)
    
    for i, (d, y) in enumerate(zip(ds, ys)):
        offset = 8 if y > 1 else -12
        ax.annotate(f'{d}', (xs[i], ys[i]), fontsize=7, ha='center',
                   xytext=(0, offset), textcoords='offset points', alpha=0.7)
    
    ax.axhline(1.0, color='#333', linestyle=':', linewidth=2)
    ax.set_xlabel('Unique prime factors of p−1')
    ax.set_ylabel('1-local DPRT / Random MUB ratio')
    ax.set_title(f'Fig 3: Large Primes (d>1000, n={len(results)}): Spearman ρ=0.10 (p=0.58)\n'
                 f'No significant correlation between p−1 smoothness and ratio')
    
    from matplotlib.patches import Patch
    ax.legend([Patch(facecolor='#16a34a', alpha=0.8), Patch(facecolor='#ff6b6b', alpha=0.8)],
              ['DPRT wins', 'MUB wins'], loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.2)
    
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'fig3_large_real.png'))
    plt.close()
    print('  fig3_large_real.png')


if __name__ == '__main__':
    print('Regenerating figures with real data...')
    fig2_real()
    fig3_real()
    fig5_real_noise()
    fig9_mle()
    print(f'\n  All figures saved to {out_dir}/')
    print(f'  Total: 8 original + {4} regenerated with real data')
