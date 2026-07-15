#!/usr/bin/env python3
"""保守替换 — 只处理最简单的独立数学符号，不碰 bra-ket 和复杂公式"""
import os, re

SIMPLE_FIXES = [
    # 独立希腊字母（被空格或标点包围）
    (' ρ ', ' $\\rho$ '), (' ρ.', ' $\\rho$.'), (' ρ,', ' $\\rho$,'),
    (' ρ;', ' $\\rho$;'), (' ρ)', ' $\\rho$)'), ('(ρ)', '($\\rho$)'),
    ('ρ 的', '$\\rho$ 的'), ('将 ρ', '将 $\\rho$'), ('对 ρ', '对 $\\rho$'),
    
    # d×d, p×p（不在$内的）
    (' d×d ', ' $d \\times d$ '), (' p×p ', ' $p \\times p$ '),
    (' d×d，', ' $d \\times d$，'), (' d×d。', ' $d \\times d$。'),
    
    # p+1 组 → 保持不碰（这是中文数字语境）
    
    # Σ 求和符号
    (' Σ_i ', ' $\\sum_i$ '), (' Σ_{i,', ' $\\sum_{i,}$ '),
    (' Σ_{i=', ' $\\sum_{i='),
    
    # ≥ 和特殊符号  
    (' ≥ ', ' $\\ge$ '), (' ≥0', ' $\\ge 0$'), (' ≥1', ' $\\ge 1$'),
    (' ≈ ', ' $\\approx$ '), (' ≠ ', ' $\\neq$ '),
    (' ∈ ', ' $\\in$ '),
    ('∞', '$\\infty$'),
    
    # 常见维度
    (' d+1 组', ' $d+1$ 组'), (' d+1 个', ' $d+1$ 个'),
    
    # σ²
    (' σ²', ' $\\sigma^2$'), ('σ²_', '$\\sigma^2_'),
    
    # F_p
    (' F_p', ' $\\mathbb{F}_p$'), ('F_p*', '$\\mathbb{F}_p^*$'),
    
    # √d
    (' 1/√d', ' $1/\\sqrt{d}$'), (' 1/√p', ' $1/\\sqrt{p}$'),
    
    # O(...)
    (' O(d²', ' $O(d^2$'), (' O(p²', ' $O(p^2$'),
    (' O(d³', ' $O(d^3$'), (' O(log', ' $O(\\log'),
]

def fix_file_safe(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        text = f.read()
    
    original = text
    
    # 保护 $$...$$ 
    display_blocks = []
    def save_display(m):
        display_blocks.append(m.group(0))
        return f'<<<DD{len(display_blocks)-1}>>>'
    text = re.sub(r'\$\$[^$]+\$\$', save_display, text, flags=re.DOTALL)
    
    # 保护 $...$（不与$$重叠，非货币）
    inline_blocks = []
    def save_inline(m):
        if '<<<DD' in m.group(0):
            return m.group(0)
        inline_blocks.append(m.group(0))
        return f'<<<II{len(inline_blocks)-1}>>>'
    text = re.sub(r'(?<!\$)\$[^$\n]+?\$(?!\$)', save_inline, text)
    
    changes = 0
    for old, new in SIMPLE_FIXES:
        if old in text:
            count = text.count(old)
            text = text.replace(old, new)
            changes += count
    
    # 恢复
    for i, block in enumerate(display_blocks):
        text = text.replace(f'<<<DD{i}>>>', block)
    for i, block in enumerate(inline_blocks):
        text = text.replace(f'<<<II{i}>>>', block)
    
    if text != original:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(text)
    
    return changes

files = [
    '/Users/arthas/git/quantum/reports/layperson_guide.md',
    '/Users/arthas/git/quantum/reports/experiment_complete_report.md', 
]
for fp in files:
    fn = os.path.basename(fp)
    print(f'{fn}: {fix_file_safe(fp)} replacements')
print('Done.')
