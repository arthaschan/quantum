#!/usr/bin/env python3
"""批量替换数学符号 — 使用简单字符串替换，避免 regex 转义问题"""
import os

REPLACEMENTS = [
    # 希腊字母（空格或标点包围的独立出现）
    ( ' ρ ',  ' $\\rho$ '),   ( ' ρ,', ' $\\rho$,'),   ( ' ρ.', ' $\\rho$.'),
    ( ' ρ;', ' $\\rho$;'),   ( ' ρ)', ' $\\rho$)'),   ( ' ρ]', ' $\\rho$]'),
    ( '(ρ ', '($\\rho$ '),   ( 'd=ρ', 'd=$\\rho$'),
    
    ( ' ψ ',  ' $\\psi$ '),   ( ' ψ,', ' $\\psi$,'),   ( ' ψ.', ' $\\psi$.'),
    ( ' ψ;', ' $\\psi$;'),   ( ' ψ)', ' $\\psi$)'),
    ( '(ψ ', '($\\psi$ '),
    
    ( ' ω ',  ' $\\omega$ '),  ( ' ω,', ' $\\omega$,'),  ( ' ω.', ' $\\omega$.'),
    ( ' ω;', ' $\\omega$;'),  ( ' ω)', ' $\\omega$)'),
    
    ( ' σ ',  ' $\\sigma$ '),  ( ' σ,', ' $\\sigma$,'),  ( ' σ.', ' $\\sigma$.'),
    
    # 上标
    ('σ²', '$\\sigma^2$'), ('ρ²', '$\\rho^2$'),
    ('O(p²', '$O(p^2$'), ('O(d²', '$O(d^2$'), ('p² ', '$p^2$ '),
    ('d² ', '$d^2$ '),
    
    # 特殊字符  
    ('d×d', '$d \\times d$'), ('p×p', '$p \\times p$'),
    
    # 求和符号（不在 $ 内的）
    ('Σ_{', '$\\sum_{'), ('Σ_i', '$\\sum_i$'),
    
    # Bra-ket（不在 $ 内的 — 用字符上下文判断）
    ('|ψ⟩', '$|\\psi\\rangle$'), ('⟨ψ|', '$\\langle\\psi|$'),
    ('|e', '$|e'), ('⟩⟨', '\\rangle\\langle'),
    
    # 常见公式
    ('ratio = σ²_DPRT / σ²_MUB', '$\\text{ratio} = \\sigma^2_{\\text{DPRT}} / \\sigma^2_{\\text{MUB}}$'),
    ('||ρ̂ − ρ||_F', '$\\|\\hat{\\rho} - \\rho\\|_F$'),
    ('tr(Oρ)', '$\\mathrm{tr}(O\\rho)$'),
    
    # ≥ 和 ≈ 等符号
    (' ≥ ', ' $\\ge$ '), (' ≥0', ' $\\ge 0$'),
    (' ≈ ', ' $\\approx$ '), (' ≠ ', ' $\\neq$ '),
    (' ∈ ', ' $\\in$ '),   (' → ', ' $\\to$ '),
    ('∞ ', '$\\infty$ '),
    
    # 分式
    (' 1/d', ' $1/d$'), (' I_d/d', ' $I_d / d$'),
    
    # 指数和排序
    ('2^k', '$2^k$'), ('O(log M)', '$O(\\log M)$'),
    ('O(log d)', '$O(\\log d)$'), ('O(d³/N)', '$O(d^3/N)$'),
    
    # 数域
    ('F_p', '$\\mathbb{F}_p$'), ('F_p*', '$\\mathbb{F}_p^*$'), 
    ('Z_p', '$\\mathbb{Z}_p$'),
    
    # 维度表示
    ('d ≤ 1000', '$d \\le 1000$'), ('d ≤ 200', '$d \\le 200$'),
    ('d ≤ 30', '$d \\le 30$'),
]

def fix_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        text = f.read()
    
    # 保护已有 $$...$$ 和 $...$ 块
    protected_blocks = []
    
    # 保护 $$...$$ 
    import re
    def save_display(m):
        protected_blocks.append(m.group(0))
        return f'<<PROTECTED_DISPLAY_{len(protected_blocks)-1}>>'
    text = re.sub(r'\$\$[^$]+\$\$', save_display, text, flags=re.DOTALL)
    
    # 保护 $...$（跳过已经是保护标记的）
    def save_inline(m):
        if '<<PROTECTED' in m.group(0):
            return m.group(0)
        # 跳过货币符号
        if re.match(r'\$\d+', m.group(0)):
            return m.group(0)
        protected_blocks.append(m.group(0))
        return f'<<PROTECTED_INLINE_{len(protected_blocks)-1}>>'
    text = re.sub(r'\$[^$]+\$', save_inline, text)
    
    changes = 0
    for old, new in REPLACEMENTS:
        count = text.count(old)
        if count > 0:
            text = text.replace(old, new)
            changes += count
    
    # 恢复保护块
    for i, block in enumerate(protected_blocks):
        text = text.replace(f'<<PROTECTED_DISPLAY_{i}>>', block)
        text = text.replace(f'<<PROTECTED_INLINE_{i}>>', block)
    
    if changes > 0:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(text)
    
    return changes


if __name__ == '__main__':
    files = [
        '/Users/arthas/git/quantum/reports/layperson_guide.md',
        '/Users/arthas/git/quantum/reports/experiment_complete_report.md',
    ]
    
    for fp in files:
        fn = os.path.basename(fp)
        print(f'Processing {fn}...')
        n = fix_file(fp)
        print(f'  {n} replacements made')
    
    print('Done.')
