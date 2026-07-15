#!/usr/bin/env python3
"""
精确修复 Markdown 数学符号 — 保守方法
只替换明确需要修复的模式，不修改已有 $$...$$ 和 $...$ 块
"""
import re, os

def fix_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    new_lines = []
    changes = 0
    
    for i, line in enumerate(lines):
        original = line
        
        # 跳过代码块
        if line.strip().startswith('```') or line.strip().startswith('~~~'):
            new_lines.append(line)
            continue
        
        # 跳过已经是纯数学公式的行（$$...$$）
        if line.strip().startswith('$$') and line.strip().endswith('$$'):
            new_lines.append(line)
            continue
        
        # 只处理非代码、非纯公式的行
        # 用占位符保护已有的 $...$ 块
        protected = []
        def save_math(m):
            protected.append(m.group(0))
            return f'\x00MATH{len(protected)-1}\x00'
        line = re.sub(r'\$[^$]+\$', save_math, line)
        
        # ── 简单替换（独立出现的 Unicode 字符）──
        # 只在不在单词内部时替换
        
        # Greek letters (standalone, not inside words)
        line = re.sub(r'(?<=[ (［【〈"\'“])ρ(?=[ ,;.）］】〉"\'”\n])', r'$\rho$', line)
        line = re.sub(r'(?<=[ (［【〈"\'“])ψ(?=[ ,;.）］】〉"\'”\n])', r'$\psi$', line)
        line = re.sub(r'(?<=[ (［【〈"\'“])ω(?=[ ,;.）］】〉"\'”\n])', r'$\omega$', line)
        line = re.sub(r'(?<=[ (［【〈"\'“])σ(?=[ ,;.）］】〉"\'”\n²³])', r'$\sigma$', line)
        line = re.sub(r'(?<=[ (［【〈"\'“])φ(?=[ ,;.）］】〉"\'”\n])', r'$\varphi$', line)
        line = re.sub(r'(?<=[ (［【〈"\'“])λ(?=[ ,;.）］】〉"\'”\n])', r'$\lambda$', line)
        
        # Also handle ρ at end of sentence or followed by ²
        line = re.sub(r'ρ²', r'$\rho^2$', line)
        line = re.sub(r'σ²', r'$\sigma^2$', line)
        
        # × → $\times$ (but not in tables or pure text)
        line = re.sub(r'(\d)\s*×\s*(\d)', r'$\1 \\times \2$', line)
        
        # Restore protected math
        for j, ph in enumerate(protected):
            line = line.replace(f'\x00MATH{j}\x00', ph)
        
        if line != original:
            changes += 1
        
        new_lines.append(line)
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(new_lines)
    
    return changes

files = [
    '/Users/arthas/git/quantum/reports/layperson_guide.md',
    '/Users/arthas/git/quantum/reports/experiment_complete_report.md',
]

for fp in files:
    print(f'Processing {os.path.basename(fp)}...')
    c = fix_file(fp)
    print(f'  Changed {c} lines')
