#!/usr/bin/env python3
"""
批量修复 Markdown 文档中的数学符号。
将 Unicode 数学符号替换为 LaTeX $...$ 格式。
保护已有 $$...$$ 和 $...$ 块不被破坏。
"""

import re
import sys

# ============================================================
# 替换规则（按顺序应用）
# ============================================================

# 规则格式: (regex_pattern, replacement, description)
# 注意：这些替换只应应用于 NOT inside existing math blocks

RULES = [
    # --- Greek letters (lowercase) ---
    # These are the most common offenders
    (r'(?<!\$)(?<!\\)\bρ\b(?![\$}])', r'$\rho$', 'Unicode ρ → $\\rho$'),
    (r'(?<!\$)(?<!\\)\bψ\b(?![\$}])', r'$\psi$', 'Unicode ψ → $\\psi$'),
    (r'(?<!\$)(?<!\\)\bω\b(?![\$}])', r'$\omega$', 'Unicode ω → $\\omega$'),
    (r'(?<!\$)(?<!\\)\bσ\b(?![\$}])', r'$\sigma$', 'Unicode σ → $\\sigma$'),
    (r'(?<!\$)(?<!\\)\bφ\b(?![\$}])', r'$\varphi$', 'Unicode φ → $\\varphi$'),
    (r'(?<!\$)(?<!\\)\bλ\b(?![\$}])', r'$\lambda$', 'Unicode λ → $\\lambda$'),
    (r'(?<!\$)(?<!\\)\bε\b(?![\$}])', r'$\varepsilon$', 'Unicode ε → $\\varepsilon$'),
    (r'(?<!\$)(?<!\\)\bχ\b(?![\$}])', r'$\chi$', 'Unicode χ → $\\chi$'),
    (r'(?<!\$)(?<!\\)\bτ\b(?![\$}])', r'$\tau$', 'Unicode τ → $\\tau$'),
    (r'(?<!\$)(?<!\\)\bπ\b(?![\$}])', r'$\pi$', 'Unicode π → $\\pi$'),
    (r'(?<!\$)(?<!\\)\bθ\b(?![\$}])', r'$\theta$', 'Unicode θ → $\\theta$'),
    (r'(?<!\$)(?<!\\)\bδ\b(?![\$}])', r'$\delta$', 'Unicode δ → $\\delta$'),
    
    # --- Bra-ket notation ---
    # ⟨x| → $\langle x|$ (but careful with existing math)
    (r'(?<!\$)(?<!\\)⟨([^⟩]*)⟩', lambda m: f'$\\langle {m.group(1)}\\rangle$', '⟨...⟩ → $\\langle...\\rangle$'),
    
    # |ψ⟩⟨ψ| patterns
    (r'(?<!\$)(?<!\\)\|([^|⟩]+)⟩', lambda m: f'$|{m.group(1)}\\rangle$', '|...⟩ → $|...\\rangle$'),
    
    # ⟨ψ| patterns (standalone bra)
    (r'(?<!\$)(?<!\\)⟨([^|]+)\|(?![\$}])', lambda m: f'$\\langle {m.group(1)}|$', '⟨...| → $\\langle...|$'),
    
    # --- Mathematical operators ---
    # Σ with subscript → $\sum_{...}$
    (r'(?<!\$)(?<!\\)Σ_\{([^}]+)\}', lambda m: f'$\\sum_{{{m.group(1)}}}$', 'Σ_{...} → $\\sum_{...}$'),
    (r'(?<!\$)(?<!\\)Σ_(\w+)', lambda m: f'$\\sum_{{{m.group(1)}}}$', 'Σ_x → $\\sum_x$'),
    
    # Π with subscript → $\prod_{...}$
    (r'(?<!\$)(?<!\\)Π_\{([^}]+)\}', lambda m: f'$\\prod_{{{m.group(1)}}}$', 'Π_{...} → $\\prod_{...}$'),
    
    # --- Special characters ---
    (r'(?<!\$)†(?![\$}])', r'$^\dagger$', '† → $^\\dagger$'),
    (r'(?<!\$)×(?![\$}])', r'$\times$', '× → $\\times$'),
    (r'(?<!\$)≈(?![\$}])', r'$\approx$', '≈ → $\\approx$'),
    (r'(?<!\$)≥(?![\$}])', r'$\ge$', '≥ → $\\ge$'),
    (r'(?<!\$)≤(?![\$}])', r'$\le$', '≤ → $\\le$'),
    (r'(?<!\$)≠(?![\$}])', r'$\neq$', '≠ → $\\neq$'),
    (r'(?<!\$)±(?![\$}])', r'$\pm$', '± → $\\pm$'),
    (r'(?<!\$)√(\w+)', lambda m: f'$\\sqrt{{{m.group(1)}}}$', '√x → $\\sqrt{x}$'),
    (r'(?<!\$)∥([^∥]+)∥', lambda m: f'$\\lVert {m.group(1)}\\rVert$', '||x|| → $\\lVert x\\rVert$'),
    
    # --- Superscripts (digits only, not in math) ---
    (r'(?<!\$)([a-zA-Zα-ω])²(?![\$}])', r'$\1^2$', 'x² → $x^2$'),
    (r'(?<!\$)([a-zA-Zα-ω])³(?![\$}])', r'$\1^3$', 'x³ → $x^3$'),
    (r'(?<!\$)2\^k(?![\$}])', r'$2^k$', '2^k → $2^k$'),
    (r'(?<!\$)p\^k(?![\$}])', r'$p^k$', 'p^k → $p^k$'),
    
    # --- Subscripts (digits and simple patterns) ---
    # v_i → $v_i$ (but not inside existing math or code)
    (r'(?<!\$)([a-zA-Z])_(\d)(?![\$}\w])', r'$\1_\2$', 'x_0 → $x_0$'),
    (r'(?<!\$)([a-zA-Z])_\{([^}]+)\}(?![\$])', lambda m: f'${m.group(1)}_{{{m.group(2)}}}$', 'x_{...} → $x_{...}$'),
    
    # --- Matrix sizes ---
    (r'(?<!\$)(\d+)×(\d+)(?![\$])', r'$\1 \\times \2$', 'd×d → $d \\times d$'),
    
    # --- Arrows ---
    (r'(?<!\$)↔(?![\$}])', r'$\leftrightarrow$', '↔ → ↔'),
    (r'(?<!\$)⇒(?![\$}])', r'$\Rightarrow$', '⇒ → ⇒'),
    
    # --- Misc ---
    (r'(?<!\$)∈(?![\$}])', r'$\in$', '∈ → $\\in$'),
    (r'(?<!\$)∝(?![\$}])', r'$\propto$', '∝ → $\\propto$'),
    (r'(?<!\$)⋯(?![\$}])', r'$\cdots$', '⋯ → $\\cdots$'),
    (r'(?<!\$)…(?![\$}])', r'$\ldots$', '… → $\\ldots$'),
    (r'(?<!\$)⨂(?![\$}])', r'$\otimes$', '⨂ → $\\otimes$'),
    (r'(?<!\$)⊗(?![\$}])', r'$\otimes$', '⊗ → $\\otimes$'),
    (r'(?<!\$)⊕(?![\$}])', r'$\oplus$', '⊕ → $\\oplus$'),
    (r'(?<!\$)∤(?![\$}])', r'$\\nmid$', '∤ → $\\nmid$'),
    (r'(?<!\$)∣\s*\(([^)]+)\)', lambda m: f'$\\mid ({m.group(1)})$', '∣(p−1) → $\\mid (p-1)$'),
    (r'(?<!\$)∣\s*(\w+)', lambda m: f'$\\mid {m.group(1)}$', '∣(p-1) → $\\mid (p-1)$'),
    (r'(?<!\$)⋯(?![\$}])', r'$\cdots$', '⋯ → $\\cdots$'),
    (r'(?<!\$)∞(?![\$}])', r'$\infty$', '∞ → $\\infty$'),
    
    # --- Special: fix double-dollar wrapping ---
    # If we accidentally created $$...$$ inside existing math, fix it
    (r'\$\$([^$]+)\$\$', lambda m: f'$${m.group(1)}$$', 'Fix potential double-wrapping'),
]

# ============================================================
# Protection: Don't modify text inside existing math blocks
# ============================================================

def protect_math_blocks(text):
    """Replace $$...$$ and $...$ blocks with placeholders."""
    placeholders = []
    
    # Protect display math $$...$$
    def save_display(m):
        placeholders.append(m.group(0))
        return f'<<<MATH_DISPLAY_{len(placeholders)-1}>>>'
    text = re.sub(r'\$\$[^$]+\$\$', save_display, text, flags=re.DOTALL)
    
    # Protect inline math $...$ (but not $$)
    def save_inline(m):
        # Skip if it's a display math placeholder
        if '<<<MATH_DISPLAY_' in m.group(0):
            return m.group(0)
        # Skip if it's currency (e.g., $100)
        if re.match(r'\$\d', m.group(0)):
            return m.group(0)
        placeholders.append(m.group(0))
        return f'<<<MATH_INLINE_{len(placeholders)-1}>>>'
    text = re.sub(r'\$[^$]+\$', save_inline, text)
    
    return text, placeholders

def restore_math_blocks(text, placeholders):
    """Restore math block placeholders."""
    for i, ph in enumerate(placeholders):
        text = text.replace(f'<<<MATH_DISPLAY_{i}>>>', ph)
        text = text.replace(f'<<<MATH_INLINE_{i}>>>', ph)
    return text

# ============================================================
# Main
# ============================================================

def fix_math_in_file(filepath):
    """Fix all math notation in a single file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        text = f.read()
    
    original = text
    
    # Step 1: Protect existing math blocks
    text, placeholders = protect_math_blocks(text)
    
    # Step 2: Apply all replacement rules (only on non-math text)
    for pattern, replacement, desc in RULES:
        old_text = text
        try:
            text = re.sub(pattern, replacement, text)
        except Exception as e:
            print(f'  WARNING: Rule "{desc}" failed: {e}')
            continue
    
    # Step 3: Restore protected math blocks
    text = restore_math_blocks(text, placeholders)
    
    # Step 4: Clean up: remove empty $ $
    text = re.sub(r'\$\s+\$', '', text)
    
    # Count changes
    changes = sum(1 for a, b in zip(original.split('\n'), text.split('\n')) if a != b)
    
    # Write back
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(text)
    
    return changes

if __name__ == '__main__':
    import os
    
    files = [
        '/Users/arthas/git/quantum/reports/layperson_guide.md',
        '/Users/arthas/git/quantum/reports/paper_zh_v2.md',
        '/Users/arthas/git/quantum/reports/experiment_complete_report.md',
    ]
    
    for fp in files:
        if os.path.exists(fp):
            print(f'\nProcessing: {fp}')
            changes = fix_math_in_file(fp)
            print(f'  Changed {changes} lines')
        else:
            print(f'  SKIP (not found): {fp}')
    
    print('\nDone!')
