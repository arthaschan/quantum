#!/usr/bin/env python3
"""
将 qa2.md 的独特内容合并到 qa.md。
策略：
- 噪声信道 → 追加到 Part 2 (量子力学基础)
- 硬件路线 + 纠缠分析 → 追加到 Part 5 (RadonShadow)
- 线性反演 vs MLE 细节 → 追加到 Part 3 (QST)
- 连续变量 Radon 变换 → 追加到 Part 1 (信号处理)
- 密度矩阵纯/混态 → 追加到 Part 2
"""
import re

with open('/Users/arthas/git/quantum/reports/qa2.md', 'r', encoding='utf-8') as f:
    qa2 = f.read()

with open('/Users/arthas/git/quantum/reports/qa.md', 'r', encoding='utf-8') as f:
    qa = f.read()

# 清理 qa2 中的多余标题——将多篇文章的 # 标题改为 ## 或 ###
# 已有的 ## 标题保留，新的 # 标题需要改为 ##
sections = qa2.split('\n# ')
header = sections[0]  # 第一篇的第一行

# 处理每一篇独立文章
clean_sections = []
clean_sections.append(header.split('\n', 1)[1] if '\n' in header else header)  # 去掉第一个 # 标题但保留内容

for sec in sections[1:]:
    if not sec.strip():
        continue
    lines = sec.split('\n')
    title = lines[0].strip()
    content = '\n'.join(lines[1:])
    
    # 映射到 qa.md 的标准编号
    # 这些是 qa2 独有的主题
    if '噪声信道' in title or 'Kraus' in title:
        num = "31"  # 追加到 Part 2 后面
        clean_sec = f"## {num}. {title}\n\n{content}"
    elif '线性反演' in title or 'MLE' in title:
        num = "32"  
        clean_sec = f"## {num}. {title}\n\n{content}"
    elif '纠缠' in title:
        num = "33"
        clean_sec = f"## {num}. {title}\n\n{content}"
    elif 'Transmon' in title or '超导' in title or '离子阱' in title:
        num = "34"
        clean_sec = f"## {num}. {title}\n\n{content}"
    elif '合数维扩展' in title:
        num = "35"
        clean_sec = f"## {num}. {title}\n\n{content}"
    elif '密度矩阵' in title or '纯态' in title or '混态' in title:
        num = "8.1"  
        clean_sec = f"### {num} {title}\n\n{content}"
    elif '连续变量' in title or 'Wigner.*相空间' in title:
        num = "4.1"
        clean_sec = f"### {num} {title}\n\n{content}"
    elif '离散变量' in title or '多体系统' in title:
        num = "9.1"
        clean_sec = f"### {num} {title}\n\n{content}"
    elif '总结' in title or '互通性' in title or 'AI 视角' in title:
        # 跳过小结
        continue
    elif '映射总览' in title:
        continue  # qa already has this
    elif '态空间' in title or '酉演化' in title or '测量公理' in title or '状态重构' in title:
        continue  # sub-sections preserved in their parent section
    else:
        # 保留为子节
        num = len([l for l in clean_sections if l.startswith('## ')]) + 1
        clean_sec = f"## {num}. {title}\n\n{content}"
    
    clean_sections.append(clean_sec)

# 合并到 qa.md
new_content = '\n\n'.join(clean_sections)

# 将新内容追加到 qa.md 末尾
qa_merged = qa.rstrip() + '\n\n---\n\n## 补充内容（来自 qa2.md）\n\n' + new_content

with open('/Users/arthas/git/quantum/reports/qa.md', 'w', encoding='utf-8') as f:
    f.write(qa_merged)

print(f"Merged: qa.md now {qa_merged.count(chr(10))} lines")
print(f"Added {len(clean_sections)} new sections from qa2.md")
