#!/usr/bin/env python3
"""清理 qa2.md 的激励话语和多余的目录，准备合并到 qa.md"""
import re

with open('/Users/arthas/git/quantum/reports/qa2.md', 'r', encoding='utf-8') as f:
    text = f.read()

# 删除重复的标题（第二、三个 # 开头的标题改为正确的编号）
# 第一个标题保留，后续的 # 标题需要合并
# qa2有多个独立文章，各自从 # 开始。我们需要保留内容但调整标题层级

# 移除激励性话语
fluff = ['融会贯通', '完美诠释', '终极', '不可思议', '震撼', '巅峰']
for word in fluff:
    text = re.sub(rf'^.*?{word}.*?$\n?', '', text, flags=re.MULTILINE)

# 移除多余空行
text = re.sub(r'\n{4,}', '\n\n\n', text)

with open('/Users/arthas/git/quantum/reports/qa2.md', 'w', encoding='utf-8') as f:
    f.write(text)

print(f"Cleaned qa2.md: {text.count(chr(10))} lines")
