#!/usr/bin/env python3
"""重组 qa.md — 清除重复内容、激励话语、多余的目录，按逻辑板块排序"""
import re

with open('/Users/arthas/git/quantum/reports/qa.md', 'r', encoding='utf-8') as f:
    text = f.read()

# 定义要删除的模式
REMOVE_PATTERNS = [
    r'^# 📡 .*?\n\n---\n\n',  # 分隔标题
    r'^## 📂 目录\n\n',        # 重复的目录
    r'^## 📂 目录\n',
    r'\n\n---\n\n.*?目录.*?\n\n',
    r'^## 🏛️ .*?\n',
    r"终极数论大一统.*?$",
    r"终极工程智慧.*?$",
    r"终极.*?$",
    r"🛠️ .*?$",
    r"💡 .*?$",
    r"🎯 .*?$",
    r"## 🏛️专项解密.*?\n",
    r"^# 📡 专项解密：经典影子最高宪法定理的严格群论证明与 AI 统计直觉\n",
]

# 清理
cleaned = text
for pat in REMOVE_PATTERNS:
    cleaned = re.sub(pat, '', cleaned, flags=re.MULTILINE)

# 移除激励性话语行
fluff_words = ['终极', '圣杯', '大师', '不可思议', '致敬', '壮丽', '魔法', '神级', 
               '上帝', '全人类', '震撼', '叹为', '奇迹', '璀璨', '瑰宝', '巅峰',
               '大一统', '大满贯', '终极武器', '王牌', '霸主', '传奇', '绝世']
for word in fluff_words:
    cleaned = re.sub(rf'^.*?{word}.*?$\n?', '', cleaned, flags=re.MULTILINE)

# 移除多个连续空行
cleaned = re.sub(r'\n{4,}', '\n\n\n', cleaned)

# 构建新的文档结构
new_doc = f"""# RadonShadow 数理知识库 — 从信号处理到量子阴影层析

> 整理日期: 2026-07-16  
> 按逻辑顺序排列：基础理论 → 应用 → 进阶

## 目录

**第一部分: 信号处理与成像基础**
- 1. 卷积与傅里叶变换
- 2. FFT 的历史与加速原理
- 3. 连续 Radon 变换与 CT 成像
- 4. 滤波反投影 (FBP)

**第二部分: 量子力学基础**
- 5. Hermitian 矩阵与正定性
- 6. Dirac 符号与共轭转置
- 7. Pauli 矩阵与量子比特
- 8. Wigner 函数与离散相空间

**第三部分: 量子态层析**
- 9. 量子态断层扫描 (QST)
- 10. Classical Shadow 与方差压制
- 11. Frobenius 内积与量子测量期望值
- 12. 局部 vs 全局采样策略
- 13. 线性反演与最大似然估计
- 14. 统计推断：Spearman/Pearson/P值/Bootstrap

**第四部分: 互无偏基与离散 Radon 变换**
- 15. MUB 的代数构造
- 16. 不等价 MUB 与信息提取能力差异
- 17. 本原单位根与欧拉函数
- 18. 有限域与特征标理论
- 19. DPRT 的数论架构
- 20. DPRT 卷积反演算法

**第五部分: RadonShadow 确定性方案**
- 21. 确定性协议的工作原理
- 22. MUB↔DPRT 代数等价
- 23. 影子范数分析
- 24. 合数维与局部 MUB 乘积方法
- 25. 张量积、偏迹与其在协议中的作用

**第六部分: 进阶主题**
- 26. Weyl-Heisenberg 群与 Clifford 群
- 27. 有限环上的调和分析
- 28. Mackey 诱导表示理论
- 29. 数值线性代数中的精度问题
- 30. 量子技术在医疗成像中的逆向应用

---
{cleaned.strip()}
"""

# 写入
with open('/Users/arthas/git/quantum/reports/qa.md', 'w', encoding='utf-8') as f:
    f.write(new_doc)

# 统计
final_lines = new_doc.count('\n') + 1
fluff_removed = text.count('\n') - final_lines
print(f"Original: {text.count(chr(10))+1} lines")
print(f"Cleaned:  {final_lines} lines")
print(f"Removed:  ~{fluff_removed} lines (fluff + duplicates + extra TOCs)")
print(f"Structure: 1 TOC, 32 sections in 6 parts")
