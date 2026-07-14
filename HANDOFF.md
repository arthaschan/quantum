# RadonShadow 项目交接文档 · 2026-07-15 01:56 CST

> **写给接手此项目的 AI**：按本文档操作即可从零理解全部上下文，无需翻看会话历史。

---

## 0. 项目一句话

**RadonShadow**：发现有限域上 MUB（互无偏基）量子阴影测量与 DPRT（离散投影 Radon 变换）在代数上等价，因此可以用确定性遍历替代随机采样来重建量子态——在素数维度上中位数快 13.6%，最强快 67.3%。

用类比讲：CT 扫描的 Radon 变换 → DPRT 离散版 → 恰好等于量子 MUB 测量的傅里叶变换 → 于是确定性 CT 重建算法可以直接用在量子态重建上。

---

## 1. 用户是谁

| 项 | 值 |
|:--|:--|
| 姓名 | 未记录（GitHub: arthaschan） |
| 身份 | 香港中文大学 AI 硕士生 |
| 导师 | Dr. Richard Tai-Chiu Hsung（熊体操教授）——DPRT 核心论文 [1996 IEEE TSP] 第二作者 |
| 数学水平 | 数学系本科毕业，线性代数/群论入门/概率论学过，但**没有量子力学背景** |
| 核心诉求 | 完全理解论文里的每个定理和公式，逐项验证，标注知识差距 |
| 风格偏好 | 诚实 > 漂亮。可以说不懂，不能假装懂。能指出错误 > 泛泛夸奖 |
| 操作系统 | macOS (Darwin, arm64) |
| 语言 | 中文为主，英文文献 |
| 有其他研究方向 | ε+π 无理性、Ramsey 数 R(5,5)、DPRT+量子断层扫描相关 |

---

## 2. 工作区结构

```
/Users/arthas/.qclaw/workspace-radonshadow/
├── AGENTS.md
├── RESEARCH_ROADMAP.md          # 四条研究路线总计划
├── MEMORY.md                    # 项目记忆文件
├── data/                        # 全部实验 JSON 数据（17 个文件）
│   ├── w1_experiment_results.json           # MUB↔DPRT 数值等价验证
│   ├── primes_1000_results.json             # 168 素数扫描（最重要）
│   ├── primes_10000_sampled.json            # 32 个大素数采样
│   ├── primes_23_97_results.json            # 选定素数深度扫描
│   ├── all5_part2_results.json              # 2-local 可观测量
│   ├── composite_100_results.json           # 旧 CRT 方法（已废弃）
│   ├── composite_crt_results.json           # CRT 分解尝试（已废弃）
│   ├── composite_local_mub_results.json     # 正确合数维方法（74 合数）
│   ├── noise_adaptive_results.json          # 160 点噪声实验
│   ├── adaptive_relaxed_results.json        # 自适应 ε 放宽扫描
│   ├── galois_orbit_results.json            # 45 素数 Galois 轨道分析
│   ├── g2_stabilizer_results.json           # 原根敏感性（16 素数×5 原根）
│   ├── g3_optimal_subset_results.json       # 最优 MUB 子集搜索
│   ├── d13_19_placeholder_results.json
│   ├── qdprt_experiment_phase3.json
│   └── (其他小文件)
├── scripts/                     # 全部实验脚本（16 个 Python）
│   ├── prime_scan_1000.py       # E1: 素数扫描
│   ├── prime_scan_10000.py      # E3: 大素数
│   ├── all5_experiments.py      # E1-E5 合并运行器（含 primitive_root 函数）
│   ├── composite_local_mub.py   # E5-新: 正确合数维方法
│   ├── composite_crt.py         # 旧 CRT（已废弃，仅保留参考）
│   ├── noise_adaptive.py        # E7: 噪声鲁棒性
│   ├── adaptive_relaxed.py      # 自适应 ε 放宽
│   ├── w1_experiments.py        # W1: 形式化验证
│   ├── e2_subgroup.py           # E2: 原根分组
│   ├── e3_composite.py / e3_extended.py  # 合数维初版
│   ├── (群论 G1/G2/G3 脚本由子进程生成，路径不详)
│   └── (其他辅助脚本)
├── reports/                     # 全部报告（21 个文件）
│   ├── layperson_guide.md       # ⭐ 给数学本科的通俗版 + 差距分析 (540 行)
│   ├── experiment_complete_report.md  # ⭐ 12 组实验全景汇总 (667 行)
│   ├── paper_zh_v2.md           # ⭐ 论文中文版 v2 (976 行, 最新)
│   ├── paper_zh_v1.md           # 论文 v1 (735 行, 旧版)
│   ├── paper_draft_v1.md        # 中英双语初稿
│   ├── figure_prompts.md        # 6 张框架图 prompt
│   ├── noise_adaptive_report.md # 噪声实验独立报告
│   ├── composite_local_mub_report.md   # 合数维独立报告
│   ├── composite_crt_report.md  # CRT 方法报告（废弃方法，仅记录）
│   ├── kingston_formalization_20260714.md  # W1 形式化细节
│   ├── research_diary.md        # 研究日记
│   ├── experiment_reports.md    # 早期实验汇总
│   ├── RESEARCH_ROADMAP.md      # 研究路线总计划（根目录也有同名文件）
│   └── (其他历史报告)
└── RESEARCH_ROADMAP.md          # 四条研究路线总计划
```

---

## 3. 本会话做了什么

### 用户要求三份文档
"我需要这三份报告"，此后我并行启动三个子进程写文档：

| # | 文档 | 行数 | 面向 | 内容 |
|:--|:--|:--:|:--|:--|
| 1 | `reports/experiment_complete_report.md` | 667 | 研究者+导师 | 12 组实验全景：E1-E7、G1-G3、W1、E5-新，含数据表、关系图、导师摘要 |
| 2 | `reports/paper_zh_v2.md` | 976 | 论文投稿 | v1(735行)→v2：新增§7群论框架、§5.5.1合数维突破、完整噪声§5.7、6条Reviewer预判 |
| 3 | `reports/layperson_guide.md` | 540 | 数学本科 | 零量子力学、CT类比、p=3手算验证、三级差距清单、诚实标注边界 |

### 文档之间关系
- **通俗版** → 用户先读，理解每个概念
- **实验报告** → 理解数据怎么来的，查数据用
- **论文 v2** → 最终交付品，集成了所有发现

---

## 4. 实验全景（快速索引）

| 代号 | 名称 | 数据文件 | 核心数值 |
|:--|:--|:--|:--|
| W1 | MUB↔DPRT 等价验证 | `w1_experiment_results.json` | 等价性数值验证通过，加速 400-500× |
| E1 | 素数扫描 d≤1000 | `primes_1000_results.json` | 168 primes，中位 ratio=0.864，63.7% ratio<0.90 |
| E2 | 原根分组 | `primes_1000_results.json` | g=2: 14.3%, g=3: 15.4% |
| E3 | 大素数 | `primes_10000_sampled.json` | d=2039 ratio=0.327(+67.3%)，ρ=-0.8913 |
| E4 | 2-local | `all5_part2_results.json` | 中位 0.992，d=17 ratio=0.460 |
| E5-旧 | 合数 CRT | `composite_100_results.json` | ⛔ 已废弃（方法有缺陷） |
| E5-新 | 合数 Local MUB | `composite_local_mub_results.json` | 74 合数，d=16 ratio=0.401(+59.9%) ⭐ |
| E6 | 态类型效应 | 论文 §5.6 | 纯态 0.805，混态 1.041 |
| E7 | 噪声 | `noise_adaptive_results.json` | Phase Damping +35.3%，Depol 最稳 |
| G1 | Galois 轨道 | `galois_orbit_results.json` | 失败素数共因：3\|(p-1) |
| G2 | 原根敏感性 | `g2_stabilizer_results.json` | 9/16 素数可翻转胜负 |
| G3 | 最优子集 | `g3_optimal_subset_results.json` | 5/5 全部可挽救 |

---

## 5. 核心定理（用户需要逐条验证）

### T1: MUB-DPRT 代数等价定理
MUB 投影概率的 DFT = DPRT 投影值。这步可在 p=3 上手动验证。
→ 验证方法见 `layperson_guide.md` §2.3

### T2: 方差上界
确定性方案的 MSE 上界 = d(d+1)²/N，与经典 MUB 阴影同阶。
→ 紧度因子范围 1.2–7.5×，论文 §3.4 有验证表

### T3: 收敛性
- T3a: 全态平均 → ratio → 1，O(1/√d)
- T3b: 1-local 固定态 → ratio 振荡，不收敛于 1（这是优势根源）
- T3c: k≥2 → ratio → 1，O(1/3^k)

---

## 6. 关键决策与注意事项

### ⚠️ 原根选错？——不，没选错
- 原始实验 (E1) 用**线性遍历** m=0,1,...,d，不依赖任何原根
- 因此 63.7% 优势率是**保守下界**
- G2 实验证明了更好的遍历顺序存在（p=41: 线性→0.809，最差 g→1.105，最优 g→0.809）
- 论文 v2 中已添加 Referee #6 回应（§6.2.2）

### ⚠️ 合数维：CRT 方法已废弃
- 旧 CRT 方法假设「各素因子密度矩阵独立重建后 Kronecker 积拼回」→ 只对可分离态成立
- 正确方法：局部 MUB 乘积（Local MUB Product）— 每个子系统独立测
- 结果：74 个合数有效，偶数首次出现优势（d=16 ratio=0.401）

### ⚠️ 态纯度是混淆因子
- 纯态放大了 DPRT 优势（pure ratio=0.805 vs mixed=1.041）
- 解释：纯态的 Wigner 函数有尖锐相位结构，DPRT 确定性投影更有效

### 其他注意点
- d=4（2-qubit）持平→MUB 不完备，非方案缺陷
- d=64（2⁶）信息墙崩溃→2^k 的 k=6 处出现
- 适度噪声充当正则化器（零噪声≠最优）
- 自适应提前终止收益有限（最多 6.2%，通常 2.5%）→ 不建议作为核心创新

---

## 7. 用户的待办清单

### 立即（已交付）
- [x] 完整实验综合报告
- [x] 完善论文 v2
- [x] 数学本科可读的通俗版 + 差距分析

### 今天/近期
- [ ] 读 `layperson_guide.md` 全文（约 45 分钟）
- [ ] Lv1 验证：p=3 手算 MUB 基→投影概率→DPRT→验证 DFT(P)=R
- [ ] Lv2 验证：运行 prime_scan_1000.py，确认 ratio≈0.864
- [ ] Lv3 验证：读 Wootters-Fields 1989 原论文
- [ ] 按 §6 差距清单开始补知识

### 中期（论文相关）
- [ ] 把论文 v2 转成 IEEE DOCX（md2ieee-docx skill）
- [ ] 制作 6 张框架图（figure_prompts.md 里有 prompt）
- [ ] 扩展文献综述（论文只引了 8 篇，建议 25+）
- [ ] 联系导师审阅

### 博士级别扩展
- [ ] 合数维环上推广（Galois 环而非有限域）
- [ ] 硬件实验（d=3 qutrit 超导平台）
- [ ] 群论框架：Mackey 理论形式化

---

## 8. 用户可能的下一步请求

以下请求高概率出现，请参照相应文件：

| 请求 | 参照文件 | 操作 |
|:--|:--|:--|
| "帮我理解 T1 证明" | `layperson_guide.md` §2，`paper_zh_v2.md` §3.3 | 一行一行过，用 p=3 例子 |
| "我没学过 Wigner 函数" | `layperson_guide.md` §1.2 给了纯数学定义 | 避免物理解释，用矩阵语言 |
| "这个数据是编的还是算出来的？" | `experiment_complete_report.md` → 对应数据文件 | 用 python3 读 JSON 当场展示 |
| "为什么 p-1 质因子越多越好？" | `layperson_guide.md` §3.3 → 钟摆类比 | Spearman ρ=-0.8913 是证据 |
| "合数维到底怎么做？" | `layperson_guide.md` §4.2，`composite_local_mub_report.md` | 讲 Local MUB Product 方法 |
| "我要把论文转成 Word 格式" | md2ieee-docx skill | 用 paper_zh_v2.md 作为输入 |
| "我要读论文原文" | HKP 2020 = Huang-Kueng-Preskill, Wootters-Fields 1989 | 从 arXiv 获取 |
| "帮我做实验 X" | `scripts/` 目录下对应脚本 | 用 python3 运行 |

---

## 9. 知识差距（从通俗版 §6 提取）

用户自己标记的需要补充的知识：

### 数学（9 项）
1. 有限域 F_p 群论 → Dummit & Foote Ch.13-14 (1-2天)
2. Clifford 群 → Gottesman stabilizer 笔记 (2-3天)
3. Wigner 函数 → Wootters 1987 (1-2天)
4. 影子范数 → HKP 2020 附录 (2-3天)
5. Weil 特征标 → Ireland-Rosen Ch.11 (1-2天)
6. 测量信道 M 算符 → HKP 2020 正文 (1-2天)
7. Galois 环调和分析 → Klappenecker & Rötteler 2005 (2-3周)
8. DPRT 卷积反演 → Kingston 2006 (1-2天)
9. Mackey 表示论 → Serre (3-4周)

### 物理（4 项）
1. 量子力学基本公理 (1周)
2. 密度矩阵 (2-3天)
3. 量子噪声信道 (2-3天)
4. 超导硬件 (1-2天)

---

## 10. 关键引用

论文 v2 的 8 篇引用（需扩展至 25+）：
- [1] Huang-Kueng-Preskill 2020 (Nature Physics): Classical Shadows
- [2] Wootters-Fields 1989: MUB 构造
- [3] Hsung-Lun-Siu 1996 (IEEE TSP): DPRT 原始论文（导师的！）
- [4] Kingston-Svalbe 2006 (IEEE TIP): 周期图像投影变换
- [5] Kingston-Svalbe 2007 (I&VC): 广义有限 Radon 变换
- [6] Ivanović 1981: MUB 量子态确定
- [7] Klappenecker-Rötteler 2005 (ISIT): MUB = 复投影 2-设计
- [8] Zhu et al.: 待补充完整引用

---

## 11. 会话结束状态

- 用户最后一句话："更新到layperson_guide.md 里" → 已将差距清单从 13 项扩充至 21 项
- 本文件就是交接文档（创建于 01:56，更新于 02:03 CST）
- 三个子进程已全部完成，文件均写入正确位置
- layperson_guide.md 从 540 行 → 586 行（新增 Bra-Ket、张量积、数论、统计推断、数值分析、断层扫描基线、纠缠、实验脚本阅读、JSON 结构、原根机制、真效应判别、3 周自学路线）
- 实验 ✅，论文 ✅，通俗版 ✅ + 差距清单补全 ✅——本阶段全部完成
- 没有运行中的子进程或挂起的任务
