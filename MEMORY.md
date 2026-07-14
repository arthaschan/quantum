- **2026-04-14**：记忆系统启用
- **2026-07-14**：RadonShadow 项目启动，工作空间切换到 /Users/arthas/.qclaw/workspace-radonshadow
- **2026-07-14 22:30**：W1/W2/素数全扫描/All-5实验/合数维全覆盖全部完成
- **2026-07-15 00:20**：文档合并（11→3）+ 文献空白确证 + 论文初稿子进程启动

## 用户身份与偏好

- 论文导师Dr. Richard Tai-Chiu Hsung（熊体操教授），在香港中文大学
- GitHub用户名为arthaschan
- 有一个 DeepSeek API Key
- 在香港中文大学读人工智能硕士
- 研究领域：量子层析（RadonShadow）+ AI医学 + R(5,5)
- 引用格式偏好：in ICLR 年份 简写格式

## RadonShadow 项目状态（2026-07-15 00:20）

### 已完成（Phase 1 + 2）
- ✅ 混合态 vs 纯态效应（d=3,4,8 全覆盖）
- ✅ 大 d 扫描（d=2..97 全素数，d≤1000 全素数 168p，1000-10000 采样 32p）
- ✅ 合数维度全覆盖（d≤100, 25 valid）
- ✅ E1-E5 五项实验（原根、子群、合数维、2-local、噪声）
- ✅ W2 定理形式化（T1-T3 + G1-G3 Kingston）
- ✅ 后处理复杂度 benchmark（400-500× 加速）
- ✅ 文档合并（prime_scan_reports, phase1_2_reports, experiment_reports）
- ✅ 文献空白确证（literature_gap_20260715.md）

### 论文核心数据
- **素数 median 1loc ratio = 0.864**（DPRT 优 13.6%），63.7% 维度 ratio<0.90
- 最佳 d=607: 42.4% 优势，d=2039: 67.3%
- 原根=2/3: median 优势 14.3%/15.4%
- 合数维 20% DPRT 优，d=63 ratio=0.363
- 去极化噪声 λ=0.05 甜区: DPRT 优 6.4%

### 进行中
- 🔄 **论文初稿 v1** — 子进程 paper_draft_v1 运行中（session: ca74c8e6）

### Novelty 确认
- DPRT↔MUB 代数等价 + 确定性 Classical Shadow 协议 = **空白地带**
- F_d 原根→测量基优化 的因果链 = **全新**
- Klappenecker & Rötteler 2004 不冲突（关注 MUB 构造存在性，非使用效率）

## QClaw 积分状态 (2026-07-15)

- 当前剩余：约 780 积分
- 当前模型：qclaw/pool-deepseek-v4-pro

## 技术规范偏好

- 用户使用 VPN 代理软件 Clash（Clash Party）
- 终端代理已配置：http_proxy/https_proxy/ALL_PROXY 走 Clash 7890/7892 端口
