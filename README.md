# Bat1K 蝙蝠免疫基因组论文：方法审计与可复现优化

本仓库是对以下论文的计算流程进行**方法审计、代码重写与可复现化**的独立工作，不是作者官方仓库：

> Morales, A. E. et al. **Bat genomes illuminate adaptations to viral tolerance and disease resistance.** *Nature* 638, 449–458 (2025). DOI: [10.1038/s41586-024-08471-0](https://doi.org/10.1038/s41586-024-08471-0)

原始作者代码与数据入口：

- 作者流程说明仓库：[ariadnamorales/2023_Bat1Kimmunity](https://github.com/ariadnamorales/2023_Bat1Kimmunity)
- 作者回归代码：[lmdavalos/count2branches](https://github.com/lmdavalos/count2branches)
- 公共数据下载：[genome.senckenberg.de/download/Bat1KImmune](http://genome.senckenberg.de/download/Bat1KImmune/)

> **状态声明（重要）**
> - 本仓库中的流程**尚未在论文的 17,130 个基因或 228 个分支上运行**，因此不包含、也不宣称任何复现结果。
> - `examples/toy/` 下的序列、树、分支计数和 ZOVER 记录全部是**人工合成的测试数据**，只用于检查文件格式和代码管线，不能用于生物学推断。
> - 论文的公共仓库主要是“命令清单”，不是可直接执行的流水线；本仓库把缺失的逐基因调度、aBSREL JSON 解析、分支回溯、富集和回归诊断补齐。
> - 论文数据（尤其是 495 MB 的 codon alignment 包）没有纳入本仓库。

## 这个仓库包含什么

```
bat1k-repro/
├── workflow/
│   ├── Snakefile                 # Snakemake 8 逐基因 aBSREL 流水线
│   ├── config/config.yaml        # 真实数据配置
│   ├── config/config.toy.yaml    # 合成玩具配置
│   ├── envs/                     # conda 环境（Python / HyPhy / R）
│   ├── profiles/slurm/           # SLURM 配置
│   └── scripts/
│       ├── make_branch_metadata.py   # 树节点命名、分支元数据
│       ├── screen_gene.py            # QC → 剪枝 → aBSREL → JSON 解析
│       ├── aggregate_absrel.py       # 分支/基因/运行汇总
│       ├── make_branch_counts.py     # 每个分支的 GO 入选基因计数
│       └── common.py                 # FASTA、QC、BH、aBSREL 辅助函数
├── R/
│   ├── count2branches_optimized.R        # 负二项回归 + 诊断 + 敏感性分析
│   ├── enrichment_gprofiler_optimized.R  # 自定义背景 vs 默认背景富集
│   └── fisher_zover_optimized.R          # Fisher + 研究聚类 GEE/bootstrap
├── examples/toy/                 # 合成数据与生成脚本
├── tests/                        # 纯函数单元测试（不调用 HyPhy/R）
└── docs/
    ├── 01_methods_audit_zh.md    # 论文方法与证据链审计
    ├── 02_workflow_zh.md         # 原代码流程与“为什么这样做”
    ├── 03_findings_zh.md         # 代码/统计问题清单与修复
    └── 04_tutorial_zh.md         # 安装、输入格式、运行、排错、验证
```

## 快速开始

### 1. 安装

需要 conda/mamba、Snakemake ≥ 8。建议先创建一个小环境运行 Snakemake：

```bash
mamba create -n bat1k-smk -c conda-forge -c bioconda snakemake-minimal=8
mamba activate bat1k-smk
```

各规则会自动创建 `workflow/envs/*.yaml` 中的环境（首次运行需要网络）。

### 2. 玩具数据（只验证管线，不做推断）

```bash
# 仅重新生成玩具 fixture 时需要；仓库已包含生成结果
python examples/toy/make_toy_data.py

# 离线核心：QC、剪枝、aBSREL、汇总
snakemake -s workflow/Snakefile \
  --configfile workflow/config/config.toy.yaml \
  --cores 2 --use-conda all

# 下游：GO 计数、g:Profiler（需网络）、贝叶斯回归
snakemake -s workflow/Snakefile \
  --configfile workflow/config/config.toy.yaml \
  --cores 2 --use-conda downstream
```

玩具数据只有 6 个物种、3 个基因、9 个分支；回归使用 `min_branches: 8` 且链很短，**结果没有统计意义**。

### 3. 真实数据

把数据放到配置指定的路径（默认在仓库根目录下）：

```
alignments/{gene}.fa              # 密码子对齐，长度是 3 的倍数
resources/species_tree.nwk        # 分支长度必须是 substitutions/site
resources/species_metadata.tsv    # taxon, order, is_bat
resources/gene_map.tsv            # gene_id, ensembl_gene_id
resources/gene2go.tsv             # gene_id, go_id
```

```bash
# 先 dry-run 检查 DAG
snakemake -s workflow/Snakefile --cores 8 --use-conda -n all

# 正式运行
snakemake -s workflow/Snakefile --cores 16 --use-conda all
snakemake -s workflow/Snakefile --cores 16 --use-conda downstream

# SLURM
snakemake -s workflow/Snakefile --profile workflow/profiles/slurm --use-conda all
```

详细输入格式、输出解释、运行时间和排错见 [docs/04_tutorial_zh.md](docs/04_tutorial_zh.md)。

## 关键优化（相对作者公开代码）

| 主题 | 原代码 | 本仓库 |
|---|---|---|
| 运行单位 | 一段 README 命令，无逐基因调度 | Snakemake 规则，单基因失败不影响全局，状态 JSON 可审计 |
| 分支映射 | aBSREL 节点标签直接使用，内部节点可能无标签/重名 | 先给全树内部节点唯一命名，剪枝后用“后代叶集合 → 全树 MRCA”回溯 |
| 输入检查 | 未显式检查树长单位、内部终止密码子、非法字符 | QC + 强制 substitutions/site + 状态记录 |
| 多重检验 | 仅在基因内对分支做 BH | 保留原定义，另加**可选**的基因级 BH 敏感性列 |
| 富集背景 | 默认全部人类注释基因 | 同时跑“该目实际可检验基因”自定义背景和默认背景，直接比较 |
| 回归模型标签 | `x:taxon` 被注释成“双截距双斜率” | 明确 `x:taxon` 是单截距双斜率，`x*taxon` 才是双截距双斜率 |
| 回归诊断 | 无 seed、无 Rhat/ESS、无后验预测检查 | 固定 seed、先验、Rhat/ESS/divergence、WAIC/LOO、`pp_check`、敏感性模型 |
| 机会量 | 只用原始计数 | 额外提供 `n_genes_tested` 及 `offset(log(n_genes_tested))` 敏感性模型 |
| ZOVER | 每条病毒记录独立做 Fisher | Fisher + 以 study 为聚类的 GEE + study 聚类 bootstrap |
| 可复现性 | `save.image()` 反复覆盖，无 session info | 每步落盘、`sessionInfo()`、版本记录、conda 环境 |

逐条证据见 [docs/03_findings_zh.md](docs/03_findings_zh.md)。

## 文档

- [01 方法与逻辑审计](docs/01_methods_audit_zh.md)：论文的样本、系统发育、选择压力、富集、回归和湿实验证据链，哪些结论强、哪些是探索性的。
- [02 代码流程](docs/02_workflow_zh.md)：从组装到 ISG15 实验，每一步“为什么做、输入输出是什么、原代码怎么实现”。
- [03 问题与修复](docs/03_findings_zh.md)：按严重程度列出代码/统计问题，附原代码片段与修复策略。
- [04 教程](docs/04_tutorial_zh.md)：安装、输入 schema、玩具/真实/SLURM 运行、输出解读、验证与排错。

## 引用与许可

如果使用本仓库的代码，请同时引用原论文和原始代码仓库。论文数据与原始代码的版权归原作者；本仓库目前**未声明许可证**，在复用前请先与仓库维护者确认。论文中的生物学结论应以原文和补充材料为准，本仓库的审计意见不构成对论文结论的否定。
