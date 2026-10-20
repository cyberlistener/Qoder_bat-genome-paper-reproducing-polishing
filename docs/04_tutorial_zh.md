# 04 · 教程：安装、输入、运行与排错

> 本教程中的所有命令都**尚未在论文数据上执行**。玩具数据命令只验证管线，不产生可解释的生物学结果。
> 真实运行前，请先读 [01 方法审计](01_methods_audit_zh.md) 和 [03 问题清单](03_findings_zh.md)。

## 1. 环境

### 1.1 需要什么

- Linux 或 macOS（Windows 建议 WSL2）；
- conda 或 mamba；
- Snakemake ≥ 8；
- 首次运行需要网络（下载 conda 包；富集需要访问 g:Profiler）；
- 真实数据：磁盘至少几十 GB（alignment + 中间文件 + 结果）；建议 ≥ 16 GB 内存。

### 1.2 创建 Snakemake 环境

```bash
mamba create -n bat1k-smk -c conda-forge -c bioconda snakemake-minimal=8
mamba activate bat1k-smk
snakemake --version
```

各规则会通过 `--use-conda` 自动创建：

| 环境 | 内容 |
|---|---|
| `workflow/envs/screen.yaml` | Python 3.10、HyPhy 2.5.8、DendroPy、Biopython、pandas |
| `workflow/envs/python.yaml` | Python 3.10、pandas、numpy |
| `workflow/envs/r.yaml` | R 4.3、brms、loo、gprofiler2、geepack、ggplot2 等 |
| `workflow/envs/snakemake.yaml` | Snakemake + SLURM executor（可选） |

如果集群不允许 conda，可先手动创建环境，然后去掉 Snakefile 中的 `conda:` 指令，或使用 `--sdm none` 并确保 `hyphy`、`Rscript` 在 PATH 中。

### 1.3 检查 HyPhy

```bash
hyphy --version   # 期望 2.5.8
```

aBSREL 的 JSON schema 在 HyPhy 版本间可能变化。本仓库按 2.5.8 的输出解析；若你用其他版本，先用一个基因检查 `branch attributes`/`tested` 字段，必要时改 `screen_gene.py` 的 `parse_absrel()`。

## 2. 目录与输入格式

### 2.1 目录

```
your-project/
├── alignments/                 # {gene}.fa
├── resources/
│   ├── species_tree.nwk
│   ├── species_metadata.tsv
│   ├── gene_map.tsv
│   └── gene2go.tsv
└── results/                    # 自动生成
```

配置里的路径相对于**运行 snakemake 的目录**（仓库根目录）。不要用绝对路径写死到别人机器上。

### 2.2 密码子对齐 `alignments/{gene}.fa`

要求：

- FASTA，每个基因一个文件，文件名（不含 `.fa`）就是 `gene_id`；
- 所有序列等长；
- 长度必须是 3 的倍数；
- 翻译后无内部终止密码子（除非显式 `allow_internal_stops: true`）；
- 序列 ID 必须与树末端标签、`species_metadata.tsv` 的 `taxon` 完全一致；
- 允许的字符：`ACGTNRYKMSWBDHV-?.`（大小写均可，读入时统一大写）。

示例：

```
>bat1
ATGGCTGATAAATTCCAAGAACTGATTGTTTCTACTCCTCGTGGTCATTATTGTTGGGCA
GCCGCGAAGGACTTTCAGGAGCTCATCGTC
>bat2
...
```

**从论文数据准备**：下载 `codon_alignments_all/ALNs_w1branch_PosSel.tar.gz`（约 495 MB），解压后把每个基因的 `.filtered.fa` 重命名为 `{gene}.fa`。如果压缩包内是嵌套目录，用 `find`/`rename` 展平；不要改变序列内容。论文的 `{file}.raw.fa` 是未清理版本，`{file}.filtered.fa` 才是选择筛查输入。

### 2.3 物种树 `resources/species_tree.nwk`

要求：

- Newick 格式，带分支长度；
- 分支长度单位必须是 **substitutions/site**；
- 末端标签与 alignment/元数据一致；
- 可以 rooted 或 unrooted；脚本用 `rooting="default-rooted"` 读取；
- 不允许负长度；零长度分支会警告，回归时会被 `bl > 0` 过滤。

**最常见错误**：直接把 ASTRAL 树（coalescent units）或 treePL 时间树（Mya）传给 aBSREL。两者都不是 substitutions/site。正确做法是先固定拓扑、重新估计替换数，例如：

```bash
# 仅示意；拓扑来自 ASTRAL，分支长度由 IQ-TREE 在串联/编码区比对下估计
iqtree2 -s concatenated_coding.fasta \
  -te astral_topology.nwk \
  -m GTR+G -nt AUTO --prefix coding_tree
# 输出 coding_tree.treefile 的分支长度是 substitutions/site
```

如果你只有时间树而想复现论文的“Mya”回归，**不要**把时间树用于 aBSREL；aBSREL 和回归可以用两棵树（同一拓扑、不同长度），但必须在结果里写清楚。

### 2.4 物种元数据 `resources/species_metadata.tsv`

```tsv
taxon	order	is_bat
HLrhiAff2	Chiroptera	TRUE
HLhipLar2	Chiroptera	TRUE
HLmopCon1	Chiroptera	TRUE
HLmanPen2	Primates	FALSE
...
```

- `taxon`：与树/FASTA 完全一致；
- `order`：目名；用于富集和分支元数据；
- `is_bat`：`TRUE/FALSE`（也接受 `1/0`, `yes/no`）。

脚本会检查：树末端是否都在元数据中；元数据中多出的物种只警告，不报错。

### 2.5 基因映射 `resources/gene_map.tsv`

```tsv
gene_id	ensembl_gene_id
ENSG00000141510	ENSG00000141510
my_transcript_1	ENSG00000012048
```

- `gene_id`：alignment 文件名/表里的 ID；
- `ensembl_gene_id`：**人类 Ensembl gene ID**（g:Profiler 的 `hsapiens` 查询需要）；
- 如果 alignment 已经用人类 Ensembl gene ID 命名，两列可以相同；
- 脚本会去掉版本号（`ENSG... .12` → `ENSG...`）；
- 一个 gene_id 只保留第一条映射；重复会被静默去重，建议自己先检查。

**怎么得到**：

- TOGA 输出通常带 human gene/transcript 对应关系；把 `gene_id` 映射到 human Ensembl gene ID；
- 或用 Ensembl BioMart / g:Profiler `gconvert` 做 ID 转换；
- 如果使用论文的 alignment 文件名（可能是 transcript ID），必须建这张表，否则富集查询会用 transcript ID，g:Profiler 可能无法识别。

### 2.6 GO 注释 `resources/gene2go.tsv`

```tsv
gene_id	go_id
ENSG00000141510	GO:0002376
ENSG00000012048	GO:0002376
ENSG00000157764	GO:0006955
```

- 长格式；
- `gene_id` 与 alignment 文件名一致；
- `go_id` 支持 `GO:0002376` 或 `GO_0002376`（脚本会统一）；
- 只用于“每个分支上 GO 术语的入选基因计数”（回归输入）；富集本身由 g:Profiler 在线完成，不依赖这张表。

**注意**：论文的回归只用了 GO:0002376（immune system process）。如果你要复现 42 vs 22，必须确保 `gene2go.tsv` 中 GO:0002376 的基因集合与论文定义一致（人类 GO 注释、版本、是否包含 IEA 等）。本仓库不附带 GO 注释；请从 Ensembl/GO 官方下载并记录版本。

### 2.7 配置文件

复制 `workflow/config/config.yaml`，改路径和阈值。关键项：

```yaml
alignments: "alignments"
tree: "resources/species_tree.nwk"
species_metadata: "resources/species_metadata.tsv"
gene_map: "resources/gene_map.tsv"
gene2go: "resources/gene2go.tsv"

qc:
  min_taxa_fraction: 0.6      # 论文的 ≥60% 115 种
  min_taxa_absrel: 4          # HyPhy 至少 4 个 taxa
  allow_internal_stops: false
  max_invalid_fraction: 0.05

absrel:
  branches: "All"
  genetic_code: "Universal"
  p_threshold: 0.05
  timeout_sec: 7200

go_terms: ["GO:0002376"]

enrichment:
  organism: "hsapiens"
  fdr: 0.05
  min_foreground: 10

regression:
  seed: 20250101
  chains: 4
  iter: 4000
  warmup: 2000
  min_branches: 20
```

## 3. 运行

### 3.1 先 dry-run

```bash
snakemake -s workflow/Snakefile --cores 8 --use-conda -n all
```

`-n` 只打印将执行的 job，不运行。检查：

- 基因数是否等于 `alignments/*.fa` 数量；
- 是否有意外的路径（相对路径错误会导致“No alignments found”）；
- 规则图是否符合预期：`snakemake -s workflow/Snakefile --dag | dot -Tpdf > dag.pdf`。

### 3.2 玩具数据（只验证管线）

```bash
python examples/toy/make_toy_data.py   # 可选；仓库已包含生成结果

snakemake -s workflow/Snakefile \
  --configfile workflow/config/config.toy.yaml \
  --cores 2 --use-conda all

snakemake -s workflow/Snakefile \
  --configfile workflow/config/config.toy.yaml \
  --cores 2 --use-conda downstream
```

玩具配置：

- 6 个物种，3 个基因；
- 分支长度是 substitutions/site；
- `min_branches: 8`，短链；
- `min_foreground: 1`，否则 2 个基因会被富集跳过；
- g:Profiler 可能对 2 个基因返回空结果，这是正常的，不代表管线坏了。

**玩具结果不能用于任何生物学解释**。它的用途是确认：

1. Snakemake 能找到基因并并行；
2. 树能被标注、剪枝、映射；
3. HyPhy 能跑并产出可解析 JSON；
4. 聚合、GO 计数、回归、富集脚本能读写正确列。

### 3.3 真实数据：分阶段

```bash
# 阶段 1：树标注 + 逐基因 aBSREL + 汇总（离线）
snakemake -s workflow/Snakefile --cores 16 --use-conda all

# 阶段 2：GO 分支计数 + g:Profiler + 回归（需要 R/网络）
snakemake -s workflow/Snakefile --cores 16 --use-conda downstream
```

建议先跑 10–20 个基因做基准：

```bash
# 临时把 alignments 指向一个子集目录，或复制 10 个文件到 alignments_bench/
```

不要直接改 `results/` 里的文件来“跳过”失败基因；正确做法是看 `errors.tsv` 和 `status.json`，修数据或调 QC。

### 3.4 SLURM

`workflow/profiles/slurm/config.yaml` 是示例：

```yaml
executor: slurm
jobs: 200
latency-wait: 120
retries: 1
default-resources:
  - mem_mb=4000
  - runtime=120
```

运行：

```bash
snakemake -s workflow/Snakefile \
  --profile workflow/profiles/slurm \
  --use-conda all
```

需要安装 Snakemake 的 SLURM executor plugin（`snakemake-executor-plugin-slurm`），并在集群上配置好 `sbatch`。不同集群的 partition/account 参数不同，按需在 profile 里加 `slurm-partition`、`slurm-account` 等。

## 4. 输出解读

### 4.1 逐基因输出

`results/genes/{gene}.branch.tsv`：

| 列 | 含义 |
|---|---|
| `gene_id` | 基因 |
| `branch_id` | 全树分支 ID（内部节点 `n00001`…，末端为物种名） |
| `node_label` | aBSREL 输入树中的节点标签 |
| `tested` | 该分支是否被 aBSREL 检验 |
| `corrected_p` | aBSREL 基因内 BH 校正 P |
| `uncorrected_p` | 未校正 P |
| `lrt` | 似然比统计量 |
| `baseline_omega` | baseline MG94xREV ω |
| `omega_max` | 该分支 rate classes 中最大的 ω |
| `prop_sites_omega_gt1` | ω>1 的位点比例（按 rate class 权重） |
| `selected` | `tested & corrected_p < threshold & omega_max > 1` |
| `mapping_ok` | 是否能映射回全树分支 |

`results/genes/{gene}.status.json`：

- `status`: `ok` / `skipped` / `error`；
- `qc`: taxa 数、codons、gap/N/invalid fraction、internal/terminal stops、variable sites；
- `n_branches_tested/selected`；
- `mapping_failures`；
- `hyphy_version`、时间戳、失败原因。

### 4.2 汇总

`results/aggregate/branch_tests.tsv.gz`：所有基因 × 所有分支的明细（最大表）。

`results/aggregate/gene_summary.tsv.gz`：

- `n_branches_total/tested/selected`；
- `min_corrected_p`、`max_omega`；
- `selected_any`（论文式）；
- `selected_orders`、`selected_bat/non_bat/mixed`；
- `q_gene`、`selected_any_fdr`（**可选敏感性**，不是论文主结果）。

`results/aggregate/branch_summary.tsv`：

- `n_genes_tested`、`n_genes_selected`、`selected_fraction`；
- 分支元数据：`branch_type`、`branch_length`、`order`、`taxon_group`、`is_order_stem`。

`run_summary.json`：状态计数、HyPhy 版本、分支/基因数、映射失败数。

`errors.tsv`：被跳过/失败的基因及原因。**先看这张表**，不要只看“成功”的基因。

### 4.3 回归

`results/regression/model_comparison_main.tsv`：

- `m0`–`m3` 的 WAIC、LOO、Bayes R²、Rhat、ESS、divergences、rank；
- `formula` 列明确写“单/双截距”。

`model_comparison_sensitivity.tsv`：

- `m1_offset`、`m1_branch_type`、`m1_no_order_stem`、`m1_terminal_only`、`m1_exclude_mixed`、`m1_log_bl`；
- 这些模型用的行不同，**不要直接比较 LOO**；看系数方向、CI 是否重叠。

`hypothesis_taxon.tsv`：`taxon_modelbat > 0` 的后验概率和 BF（brms `hypothesis`）。

`branch_observed_expected.tsv`：每个分支的 observed、expected、95% CI、residual；ancestral 分支在 `is_order_stem`/`order` 列中。

`fitted_curve.tsv` + `regression_curve.png`：sqrt(length) 与期望计数的曲线。

`pp_check_m1.png`：后验预测检查；如果观测计数的尾部/零计数与模拟差很多，负二项可能不够。

`diagnostics.txt`：seed、链设置、Rhat/ESS/divergence 摘要、公式标签说明。

### 4.4 富集

- `gprofiler_custom.tsv`：背景 = 该目所有 tested 基因；
- `gprofiler_default.tsv`：背景 = g:Profiler 默认注释域（更接近论文）；
- `foreground_sizes.tsv`：每个目的 foreground/background 大小、返回术语数、状态；
- `sessionInfo.txt`：R/gprofiler2 版本。

比较两份结果：

- 若同一术语在两种背景下都显著，结论更稳；
- 若只在 default 背景显著、custom 背景不显著，说明它可能只是“可检验基因集中常见”；
- 不要直接拿 FDR 数值和论文的 g:SCS 数值逐位比较。

### 4.5 ZOVER

```bash
Rscript R/fisher_zover_optimized.R \
  --input examples/toy/resources/zover_records.tsv \
  --outdir results/zover
```

输出：

- `zover_fisher.tsv`：论文式 Fisher OR/CI/P + BH；
- `zover_gee.tsv`：study-clustered GEE OR/CI/P；
- `zover_bootstrap.tsv`：study 聚类 bootstrap CI；
- `zover_study_summary.tsv`：每个 study 的记录/病例/宿主构成。

如果 GEE 报“not estimable”，先检查是否只有一个 study、或暴露变量没有变异；小数据/完全分离时 GEE 不稳定，应以 Fisher + bootstrap 为主并明确标注。

## 5. 验证清单

真实运行后，按顺序检查：

1. **状态计数**：`run_summary.json` 的 `status_counts`；skipped/error 比例是否异常（例如 >20% 说明 QC/物种名有问题）。
2. **映射失败**：`n_mapping_failures` 应为 0；若有，查 `unmapped_branches.tsv`，不要忽略。
3. **基因数**：`n_genes_in_status` 是否等于 alignment 文件数；`genes` 是否等于成功解析的基因数。
4. **分支数**：每个成功基因的 `n_branches_total` 应接近 228（少 taxa 的基因会少一些）。
5. **selected 比例**：与论文 8,608/17,130 ≈ 50% 比较；差异大不一定是错，可能来自树/阈值/ω 判定。
6. **回归诊断**：Rhat < 1.01，bulk ESS > 400，divergences = 0；否则增加 `iter`、`adapt_delta`，或简化模型。
7. **后验预测**：`pp_check` 是否大致覆盖观测；若零计数过多，考虑零膨胀负二项（本仓库未实现）。
8. **敏感性方向**：offset/branch_type/terminal-only/log 模型中，蝙蝠截距方向是否一致；若不一致，报告中要说明。
9. **富集背景**：custom vs default 的术语重叠；只在一方显著的术语要谨慎解释。
10. **ZOVER**：Fisher 与 GEE 的 OR 方向是否一致；若 GEE CI 跨 1 而 Fisher 极小，说明聚类影响很大。

**不要**在没有跑完整数据的情况下声称复现了论文的 17,130/8,608/42 vs 22。本仓库目前没有这些结果。

## 6. 常见错误

### `No alignments found under 'alignments'`

- 配置路径相对于运行目录；
- 文件必须命名为 `{gene}.fa`，不能是 `.fasta`/`.fa.gz`；
- 目录为空或只有隐藏文件。

### `tree has branches without lengths`

树里有无长度分支。aBSREL 需要 substitutions/site；补长度或换树。不要用 1.0 之类占位值。

### `tree leaves missing from species metadata`

叶名与 `taxon` 不一致（大小写、空格、下划线）。Newick 中下划线在 `preserve_underscores=True` 下保留；确保两边完全一致。

### `only N taxa after intersecting alignment and tree; need >= M`

- alignment 的序列 ID 与树叶名不一致；
- 基因确实缺少太多物种；
- `min_taxa_fraction` 太高。先看 `status.json` 的 `qc`。

### `aBSREL executable not found: hyphy`

- 没有 `--use-conda`；
- conda 环境创建失败；
- 手动环境里 `hyphy` 不在 PATH。

### aBSREL 超时

- 提高 `absrel.timeout_sec`；
- 检查比对是否过长/过多 taxa；
- 单基因重跑：`snakemake ... results/genes/GENE.branch.tsv --forcerun screen_gene`。

### g:Profiler 失败/超时

- 需要网络；公司/集群代理可能拦截；
- 可先跑离线核心 `all`，再单独跑 `Rscript R/enrichment_gprofiler_optimized.R`；
- `gene_map.tsv` 缺失/ID 不是 Ensembl 会导致查询为空；
- g:Profiler 服务版本变化会让结果不可逐位复现，记录 `sessionInfo()` 和查询日期。

### `fewer than N usable branches; refusing to fit`

- 分支表太小（玩具数据用 `min_branches: 8`）；
- 零长度分支被过滤太多；
- `taxon` 只有一种水平（脚本会直接报错）。

### R 包缺失

用 `workflow/envs/r.yaml` 创建环境；或手动安装 `brms`, `loo`, `gprofiler2`, `geepack`, `optparse`, `readr`, `dplyr`, `ggplot2`, `posterior`。版本差异可能导致 `brms` 参数名变化，建议用 conda 环境。

## 7. 单元测试

```bash
# 需要 numpy/pandas；可直接使用 workflow/envs/python.yaml 或 screen.yaml
conda env create -f workflow/envs/python.yaml
conda activate bat1k-post
python -m unittest discover -s tests -v
```

测试只覆盖：

- BH q 值；
- aBSREL rate distribution 解析；
- codon QC；
- aBSREL JSON 解析（需要 dendropy）。

不调用 HyPhy/R/Snakemake/g:Profiler，不会产生分析结果。

## 8. 运行时间与并行

- **aBSREL**：单基因从几秒到几分钟不等；17,130 个基因是主要成本。用 `--cores` 并行；每 job 1 thread（`threads: 1`），内存默认 4 GB。
- **聚合**：几分钟；`branch_tests.tsv.gz` 可能很大（数百万行），确保磁盘足够。
- **富集**：每个目一次网络查询；几分钟。
- **回归**：228 行 × 4 模型 × 4 链 × 4000 iter；笔记本几十分钟到几小时；可用 `--quick` 先验证管线（**不要用于推断**）。
- **玩具**：分钟级。

## 9. 复现论文数字的建议顺序

1. 先只跑 `all`，确认 17,130 个基因的状态、分支数、映射失败为 0。
2. 用 `gene_summary` 统计 `selected_any` 的基因数，与 8,608 比较；若不一致，先查：
   - 树/分支长度；
   - aBSREL 版本与 `--branches`；
   - `selected` 是否要求 `omega_max > 1`（论文可能只要求 corrected P）；
   - 60% 阈值和 alignment 版本（raw vs filtered）。
3. 用 `branch_counts` 生成 GO:0002376 的 228 行计数，与作者 `summary.fittedValues_ed.csv` 的 `count` 比较；若逐分支不一致，先解决选择判定/基因集差异。
4. 再跑回归，比较模型排序和蝙蝠截距方向；不要期望数值逐位一致（先验、Stan 版本、seed 不同）。
5. 最后跑富集，比较 custom/default 与论文的术语排序；P 值不会完全一致（g:SCS vs FDR、数据库版本）。

## 10. 重要提醒

- 本仓库**没有**运行论文数据，因此没有“复现成功”的声明。
- 玩具数据是合成的，不能用于生物学结论。
- 论文的原始数据/代码版权归原作者；使用前请引用原文和原仓库。
- 如果你的目标是发表级复现，建议把本仓库的 `config.yaml`、`sessionInfo()`、`run_summary.json`、`errors.tsv` 和最终模型一起归档，并记录所有数据库版本（Ensembl、GO、g:Profiler）。
