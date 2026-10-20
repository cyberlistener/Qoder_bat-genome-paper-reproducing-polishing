# 02 · 代码运行流程与“为什么这样做”

这篇论文的计算部分不是一个单一脚本，而是“基因组组装 → 注释/同源 → 密码子比对 → 系统发育 → 全分支选择筛查 → 富集 → 回归 → 湿实验”的串联。公开仓库 [2023_Bat1Kimmunity](https://github.com/ariadnamorales/2023_Bat1Kimmunity) 只有 README 形式的命令说明；真正可执行的统计代码主要在 [count2branches](https://github.com/lmdavalos/count2branches)。下面把每一段的目的、输入输出和逻辑讲清楚。

## 0. 全局流程

```
10 个新蝙蝠组装 + 105 个公共哺乳动物基因组
        │
        ├─ TOGA：orthology / gene annotation（人类参考）
        │
        ├─ 1:1 orthologs → MACSE 密码子比对 → HmmCleaner
        │
        ├─ RAxML（串联）/ ASTRAL（溯祖）→ 拓扑
        │        └─ treePL + 17 个化石 → 时间树
        │
        ├─ 17,130 个基因 × 228 个分支 → aBSREL（HyPhy）
        │        └─ 每基因内 BH → 8,608 个“至少一个分支入选”的基因
        │
        ├─ 按目 g:Profiler 富集 → 免疫系统过程等
        │
        ├─ 每分支免疫基因计数 + 分支长度 → 负二项回归
        │        └─ 蝙蝠截距、observed − expected、ancestral 分支残差
        │
        ├─ ZOVER：蝙蝠/啮齿类/蝙蝠科的冠状病毒检出比例 → Fisher
        │
        └─ ISG15：结构、二聚化、抗病毒、分泌、去 ISGylation
```

**为什么顺序不能随意换**：

- 选择筛查的分支必须在**同一棵带分支长度的树**上定义，否则不同基因的分支不可比。
- 富集的 foreground 必须来自**同一套选择判定规则**；否则“蝙蝠富集更强”可能只是判定阈值不同。
- 回归的分支长度必须与 aBSREL 使用的树/单位协调；拿时间树长度解释“替换机会”会混淆时间与突变率。
- 湿实验是独立证据，不依赖上述统计，但它是“这些基因是否真的影响病毒”的关键验证。

## 1. 组装、注释、同源与比对

### 1.1 组装

- 9 个物种：PacBio HiFi + Hi-C；*M. spasma*：ONT + Bionano + 10× + Hi-C。
- 目的：获得染色体级、低错误率的参考基因组，使 TOGA 能准确推断完整 ORF，减少“基因为假缺失/假突变”的误判。
- 为什么重要：选择筛查的输入是 codon alignment；如果组装有 frameshift/错拼，密码子比对会错位，aBSREL 可能把错位当成正选择。

### 1.2 TOGA 注释与 orthology

- TOGA 以人类 GENCODE V38/Ensembl 104 为参考，把人类基因投射到每个组装，输出 intact/missing 等分类。
- 目的：得到跨 115 种的 1:1 ortholog 集合，并区分“真缺失”和“测序/组装失败”。
- 为什么用参考引导：跨 115 个物种做 de novo 注释不现实；TOGA 能利用外显子结构和剪接位点，比纯序列比对更适合远缘物种。

### 1.3 密码子比对与清理

- MACSE v2：逐外显子做密码子感知比对，保留阅读框。
- HmmCleaner：识别并清理疑似比对错误区域。
- 目的：aBSREL 需要密码子对齐；错配会产生假替换，尤其是免疫基因这种快速演化、含重复/低复杂度区域的基因。
- 风险：清理过度会删掉真实变异；清理不足会保留错配。论文对感兴趣基因做了人工检查，但没有对全部基因做统一敏感性分析。

### 1.4 系统发育

- 串联树：RAxML 8.1.16，1:1 orthologs 共 30,354,372 bp。
- 溯祖树：ASTRAL 5.5.9。
- 时间树：treePL + 17 个化石校准点。
- 编码区分支长度：IQ-TREE 2.1.3。
- 4D 中性位点分支长度：phyloFit。

**为什么同时要多种树/长度**：

- 拓扑用于定义分支和同源关系；串联树和溯祖树一致，说明拓扑稳健。
- 选择筛查需要 substitutions/site；时间树用于解释“每百万年”，4D 位点用于“中性替换”，编码区用于“蛋白编码替换”。
- 三种长度相关性高但不相同；论文在回归中分别拟合，避免结论依赖某一种长度。

## 2. 选择压力筛查：aBSREL

### 2.1 为什么用 aBSREL

aBSREL 是 HyPhy 的 **adaptive branch-site random effects likelihood** 方法：

- 对每个分支，比较“允许一部分位点 ω>1”的模型与“更简单的 ω 类别”模型；
- 用 AIC 惩罚模型复杂度；
- 输出每个分支的 LRT、uncorrected P、corrected P、rate classes、ω 分布。

它适合这篇论文的问题，因为：

1. 作者不想预设“哪些分支应该受选择”，而是想做**全分支探索**；
2. 正选择可能是**episodic**（只在某些分支、某些位点），不是全树平均；
3. 228 个分支可以逐分支比较，便于后续“分支长度 vs 入选基因数”的回归。

### 2.2 为什么是 17,130 个基因、228 个分支

- 115 个叶子的 rooted binary tree 有 `2×115−2 = 228` 个分支。
- 基因集要求至少在 60% 的物种（约 69 种）中有序列；最终 17,130 个基因通过。
- 每个基因内部，对全部 tested branches 做 BH；输出 corrected P。
- “selected” 的判定通常是：`tested == TRUE & corrected P < 0.05 & ω_max > 1`。不同实现可能略有差异；本仓库把这个规则显式写在 `screen_gene.py`。

**为什么在基因内做 BH，而不是全局**：

- aBSREL 一次运行会测 228 个分支；基因内 BH 控制“这个基因的 228 次检验”的 FDR。
- 但跨 17,130 个基因仍然有大量检验；因此“8,608 个基因”不是基因级 FDR 控制后的数字。论文的用法是**探索性筛查 + 相对比较**，不是“8,608 个基因都确定受选择”。
- 本仓库保留基因内 BH 作为主定义，另加可选的基因级 BH（对每个基因的最小 corrected P 做 BH），只作为敏感性分析，避免偷换论文定义。

### 2.3 原代码如何跑

公开 README 给出的核心命令类似：

```bash
hyphy absrel \
  --alignment gene.codon.fa \
  --tree species_tree.nwk \
  --branches All \
  --output gene.ABSREL.json
```

然后需要从 JSON 中提取每个分支的 `Corrected P-value`、`Rate Distributions` 等。**公开仓库没有提供这一步的解析/汇总代码**，这是复现的最大缺口之一。本仓库的 `screen_gene.py` 补上了：

1. 读 FASTA、树、物种元数据；
2. 取 alignment 与 tree 的交集；
3. 检查长度、内部终止密码子、非法字符、可变位点；
4. 按交集剪枝树；
5. 给剪枝树的内部节点重新命名，并用“后代叶集合 → 全树 MRCA”映射回全树分支；
6. 调 HyPhy；
7. 解析 JSON，输出逐分支 TSV 和状态 JSON。

### 2.4 为什么分支映射要这么麻烦

aBSREL 的输出以**节点标签**为键。如果：

- 内部节点没有标签，或
- 两个内部节点标签相同，或
- 基因缺少某些物种，剪枝后节点被压缩/重建，

那么直接把 JSON 的 key 当成“全树分支 ID”就会错位。错误不会让 HyPhy 报错，但会让后续“哪个分支入选”完全错乱。

本仓库的做法：

- 先给全树每个内部节点一个唯一、确定的标签（`n00001`…），叶节点保留物种名；
- 对每个基因，按 alignment 的实际物种剪枝；
- 对剪枝树的每个节点，取其后代叶集合，在全树中找 MRCA；
- 用 MRCA 的标签作为全树分支 ID；
- 如果两个剪枝节点映射到同一全树分支，直接报错，而不是静默合并。

这样即使某个基因缺少 20% 的物种，剩余分支仍然对应到正确的全树分支。

## 3. 富集分析

### 3.1 论文的逻辑

- 对每个目，取“在该目分支上入选”的基因作为 foreground；
- 用 g:Profiler 查 GO/KEGG/Reactome/WikiPathways 等；
- 用全部人类注释基因作为 background；
- 重点比较 5 个各 20 种的目，报告 corrected P 和 precision；
- 用 10 物种子抽样、替代组装做稳健性。

### 3.2 为什么背景很重要

富集检验问的是：

> 在某个“机会集合”里，foreground 中带术语 T 的比例，是否高于 background 中带 T 的比例？

如果 background 是“所有人类注释基因”，而实际只有 17,130 个基因可能被选为 foreground，那么：

- 一些在可检验基因中很常见的术语，会因为 background 里比例较低而显得“富集”；
- 不同目的可检验基因集不同，P 值不完全可比；
- “免疫系统过程”这种大术语尤其容易受影响。

论文的选择在**描述性**上是合理的（它想回答“在人类注释空间里是否富集”），但不等于“在可检验基因中是否富集”。本仓库同时跑两种背景：

- `custom_tested_in_order`：该目所有有 tested branch 的基因；
- `gprofiler_annotated`：g:Profiler 默认注释域。

两者都输出，方便判断结论是否依赖背景。

### 3.3 原代码的富集片段

公开 README/代码中的 g:Profiler 调用大致是：

```r
gost(
  query = selected_genes,
  organism = "hsapiens",
  significant = TRUE
)
```

问题：

- 没有 `custom_bg`，默认背景是全部注释基因；
- `significant = TRUE` 只返回显著术语，无法做“所有术语”的比较或后续校正；
- 没有记录 g:Profiler 数据库版本/校正方法；
- 没有把“每个目的可检验基因”作为背景。

本仓库的 `R/enrichment_gprofiler_optimized.R` 修正为：

- `significant = FALSE`, `user_threshold = 1` 返回全部术语；
- `correction_method = "fdr"` 显式化；
- 同时跑 custom/default 两种背景；
- 输出 `p_value`、`fold_enrichment`、`precision` 相关列、foreground/background 大小；
- 记录 `sessionInfo()`。

注意：论文可能使用的是 g:Profiler 默认的 g:SCS 校正，而不是 FDR；本仓库用 FDR 是为了跨目统一和可解释，不能声称与论文数值完全一致。

## 4. 分支长度回归：原代码逐行说明

原仓库有三个 R 脚本：

### 4.1 `count2branches_v2.R`：拟合模型

1. **读数据**：`summary.fittedValues_ed.csv`，228 行，每行一个分支；列包括 `label, blma, blsr1, blsr2, count, taxon, taxon2, taxon3, lab`。
   - `count`：该分支上 GO:0002376 的入选基因数；
   - `blma`：时间树长度（Mya）；
   - `blsr1`：4D 中性位点替换数；
   - `blsr2`：编码区替换数；
   - `taxon`：bat/not；
   - `taxon2`：bat/ancestral/not（用于画图颜色）；
   - `taxon3`：lg/sm（用于点大小；不是模型变量）。
2. **设置因子水平**：`taxon` 的 reference 设为 `not`，这样 `taxonbat` 的系数就是“蝙蝠相对非蝙蝠的截距差”。
3. **对每种分支长度拟合 4 个模型**：
   - `count ~ sqrt(bl)`：单截距单斜率；
   - `count ~ sqrt(bl) + taxon`：双截距单斜率（ANCOVA）；
   - `count ~ sqrt(bl) : taxon`：单截距双斜率；
   - `count ~ sqrt(bl) * taxon`：双截距双斜率。
4. **family = negbinomial(link='log')**：计数过离散；log link 保证均值正。
5. **对“最佳”模型再拟合 Poisson 和“不变换”的负二项**，作为分布/变换的替代。
6. **`add_criterion(..., "waic")`**：为每个模型加 WAIC。
7. **`save.image("count2branches.RData")`**：保存所有模型。

原代码的关键问题：

- 注释把第 3 个模型（`:`）写成“two intercepts, two slopes”，但它只有一个截距；第 4 个（`*`）才是双截距双斜率。这会影响补充表的解读。
- 没有 `set.seed()`；brms 默认随机数不同，后验抽样的 Monte Carlo 误差不可复现。
- 没有 Rhat/ESS/divergence 检查；没有先验说明；没有后验预测检查。
- 每个 block 末尾都 `save.image`，后一个 block 会覆盖前一个；如果中途失败，之前的结果可能丢失。
- 没有 `sessionInfo()`，也没有记录 brms/Stan 版本。
- `sqrt()` 与论文 Box–Cox λ≈0.05 不一致；原代码没有解释。

### 4.2 `print_count2branches_v1.R`：输出比较

1. `load("count2branches.RData")`；
2. 对每种长度，`loo_compare(fit0, fit.p, fit1, fit2, fit3, fit4, criterion="waic")`；
3. 打最佳模型 `summary(fit2)` 和 `bayes_R2(fit2)`；
4. 跨长度比较 `fit2, fit2i, fit2m, fit4m`；
5. 用 `predict(fit, robust=TRUE)` 生成每分支预测，导出 `predicted.csv`。

**为什么用 WAIC**：WAIC 是贝叶斯模型比较指标，惩罚有效参数数；论文还报告 R²。问题是：

- `loo_compare(..., criterion="waic")` 是合法用法，但 `loo_compare` 更常用于 LOO；两种指标不应混着解释。
- 不同模型如果用了不同的行（例如缺失值处理不同），WAIC 不可直接比较；原代码没有检查。
- `predict(robust=TRUE)` 对负二项模型返回的是响应尺度还是链接尺度，需要明确；不同版本 brms 行为可能变化。

### 4.3 `plot_count2branches_v2.R`：画图

1. `load("count2branches.RData")`；
2. 对每种长度，用 `data_grid()` 生成预测网格；
3. `add_epred_draws(fit)` 得到后验期望；
4. `stat_lineribbon(.width=0.95)` 画 95% 区间；
5. `geom_point()` 画原始分支，颜色 `taxon2`，大小 `taxon3`；
6. `geom_text()` 标注 `lab` 中的重点分支；
7. `ggsave()` 输出 PDF。

原代码的问题：

- `blsr2` 的 x 轴标签写成 “Branch Lengths (coalescent units)”；`blsr2` 是 IQ-TREE 编码区替换数，单位是 substitutions/site，不是 coalescent units。ASTRAL 的 coalescent units 是另一回事。
- `data_grid()` 里同时放 `taxon2, taxon3, lab`，会生成大量组合；虽然 `add_epred_draws` 只用模型变量，但网格规模不必要地大。
- `scale_color_discrete(type=..., breaks=...)` 给了 4 个颜色但只有 3 个 break；可能有未使用的颜色/图例问题。
- `geom_text(nudge_x=0.01)` 对 blsr1/blsr2 合适，对 blma 用 3；这反映作者知道单位不同，但标签仍然写错。
- 图注/图例没有说明 `taxon3` 的大小编码是什么（lg/sm 不是模型变量）。

## 5. ZOVER 分析

论文：

1. 从 ZOVER 取 metagenomic studies；
2. 统一蝙蝠分类；
3. 构造 2×2 表（蝙蝠/啮齿类 × 冠状病毒/其他病毒）；
4. 做 two-tailed Fisher exact test；
5. 对 Rhinolophidae、Hipposideridae 分别与其他蝙蝠比较。

**为什么用 Fisher**：2×2、样本量大、期望频数低时，Fisher 比卡方更合适。

**逻辑问题**：

- 每条病毒记录作为独立观测；同一 study/物种/样本的记录不独立；
- 检测方法、采样地点、研究目的不同，比例不可直接比较；
- P 值极小（5×10⁻⁵⁰）对聚类非常敏感；
- “检出比例”不是“感染率”，也不是“宿主易感性”。

本仓库的 `R/fisher_zover_optimized.R`：

- 保留论文式 Fisher；
- 加 `geepack::geeglm`，以 `study_id` 为 cluster，exchangeable working correlation；
- 加 study-level cluster bootstrap，给 OR 的 percentile CI；
- 对 3 个主要比较做 BH；
- 输出 study-level 汇总，便于检查是否有单个研究驱动结果。

## 6. 本仓库的 Snakemake 流程

### 6.1 规则依赖

```
label_tree
   ├─ tree.labeled.nwk
   └─ branch_metadata.tsv
        │
screen_gene (每个基因一个 job)
   ├─ results/genes/{gene}.branch.tsv
   └─ results/genes/{gene}.status.json
        │
aggregate
   ├─ branch_tests.tsv.gz
   ├─ gene_summary.tsv.gz
   ├─ branch_summary.tsv
   └─ run_summary.json
        │
branch_counts (每个 GO term)
   └─ branch_counts_GO_xxxxxxx.tsv
        │
regression
   └─ model_comparison_main.tsv, ... diagnostics.txt
```

`enrichment` 直接读 `branch_tests.tsv.gz` 和 `gene_map.tsv`；它不需要等 `branch_counts`。

### 6.2 每个脚本的输入输出

| 脚本 | 输入 | 输出 | 关键点 |
|---|---|---|---|
| `make_branch_metadata.py` | 树、物种元数据 | 带唯一节点标签的树、分支元数据 | 检查分支长度、叶名、重复 taxon；标注 order/taxon_group/is_order_stem |
| `screen_gene.py` | 单基因 FASTA、带标签树、物种元数据 | 分支 TSV、状态 JSON | QC、剪枝、MRCA 映射、HyPhy、JSON 解析 |
| `aggregate_absrel.py` | 所有分支 TSV、状态 JSON、分支元数据 | 分支/基因/运行汇总 | 合并元数据、可选基因级 BH、错误/未映射报告 |
| `make_branch_counts.py` | 分支测试表、gene2go | 每分支 GO 计数 | 去重基因、`n_genes_tested` 作为 offset |
| `enrichment_gprofiler_optimized.R` | 分支测试表、gene map | custom/default 富集表 | 两种背景、FDR、foreground 大小 |
| `count2branches_optimized.R` | 分支计数表 | 模型比较、诊断、曲线、observed/expected | 修正公式标签、先验、seed、Rhat/ESS、offset 敏感性 |
| `fisher_zover_optimized.R` | record-level ZOVER 表 | Fisher/GEE/bootstrap | study 聚类 |

### 6.3 为什么用“每基因一个 job”

- 17,130 个基因，单个 aBSREL 可能几十秒到几分钟；整段脚本跑一旦中途失败，重跑成本高。
- Snakemake 以输出文件为 checkpoint；失败基因重试一次，仍失败则写状态 JSON，不阻塞其他基因。
- 每个基因的临时目录独立，避免并发写冲突。
- `status.json` 记录 QC、taxa 数、tested/selected 数、HyPhy 版本、失败原因，便于审计“为什么这个基因没有结果”。

### 6.4 为什么状态 JSON 也要作为 output

如果只用 `branch.tsv` 作为输出，失败的基因不会生成文件，Snakemake 会一直重试或报错；把 `status.json` 也声明为输出，可以让“失败”成为一个**已完成的、可审计的状态**。聚合脚本会统计 `ok/skipped/error`，并把错误写到 `errors.tsv`。

## 7. 运行时间与资源

- 17,130 个基因 × aBSREL：主要成本。单基因取决于比对长度、taxa 数、CPU；用 `--cores` 并行。
- 内存：每基因 aBSREL 通常 < 4 GB；config 默认 4000 MB。
- 富集：g:Profiler 需要网络；每个目一次查询，通常几秒到几分钟。
- 回归：228 个分支 × 4 个主模型 × 4 条链 × 4000 iter；普通笔记本几十分钟到几小时，取决于 CPU。
- 玩具数据：6 taxa、3 基因，分钟级。

## 8. 一句话总结

原论文的统计流程是“**全分支选择筛查 → 按目富集 → 分支长度负二项回归 → 功能实验**”；原公开代码只覆盖了最后两步中的回归和画图，前面的 aBSREL 解析/汇总、分支映射、富集参数和诊断都缺失。本仓库把这些环节拆成可恢复、可审计的规则，并把统计上最容易被质疑的地方（基因级 FDR、富集背景、分支非独立性、检验机会 offset、树长单位、模型标签）做成显式的敏感性分析。
