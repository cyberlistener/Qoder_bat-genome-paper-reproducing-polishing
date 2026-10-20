# 03 · 代码与统计问题清单

> 范围：公开的 `2023_Bat1Kimmunity` README、`count2branches` 三个 R 脚本，以及论文 Methods。
> 标记：**确认** = 代码/公式可直接核对；**关注** = 方法学上需要补充说明或敏感性分析；**复现缺口** = 公开材料不足以独立复核。
> 本仓库的修复不改变论文的原始定义；所有“改进”都以**主分析 + 敏感性分析**的方式并列，避免把优化偷换成新的主结论。

## 总表

| ID | 严重度 | 类型 | 位置 | 问题 | 影响 | 本仓库处理 |
|---|---|---|---|---|---|---|
| F1 | 高 | 确认 | `count2branches_v2.R` 模型注释 | `x:taxon` 与 `x*taxon` 的“单/双截距”注释互换 | 补充表/文字可能把单截距模型当双截距解释 | 修正标签；`m2`/`m3` 分开拟合并明确写出公式 |
| F2 | 高 | 复现缺口 | aBSREL 输入树 | 公开材料未给出 aBSREL 树的单位；数据门户只有 dated tree（Mya） | 若用 Mya/coalescent units，aBSREL 的分支长度尺度不对 | 强制 substitutions/site；检查零长/缺失长度；文档要求作者确认 |
| F3 | 高 | 关注 | 选择筛查 | 只在基因内 BH，无基因级 FDR | 8,608 个“selected genes”是上限；富集基因池偏大 | 保留主定义；加 `q_gene`/`selected_any_fdr` 敏感性列 |
| F4 | 高 | 关注 | 回归 | 228 个分支当独立观测；无系统发育相关性 | 截距/残差的 CI 偏窄；outlier 判断过强 | 加 branch_type、去 order-stem、terminal-only、去 mixed 等敏感性；文档标注 |
| F5 | 高 | 关注 | 回归 | 计数未除以该分支可检验基因数 | 分支长度可能混入“检验机会” | 输出 `n_genes_tested`；加 `offset(log(n_genes_tested))` 敏感性模型 |
| F6 | 高 | 关注 | ZOVER | 每条病毒记录独立做 Fisher | P=5×10⁻⁵⁰ 对 study 聚类敏感 | 加 study-clustered GEE + cluster bootstrap + BH |
| F7 | 中 | 确认 | g:Profiler 调用 | `significant=TRUE`、无 `custom_bg`（变量名 `noBG`） | 只返回显著项；背景是全部注释基因，非 17,130 可检验基因 | 两种背景并列；`significant=FALSE`；FDR 显式化 |
| F8 | 中 | 关注 | 富集比较 | “蝙蝠 P 最小”是跨目描述性排序 | P 值受 query size 影响；非正式交互检验 | 输出 precision/fold enrichment；文档建议正式交互/置换检验 |
| F9 | 中 | 确认 | `plot_count2branches_v2.R` | `blsr2` 标成 “coalescent units” | 单位解释错误；读者会误以为 ASTRAL 长度 | 文档纠正；图中写 substitutions/site |
| F10 | 中 | 确认 | Methods + 代码 | Box–Cox λ≈0.05，却用 sqrt（λ=0.5） | 变换选择与报告不一致 | 主模型保留 sqrt；加 log(bl) 敏感性模型 |
| F11 | 中 | 确认 | `count2branches_v2.R` | 无 seed、无 Rhat/ESS/divergence、无 pp_check | 后验不可复现；收敛问题不可见 | 固定 seed、先验、诊断、`pp_check`、`sessionInfo` |
| F12 | 中 | 确认 | `count2branches_v2.R` | 多个 `save.image` 覆盖同一 RData | 中途失败会丢失之前结果 | 每模型 `saveRDS`；输出分开 |
| F13 | 中 | 确认 | `print_count2branches_v1.R` | `loo_compare(..., criterion="waic")` 与 LOO 混用 | 指标解释混乱；不同行数的模型不可比 | 主表同时给 WAIC/LOO；敏感性表标注“不同行，不可头对头比较” |
| F14 | 中 | 关注 | 分支定义 | `taxon2=ancestral` 但模型只用 bat/not；mixed 内部支归为 non-bat | ancestral 非蝙蝠分支被当作普通 non-bat；mixed 支可能稀释信号 | 输出 `taxon_group`；加去 mixed 敏感性；文档明确 |
| F15 | 中 | 复现缺口 | 60% 过滤 | 未报告被排除基因的功能/数量明细 | 免疫基因可能被系统性排除；选择集不可复现 | 状态 JSON 记录 retained/dropped taxa；文档要求报告排除清单 |
| F16 | 中 | 关注 | 湿实验 | pairwise t-test，无统一多重校正/批次模型 | 多构建、多病毒、多读数的假阳性风险 | 未重写（原始数据未公开）；审计标注 |
| F17 | 低 | 确认 | `plot_count2branches_v2.R` | `data_grid` 展开 `taxon2/taxon3/lab`；颜色数/breaks 不一致 | 网格过大；图例可能有多余颜色 | 优化脚本只生成模型需要的网格 |
| F18 | 低 | 确认 | 两个原 R 脚本 | 加载大量未使用包、`rm(list=ls())` | 环境冲突/可读性差 | 只加载必要包；不全局清空 |
| F19 | 高 | 复现缺口 | 公共仓库 | README 命令多、占位符多；无 per-gene 解析/汇总 | 无法一键复现 17,130 基因筛查 | Snakemake + 逐基因状态 + 聚合 |
| F20 | 中 | 复现缺口 | 公共数据 | 只有入选分支的 alignment 包；无逐分支 P/ω 表 | 无法复核 8,608/42 等数字 | 文档列出所需数据；流水线可直接产出这些表 |

---

## F1（确认）回归公式的“单/双截距”注释互换

原代码：

```r
##two intercepts, two slopes
fit3 <- brm(bf(count ~ sqrt(blsr1) : taxon), ...)

##one intercept, two slopes
fit4 <- brm(bf(count ~ sqrt(blsr1) * taxon), ...)
```

R 公式的实际含义：

- `count ~ x:taxon`：保留全局截距，按 taxon 估计两条斜率；**一个截距、两个斜率**。
- `count ~ x*taxon`：展开为 `x + taxon + x:taxon`；**两个截距、两个斜率**。
- `count ~ x + taxon`：两个截距、一个公共斜率（ANCOVA）。

所以 `fit3` 和 `fit4` 的注释正好互换。WAIC 仍然对应实际拟合的公式，不会因为注释错而变；但如果 Supplementary Table 8/9 或正文按注释解释“哪个模型是双截距”，就会把模型含义说反。

本仓库的 `m0`–`m3` 明确写出：

```r
m2 = count ~ sqrt(bl) * taxon_model   # two intercepts, two slopes
m3 = count ~ sqrt(bl):taxon_model     # one intercept, two slopes
```

并用 `model_comparison_main.tsv` 的 `formula` 列把公式和标签一起落盘。

## F2（复现缺口）aBSREL 的树长单位

aBSREL 的 branch-site 模型把分支长度当作**每个位点的期望替换数**（substitutions/site）。论文 Methods 说“使用 ASTRAL topology”，但：

- ASTRAL 的输出分支长度默认是 **coalescent units**；
- 数据门户只提供 `datedTree.50bats.nh` 和 `datedTree.115mammals.nh`，都是 **Mya**；
- 公开 README 的 aBSREL 命令用了一个未说明来源的 `${tree}`；
- 论文没有给 aBSREL 输入树的文件或单位。

如果实际输入是 Mya 或 coalescent units，aBSREL 的似然尺度与“替换数”不一致，分支级 P 值和 ω 可能改变。这不是说结果一定错，而是**必须由作者澄清**，并做“用 substitutions/site 树 vs 原树”的敏感性比较。

本仓库的处理：

- `make_branch_metadata.py` 对缺失分支长度直接报错；
- 文档反复声明“必须是 substitutions/site”；
- `screen_gene.py` 把树剪枝后原样传给 HyPhy，不偷偷重标定；
- 教程要求用户明确记录树来源和单位。

## F3（关注）基因级 FDR 缺失

论文的 BH 是“每个基因内部，对 228 个分支”的校正。跨 17,130 个基因没有再做 FDR。后果：

- “8,608 个基因至少一个分支入选”不是基因级 FDR<0.05 的基因数；
- 富集 foreground 包含大量可能由机会产生的基因，尤其是长基因/多物种/高分化基因；
- 但蝙蝠 vs 非蝙蝠使用同一规则，**相对比较**比绝对数量稳健。

本仓库：

- 主表保留论文定义（`selected` 基于基因内 corrected P）；
- `aggregate_absrel.py` 额外计算 `q_gene = BH(min_corrected_p)` 和 `selected_any_fdr`；
- 文档明确：这是敏感性列，不是论文主结果。

## F4/F5（关注）分支非独立与检验机会

228 个分支来自一棵树，亲缘分支共享历史；普通 NB 回归假设观测独立。更严重的是，计数没有除以“该分支上 GO 基因的可检验数”：

```
count_b ~ NegBin(mu_b, shape)
log(mu_b) = b0 + b1*sqrt(length_b) + b_taxon*I(bat)
```

如果分支 b 有 100 个 GO 基因可检验，分支 c 只有 20 个，b 的 count 自然更大。正确写法之一：

```
log(mu_b) = b0 + b1*sqrt(length_b) + b_taxon*I(bat) + log(n_genes_tested_b)
```

本仓库不把 offset 模型当作主模型（因为论文没用），而是作为 `m1_offset` 敏感性模型；同时输出 `n_genes_tested` 和 `n_genes_tested_all`，让读者自行判断。

分支非独立性目前只能部分处理：

- `m1_branch_type`：区分 terminal/internal；
- `m1_no_order_stem`：去掉目干分支；
- `m1_terminal_only`：只用末端分支；
- `m1_exclude_mixed`：去掉混合后代分支；
- 文档建议：若要正式 inference，应用 phylogenetic mixed model / GEE / phylolm，或至少按目做 cluster-robust SE。

## F6（关注）ZOVER 的伪重复

原分析：

```r
fisher.test(table(host_group, is_coronavirus))
```

每条病毒记录一行，同一 study 的记录被当作独立。若一个 study 贡献 500 条记录，它对 P 值的影响远大于只贡献 5 条的 study。P=5×10⁻⁵⁰ 这种量级尤其需要聚类稳健性。

本仓库 `R/fisher_zover_optimized.R`：

- 保留 Fisher（论文式）；
- `geeglm(coronavirus ~ exposed, id=study_id, family=binomial, corstr="exchangeable")`；
- study-level cluster bootstrap（Haldane 校正 OR，percentile CI）；
- 对 3 个主要比较做 BH；
- 输出 study-level 表，检查是否由单一研究驱动。

注意：GEE 需要 ≥2 个 study 且有 exposure 变异；玩具数据只有 8 个 study，GEE 结果只用于演示。

## F7（确认）g:Profiler 背景与 `significant=TRUE`

原 README：

```r
multi_sig_noBG_tmp_gostres = gost(
  listGenes_groups,
  organism = "hsapiens",
  significant = TRUE
)
```

变量名 `noBG` 直接说明“没有 background”。问题：

1. g:Profiler 默认 background 是全部注释基因，不是 17,130 个可检验基因；
2. `significant=TRUE` 只返回显著项，无法检查未显著项、无法做统一 FDR/比较；
3. 没有记录 `correction_method`（默认 g:SCS）和数据库版本；
4. `gost` 接受 gene list 的 list；返回也是 list；直接 `sapply(..., paste)` 会把结果压成字符串，不利于复算。

本仓库：

- `significant=FALSE`, `user_threshold=1`, `correction_method="fdr"`；
- 每个目跑两次：custom background（该目 tested 基因）和 default；
- 输出标准表 + `foreground_sizes.tsv` + `sessionInfo.txt`；
- 明确说明：FDR 与论文可能的 g:SCS 不同，数值不应逐位对比。

## F8（关注）跨目 P 值比较不是正式检验

论文说蝙蝠的 immune system process P 值最小、precision 最高，并比较 5 个 20 物种目。这是**描述性排序**，不是“蝙蝠的富集程度显著高于啮齿类”的检验。P 值依赖 query size；precision 更稳定，但仍有依赖。

更正式的做法：

- 以“目”为协变量，对“基因是否带术语 T”做 logistic 模型，加入目 × 选择状态的交互；
- 或用置换：在目间打乱“selected”标签，保持每个目的基因数和 GO 注释结构，比较 P/precision 的零分布；
- 或对每个术语拟合“selected ~ order”，报告 order 的整体效应和蝙蝠 vs 其他目的对比。

本仓库目前输出 precision/fold enrichment，便于人工比较；没有实现置换检验，文档把它列为下一步。

## F9（确认）`blsr2` 被标成 coalescent units

原图：

```r
labs(x='Branch Lengths (coalescent units)', ...)
```

`blsr2` 来自 IQ-TREE 编码区替换数，单位是 substitutions/site。coalescent units 是 ASTRAL 的长度。这个标签会误导读者把“编码区替换数”当成“溯祖单位”。本仓库文档和优化图的 x 轴统一写 “substitutions/site”，并要求在图中注明树来源。

## F10（确认）Box–Cox λ≈0.05 与 sqrt 不一致

论文 Methods 写“Box–Cox 最大似然在 λ≈0.05，因此使用平方根变换”。λ=0.05 接近对数变换；平方根对应 λ=0.5。两者不一致，可能是：

- 报告笔误（实际 λ≈0.5）；
- 或作者在 AIC 比较后选择了 sqrt，但文字写成 Box–Cox 结果。

无论哪种，都应报告 log 变换的模型比较。本仓库主模型保留 sqrt（匹配论文），加 `m1_log_bl` 敏感性；如果 log 明显更好，读者可以看到。

## F11/F12（确认）可复现性与诊断

原代码：

- 无 `set.seed()`：brms/Stan 的后验抽样不可逐位复现；
- 无 `prior`：使用默认先验，读者不知道；
- 无 Rhat/ESS/divergence：无法判断链是否收敛；
- 无 `pp_check`：无法检查负二项是否拟合计数分布；
- 多个 `save.image("count2branches.RData")`：后一个 block 覆盖前一个；
- 无 `sessionInfo()`：brms/Stan/R 版本不可追溯。

本仓库：

- `seed`、`chains`、`iter`、`warmup` 全部 CLI/config 化；
- 显式先验：`normal(0,2)` intercept、`normal(0,1)` slope、`gamma(2,0.1)` shape；
- `adapt_delta=0.95`, `max_treedepth=12`；
- `rhat_max`, `ess_bulk_min`, `divergences` 写入 `model_comparison_main.tsv`；
- `pp_check_m1.png`、`branch_observed_expected.tsv`、`fitted_curve.tsv`；
- 每个模型单独 `saveRDS`；
- `diagnostics.txt` + `sessionInfo.txt`。

## F13（确认）WAIC/LOO 混用

原 `print_count2branches_v1.R`：

```r
loo_compare(fit0, fit.p, fit1, fit2, fit3, fit4, criterion = "waic")
```

`loo_compare` 支持 `criterion="waic"`，但这会把 WAIC 的差异当作 LOO 式的比较；两者定义不同，不能混着解释。另外，如果模型用了不同的数据行（缺失值/taxon 水平），比较不成立。

本仓库：

- 主表同时给 `waic`, `looic`, `elpd_loo`, `p_loo`；
- `rank` 明确按 `elpd_loo`；
- 敏感性表标注“不同行，LOO/WAIC 仅描述性”。

## F14（关注）ancestral 与 mixed 分支的归类

原数据：

- `taxon2` 有 `bat/ancestral/not`；
- 模型只用 `taxon`（bat/not）；
- `taxon2=ancestral` 的 15 个分支在模型中都是 `not`（非蝙蝠 ancestral），只有 Chiroptera ancestral 是 `bat`；
- 内部混合分支（如 Boreoeutheria）也被归为 `not`。

这在“蝙蝠 vs 非蝙蝠”的二分类下是合理的，但会带来两个问题：

1. ancestral 非蝙蝠分支和 terminal 非蝙蝠分支共享一个截距，可能掩盖深部分支的差异；
2. mixed 分支的 `taxon=not` 是一种简化，可能稀释信号。

本仓库输出 `taxon_group`（bat/non-bat/mixed），并加 `m1_exclude_mixed` 敏感性；主模型仍按 `taxon` 二分类，保持与论文可比。

## F15（复现缺口）60% 过滤的排除清单

论文只报告 17,130 个通过过滤的基因，没有报告：

- 每个基因在 115 种中的覆盖数分布；
- 被排除基因的 GO/免疫注释比例；
- 排除是否与基因组质量/基因家族相关。

如果免疫基因（尤其快速演化、含重复的基因）更容易低于 60%，那么“免疫基因富集”可能部分来自**可检验性偏差**。本仓库的 `screen_gene.py` 会记录 `alignment_taxa`, `tree_taxa`, `retained_taxa`, `dropped_alignment_taxa`, `tree_taxa_absent_from_alignment`，聚合后可做这个审计；但论文没有公开这层日志。

## F16（关注）湿实验统计

论文图注多处写 two-tailed t-test，n=4–7。潜在问题：

- 多构建、多病毒、多时间点、多读数 → 多重比较；
- 没有报告批次/平板随机效应；
- 没有报告正态性/方差齐性；
- 有些比较是计划好的（planned contrasts），有些是事后选的，应区分。

建议：混合效应模型（batch/plate random intercept）或至少 BH 校正；报告效应量和 CI，而不是只报 P。由于原始数据未公开，本仓库没有重写。

## F17/F18（确认）绘图与代码卫生

低严重度，但会影响可读性：

- `data_grid()` 同时展开 `taxon2/taxon3/lab`，生成不必要的组合；
- `scale_color_discrete(type=4 colors, breaks=3)` 可能有多余颜色；
- `library()` 加载大量未使用包（MCMCglmm、rstanarm、cowplot 等）；
- `rm(list=ls())` 在脚本中会清空用户环境，不适合被 source；
- `save.image()` 把整个全局环境写入 RData，体积大且不可审计。

本仓库的 R 脚本只加载必要包，不 `rm(list=ls())`，输出明确文件。

## F19/F20（复现缺口）公共仓库与数据

- `2023_Bat1Kimmunity` 是 README 命令集合：TOGA、MACSE、HmmCleaner、RAxML、ASTRAL、treePL、HyPhy 的命令都有，但没有逐基因循环、错误处理、JSON 解析、汇总表。
- `count2branches` 只有回归/画图，输入 `summary.fittedValues_ed.csv` 是作者整理后的 228 行表；生成这张表的代码不在仓库里。
- 公共数据只有“至少一个分支入选”的 alignment 压缩包，没有逐分支 P/ω 表；因此无法从原始 alignment 独立复核 8,608、42 vs 22 等数字，除非重跑 17,130 个 aBSREL。

本仓库能补上前两步（逐基因筛查 + 汇总 + 回归/富集），但仍需要作者提供：

1. aBSREL 输入树及单位；
2. 完整基因列表/过滤日志；
3. 逐分支 aBSREL 结果（或完整 JSON）；
4. g:Profiler 参数与原始输出；
5. ZOVER 记录级数据和 study 标识；
6. 湿实验原始读数。

---

## 修复原则：不改主定义，只加“可检验的旁边”

本仓库刻意避免以下“优化”：

- **不**把基因级 FDR 当成新的主筛选阈值；
- **不**用 custom background 替换论文的 default background，而是并列；
- **不**用 offset 模型替换主回归，而是作为敏感性；
- **不**用 GEE 替换 Fisher，而是并列；
- **不**改 aBSREL 的 P 阈值或 ω 判定。

原因是：审计的目标是让读者看到“结论对分析选择有多敏感”，而不是用一套新的偏好把原结论推翻。若敏感性结果与主结果方向一致，结论更可信；若不一致，应报告并解释，而不是只留一个“更好看”的版本。
