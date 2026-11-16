# TBK 年代扰动实验：结果记录

2026 年 11 月 11 日完成全部 **123 个输入、492 次方法调用和 984 条方向记录**。年代扰动改变了南北端点的统计方向证据；解释这些变化时，准确年代下的方法局限同样重要。

最明确的转换来自 GC：双向证据由准确年代的 3/3 降至最高误差档的 9/30，但固定 VAR(1) 的残差在全部输入上都拒绝白噪声。LKIF 在准确年代已呈单向，低误差档反而全部呈双向。PCMCI+ 和经过 lag screen 的 CCM 在准确年代均有未决证据，因此四种方法并未共同呈现同一种“双向关系被年代误差破坏”的过程。

## 比较设置与图件

采用原 notebook 零起算 cell 12 的四格模型：南端扩散系数 0.005、原有过程噪声 SD 0.005，无额外测量误差。年步长积分 39 kyr，保留模型时刻 9000–39000 年的 30 kyr；外部输送事件沿原有 1500 年周期延长。T1 为 South/X，T4 为 North/Y；每 20 年同步采样，在固定 29–1 kyr BP 网格分析，共 1401 点。

R1/R2/R3 对应种子 20261009、20261301、20261302。每条历史有自己的南北不确定年代中心；各非零档均对两条记录采用 C 型误差，每层 proposal SD 为 20、50、100、200 年，相关长度为 1000 年，每档十个成员围绕中心生成。准确年代分母为 3；各非零档为 **3 条历史 × 10 个成员 = 30**，不是三十条独立历史。温度值和因果设置在年代比较中保持不变。

PCMCI+/CPU GPDC 与 LKIF 最大滞后为 5 步，GC 阶数固定为 1。CCM 使用 E=3、tau=-1、exclusion=10、L=300，原始技能平均十个库，使用 100 个 swap surrogates；lag screen 使用 100 个配对库和 61 个位移（±600 年）。详细设计见[方案](04_tbk_chronology_plan.md)。

依次阅读 [00 模型与记录](../figures/04_check_tbk_chronology/00_model_and_records.png)、[01 方向证据](../figures/04_check_tbk_chronology/01_evidence_vs_uncertainty.png)、[02 配对分数变化](../figures/04_check_tbk_chronology/02_paired_score_changes.png)、[03 CCM 原始与滞后证据](../figures/04_check_tbk_chronology/03_ccm_evidence_and_timing.png)。每图有 PDF，说明见[英文 captions](figure_review/09_tbk_chronology.md)。

## 准确年代已经存在的差异

| 方法 | 三条准确年代历史的结果 |
|---|---|
| PCMCI+ | 三条均保留两个滞后方向，同时具有 native contemporaneous o-o，完整状态均为未决。 |
| LKIF | 三条均只保留 North→South。 |
| GC | 三条均保留双向，无未决标记。 |
| CCM | 原始 Tp=0 检验三条均为双向显著；lag screen 后三条均未决。R1/R3 仍保留 North→South，R2 两方向均未保留；South→North 三条均因峰符号不稳定而未决。 |

参考状态不是物理真值的替代。准确年代的两方向 CCM rho 均很高，范围为 **0.99883–0.99984**，但高重构技能未带来同样明确的时间方向。South→North 的库峰平均有 76% 位于零位移；此处汇总只作描述，90% 筛选仍对每个输入分别执行。

## 双向证据随误差变化

下表为**两个方向均保留且无未决证据**的输入数，分母包含全部请求输入。

| 方法 | Accurate | 20 年 | 50 年 | 100 年 | 200 年 |
|---|---:|---:|---:|---:|---:|
| PCMCI+ | 0/3 | 0/30 | 0/30 | 0/30 | 0/30 |
| LKIF | 0/3 | 30/30 | 17/30 | 2/30 | 1/30 |
| GC | 3/3 | 30/30 | 30/30 | 22/30 | 9/30 |
| CCM，筛选后 | 0/3 | 0/30 | 0/30 | 0/30 | 0/30 |

**PCMCI+ 的零值来自持续未决，不能读成没有方向证据。** 全部 123 个输入均保留 contemporaneous o-o；South→North 和 North→South 各保留 **122/123**。200 年档每方向为 29/30，其余档位全数保留。每条历史各非零档的两方向 q 均为 0.9–1.0，但未决证据同时存在。只报告高 q 会遗漏这一状态。

**LKIF 并非从准确年代双向开始单调退化。** 三条准确年代历史均为仅 North→South；20 年档变为 30/30 双向。50、100、200 年档的“均未保留”分别为 3/30、15/30、24/30，其余包括双向和两个单向状态。低档位增加了 South→North 的显著保留，不能据此认定年代误差提高了因果正确性。

**GC 可以计算严格定义的 both-to-one 转换。** 三条准确年代历史均为无未决的双向，因此各非零档的条件分母为 3 × 10 = 30：

| 每层 SD | 仍为双向 | 仅 South→North | 仅 North→South | 均未保留 | Both→one |
|---|---:|---:|---:|---:|---:|
| 20 年 | 30 | 0 | 0 | 0 | 0/30 |
| 50 年 | 30 | 0 | 0 | 0 | 0/30 |
| 100 年 | 22 | 0 | 7 | 1 | 7/30 |
| 200 年 | 9 | 1 | 7 | 13 | 8/30 |

GC 的双向保留数按 R1/R2/R3，在 100 年档为 10/10、6/10、6/10，在 200 年档为 5/10、2/10、2/10。其他三种方法的准确年代均没有符合完整状态定义的双向历史，因此它们的同类 both-to-one 比例为**不适用**，不是 0%。个体状态与基线转换见 [case_states.csv](../data/processed/04_check_tbk_chronology/case_states.csv) 和 [baseline_transitions.csv](../data/processed/04_check_tbk_chronology/baseline_transitions.csv)。

## GC 的适用性诊断

全部 **123/123** 个 VAR(1) 拟合满足软件的稳定性诊断，但全部 **123/123** 个残差白噪声检验均有 p<0.05。准确年代三条历史的 p 为约 2.0×10⁻¹⁴¹ 至 2.7×10⁻¹³⁵；所有扰动输入的软件结果为 0，不把这一浮点数解释为数学上的严格零概率。该诊断使用十阶残差检验。

因此，GC 的转换说明固定 VAR(1) 流程对年代敏感；残差相关限制了其 F 检验和物理方向解释。我们保留全部结果，没有按成员调阶、筛除基线或事后重拟合。增加历史阶数或处理规定强迫应作为单独的适用性诊断。完整字段位于各 [GC native archives](../data/processed/04_check_tbk_chronology/native/)。

## CCM：原始显著性与 lag screen

下表每个单元格为 **South→North / North→South** 的输入数。

| 每层 SD | 每方向请求数 | 原始 Tp=0 显著 | 筛选后方向保留 | 原始显著且 lag 未决 |
|---|---:|---:|---:|---:|
| Accurate | 3 | 3 / 3 | 0 / 2 | 3 / 1 |
| 20 年 | 30 | 30 / 30 | 1 / 12 | 29 / 17 |
| 50 年 | 30 | 27 / 24 | 9 / 3 | 16 / 20 |
| 100 年 | 30 | 11 / 12 | 0 / 1 | 10 / 10 |
| 200 年 | 30 | 3 / 3 | 0 / 0 | 2 / 3 |

某方向保留可与另一方向未决共存，边际计数不能当作互斥单向状态计数。50 年档仅有 2/30 输入达到完整“仅 South→North”状态；CCM 的完整未决状态在四个非零档依次为 30/30、26/30、12/30、4/30，“均未保留”依次为 0/30、2/30、18/30、26/30。高档位未决数减少主要伴随原始显著性丢失，不表示方向更确定。

两方向原始 rho 中位数由准确年代约 0.9998，降至最高档的 **0.084 和 0.059**。最高档的均值 lag profile 在 South→North 有 10/30、North→South 有 5/30 达到扫描边界；保留这些结果，没有事后扩大窗口。符号不稳定和零位移峰已出现在准确年代，较大误差又增加边界和原始不显著情形。它们反映该 lag screen 在此记录上的局限，不证明物理反馈不存在。[完整方向及峰诊断](../data/processed/04_check_tbk_chronology/ccm_lag_summary.csv)和全部 100×61×2 rho 数组均已保存。

## Native scores 与实际年代误差

全部 **960/960** 个成员—本历史准确年代的方向分数差均成功且有限，包括不显著的结果。最高档的配对差中位数如下；方法分别使用自身单位，不比较跨方法数值大小。

| 方法 | South→North 的 Δscore | North→South 的 Δscore |
|---|---:|---:|
| PCMCI+ | −0.332 | −0.549 |
| LKIF | −0.229 | −0.256 |
| GC | −0.204 | −3.564 |
| CCM，原始 Tp=0 rho | −0.916 | −0.940 |

两方向 CCM 的分数差在所有扰动成员上均为负；分数、原始显著性和筛选状态仍需分开读取。全部分布与每条历史细节见 [paired_score_changes.csv](../data/processed/04_check_tbk_chronology/paired_score_changes.csv) 和 [paired_change_summary.csv](../data/processed/04_check_tbk_chronology/paired_change_summary.csv)。

| 每层 proposal SD | 南北中心相对 RMS 中位数 | 南北成员相对 RMS 中位数 [最小，最大] |
|---|---:|---:|
| 20 年 | 26.8 年 | 39.2 [29.2，48.6] 年 |
| 50 年 | 67.1 年 | 98.0 [73.1，121.6] 年 |
| 100 年 | 134.2 年 | 196.0 [146.1，243.2] 年 |
| 200 年 | 268.3 年 | 392.4 [292.3，487.5] 年 |

中心统计使用三对观测年代；成员统计使用每档三十对成员。范围是描述性范围，不是置信区间。全部 24 个观测年代中心首次提议即被接受；240 个单记录成员来自 243 次提议，三次失序拒绝均在 200 年档，没有覆盖拒绝。接受曲线均保持顺序且覆盖固定网格。高档位的拒绝改变了部分成员与候选 ID 的对应关系，同名成员跨档位不一概视为完全配对；本报告的差值均相对于本历史准确年代参考。[相对误差诊断](../data/processed/04_check_tbk_chronology/relative_error_diagnostics.csv)和[接受计数](../data/processed/04_check_tbk_chronology/chronology_acceptance.csv)保留全部明细。

## 解释边界

1. 模型的相互温度作用来自扩散路径和全局均温恢复项；输送 q 由外部规定，没有温度反馈至 q 的方程。端点双向证据不能直接称为“恢复了 AMOC 反馈”。
2. 三条独立噪声历史共享同一规定的 1500 年强迫周期，没有同时改变强迫形态或物理机制；三十个年代成员不是三十条独立气候历史。
3. 约 2.8–154 年的温度衰减模态、受周期强迫影响的记录自相关和最佳估计滞后是不同量；20 年观测不能解析最快模态，年度方程图不是这一尺度下唯一的统计真图。
4. 共同周期强迫可能支持很高的 CCM 重构技能；当前 surrogate 和 90% lag-sign 规则没有在此模型中校准成反馈检验。拒绝、未决或不显著均不证明物理反馈消失。
5. 相对成员误差 SD=2σ 是未条件化的点位理论量，不是实际接受曲线的时间 RMS。严格跨档配对应核对两条记录的中心和成员候选 ID。
6. 改变积分长度会改变 legacy RandomState 预生成数组布局。同种子下 10 kyr 与 39 kyr 不共享全部前段噪声；本轮配对发生在固定 39 kyr 历史内部。

这轮已给出同一物理历史在年代改变后呈现不同统计方向证据的实例，也揭示了准确年代下的适用性问题。下一轮宜先解释 GC 残差相关、PCMCI+ contemporaneous 未决和 CCM 零位移/周期性峰，再决定是否开展预先规定的观测分辨率或历史长度敏感性。本轮未据结果调整科学参数。

## 验证与运行记录

[原始 notebook](../data/archives/tbk_original_20261110T214805311811Z/One_D_model_TBK.ipynb)与[来源清单](../data/archives/tbk_original_20261110T214805311811Z/source_manifest.json)保留迁移依据，[迁移检验](../data/benchmarks/tbk_migration_checks/20261110T215254006646Z/report.json)及[整体预检](../data/benchmarks/tbk_chronology_checks/20261111T165732295477Z/report.json)通过。首次预检的归档名称问题已修复，未改变计算；[旧报告](../data/benchmarks/tbk_chronology_checks/20261110T215752566281Z/report.json)保留追踪。

[正式运行目录](../data/logs/tbk_chronology_20261111T171153343045Z/)保存源码快照、日志和资源监测。[产物审计](../data/logs/tbk_chronology_20261111T171153343045Z/artifact_audit.json)确认 123 个完整提交、492 个 native archives 和七张汇总表均一致；没有失败、待计算或 CCM lag-scan failure。198 个方法—输入组合有未决证据，计算成功不能等同于方向已判定。[保全审计](../data/logs/tbk_chronology_20261111T171153343045Z/preservation_check.json)确认 36,677 个已有文件的大小、修改时间和哈希未变。[图件审阅](../data/logs/tbk_chronology_20261111T171153343045Z/figure_review.json)核对实际保存数据与输出图；绘图更新未改变计算摘要或结果 CSV。

实验脚本耗时 **1460.47 秒（24.34 分钟）**，包含启动与监测收尾的 supervisor 耗时 1466.14 秒。32 个 workers 运行期间，采样到的最低可用物理内存为 **6.89 GiB**。累计 worker 计时中，PCMCI+/GPDC 为 39054 秒（88.6%），CCM 为 4920 秒（11.2%），LKIF 与 GC 合计约 106 秒。该计时包括各方法的计算与保存，是并行任务累计耗时，不是墙钟时间分解；当前主要耗时位于 PCMCI+/GPDC 流程。

实现见 [tbk_model.py](../toolbox/tbk_model.py)、[tbk_observation.py](../toolbox/tbk_observation.py)、[检查 notebook](../04_0_inspect_tbk_records.ipynb)与[主脚本](../04_check_tbk_chronology.py)。完整结果位于 [data/processed/04_check_tbk_chronology/](../data/processed/04_check_tbk_chronology/)，研究动机见[讨论记录](research_direction_and_tbk_notes.md)。
