# 当前研究方法

研究在相同患者中比较 SHR、TyG 和 AIP。MIMIC-IV 为主队列，MIMIC-III 和 NWICU 重复验证死亡关联。各库先从完整 ICU 表选择每位患者首次 ICU stay，再应用源年龄可得、成人、时长和实验室条件；每人仅贡献一次。

入组要求年龄 ≥18 岁、ICU 时长及自 ICU 入科至医院出院均 ≥24 小时。四项暴露成分在入科前 6 小时至入科后 24 小时内取均值，单位为 glucose/TG/HDL-C：mg/dL，HbA1c：%。SHR = glucose/(28.7×HbA1c−46.7)，TyG = ln(TG×glucose/2)，AIP = log10[(TG/88.57)/(HDL-C/38.67)]。first-available 是窗口内最早有效值，closest-to-admission 是距入科时点绝对时间差最小的值，可在入科前或后。

院内死亡为主要结局，首次 ICU 死亡和 28 天死亡为次要结局。院内及 ICU 死亡采用 Fine–Gray 模型，存活出院或出 ICU 为竞争事件；MIMIC 28 天死亡采用 Cox 模型。NWICU 28 天记录死亡单列探索性结果。Model 1 未调整，Model 2 调整年龄、性别，Model 3 再加入 SOFA、Charlson、糖尿病。三库正式协变量完整，不进行 MICE 或 Rubin 合并；兼容文件名中的 M5/M20 不代表实际插补。

MIMIC-III 去标识化年龄 ≥300 岁编码为 90 岁，最终队列 65 例。NWICU 使用适配 SOFA，规则见当前补充方法及重建 SQL。

Table 1 展示填补前的观测数据及 P 值。对每个连续变量按结局分组执行 Shapiro–Wilk 检验；两组均 P≥0.05 时用均数（标准差）和 Welch t 检验，否则用中位数（四分位距）和 Mann–Whitney U 检验。当前 29 个连续变量均进入非参数分支。分类变量采用卡方或 Fisher 检验。缺失情况见 Table S2 及逐变量审计 CSV。

预测分析仅在 MIMIC-IV 内执行，比较各指标加入 SOFA、SAPS II、APS III、OASIS 的增量。采用配对 10×5 折交叉验证，200 次患者 bootstrap 重拟合整个流程；重复抽中的同一患者不跨训练/测试折，辅助评分的填补限于训练折。

当前补充材料为 10 张图、15 张表。敏感性分析保留替代实验室汇总、24 小时 landmark、合并症排除和亚组。54 个亚组效应使用共同队列 SD；27 个交互采用双侧 Wald 检验及 Benjamini–Hochberg 校正，见 Figure S10、Table S15。DCA、校准、Basic Model 及其他已删除分析仅保留聚合计算记录，不列为当前稿件结果。

运行入口见 [运行顺序](00_阅读入口与运行顺序.md)。
