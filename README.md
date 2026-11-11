# IR-critical

研究 SHR、TyG、AIP 与危重患者死亡结局的关联及预测表现。本目录是当前项目唯一入口；代码、汇总结果、图表和本地投稿材料采用英文路径。

## 现在先看

- 本地稿件：[Manuscript](deliverables/01_Manuscript_Clean.docx)、[Supplement](deliverables/03_Supplementary_Materials.docx)、[Reviewer responses](deliverables/04_Response_to_Reviewers.docx)。Word 文件保留 Zotero 动态引用，仅在本地提供。
- [图表](deliverables/figures/)、[结果与表图映射](results/analysis/10_CR_table_result_figure_mapping.csv)。
- [运行顺序](docs/RUN_ORDER.md)、[受限输入字典](docs/INPUT_DICTIONARY.md)。

|路径|用途|主要输入与输出|使用入口|
|---|---|---|---|
|`src/extraction/`|三库 SQL 提取及 NWICU 衍生表|授权数据库 → 受限 CSV|各队列 SQL；运行顺序见 docs|
|`src/analysis/`|正式分析、敏感性分析及绘图入口|受限输入 → 汇总结果 → 图|01–14 编号入口|
|`src/shared/`|统计及预测共用函数|由分析入口调用|04_prediction.R 为预测入口|
|`src/figures/`|共用图形样式|汇总结果 → 图形对象|由 12_build_figures.py 调用|
|`src/validation/`|公开文件与结果一致性检查|Git 文件集及汇总结果 → 检查结果|01_validate_public.py|
|`src/reporting/`|本地文稿、回复和动态引用更新|作者稿、汇总结果 → DOCX|本地保留，不上传|
|`data/restricted/`|授权患者级输入及模型对象|原始提取与派生数据|本地保留，不上传|
|`results/analysis/`|当前回归及统计汇总|模型、风险集、流程、敏感性结果|表图映射 CSV|
|`results/prediction/`|预测性能汇总|ROC、DeLong、交叉验证等|mimiciv 子目录|
|`deliverables/`|当前投稿材料及图形导出|DOCX、PNG/PDF/SVG 和图源 CSV|正文及 figures|
|`reports/study/`|中文方法、结果和解读|当前正式结果 → 中文说明|01_methods_zh.md|
|`reports/reproducibility/`|版本、迁移及公开文件清单|软件记录、原始路径映射、校验文件|PATH_MIGRATION.csv|
|`archive/`|本地历史资料与原工作区|两个可解压 ZIP，保留原始资料|不上传|
|`reports/internal/`|内部审计及私人引用映射|执行证据、原始审稿材料|本地保留，不上传|

## 当前状态与定义

三库共同指标队列分别为 3,056、1,133、1,556 人。院内及 ICU 死亡采用 Fine–Gray；MIMIC 28 天死亡采用 Cox。正式模型协变量完整，不执行 MICE。旧结果文件中的 M5/M20 为兼容标签，实际 `m=0, maxit=0`；原结果文件名保留以确保引用可追溯，执行入口已按运行顺序重排。当前为 6 组主图、10 组补充图和 15 张补充表。

本次仅整理路径与入口，复用已完成的统计结果和作者稿，不重新运行患者级分析。新增路径已通过脚本语法、汇总矩阵及实际绘图检查。

## 运行环境

正式结果使用 R 4.5.3；核心包为 cmprsk、survival、Hmisc、mice 和 pROC，具体版本见 `reports/reproducibility/*session.txt`。mice 是共用函数依赖，正式协变量未插补。Python 依赖列于 `requirements.txt`；提取还需 PostgreSQL 客户端及授权数据库。

## 更新和保留

按 [RUN_ORDER](docs/RUN_ORDER.md) 设置 `IR_REVISION_ROOT` 为本目录，在临时目录生成后替换固定成品。普通更新直接修改源文件；作者手改正文及动态引用应保留。

历史原始资料及旧研究方案集中保存在本地 `archive/01_historical_materials.zip`；原工作区保存在 `archive/02_previous_workspace.zip`。归档保留原内部路径以维持冻结资料，不纳入当前运行入口，也不上传。`reports/reproducibility/PATH_MIGRATION.csv` 记录原公开包到当前路径的对应关系。公开仓库只包含可复现代码、非识别性汇总结果和图形，受限数据及本地投稿工作文件由 `.gitignore` 排除。

数据资源：MIMIC-IV 3.1 (10.13026/kpb9-mt58)、MIMIC-III 1.4 (10.13026/C2XW26)、NWICU 0.1.0 (10.13026/s84w-1829)。受限输入需按各数据使用协议取得授权。

Repository: https://github.com/jialeyang0104/IR-critical
