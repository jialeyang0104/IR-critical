#!/usr/bin/env Rscript
# 从受限输入导出本轮Table 1、纳入/排除比较和与Table 1逐变量对应的缺失率。
# 只加载既有描述统计函数，不运行历史Cox分析入口。
root <- Sys.getenv('IR_REVISION_ROOT')
source(file.path(root, '08_CR_竞争风险修订/02_代码/01_CR_CORE.R'))
env <- new.env(parent = globalenv())
wanted <- c('fmt_p', 'fmt_cont', 'fmt_binary', 'group_p',
            'm4_baseline_source', 'comparison_table')
loaded <- character()
for (ex in parse(file.path(root, '00_总说明与运行顺序/19_RUN_ALL_ANALYSES.R'))) {
  if (is.call(ex) && identical(ex[[1]], as.name('<-')) &&
      is.symbol(ex[[2]]) && as.character(ex[[2]]) %in% wanted) {
    eval(ex, env)
    loaded <- c(loaded, as.character(ex[[2]]))
  }
}
stopifnot(setequal(loaded, wanted))
raw <- env$m4_baseline_source(ir_input_path('MIMIC-IV'))
selected <- complete.cases(raw[, c('Glucose', 'HbA1c', 'TG', 'HDL', 'SHR', 'TyG', 'AIP')])
inc <- raw[selected, , drop = FALSE]
group <- ifelse(inc$hospital_outcome == 1, 2L, 1L)
table1 <- env$comparison_table(inc, group, c('Survivor', 'Non-survivor'))
# Table 1 selects its summary and two-group test from observed groupwise
# normality. Both outcome groups contain fewer than 5,000 patients, so the
# Shapiro-Wilk test uses every observed value without subsampling.
normality <- list()
for (i in which(table1$type == 'continuous')) {
  v <- table1$variable[i]
  values <- lapply(1:2, function(g) {
    a <- inc[[v]][group == g]
    a[is.finite(a)]
  })
  stopifnot(all(lengths(values) <= 5000))
  ps <- vapply(values, function(a) {
    if (length(a) < 3 || length(unique(a)) < 2) return(NA_real_)
    shapiro.test(a)$p.value
  }, numeric(1))
  use_normal <- all(is.finite(ps)) && all(ps >= .05)
  if (use_normal) {
    digits <- if (v == 'AIP') 2L else 1L
    display <- function(a) {
      a <- a[is.finite(a)]
      sprintf(paste0('%.', digits, 'f (%.', digits, 'f)'), mean(a), sd(a))
    }
    table1$overall[i] <- display(inc[[v]])
    table1$group_1[i] <- display(values[[1]])
    table1$group_2[i] <- display(values[[2]])
    table1$p_value[i] <- t.test(values[[1]], values[[2]], var.equal = FALSE)$p.value
    table1$p_display[i] <- env$fmt_p(table1$p_value[i])
  }
  normality[[length(normality) + 1L]] <- data.frame(
    variable = v, label = table1$label[i],
    survivor_n = length(values[[1]]), non_survivor_n = length(values[[2]]),
    normality_test = 'Shapiro-Wilk', survivor_normality_p = ps[1],
    non_survivor_normality_p = ps[2], normality_alpha = .05,
    both_groups_pass = use_normal,
    summary = if (use_normal) 'Mean (SD)' else 'Median (IQR)',
    comparison_test = if (use_normal) 'Welch t-test' else 'Mann-Whitney U',
    comparison_p = table1$p_value[i])
}
normality <- do.call(rbind, normality)
table1$total_n <- nrow(inc)
table1$group_1_n <- sum(group == 1)
table1$group_2_n <- sum(group == 2)
comparison <- env$comparison_table(raw, ifelse(selected, 1L, 2L),
  c('Included', 'Excluded: missing ≥1 exposure laboratory'), extended = TRUE)
comparison$total_n <- nrow(raw)
comparison$group_1_n <- sum(selected)
comparison$group_2_n <- sum(!selected)
handling <- vapply(table1$variable, function(v) {
  if (v %in% c('age', 'gender', 'sofa', 'charlson', 'Diabetes'))
    return('Formal adjustment covariate; no missing values in the MIMIC-IV common-index cohort')
  if (v %in% c('Glucose', 'HbA1c', 'TG', 'HDL', 'SHR', 'TyG', 'AIP'))
    return('Required exposure component or index; not imputed')
  'Observed descriptive or supplementary-analysis variable; not imputed for the formal regression hierarchy'
}, character(1))
missingness <- data.frame(variable = table1$variable, label = table1$label,
  total_n = nrow(inc), nonmissing_n = table1$overall_nonmissing_n,
  missing_n = nrow(inc) - table1$overall_nonmissing_n,
  missing_percent = 100 * (nrow(inc) - table1$overall_nonmissing_n) / nrow(inc),
  handling = unname(handling))
cr_write(table1, '00_CR_M4_table1_observed_baseline.csv')
cr_write(normality, '00_CR_M4_Table1_distribution_audit.csv')
cr_write(comparison, '00_CR_M4_included_excluded_comparison.csv')
cr_write(missingness, '00_CR_M4_Table1_variable_missingness.csv')
cr_write(data.frame(eligible_n = nrow(raw), included_n = sum(selected),
  excluded_n = sum(!selected), included_percent = 100 * mean(selected)),
  '00_CR_M4_included_excluded_counts.csv')
message('描述统计已导出；Table 1与Table S1变量数：', nrow(table1))
