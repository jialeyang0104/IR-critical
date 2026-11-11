#!/usr/bin/env Rscript

# revision v3高级敏感性分析。所有正式结果只读取本轮数据库提取生成的
# .restricted_input文件；不读取revision v2结果、旧RDS或项目根目录旧CSV。

options(stringsAsFactors = FALSE, warn = 1)
root <- Sys.getenv("IR_REVISION_ROOT", unset = "")
if (!nzchar(root)) stop("必须设置IR_REVISION_ROOT。", call. = FALSE)
root <- normalizePath(root, mustWork = TRUE)
source(file.path(root, "00_总说明与运行顺序/18_ANALYSIS_CORE.R"), chdir = TRUE)
ir_check_required_packages(c("mice", "survival", "cmprsk"))

cohorts_all <- c("MIMIC-IV", "MIMIC-III", "NWICU")
requested <- toupper(trimws(strsplit(Sys.getenv("IR_ANALYSIS_COHORTS", "MIMIC-IV,MIMIC-III,NWICU"), ",", fixed = TRUE)[[1]]))
aliases <- c("MIMIC-IV"="MIMIC-IV", "M4"="MIMIC-IV", "MIMICIV"="MIMIC-IV",
             "MIMIC-III"="MIMIC-III", "M3"="MIMIC-III", "MIMICIII"="MIMIC-III",
             "NWICU"="NWICU", "NW"="NWICU")
cohorts <- unique(unname(aliases[requested])); cohorts <- cohorts[!is.na(cohorts)]
if (!length(cohorts)) stop("IR_ANALYSIS_COHORTS未包含可识别队列。", call. = FALSE)
codes <- c("MIMIC-IV" = "M4", "MIMIC-III" = "M3", "NWICU" = "NW")
dirs <- c("MIMIC-IV" = "01_MIMICIV", "MIMIC-III" = "02_MIMICIII外部验证", "NWICU" = "03_NWICU外部验证")
seeds <- c("MIMIC-IV" = 41001L, "MIMIC-III" = 42001L, "NWICU" = 43001L)
order_prefix <- function(cohort) sprintf("%02d", match(cohort, cohorts_all))

write_block <- function(section, cohort, reason) {
  p <- file.path(root, section, dirs[[cohort]], "03_汇总结果",
                 paste0("99_BLOCKED_", codes[[cohort]], ".md"))
  ir_write_block(p, paste0(cohort, "高级敏感性分析阻塞"), reason)
}

pool_from_fits <- function(fits, df_n) {
  stopifnot(length(unique(vapply(fits, `[[`, numeric(1), "n"))) == 1,
            length(unique(vapply(fits, `[[`, numeric(1), "events"))) == 1)
  ir_pool_scalar(vapply(fits, `[[`, numeric(1), "q"),
                 vapply(fits, `[[`, numeric(1), "u"), df_n,
                 vapply(fits, `[[`, numeric(1), "k")[1])
}

result_row <- function(cohort, exposure, scale, outcome, model, measure, fit, pooled,
                       analysis_type, sample, note, extra = list()) {
  row <- data.frame(
    analysis_level = "SA", cohort = cohort, time_origin = "AdmissionOrigin",
    imputation = "M5 Rubin pooled", exposure = exposure, exposure_scale = scale,
    outcome = outcome, model = model,
    term = if (scale == "per_sd") "每升高1个标准差" else "每升高1单位",
    effect_measure = measure, n = fit$n, events = fit$events,
    estimate = exp(pooled["beta"]), conf_low = exp(pooled["low"]),
    conf_high = exp(pooled["high"]), p_value = pooled["p"],
    rubin_df = pooled["df"], fraction_missing_information = pooled["fmi"],
    m = 5, maxit = 10, sample_definition = sample, denominator_note = note,
    analysis_type = analysis_type, stringsAsFactors = FALSE
  )
  for (nm in names(extra)) row[[nm]] <- extra[[nm]]
  row
}

run_logistic_finegray <- function(d, cohort) {
  common <- ir_common_exposure_sample(d)
  sets <- ir_make_imputations(common, 5, seeds[[cohort]])
  logistic <- finegray <- list()
  for (x in c("SHR", "TyG", "AIP")) for (y in names(ir_outcomes)) for (model in names(ir_model_covariates)) {
    covars <- ir_model_covariates[[model]]
    fits <- lapply(sets, function(dat) {
      z <- ir_risk_set(dat, y, 0); z$.x <- z[[x]]
      f <- glm(as.formula(paste(".event ~", paste(c(".x", covars), collapse = " + "))),
               z, family = binomial())
      list(q = unname(coef(f)[".x"]), u = unname(vcov(f)[".x", ".x"]),
           n = nobs(f), events = sum(model.response(model.frame(f))), k = length(coef(f)))
    })
    p <- pool_from_fits(fits, fits[[1]]$n)
    logistic[[length(logistic) + 1]] <- result_row(
      cohort, x, "per_unit", y, model, "OR", fits[[1]], p,
      "fixed_outcome_logistic", "SHR、TyG和AIP共同完整样本",
      "固定二分类结局；协变量缺失在5套完成数据中处理后按Rubin原则合并。")
  }
  for (x in c("SHR", "TyG", "AIP")) for (y in c("hospital_mortality", "icu_mortality")) for (model in names(ir_model_covariates)) {
    covars <- ir_model_covariates[[model]]
    fits <- lapply(sets, function(dat) {
      z <- ir_risk_set(dat, y, 0)
      form <- as.formula(paste("~", paste(c(x, covars), collapse = " + ")))
      mm <- model.matrix(form, z)[, -1, drop = FALSE]
      f <- cmprsk::crr(ftime = z$.time, fstatus = ifelse(z$.event == 1, 1L, 2L),
                       cov1 = mm, failcode = 1, cencode = 0)
      idx <- match(x, colnames(mm))
      list(q = unname(f$coef[idx]), u = unname(f$var[idx, idx]), n = nrow(z),
           events = sum(z$.event == 1), k = ncol(mm))
    })
    p <- pool_from_fits(fits, fits[[1]]$events)
    finegray[[length(finegray) + 1]] <- result_row(
      cohort, x, "per_unit", y, model, "subdistribution HR", fits[[1]], p,
      "fine_gray", "SHR、TyG和AIP共同完整样本",
      "死亡为目标事件，存活出院为竞争事件；5套模型按Rubin原则合并。")
  }
  out <- file.path(root, "04_SA_敏感性分析/06_FineGray与固定结局Logistic", dirs[[cohort]], "03_汇总结果")
  ir_write_csv(do.call(rbind, logistic), file.path(out, paste0(order_prefix(cohort), "_SA_", codes[[cohort]], "_AO_M5P_ALL_logistic_results.csv")))
  ir_write_csv(do.call(rbind, finegray), file.path(out, paste0(order_prefix(cohort), "_SA_", codes[[cohort]], "_AO_M5P_ALL_finegray_results.csv")))
}

run_exclusions <- function(d, cohort) {
  flags <- c(exclude_diabetes = "diabetes", exclude_liver_failure = "liver_failure",
             exclude_renal_failure = "renal_failure", exclude_longterm_insulin = "longterm_insulin",
             exclude_longterm_antidiabetics = "longterm_antidiabetics")
  absent <- names(flags)[vapply(flags, function(v) all(is.na(d[[v]])), logical(1))]
  available <- setdiff(names(flags), absent)
  if (!length(available)) stop("全部排除标志均缺失。", call. = FALSE)
  common <- ir_common_exposure_sample(d); sets <- ir_make_imputations(common, 5, seeds[[cohort]] + 1L)
  rows <- list()
  for (exc in available) for (x in c("SHR", "TyG", "AIP")) for (y in names(ir_outcomes)) for (scale in c("per_unit", "per_sd")) {
    flag <- flags[[exc]]
    xmean <- mean(common[[x]], na.rm = TRUE); xsd <- sd(common[[x]], na.rm = TRUE)
    covars <- setdiff(ir_model_covariates[["Model 3"]], if (flag == "diabetes") "diabetes" else character())
    fits <- lapply(sets, function(dat) {
      keep <- ir_event01(dat[[flag]]) == 0
      z <- ir_risk_set(dat[!is.na(keep) & keep, , drop = FALSE], y, 0)
      z$.x <- if (scale == "per_sd") (z[[x]] - xmean) / xsd else z[[x]]
      f <- survival::coxph(as.formula(paste("survival::Surv(.time,.event) ~", paste(c(".x", covars), collapse = " + "))),
                           z, ties = "efron")
      list(q = unname(coef(f)[".x"]), u = unname(vcov(f)[".x", ".x"]), n = f$n,
           events = f$nevent, k = length(coef(f)))
    })
    p <- pool_from_fits(fits, fits[[1]]$events)
    rows[[length(rows) + 1]] <- result_row(
      cohort, x, scale, y, "Model 3", "HR", fits[[1]], p, "comorbidity_exclusion",
      paste0("排除", flag, "阳性患者后的三指标共同样本"),
      if (flag == "diabetes") "排除糖尿病后Model 3不再加入糖尿病常量项。" else "Model 3定义不变。",
      list(exclusion = exc))
  }
  out <- file.path(root, "04_SA_敏感性分析/08_排除特定合并症", dirs[[cohort]], "03_汇总结果")
  ir_write_csv(do.call(rbind, rows), file.path(out, paste0(order_prefix(cohort), "_SA_", codes[[cohort]], "_AO_M5P_ALL_exclusion_results.csv")))
  ir_write_csv(data.frame(cohort = cohort, exclusion = names(flags), source_field = unname(flags),
                          executable = names(flags) %in% available,
                          reason = ifelse(names(flags) %in% available, "字段可用并已执行",
                                          "本轮受限提取缺少该标志；未以其他疾病定义静默替代")),
               file.path(out, paste0(order_prefix(cohort), "_SA_", codes[[cohort]], "_AO_ALL_ALL_exclusion_field_audit.csv")))
}

mode_value <- function(x) names(which.max(table(x, useNA = "no")))[1]
run_selection_weighting <- function(d, cohort) {
  d$selected <- as.integer(complete.cases(d[, c("SHR", "TyG", "AIP")]))
  z <- d
  for (v in c("age", "sofa", "charlson")) z[[v]][is.na(z[[v]])] <- median(z[[v]], na.rm = TRUE)
  for (v in c("gender", "diabetes")) {
    fill <- mode_value(z[[v]]); z[[v]][is.na(z[[v]])] <- fill; z[[v]] <- droplevels(z[[v]])
  }
  psfit <- glm(selected ~ age + gender + sofa + charlson + diabetes, z, family = binomial())
  z$ps <- pmin(pmax(predict(psfit, type = "response"), 1e-4), .9999)
  sel <- z[z$selected == 1, , drop = FALSE]
  sel$w_untrimmed <- mean(z$selected) / sel$ps
  trims <- list(unweighted = rep(1, nrow(sel)), untrimmed = sel$w_untrimmed)
  add_trim <- function(name, probs) {
    b <- quantile(sel$w_untrimmed, probs, na.rm = TRUE)
    trims[[name]] <<- pmin(pmax(sel$w_untrimmed, b[1]), b[2])
  }
  add_trim("p01_p99", c(.01, .99)); add_trim("p025_p975", c(.025, .975)); add_trim("p05_p95", c(.05, .95))
  common <- ir_common_exposure_sample(d); sets <- ir_make_imputations(common, 5, seeds[[cohort]] + 2L)
  weight_by_row <- lapply(trims, function(w) setNames(w, sel$source_row))
  rows <- list()
  for (tn in names(trims)) for (x in c("SHR", "TyG", "AIP")) for (y in names(ir_outcomes)) {
    xm <- mean(common[[x]], na.rm = TRUE); xs <- sd(common[[x]], na.rm = TRUE)
    fits <- lapply(sets, function(dat) {
      zz <- ir_risk_set(dat, y, 0); zz$.x <- (zz[[x]] - xm) / xs
      zz$.w <- unname(weight_by_row[[tn]][as.character(zz$source_row)])
      f <- survival::coxph(survival::Surv(.time,.event) ~ .x + age + gender + sofa + charlson + diabetes,
                           zz, weights = .w, robust = TRUE, ties = "efron")
      list(q = unname(coef(f)[".x"]), u = unname(vcov(f)[".x", ".x"]), n = f$n,
           events = f$nevent, k = length(coef(f)))
    })
    p <- pool_from_fits(fits, fits[[1]]$events)
    w <- trims[[tn]]; ess <- sum(w)^2 / sum(w^2)
    rows[[length(rows) + 1]] <- result_row(
      cohort, x, "per_sd", y, "Model 3", "weighted HR", fits[[1]], p,
      "selection_weighting", "具有全部暴露定义实验室检查的共同样本",
      "基于年龄、性别、SOFA、Charlson指数和糖尿病的稳定逆纳入概率权重。",
      list(trimming = tn, effective_sample_size = ess))
  }
  balance <- lapply(c("age", "sofa", "charlson"), function(v) {
    target <- z[[v]]; observed <- sel[[v]]
    data.frame(variable = v, target_mean = mean(target), selected_mean = mean(observed),
      smd_unweighted = (mean(observed) - mean(target)) / sqrt((var(observed) + var(target)) / 2),
      smd_p01_p99 = (weighted.mean(observed, trims$p01_p99) - mean(target)) / sqrt((var(observed) + var(target)) / 2))
  })
  out <- file.path(root, "04_SA_敏感性分析/07_选择加权", dirs[[cohort]], "03_汇总结果")
  ir_write_csv(do.call(rbind, rows), file.path(out, paste0(order_prefix(cohort), "_SA_", codes[[cohort]], "_AO_M5P_ALL_selection_weighted_results.csv")))
  ir_write_csv(do.call(rbind, balance), file.path(out, paste0(order_prefix(cohort), "_SA_", codes[[cohort]], "_AO_M5P_ALL_covariate_balance.csv")))
  ir_write_csv(data.frame(cohort = cohort, eligible_n = nrow(d), selected_n = sum(d$selected),
                          selected_percent = 100 * mean(d$selected), ps_min = min(sel$ps),
                          ps_median = median(sel$ps), ps_max = max(sel$ps),
                          ess_p01_p99 = sum(trims$p01_p99)^2 / sum(trims$p01_p99^2)),
               file.path(out, paste0(order_prefix(cohort), "_SA_", codes[[cohort]], "_AO_M5P_ALL_weight_summary.csv")))
}

first_existing <- function(d, candidates) {
  hit <- candidates[candidates %in% names(d)]
  if (length(hit)) hit[1] else NA_character_
}
run_lab_timing <- function(input, cohort) {
  raw <- read.csv(input, check.names = FALSE)
  prefix <- c("MIMIC-IV"="mimiciv", "MIMIC-III"="mimiciii", "NWICU"="nwicu")[[cohort]]
  definition_path <- file.path(root, ".restricted_input", paste0(prefix, "_lab_timing_definitions.csv"))
  if (file.exists(definition_path)) {
    definitions <- read.csv(definition_path, check.names = FALSE)
  } else {
    embedded <- list(
      glucose = c("glucose", "glucose_first_value", "glucose_closest_value"),
      hba1c = c("hba1c", "hba1c_first_value", "hba1c_closest_value"),
      tg = c(if ("tg" %in% names(raw)) "tg" else "triglycerides", "tg_first_value", "tg_closest_value"),
      hdl = c("hdl", "hdl_first_value", "hdl_closest_value")
    )
    required_embedded <- unique(c("stay_id", unlist(embedded)))
    if (!all(required_embedded %in% names(raw)))
      stop("缺少标准化实验室定义文件，且主提取未包含全部first/closest/mean字段。", call. = FALSE)
    definitions <- do.call(rbind, lapply(names(embedded), function(a) data.frame(
      stay_id = raw$stay_id, analyte = a,
      first_available_value = raw[[embedded[[a]][2]]],
      closest_to_icu_value = raw[[embedded[[a]][3]]],
      early_window_mean = raw[[embedded[[a]][1]]], stringsAsFactors = FALSE)))
  }
  if (!"stay_id" %in% names(raw) && "icustay_id" %in% names(raw)) raw$stay_id <- raw$icustay_id
  if (!"stay_id" %in% names(definitions) && "icustay_id" %in% names(definitions)) definitions$stay_id <- definitions$icustay_id
  required <- c("stay_id", "analyte", "first_available_value", "closest_to_icu_value", "early_window_mean")
  if (!all(required %in% names(definitions))) stop("实验室取值定义文件缺少字段：", paste(setdiff(required, names(definitions)), collapse = "、"), call. = FALSE)
  if (!"stay_id" %in% names(raw)) stop("主提取缺少stay_id，无法合并实验室时间定义。", call. = FALSE)
  definitions$analyte <- tolower(as.character(definitions$analyte))
  definitions$analyte[definitions$analyte %in% c("triglycerides", "triglyceride")] <- "tg"
  definitions <- definitions[definitions$analyte %in% c("glucose", "hba1c", "tg", "hdl"), , drop = FALSE]
  value_columns <- c(first_available="first_available_value",
                     closest_to_admission="closest_to_icu_value",
                     early_window_mean="early_window_mean")
  sets <- list()
  temp_paths <- character()
  on.exit(unlink(temp_paths), add = TRUE)
  for (method in names(value_columns)) {
    dd <- definitions[, c("stay_id", "analyte", value_columns[[method]]), drop = FALSE]
    names(dd)[3] <- "value"
    if (anyDuplicated(dd[, c("stay_id", "analyte")])) stop("实验室定义文件中stay_id-analyte不唯一。", call. = FALSE)
    wide <- reshape(dd, idvar="stay_id", timevar="analyte", direction="wide")
    names(wide) <- sub("^value\\.", "", names(wide))
    raw$.revision_v3_row_order <- seq_len(nrow(raw))
    tmp <- merge(raw, wide, by="stay_id", all.x=TRUE, sort=FALSE, suffixes=c(".original", ""))
    tmp <- tmp[order(tmp$.revision_v3_row_order), , drop = FALSE]
    for (a in c("glucose", "hba1c", "tg", "hdl")) if (!a %in% names(tmp)) stop(method, "缺少", a, "定义值。", call. = FALSE)
    p <- file.path(root, ".restricted_input", paste0("_", codes[[cohort]], "_", method, ".csv"))
    temp_paths <- c(temp_paths, p)
    write.csv(tmp, p, row.names = FALSE, na = "")
    sets[[method]] <- ir_prepare_cohort(cohort, p)
  }
  common_rows <- Reduce(intersect, lapply(sets, function(z) ir_common_exposure_sample(z)$source_row))
  rows <- list()
  for (i in seq_along(sets)) {
    method <- names(sets)[i]; z <- sets[[i]][sets[[i]]$source_row %in% common_rows, , drop = FALSE]
    ans <- ir_run_cox_matrix(z, 5, "AdmissionOrigin", "per_unit", "SA", seeds[[cohort]] + 10L + i)
    ans$lab_definition <- method
    ans$sample_definition <- "三种实验室取值定义下均可计算SHR、TyG和AIP的共同患者样本"
    rows[[length(rows) + 1]] <- ans
  }
  out <- file.path(root, "04_SA_敏感性分析/05_实验室取值时间", dirs[[cohort]], "03_汇总结果")
  ir_write_csv(do.call(rbind, rows), file.path(out, paste0(order_prefix(cohort), "_SA_", codes[[cohort]], "_AO_M5P_ALL_lab_timing_results.csv")))
  ir_write_csv(data.frame(cohort = cohort, common_across_definitions_n = length(common_rows)),
               file.path(out, paste0(order_prefix(cohort), "_SA_", codes[[cohort]], "_AO_M5P_ALL_lab_timing_counts.csv")))
}

audit <- list()
for (cohort in cohorts) {
  input <- tryCatch(ir_input_path(cohort), error = function(e) e)
  if (inherits(input, "error")) next
  d <- ir_prepare_cohort(cohort, input)
  tasks <- list(
    logistic_finegray = function() run_logistic_finegray(d, cohort),
    exclusions = function() run_exclusions(d, cohort),
    selection_weighting = function() run_selection_weighting(d, cohort),
    lab_timing = function() run_lab_timing(input, cohort)
  )
  for (nm in names(tasks)) {
    status <- tryCatch({ tasks[[nm]](); "PASS" }, error = function(e) paste0("BLOCKED: ", conditionMessage(e)))
    section <- switch(nm, logistic_finegray = "04_SA_敏感性分析/06_FineGray与固定结局Logistic",
                      exclusions = "04_SA_敏感性分析/08_排除特定合并症",
                      selection_weighting = "04_SA_敏感性分析/07_选择加权",
                      lab_timing = "04_SA_敏感性分析/05_实验室取值时间")
    if (startsWith(status, "BLOCKED")) write_block(section, cohort, status)
    unit <- file.path(root, section, dirs[[cohort]])
    method_text <- switch(nm,
      logistic_finegray = "固定二分类结局Logistic模型覆盖三指标、三结局和Model 1–3；院内死亡和ICU死亡另采用Fine–Gray模型并将存活出院视为竞争事件。协变量缺失在5套完成数据中处理，系数与方差按Rubin原则合并。",
      exclusions = "分别评估排除糖尿病、肝衰竭、肾衰竭、长期胰岛素和长期降糖药患者后的Model 3关联，同时报告per-unit和per-SD结果。字段不可用的排除条件仅记录为不可执行，不使用近似疾病标志替代。",
      selection_weighting = "以是否同时具有四项有效暴露定义实验室检查为选择结局，依据年龄、性别、SOFA、Charlson指数和糖尿病估计稳定逆纳入概率权重，并报告多种截尾、有效样本量和协变量平衡。",
      lab_timing = "在同一早期窗口内分别使用首次可用值、最接近ICU入科值和窗口均值重算SHR、TyG和AIP；比较限制在三种定义均可计算三个指标的共同患者样本，采用M5 Rubin合并Model 1–3。")
    interpretation <- switch(nm,
      logistic_finegray = "方向一致可说明主要关联不完全依赖Cox设定；OR、原因特异HR和亚分布HR回答的问题不同，不能直接按数值大小互换。",
      exclusions = "排除分析用于判断结果是否由特定代谢或器官功能亚群驱动；样本量和事件数下降导致区间变宽，不等同于不存在关联。",
      selection_weighting = "该分析只能处理由已测变量解释的选择机制；极端权重、有效样本量下降或残余不平衡均提示仍有选择偏倚。",
      lab_timing = "不同取值定义下方向一致支持早期取值规则的稳健性，但仍不能恢复治疗前暴露或排除输注、营养和胰岛素治疗影响。")
    writeLines(c(paste0("# ", cohort, "方法说明"), "", method_text), file.path(unit, "04_中文方法说明.md"), useBytes = TRUE)
    writeLines(c(paste0("# ", cohort, "结果说明"), "", paste0("本轮revision v3运行状态：", status),
                 "所有数值均来自本轮数据库重提取后的受限输入；未沿用revision v2结果。"),
               file.path(unit, "05_中文结果说明.md"), useBytes = TRUE)
    writeLines(c(paste0("# ", cohort, "结果解读"), "", interpretation),
               file.path(unit, "06_中文结果解读.md"), useBytes = TRUE)
    audit[[length(audit) + 1]] <- data.frame(cohort = cohort, analysis = nm, status = status)
  }
}
ir_write_csv(do.call(rbind, audit), file.path(root, "05_QA_内部统计审计/21_QA_ALL_ALL_ALL_ALL_高级敏感性运行审计.csv"))
if (any(startsWith(vapply(audit, `[[`, character(1), "status"), "BLOCKED"))) quit(status = 2)
