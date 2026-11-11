#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
root <- Sys.getenv("IR_REVISION_ROOT", unset = "")
if (!nzchar(root)) stop("必须设置IR_REVISION_ROOT。", call. = FALSE)
root <- normalizePath(root, mustWork = TRUE)
source(file.path(root, "src/shared/01_analysis_core.R"), chdir = TRUE)
ir_check_required_packages(c("mice", "survival", "pROC", "Hmisc"))

cohort_codes <- c("MIMIC-IV" = "M4", "MIMIC-III" = "M3", "NWICU" = "NW")
cohort_dirs <- c("MIMIC-IV" = "mimiciv", "MIMIC-III" = "mimiciii", "NWICU" = "nwicu")
seed_base <- c("MIMIC-IV" = 31001L, "MIMIC-III" = 32001L, "NWICU" = 33001L)
requested <- toupper(trimws(strsplit(Sys.getenv("IR_ANALYSIS_COHORTS", "MIMIC-IV,MIMIC-III,NWICU"), ",", fixed = TRUE)[[1]]))
cohort_alias <- c("MIMIC-IV" = "MIMIC-IV", "M4" = "MIMIC-IV", "MIMICIV" = "MIMIC-IV",
                  "MIMIC-III" = "MIMIC-III", "M3" = "MIMIC-III", "MIMICIII" = "MIMIC-III",
                  "NWICU" = "NWICU", "NW" = "NWICU")
run_cohorts <- unique(unname(cohort_alias[requested]))
run_cohorts <- run_cohorts[!is.na(run_cohorts)]
if (!length(run_cohorts)) stop("IR_ANALYSIS_COHORTS未包含可识别队列。", call. = FALSE)

result_path <- function(cohort, job) {
  code <- cohort_codes[[cohort]]; cdir <- cohort_dirs[[cohort]]
  sa_subdir <- cdir
  switch(job,
    AO_M5 = file.path(root, "src/extraction", cdir, "results/analysis",
                      sprintf("%02d_PA_%s_AO_M5P_ALL_cox_results.csv", match(cohort, names(cohort_codes)), code)),
    AO_M20 = file.path(root, "results/historical/01_AdmissionOrigin_M20P", sa_subdir, "results/analysis",
                       sprintf("%02d_SA_%s_AO_M20P_ALL_cox_results.csv", match(cohort, names(cohort_codes)), code)),
    L24_M5 = file.path(root, "results/historical/02_Landmark24h_M5P", sa_subdir, "results/analysis",
                       sprintf("%02d_SA_%s_L24_M5P_ALL_cox_results.csv", match(cohort, names(cohort_codes)), code)),
    L24_M20 = file.path(root, "results/historical/03_Landmark24h_M20P", sa_subdir, "results/analysis",
                        sprintf("%02d_SA_%s_L24_M20P_ALL_cox_results.csv", match(cohort, names(cohort_codes)), code)),
    PER_SD = file.path(root, "results/prediction", cdir, "results/analysis",
                       sprintf("%02d_PSA_%s_AO_M5P_ALL_per_sd_cox_results.csv", match(cohort, names(cohort_codes)), code))
  )
}

run_core_jobs <- function(d, cohort) {
  jobs <- list(
    AO_M5 = list(m = 5, origin = "AdmissionOrigin", level = "PA", scales = c("per_unit", "tertile", "p_trend"), seed = seed_base[[cohort]]),
    AO_M20 = list(m = 20, origin = "AdmissionOrigin", level = "SA", scales = c("per_unit", "tertile", "p_trend"), seed = seed_base[[cohort]] + 1L),
    L24_M5 = list(m = 5, origin = "Landmark24h", level = "SA", scales = c("per_unit", "tertile", "p_trend"), seed = seed_base[[cohort]] + 2L),
    L24_M20 = list(m = 20, origin = "Landmark24h", level = "SA", scales = c("per_unit", "tertile", "p_trend"), seed = seed_base[[cohort]] + 3L),
    PER_SD = list(m = 5, origin = "AdmissionOrigin", level = "PSA", scales = "per_sd", seed = seed_base[[cohort]])
  )
  for (nm in names(jobs)) {
    j <- jobs[[nm]]
    ans <- ir_run_cox_matrix(d, j$m, j$origin, j$scales, j$level, j$seed)
    output_path <- result_path(cohort, nm)
    ir_write_csv(ans, output_path)
    unit <- dirname(dirname(output_path))
    key <- ans[ans$exposure_scale %in% c("per_unit", "per_sd") & ans$model == "Model 3", , drop = FALSE]
    effect_lines <- if (nrow(key)) vapply(seq_len(nrow(key)), function(i) sprintf(
      "- %s / %s：HR %.2f（95%% CI %.2f–%.2f），P=%s，N=%d，事件=%d。",
      key$exposure[i], key$outcome[i], key$estimate[i], key$conf_low[i], key$conf_high[i],
      ifelse(key$p_value[i] < .001, "<0.001", sprintf("%.3f", key$p_value[i])), key$n[i], key$events[i]), character(1)) else "- 无可报告结果。"
    method_line <- paste0("时间起点：", j$origin, "；插补：m=", j$m,
                          "、maxit=10，每套完成数据分别拟合Cox模型后按Rubin原则合并。Model 1未调整；Model 2调整年龄和性别；Model 3调整年龄、性别、SOFA、Charlson合并症指数和糖尿病。暴露定义实验室检查不插补。")
    writeLines(c(paste0("# ", cohort, "方法说明"), "", method_line), file.path(unit, "04_中文方法说明.md"), useBytes = TRUE)
    writeLines(c(paste0("# ", cohort, "结果说明"), "", effect_lines), file.path(unit, "05_中文结果说明.md"), useBytes = TRUE)
    writeLines(c(paste0("# ", cohort, "结果解读"), "",
                 "结果用于描述所选共同三指标样本中的预后关联，不作因果、治疗触发或临床优越性解释。不同时间起点和插补数的结果应与主分析并列判断稳健性。"),
               file.path(unit, "06_中文结果解读.md"), useBytes = TRUE)
    message(cohort, " ", nm, "：", nrow(ans), "行")
  }
}

pool_rcs_one <- function(d, cohort, exposure, outcome, seed) {
  datasets <- ir_make_imputations(d, 5, seed)
  x0 <- d[[exposure]]
  knots <- as.numeric(quantile(x0, c(.05, .35, .65, .95), na.rm = TRUE))
  grid <- seq(quantile(x0, .01, na.rm = TRUE), quantile(x0, .99, na.rm = TRUE), length.out = 100)
  reference <- median(x0, na.rm = TRUE)
  fits <- lapply(datasets, function(dat) {
    z <- ir_risk_set(dat, outcome, 0)
    basis <- Hmisc::rcspline.eval(z[[exposure]], knots = knots, inclx = TRUE)
    colnames(basis) <- paste0(".rcs", seq_len(ncol(basis)))
    z <- cbind(z, basis)
    survival::coxph(as.formula(paste("survival::Surv(.time,.event) ~",
      paste(c(colnames(basis), ir_model_covariates[["Model 3"]]), collapse = " + "))), z, ties = "efron")
  })
  terms <- grep("^\\.rcs", names(coef(fits[[1]])), value = TRUE)
  qmat <- do.call(cbind, lapply(fits, function(f) coef(f)[terms]))
  qbar <- rowMeans(qmat)
  ubar <- Reduce("+", lapply(fits, function(f) vcov(f)[terms, terms, drop = FALSE])) / length(fits)
  between <- if (length(fits) > 1) cov(t(qmat)) else matrix(0, length(terms), length(terms))
  total <- ubar + (1 + 1 / length(fits)) * between
  bg <- Hmisc::rcspline.eval(grid, knots = knots, inclx = TRUE)
  br <- as.numeric(Hmisc::rcspline.eval(reference, knots = knots, inclx = TRUE))
  delta <- sweep(bg, 2, br, "-")
  eta <- as.numeric(delta %*% qbar); se <- sqrt(rowSums((delta %*% total) * delta))
  wald <- function(beta, variance) tryCatch(as.numeric(crossprod(beta, solve(variance, beta))), error = function(e) NA_real_)
  w_all <- wald(qbar, total)
  w_non <- if (length(qbar) > 1) wald(qbar[-1], total[-1, -1, drop = FALSE]) else NA_real_
  data.frame(
    analysis_level = "PSA", cohort = cohort, time_origin = "AdmissionOrigin", imputation = "M5 Rubin pooled",
    exposure = exposure, outcome = outcome, model = "Model 3", x = grid, reference = reference,
    estimate = exp(eta), conf_low = exp(eta - 1.96 * se), conf_high = exp(eta + 1.96 * se),
    n = fits[[1]]$n, events = fits[[1]]$nevent, knots = paste(signif(knots, 8), collapse = ";"),
    p_overall = pchisq(w_all, length(qbar), lower.tail = FALSE),
    p_nonlinear = pchisq(w_non, length(qbar) - 1, lower.tail = FALSE), stringsAsFactors = FALSE
  )
}

km_roc_one <- function(d, cohort) {
  km <- roc_sum <- roc_coords <- delong <- list()
  for (outcome in names(ir_outcomes)) {
    z <- ir_risk_set(d, outcome, 0); rocs <- list()
    for (exposure in c("SHR", "TyG", "AIP")) {
      cuts <- as.numeric(quantile(d[[exposure]], c(0, 1/3, 2/3, 1), na.rm = TRUE)); cuts[1] <- -Inf; cuts[4] <- Inf
      z$.tertile <- cut(z[[exposure]], cuts, labels = c("T1", "T2", "T3"), include.lowest = TRUE)
      sf <- survival::survfit(survival::Surv(.time, .event) ~ .tertile, z)
      ss <- summary(sf)
      km[[length(km) + 1]] <- data.frame(cohort = cohort, exposure = exposure, outcome = outcome,
        tertile = sub(".tertile=", "", ss$strata, fixed = TRUE), time = ss$time,
        n_risk = ss$n.risk, n_event = ss$n.event, survival = ss$surv,
        conf_low = ss$lower, conf_high = ss$upper)
      # 固定“指标越高=风险越高”的方向，避免AIP等弱/反向指标被auto静默翻转，
      # 并保证三条ROC可进行有效的配对DeLong比较。
      ro <- pROC::roc(z$.event, z[[exposure]], levels = c(0, 1), direction = "<", quiet = TRUE, ci = TRUE)
      rocs[[exposure]] <- ro
      roc_sum[[length(roc_sum) + 1]] <- data.frame(cohort = cohort, outcome = outcome, exposure = exposure,
        n = nrow(z), events = sum(z$.event), auc = as.numeric(ro$auc),
        conf_low = as.numeric(ro$ci[1]), conf_high = as.numeric(ro$ci[3]))
      cc <- pROC::coords(ro, x = "all", ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE)
      roc_coords[[length(roc_coords) + 1]] <- data.frame(cohort = cohort, outcome = outcome, exposure = exposure, cc)
    }
    for (pair in combn(c("SHR", "TyG", "AIP"), 2, simplify = FALSE)) {
      tt <- pROC::roc.test(rocs[[pair[1]]], rocs[[pair[2]]], method = "delong", paired = TRUE)
      delong[[length(delong) + 1]] <- data.frame(cohort = cohort, outcome = outcome,
        exposure_1 = pair[1], exposure_2 = pair[2], method = "paired DeLong", p_value = tt$p.value)
    }
  }
  list(km = do.call(rbind, km), roc = do.call(rbind, roc_sum), coords = do.call(rbind, roc_coords), delong = do.call(rbind, delong))
}

assign_folds <- function(y, k) {
  fold <- integer(length(y))
  for (v in sort(unique(y))) { ii <- sample(which(y == v)); fold[ii] <- rep(seq_len(k), length.out = length(ii)) }
  fold
}

prediction_metrics <- function(y, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  ro <- pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
  lp <- qlogis(p)
  c(auc = as.numeric(ro$auc), brier = mean((y - p)^2),
    calibration_intercept = unname(coef(glm(y ~ 1 + offset(lp), family = binomial()))[1]),
    calibration_slope = unname(coef(glm(y ~ lp, family = binomial()))[2]))
}

cv_performance <- function(d, cohort, repeats = 20, folds = 5, seed = 34001L) {
  d <- ir_common_exposure_sample(d)
  datasets <- ir_make_imputations(d, 5, seed)
  rows <- list()
  for (outcome in names(ir_outcomes)) for (exposure in c("SHR", "TyG", "AIP")) for (model in names(ir_model_covariates)) {
    z0 <- ir_risk_set(datasets[[1]], outcome, 0); y <- z0$.event
    vars <- c(exposure, ir_model_covariates[[model]])
    metrics <- matrix(NA_real_, repeats, 4, dimnames = list(NULL, c("auc", "brier", "calibration_intercept", "calibration_slope")))
    for (r in seq_len(repeats)) {
      set.seed(seed + r + 100 * match(exposure, c("SHR", "TyG", "AIP")) + 1000 * match(outcome, names(ir_outcomes)))
      fold <- assign_folds(y, folds); p_imp <- matrix(NA_real_, nrow(z0), length(datasets))
      for (i in seq_along(datasets)) {
        z <- ir_risk_set(datasets[[i]], outcome, 0); z$.event <- y
        form <- as.formula(paste(".event ~", paste(vars, collapse = " + ")))
        for (f in seq_len(folds)) {
          fit <- suppressWarnings(glm(form, z[fold != f, , drop = FALSE], family = binomial()))
          p_imp[fold == f, i] <- suppressWarnings(predict(fit, z[fold == f, , drop = FALSE], type = "response"))
        }
      }
      metrics[r, ] <- prediction_metrics(y, rowMeans(p_imp))
    }
    rows[[length(rows) + 1]] <- data.frame(
      analysis_level = "PSA", cohort = cohort, validation_method = "5-fold CV repeated 20 times",
      exposure = exposure, outcome = outcome, model = model, n = length(y), events = sum(y),
      m = 5, repeats = repeats, folds = folds,
      auc_mean = mean(metrics[, "auc"]), auc_sd = sd(metrics[, "auc"]),
      auc_p025 = quantile(metrics[, "auc"], .025), auc_p975 = quantile(metrics[, "auc"], .975),
      brier_mean = mean(metrics[, "brier"]), brier_sd = sd(metrics[, "brier"]),
      calibration_intercept_mean = mean(metrics[, "calibration_intercept"]),
      calibration_slope_mean = mean(metrics[, "calibration_slope"]),
      pooling_note = "每个重复内对5套完成数据的out-of-fold预测概率取患者级均值；AUC不使用Rubin合并。",
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "<0.001", sprintf("%.3f", p)))
fmt_cont <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("")
  sprintf("%.1f (%.1f–%.1f)", median(x), quantile(x, .25), quantile(x, .75))
}
fmt_binary <- function(x) {
  x <- ir_event01(x); den <- sum(!is.na(x)); num <- sum(x == 1, na.rm = TRUE)
  if (!den) return("")
  sprintf("%d (%.1f%%)", num, 100 * num / den)
}
group_p <- function(x, g, type) {
  ok <- !is.na(x) & !is.na(g); x <- x[ok]; g <- g[ok]
  if (length(unique(g)) < 2 || !length(x)) return(NA_real_)
  if (type == "continuous") return(tryCatch(wilcox.test(x ~ g, exact = FALSE)$p.value, error = function(e) NA_real_))
  tab <- table(ir_event01(x), g)
  if (any(dim(tab) < 2)) return(NA_real_)
  tryCatch(if (any(tab < 5)) fisher.test(tab)$p.value else suppressWarnings(chisq.test(tab, correct = FALSE)$p.value),
           error = function(e) NA_real_)
}

m4_baseline_source <- function(input) {
  raw <- read.csv(input, check.names = FALSE)
  getn <- function(cands, lo = -Inf, hi = Inf) ir_num(ir_get(raw, cands, required = FALSE), lo, hi)
  getb <- function(cands) ir_event01(ir_get(raw, cands, required = FALSE))
  intime <- ir_time(ir_get(raw, c("icu_intime", "intime")))
  outtime <- ir_time(ir_get(raw, c("icu_outtime", "outtime")))
  dischtime <- ir_time(ir_get(raw, c("dischtime", "hospital_dischtime")))
  eligible <- is.finite(as.numeric(intime)) & is.finite(as.numeric(outtime)) & is.finite(as.numeric(dischtime)) &
    as.numeric(difftime(outtime, intime, units = "hours")) >= 24 &
    as.numeric(difftime(dischtime, intime, units = "hours")) >= 24
  raw <- raw[eligible, , drop = FALSE]
  getn <- function(cands, lo = -Inf, hi = Inf) ir_num(ir_get(raw, cands, required = FALSE), lo, hi)
  getb <- function(cands) ir_event01(ir_get(raw, cands, required = FALSE))
  glucose <- getn(c("glucose", "first_day_glucose_mean", "glucose_mean"), 10, 2000)
  hba1c <- getn(c("hba1c", "first_day_hba1c_mean", "hba1c_mean"), 2, 25)
  tg <- getn(c("tg", "triglycerides", "first_day_triglycerides_mean"), 1, 10000)
  hdl <- getn(c("hdl", "first_day_hdl_mean"), 1, 500)
  idx <- ir_derive_indices(glucose, hba1c, tg, hdl)
  x <- data.frame(
    hospital_outcome = getb(c("hospital_outcome", "hospital_mortality")),
    outcome_28d = getb(c("outcome_28d", "mortality_28d")),
    icu_outcome = getb(c("icu_outcome", "icu_mortality")),
    gender = getb(c("sex_male")),
    age = getn(c("admission_age", "age"), 18, 120),
    los_hospital = getn(c("los_hospital", "hospital_los_days", "hospital_time_days"), 0, Inf),
    los_icu = getn(c("los_icu", "icu_los_days", "icu_time_days"), 0, Inf),
    sofa = getn(c("sofa_score", "sofa"), 0, 30), sapsii = getn(c("sapsii", "saps_ii")),
    apsiii = getn(c("apsiii", "aps_iii")), oasis = getn(c("oasis")),
    gcs = getn(c("gcs_min", "gcs"), 0, 15),
    charlson = getn(c("charlson_comorbidity_index", "charlson", "charlson_score"), 0, 40),
    sirs = getn(c("sirs"), 0, 4), HR = getn(c("heart_rate_mean", "HR", "hr_mean")),
    RR = getn(c("resp_rate_mean", "RR", "rr_mean")), MBP = getn(c("mbp_mean", "MBP")),
    Tm = getn(c("temperature_mean", "Tm")), WBC = getn(c("wbc_mean", "WBC")),
    PLT = getn(c("platelets_mean", "PLT")), Hb = getn(c("hemoglobin_mean", "Hb")),
    BUN = getn(c("bun_mean", "BUN")), Scr = getn(c("creatinine_mean", "Scr")),
    Na = getn(c("sodium_mean", "Na")), K = getn(c("potassium_mean", "K")),
    INR = getn(c("inr_mean", "INR")), PT = getn(c("pt_mean", "PT")),
    ALT = getn(c("alt_mean", "ALT")), AST = getn(c("ast_mean", "AST")),
    TG = tg, Glucose = glucose, HbA1c = hba1c, HDL = hdl,
    Glucocorticoids = getb(c("glucocorticoids_1d", "glucocorticoids")),
    Vasopressors = getb(c("vasopressor", "vasopressors")), Ventilation = getb(c("ventilation", "mechanical_ventilation")),
    Hypertension = getb(c("hypertension")), MI = getb(c("myocardial_infarct", "MI")),
    CHF = getb(c("congestive_heart_failure", "CHF")), CPD = getb(c("chronic_pulmonary_disease", "CPD")),
    Diabetes = getb(c("diabetes", "Diabetes")), liver_failure = getb(c("liver_failure")),
    renal_failure = getb(c("renal_failure")), longterm_insulin = getb(c("longterm_insulin", "long_term_insulin")),
    longterm_antidiabetics = getb(c("longterm_antidiabetics", "long_term_antidiabetics")),
    TyG = idx$TyG, SHR = idx$SHR, AIP = idx$AIP, stringsAsFactors = FALSE)
  sex_raw <- toupper(trimws(as.character(ir_get(raw, c("gender", "sex"), required = FALSE))))
  x$gender[is.na(x$gender)] <- as.integer(sex_raw[is.na(x$gender)] %in% c("M", "MALE", "1"))
  x
}

comparison_table <- function(d, group, group_labels, extended = FALSE) {
  labels <- c(gender = "Male sex, n (%)", age = "Age (years)", los_hospital = "Hospital length of stay (days)",
    los_icu = "ICU length of stay (days)", sofa = "SOFA", gcs = "GCS", charlson = "Charlson comorbidity index",
    sirs = "SIRS", HR = "Heart rate (beats/min)", RR = "Respiratory rate (breaths/min)", MBP = "Mean blood pressure (mmHg)", Tm = "Temperature (°C)",
    WBC = "WBC (10^9/L)", PLT = "PLT (10^9/L)", Hb = "Hb (g/dL)", BUN = "BUN (mg/dL)",
    Scr = "Scr (mg/dL)", Na = "Na (mmol/L)", K = "K (mmol/L)", INR = "INR", PT = "PT (s)",
    ALT = "ALT (IU/L)", AST = "AST (IU/L)", TG = "TG (mg/dL)", Glucose = "Glucose (mg/dL)",
    HbA1c = "HbA1c (%)", HDL = "HDL-C (mg/dL)", Glucocorticoids = "Glucocorticoids, n (%)",
    Vasopressors = "Vasopressors, n (%)", Ventilation = "Mechanical ventilation, n (%)",
    Hypertension = "Hypertension, n (%)", MI = "MI, n (%)", CHF = "CHF, n (%)", CPD = "CPD, n (%)",
    Diabetes = "Diabetes, n (%)", liver_failure = "Liver failure, n (%)", renal_failure = "Renal failure, n (%)",
    longterm_insulin = "Long-term use of insulin, n (%)", longterm_antidiabetics = "Long-term use of antidiabetics, n (%)",
    TyG = "TyG", SHR = "SHR", AIP = "AIP")
  if (extended) labels <- c(labels,
    sapsii = "SAPS II", apsiii = "APS III", oasis = "OASIS",
    hospital_outcome = "In-hospital mortality, n (%)", outcome_28d = "28-day mortality, n (%)",
    icu_outcome = "ICU mortality, n (%)")
  binary <- c("gender", "Glucocorticoids", "Vasopressors", "Ventilation", "Hypertension", "MI", "CHF", "CPD",
              "Diabetes", "liver_failure", "renal_failure", "longterm_insulin", "longterm_antidiabetics",
              "hospital_outcome", "outcome_28d", "icu_outcome")
  rows <- lapply(names(labels), function(v) {
    type <- if (v %in% binary) "binary" else "continuous"; fun <- if (type == "binary") fmt_binary else fmt_cont
    if(v=='AIP') fun <- function(a) {q<-quantile(a,c(.5,.25,.75),na.rm=TRUE);q[abs(q)<5e-3]<-0;sprintf('%.2f (%.2f–%.2f)',q[1],q[2],q[3])}
    vals <- lapply(seq_along(group_labels), function(i) fun(d[[v]][group == i]))
    p <- group_p(d[[v]], group, type)
    data.frame(variable = v, label = labels[[v]], type = type,
      overall = fun(d[[v]]), group_1 = vals[[1]], group_2 = vals[[2]], p_value = p, p_display = fmt_p(p),
      overall_nonmissing_n = sum(!is.na(d[[v]])), group_1_nonmissing_n = sum(group == 1 & !is.na(d[[v]])),
      group_2_nonmissing_n = sum(group == 2 & !is.na(d[[v]])), group_1_label = group_labels[1], group_2_label = group_labels[2],
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

run_m4_baseline <- function(input) {
  source_raw <- read.csv(input, check.names = FALSE)
  raw <- m4_baseline_source(input)
  selected <- complete.cases(raw[, c("Glucose", "HbA1c", "TG", "HDL")]) &
    complete.cases(raw[, c("SHR", "TyG", "AIP")])
  inc <- raw[selected, , drop = FALSE]
  group_death <- ifelse(inc$hospital_outcome == 1, 2L, 1L)
  table1 <- comparison_table(inc, group_death, c("Survivor", "Non-survivor"))
  table1$total_n <- nrow(inc); table1$group_1_n <- sum(group_death == 1); table1$group_2_n <- sum(group_death == 2)
  inc_exc <- comparison_table(raw, ifelse(selected, 1L, 2L), c("Included", "Excluded: missing ≥1 exposure laboratory"), extended = TRUE)
  inc_exc$total_n <- nrow(raw); inc_exc$group_1_n <- sum(selected); inc_exc$group_2_n <- sum(!selected)
  out <- file.path(root, "results/prediction/mimiciv")
  ir_write_csv(table1, file.path(out, "01_PSA_M4_AO_M5P_ALL_table1_complete_baseline.csv"))
  ir_write_csv(inc_exc, file.path(out, "07_PSA_M4_AO_NA_ALL_included_excluded_complete_comparison.csv"))
  ir_write_csv(data.frame(eligible_n = nrow(raw), included_n = sum(selected), excluded_n = sum(!selected),
                          included_percent = 100 * mean(selected)),
               file.path(out, "08_PSA_M4_AO_NA_ALL_included_excluded_counts.csv"))
  id_time <- c("subject_id","hadm_id","stay_id","admittime","dischtime","deathtime","icu_intime","icu_outtime")
  vars <- setdiff(names(source_raw), id_time)
  miss <- do.call(rbind, lapply(vars, function(v) {
    x <- source_raw[[v]]
    missing <- if (is.character(x)) is.na(x) | trimws(x) == "" else is.na(x)
    method <- if (tolower(v) %in% c("glucose","hba1c","tg","triglycerides","hdl","shr","tyg","aip"))
      "Not imputed (exposure-defining variable)"
    else if (tolower(v) %in% c("admission_age","age","gender","sex","sofa","sofa_score","charlson","charlson_comorbidity_index","diabetes"))
      "MICE when missing; pmm for numeric and logistic regression for binary"
    else "Not used as an imputed covariate in the formal Model 1–3 Cox hierarchy"
    data.frame(variable = v, total_n = length(x), missing_n = sum(missing),
               missing_percent = 100 * mean(missing), imputation_method = method)
  }))
  ir_write_csv(miss, file.path(out, "03_PSA_M4_AO_M5P_ALL_missingness.csv"))
}

run_supplement <- function(d, cohort) {
  cdir <- cohort_dirs[[cohort]]; code <- cohort_codes[[cohort]]
  out <- file.path(root, "results/prediction", cdir, "results/analysis")
  common <- ir_common_exposure_sample(d)
  rcs <- do.call(rbind, unlist(lapply(c("SHR", "TyG", "AIP"), function(x)
    lapply(names(ir_outcomes), function(y) pool_rcs_one(common, cohort, x, y, seed_base[[cohort]]))), recursive = FALSE))
  sr <- km_roc_one(common, cohort)
  perf <- cv_performance(common, cohort, seed = seed_base[[cohort]] + 10L)
  ir_write_csv(rcs, file.path(out, sprintf("%02d_PSA_%s_AO_M5P_ALL_rcs_source_data.csv", match(cohort, names(cohort_codes)), code)))
  ir_write_csv(sr$km, file.path(out, sprintf("%02d_PSA_%s_AO_NA_ALL_km_source_data.csv", match(cohort, names(cohort_codes)), code)))
  ir_write_csv(sr$roc, file.path(out, sprintf("%02d_PSA_%s_AO_NA_ALL_roc_auc.csv", match(cohort, names(cohort_codes)), code)))
  ir_write_csv(sr$coords, file.path(out, sprintf("%02d_PSA_%s_AO_NA_ALL_roc_coordinates.csv", match(cohort, names(cohort_codes)), code)))
  ir_write_csv(sr$delong, file.path(out, sprintf("%02d_PSA_%s_AO_NA_ALL_delong.csv", match(cohort, names(cohort_codes)), code)))
  ir_write_csv(perf, file.path(out, sprintf("%02d_PSA_%s_AO_M5P_ALL_cv_auc_calibration_brier.csv", match(cohort, names(cohort_codes)), code)))
}

run_cc_qa <- function(d, cohort) {
  common <- ir_common_exposure_sample(d)
  needed <- unique(c(ir_model_covariates[["Model 3"]], unlist(ir_outcomes)))
  cc <- common[complete.cases(common[, needed, drop = FALSE]), , drop = FALSE]
  if (!nrow(cc)) stop("完整病例QA样本为空。", call. = FALSE)
  ans <- ir_run_cox_matrix(cc, 1, "AdmissionOrigin", c("per_unit", "per_sd"), "QA", seed_base[[cohort]] + 20L)
  ans$imputation <- "Complete case; internal QA only"
  ir_write_csv(ans, file.path(root, "05_QA_内部统计审计", sprintf("20_QA_%s_AO_CC_ALL_complete_case_results.csv", cohort_codes[[cohort]])))
}

write_imputation_audit <- function(d, cohort) {
  common <- ir_common_exposure_sample(d)
  vars <- ir_model_covariates[["Model 3"]]
  fixed <- c("SHR", "TyG", "AIP", unlist(ir_outcomes))
  frame <- common[, unique(c(vars, fixed)), drop = FALSE]
  method <- mice::make.method(frame); method[] <- ""
  for (v in vars) if (anyNA(frame[[v]])) method[v] <- if (is.factor(frame[[v]])) "logreg" else "pmm"
  pred <- mice::make.predictorMatrix(frame); diag(pred) <- 0; pred[setdiff(names(frame), vars), ] <- 0
  unit <- file.path(root, "src/extraction", cohort_dirs[[cohort]], "data_dictionary")
  ir_write_csv(data.frame(variable = vars, total_n = nrow(frame), missing_n = vapply(frame[vars], function(x) sum(is.na(x)), numeric(1)),
                          missing_percent = 100 * vapply(frame[vars], function(x) mean(is.na(x)), numeric(1)),
                          method = unname(method[vars])), file.path(unit, paste0("01_PA_", cohort_codes[[cohort]], "_AO_M5P_ALL_covariate_missingness.csv")))
  ir_write_csv(data.frame(target = rep(rownames(pred), each = ncol(pred)), predictor = rep(colnames(pred), times = nrow(pred)),
                          include = as.integer(as.vector(t(pred)))),
               file.path(unit, paste0("02_PA_", cohort_codes[[cohort]], "_AO_M5P_ALL_predictor_matrix.csv")))
  ir_write_csv(data.frame(cohort = cohort, common_sample_n = nrow(frame), exposure_variables_imputed = FALSE,
                          m_primary = 5, m_sensitivity = 20, maxit = 10,
                          seed_m5 = seed_base[[cohort]], seed_m20 = seed_base[[cohort]] + 1L,
                          pooling = "Rubin 1987 scalar pooling on log-effect and model-based variance"),
               file.path(unit, paste0("03_PA_", cohort_codes[[cohort]], "_AO_M5P_ALL_imputation_run_metadata.csv")))
}

audit_rows <- blocks <- list()
for (cohort in run_cohorts) {
  input <- tryCatch(ir_input_path(cohort), error = function(e) e)
  if (inherits(input, "error")) {
    p <- file.path(root, "05_QA_内部统计审计", paste0("19_QA_", cohort_codes[[cohort]], "_ALL_ALL_ALL_revision_v3阻塞.md"))
    ir_write_block(p, paste0(cohort, " revision v3重跑阻塞"), conditionMessage(input))
    blocks[[cohort]] <- conditionMessage(input)
    next
  }
  d <- tryCatch(ir_prepare_cohort(cohort, input), error = function(e) e)
  if (inherits(d, "error")) {
    p <- file.path(root, "05_QA_内部统计审计", paste0("19_QA_", cohort_codes[[cohort]], "_ALL_ALL_ALL_revision_v3阻塞.md"))
    ir_write_block(p, paste0(cohort, " revision v3字段阻塞"), conditionMessage(d))
    blocks[[cohort]] <- conditionMessage(d)
    next
  }
  common <- ir_common_exposure_sample(d)
  audit_rows[[length(audit_rows) + 1]] <- data.frame(
    cohort = cohort, source_file = attr(d, "source_path"), extracted_n = nrow(d), common_exposure_n = nrow(common),
    hospital_risk_n = sum(is.finite(common$hospital_time_days) & common$hospital_time_days > 0),
    landmark_hospital_risk_n = sum(is.finite(common$hospital_time_days) & common$hospital_time_days > 1),
    mortality28_risk_n = sum(is.finite(common$time_28d_days) & common$time_28d_days > 0),
    landmark_mortality28_risk_n = sum(is.finite(common$time_28d_days) & common$time_28d_days > 1),
    icu_risk_n = sum(is.finite(common$icu_time_days) & common$icu_time_days > 0),
    landmark_icu_risk_n = sum(is.finite(common$icu_time_days) & common$icu_time_days > 1),
    model3_complete_n = sum(complete.cases(common[, ir_model_covariates[["Model 3"]], drop = FALSE])),
    stringsAsFactors = FALSE)
  write_imputation_audit(d, cohort)
  run_core_jobs(d, cohort)
  run_supplement(d, cohort)
  run_cc_qa(d, cohort)
  if (cohort == "MIMIC-IV") run_m4_baseline(input)
}

if (length(audit_rows)) ir_write_csv(do.call(rbind, audit_rows),
  file.path(root, "05_QA_内部统计审计/19_QA_ALL_ALL_ALL_ALL_revision_v3输入分母审计.csv"))
writeLines(c(paste0("运行时间：", Sys.time()), paste0("完成队列：", paste(vapply(audit_rows, `[[`, character(1), "cohort"), collapse = "、")),
             if (length(blocks)) paste0("阻塞队列：", paste(names(blocks), collapse = "、")) else "阻塞队列：无",
             capture.output(sessionInfo())),
           file.path(root, "05_QA_内部统计审计/19_QA_ALL_ALL_ALL_ALL_revision_v3运行日志.txt"), useBytes = TRUE)

if (length(blocks)) quit(status = 2)
