#!/usr/bin/env Rscript

# revision v3预测性能、临床效用和探索性亚组。Basic Model变量集沿用原稿
# 已披露的投票筛选结果，但本轮全部系数和性能均在新提取数据上重新估计。

options(stringsAsFactors = FALSE, warn = 1)
root <- Sys.getenv("IR_REVISION_ROOT", unset = "")
if (!nzchar(root)) stop("必须设置IR_REVISION_ROOT。", call. = FALSE)
root <- normalizePath(root, mustWork = TRUE)
source(file.path(root, "00_总说明与运行顺序/18_ANALYSIS_CORE.R"), chdir = TRUE)
ir_check_required_packages(c("mice", "survival", "pROC"))

out <- file.path(root, "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果")
logdir <- file.path(root, "03_PSA_主分析补充分析/01_MIMICIV/07_运行日志")
dir.create(out, recursive = TRUE, showWarnings = FALSE); dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
input <- ir_input_path("MIMIC-IV")
prepared <- ir_prepare_cohort("MIMIC-IV", input)
common <- ir_common_exposure_sample(prepared)
raw <- read.csv(input, check.names = FALSE)
raw <- raw[common$source_row, , drop = FALSE]

getn <- function(cands) ir_num(ir_get(raw, cands, required = FALSE))
getb <- function(cands) ir_event01(ir_get(raw, cands, required = FALSE))
dat <- data.frame(
  source_row = common$source_row,
  hospital_mortality = common$hospital_outcome,
  mortality_28d = common$outcome_28d,
  icu_mortality = common$icu_outcome,
  age = common$age, gender = common$gender, sofa = common$sofa,
  charlson = common$charlson, diabetes = common$diabetes,
  SHR = common$SHR, TyG = common$TyG, AIP = common$AIP,
  glucose = common$glucose, hba1c = common$hba1c,
  sapsii = getn(c("sapsii", "saps_ii")), apsiii = getn(c("apsiii", "aps_iii")), oasis = getn("oasis"),
  sirs = getn("sirs"), RR = getn(c("resp_rate_mean", "RR")), Tm = getn(c("temperature_mean", "Tm")),
  WBC = getn(c("wbc_mean", "WBC")), Hb = getn(c("hemoglobin_mean", "Hb")), BUN = getn(c("bun_mean", "BUN")),
  Na = getn(c("sodium_mean", "Na")), K = getn(c("potassium_mean", "K")), MI = getb(c("myocardial_infarct", "MI")),
  ventilation = getb("ventilation"), liver_failure = getb("liver_failure"), HR = getn(c("heart_rate_mean", "HR")),
  MBP = getn(c("mbp_mean", "MBP")), PLT = getn(c("platelets_mean", "PLT")), Scr = getn(c("creatinine_mean", "Scr")),
  PT = getn(c("pt_mean", "PT")), ALT = getn(c("alt_mean", "ALT")), AST = getn(c("ast_mean", "AST")),
  GCS = getn(c("gcs_min", "gcs")), stringsAsFactors = FALSE
)

# 预测分析使用训练折内可重现的简单中位数/众数填补。暴露指标本身不填补。
fill_prediction <- function(z, vars, reference = NULL) {
  if (is.null(reference)) reference <- z
  for (v in vars) {
    if (is.factor(z[[v]]) || is.character(z[[v]])) {
      tab <- table(reference[[v]], useNA = "no"); value <- names(tab)[which.max(tab)]
      z[[v]][is.na(z[[v]])] <- value; z[[v]] <- droplevels(factor(z[[v]]))
    } else {
      value <- median(reference[[v]], na.rm = TRUE); z[[v]][!is.finite(z[[v]])] <- value
    }
  }
  z
}

folds_for <- function(y, k, seed) {
  set.seed(seed); fold <- integer(length(y))
  for (v in sort(unique(y))) { ii <- sample(which(y == v)); fold[ii] <- rep(seq_len(k), length.out = length(ii)) }
  fold
}

cv_pair <- function(z, outcome, base, added = character(), repeats = 10, k = 5, seed = 51001L) {
  vars <- unique(c(base, added)); zz <- z[, c(outcome, vars), drop = FALSE]
  y <- ir_event01(zz[[outcome]]); base_auc <- ext_auc <- numeric(repeats)
  for (r in seq_len(repeats)) {
    fold <- folds_for(y, k, seed + r); pb <- pe <- rep(NA_real_, length(y))
    for (f in seq_len(k)) {
      train <- fold != f; test <- !train
      train_dat <- fill_prediction(zz[train, , drop = FALSE], vars)
      test_dat <- fill_prediction(zz[test, , drop = FALSE], vars, reference = zz[train, , drop = FALSE])
      fb <- glm(as.formula(paste(outcome, "~", paste(base, collapse = " + "))), train_dat, family = binomial())
      fe <- glm(as.formula(paste(outcome, "~", paste(vars, collapse = " + "))), train_dat, family = binomial())
      pb[test] <- predict(fb, test_dat, type = "response"); pe[test] <- predict(fe, test_dat, type = "response")
    }
    base_auc[r] <- as.numeric(pROC::roc(y, pb, levels = c(0, 1), direction = "<", quiet = TRUE)$auc)
    ext_auc[r] <- as.numeric(pROC::roc(y, pe, levels = c(0, 1), direction = "<", quiet = TRUE)$auc)
  }
  c(n = length(y), events = sum(y), baseline_auc = mean(base_auc), baseline_sd = sd(base_auc),
    extended_auc = mean(ext_auc), extended_sd = sd(ext_auc), delta_auc = mean(ext_auc - base_auc),
    delta_sd = sd(ext_auc - base_auc), repeats = repeats, folds = k)
}

score_vars <- c(SOFA = "sofa", `SAPS II` = "sapsii", `APS III` = "apsiii", OASIS = "oasis")
score_rows <- list()
for (y in c("hospital_mortality", "mortality_28d", "icu_mortality")) for (sn in names(score_vars)) for (x in c("SHR", "TyG", "AIP")) {
  a <- cv_pair(dat, y, score_vars[[sn]], x, seed = 52000 + match(y, c("hospital_mortality", "mortality_28d", "icu_mortality")) * 100)
  score_rows[[length(score_rows) + 1]] <- data.frame(cohort = "MIMIC-IV", outcome = y, baseline = sn, exposure = x,
    n = a["n"], events = a["events"], baseline_auc = a["baseline_auc"], baseline_auc_sd = a["baseline_sd"],
    extended_auc = a["extended_auc"], extended_auc_sd = a["extended_sd"], delta_auc = a["delta_auc"],
    delta_auc_sd = a["delta_sd"], repeats = a["repeats"], folds = a["folds"])
}
score_perf <- do.call(rbind, score_rows)

component_rows <- list()
components <- list(SHR = "SHR", Glucose = "glucose", HbA1c = "hba1c")
for (y in c("hospital_mortality", "mortality_28d", "icu_mortality")) for (nm in names(components)) {
  a <- cv_pair(dat, y, components[[nm]], character(), seed = 53000 + match(y, c("hospital_mortality", "mortality_28d", "icu_mortality")) * 100)
  component_rows[[length(component_rows) + 1]] <- data.frame(cohort = "MIMIC-IV", outcome = y, model = nm,
    n = a["n"], events = a["events"], cv_auc_mean = a["baseline_auc"], cv_auc_sd = a["baseline_sd"], repeats = 10, folds = 5)
}
component <- do.call(rbind, component_rows)

basic_vars <- list(
  hospital_mortality = c("age","sofa","charlson","sirs","RR","Tm","WBC","Hb","BUN","Na","K","MI","ventilation","liver_failure","HR","MBP","PLT","Scr","PT","ALT","AST","GCS"),
  mortality_28d = c("age","sofa","charlson","sirs","RR","Tm","WBC","Hb","Na","K","MI","ventilation","liver_failure","GCS","HR","MBP","PLT","BUN","Scr","PT","ALT"),
  icu_mortality = c("sofa","charlson","sirs","RR","WBC","Na","ventilation","liver_failure","age","HR","MBP","Tm","PLT","Hb","BUN","Scr","K","PT","ALT","AST","GCS")
)
basic_rows <- list()
for (y in names(basic_vars)) for (x in c("SHR", "TyG", "AIP")) {
  a <- cv_pair(dat, y, basic_vars[[y]], x, seed = 54000 + match(y, names(basic_vars)) * 100)
  basic_rows[[length(basic_rows) + 1]] <- data.frame(cohort = "MIMIC-IV", outcome = y, exposure = x,
    n = a["n"], events = a["events"], baseline_auc = a["baseline_auc"], baseline_auc_sd = a["baseline_sd"],
    extended_auc = a["extended_auc"], extended_auc_sd = a["extended_sd"], delta_auc = a["delta_auc"],
    delta_auc_sd = a["delta_sd"], repeats = 10, folds = 5,
    selection_note = "原稿投票筛选变量集在交叉验证前固定；变量筛选未嵌套于重采样，结果仅作探索性分析。")
}
basic <- do.call(rbind, basic_rows)

reclassification_once <- function(z, outcome, base, exposure, index = seq_len(nrow(z))) {
  zz <- z[index, c(outcome, base, exposure), drop = FALSE]
  zz <- fill_prediction(zz, unique(c(base, exposure))); y <- ir_event01(zz[[outcome]])
  fb <- glm(as.formula(paste(outcome, "~", paste(base, collapse = " + "))), zz, family = binomial())
  fe <- glm(as.formula(paste(outcome, "~", paste(c(base, exposure), collapse = " + "))), zz, family = binomial())
  pb <- predict(fb, type = "response"); pe <- predict(fe, type = "response"); up <- pe > pb; down <- pe < pb
  nri <- (mean(up[y == 1]) - mean(down[y == 1])) + (mean(down[y == 0]) - mean(up[y == 0]))
  idi <- (mean(pe[y == 1]) - mean(pe[y == 0])) - (mean(pb[y == 1]) - mean(pb[y == 0]))
  c(nri_continuous = nri, idi = idi)
}
reclass_rows <- list(); set.seed(56001)
for (y in names(basic_vars)) for (x in c("SHR", "TyG", "AIP")) {
  point <- reclassification_once(dat, y, basic_vars[[y]], x)
  yy <- ir_event01(dat[[y]]); boots <- matrix(NA_real_, 200, 2, dimnames = list(NULL, names(point)))
  for (b in 1:200) {
    idx <- c(sample(which(yy == 1), sum(yy == 1), replace = TRUE),
             sample(which(yy == 0), sum(yy == 0), replace = TRUE))
    boots[b, ] <- tryCatch(reclassification_once(dat, y, basic_vars[[y]], x, idx),
                           error = function(e) c(NA_real_, NA_real_))
  }
  reclass_rows[[length(reclass_rows) + 1]] <- data.frame(
    cohort = "MIMIC-IV", outcome = y, comparison = paste0("Basic+", x, " vs Basic"), exposure = x,
    nri_continuous = point["nri_continuous"], nri_low = quantile(boots[,1], .025, na.rm=TRUE),
    nri_high = quantile(boots[,1], .975, na.rm=TRUE), idi = point["idi"],
    idi_low = quantile(boots[,2], .025, na.rm=TRUE), idi_high = quantile(boots[,2], .975, na.rm=TRUE),
    bootstrap_replicates = 200,
    note = "表观、post-selection探索性重分类；变量选择未嵌套于重采样。")
}
reclassification <- do.call(rbind, reclass_rows)

# M5完成数据的固定5折out-of-fold预测按患者求平均，用于校准、Brier和DCA。
sets <- ir_make_imputations(common, 5, 55001L)
prediction <- calibration <- dca <- list()
for (y in names(ir_outcomes)) for (x in c("None", "SHR", "TyG", "AIP")) {
  z0 <- ir_risk_set(sets[[1]], y, 0); yy <- z0$.event; fold <- folds_for(yy, 5, 55100 + match(y, names(ir_outcomes)))
  pm <- matrix(NA_real_, nrow(z0), length(sets))
  for (i in seq_along(sets)) {
    z <- ir_risk_set(sets[[i]], y, 0)
    form <- as.formula(paste(".event ~", paste(c(if(x=='None')character()else x, ir_model_covariates[["Model 3"]]), collapse = " + ")))
    for (f in 1:5) {
      fit <- glm(form, z[fold != f, ], family = binomial())
      pm[fold == f, i] <- predict(fit, z[fold == f, ], type = "response")
    }
  }
  p <- rowMeans(pm); ro <- pROC::roc(yy, p, levels = c(0,1), direction = "<", quiet = TRUE, ci = TRUE)
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  prediction[[length(prediction) + 1]] <- data.frame(cohort = "MIMIC-IV", outcome = y, exposure = x,
    model = if(x=='None')'Model 3'else paste0('Model 3 + ',x), n = length(yy), events = sum(yy), auc = as.numeric(ro$auc),
    auc_low = as.numeric(ro$ci[1]), auc_high = as.numeric(ro$ci[3]), brier = mean((yy-p)^2),
    calibration_intercept = unname(coef(glm(yy ~ 1 + offset(lp), family = binomial()))[1]),
    calibration_slope = unname(coef(glm(yy ~ lp, family = binomial()))[2]),
    prediction_pooling = "Five-fold out-of-fold probabilities; complete covariates")
  br <- unique(quantile(p, seq(0,1,.1), na.rm = TRUE)); dec <- cut(p, br, include.lowest = TRUE)
  cs <- aggregate(cbind(predicted = p, observed = yy) ~ dec, FUN = mean); cs$n_bin <- as.integer(table(dec)[as.character(cs$dec)])
  cs$cohort <- "MIMIC-IV"; cs$outcome <- y; cs$exposure <- x; calibration[[length(calibration)+1]] <- cs
  for (thr in seq(.05,.30,.01)) {
    pos <- p >= thr; n <- length(yy); tp <- sum(pos & yy==1); fp <- sum(pos & yy==0); prev <- mean(yy)
    dca[[length(dca)+1]] <- data.frame(cohort="MIMIC-IV", outcome=y, exposure=x, threshold=thr,
      net_benefit_model=tp/n-fp/n*thr/(1-thr), net_benefit_treat_all=prev-(1-prev)*thr/(1-thr),
      net_benefit_treat_none=0, n=n, events=sum(yy))
  }
}

# 探索性亚组已统一由 08_CR_竞争风险修订/02_代码/02_CR_RUN_ANALYSES.R 按竞争风险/Cox 导出。
ir_write_csv(component, file.path(out,"08_PSA_M4_AO_NA_SHR_component_cv_auc.csv"))
ir_write_csv(score_perf, file.path(out,"09_PSA_M4_AO_NA_ALL_score_comparison_cv_auc.csv"))
ir_write_csv(do.call(rbind,prediction), file.path(out,"10_PSA_M4_AO_M5P_ALL_calibration_brier_auc.csv"))
ir_write_csv(do.call(rbind,calibration), file.path(out,"11_PSA_M4_AO_M5P_ALL_calibration_deciles.csv"))
ir_write_csv(do.call(rbind,dca), file.path(out,"12_PSA_M4_AO_M5P_ALL_dca.csv"))
ir_write_csv(basic, file.path(out,"14_PSA_M4_AO_NA_ALL_basic_model_cv_auc.csv"))
ir_write_csv(data.frame(outcome=names(basic_vars), variables=vapply(basic_vars,paste,collapse=";",FUN.VALUE=character(1)),
                        selection_method="原稿LASSO、随机森林、XGBoost和Boruta投票后至少两种方法入选；本轮未重新筛选"),
             file.path(out,"15_PSA_M4_AO_NA_ALL_basic_model_variable_sets.csv"))
ir_write_csv(reclassification, file.path(out,"16_PSA_M4_AO_NA_ALL_basic_model_nri_idi.csv"))
writeLines(c("# 预测性能与探索性亚组方法说明","",
  "严重度评分和SHR组成比较采用重复10次5折交叉验证；每一比较的基线模型与扩展模型共用折分配。",
  "Basic Model沿用原稿已披露的投票筛选变量集，但全部系数和性能在revision v3新提取样本中重新估计。由于变量筛选未嵌套于交叉验证，该结果仅作探索性分析。",
  "校准、Brier和DCA基于5套完成数据固定5折out-of-fold预测概率的患者级均值；DCA阈值为0.05–0.30。",
  "亚组覆盖年龄、性别和糖尿病，报告三指标、三结局的per-SD Model 3 HR；未进行交互检验或多重比较校正。"),
  file.path(root,"03_PSA_主分析补充分析/01_MIMICIV/08_中文扩展方法说明.md"))
writeLines(c(paste0("运行时间：",Sys.time()),"状态：PASS",paste0("共同样本N：",nrow(common)),capture.output(sessionInfo())),
           file.path(logdir,"02_PSA_M4_AO_M5P_ALL_prediction_subgroups.log"))
