#!/usr/bin/env Rscript

# revision v3统一统计分析内核。
# 本文件不读取数据库，也不写公开目录；所有患者级输入和MIDS对象均限制在
# IR_REVISION_ROOT/.restricted_input。调用脚本必须显式设置IR_REVISION_ROOT。

options(stringsAsFactors = FALSE, warn = 1)

ir_revision_root <- function() {
  root <- Sys.getenv("IR_REVISION_ROOT", unset = "")
  if (!nzchar(root)) stop("必须设置IR_REVISION_ROOT。", call. = FALSE)
  normalizePath(root, mustWork = TRUE)
}

ir_event01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) return(as.integer(x))
  if (is.character(x)) return(as.integer(tolower(trimws(x)) %in% c("1", "yes", "true", "dead", "death", "死亡")))
  as.integer(suppressWarnings(as.numeric(x)) == 1)
}

ir_num <- function(x, lower = -Inf, upper = Inf) {
  z <- suppressWarnings(as.numeric(x))
  z[!is.finite(z) | z < lower | z > upper] <- NA_real_
  z
}

ir_time <- function(x) {
  x <- as.character(x)
  x[!nzchar(trimws(x)) | toupper(trimws(x)) %in% c("NA", "NULL", "NAT")] <- NA_character_
  as.POSIXct(x, tz = "UTC", tryFormats = c(
    "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M:%OS", "%Y-%m-%dT%H:%M:%S",
    "%Y-%m-%dT%H:%M:%OS", "%Y-%m-%d"
  ))
}

ir_first_name <- function(d, candidates, required = TRUE, label = candidates[1]) {
  hit <- candidates[candidates %in% names(d)]
  if (length(hit)) return(hit[1])
  if (required) stop("受限输入缺少字段：", label, "（可接受：", paste(candidates, collapse = ", "), "）。", call. = FALSE)
  NA_character_
}

ir_get <- function(d, candidates, required = TRUE, label = candidates[1]) {
  nm <- ir_first_name(d, candidates, required, label)
  if (is.na(nm)) rep(NA, nrow(d)) else d[[nm]]
}

ir_input_path <- function(cohort) {
  root <- ir_revision_root()
  candidates <- switch(cohort,
    "MIMIC-IV" = c("mimiciv_analysis_extract.csv", "mimiciv_analysis.csv", "mimiciv_primary.csv", "mimiciv_restricted_extract.csv"),
    "MIMIC-III" = c("mimiciii_external_validation.csv", "mimiciii_analysis_extract.csv", "mimiciii_analysis.csv"),
    "NWICU" = c("nwicu_external_validation.csv", "nwicu_analysis_extract.csv", "nwicu_analysis.csv")
  )
  patient_file <- paste0(c('MIMIC-IV'='mimiciv','MIMIC-III'='mimiciii','NWICU'='nwicu')[[cohort]], '_patient_level.csv')
  # Current analyses require the audited full-database first-stay derivative.
  candidates <- patient_file
  paths <- file.path(root, ".restricted_input", candidates)
  exists_nonempty <- file.exists(paths) & !is.na(file.info(paths)$size) & file.info(paths)$size > 0
  hit <- paths[exists_nonempty]
  if (!length(hit)) stop("缺少", cohort, "标准化受限提取：", paste(paths, collapse = "；"), call. = FALSE)
  hit[1]
}

ir_derive_indices <- function(glucose_mg_dl, hba1c_pct, tg_mg_dl, hdl_mg_dl) {
  denom <- 28.7 * hba1c_pct - 46.7
  data.frame(
    SHR = ifelse(is.finite(denom) & denom > 0, glucose_mg_dl / denom, NA_real_),
    TyG = ifelse(glucose_mg_dl > 0 & tg_mg_dl > 0, log(tg_mg_dl * glucose_mg_dl / 2), NA_real_),
    # TG和HDL-C原始单位均为mg/dL，先分别换算为mmol/L后计算AIP。
    AIP = ifelse(tg_mg_dl > 0 & hdl_mg_dl > 0,
                 log10((tg_mg_dl / 88.57) / (hdl_mg_dl / 38.67)), NA_real_)
  )
}

ir_duration_from_icu <- function(start, stop, event, event_time = NULL, cap = Inf) {
  start <- ir_time(start); stop <- ir_time(stop)
  event <- ir_event01(event)
  endpoint <- stop
  if (!is.null(event_time)) {
    death <- ir_time(event_time)
    use_death <- event == 1 & !is.na(death)
    endpoint[use_death] <- death[use_death]
  }
  duration <- as.numeric(difftime(endpoint, start, units = "days"))
  duration[!is.finite(duration) | duration <= 0] <- NA_real_
  if (is.finite(cap)) duration <- pmin(duration, cap)
  duration
}

ir_prepare_cohort <- function(cohort, path = ir_input_path(cohort)) {
  d <- read.csv(path, check.names = FALSE)
  if (!nrow(d)) stop(cohort, "受限输入为空。", call. = FALSE)

  age_raw <- ir_num(ir_get(d, c("age", "admission_age")))
  elderly <- if ('age_deidentified' %in% names(d)) d$age_deidentified == 1 else cohort == 'MIMIC-III' & age_raw >= 300
  age_raw[!is.na(elderly) & elderly] <- 90
  age <- ir_num(age_raw, 18, 120)
  sex0 <- toupper(trimws(as.character(ir_get(d, c("sex", "gender")))))
  sex <- factor(ifelse(sex0 %in% c("M", "MALE", "1"), "Male",
                       ifelse(sex0 %in% c("F", "FEMALE", "0"), "Female", NA_character_)),
                levels = c("Female", "Male"))
  sofa_candidates <- if (cohort == "NWICU") c("sofa_day1", "sofa_day1_adapted", "sofa", "sofa_score") else c("sofa", "sofa_score")
  sofa <- ir_num(ir_get(d, sofa_candidates), 0, 30)
  charlson <- ir_num(ir_get(d, c("charlson_comorbidity_index", "charlson", "charlson_score")), 0, 40)
  diabetes <- factor(ir_event01(ir_get(d, c("diabetes", "icd_diabetes", "Diabetes"))), levels = c(0, 1), labels = c("No", "Yes"))

  # 宽松的数据质量边界：用于清除非有限值、9999999类哨兵值及明显录入错误，
  # 不是临床正常参考区间，也不是依据结局或本研究效应估计选择的截点。
  glucose <- ir_num(ir_get(d, c("glucose", "first_day_glucose_mean", "glucose_mean")), 10, 2000)
  hba1c <- ir_num(ir_get(d, c("hba1c", "first_day_hba1c_mean", "hba1c_mean")), 2, 25)
  tg <- ir_num(ir_get(d, c("tg", "triglycerides", "first_day_triglycerides_mean")), 1, 10000)
  hdl <- ir_num(ir_get(d, c("hdl", "first_day_hdl_mean")), 1, 500)
  idx <- ir_derive_indices(glucose, hba1c, tg, hdl)

  hosp_event <- ir_event01(ir_get(d, c("hospital_outcome", "hospital_mortality")))
  icu_event <- ir_event01(ir_get(d, c("icu_outcome", "icu_mortality")))
  intime <- ir_get(d, c("icu_intime", "intime"))
  outtime <- ir_get(d, c("icu_outtime", "outtime"))
  dischtime <- ir_get(d, c("dischtime", "hospital_dischtime"))
  death_raw <- ir_get(d, c("deathtime", "dod", "death_time"), required = FALSE)

  hospital_time <- ir_duration_from_icu(intime, dischtime, hosp_event, death_raw)
  icu_time <- ir_duration_from_icu(intime, outtime, icu_event, death_raw)
  death_time <- ir_time(death_raw)
  death_days <- as.numeric(difftime(death_time, ir_time(intime), units = "days"))
  death28 <- as.integer(is.finite(death_days) & death_days >= 0 & death_days <= 28)
  time28 <- ifelse(death28 == 1, pmax(death_days, 1 / 86400), 28)

  # 若标准化提取已经提供从ICU入科起算的核验字段，以时间戳重建值为主，
  # 仅在源时间戳无法解析时回退到明确命名的duration字段。
  fallback_h <- ir_num(ir_get(d, c("hospital_time_days", "time_hospital_from_icu_days"), required = FALSE), 0, Inf)
  fallback_i <- ir_num(ir_get(d, c("icu_time_days", "time_icu_from_icu_days"), required = FALSE), 0, Inf)
  fallback_28 <- ir_num(ir_get(d, c("time_28d_days"), required = FALSE), 0, 28)
  hospital_time[is.na(hospital_time)] <- fallback_h[is.na(hospital_time)]
  icu_time[is.na(icu_time)] <- fallback_i[is.na(icu_time)]
  time28[is.na(time28)] <- fallback_28[is.na(time28)]
  provided28_raw <- ir_get(d, c("outcome_28d", "mortality_28d"), required = FALSE)
  provided28 <- ir_event01(provided28_raw)
  if (all(is.na(death_time))) {
    if (all(is.na(provided28_raw)) || all(is.na(fallback_28))) {
      stop(cohort, "缺少从ICU入科起算的精确死亡时间，且未同时提供outcome_28d与time_28d_days；不得用capped survival_time<28近似28天死亡。",
           call. = FALSE)
    }
    death28 <- provided28
    time28 <- fallback_28
  } else {
    use_provided <- is.na(death_time) & !is.na(provided28_raw) & is.finite(fallback_28)
    death28[use_provided] <- provided28[use_provided]
    time28[use_provided] <- fallback_28[use_provided]
  }

  out <- data.frame(
    cohort = cohort, source_row = seq_len(nrow(d)),
    subject_id = ir_get(d, 'subject_id'), stay_id = ir_get(d, c('stay_id','icustay_id')),
    age_deidentified = as.integer(elderly), age = age, gender = sex,
    sofa = sofa, charlson = charlson, diabetes = diabetes,
    glucose = glucose, hba1c = hba1c, tg = tg, hdl = hdl,
    SHR = idx$SHR, TyG = idx$TyG, AIP = idx$AIP,
    hospital_time_days = hospital_time, hospital_outcome = hosp_event,
    icu_time_days = icu_time, icu_outcome = icu_event,
    time_28d_days = time28, outcome_28d = death28,
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(out$subject_id)) stop('Current input must contain one first ICU stay per patient.')
  for (v in c('admission_type','first_careunit','diagnosis_category','dbsource'))
    out[[v]] <- as.character(ir_get(d,v,required=FALSE))
  for (v in c('covid_u071','coronavirus_b9729','diagnosis_available','steroid_24h','coded_sepsis','coded_head_trauma'))
    out[[v]] <- ir_num(ir_get(d,v,required=FALSE))

  # 非暴露敏感性分析所需字段若存在则原样映射；缺失字段保留NA并由接口审计。
  optional_map <- list(
    liver_failure = c("liver_failure"), renal_failure = c("renal_failure", "charlson_renal_disease"),
    longterm_insulin = c("longterm_insulin", "long_term_insulin"),
    longterm_antidiabetics = c("longterm_antidiabetics", "long_term_antidiabetics"),
    sapsii = c("sapsii", "saps_ii"), apsiii = c("apsiii", "aps_iii"),
    oasis = c("oasis"), gcs = c("gcs", "gcs_min"), sirs = c("sirs"),
    ventilation = c("ventilation", "ventilation_icu"), vasopressor = c("vasopressor", "vasopressor_icu")
  )
  for (nm in names(optional_map)) out[[nm]] <- ir_num(ir_get(d, optional_map[[nm]], required = FALSE))
  if (all(is.na(out$liver_failure)) && any(c("charlson_mild_liver_disease", "charlson_severe_liver_disease") %in% names(d))) {
    mild <- ir_event01(ir_get(d, "charlson_mild_liver_disease", required = FALSE))
    severe <- ir_event01(ir_get(d, "charlson_severe_liver_disease", required = FALSE))
    out$liver_failure <- as.integer(mild == 1 | severe == 1)
    out$liver_failure[is.na(mild) & is.na(severe)] <- NA_integer_
  }
  attr(out, "source_path") <- normalizePath(path, mustWork = TRUE)
  out
}

ir_common_exposure_sample <- function(d) {
  ok <- complete.cases(d[, c("SHR", "TyG", "AIP")]) &
    is.finite(d$SHR) & is.finite(d$TyG) & is.finite(d$AIP)
  d[ok, , drop = FALSE]
}

ir_model_covariates <- list(
  "Model 1" = character(),
  "Model 2" = c("age", "gender"),
  "Model 3" = c("age", "gender", "sofa", "charlson", "diabetes")
)

ir_outcomes <- list(
  hospital_mortality = c("hospital_time_days", "hospital_outcome"),
  mortality_28d = c("time_28d_days", "outcome_28d"),
  icu_mortality = c("icu_time_days", "icu_outcome")
)

ir_make_imputations <- function(d, m, seed) {
  vars <- ir_model_covariates[["Model 3"]]
  fixed <- c("SHR", "TyG", "AIP", unlist(ir_outcomes))
  frame <- d[, unique(c(vars, fixed)), drop = FALSE]
  if (!anyNA(frame[, vars, drop = FALSE])) return(replicate(m, d, simplify = FALSE))
  if (!requireNamespace("mice", quietly = TRUE)) stop("缺少R包mice。", call. = FALSE)
  method <- mice::make.method(frame); method[] <- ""
  for (v in vars) if (anyNA(frame[[v]])) {
    method[v] <- if (is.factor(frame[[v]])) "logreg" else "pmm"
  }
  pred <- mice::make.predictorMatrix(frame); diag(pred) <- 0
  pred[setdiff(names(frame), vars), ] <- 0
  imp <- mice::mice(frame, m = m, maxit = 10, method = method,
                    predictorMatrix = pred, seed = seed, printFlag = FALSE)
  completed <- lapply(seq_len(m), function(i) {
    z <- d; ci <- mice::complete(imp, i); z[, vars] <- ci[, vars]; z
  })
  attr(completed, "mids") <- imp
  completed
}

ir_risk_set <- function(d, outcome, landmark_days = 0) {
  def <- ir_outcomes[[outcome]]
  time <- ir_num(d[[def[1]]], 0, Inf); event <- ir_event01(d[[def[2]]])
  keep <- is.finite(time) & time > landmark_days & !is.na(event)
  z <- d[keep, , drop = FALSE]
  z$.time <- time[keep] - landmark_days
  z$.event <- event[keep]
  z
}

ir_pool_scalar <- function(q, u, dfcom, k) {
  if (length(q)==1L || (all(q==q[1]) && all(u==u[1]))) {
    se <- sqrt(u[1]); b <- q[1]
    return(c(beta=b,low=b-qnorm(.975)*se,high=b+qnorm(.975)*se,
             p=2*pnorm(abs(b/se),lower.tail=FALSE),df=Inf,fmi=0))
  }
  if (!requireNamespace("mice", quietly = TRUE)) stop("缺少R包mice。", call. = FALSE)
  p <- mice::pool.scalar(Q = q, U = u, n = dfcom, k = k, rule = "rubin1987")
  se <- sqrt(p$t); crit <- stats::qt(.975, df = p$df)
  c(beta = p$qbar, low = p$qbar - crit * se, high = p$qbar + crit * se,
    p = 2 * stats::pt(abs(p$qbar / se), df = p$df, lower.tail = FALSE),
    df = p$df, fmi = p$fmi)
}

ir_fit_cox_term <- function(datasets, exposure, outcome, model, scale,
                            landmark_days, cuts, x_mean, x_sd) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("缺少R包survival。", call. = FALSE)
  fits <- lapply(datasets, function(dat) {
    z <- ir_risk_set(dat, outcome, landmark_days)
    if (scale == "per_unit") { z$.x <- z[[exposure]]; rhs <- ".x" }
    if (scale == "per_sd") { z$.x <- (z[[exposure]] - x_mean) / x_sd; rhs <- ".x" }
    if (scale == "tertile") {
      z$.tertile <- cut(z[[exposure]], breaks = cuts, labels = c("T1", "T2", "T3"), include.lowest = TRUE)
      rhs <- ".tertile"
    }
    if (scale == "p_trend") {
      z$.trend <- as.integer(cut(z[[exposure]], breaks = cuts, labels = FALSE, include.lowest = TRUE))
      rhs <- ".trend"
    }
    covars <- ir_model_covariates[[model]]
    survival::coxph(stats::as.formula(paste("survival::Surv(.time,.event) ~", paste(c(rhs, covars), collapse = " + "))),
                    data = z, ties = "efron", x = TRUE, model = TRUE)
  })
  terms <- switch(scale, per_unit = ".x", per_sd = ".x",
                  tertile = c(".tertileT2", ".tertileT3"), p_trend = ".trend")
  if (length(unique(vapply(fits, function(f) f$n, numeric(1)))) != 1 ||
      length(unique(vapply(fits, function(f) f$nevent, numeric(1)))) != 1) {
    stop("插补间分母或事件数不一致：", exposure, "/", outcome, "/", model, "/", scale, call. = FALSE)
  }
  lapply(terms, function(term) {
    q <- vapply(fits, function(f) unname(stats::coef(f)[term]), numeric(1))
    u <- vapply(fits, function(f) unname(stats::vcov(f)[term, term]), numeric(1))
    pooled <- ir_pool_scalar(q, u, fits[[1]]$nevent, length(stats::coef(fits[[1]])))
    list(term = term, pooled = pooled, n = fits[[1]]$n, events = fits[[1]]$nevent)
  })
}

ir_run_cox_matrix <- function(d, m, time_origin = c("AdmissionOrigin", "Landmark24h"),
                              scales = c("per_unit", "tertile", "p_trend"),
                              analysis_level = "PA", seed = 20260725) {
  time_origin <- match.arg(time_origin)
  d <- ir_common_exposure_sample(d)
  if (!nrow(d)) stop("三指标共同完整样本为空。", call. = FALSE)
  datasets <- ir_make_imputations(d, m, seed)
  landmark <- if (time_origin == "Landmark24h") 1 else 0
  rows <- list()
  for (exposure in c("SHR", "TyG", "AIP")) {
    x <- d[[exposure]]; cuts <- as.numeric(stats::quantile(x, c(0, 1/3, 2/3, 1), na.rm = TRUE))
    cuts[1] <- -Inf; cuts[4] <- Inf
    x_sd <- stats::sd(x, na.rm = TRUE); x_mean <- mean(x, na.rm = TRUE)
    if (!is.finite(x_sd) || x_sd <= 0) stop("指标标准差无效：", exposure, call. = FALSE)
    for (outcome in names(ir_outcomes)) for (model in names(ir_model_covariates)) for (scale in scales) {
      fit <- ir_fit_cox_term(datasets, exposure, outcome, model, scale, landmark, cuts, x_mean, x_sd)
      for (z in fit) {
        term_label <- switch(z$term, .x = if (scale == "per_sd") "每升高1个标准差" else "每升高1单位",
                             .tertileT2 = "T2相对T1", .tertileT3 = "T3相对T1",
                             .trend = "每升高1个三分位等级", z$term)
        rows[[length(rows) + 1]] <- data.frame(
          analysis_level = analysis_level, cohort = unique(d$cohort), time_origin = time_origin,
          imputation = paste0("M", m, " Rubin pooled"), exposure = exposure,
          exposure_scale = scale, outcome = outcome, model = model, term = term_label,
          effect_measure = "HR", n = z$n, events = z$events,
          estimate = exp(z$pooled["beta"]), conf_low = exp(z$pooled["low"]),
          conf_high = exp(z$pooled["high"]), p_value = z$pooled["p"],
          rubin_df = z$pooled["df"], fraction_missing_information = z$pooled["fmi"],
          m = m, maxit = 10, sample_definition = "SHR、TyG和AIP共同完整样本",
          denominator_note = if (landmark == 0) "从ICU入科起算的结局特异风险集。" else "仅纳入该结局在ICU入科后24小时仍处于风险集者，随访时间减去24小时。",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

ir_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

ir_write_block <- function(path, title, reason, required_fields = character()) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(paste0("# ", title), "", paste0("状态：BLOCKED（", Sys.time(), "）"), "",
               paste0("原因：", reason),
               if (length(required_fields)) paste0("所需字段：", paste(required_fields, collapse = "、")) else NULL,
               "", "未生成或沿用任何数值结果；旧结果不得进入revision v3稿件。"), path, useBytes = TRUE)
}

ir_check_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("缺少R包：", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}
