#!/usr/bin/env Rscript
# 当前 MIMIC-IV 亚组效应和交互作用；使用共同队列 SD，不输出患者明细。
options(stringsAsFactors = FALSE, warn = 1)
root <- Sys.getenv('IR_REVISION_ROOT')
source(file.path(root, 'src/analysis/00_core.R'))
d <- ir_common_exposure_sample(cr_prepare('MIMIC-IV'))
stopifnot(nrow(d) == 3056, !anyNA(d[, ir_model_covariates[['Model 3']]]))
groups <- list(
  age = list(variable = 'age', reference = '<65', comparison = '>=65',
             value = factor(ifelse(d$age < 65, '<65', '>=65'), levels = c('<65', '>=65'))),
  sex = list(variable = 'gender', reference = 'Female', comparison = 'Male', value = d$gender),
  diabetes = list(variable = 'diabetes', reference = 'No', comparison = 'Yes', value = d$diabetes))
subgroup_rows <- interaction_rows <- list()
for (x in c('SHR', 'TyG', 'AIP')) {
  index_mean <- mean(d[[x]])
  index_sd <- sd(d[[x]])
  for (sg in names(groups)) {
    definition <- groups[[sg]]
    dat <- d
    dat$.x <- (dat[[x]] - index_mean) / index_sd
    dat$.modifier <- as.integer(definition$value == definition$comparison)
    dat$.interaction <- dat$.x * dat$.modifier
    # For age strata, retain continuous age adjustment. Sex/diabetes main
    # effects are represented by .modifier in the interaction model.
    covars <- setdiff(ir_model_covariates[['Model 3']],
                     if (sg == 'age') character() else definition$variable)
    for (y in cr_outcomes) {
      z <- cr_risk(dat, y)
      fit <- cr_fit(z, y, c('.x', '.modifier', '.interaction'), covars)
      beta <- unname(fit$b['.interaction'])
      se <- sqrt(fit$v['.interaction', '.interaction'])
      wald <- (beta / se)^2
      p <- pchisq(wald, df = 1, lower.tail = FALSE)
      # Independent calculation checks the 1-df two-sided Wald p value.
      stopifnot(is.finite(se), se > 0,
                isTRUE(all.equal(p, 2 * pnorm(-abs(beta / se)), tolerance = 1e-12)))
      interaction_rows[[length(interaction_rows) + 1L]] <- data.frame(
        cohort = 'MIMIC-IV', exposure = x, outcome = y, subgroup = sg,
        reference_level = definition$reference, comparison_level = definition$comparison,
        model = 'Model 3', n = fit$n, events = fit$events,
        reference_n = sum(z$.modifier == 0), comparison_n = sum(z$.modifier == 1),
        reference_events = sum(z$.modifier == 0 & z$.event == 1),
        comparison_events = sum(z$.modifier == 1 & z$.event == 1),
        exposure_scale = 'per_common_cohort_sd', index_mean = index_mean, index_sd = index_sd,
        beta_interaction = beta, se_interaction = se,
        interaction_ratio = exp(beta), conf_low = exp(beta - qnorm(.975) * se),
        conf_high = exp(beta + qnorm(.975) * se), wald_chisq = wald, df = 1L,
        p_interaction = p, regression_method = fit$type,
        test = 'Two-sided 1-df Wald test of index-by-subgroup product term',
        adjustment_covariates = paste(covars, collapse = ';'),
        convergence = if (y == 'mortality_28d') all(is.finite(fit$b)) else fit$fit$converged)
      for (lev in 0:1) {
        zz <- z[z$.modifier == lev, , drop = FALSE]
        subfit <- cr_fit(zz, y, '.x', covars)
        subgroup_rows[[length(subgroup_rows) + 1L]] <- cr_row(
          list(subfit), '.x', 'MIMIC-IV', x, y, 'Model 3', 'per_sd',
          'AdmissionOrigin', 0, 'PSA', list(
            subgroup = sg, level = if (lev == 0) definition$reference else definition$comparison,
            standardization = 'SD of the full MIMIC-IV common-index cohort',
            index_mean = index_mean, index_sd = index_sd))
      }
    }
  }
}
interactions <- do.call(rbind, interaction_rows)
interactions$p_interaction_bh <- p.adjust(interactions$p_interaction, method = 'BH')
interactions$multiplicity_family <- 'All 27 index-by-subgroup tests (3 indices x 3 outcomes x 3 modifiers)'
subgroups <- do.call(rbind, subgroup_rows)
subgroups <- merge(subgroups, interactions[, c('cohort', 'exposure', 'outcome', 'subgroup',
  'p_interaction', 'p_interaction_bh')], by = c('cohort', 'exposure', 'outcome', 'subgroup'), sort = FALSE)
stopifnot(nrow(interactions) == 27, nrow(subgroups) == 54,
          all(interactions$convergence), all(interactions$n == 3056),
          all(interactions$reference_n + interactions$comparison_n == interactions$n))
cr_write(subgroups, '05_CR_M4_AO_M5_exploratory_subgroups.csv')
cr_write(interactions, '05_CR_M4_AO_subgroup_interactions.csv')
writeLines(c('Completed: 54 common-SD subgroup estimates and 27 interaction Wald tests.',
             'BH adjustment family: 27 interaction tests.', capture.output(sessionInfo())),
           file.path(crdir, 'reports/internal/02b_CR_subgroup_session.txt'))
message('Completed 54 subgroup estimates and 27 interaction tests; nominal P<0.05: ',
        sum(interactions$p_interaction < .05), '; BH P<0.05: ', sum(interactions$p_interaction_bh < .05))
