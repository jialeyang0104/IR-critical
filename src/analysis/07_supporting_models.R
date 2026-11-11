#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
root<-Sys.getenv('IR_REVISION_ROOT');source(file.path(root,'src/analysis/00_core.R'))
cox<-logistic<-discharge<-list()
for(cohort in names(cr_codes)) {
  d<-ir_common_exposure_sample(cr_prepare(cohort));sets<-ir_make_imputations(d,5,cr_seeds[[cohort]])
  for(x in c('SHR','TyG','AIP'))for(y in c('hospital_mortality','icu_mortality'))for(event in c('death','alive_discharge')) {
    z<-cr_risk(d,y);z$.x<-as.numeric(scale(z[[x]]))
    z$.event<-as.integer(z$.status==if(event=='death')1 else 2)
    f<-cr_fit(z,y,'.x',ir_model_covariates[['Model 3']],force_cause_specific=TRUE)
    r<-cr_row(list(f),'.x',cohort,x,y,'Model 3','per_sd','AdmissionOrigin',0,'SA',list(cause=event));r$effect_measure<-'HR'
    discharge[[length(discharge)+1]]<-r
  }
  for(x in c('SHR','TyG','AIP'))for(y in cr_outcomes)for(model in names(ir_model_covariates)) {
    covars<-ir_model_covariates[[model]]
    f<-cr_apply_mi(sets,function(dat){z<-cr_risk(dat,y);z$.x<-z[[x]];cr_fit(z,y,'.x',covars,force_cause_specific=TRUE)},covars)
    r<-cr_row(f,'.x',cohort,x,y,model,'per_unit','AdmissionOrigin',5,'SA');r$effect_measure<-'HR'
    cox[[length(cox)+1]]<-r
    lf<-cr_apply_mi(sets,function(dat){z<-cr_risk(dat,y);z$.x<-z[[x]]
      fit<-glm(as.formula(paste('.event ~',paste(c('.x',covars),collapse=' + '))),z,family=binomial())
      list(b=coef(fit),v=vcov(fit),n=nrow(z),events=sum(z$.event),competing=sum(z$.status==2),censored=sum(z$.status==0),k=length(coef(fit)),type='fixed_binary_logistic')},covars)
    r<-cr_row(lf,'.x',cohort,x,y,model,'per_unit','AdmissionOrigin',5,'SA');r$effect_measure<-'OR'
    p<-ir_pool_scalar(vapply(lf,function(a)unname(a$b['.x']),numeric(1)),vapply(lf,function(a)unname(a$v['.x','.x']),numeric(1)),lf[[1]]$n,lf[[1]]$k)
    r$estimate<-exp(p['beta']);r$conf_low<-exp(p['low']);r$conf_high<-exp(p['high']);r$p_value<-p['p'];r$rubin_df<-p['df'];r$fraction_missing_information<-p['fmi']
    logistic[[length(logistic)+1]]<-r
  }
}
cr_write(do.call(rbind,cox),'10a_CR_ALL_AO_M5_cause_specific_cox_results.csv')
cr_write(do.call(rbind,logistic),'10b_CR_ALL_AO_M5_logistic_results.csv')
cr_write(do.call(rbind,discharge),'10c_CR_ALL_death_discharge_cause_specific.csv')
writeLines(c('原因别Cox与固定二分类Logistic支持模型运行完成。',capture.output(sessionInfo())),file.path(crdir,'reports/internal/03b_CR_support_session.txt'))
