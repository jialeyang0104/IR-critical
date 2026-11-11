#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
root<-Sys.getenv('IR_REVISION_ROOT');source(file.path(root,'08_CR_竞争风险修订/02_代码/01_CR_CORE.R'))
seeds<-c('MIMIC-IV'=41001L,'MIMIC-III'=42001L,'NWICU'=43001L)
# 复用经过核对的时间取值构建函数；只重定向模型与聚合输出，不运行旧脚本顶层代码。
labenv<-new.env(parent=globalenv())
exprs<-parse(file.path(root,'00_总说明与运行顺序/20_RUN_ADVANCED_SENSITIVITIES.R'))
for(ex in exprs) if(is.call(ex)&&identical(ex[[1]],as.name('<-'))&&as.character(ex[[2]])=='run_lab_timing') eval(ex,labenv)
labenv$codes<-cr_codes;labenv$dirs<-c('MIMIC-IV'='01_MIMICIV','MIMIC-III'='02_MIMICIII外部验证','NWICU'='03_NWICU外部验证')
labenv$seeds<-seeds;labenv$order_prefix<-function(c) '06';labenv$ir_prepare_cohort<-cr_prepare
labenv$ir_run_cox_matrix<-function(d,m,origin,scales,level,seed) cr_matrix(d,m,origin,c('per_unit','per_sd'),level,seed)
labenv$ir_write_csv<-function(x,path) {
  nm<-basename(path);nm<-sub('^06_SA_','06_CR_',nm)
  cr_write(x,nm)
}
all_exclusions<-all_weights<-all_balance<-all_wsummary<-field_audit<-weighted_check<-ph<-list()
for(cohort in names(cr_codes)) {
  code<-cr_codes[[cohort]];message('START ',cohort);d<-cr_prepare(cohort);common<-ir_common_exposure_sample(d)
  labenv$run_lab_timing(ir_input_path(cohort),cohort);message('LAB ',cohort,' complete')
  sets<-ir_make_imputations(common,5,seeds[[cohort]]+1L)
  flags<-c(exclude_diabetes='diabetes',exclude_liver_failure='liver_failure',exclude_renal_failure='renal_failure',exclude_longterm_insulin='longterm_insulin',exclude_longterm_antidiabetics='longterm_antidiabetics')
  for(exc in names(flags)) {
    flag<-flags[[exc]];available<-any(!is.na(d[[flag]]))
    field_audit[[length(field_audit)+1]]<-data.frame(cohort=cohort,exclusion=exc,executable=available,reason=if(available)'字段可用' else '受限提取缺少可靠长期用药史；未以阴性诊断标志替代')
    if(!available)next
    covars<-setdiff(ir_model_covariates[['Model 3']],if(flag=='diabetes')'diabetes'else character())
    for(x in c('SHR','TyG','AIP'))for(y in cr_outcomes)for(sc in c('per_unit','per_sd')) {
      xm<-mean(common[[x]]);xs<-sd(common[[x]])
      fits<-cr_apply_mi(sets,function(dat){keep<-ir_event01(dat[[flag]])==0;z<-cr_risk(dat[!is.na(keep)&keep,],y);z$.x<-if(sc=='per_sd')(z[[x]]-xm)/xs else z[[x]];cr_fit(z,y,'.x',covars)},covars)
      all_exclusions[[length(all_exclusions)+1]]<-cr_row(fits,'.x',cohort,x,y,'Model 3',sc,'AdmissionOrigin',5,'SA',list(exclusion=exc))
    }
  }
  message('EXCLUSIONS ',cohort,' complete')
  z<-d;z$selected<-as.integer(complete.cases(z[,c('SHR','TyG','AIP')]))
  for(v in c('age','sofa','charlson'))z[[v]][is.na(z[[v]])]<-median(z[[v]],na.rm=TRUE)
  for(v in c('gender','diabetes')) {fill<-names(which.max(table(z[[v]])));z[[v]][is.na(z[[v]])]<-fill;z[[v]]<-droplevels(z[[v]])}
  for(v in c('admission_type','first_careunit','diagnosis_category')) {
    z[[v]][is.na(z[[v]])|z[[v]]==''] <- 'Unknown'
    rare <- names(which(table(z[[v]])<100))
    z[[v]][z[[v]]%in%rare] <- 'Other (fewer than 100 eligible patients)'
    z[[v]] <- factor(z[[v]])
  }
  psfit<-glm(selected~age+gender+sofa+charlson+diabetes+admission_type+first_careunit+diagnosis_category,z,family=binomial())
  z$ps<-pmin(pmax(predict(psfit,type='response'),1e-4),.9999);sel<-z[z$selected==1,];w<-mean(z$selected)/sel$ps
  trims<-list(unweighted=rep(1,nrow(sel)),untrimmed=w)
  for(nm in c('p01_p99','p025_p975','p05_p95')) {pr<-switch(nm,p01_p99=c(.01,.99),p025_p975=c(.025,.975),p05_p95=c(.05,.95));b<-quantile(w,pr);trims[[nm]]<-pmin(pmax(w,b[1]),b[2])}
  sets<-ir_make_imputations(common,5,seeds[[cohort]]+2L)
  for(tn in names(trims)) {
    wb<-setNames(trims[[tn]],sel$source_row);ess<-sum(trims[[tn]])^2/sum(trims[[tn]]^2)
    for(x in c('SHR','TyG','AIP'))for(y in cr_outcomes) {
      xm<-mean(common[[x]]);xs<-sd(common[[x]])
      fits<-cr_apply_mi(sets,function(dat){zz<-cr_risk(dat,y);zz$.x<-(zz[[x]]-xm)/xs;zz$.w<-unname(wb[as.character(zz$source_row)]);cr_fit(zz,y,'.x',ir_model_covariates[['Model 3']],weighted=TRUE)},ir_model_covariates[['Model 3']])
      all_weights[[length(all_weights)+1]]<-cr_row(fits,'.x',cohort,x,y,'Model 3','per_sd','AdmissionOrigin',5,'SA',list(trimming=tn,effective_sample_size=ess))
      if(tn=='unweighted'&&y!='mortality_28d') {
        zz<-cr_risk(sets[[1]],y);zz$.x<-(zz[[x]]-xm)/xs;ff<-cr_fit(zz,y,'.x',ir_model_covariates[['Model 3']])
        weighted_check[[length(weighted_check)+1]]<-data.frame(cohort=cohort,exposure=x,outcome=y,crr_beta=ff$b['.x'],finegray_expanded_beta=fits[[1]]$b['.x'],difference=abs(ff$b['.x']-fits[[1]]$b['.x']))
      }
    }
  }
  balance_matrix <- model.matrix(~age+sofa+charlson+gender+diabetes+admission_type+first_careunit+diagnosis_category,z)[,-1,drop=FALSE]
  for(v in colnames(balance_matrix)) {
    target<-balance_matrix[,v]
    obs<-target[z$selected==1];den<-sqrt((var(obs)+var(target))/2)
    for(tn in names(trims))all_balance[[length(all_balance)+1]]<-data.frame(cohort=cohort,variable=v,trimming=tn,target_mean=mean(target),weighted_mean=weighted.mean(obs,trims[[tn]]),standardized_mean_difference=(weighted.mean(obs,trims[[tn]])-mean(target))/den)
  }
  for(tn in names(trims))all_wsummary[[length(all_wsummary)+1]]<-data.frame(cohort=cohort,trimming=tn,eligible_n=nrow(d),selected_n=nrow(sel),effective_sample_size=sum(trims[[tn]])^2/sum(trims[[tn]]^2),weight_min=min(trims[[tn]]),weight_max=max(trims[[tn]]))
  message('WEIGHTS ',cohort,' complete')
  # 审计：在Fine-Gray模型加入暴露×标准化log(t)交互，检验暴露亚分布效应的时间变化。
  for(x in c('SHR','TyG','AIP'))for(y in c('hospital_mortality','icu_mortality')) {
    covars<-ir_model_covariates[['Model 3']]
    fits<-cr_apply_mi(sets,function(dat){zz<-cr_risk(dat,y);zz$.x<-as.numeric(scale(zz[[x]]));mm<-model.matrix(~.x+age+gender+sofa+charlson+diabetes,zz)[,-1,drop=FALSE]
      f<-cmprsk::crr(zz$.time,zz$.status,cov1=mm,cov2=as.matrix(zz$.x),tf=function(t)matrix(log1p(t)-mean(log1p(t)),ncol=1),failcode=1,cencode=0)
      idx<-length(f$coef);b<-f$coef;names(b)[idx]<-'.time_effect';v<-f$var;dimnames(v)<-list(names(b),names(b));list(b=b,v=v,n=nrow(zz),events=sum(zz$.event),k=length(b))},covars)
    p<-cr_pool(fits,'.time_effect');ph[[length(ph)+1]]<-data.frame(cohort=cohort,exposure=x,outcome=y,beta_time=p['beta'],p_time=p['p'],diagnostic='per-SD exposure × centered log(1+days); Model 3; complete observed covariates')
  }
}
cr_write(do.call(rbind,all_exclusions),'07_CR_ALL_AO_M5_exclusion_results.csv')
cr_write(do.call(rbind,field_audit),'07_CR_ALL_exclusion_field_audit.csv')
cr_write(do.call(rbind,all_weights),'08_CR_ALL_AO_M5_selection_weighted_results.csv')
cr_write(do.call(rbind,all_balance),'08_CR_ALL_covariate_balance.csv')
cr_write(do.call(rbind,all_wsummary),'08_CR_ALL_weight_summary.csv')
cr_write(do.call(rbind,weighted_check),'09_CR_weighted_implementation_equivalence.csv')
cr_write(do.call(rbind,ph),'09_CR_FineGray_time_effect_audit.csv')
writeLines(c('敏感性运行完成。',capture.output(sessionInfo())),file.path(crdir,'07_运行日志/03_CR_sensitivity_session.txt'))
