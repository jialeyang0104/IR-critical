#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
root<-Sys.getenv('IR_REVISION_ROOT');source(file.path(root,'08_CR_竞争风险修订/02_代码/01_CR_CORE.R'))
requested<-strsplit(Sys.getenv('IR_CR_COHORTS','MIMIC-IV,MIMIC-III,NWICU'),',',fixed=TRUE)[[1]]
cr_audit<-list();missing_audit<-list();rcs_all<-curves_all<-tests_all<-risk_all<-cif_all<-list()
cr_spline <- function(d,sets,x,y) {
  knots<-as.numeric(quantile(d[[x]],c(.05,.35,.65,.95)));ref<-median(d[[x]])
  grid<-seq(quantile(d[[x]],.01),quantile(d[[x]],.99),length.out=100)
  fn<-function(dat) {
    z<-cr_risk(dat,y);b<-Hmisc::rcspline.eval(z[[x]],knots=knots,inclx=TRUE)
    colnames(b)<-paste0('.rcs',seq_len(ncol(b)));z<-cbind(z,b)
    cr_fit(z,y,colnames(b),ir_model_covariates[['Model 3']])
  }
  fits<-cr_apply_mi(sets,fn,ir_model_covariates[['Model 3']]);tt<-grep('^\\.rcs',names(fits[[1]]$b),value=TRUE)
  q<-do.call(cbind,lapply(fits,function(f) f$b[tt]));qb<-rowMeans(q)
  u<-Reduce('+',lapply(fits,function(f) f$v[tt,tt,drop=FALSE]))/length(fits)
  total<-u+(1+1/length(fits))*cov(t(q))
  bg<-Hmisc::rcspline.eval(grid,knots=knots,inclx=TRUE);br<-as.numeric(Hmisc::rcspline.eval(ref,knots=knots,inclx=TRUE))
  delta<-sweep(bg,2,br,'-');eta<-as.numeric(delta%*%qb);se<-sqrt(rowSums((delta%*%total)*delta))
  wald<-function(b,v) as.numeric(crossprod(b,solve(v,b)))
  data.frame(cohort=unique(d$cohort),exposure=x,outcome=y,model='Model 3',effect_measure=cr_measure(y),
    x=grid,reference=ref,knots=paste(knots,collapse=';'),estimate=exp(eta),conf_low=exp(eta-1.96*se),conf_high=exp(eta+1.96*se),
    p_overall=pchisq(wald(qb,total),length(qb),lower.tail=FALSE),
    p_nonlinear=pchisq(wald(qb[-1],total[-1,-1,drop=FALSE]),length(qb)-1,lower.tail=FALSE),
    n=fits[[1]]$n,events=fits[[1]]$events)
}
cr_curves <- function(d,x,y) {
  z<-cr_risk(d,y);cuts<-quantile(d[[x]],c(0,1/3,2/3,1));cuts[c(1,4)]<-c(-Inf,Inf)
  z$.group<-cut(z[[x]],cuts,labels=c('T1','T2','T3'),include.lowest=TRUE)
  curves<-risks<-cifs<-list()
  if(y!='mortality_28d') {
    ci<-cmprsk::cuminc(z$.time,z$.status,z$.group,cencode=0)
    pv<-ci$Tests['1','pv'];test<-'Gray'
    for(g in levels(z$.group)) {
      a<-ci[[paste(g,'1')]]; sdv<-sqrt(pmax(a$var,0))
      curves[[g]]<-data.frame(tertile=g,time=a$time,probability=a$est,conf_low=pmax(0,a$est-1.96*sdv),conf_high=pmin(1,a$est+1.96*sdv),curve_type='CIF')
      for(t in c(7,14,28)) {
        j<-max(which(a$time<=t));p<-a$est[j];ss<-sqrt(a$var[j])
        cifs[[length(cifs)+1]]<-data.frame(tertile=g,time=t,estimate=p,conf_low=max(0,p-1.96*ss),conf_high=min(1,p+1.96*ss))
      }
    }
  } else {
    sf<-survival::survfit(survival::Surv(.time,.event)~.group,z)
    lr<-survival::survdiff(survival::Surv(.time,.event)~.group,z);pv<-pchisq(lr$chisq,2,lower.tail=FALSE);test<-'Log-rank'
    for(g in levels(z$.group)) {
      f<-sf[paste0('.group=',g)]
      curves[[g]]<-data.frame(tertile=g,time=c(0,f$time),probability=c(1,f$surv),conf_low=c(1,f$lower),conf_high=c(1,f$upper),curve_type='KM')
    }
  }
  for(g in levels(z$.group)) for(t in c(0,7,14,21,28)) risks[[length(risks)+1]]<-data.frame(tertile=g,time=t,n_risk=sum(z$.group==g & z$.time>=t))
  meta<-function(a) {a$cohort<-unique(d$cohort);a$exposure<-x;a$outcome<-y;a}
  list(curves=meta(do.call(rbind,curves)),risk=meta(do.call(rbind,risks)),
       tests=data.frame(cohort=unique(d$cohort),exposure=x,outcome=y,test=test,p_value=pv),
       cif=if(length(cifs))meta(do.call(rbind,cifs))else NULL)
}
for(cohort in requested) {
  message('START ',cohort);code<-cr_codes[[cohort]];d<-cr_prepare(cohort);common<-ir_common_exposure_sample(d)
  for(y in cr_outcomes) {
    z<-cr_risk(common,y)
    cr_audit[[length(cr_audit)+1]]<-data.frame(cohort=cohort,outcome=y,n=nrow(z),events=sum(z$.event),competing_events=sum(z$.status==2),censored=sum(z$.status==0),missing_state=sum(is.na(z$.status)),
      endpoint=if(y=='hospital_mortality')'本次住院死亡/存活出院' else if(y=='icu_mortality')'首次ICU住院死亡/存活离开首次ICU' else '入科后28天内记录死亡；其余第28天行政删失')
  }
  missing_audit[[cohort]]<-data.frame(cohort=cohort,variable=ir_model_covariates[['Model 3']],n=nrow(common),missing=vapply(common[,ir_model_covariates[['Model 3']]],function(a)sum(is.na(a)),integer(1)))
  jobs<-list(AO_M5=list(5,'AdmissionOrigin',0L,'PA'),L24_M5=list(5,'Landmark24h',2L,'SA'),AO_M20=list(20,'AdmissionOrigin',1L,'SA'),L24_M20=list(20,'Landmark24h',3L,'SA'))
  for(job in names(jobs)) {
    j<-jobs[[job]];a<-cr_matrix(d,j[[1]],j[[2]],c('per_unit','tertile','p_trend'),j[[4]],cr_seeds[[cohort]]+j[[3]])
    cr_write(a,paste0('01_CR_',code,'_',job,'_regression_results.csv'));message(cohort,' ',job,' ',nrow(a),' rows')
  }
  cr_write(cr_matrix(d,5,'AdmissionOrigin','per_sd','PSA',cr_seeds[[cohort]]),paste0('02_CR_',code,'_AO_M5_per_sd_results.csv'))
  sets<-ir_make_imputations(common,5,cr_seeds[[cohort]])
  for(x in c('SHR','TyG','AIP')) for(y in cr_outcomes) {
    rcs_all[[length(rcs_all)+1]]<-cr_spline(common,sets,x,y)
    cv<-cr_curves(common,x,y)
    curves_all[[length(curves_all)+1]]<-cv$curves;tests_all[[length(tests_all)+1]]<-cv$tests;risk_all[[length(risk_all)+1]]<-cv$risk
    if(!is.null(cv$cif))cif_all[[length(cif_all)+1]]<-cv$cif
  }
  message('CURVES ',cohort,' complete')

}
cr_write(do.call(rbind,cr_audit),'00_CR_endpoint_state_audit.csv')
cr_write(do.call(rbind,missing_audit),'00_CR_covariate_missingness.csv')
cr_write(do.call(rbind,rcs_all),'03_CR_ALL_AO_M5_rcs_source_data.csv')
cr_write(do.call(rbind,curves_all),'04_CR_ALL_AO_tertile_curve_coordinates.csv')
cr_write(do.call(rbind,tests_all),'04_CR_ALL_AO_Gray_logrank_tests.csv')
cr_write(do.call(rbind,risk_all),'04_CR_ALL_AO_fixed_time_risk_counts.csv')
cr_write(do.call(rbind,cif_all),'04_CR_ALL_AO_CIF_7_14_28_days.csv')
if ('MIMIC-IV' %in% requested) source(file.path(root, '08_CR_竞争风险修订/02_代码/02b_CR_RUN_SUBGROUPS.R'))
writeLines(c('主分析与曲线运行完成。',capture.output(sessionInfo())),file.path(crdir,'07_运行日志/02_CR_main_session.txt'))
