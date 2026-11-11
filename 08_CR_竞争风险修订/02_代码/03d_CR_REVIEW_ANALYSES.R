#!/usr/bin/env Rscript
# Patient-level review analyses. All exports are aggregates.
root <- Sys.getenv('IR_REVISION_ROOT')
source(file.path(root,'08_CR_竞争风险修订/02_代码/01_CR_CORE.R'))
audit<-scales<-selection<-special<-diagnostics<-list()
flow<-read.csv(file.path(resdir,'00_CR_patient_flow.csv'))
for(cohort in names(cr_codes)) {
  all<-cr_prepare(cohort);d<-ir_common_exposure_sample(all)
  stopifnot(!anyDuplicated(d$subject_id),!anyNA(d[,ir_model_covariates[['Model 3']]]))
  flow$complete_indices_n[flow$cohort==cohort]<-nrow(d)
  flow$excluded_laboratory[flow$cohort==cohort]<-nrow(all)-nrow(d)
  audit[[cohort]]<-data.frame(cohort=cohort,n=nrow(d),unique_patients=length(unique(d$subject_id)),
    elderly_recoded=sum(d$age_deidentified),age_missing=sum(is.na(d$age)),
    main_age=90,sensitivity_age=91.4,icu_duration_hours_min=min(d$icu_time_days)*24,
    coded_sepsis_n=sum(d$coded_sepsis,na.rm=TRUE),coded_head_trauma_n=sum(d$coded_head_trauma,na.rm=TRUE),
    diagnostic_flags_missing=sum(is.na(d$coded_sepsis)))
  for(x in c('SHR','TyG','AIP')) {
    q<-quantile(d[[x]],c(1/3,2/3))
    scales[[length(scales)+1]]<-data.frame(cohort=cohort,exposure=x,n=nrow(d),mean=mean(d[[x]]),sd=sd(d[[x]]),
      tertile_cut1=q[1],tertile_cut2=q[2])
  }
  for(v in c('admission_type','first_careunit','diagnosis_category')) {
    a<-all[[v]];a[is.na(a)|a=='']<-'Unknown';inc<-all$source_row%in%d$source_row
    for(g in sort(unique(a))) selection[[length(selection)+1]]<-data.frame(cohort=cohort,variable=v,category=g,
      included_n=sum(a==g & inc),included_percent=100*sum(a==g & inc)/sum(inc),
      excluded_n=sum(a==g & !inc),excluded_percent=100*sum(a==g & !inc)/sum(!inc))
  }
  variants<-list()
  if(cohort=='MIMIC-III') {
    a<-d;a$age[a$age_deidentified==1]<-91.4;variants$elderly_age_91.4<-a
    variants$CareVue_only<-d[d$dbsource=='carevue',]
  }
  if(cohort=='NWICU') {
    variants$exclude_COVID_U071<-d[d$covid_u071==0 & d$diagnosis_available==1,]
    variants$exclude_U071_or_B9729<-d[d$covid_u071==0 & d$coronavirus_b9729==0 & d$diagnosis_available==1,]
    variants$exclude_steroid_24h<-d[d$steroid_24h==0,]
    variants$adjust_COVID_and_steroid<-d
    cr_write(data.frame(n=nrow(d),diagnosis_available=sum(d$diagnosis_available),
      COVID_U071=sum(d$covid_u071),coronavirus_B9729=sum(d$coronavirus_b9729),
      U071_or_B9729=sum(d$covid_u071==1|d$coronavirus_b9729==1),
      documented_systemic_steroid_24h=sum(d$steroid_24h)), '11_CR_NW_COVID_steroid_description.csv')
  }
  for(v in names(variants)) for(x in c('SHR','TyG','AIP')) for(y in cr_outcomes) {
    zz<-cr_risk(variants[[v]],y);zz$.x<-(zz[[x]]-mean(d[[x]]))/sd(d[[x]])
    covars<-ir_model_covariates[['Model 3']]
    if(v=='adjust_COVID_and_steroid')covars<-c(covars,'covid_u071','steroid_24h')
    f<-cr_fit(zz,y,'.x',covars)
    special[[length(special)+1]]<-cr_row(list(f),'.x',cohort,x,y,'Model 3','per_sd','AdmissionOrigin',0,'SA',list(sensitivity=v))
  }
  for(x in c('SHR','TyG','AIP'))for(y in cr_outcomes) {
    z<-cr_risk(d,y);z$.x<-as.numeric(scale(z[[x]]))
    f<-cr_fit(z,y,'.x',ir_model_covariates[['Model 3']],force_cause_specific=TRUE)
    ph<-survival::cox.zph(f$fit,transform='log')$table
    for(term in c('.x','GLOBAL'))diagnostics[[length(diagnostics)+1]]<-data.frame(
      cohort=cohort,exposure=x,outcome=y,diagnostic='Cause-specific Cox Schoenfeld test',term=term,
      chisq=ph[term,'chisq'],df=ph[term,'df'],p_value=ph[term,'p'])
  }
}
cr_write(flow,'00_CR_patient_flow.csv')
cr_write(do.call(rbind,audit),'00_CR_common_patient_audit.csv')
cr_write(do.call(rbind,scales),'11_CR_index_scales_tertiles.csv')
cr_write(do.call(rbind,selection),'11_CR_admission_unit_diagnosis_selection.csv')
cr_write(do.call(rbind,special),'11_CR_age_CareVue_COVID_steroid_sensitivities.csv')
cr_write(do.call(rbind,diagnostics),'11_CR_cause_specific_diagnostics.csv')

# Paired patient bootstrap: standardize each index within each resampled cohort.
# The primary comparison is MIMIC-IV in-hospital mortality; other rows are exploratory.
B<-200L
contrast<-list();failures<-list()
for(cohort in names(cr_codes)) {
  d<-ir_common_exposure_sample(cr_prepare(cohort))
  for(y in cr_outcomes) {
    z<-cr_risk(d,y)
    once<-function(a) {
      b<-vapply(c('SHR','TyG','AIP'),function(x){a$.x<-as.numeric(scale(a[[x]]));unname(cr_fit(a,y,'.x',ir_model_covariates[['Model 3']])$b['.x'])},numeric(1))
      c('SHR minus TyG'=unname(b[1]-b[2]),'SHR minus AIP'=unname(b[1]-b[3]))
    }
    point<-once(z);seed<-61000L+match(cohort,names(cr_codes))*100L+match(y,cr_outcomes)
    set.seed(seed);indices<-replicate(B,sample.int(nrow(z),nrow(z),replace=TRUE),simplify=FALSE)
    result<-parallel::mclapply(indices,function(ii)tryCatch(list(value=once(z[ii,]),error=''),error=function(e)list(value=c(NA_real_,NA_real_),error=conditionMessage(e))),mc.cores=4L,mc.set.seed=FALSE)
    bb<-do.call(rbind,lapply(result,`[[`,'value'));bad<-vapply(result,function(a)nzchar(a$error),logical(1))
    for(j in seq_along(point))contrast[[length(contrast)+1]]<-data.frame(cohort=cohort,outcome=y,
      contrast=names(point)[j],log_effect_difference=point[j],conf_low=quantile(bb[,j],.025,na.rm=TRUE),
      conf_high=quantile(bb[,j],.975,na.rm=TRUE),bootstrap_planned=B,bootstrap_success=sum(!bad),seed=seed,
      scale='Log effect per within-resample SD',inference_family=if(cohort=='MIMIC-IV'&&y=='hospital_mortality')'Primary comparison'else'Exploratory')
    if(any(bad))for(i in which(bad))failures[[length(failures)+1]]<-data.frame(cohort=cohort,outcome=y,replicate=i,error=result[[i]]$error)
    message('Paired bootstrap complete: ',cohort,' ',y,' ',sum(!bad),'/',B)
  }
}
cr_write(do.call(rbind,contrast),'11_CR_paired_perSD_effect_contrasts.csv')
cr_write(if(length(failures))do.call(rbind,failures)else data.frame(cohort=character(),outcome=character(),replicate=integer(),error=character()),'11_CR_effect_bootstrap_failures.csv')
writeLines(c('Review analyses complete',capture.output(sessionInfo())),file.path(crdir,'07_运行日志/03d_CR_review_session.txt'))
