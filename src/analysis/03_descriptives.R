#!/usr/bin/env Rscript
# Observed baseline, missingness, selection, and MIMIC-IV survival-status tables.
# Base R only: no regression, imputation, or prediction packages are needed.
root <- Sys.getenv('IR_REVISION_ROOT')
source(file.path(root, 'src/shared/01_analysis_core.R'))
outdir <- file.path(root, 'results/analysis')
codes <- c('MIMIC-IV'='M4', 'MIMIC-III'='M3', 'NWICU'='NW')
labels <- c(gender='Male sex, n (%)', age='Age (years)',
  los_hospital='Hospital length of stay (days)', los_icu='ICU length of stay (days)',
  sofa='SOFA', gcs='GCS', charlson='Charlson comorbidity index', sirs='SIRS',
  HR='Heart rate (beats/min)', RR='Respiratory rate (breaths/min)',
  MBP='Mean blood pressure (mmHg)', Tm='Temperature (°C)', WBC='WBC (10^9/L)',
  PLT='PLT (10^9/L)', Hb='Hb (g/dL)', BUN='BUN (mg/dL)', Scr='Scr (mg/dL)',
  Na='Na (mmol/L)', K='K (mmol/L)', INR='INR', PT='PT (s)', ALT='ALT (IU/L)',
  AST='AST (IU/L)', TG='TG (mg/dL)', Glucose='Glucose (mg/dL)',
  HbA1c='HbA1c (%)', HDL='HDL-C (mg/dL)', Glucocorticoids='Glucocorticoids, n (%)',
  Vasopressors='Vasopressors, n (%)', Ventilation='Mechanical ventilation, n (%)',
  Hypertension='Hypertension, n (%)', MI='MI, n (%)', CHF='CHF, n (%)',
  CPD='CPD, n (%)', Diabetes='Diabetes, n (%)', liver_failure='Liver failure, n (%)',
  renal_failure='Renal failure, n (%)', longterm_insulin='Long-term use of insulin, n (%)',
  longterm_antidiabetics='Long-term use of antidiabetics, n (%)', TyG='TyG', SHR='SHR', AIP='AIP')
extended <- c(labels, sapsii='SAPS II', apsiii='APS III', oasis='OASIS',
  hospital_outcome='In-hospital mortality, n (%)', outcome_28d='28-day mortality, n (%)',
  icu_outcome='ICU mortality, n (%)')
binary <- c('gender','Glucocorticoids','Vasopressors','Ventilation','Hypertension','MI','CHF','CPD',
  'Diabetes','liver_failure','renal_failure','longterm_insulin','longterm_antidiabetics',
  'hospital_outcome','outcome_28d','icu_outcome')
write_result <- function(x, name) ir_write_csv(x, file.path(outdir, name))
bool <- function(x) {
  z <- rep(NA_integer_, length(x))
  z[tolower(as.character(x)) %in% c('1','yes','true')] <- 1L
  z[tolower(as.character(x)) %in% c('0','no','false')] <- 0L
  z
}
sources <- list(); field_audit <- list()
for (co in names(codes)) {
  raw <- read.csv(ir_input_path(co), check.names=FALSE)
  core <- ir_prepare_cohort(co)
  if (co!='MIMIC-IV') {
    db <- if (co=='MIMIC-III') 'mimiciii' else 'nwicu'
    extra <- read.csv(file.path(root, 'data/restricted', paste0(db, '_baseline_supplement.csv')))
    key <- if (co=='MIMIC-III') 'icustay_id' else 'stay_id'
    stopifnot(!anyDuplicated(extra[[key]]), setequal(raw[[key]], extra[[key]]))
    extra <- extra[match(raw[[key]],extra[[key]]),,drop=FALSE]
    stopifnot(identical(raw$subject_id, extra$subject_id), identical(raw$hadm_id, extra$hadm_id))
    for (v in setdiff(names(extra), c('subject_id','hadm_id',key))) raw[[v]] <- extra[[v]]
  }
  field_map <- list(gcs=c('gcs_min','mingcs','gcs'), sirs='sirs',
    HR=c('heart_rate_mean','HR'), RR=c('resp_rate_mean','RR'), MBP=c('mbp_mean','MBP'),
    Tm=c('temperature_mean','Tm'), WBC=c('wbc_mean','WBC'), PLT=c('platelets_mean','PLT'),
    Hb=c('hemoglobin_mean','Hb'), BUN=c('bun_mean','BUN'), Scr=c('creatinine_mean','Scr'),
    Na=c('sodium_mean','Na'), K=c('potassium_mean','K'), INR=c('inr_mean','INR'),
    PT=c('pt_mean','PT'), ALT=c('alt_mean','ALT'), AST=c('ast_mean','AST'),
    Glucocorticoids=c('glucocorticoids_1d','glucocorticoids'),
    Vasopressors=c('vasopressor_day1','vasopressor','vasopressor_icu'),
    Ventilation=c('ventilation_day1','ventilation','ventilation_icu'),
    Hypertension=c('hypertension','icd_hypertension'), MI='myocardial_infarct',
    CHF='congestive_heart_failure', CPD='chronic_pulmonary_disease',
    liver_failure='liver_failure', renal_failure=c('renal_failure','charlson_renal_disease'),
    longterm_insulin='longterm_insulin', longterm_antidiabetics='longterm_antidiabetics',
    sapsii='sapsii', apsiii='apsiii', oasis='oasis', hospital_outcome='hospital_outcome',
    outcome_28d=c('outcome_28d','mortality_28d'), icu_outcome='icu_outcome')
  d <- data.frame(gender=as.integer(core$gender=='Male'),age=core$age,
    los_hospital=as.numeric(difftime(ir_time(raw$dischtime),ir_time(raw$admittime),units='days')),
    los_icu=as.numeric(difftime(ir_time(raw$icu_outtime),ir_time(raw$icu_intime),units='days')),
    sofa=core$sofa,charlson=core$charlson,Diabetes=as.integer(core$diabetes=='Yes'),
    TG=core$tg,Glucose=core$glucose,HbA1c=core$hba1c,HDL=core$hdl,
    TyG=core$TyG,SHR=core$SHR,AIP=core$AIP)
  stopifnot(all(d$los_icu>=1),all(as.numeric(difftime(ir_time(raw$dischtime),ir_time(raw$icu_intime),units='days'))>=1))
  for (v in names(field_map)) {
    src <- ir_first_name(raw,field_map[[v]],required=FALSE)
    z <- if (is.na(src)) rep(NA_real_,nrow(raw)) else raw[[src]]
    d[[v]] <- if (v %in% binary) bool(z) else ir_num(z)
    field_audit[[length(field_audit)+1L]] <- data.frame(cohort=co,variable=v,
      source_field=if(is.na(src)) '' else src, source_available=!is.na(src))
  }
  if (co=='NWICU') {
    mild <- bool(raw$charlson_mild_liver_disease); severe <- bool(raw$charlson_severe_liver_disease)
    d$liver_failure <- ifelse(is.na(mild)&is.na(severe),NA_integer_,as.integer(mild==1|severe==1))
  }
  derived_fields <- c(gender='gender or sex; mapped to male indicator',
    age='age or admission_age; deidentified elderly ages coded as 90',
    los_hospital='dischtime minus admittime; exact days',
    los_icu='icu_outtime minus icu_intime; exact days',
    sofa='sofa or sofa_score; unchanged formal cohort value',
    charlson='charlson_comorbidity_index',Diabetes='diabetes or icd_diabetes',
    TG='tg, triglycerides, or first_day_triglycerides_mean',
    Glucose='glucose or first_day_glucose_mean',HbA1c='hba1c or first_day_hba1c_mean',
    HDL='hdl or first_day_hdl_mean',TyG='ln(TG * glucose / 2)',
    SHR='glucose / (28.7 * HbA1c - 46.7)',
    AIP='log10((TG / 88.57) / (HDL / 38.67))')
  for (v in names(derived_fields)) field_audit[[length(field_audit)+1L]] <- data.frame(
    cohort=co,variable=v,source_field=derived_fields[[v]],source_available=TRUE)
  d <- d[,names(extended),drop=FALSE]
  d$cohort <- co
  d$included <- complete.cases(d[,c('Glucose','HbA1c','TG','HDL','SHR','TyG','AIP')])
  stopifnot(identical(d$included,complete.cases(core[,c('SHR','TyG','AIP')])))
  sources[[co]] <- d
}
sources <- do.call(rbind,sources)
availability <- expand.grid(cohort=names(codes),variable=names(extended),stringsAsFactors=FALSE)
availability$status <- 'Observed source available'
availability$status[availability$cohort=='MIMIC-III' & availability$variable=='longterm_antidiabetics'] <-
  'Not available: no comparable dedicated ICD-9 long-term oral antidiabetic code'
availability$status[availability$cohort=='NWICU' & availability$variable %in% c('sapsii','apsiii','oasis')] <-
  'Not available: score not present in the NWICU extraction'
available <- function(co,v) availability$status[availability$cohort==co & availability$variable==v]=='Observed source available'
format_p <- function(p) if(!is.finite(p)) 'NA' else if(p<.001) '<0.001' else sprintf('%.3f',p)
summary_value <- function(z,v,ok=TRUE) {
  if(!ok) return('NA')
  z <- z[is.finite(z)]
  if(!length(z)) return('NA')
  if(v %in% binary) return(sprintf('%d (%.1f%%)',sum(z==1),100*mean(z==1)))
  digits <- if(v=='AIP') 2L else 1L
  q <- quantile(z,c(.5,.25,.75));q[abs(q)<.5*10^(-digits)] <- 0
  sprintf(paste0('%.',digits,'f (%.',digits,'f–%.',digits,'f)'),q[1],q[2],q[3])
}
compare <- function(z,g,v) {
  keep <- is.finite(z)&!is.na(g);z<-z[keep];g<-droplevels(factor(g[keep]))
  if(length(levels(g))<2 || length(unique(z))<2) return(list(p=NA_real_,test='Not estimable: insufficient variation'))
  if(!v %in% binary) {
    if(nlevels(g)==3) return(list(p=kruskal.test(z~g)$p.value,test='Kruskal-Wallis'))
    return(list(p=wilcox.test(z~g,exact=FALSE)$p.value,test='Mann-Whitney U'))
  }
  tab <- table(z,g); chi <- suppressWarnings(chisq.test(tab,correct=FALSE))
  if(any(chi$expected<5)) return(list(p=fisher.test(tab,workspace=2e7)$p.value,test='Fisher exact'))
  list(p=chi$p.value,test='Pearson chi-square')
}
inc <- sources[sources$included,,drop=FALSE]
m4 <- inc[inc$cohort=='MIMIC-IV',,drop=FALSE]
stopifnot(nrow(m4)==3056,!anyNA(m4$hospital_outcome),
  sum(m4$hospital_outcome==0)==2735,sum(m4$hospital_outcome==1)==321)
survival_baseline <- do.call(rbind,lapply(names(labels),function(v) {
  alive<-m4$hospital_outcome==0;dead<-m4$hospital_outcome==1
  test<-compare(m4[[v]],m4$hospital_outcome,v)
  data.frame(cohort='MIMIC-IV',variable=v,label=labels[[v]],
    type=if(v%in%binary)'binary' else 'continuous',
    overall=summary_value(m4[[v]],v),group_1=summary_value(m4[alive,v],v),
    group_2=summary_value(m4[dead,v],v),p_value=test$p,p_display=format_p(test$p),
    comparison_test=test$test,total_n=nrow(m4),group_1_n=sum(alive),group_2_n=sum(dead),
    overall_nonmissing_n=sum(!is.na(m4[[v]])),
    group_1_nonmissing_n=sum(alive&!is.na(m4[[v]])),
    group_2_nonmissing_n=sum(dead&!is.na(m4[[v]])),
    group_1_label='Hospital survivor',group_2_label='Hospital non-survivor')
}))
table1 <- do.call(rbind,lapply(names(labels),function(v) {
  a <- lapply(names(codes),function(co) summary_value(inc[inc$cohort==co,v],v,available(co,v)))
  ok <- vapply(inc$cohort,function(co) available(co,v),logical(1))
  test <- compare(inc[[v]][ok],inc$cohort[ok],v)
  data.frame(variable=v,label=labels[[v]],type=if(v%in%binary)'binary' else 'continuous',
    MIMICIV=a[[1]],MIMICIII=a[[2]],NWICU=a[[3]],p_value=test$p,p_display=format_p(test$p),
    comparison_test=test$test,comparison_scope=paste(names(codes)[vapply(names(codes),function(co)
      available(co,v)&&any(is.finite(inc[inc$cohort==co,v])),logical(1))],collapse=';'),
    MIMICIV_n=sum(inc$cohort=='MIMIC-IV'),MIMICIII_n=sum(inc$cohort=='MIMIC-III'),NWICU_n=sum(inc$cohort=='NWICU'))
}))
missingness <- do.call(rbind,lapply(names(codes),function(co) do.call(rbind,lapply(names(labels),function(v) {
  z <- inc[inc$cohort==co,v];ok<-available(co,v);n<-length(z);miss<-if(ok)sum(is.na(z))else NA_integer_
  data.frame(cohort=co,variable=v,label=labels[[v]],total_n=n,
    nonmissing_n=if(ok)n-miss else NA_integer_,missing_n=miss,
    missing_percent=if(ok)100*miss/n else NA_real_,
    display=if(ok)sprintf('%d (%.1f%%)',miss,100*miss/n)else 'NA',
    availability=availability$status[availability$cohort==co & availability$variable==v])
}))))
comparison <- do.call(rbind,lapply(names(codes),function(co) do.call(rbind,lapply(names(extended),function(v) {
  d<-sources[sources$cohort==co,,drop=FALSE];ok<-available(co,v)
  test<-if(ok)compare(d[[v]],d$included,v)else list(p=NA_real_,test='Not available')
  data.frame(cohort=co,variable=v,label=extended[[v]],type=if(v%in%binary)'binary' else 'continuous',
    group_1=summary_value(d[d$included,v],v,ok),group_2=summary_value(d[!d$included,v],v,ok),
    p_value=test$p,p_display=format_p(test$p),comparison_test=test$test,
    group_1_n=sum(d$included),group_2_n=sum(!d$included),
    group_1_nonmissing_n=if(ok)sum(d$included&!is.na(d[[v]]))else NA_integer_,
    group_2_nonmissing_n=if(ok)sum(!d$included&!is.na(d[[v]]))else NA_integer_)
}))))
counts<-do.call(rbind,lapply(names(codes),function(co) {
  d<-sources[sources$cohort==co,];data.frame(cohort=co,eligible_n=nrow(d),
    included_n=sum(d$included),excluded_n=sum(!d$included),included_percent=100*mean(d$included))
}))
stopifnot(counts$eligible_n==c(51564,32122,17351),counts$included_n==c(3056,1133,1556))
stopifnot(all(missingness$missing_n[missingness$variable%in%c('age','gender','sofa','charlson','Diabetes')]==0))
write_result(table1,'00_CR_ALL_table1_observed_baseline.csv')
write_result(survival_baseline,'00_CR_M4_survivor_nonsurvivor_baseline.csv')
write_result(missingness,'00_CR_ALL_Table1_variable_missingness.csv')
write_result(comparison,'00_CR_ALL_included_excluded_comparison.csv')
write_result(counts,'00_CR_ALL_included_excluded_counts.csv')
write_result(availability,'00_CR_ALL_baseline_field_availability.csv')
fields<-do.call(rbind,field_audit)
nw_liver<-fields$cohort=='NWICU' & fields$variable=='liver_failure'
fields$source_field[nw_liver]<-'charlson_mild_liver_disease or charlson_severe_liver_disease'
fields$source_available[nw_liver]<-TRUE
stopifnot(nrow(fields)==3*length(extended),!anyDuplicated(fields[,c('cohort','variable')]))
write_result(fields,'00_CR_ALL_baseline_source_fields.csv')
message('Observed baseline tables exported: 42 Table 1 rows; 126 missingness rows; 144 selection-comparison rows; 42 MIMIC-IV survival-status rows.')
