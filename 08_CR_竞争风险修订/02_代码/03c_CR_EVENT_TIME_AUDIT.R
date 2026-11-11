#!/usr/bin/env Rscript
root<-Sys.getenv('IR_REVISION_ROOT');source(file.path(root,'08_CR_竞争风险修订/02_代码/01_CR_CORE.R'))
ans<-list()
for(cohort in names(cr_codes)) {
 p<-ir_input_path(cohort);raw<-read.csv(p,check.names=FALSE);d<-ir_common_exposure_sample(cr_prepare(cohort,p));r<-raw[d$source_row,]
 dt<-ir_time(ir_get(r,'deathtime',FALSE));src<-ir_get(r,c('mortality_28d_time_source','death_time_source'),FALSE)
 ans[[length(ans)+1]]<-data.frame(cohort=cohort,n=nrow(d),hospital_deaths=sum(d$hospital_outcome==1),
   hospital_deaths_without_exact_timestamp=sum(d$hospital_outcome==1&is.na(dt)),icu_deaths=sum(d$icu_outcome==1),
   icu_deaths_without_exact_timestamp=sum(d$icu_outcome==1&is.na(dt)),deaths28=sum(d$outcome_28d==1),
   deaths28_date_only=sum(d$outcome_28d==1&grepl('dod_date',src)),
   deaths28_hospital_timestamp=sum(d$outcome_28d==1&grepl('deathtime',src)))
}
ir_write_csv(do.call(rbind,ans),file.path(resdir,'00_CR_event_time_precision_audit.csv'))
risks<-list()
for(cohort in names(cr_codes)) {
 p<-ir_input_path(cohort);raw<-read.csv(p,check.names=FALSE);d<-ir_common_exposure_sample(cr_prepare(cohort,p));r<-raw[d$source_row,]
 src<-ir_get(r,c('mortality_28d_time_source','death_time_source'),FALSE)
 for(y in cr_outcomes)for(lm in c(0,1)) {
  keys<-ir_outcomes[[y]];tt<-d[[keys[1]]];ev<-d[[keys[2]]];included<-is.finite(tt)&tt>lm&!is.na(ev)
  risks[[length(risks)+1]]<-data.frame(cohort=cohort,outcome=y,time_origin=if(lm==0)'AdmissionOrigin'else'Landmark24h',
    common_n=nrow(d),included_n=sum(included),excluded_n=sum(!included),excluded_events=sum(ev[!included]==1,na.rm=TRUE),
    excluded_date_only_deaths=if(y=='mortality_28d')sum(!included&ev==1&grepl('dod_date',src),na.rm=TRUE)else 0)
 }
}
ir_write_csv(do.call(rbind,risks),file.path(resdir,'00_CR_risk_set_audit.csv'))
