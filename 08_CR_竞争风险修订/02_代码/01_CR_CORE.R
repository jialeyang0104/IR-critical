# 竞争风险修订内核；旧Cox入口和结果原样保留。
root <- Sys.getenv('IR_REVISION_ROOT')
if (!nzchar(root)) stop('必须设置IR_REVISION_ROOT。')
source(file.path(root,'00_总说明与运行顺序/18_ANALYSIS_CORE.R'))
ir_check_required_packages(c('cmprsk','survival','mice','Hmisc'))
crdir <- file.path(root,'08_CR_竞争风险修订')
resdir <- file.path(crdir,'03_汇总结果')
dir.create(resdir,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(crdir,'07_运行日志'),recursive=TRUE,showWarnings=FALSE)
cr_outcomes <- c('hospital_mortality','mortality_28d','icu_mortality')
cr_codes <- c('MIMIC-IV'='M4','MIMIC-III'='M3','NWICU'='NW')
cr_seeds <- c('MIMIC-IV'=31001L,'MIMIC-III'=32001L,'NWICU'=33001L)
cr_write <- function(x,name) {
  if(all(c('cohort','outcome')%in%names(x))) {
    limited<-x$cohort=='NWICU' & x$outcome=='mortality_28d'
    x$endpoint_scope<-ifelse(limited,'28-day deaths recorded in NM hospital data',
      ifelse(x$outcome=='mortality_28d','28-day all-cause mortality',
        ifelse(x$outcome=='hospital_mortality','Index-hospitalization mortality','First-ICU-stay mortality')))
    x$formal_status<-ifelse(limited,'exploratory_limited_ascertainment','formal')
  }
  ir_write_csv(x,file.path(resdir,name))
}
cr_measure <- function(y) if(y=='mortality_28d') 'HR' else 'sHR'
cr_attach_states <- function(d,raw) {
  stopifnot(nrow(d)==nrow(raw))
  # 0是真正右删失，1是目标死亡，2是存活离开本次住院/首次ICU。
  # 当前提取没有独立的最后联系时点；非死亡且缺少结束时点时停止，避免伪造竞争事件。
  for (y in c('hospital_mortality','icu_mortality')) {
    def <- ir_outcomes[[y]]; ek <- if(y=='hospital_mortality') 'dischtime' else 'icu_outtime'
    ek <- ir_first_name(raw,if(ek=='icu_outtime') c('icu_outtime','outtime') else c('dischtime','hospital_dischtime'))
    flags <- ir_num(ir_get(raw,if(y=='hospital_mortality') c('hospital_outcome','hospital_mortality') else c('icu_outcome','icu_mortality')))
    if(any(!is.na(flags)&!flags%in%c(0,1))) stop('非法结局编码。')
    end <- ir_time(raw[[ek]])
    state <- rep(NA_integer_,nrow(d)); state[!is.na(flags)&flags==1] <- 1L
    state[!is.na(flags)&flags==0 & !is.na(end)] <- 2L
    d[[paste0(y,'_status')]] <- state
  }
  d
}
cr_prepare <- function(cohort,path=ir_input_path(cohort)) {
  cr_attach_states(ir_prepare_cohort(cohort,path),read.csv(path,check.names=FALSE))
}
cr_risk <- function(d,y,landmark=0) {
  z <- ir_risk_set(d,y,landmark)
  if(y!='mortality_28d') {
    z$.status <- z[[paste0(y,'_status')]]
    if(anyNA(z$.status)) stop('竞争事件/删失状态未能从源字段确定：',unique(d$cohort),'/',y)
    stopifnot(all((z$.status==1)==(z$.event==1)))
  } else z$.status <- z$.event
  z
}
cr_apply_mi <- function(sets,fn,vars) {
  same <- all(vapply(sets,function(z) identical(z[,vars,drop=FALSE],sets[[1]][,vars,drop=FALSE]),logical(1)))
  if(same) rep(list(fn(sets[[1]])),length(sets)) else lapply(sets,fn)
}
cr_fit <- function(z,y,rhs,covars,weighted=FALSE,force_cause_specific=FALSE) {
  form <- as.formula(paste('~',paste(c(rhs,covars),collapse=' + ')))
  if(y=='mortality_28d'||force_cause_specific) {
    f <- survival::coxph(as.formula(paste('survival::Surv(.time,.event)',paste(deparse(form),collapse=' '))),
      z,ties='efron',weights=if(weighted) .w else rep(1,nrow(z)),robust=weighted,x=TRUE)
    b<-coef(f);v<-vcov(f); fit_type<-'cause_specific_cox'
  } else if(!weighted) {
    mm <- model.matrix(form,z)[,-1,drop=FALSE]
    f <- cmprsk::crr(z$.time,z$.status,cov1=mm,failcode=1,cencode=0)
    if(!isTRUE(f$converged)) stop('Fine-Gray模型未收敛。')
    b<-f$coef;v<-f$var;dimnames(v)<-list(names(b),names(b)); fit_type<-'fine_gray'
  } else {
    z$.id <- z$source_row
    z$.state <- factor(z$.status,levels=0:2,labels=c('censor','death','alive_exit'))
    fg <- survival::finegray(as.formula(paste('survival::Surv(.time,.state) ~',paste(c(rhs,covars,'.id'),collapse=' + '))),
                             z,weights=z$.w,etype='death',id=z$.id)
    f <- survival::coxph(as.formula(paste('survival::Surv(fgstart,fgstop,fgstatus) ~',paste(c(rhs,covars),collapse=' + '))),
                        fg,weights=fgwt,cluster=.id,robust=TRUE,ties='breslow')
    b<-coef(f);v<-vcov(f);fit_type<-'weighted_fine_gray'
  }
  if(any(!is.finite(b))||any(!is.finite(v))) stop('模型产生非有限系数/方差。')
  list(b=b,v=v,n=nrow(z),events=sum(z$.event),competing=sum(z$.status==2),censored=sum(z$.status==0),
       k=length(b),type=fit_type,fit=f)
}
cr_pool <- function(fits,term) {
  stopifnot(length(unique(vapply(fits,`[[`,numeric(1),'n')))==1,
            length(unique(vapply(fits,`[[`,numeric(1),'events')))==1)
  ir_pool_scalar(vapply(fits,function(f) unname(f$b[term]),numeric(1)),
                 vapply(fits,function(f) unname(f$v[term,term]),numeric(1)),fits[[1]]$events,fits[[1]]$k)
}
cr_row <- function(fits,term,cohort,x,y,model,scale,origin,m,level='PA',extra=list()) {
  p<-cr_pool(fits,term);f<-fits[[1]]
  complete <- is.infinite(p['df'])
  row<-data.frame(analysis_level=level,cohort=cohort,time_origin=origin,imputation=if(complete)'Complete covariates; no imputation'else paste0('M',m,' Rubin pooled'),
    exposure=x,exposure_scale=scale,outcome=y,model=model,term=term,effect_measure=cr_measure(y),
    n=f$n,events=f$events,competing_events=f$competing,censored=f$censored,
    estimate=exp(p['beta']),conf_low=exp(p['low']),conf_high=exp(p['high']),p_value=p['p'],
    regression_method=f$type,rubin_df=p['df'],fraction_missing_information=p['fmi'],m=if(complete)0 else m,maxit=if(complete)0 else 10,
    sample_definition='SHR、TyG和AIP共同完整样本',stringsAsFactors=FALSE)
  for(nm in names(extra)) row[[nm]]<-extra[[nm]]
  row
}
cr_matrix <- function(d,m,origin,scales,level,seed) {
  d<-ir_common_exposure_sample(d);sets<-ir_make_imputations(d,m,seed)
  lm<-if(origin=='Landmark24h')1 else 0;rows<-list()
  for(x in c('SHR','TyG','AIP')) {
    cuts<-quantile(d[[x]],c(0,1/3,2/3,1));cuts[c(1,4)]<-c(-Inf,Inf);xm<-mean(d[[x]]);xs<-sd(d[[x]])
    for(y in cr_outcomes) for(model in names(ir_model_covariates)) for(sc in scales) {
      covars<-ir_model_covariates[[model]]
      fn<-function(dat) {
        z<-cr_risk(dat,y,lm)
        if(sc=='per_unit') {z$.x<-z[[x]];rhs<-'.x'}
        if(sc=='per_sd') {z$.x<-(z[[x]]-xm)/xs;rhs<-'.x'}
        if(sc=='tertile') {z$.tertile<-cut(z[[x]],cuts,labels=c('T1','T2','T3'),include.lowest=TRUE);rhs<-'.tertile'}
        if(sc=='p_trend') {z$.trend<-as.integer(cut(z[[x]],cuts,include.lowest=TRUE));rhs<-'.trend'}
        cr_fit(z,y,rhs,covars)
      }
      fits<-cr_apply_mi(sets,fn,covars)
      tt<-switch(sc,per_unit='.x',per_sd='.x',tertile=c('.tertileT2','.tertileT3'),p_trend='.trend')
      for(t in tt) rows[[length(rows)+1]]<-cr_row(fits,t,unique(d$cohort),x,y,model,sc,origin,m,level)
    }
  }
  do.call(rbind,rows)
}
