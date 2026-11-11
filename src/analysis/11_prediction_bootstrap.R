#!/usr/bin/env Rscript
# Refit paired repeated CV in each patient bootstrap; repeated copies share folds.
root<-Sys.getenv('IR_REVISION_ROOT')
source(file.path(root,'src/analysis/00_core.R'))
ir_check_required_packages('pROC')
# Load only data preparation, named variable sets and functions from shared prediction code.
e<-new.env(parent=globalenv())
wanted<-c('input','prepared','common','raw','getn','getb','dat','fill_prediction','folds_for','cv_pair','basic_vars','score_vars')
for(ex in parse(file.path(root,'src/shared/04_prediction.R'))) {
 if(is.call(ex)&&identical(ex[[1]],as.name('<-'))&&is.symbol(ex[[2]])&&as.character(ex[[2]])%in%wanted)eval(ex,e)
}
d<-e$dat
out<-file.path(root,'results/prediction/mimiciv')
roc<-coords<-delong<-list()
for(y in cr_outcomes) {
 rr<-list()
 for(x in c('SHR','TyG','AIP')) {
  r<-pROC::roc(d[[y]],d[[x]],direction='<',levels=0:1,quiet=TRUE,ci=TRUE);rr[[x]]<-r
  roc[[length(roc)+1]]<-data.frame(cohort='MIMIC-IV',outcome=y,exposure=x,n=nrow(d),events=sum(d[[y]]),auc=as.numeric(r$auc),conf_low=r$ci[1],conf_high=r$ci[3])
  a<-pROC::coords(r,x='all',ret=c('threshold','specificity','sensitivity'),transpose=FALSE)
  coords[[length(coords)+1]]<-data.frame(cohort='MIMIC-IV',outcome=y,exposure=x,a)
 }
 for(pair in combn(names(rr),2,simplify=FALSE))delong[[length(delong)+1]]<-data.frame(cohort='MIMIC-IV',outcome=y,exposure_1=pair[1],exposure_2=pair[2],p_value=pROC::roc.test(rr[[pair[1]]],rr[[pair[2]]],paired=TRUE,method='delong')$p.value)
}
ir_write_csv(do.call(rbind,roc),file.path(out,'01_PSA_M4_AO_NA_ALL_roc_auc.csv'))
ir_write_csv(do.call(rbind,coords),file.path(out,'01_PSA_M4_AO_NA_ALL_roc_coordinates.csv'))
ir_write_csv(do.call(rbind,delong),file.path(out,'01_PSA_M4_AO_NA_ALL_delong.csv'))

cv_all<-function(z,y,vars,index,seed) {
 # Compress bootstrap copies for fitting but expand for AUC evaluation.
 freq<-tabulate(index,nbins=nrow(z));ids<-which(freq>0);a<-z[ids,];w<-freq[ids]
 yy<-a[[y]];xs<-as.matrix(a[,vars,drop=FALSE]);xx<-as.matrix(a[,c('SHR','TyG','AIP')])
 au<-matrix(NA_real_,10,4,dimnames=list(NULL,c('Baseline','SHR','TyG','AIP')))
 for(r in 1:10) {
  fold<-e$folds_for(yy,5,seed+r);pred<-matrix(NA_real_,nrow(a),4)
  for(f in 1:5) {
   train<-fold!=f;test<-!train;filled<-xs
   for(v in seq_len(ncol(xs)))if(anyNA(xs[,v])) {
    vv<-rep(xs[train,v],w[train]);med<-median(vv,na.rm=TRUE)
    filled[is.na(filled[,v]),v]<-med
   }
   xb<-cbind(Intercept=1,filled)
   for(j in 1:4) {
    xm<-if(j==1)xb else cbind(xb,xx[,j-1])
    fit<-glm.fit(xm[train,,drop=FALSE],yy[train],weights=w[train],family=binomial())
    if(!fit$converged || any(!is.finite(fit$coefficients)))stop('CV fit failed or singular')
    pred[test,j]<-plogis(xm[test,,drop=FALSE]%*%fit$coefficients)
   }
  }
  expanded<-rep(seq_len(nrow(a)),w)
  for(j in 1:4)au[r,j]<-as.numeric(pROC::roc(yy[expanded],pred[expanded,j],direction='<',levels=0:1,quiet=TRUE)$auc)
 }
 delta<-sweep(au[,-1,drop=FALSE],1,au[,1],'-')
 list(base=mean(au[,1]),base_sd=sd(au[,1]),ext=colMeans(au[,-1]),ext_sd=apply(au[,-1],2,sd),
      delta=colMeans(delta),delta_sd=apply(delta,2,sd))
}
B<-200L;rows<-contrasts<-fails<-list()
for(y in cr_outcomes)for(base in c(names(e$score_vars),'Basic')) {
 vars<-if(base=='Basic')e$basic_vars[[y]]else e$score_vars[[base]]
 seed<-71000L+match(y,cr_outcomes)*100L+match(base,c(names(e$score_vars),'Basic'))
 point<-cv_all(d,y,vars,seq_len(nrow(d)),seed)
 set.seed(seed);samples<-replicate(B,sample.int(nrow(d),nrow(d),replace=TRUE),simplify=FALSE)
 boot<-parallel::mclapply(seq_len(B),function(b)tryCatch(list(value=cv_all(d,y,vars,samples[[b]],seed+b*1000L),error=''),error=function(e)list(value=NULL,error=conditionMessage(e))),mc.cores=4L,mc.set.seed=FALSE)
 ok<-!vapply(boot,function(a)nzchar(a$error),logical(1));bb<-do.call(rbind,lapply(boot[ok],function(a)a$value$delta))
 for(j in 1:3)rows[[length(rows)+1]]<-data.frame(cohort='MIMIC-IV',outcome=y,baseline=base,exposure=c('SHR','TyG','AIP')[j],
  n=nrow(d),events=sum(d[[y]]),baseline_auc=point$base,baseline_auc_sd=point$base_sd,
  extended_auc=point$ext[j],extended_auc_sd=point$ext_sd[j],delta_auc=point$delta[j],delta_auc_sd=point$delta_sd[j],
  delta_auc_low=quantile(bb[,j],.025,na.rm=TRUE),delta_auc_high=quantile(bb[,j],.975,na.rm=TRUE),
  bootstrap_planned=B,bootstrap_success=sum(ok),seed=seed,repeats=10,folds=5,
  selection_note=if(base=='Basic')'Variable selection preceded resampling; exploratory'else'Prespecified severity score')
 for(j in 2:3)contrasts[[length(contrasts)+1]]<-data.frame(cohort='MIMIC-IV',outcome=y,baseline=base,
  contrast=paste('SHR minus',c('SHR','TyG','AIP')[j]),estimate=point$delta[1]-point$delta[j],
  conf_low=quantile(bb[,1]-bb[,j],.025,na.rm=TRUE),conf_high=quantile(bb[,1]-bb[,j],.975,na.rm=TRUE),
  bootstrap_planned=B,bootstrap_success=sum(ok),seed=seed)
 if(any(!ok))for(b in which(!ok))fails[[length(fails)+1]]<-data.frame(outcome=y,baseline=base,replicate=b,error=boot[[b]]$error)
 message('CV bootstrap: ',y,' ',base,' ',sum(ok),'/',B)
 cr_write(do.call(rbind,rows),'12_CR_patient_bootstrap_CV_increments.csv')
 cr_write(do.call(rbind,contrasts),'12_CR_paired_CV_increment_contrasts.csv')
}
a<-do.call(rbind,rows)
ir_write_csv(a[a$baseline!='Basic',],file.path(out,'09_PSA_M4_AO_NA_ALL_score_comparison_cv_auc.csv'))
ir_write_csv(a[a$baseline=='Basic',],file.path(out,'14_PSA_M4_AO_NA_ALL_basic_model_cv_auc.csv'))
cr_write(if(length(fails))do.call(rbind,fails)else data.frame(outcome=character(),baseline=character(),replicate=integer(),error=character()),'12_CR_CV_bootstrap_failures.csv')
writeLines(c('Paired 200-patient-bootstrap repeated 10x5 CV complete',capture.output(sessionInfo())),file.path(crdir,'reports/internal/03e_CR_prediction_session.txt'))
