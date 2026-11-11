#!/usr/bin/env python3
"""竞争风险修订图形。沿用项目Python样式；CIF/KM的统计坐标由R生成。"""
from pathlib import Path
import os,shutil,importlib.util,runpy
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(os.environ['IR_REVISION_ROOT']);C=R/'08_CR_竞争风险修订';D=C/'03_汇总结果';O=Path(os.environ['IR_CR_STAGE']);F=O/'figures';S=F/'source_data'
F.mkdir(parents=True,exist_ok=True);S.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('oldfig',R/'01_最终投稿文件/00_build_submission_figures.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
old.FIG=F;old.SRC=S
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
                     'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':300})
def save_cns_figure(fig,filename):
 fig.savefig(f'{filename}.pdf',bbox_inches='tight',dpi=300)
 fig.savefig(f'{filename}.svg',bbox_inches='tight',dpi=300)
 fig.savefig(f'{filename}.png',bbox_inches='tight',dpi=300)
old.save_cns_figure=save_cns_figure
_original_ticks=old.set_readable_log_ticks
def ticks_including_null(ax,low,high,**kwargs):
    return _original_ticks(ax,min(low,1.0),max(high,1.0),**kwargs)
old.set_readable_log_ticks=ticks_including_null
plt.rcParams.update({'font.size':7,'axes.titlesize':7,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':7,'figure.titlesize':9})
EXP=['SHR','TyG','AIP'];COH=['MIMIC-IV','MIMIC-III','NWICU'];CODE={'MIMIC-IV':'M4','MIMIC-III':'M3','NWICU':'NW'}
Y=['hospital_mortality','mortality_28d','icu_mortality'];YN={'hospital_mortality':'In-hospital death · sHR','mortality_28d':'28-day death · HR','icu_mortality':'First ICU death · sHR'}
manifest=[]
def save(fig,name):
 old.save_cns_figure(fig,F/name);plt.close(fig);manifest.append(name)
def load(name,scale='per_unit',formal=True):
 a=pd.read_csv(D/name)
 if formal and 'formal_status' in a:a=a[a.formal_status=='formal']
 if scale and 'exposure_scale' in a:a=a[a.exposure_scale==scale]
 a=a.copy();a['outcome_label']=a.outcome.map(YN);a['cohort_label']=a.cohort
 return a

def forest(a,name,title,cols=('outcome_label',),height=110,model=True,catcol='model',ylabel='Effect ratio (95% CI)'):
 a=a.copy();a.to_csv(S/(name+'_source.csv'),index=False)
 old.forest_by_exposure(a,name,title,label_cols=cols,model_colors=model,model_col=catcol,effect_label=ylabel,width=183,height=height)
 manifest.append(name)

runpy.run_path(str(C/'02_代码/04b_CR_BUILD_FLOWCHART.py'), run_name='__main__')
manifest.append('Figure1_Cohort_flowchart')
pred=R/'03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果'
rc=pd.read_csv(pred/'01_PSA_M4_AO_NA_ALL_roc_coordinates.csv')
rc['false_positive_rate']=1-rc.specificity;rc['true_positive_rate']=rc.sensitivity
rc.to_csv(S/'Figure5_ROC_curve_coordinates.csv',index=False)
pd.read_csv(pred/'01_PSA_M4_AO_NA_ALL_roc_auc.csv').to_csv(S/'Figure5_ROC_AUC_source.csv',index=False)
old.build_roc();manifest.append('Figure5_ROC')
inc=pd.read_csv(D/'12_CR_patient_bootstrap_CV_increments.csv')
inc.to_csv(S/'Figure6_Increment_source.csv',index=False)
fig,axs=plt.subplots(1,3,figsize=(183/25.4,100/25.4))
bases=['SOFA','SAPS II','APS III','OASIS']
for ax,y in zip(axs,Y):
 for j,x in enumerate(EXP):
  a=inc[(inc.outcome==y)&(inc.exposure==x)&inc.baseline.isin(bases)].set_index('baseline').loc[bases]
  pos=np.arange(4)+(j-1)*.18
  ax.hlines(pos,a.delta_auc_low,a.delta_auc_high,color=old.EXPOSURE_COLORS[x],lw=1)
  ax.scatter(a.delta_auc,pos,color=old.EXPOSURE_COLORS[x],s=13,label=x)
 ax.axvline(0,color='#888888',ls='--',lw=.7);ax.set_yticks(range(4),bases);ax.invert_yaxis()
 ax.set_xlabel('ΔAUC (patient bootstrap 95% CI)');ax.set_title(YN[y].split(' · ')[0],fontsize=7)
axs[0].legend(fontsize=6);fig.tight_layout();save(fig,'Figure6_Incremental_performance')

primary=load('01_CR_M4_AO_M5_regression_results.csv')
forest(primary,'FigureS01_Primary_MIMICIV','MIMIC-IV mortality associations')
for i,job,title in [(1,'AO_M5','Primary admission-origin analysis'),(2,'L24_M5','24-hour landmark analysis')]:
 a=pd.concat([load(f'01_CR_{code}_{job}_regression_results.csv') for code in CODE.values()])
 forest(a,'Figure2_All_cohorts' if i==1 else 'FigureS02_L24_M5_all_cohorts',title,cols=('cohort_label','outcome_label'),height=185)
lab=pd.concat([load(f'06_CR_{c}_AO_M5P_ALL_lab_timing_results.csv',scale='per_sd') for c in CODE.values()]);lab['definition']=lab.lab_definition.map({'early_window_mean':'Window mean','first_available':'First available','closest_to_admission':'Closest to ICU'})
for i,co in enumerate(COH):
 forest(lab[lab.cohort==co],f'FigureS{i+3:02d}_Alternative_lab_{CODE[co]}',co+' · alternative laboratory timing',cols=('definition','outcome_label'),height=185)
ex=load('07_CR_ALL_AO_M5_exclusion_results.csv');ex['exclusion_label']=ex.exclusion.str.replace('exclude_','Exclude ').str.replace('_',' ')
for i,co in enumerate(COH):
 forest(ex[ex.cohort==co],f'FigureS{i+6:02d}_Exclusions_{CODE[co]}',co+' · comorbidity exclusions',cols=('exclusion_label','outcome_label'),height=210,model=False,catcol='exclusion_label')
# NWICU缺少院外死亡来源，28天记录死亡单列探索性图，不混入正式验证。
nw=[]
for job in ['AO_M5','L24_M5']:
 a=load(f'01_CR_NW_{job}_regression_results.csv',formal=False);a=a[a.outcome=='mortality_28d'].copy();a['setting']=job;nw.append(a)
forest(pd.concat(nw),'FigureS09_NW_recorded_28d_exploratory','NWICU · recorded 28-day mortality · exploratory',cols=('setting',),height=135,ylabel='Recorded-death HR (95% CI)')

# RCS 3×3：两列sHR，一列全因死亡HR，逐列明确效应名称。
rcs=load('03_CR_ALL_AO_M5_rcs_source_data.csv',scale=None);rcs=rcs[rcs.cohort=='MIMIC-IV'];rcs.to_csv(S/'Figure3_RCS_source.csv',index=False)
old.validate_effects(rcs,'RCS')
assert (rcs[['estimate','conf_low','conf_high']] > 0).all().all(), 'Log effect axes require positive values'
fig,axs=plt.subplots(3,3,figsize=(183/25.4,220/25.4))
for i,x in enumerate(EXP):
 for j,y in enumerate(Y):
  ax=axs[i,j];a=rcs[(rcs.exposure==x)&(rcs.outcome==y)].sort_values('x');color=old.EXPOSURE_COLORS[x]
  ax.fill_between(a.x,a.conf_low,a.conf_high,color=color,alpha=.18,lw=0);ax.plot(a.x,a.estimate,color=color,lw=1.3)
  ax.axhline(1,ls='--',color='#777777',lw=.6);ax.axvline(a.reference.iloc[0],ls=':',color='#888888',lw=.6)
  ax.set_yscale('log')
  low,high=ax.get_ylim()
  ticks=[t for t in [.1,.2,.3,.5,.75,1,1.5,2,3,5,7.5,10]if low<=t<=high]
  if len(ticks)>5:ticks=ticks[::2]if 1 in ticks[::2]else sorted(set(ticks[::2]+[1]))
  ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
  ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,pos:f'{v:g}'))
  ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
  ax.set_title(f'{chr(97+i*3+j)}  {x} · '+YN[y].split(' · ')[0],loc='left',fontweight='bold')
  ax.set_xlabel(x);ax.set_ylabel('HR' if y=='mortality_28d' else 'Subdistribution HR')
  ax.text(.03,.97,'Overall P '+old.p_text(a.p_overall.iloc[0])+'\nNonlinear P '+old.p_text(a.p_nonlinear.iloc[0]),transform=ax.transAxes,va='top',fontsize=6)
fig.suptitle('Dose–response associations · Model 3',fontweight='bold',y=.995);fig.tight_layout(rect=(0,0,1,.97));save(fig,'Figure3_RCS')
# CIF/KM不混用概率名称；风险表置于独立子网格，避免与坐标和曲线重叠。
coords=load('04_CR_ALL_AO_tertile_curve_coordinates.csv',None);tests=load('04_CR_ALL_AO_Gray_logrank_tests.csv',None);risk=load('04_CR_ALL_AO_fixed_time_risk_counts.csv',None)
coords=coords[coords.cohort=='MIMIC-IV'];tests=tests[tests.cohort=='MIMIC-IV'];risk=risk[risk.cohort=='MIMIC-IV']
for a,n in [(coords,'curve_coordinates'),(tests,'Gray_logrank_tests'),(risk,'risk_counts')]:a.to_csv(S/('Figure4_'+n+'.csv'),index=False)
fig=plt.figure(figsize=(183/25.4,238/25.4));outer=fig.add_gridspec(3,3,hspace=.43,wspace=.45,left=.09,right=.98,bottom=.035,top=.93)
for i,x in enumerate(EXP):
 for j,y in enumerate(Y):
  sub=outer[i,j].subgridspec(2,1,height_ratios=[3.4,1.4],hspace=.62);ax=fig.add_subplot(sub[0]);rtax=fig.add_subplot(sub[1]);rtax.axis('off')
  a=coords[(coords.exposure==x)&(coords.outcome==y)];test=tests[(tests.exposure==x)&(tests.outcome==y)].iloc[0]
  for g in ['T1','T2','T3']:
   b=a[a.tertile==g].sort_values('time',kind='stable');color=old.TERTILE_COLORS[g]
   ax.step(b.time,b.probability,where='post',color=color,lw=1.05)
   ax.fill_between(b.time,b.conf_low,b.conf_high,step='post',color=color,alpha=.09,lw=0)
  ax.set_xlim(0,28);ax.set_xticks([0,7,14,21,28]);ax.set_xlabel('Days from ICU admission')
  ax.set_title(f'{chr(97+i*3+j)}  {x} · '+YN[y].split(' · ')[0],loc='left',fontsize=6.5,fontweight='bold')
  ax.set_ylabel('Survival probability' if y=='mortality_28d' else 'Cumulative incidence')
  if y=='mortality_28d':ax.set_ylim(.65,1.005)
  else:ax.set_ylim(0,.27 if y=='hospital_mortality' else .18)
  ax.text(.03,.06 if y=='mortality_28d' else .94,test['test']+' P '+old.p_text(test.p_value),transform=ax.transAxes,va='bottom' if y=='mortality_28d' else 'top',fontsize=6)
  r=risk[(risk.exposure==x)&(risk.outcome==y)]
  rtax.text(-.08,.97,'At risk',ha='right',fontsize=5.5,transform=rtax.transAxes)
  for t in [0,7,14,21,28]:rtax.text(t/28,.97,str(t),ha='center',fontsize=5.5,transform=rtax.transAxes)
  for k,g in enumerate(['T1','T2','T3']):
   yy=.68-k*.29;rtax.text(-.08,yy,g,ha='right',color=old.TERTILE_COLORS[g],fontsize=5.5,transform=rtax.transAxes)
   for t in [0,7,14,21,28]:rtax.text(t/28,yy,str(int(r[(r.tertile==g)&(r.time==t)].n_risk.iloc[0])),ha='center',fontsize=5.5,transform=rtax.transAxes)
fig.legend(handles=[Line2D([0],[0],color=old.TERTILE_COLORS[g],label=g)for g in ['T1','T2','T3']],ncol=3,loc='upper center',bbox_to_anchor=(.5,.977))
fig.suptitle('Cumulative incidence and 28-day survival by index tertile',fontweight='bold',y=.998);save(fig,'Figure4_CIF_KM')
# Use the same dedicated entry for full builds and focused subgroup updates.
runpy.run_path(str(C/'02_代码/04a_CR_BUILD_SUBGROUP_FIGURE.py'))
manifest.append('FigureS10_Subgroup_forest')
pd.DataFrame({'figure':manifest,'png_exists':[(F/(n+'.png')).exists()for n in manifest],'pdf_exists':[(F/(n+'.pdf')).exists()for n in manifest]}).to_csv(F/'figure_export_manifest.csv',index=False)
print('Created',len(manifest),'figure groups')

# The staging directory holds only the current, numbered figure set.
allowed=set(manifest)
for path in F.glob('Figure*'):
 if path.suffix in {'.pdf','.svg','.png'} and path.stem not in allowed:
  path.unlink()
