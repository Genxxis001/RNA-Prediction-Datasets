"""Writes the LaTeX code of Tables 3, 4, 5, 7, S1, S3 and S4 (longtable format used in main.tex).

Needs the CSVs from compute_test_metrics.py and compute_extra_metrics.py.
Tables 1, 2, 6 and S2 are descriptive and were typed by hand; their numbers come from
dataset_summary.py (Tables 1/S2) and fig4_to_fig10_sequence_logos (Table 6).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
import pandas as pd,numpy as np
R=pd.read_csv(config.OUT_RES/'final_test05.csv');T=pd.read_csv(config.OUT_RES/'human_train_val_metrics.csv');DL=pd.read_csv(config.OUT_RES/'delong.csv')
M=pd.read_csv(config.OUT_RES/'mouse_metrics.csv');F=pd.read_csv(config.OUT_RES/'human_full05.csv')
ORD=['m6A','m1A','m5C','m7G','Pseudo','Nm']
MOD={'m6A':'\\mSixA{}','m1A':'\\mOneA{}','m5C':'\\mFiveC{}','m7G':'\\mSevenG{}','Pseudo':'\\Psu{}','Nm':'Nm'}  # macros defined in main.tex preamble
NM={'deepSRAMP-seqonly41':'deepSRAMP','CLSM6A-41nt':'CLSM6A','RNA_m5C_predict':'RNA\\_m5C\\_\\allowbreak{}predict'}
REF={'MKE-ResNet':'mke','m6A-IIN':'iin','BiRNA-BERT':'birna','CLSM6A-41nt':'clsm6a','deepSRAMP-seqonly41':'deepsramp','EvoRMD':'evormd','m1A-Ensem':'m1aensem','BERT-RNA':'bertrna','MRM-BERT':'mrmbert','FourierAN-m1A':'fourieran','m5C-pred':'m5cpred','RNA_m5C_predict':'rnam5c','Trans-m5C':'transm5c','m5C-Seq':'m5cseq','FusDRM-m5C':'fusdrm','Moss-m7G':'moss','CAP-m7G':'capm7g','BF-m7GPred':'bfm7g','MCAMEF-BERT':'mcamef','autoBioSeqpy':'autobio','Definer-PseU':'definer','MetaPseU':'metapseu','PseU-FKeERF':'fkeerf','PseU-ST':'pseust','RSCNN-PseU':'rscnn','Nmix':'nmix','2OMe-LM':'omelm','BERT2OME':'bert2ome','Meta-2OM':'meta2om','2OM-Pred':'ompred'}
n=lambda m:NM.get(m,m)
f=lambda v,d=3:(f"{v:.{d}f}").replace('-','$-$') if v<0 else f"{v:.{d}f}"
def lt(spec,cap,lab,head,rows,size='\\scriptsize\\setlength{\\tabcolsep}{3pt}'):
    nc=head.count('&')+1
    s=f"{{{size}\n\\begin{{longtable}}[]{{{spec}}}\n\\caption{{{cap}}}\\label{{{lab}}}\\\\\n\\toprule\n{head}\\\\\n\\midrule\n\\endfirsthead\n"
    s+=f"\\multicolumn{{{nc}}}{{@{{}}l}}{{\\small\\textit{{{lab.replace('tab:','Table ')} (continued)}}}}\\\\\n\\toprule\n{head}\\\\\n\\midrule\n\\endhead\n\\bottomrule\n\\endlastfoot\n"
    return s+"\n".join(rows)+"\n\\end{longtable}\n}\n"
P=lambda w:'>{\\raggedright\\arraybackslash}p{'+w+'}'
out={}
# Table 3 human test
H=R[R.sp=='human'];rows=[];prev=None
for mod in ORD:
    g=H[H['mod']==mod].sort_values('MCC',ascending=False);best={c:g[c].max() for c in ['Acc','Pre','Sn','Sp','F1','MCC','AUC']}
    if prev: rows.append('\\midrule')
    prev=1
    for _,r in g.iterrows():
        B=lambda c,txt: f"\\textbf{{{txt}}}" if abs(r[c]-best[c])<1e-12 else txt
        rows.append(f"{MOD[mod]} & {n(r.method)} & {B('Acc',f(r.Acc,2))} & {B('Pre',f(r.Pre,2))} & {B('Sn',f(r.Sn,2))} & {B('Sp',f(r.Sp,2))} & {B('F1',f(r.F1,2))} & {B('MCC',f(r.MCC))} [{f(r.MCC_lo)}, {f(r.MCC_hi)}] & {B('AUC',f(r.AUC))} \\\\")
out['t3']=lt('@{}l'+P('2.4cm')+'rrrrrlr@{}','\\textbf{Table 3.} Held-out test performance of the 30 retrained predictors (threshold 0.5). Rows are sorted by MCC within each modification; brackets give 95\\% bootstrap confidence intervals. Acc, Pre, Sn, Sp and F1 in \\%. Bold marks the best value per column within each modification.','tab:3','Modification & Method & Acc & Pre & Sn & Sp & F1 & MCC [95\\% CI] & AUROC',rows)
# Table 4 summary
H2=H.merge(T[['mod','method','train_AUC']],on=['mod','method']).merge(F[['mod','method','MCC']].rename(columns={'MCC':'fMCC'}),on=['mod','method'])
rows=[]
for mod in ORD:
    g=H2[H2['mod']==mod];tb=g.loc[g.MCC.idxmax()];ab=g.loc[g.AUC.idxmax()];fb=g.loc[g.fMCC.idxmax()]
    rk=int(g.MCC.rank(ascending=False,method='min')[fb.name]);d=DL[(DL.sp=='human')&(DL['mod']==mod)];ns=(d.p_holm<0.05).sum()
    N=int(tb.TP+tb.FN+tb.TN+tb.FP)
    rows.append(f"{MOD[mod]} & {N:,} & {n(tb.method)} ({f(tb.MCC)}) & {n(ab.method)} ({f(ab.AUC)}) & {f(g.MCC.mean())} $\\pm$ {f(g.MCC.std())} & {n(fb.method)} $\\rightarrow$ {rk} & {ns} of 4 & {f((g.train_AUC-g.AUC).mean())} \\\\")
out['t4']=lt('@{}lr'+P('2.2cm')+P('2.2cm')+'l'+P('2.3cm')+P('1.9cm')+P('1.5cm')+'@{}','\\textbf{Table 4.} Summary of held-out results per modification. ``Full-data best $\\rightarrow$ test rank\'\' gives the test MCC rank of the method that ranked first on the full dataset. ``Significantly below best AUROC\'\' counts competitors whose test AUROC is significantly lower than that of the best method (paired DeLong test, Holm-corrected $p<0.05$). The last column is the mean difference between training and test AUROC.','tab:4','Modification & Test records & Best MCC & Best AUROC & Mean MCC $\\pm$ SD & Full-data best $\\rightarrow$ test rank & Significantly below best AUROC & Train $-$ test AUROC',rows)
# Table 5 bias
H['bias']=H.Sn-H.Sp;b=H[abs(H.bias)>=10].sort_values('bias',ascending=False);rows=[]
for _,r in b.iterrows(): rows.append(f"{MOD[r['mod']]} & {n(r.method)} & {int(r.TP):,} & {int(r.FN):,} & {int(r.TN):,} & {int(r.FP):,} & {'+' if r.bias>0 else '$-$'}{abs(r.bias):.1f} \\\\")
out['t5']=lt('@{}llrrrrr@{}','\\textbf{Table 5.} Test-set confusion matrices for predictors whose sensitivity and specificity differ by at least 10 percentage points.','tab:5','Modification & Method & TP & FN & TN & FP & Sn $-$ Sp (points)',rows,'\\footnotesize')
# Table 7 mouse
Mo=R[R.sp=='mouse'].merge(M[['mod','method','train_AUC']],on=['mod','method']);rows=[];prev=None
for mod in ['m6A','m5C','Pseudo']:
    g=Mo[Mo['mod']==mod].sort_values('MCC',ascending=False);best={c:g[c].max() for c in ['MCC','AUC']}
    if prev: rows.append('\\midrule')
    prev=1
    for _,r in g.iterrows():
        B=lambda c,txt: f"\\textbf{{{txt}}}" if abs(r[c]-best[c])<1e-12 else txt
        rows.append(f"{MOD[mod]} & {n(r.method)} & {f(r.Sn,2)} & {f(r.Sp,2)} & {B('MCC',f(r.MCC))} [{f(r.MCC_lo)}, {f(r.MCC_hi)}] & {B('AUC',f(r.AUC))} & {f(r.train_AUC-r.AUC)} \\\\")
out['t7']=lt('@{}l'+P('2.4cm')+'rrlrr@{}','\\textbf{Table 7.} Held-out test performance on the three full-size mouse datasets (threshold 0.5; test sets of 4,000, 4,000 and 972 records). Sn and Sp in \\%; brackets give 95\\% bootstrap confidence intervals. Bold marks the best MCC and AUROC within each modification.','tab:7','Modification & Method & Sn & Sp & MCC [95\\% CI] & AUROC & Train $-$ test AUROC',rows,'\\footnotesize')
# Table S1 full data
rows=[];prev=None
for mod in ORD:
    g=F[F['mod']==mod].sort_values('MCC',ascending=False)
    if prev: rows.append('\\midrule')
    prev=1
    for _,r in g.iterrows(): rows.append(f"{MOD[mod]} & {n(r.method)} & {f(r.Acc,2)} & {f(r.Pre,2)} & {f(r.Sn,2)} & {f(r.Sp,2)} & {f(r.F1,2)} & {f(r.MCC)} & {f(r.AUC)} \\\\")
out['ts1']=lt('@{}l'+P('2.6cm')+'rrrrrrr@{}','\\textbf{Table S1.} Full-dataset performance (training, validation and test records combined; threshold 0.5). These values include records used for training and are shown only for comparison with Table 3.','tab:S1','Modification & Method & Acc & Pre & Sn & Sp & F1 & MCC & AUROC',rows)
# Table S3 pilot mouse
rows=[];prev=None
for mod in ['m1A','m7G','Nm']:
    g=M[M['mod']==mod].sort_values('full_MCC',ascending=False)
    if prev: rows.append('\\midrule')
    prev=1
    for _,r in g.iterrows(): rows.append(f"{MOD[mod]} & {n(r.method)} & {f(r.full_MCC)} & {f(r.train_AUC)} & {f(r.T_MCC)} & {f(r.T_AUC)} & {r.T_Sn:.1f} & {r.T_Sp:.1f} \\\\")
out['ts3']=lt('@{}l'+P('2.6cm')+'rrrrrr@{}','\\textbf{Table S3.} Pilot mouse datasets (m\\textsuperscript{1}A, m\\textsuperscript{7}G and Nm; test sets of 26, 6 and 26 records). Values are reported for completeness and are not used for ranking.','tab:S3','Modification & Method & Full-data MCC & Training AUROC & Test MCC & Test AUROC & Test Sn & Test Sp',rows,'\\footnotesize')
# Table S4 Psi validation thresholds (human)
P2=T[T['mod']=='Pseudo'].sort_values('Tv_MCC',ascending=False);rows=[]
for _,r in P2.iterrows(): rows.append(f"{n(r.method)} & {r.val_thr:.3f} & {f(r.Tv_Sn,2)} & {f(r.Tv_Sp,2)} & {f(r.Tv_MCC)} & {f(H[(H['mod']=='Pseudo')&(H.method==r.method)].MCC.iloc[0])} \\\\")
out['ts4']=lt('@{}lrrrrr@{}','\\textbf{Table S4.} Human $\\Psi$ test results when each predictor\'s threshold is chosen by maximum MCC on the validation partition, compared with the common threshold of 0.5 used in the main text.','tab:S4','Method & Threshold & Sn & Sp & MCC (validation threshold) & MCC (0.5)',rows,'\\footnotesize')
import json;json.dump(out,open(config.OUT_RES/'tables.json','w'))
NAMES={'t3':'table3_human_test','t4':'table4_summary','t5':'table5_bias_confusion','t7':'table7_mouse_test',
       'ts1':'tableS1_full_data','ts3':'tableS3_pilot_mouse','ts4':'tableS4_psi_val_thresholds'}
for k,v in out.items(): open(config.OUT_TAB/f"{NAMES[k]}.tex",'w').write(v)
print('LaTeX tables written:',', '.join(NAMES.values()))
