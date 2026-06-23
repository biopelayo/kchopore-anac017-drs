#!/usr/bin/env python3
# =============================================================
# K-CHOPORE | m6A epitranscriptome from ELIGOS2 (pair_diff vs m6A writers)
# v2 publication/poster grade: English, large axis labels, re-centred motif,
# COVERAGE-NORMALISED rates, gene heatmap. Design 2x2: WT/anac017-1 x C/AA(Antimycin A).
# =============================================================
import os, re, glob
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

RES=os.path.expanduser("~/eligos_results"); OUT=os.path.expanduser("~/eligos_figs_v2"); os.makedirs(OUT,exist_ok=True)
CONDS=["WT_C","WT_AA","anac017-1_C","anac017-1_AA"]
LABEL={"WT_C":"WT · Ctrl","WT_AA":"WT · AA","anac017-1_C":"anac017 · Ctrl","anac017-1_AA":"anac017 · AA"}
WRITERS=["mta","mtb","fip37","hakai","vir"]; MIN_W=3
PAL={"WT_C":"#4C78A8","WT_AA":"#F58518","anac017-1_C":"#54A24B","anac017-1_AA":"#E45756"}
USECOLS=["chrom","start_loc","strand","name","ref","kmer7","oddR","pval","ESB_test","ESB_ctrl"]
ATRE=re.compile(r"AT[1-5CM]G\d{5}",re.I)

sns.set_theme(context="talk",style="ticks")
mpl.rcParams.update({"figure.dpi":120,"savefig.dpi":350,"savefig.bbox":"tight",
 "font.family":"DejaVu Sans","font.size":16,"axes.titlesize":21,"axes.titleweight":"bold",
 "axes.labelsize":19,"axes.labelweight":"bold","xtick.labelsize":15,"ytick.labelsize":15,
 "legend.fontsize":14,"axes.spines.top":False,"axes.spines.right":False,
 "axes.grid":True,"grid.alpha":0.25,"axes.axisbelow":True,"axes.linewidth":1.3})
def save(fig,n):
    for e in ("png","pdf"): fig.savefig(f"{OUT}/{n}.{e}")
    plt.close(fig); print("  fig:",n,flush=True)
def ffile(cond,w):
    g=[x for x in glob.glob(f"{RES}/{cond}_results_{w}/{w}_sorted_merged_vs_{cond}_sorted_merged_*_baseExt0.txt") if "._" not in os.path.basename(x)]
    return g[0] if g else None

def recentre(kmers):
    """align each 7-mer on the A of the first RAC/AC so the modified A is centred."""
    out=[]
    for k in kmers:
        if not isinstance(k,str): continue
        k=k.replace("U","T")
        idx=[m.start() for m in re.finditer("AC",k)]
        if not idx: continue
        a=idx[0]                      # modified A position
        win=k[max(0,a-2):a+3]
        if len(win)==5: out.append(win)
    return out

print("== load + consensus ==",flush=True)
percond={}; A_tested={}; raw={}
for cond in CONDS:
    cnt={}; esb={}; km={}; gn={}
    a_tested=0
    for w in WRITERS:
        f=ffile(cond,w)
        if not f: continue
        df=pd.read_csv(f,sep="\t",usecols=USECOLS,dtype={"chrom":str,"strand":str},low_memory=False)
        df=df[df["ref"].astype(str).str.upper()=="A"]
        for c in ("oddR","pval","ESB_test","ESB_ctrl"): df[c]=pd.to_numeric(df[c],errors="coerce")
        df=df.dropna(subset=["oddR","pval","ESB_test","ESB_ctrl"])
        if w=="mta": a_tested=len(df)              # coverage proxy: A interrogated
        sig=df[(df.pval<0.05)&(df.oddR<1)&(df.ESB_ctrl>df.ESB_test)].copy()
        sig["site"]=sig.chrom+":"+sig.start_loc.astype(str)+":"+sig.strand
        sig["gene"]=sig.name.astype(str).str.extract("("+ATRE.pattern+")",expand=False)
        if cond not in raw: raw[cond]=df[["oddR","pval","ESB_ctrl","ESB_test"]]
        for _,r in sig[["site","ESB_ctrl","kmer7","gene"]].iterrows():
            cnt[r.site]=cnt.get(r.site,0)+1; esb.setdefault(r.site,[]).append(r.ESB_ctrl); km[r.site]=r.kmer7; gn[r.site]=r.gene
    sites=[s for s,n in cnt.items() if n>=MIN_W]
    percond[cond]=pd.DataFrame({"site":sites,"ESB":[np.mean(esb[s]) for s in sites],
                                "kmer7":[km[s] for s in sites],"gene":[gn[s] for s in sites]})
    A_tested[cond]=max(a_tested,1)
    print(f"  {cond}: {len(sites)} m6A | A tested(mta)={a_tested}",flush=True)

allg=set().union(*[set(percond[c].site) for c in CONDS])
master=pd.DataFrame({"site":sorted(allg)})
for c in CONDS: master[c]=master.site.isin(set(percond[c].site)).astype(int)
gmap={}; [gmap.update(dict(zip(percond[c].site,percond[c].gene))) for c in CONDS]
master["gene"]=master.site.map(gmap); master.to_csv(f"{OUT}/m6A_master_table.csv",index=False)

# FIG1 motif (re-centred, per condition 2x2)
try:
    import logomaker
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for ax,c in zip(axes.ravel(),CONDS):
        w5=recentre(list(percond[c]["kmer7"]))
        if not w5: ax.set_visible(False); continue
        m=logomaker.alignment_to_matrix(w5); info=logomaker.transform_matrix(m,from_type="counts",to_type="information")
        logomaker.Logo(info,ax=ax,color_scheme="classic")
        ax.set_title(f"{LABEL[c]}  (n={len(w5):,})"); ax.set_ylabel("bits")
        ax.set_xticks(range(5)); ax.set_xticklabels(["-2","-1","m6A","+1","+2"])
    fig.suptitle("m6A sequence motif (RRACH) by condition",fontsize=23,fontweight="bold")
    save(fig,"01_motif_logo_byCond")
except Exception as e: print("  [logo]",e)

# FIG2 coverage-normalised m6A rate (+ raw counts annotated)
fig,axes=plt.subplots(1,2,figsize=(15,6))
raw_n=[len(percond[c]) for c in CONDS]
rate=[1000*len(percond[c])/A_tested[c] for c in CONDS]
for ax,(vals,ylab,ttl) in zip(axes,[(raw_n,"m6A sites (count)","Raw site count"),
                                     (rate,"m6A sites per 1,000 A tested","Coverage-normalised rate")]):
    b=ax.bar([LABEL[c] for c in CONDS],vals,color=[PAL[c] for c in CONDS],edgecolor="black",lw=0.8)
    for bar,v in zip(b,vals): ax.text(bar.get_x()+bar.get_width()/2,v,(f"{v:,.0f}" if v>100 else f"{v:.1f}"),ha="center",va="bottom",fontsize=14,fontweight="bold")
    ax.set_ylabel(ylab); ax.set_title(ttl); ax.set_xticklabels([LABEL[c] for c in CONDS],rotation=18,ha="right")
fig.suptitle("m6A sites per condition — raw vs coverage-normalised",fontsize=22,fontweight="bold")
save(fig,"02_sites_raw_vs_normalised")

# FIG3 UpSet
try:
    from upsetplot import UpSet, from_contents
    data=from_contents({LABEL[c]:set(percond[c].site) for c in CONDS})
    fig=plt.figure(figsize=(12,7))
    UpSet(data,subset_size="count",show_counts=True,sort_by="cardinality").plot(fig=fig)
    fig.suptitle("Shared and condition-specific m6A sites",fontsize=22,fontweight="bold")
    save(fig,"03_upset_overlap")
except Exception as e: print("  [upset]",e)

# FIG4 volcano panels (vs mta)
fig,axes=plt.subplots(2,2,figsize=(14,12),sharex=True,sharey=True)
for ax,c in zip(axes.ravel(),CONDS):
    d=raw[c].copy(); d["nlp"]=-np.log10(d.pval.clip(lower=1e-300))
    hit=(d.pval<0.05)&(d.oddR<1)&(d.ESB_ctrl>d.ESB_test)
    ax.scatter(np.log2(d.loc[~hit,"oddR"].clip(lower=1e-3)),d.loc[~hit,"nlp"],s=4,c="#cccccc",alpha=0.3,rasterized=True)
    ax.scatter(np.log2(d.loc[hit,"oddR"].clip(lower=1e-3)),d.loc[hit,"nlp"],s=7,c=PAL[c],alpha=0.55,rasterized=True)
    ax.axhline(-np.log10(0.05),ls="--",c="k",lw=1); ax.axvline(0,ls="--",c="k",lw=1)
    ax.set_title(f"{LABEL[c]}  ·  {int(hit.sum()):,} m6A")
fig.supxlabel("log2(odds ratio)   ←  more modified in sample",fontsize=19,fontweight="bold")
fig.supylabel("-log10(p-value)",fontsize=19,fontweight="bold")
fig.suptitle("ELIGOS2 error-based modification (vs mta writer mutant)",fontsize=22,fontweight="bold")
save(fig,"04_volcano_panels")

# FIG5 gene heatmap: top genes by total m6A sites, presence across conditions
g=master[master.gene.notna()].groupby("gene")[CONDS].sum()
g["tot"]=g.sum(axis=1); top=g.sort_values("tot",ascending=False).head(30).drop(columns="tot")
fig,ax=plt.subplots(figsize=(9,12))
sns.heatmap(top,cmap="rocket_r",annot=True,fmt="d",linewidths=.5,linecolor="white",
            cbar_kws={"label":"m6A sites in gene"},ax=ax,xticklabels=[LABEL[c] for c in CONDS])
ax.set_title("Top 30 genes by m6A site number",fontweight="bold"); ax.set_ylabel("gene (TAIR ID)"); ax.set_xlabel("")
plt.setp(ax.get_xticklabels(),rotation=18,ha="right")
save(fig,"05_gene_heatmap")

# FIG6 AA effect / ANAC017 dependence (normalised deltas)
def s(c): return set(percond[c].site)
panel={"AA-gained (WT)":len(s("WT_AA")-s("WT_C")),"AA-gained (anac017)":len(s("anac017-1_AA")-s("anac017-1_C")),
       "AA-lost (WT)":len(s("WT_C")-s("WT_AA")),"AA-lost (anac017)":len(s("anac017-1_C")-s("anac017-1_AA")),
       "WT_AA-specific vs anac017_AA":len(s("WT_AA")-s("anac017-1_AA"))}
fig,ax=plt.subplots(figsize=(11,6))
ks=list(panel); vs=[panel[k] for k in ks]
ax.barh(ks,vs,color=["#F58518","#E45756","#9ecae1","#fcbba1","#756bb1"],edgecolor="black",lw=0.8)
for i,v in enumerate(vs): ax.text(v,i,f" {v:,}",va="center",fontweight="bold",fontsize=14)
ax.invert_yaxis(); ax.set_xlabel("m6A sites"); ax.set_title("m6A remodelling: Antimycin A effect & ANAC017 dependence")
save(fig,"06_AA_ANAC017")
print("== DONE v2 ==",flush=True)
