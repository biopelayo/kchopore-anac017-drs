#!/usr/bin/env python3
# =============================================================
# K-CHOPORE | m6A epitranscriptome desde ELIGOS2 (pair_diff vs writers)
# Consolida los 20 contrastes (4 condiciones x 5 writers del complejo m6A),
# llama sitios m6A por consenso de writers, y genera figuras publication-grade.
# Diseno 2x2: WT/anac017-1 x C(control)/AA(Antimycin A).
# =============================================================
import os, re, glob, sys
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

RES = os.path.expanduser("~/eligos_results")
OUT = os.path.expanduser("~/eligos_figs"); os.makedirs(OUT, exist_ok=True)
CONDS = ["WT_C", "WT_AA", "anac017-1_C", "anac017-1_AA"]
WRITERS = ["mta", "mtb", "fip37", "hakai", "vir"]
MIN_WRITERS = 3            # consenso: sitio m6A si significativo vs >=3 writers
USECOLS = ["chrom","start_loc","strand","name","ref","kmer7","oddR","pval","ESB_test","ESB_ctrl"]
PAL = {"WT_C":"#4C78A8","WT_AA":"#F58518","anac017-1_C":"#54A24B","anac017-1_AA":"#E45756"}
ATRE = re.compile(r"AT[1-5CM]G\d{5}", re.I)

# ---- estilo god-level ----
sns.set_theme(context="talk", style="ticks")
mpl.rcParams.update({
    "figure.dpi":120,"savefig.dpi":300,"savefig.bbox":"tight",
    "font.family":"DejaVu Sans","axes.titleweight":"bold","axes.spines.top":False,
    "axes.spines.right":False,"axes.grid":True,"grid.alpha":0.25,"axes.axisbelow":True})
def save(fig,name):
    for ext in ("png","pdf"): fig.savefig(f"{OUT}/{name}.{ext}")
    plt.close(fig); print("  fig:",name)

def find_file(cond,writer):
    g=glob.glob(f"{RES}/{cond}_results_{writer}/{writer}_sorted_merged_vs_{cond}_sorted_merged_*_baseExt0.txt")
    g=[x for x in g if "/._" not in x and "._" not in os.path.basename(x)]
    return g[0] if g else None

def load_sites(cond,writer):
    f=find_file(cond,writer)
    if not f: return None
    df=pd.read_csv(f,sep="\t",usecols=USECOLS,dtype={"chrom":str,"strand":str},low_memory=False)
    df=df[df["ref"].astype(str).str.upper()=="A"]
    for c in ("oddR","pval","ESB_test","ESB_ctrl"):
        df[c]=pd.to_numeric(df[c],errors="coerce")
    df=df.dropna(subset=["oddR","pval","ESB_test","ESB_ctrl"])
    # m6A: la condicion (ctrl) tiene mas error que el writer-mutante (test)
    sig=df[(df["pval"]<0.05)&(df["oddR"]<1)&(df["ESB_ctrl"]>df["ESB_test"])].copy()
    sig["site"]=sig["chrom"]+":"+sig["start_loc"].astype(str)+":"+sig["strand"]
    sig["gene"]=sig["name"].astype(str).str.extract("("+ATRE.pattern+")",expand=False)
    return sig

print("== Consolidando 20 contrastes ==",flush=True)
percond={}; raw_one={}
for cond in CONDS:
    counter={}; esb={}; kmer={}; gene={}
    nwr=0
    for w in WRITERS:
        s=load_sites(cond,w)
        if s is None: print(f"  falta {cond} vs {w}"); continue
        nwr+=1
        if cond not in raw_one: raw_one[cond]=(w,s)  # guarda 1 contraste para volcano
        for _,r in s[["site","ESB_ctrl","kmer7","gene"]].iterrows():
            counter[r.site]=counter.get(r.site,0)+1
            esb.setdefault(r.site,[]).append(r.ESB_ctrl)
            kmer[r.site]=r.kmer7; gene[r.site]=r.gene
    # consenso
    sites=[s for s,n in counter.items() if n>=MIN_WRITERS]
    percond[cond]=pd.DataFrame({"site":sites,
        "n_writers":[counter[s] for s in sites],
        "ESB":[np.mean(esb[s]) for s in sites],
        "kmer7":[kmer[s] for s in sites],
        "gene":[gene[s] for s in sites]})
    print(f"  {cond}: {len(sites)} sitios m6A consenso (>= {MIN_WRITERS}/{nwr} writers)",flush=True)

# tabla maestra
allsites=sorted(set().union(*[set(d.site) for d in percond.values()]))
master=pd.DataFrame({"site":allsites})
for c in CONDS:
    m=percond[c].set_index("site")
    master[f"m6A_{c}"]=master.site.isin(m.index).astype(int)
    master[f"ESB_{c}"]=master.site.map(m["ESB"]) if len(m) else np.nan
gmap={}
for c in CONDS:
    for _,r in percond[c].iterrows(): gmap[r.site]=r.gene
master["gene"]=master.site.map(gmap)
master.to_csv(f"{OUT}/m6A_master_table.csv",index=False)
print("  tabla maestra:",master.shape,"->",f"{OUT}/m6A_master_table.csv",flush=True)

# ---------- FIG 1: sequence logo (motivo) ----------
try:
    import logomaker
    kmers=[k for k in percond["WT_C"]["kmer7"] if isinstance(k,str) and len(k)==7]
    mat=logomaker.alignment_to_matrix([k.replace("U","T") for k in kmers])
    info=logomaker.transform_matrix(mat,from_type="counts",to_type="information")
    fig,ax=plt.subplots(figsize=(6,3))
    logomaker.Logo(info,ax=ax,color_scheme="classic")
    ax.set_title("Motivo de los sitios m6A (WT_C)  ·  centro = A modificada")
    ax.set_xlabel("posicion (k-mer de 7)"); ax.set_ylabel("bits")
    ax.set_xticks(range(7)); ax.set_xticklabels(["-3","-2","-1","A","+1","+2","+3"])
    save(fig,"01_motif_logo_WT_C")
except Exception as e: print("  [logo omitido]",e)

# ---------- FIG 2: nº sitios m6A por condicion ----------
fig,ax=plt.subplots(figsize=(7,5))
counts=[len(percond[c]) for c in CONDS]
bars=ax.bar(CONDS,counts,color=[PAL[c] for c in CONDS],edgecolor="black",linewidth=0.6)
for b,v in zip(bars,counts): ax.text(b.get_x()+b.get_width()/2,v,f"{v:,}",ha="center",va="bottom",fontsize=12,fontweight="bold")
ax.set_ylabel("nº sitios m6A (consenso ≥3 writers)"); ax.set_title("Sitios m6A por condicion")
ax.set_xticklabels(CONDS,rotation=20,ha="right")
save(fig,"02_nsites_per_condition")

# ---------- FIG 3: UpSet solapamiento ----------
try:
    from upsetplot import UpSet, from_contents
    data=from_contents({c:set(percond[c].site) for c in CONDS})
    fig=plt.figure(figsize=(10,6))
    UpSet(data,subset_size="count",show_counts=True,sort_by="cardinality").plot(fig=fig)
    fig.suptitle("Solapamiento de sitios m6A entre condiciones",fontweight="bold")
    save(fig,"03_upset_overlap")
except Exception as e: print("  [upset omitido]",e)

# ---------- FIG 4: volcano (un contraste por condicion: vs mta) ----------
fig,axes=plt.subplots(2,2,figsize=(13,11),sharex=True,sharey=True)
for ax,cond in zip(axes.ravel(),CONDS):
    w,s=raw_one.get(cond,(None,None))
    f=find_file(cond,"mta") or (find_file(cond,w) if w else None)
    if not f: ax.set_visible(False); continue
    d=pd.read_csv(f,sep="\t",usecols=["ref","oddR","pval","ESB_ctrl","ESB_test"],low_memory=False)
    d=d[d["ref"].astype(str).str.upper()=="A"]
    for c in ("oddR","pval","ESB_ctrl","ESB_test"): d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d.dropna(); d["nlp"]=-np.log10(d["pval"].clip(lower=1e-300))
    hit=(d["pval"]<0.05)&(d["oddR"]<1)&(d["ESB_ctrl"]>d["ESB_test"])
    ax.scatter(np.log2(d.loc[~hit,"oddR"].clip(lower=1e-3)),d.loc[~hit,"nlp"],s=4,c="#bbbbbb",alpha=0.3,rasterized=True)
    ax.scatter(np.log2(d.loc[hit,"oddR"].clip(lower=1e-3)),d.loc[hit,"nlp"],s=6,c=PAL[cond],alpha=0.6,rasterized=True)
    ax.axhline(-np.log10(0.05),ls="--",c="k",lw=0.8); ax.axvline(0,ls="--",c="k",lw=0.8)
    ax.set_title(f"{cond}  (vs mta)  ·  {int(hit.sum()):,} m6A");
fig.supxlabel("log2(oddR)  ← mas modificado en la condicion"); fig.supylabel("-log10(pval)")
fig.suptitle("Volcano de modificacion por error (ELIGOS2)",fontweight="bold")
save(fig,"04_volcano_panels")

# ---------- FIG 5: distribucion ESB de sitios m6A ----------
fig,ax=plt.subplots(figsize=(8,5))
dd=pd.concat([percond[c].assign(cond=c) for c in CONDS])
sns.violinplot(data=dd,x="cond",y="ESB",hue="cond",palette=PAL,legend=False,cut=0,inner="quartile",ax=ax)
ax.set_xticklabels(CONDS,rotation=20,ha="right"); ax.set_ylabel("ESB (error specific base)")
ax.set_title("Intensidad de la senal m6A por condicion"); save(fig,"05_esb_violin")

# ---------- FIG 6: efecto AA (ganancia/perdida) y dependencia ANAC017 ----------
def setof(c): return set(percond[c].site)
comp={
 "WT: AA gana": len(setof("WT_AA")-setof("WT_C")),
 "WT: AA pierde": len(setof("WT_C")-setof("WT_AA")),
 "anac017: AA gana": len(setof("anac017-1_AA")-setof("anac017-1_C")),
 "anac017: AA pierde": len(setof("anac017-1_C")-setof("anac017-1_AA")),
 "ANAC017-dep (WT_AA no anac017_AA)": len(setof("WT_AA")-setof("anac017-1_AA")),
}
fig,ax=plt.subplots(figsize=(9,5))
ks=list(comp); vs=[comp[k] for k in ks]
ax.barh(ks,vs,color=["#F58518","#9ecae1","#E45756","#c7e9c0","#756bb1"],edgecolor="black",linewidth=0.6)
for i,v in enumerate(vs): ax.text(v,i,f" {v:,}",va="center",fontweight="bold")
ax.set_xlabel("nº sitios m6A"); ax.set_title("Remodelado de m6A: efecto de AA y dependencia de ANAC017")
save(fig,"06_AA_effect_ANAC017")

print("== HECHO. Figuras en",OUT,"==")
print(open(f"{OUT}/m6A_master_table.csv").readline())
