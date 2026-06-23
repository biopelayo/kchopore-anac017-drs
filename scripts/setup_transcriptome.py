#!/usr/bin/env python3
# Monta el workdir run_transcriptome para K-CHOPORE leyendo del NAS + BAMs locales,
# con todo lo pesado IN PLACE para que Snakemake solo rehaga QC/plots.
import os, glob, subprocess, textwrap
HOME=os.path.expanduser("~")
V=f"{HOME}/nas".replace("/nas","")  # placeholder
V_KCHOP="/media/usuario2/ssd4TB1/pelamovic/K-CHOPORE-V2"
RUN="/media/usuario2/ssd4TB1/kchopore_arabidopsis/run_transcriptome"
NAS=f"{HOME}/nas/Comun/Chus"
GUPPY=f"{NAS}/Chus_DRS_guppy/guppy"
SUMS=f"{NAS}/guppy_sequencing_summaries"
REF_SRC=f"{NAS}/Chus_DRS_guppy/flair.collapse.isoforms.fa"
BAMS=f"{HOME}/bams_transcriptome"

SAMPLES=["WT_C_R1","WT_C_R2","WT_C_R3","WT_AA_R1","WT_AA_R2","WT_AA_R3",
         "anac017-1_C_R1","anac017-1_C_R2","anac017-1_C_R3",
         "anac017-1_AA_R1","anac017-1_AA_R2","anac017-1_AA_R3"]
# mapeo de summary (typo en C_R3 del NAS)
def sum_name(s): return "anac017-_C_R3_sequencing_summary.txt" if s=="anac017-1_C_R3" else f"{s}_sequencing_summary.txt"

def ln(src,dst):
    if os.path.lexists(dst): os.remove(dst)
    os.symlink(src,dst)

os.makedirs(RUN,exist_ok=True)
for d in ["config","data/raw/fastq","data/raw/summaries","data/reference/genome",
          "data/reference/annotations","results/fastq_filtered","results/sorted_bam","logs"]:
    os.makedirs(f"{RUN}/{d}",exist_ok=True)
# symlinks del pipeline
for x in ["Snakefile","scripts","envs"]:
    ln(f"{V_KCHOP}/{x}",f"{RUN}/{x}")
# referencia (transcriptoma) copiada y escribible
ref_dst=f"{RUN}/data/reference/genome/transcriptome.fa"
if not os.path.exists(ref_dst):
    subprocess.run(["cp",REF_SRC,ref_dst],check=True)
subprocess.run(f"source {HOME}/miniconda/etc/profile.d/conda.sh; conda activate kchopore; "
               f"samtools faidx {ref_dst}; minimap2 -d {ref_dst}.mmi -k14 {ref_dst}",shell=True,executable="/bin/bash")

missing=[]
for s in SAMPLES:
    # fastq crudo + filtered (symlink al fastq -> nanofilt no recomputa)
    fq=f"{GUPPY}/{s}.fastq"
    ln(fq,f"{RUN}/data/raw/fastq/{s}.fastq")
    ln(fq,f"{RUN}/results/fastq_filtered/{s}_filtered.fastq")
    # summary
    su=f"{SUMS}/{sum_name(s)}"
    if os.path.exists(su): ln(su,f"{RUN}/data/raw/summaries/{s}_sequencing_summary.txt")
    else: missing.append(f"summary {s}")
    # bam local (mas nuevo que fastq -> minimap2/sort no recomputan)
    bam=f"{BAMS}/{s}.sorted.bam"
    if os.path.exists(bam):
        ln(bam,f"{RUN}/results/sorted_bam/{s}_sorted.bam")
        if os.path.exists(bam+".bai"): ln(bam+".bai",f"{RUN}/results/sorted_bam/{s}_sorted.bam.bai")
    else: missing.append(f"bam {s}")

# config
seqs="\n".join(f'    - "data/raw/summaries/{s}_sequencing_summary.txt"' for s in SAMPLES)
fqs="\n".join(f'    - "data/raw/fastq/{s}.fastq"' for s in SAMPLES)
slist="\n".join(f'  - "{s}"' for s in SAMPLES)
conds="\n".join(f'  {s}: "{"control" if "_C_" in s else "treatment"}"' for s in SAMPLES)
cfg=textwrap.dedent(f"""\
samples:
{slist}
conditions:
{conds}
samples_with_fast5: []
input_files:
  fastq_dir: "data/raw/fastq"
  fast5_dir: "data/raw/fast5"
  pod5_dir: "data/raw/pod5"
  sequencing_summaries_dir: "data/raw/summaries"
  fastq_files:
{fqs}
  sequencing_summaries:
{seqs}
  reference_genome: "data/reference/genome/transcriptome.fa"
  reference_genome_mmi: "data/reference/genome/transcriptome.fa.mmi"
  gtf_file: "data/reference/annotations/AtRTDv2_QUASI_19April2016.gtf"
  bed_file: "results/flair/flair.collapse.isoforms.bed"
  transcriptome_fasta: "data/reference/genome/transcriptome.fa"
  reads_manifest: "data/raw/reads_manifest.tsv"
output:
  path: "results"
  logs: "logs"
""")
# anexar bloque output completo y tools/params del config original (reusar el resto)
import re
orig=open(f"{V_KCHOP}/config/config.yml").read()
# tomar de 'output:' en adelante del original (rutas de salida + tools + params)
out_block=orig[orig.index("output:"):]
# tools/params: forzar map-ont y QC on, pesados off
import yaml
base=yaml.safe_load(orig)
newc=yaml.safe_load(cfg)
base.update({k:newc[k] for k in ["samples","conditions","samples_with_fast5","input_files"]})
base["tools"]["minimap2_preset"]="map-ont"
p=base["params"]
for k,v in {"run_basecalling":False,"run_nanofilt":True,"run_nanoplot":True,"run_nanocomp":True,
            "run_pycoqc":True,"run_flair":False,"run_stringtie":False,"run_bambu":False,
            "run_eligos2":False,"run_m6anet":False,"run_xpore":False,"run_deseq2":False,"run_multiqc":True}.items():
    p[k]=v
yaml.safe_dump(base,open(f"{RUN}/config/config_transcriptome.yml","w"),sort_keys=False)
# manifest
with open(f"{RUN}/data/raw/reads_manifest.tsv","w") as fh:
    for s in SAMPLES:
        fh.write(f"{s}\t{'control' if '_C_' in s else 'treatment'}\tbatch1\tdata/raw/fastq/{s}.fastq\n")
print("WORKDIR:",RUN)
print("samples:",len(SAMPLES),"| BAMs in place:",len(glob.glob(f'{RUN}/results/sorted_bam/*_sorted.bam')),
      "| summaries:",len(glob.glob(f'{RUN}/data/raw/summaries/*.txt')))
if missing: print("FALTAN:",missing)
else: print("todo in place")
