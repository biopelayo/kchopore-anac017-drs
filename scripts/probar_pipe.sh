#!/usr/bin/env bash
# =============================================================
# K-CHOPORE · anac017 DRS — probar el pipeline end-to-end
# Idempotente: reutiliza todo lo ya generado (Snakemake --rerun-triggers mtime),
# solo rehace lo que falte. Ejecutar en EpiPower (usuario2).
#
#   bash ~/probar_pipe.sh [etapa]
#   etapa ∈ {check, dryrun, qc, deseq2, go, m6a, docker, all}   (def: check)
# =============================================================
set -uo pipefail

RUN="/media/usuario2/ssd4TB1/kchopore_arabidopsis/run_transcriptome"
CFG="config/config_transcriptome.yml"
CONDA="$HOME/miniconda/etc/profile.d/conda.sh"
IMG="kchopore-anac017-drs:latest"
ETAPA="${1:-check}"

c(){ printf "\033[1;32m==> %s\033[0m\n" "$*"; }   # verde
w(){ printf "\033[1;33m[!] %s\033[0m\n" "$*"; }   # amarillo
e(){ printf "\033[1;31m[x] %s\033[0m\n" "$*"; }   # rojo
source "$CONDA"

# ---------- check: entorno, NAS, workdir, imagen ----------
check(){
  c "Entorno"
  command -v snakemake >/dev/null || conda activate kchopore
  conda activate kchopore 2>/dev/null
  echo "  snakemake $(snakemake --version 2>/dev/null) | $(samtools --version|head -1)"
  c "NAS montado"
  if mountpoint -q "$HOME/nas/Comun" 2>/dev/null || ls "$HOME/nas/Comun/Chus" >/dev/null 2>&1; then
    echo "  NAS OK"
  else
    w "NAS no montado; montando…"
    "$HOME/bin/rclone" mount nas:Comun "$HOME/nas/Comun" --read-only --vfs-cache-mode off --dir-cache-time 720h --daemon
    sleep 5; ls "$HOME/nas/Comun/Chus" >/dev/null 2>&1 && echo "  NAS OK" || e "fallo NAS"
  fi
  c "Workdir + datos in place"
  [ -d "$RUN" ] && echo "  $RUN OK" || { e "falta $RUN — corre: python3 ~/setup_transcriptome.py"; return 1; }
  echo "  BAMs: $(ls "$RUN"/results/sorted_bam/*_sorted.bam 2>/dev/null|wc -l)/12 | fastq: $(ls "$RUN"/data/raw/fastq/*.fastq 2>/dev/null|wc -l)/12 | summaries: $(ls "$RUN"/data/raw/summaries/*.txt 2>/dev/null|wc -l)/12"
  c "Imagen Docker"
  docker image inspect "$IMG" >/dev/null 2>&1 && echo "  $IMG OK" || w "imagen no encontrada (docker build en ~/kchopore-anac017-drs)"
}

# ---------- dryrun: qué haría Snakemake (debe ser solo QC/plots) ----------
dryrun(){ c "DRY-RUN (plan; reutiliza lo hecho)"; cd "$RUN"; conda activate kchopore
  snakemake -n --snakefile Snakefile --configfile "$CFG" --rerun-triggers mtime | tail -25; }

# ---------- qc: pipeline Snakemake (QC/align/plots) ----------
qc(){ c "Pipeline QC (Snakemake)"; cd "$RUN"; conda activate kchopore
  snakemake --snakefile Snakefile --configfile "$CFG" --cores 12 --rerun-triggers mtime --keep-going --printshellcmds
  c "Resultados QC"
  echo "  NanoPlot: $(ls -d results/nanoplot/*/ 2>/dev/null|wc -l)/12 | pycoQC: $(ls results/quality_analysis/*.html 2>/dev/null|wc -l)/12 | NanoComp: $([ -f results/nanocomp/NanoComp-report.html ]&&echo si||echo no) | MultiQC: $([ -f results/multiqc/multiqc_report.html ]&&echo si||echo no)"; }

# ---------- forzar una regla para VERLA correr (demo) ----------
demo(){ c "Demo: forzar MultiQC para verlo correr"; cd "$RUN"; conda activate kchopore
  snakemake --snakefile Snakefile --configfile "$CFG" --cores 4 --rerun-triggers mtime --forcerun multiqc results/multiqc/multiqc_report.html; }

deseq2(){ c "DESeq2 2x2"; bash "$HOME/run_deseq2.sh"; ls "$HOME/deseq2/out"/*.png 2>/dev/null|wc -l|xargs echo "  figuras:"; }
go(){    c "GO/KEGG (g:Profiler)"; conda activate kchopore; Rscript "$HOME/go_enrich.R"; }
m6a(){   c "m6A figuras (ELIGOS2 downstream)"; conda activate viz; python "$HOME/eligos_m6a_paper.py"; }

docker_run(){ c "Pipeline vía Docker (reproducible FAIR)"; cd "$RUN"
  docker run --rm -v "$RUN":/workspace -v "$HOME/nas":/nas -w /workspace "$IMG" \
    snakemake --snakefile Snakefile --configfile "$CFG" --cores 12 --rerun-triggers mtime --keep-going; }

case "$ETAPA" in
  check)  check ;;
  dryrun) check && dryrun ;;
  qc)     check && qc ;;
  demo)   check && demo ;;
  deseq2) deseq2 ;;
  go)     go ;;
  m6a)    m6a ;;
  docker) check && docker_run ;;
  all)    check && dryrun && qc && deseq2 && go && m6a ;;
  *) e "etapa desconocida: $ETAPA"; echo "usa: check|dryrun|qc|demo|deseq2|go|m6a|docker|all" ;;
esac
