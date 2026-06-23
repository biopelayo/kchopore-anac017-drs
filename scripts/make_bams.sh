#!/usr/bin/env bash
source /home/usuario2/miniconda/etc/profile.d/conda.sh; conda activate kchopore
REF="$HOME/nas/Comun/Chus/Chus_DRS_guppy/flair.collapse.isoforms.fa"
GUPPY="$HOME/nas/Comun/Chus/Chus_DRS_guppy/guppy"
OUT="$HOME/bams_transcriptome"
for s in anac017-1_C_R2 anac017-1_C_R3 anac017-1_AA_R1 anac017-1_AA_R2 anac017-1_AA_R3; do
  echo "[$(date)] $s"
  [ -s "$OUT/${s}.sorted.bam" ] && { echo "  ya existe"; continue; }
  minimap2 -ax map-ont -L -t 12 "$REF" "$GUPPY/${s}.fastq" 2>"$OUT/${s}.mm2.log" \
    | samtools view -bh -F 2324 -q 10 - \
    | samtools sort -@ 4 -O bam -o "$OUT/${s}.sorted.bam" -
  samtools index "$OUT/${s}.sorted.bam"
  echo "  -> $(du -h $OUT/${s}.sorted.bam | cut -f1)  $(samtools flagstat $OUT/${s}.sorted.bam | grep -m1 mapped)"
done
echo BAMS_DONE
