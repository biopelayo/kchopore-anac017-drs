# K-CHOPORE · Runbook CLI (corrida manual end-to-end)

Servidor: `ssh usuario2@156.35.42.17` (VPN UniOvi activa). Repo: `~/kchopore-anac017-drs`.

## 0. Montar el NAS (datos, read-only, sin sudo)
```bash
~/bin/rclone mount nas:Comun ~/nas/Comun --read-only --vfs-cache-mode off --dir-cache-time 720h --daemon
~/bin/rclone mount nas:HTData_and_DBs ~/nas/HTData --read-only --vfs-cache-mode off --dir-cache-time 720h --daemon
ls ~/nas/Comun/Chus   # comprobar
```

## 1. Preparar el workdir (datos del NAS in place)
```bash
source ~/miniconda/etc/profile.d/conda.sh && conda activate kchopore
python3 ~/setup_transcriptome.py          # crea run_transcriptome con fastq/summaries/BAMs/ref enlazados
```

## 2. Pipeline QC + alineamiento (Snakemake — solo rehace lo que falta)
```bash
cd /media/usuario2/ssd4TB1/kchopore_arabidopsis/run_transcriptome
snakemake -n --configfile config/config_transcriptome.yml --rerun-triggers mtime   # DRY-RUN (ver qué hará)
snakemake    --configfile config/config_transcriptome.yml --cores 12 --rerun-triggers mtime --keep-going
# -> results/{nanoplot,nanocomp,quality_analysis,samtools_stats,multiqc}
```

## 3. Expresión diferencial (DESeq2 2×2)
```bash
bash ~/run_deseq2.sh      # idxstats de los 12 BAMs -> counts -> DESeq2 -> ~/deseq2/out (PCA, volcanos, heatmap)
```

## 4. GO/KEGG (g:Profiler)
```bash
conda activate kchopore && Rscript ~/go_enrich.R     # -> ~/deseq2/go
```

## 5. m6A (ELIGOS2 ya hecho por Chus) — figuras
```bash
conda run -n viz python ~/eligos_m6a_figures_v2.py   # -> ~/eligos_figs_v2 (motivo, normalizado, volcano, heatmap)
# (los baseExt0 ya extraidos en ~/eligos_results)
```

## 6. m6anet (validación ortogonal, subconjunto por espacio)
```bash
# bajar 1 eventalign a local y procesar con FIFO local (mas rapido que SMB):
cp "~/nas/Comun/Chus/Chus_DRS_Nanopolish_eventalign_m6Anet_NO BORRAR/anac017-1_C_R3_eventalign.txt.gz" /tmp/
bash ~/m6anet_stream.sh anac017-1_C_R3   # (ajustar para leer de /tmp)
```

## Vía Docker (reproducible, FAIR)
```bash
cd ~/kchopore-anac017-drs
docker build -t kchopore-anac017-drs:latest .          # ya construida
docker run --rm -v $PWD:/workspace -v ~/nas:/nas -w /workspace kchopore-anac017-drs:latest \
  snakemake --configfile config/config_transcriptome.yml --cores 12 --rerun-triggers mtime --keep-going
```

## Qué editar para una corrida nueva
Solo **`config/config_transcriptome.yml`**: `samples`, `conditions`, rutas en `input_files`, y los flags `run_*` (qué módulos correr). Las reglas no se tocan.
