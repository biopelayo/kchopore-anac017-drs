# K-CHOPORE anac017 DRS runbook

End-to-end manual run on the compute server. Connect with `ssh <server>` (UniOvi VPN
required). The pipeline repo lives at `<conda-prefix>/kchopore-anac017-drs`. The working
directory for a run is `<data-dir>/run_transcriptome`.

The only file you edit for a new run is `config/config_transcriptome.yml`: `samples`,
`conditions`, the `input_files` paths, and the `run_*` flags that switch modules on or off.
The `Snakefile` rules are never touched.

## 0. Mount the NAS (read-only, no sudo)

Data live on the NAS and are read in place; nothing is copied.

```bash
rclone mount nas:Comun          <data-dir>/nas/Comun  --read-only --vfs-cache-mode off --dir-cache-time 720h --daemon
rclone mount nas:HTData_and_DBs <data-dir>/nas/HTData --read-only --vfs-cache-mode off --dir-cache-time 720h --daemon
ls <data-dir>/nas/Comun/Chus   # sanity check
```

## 1. Prepare the working directory

```bash
source <conda-prefix>/etc/profile.d/conda.sh
conda activate kchopore
python3 scripts/setup_transcriptome.py   # builds run_transcriptome with fastq/summaries/BAMs/ref symlinked
```

## 2. QC and alignment (Snakemake)

Run everything from the working directory:

```bash
cd <data-dir>/run_transcriptome
```

Dry-run first to see the plan. `Nothing to be done` means everything is up to date (reuse, no
recompute). Always keep `--rerun-triggers mtime` so finished heavy steps are reused.

```bash
snakemake -n \
  --snakefile Snakefile \
  --configfile config/config_transcriptome.yml \
  --rerun-triggers mtime
```

If the plan looks right, drop `-n` and launch:

```bash
snakemake \
  --snakefile Snakefile \
  --configfile config/config_transcriptome.yml \
  --cores 12 \
  --rerun-triggers mtime \
  --keep-going
# -> results/{nanoplot,nanocomp,quality_analysis,samtools_stats,multiqc}
```

Force one rule to watch it run:

```bash
snakemake --snakefile Snakefile --configfile config/config_transcriptome.yml \
  --cores 4 --rerun-triggers mtime \
  --forcerun multiqc results/multiqc/multiqc_report.html
```

Request one specific target (Snakemake does only what is needed to produce it):

```bash
snakemake --snakefile Snakefile --configfile config/config_transcriptome.yml \
  --cores 4 --rerun-triggers mtime \
  results/nanoplot/WT_C_R1/NanoStats.txt
```

### Snakemake flag reference

| Flag | What it does |
|------|--------------|
| `-n` | dry-run: shows what it would do, changes nothing |
| `--configfile` | the file that controls everything (samples, modules, paths) |
| `--rerun-triggers mtime` | decide what to redo by file date only, so finished heavy steps are not recomputed |
| `--cores N` / `-j N` | cores in parallel (40 available; 12 to 16 is plenty) |
| `--keep-going` | if a rule fails, continue with the rest |
| `--forcerun <rule>` | re-run a rule even if already done |
| `-p` / `--printshellcmds` | print the real shell command of each step |
| `-r` / `--reason` | explain why each rule runs |
| `--unlock` | unlock the directory after an interruption, then relaunch |
| `--dag \| dot -Tpng > dag.png` | draw the dependency graph |

## 3. Differential expression (DESeq2, 2x2)

```bash
Rscript scripts/run_deseq2.R
# idxstats of the 12 BAMs -> counts -> DESeq2 -> out (PCA, volcanoes, heatmap)
```

## 4. GO / KEGG enrichment (g:Profiler)

```bash
conda activate kchopore
Rscript scripts/go_enrich.R      # -> go
```

## 5. m6A figures (ELIGOS2 already run upstream)

```bash
conda run -n viz python scripts/eligos_m6a_figures_v2.py
# -> eligos_figs_v2 (motif, normalised rate, volcano, heatmap)
# baseExt0 files already extracted in eligos_results
```

## 6. m6anet (orthogonal validation, subset by disk space)

```bash
# copy one eventalign locally and process through a local FIFO (faster than SMB):
cp "<data-dir>/nas/Comun/Chus/Chus_DRS_Nanopolish_eventalign_m6Anet_NO BORRAR/anac017-1_C_R3_eventalign.txt.gz" /tmp/
bash scripts/m6anet_stream.sh anac017-1_C_R3   # (adjust to read from /tmp)
```

## Docker (reproducible)

The container mounts the repo at `/workspace` and the NAS at `/nas`, so the same Snakemake
invocation runs inside or outside Docker.

```bash
cd <conda-prefix>/kchopore-anac017-drs
docker build -t kchopore-anac017-drs:latest .          # image already built
docker run --rm \
  -v $PWD:/workspace -v <data-dir>/nas:/nas -w /workspace \
  kchopore-anac017-drs:latest \
  snakemake --configfile config/config_transcriptome.yml \
            --cores 12 --rerun-triggers mtime --keep-going
```
