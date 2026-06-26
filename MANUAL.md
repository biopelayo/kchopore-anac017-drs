# K-CHOPORE · anac017 DRS — Consolidated Code Manual

Direct RNA sequencing (DRS, Oxford Nanopore) analysis of the *Arabidopsis thaliana*
**WT / anac017-1 × Control / Antimycin A** 2×2 design (3 replicates per group, 12 libraries).
The workflow covers two parallel axes: the **transcriptome** axis (FLAIR isoforms — the
epitranscriptomic, primary track) and the **genome** axis (TAIR10 — alternative QC/DE track).

> This manual consolidates `METHODS_TOOLS.md` and `RUN_CLI.md`. All tool settings, filters
> and file paths are taken verbatim from those sources.

---

## 0. Scope and repo map

**Design.** 2×2 factorial: genotype (WT / anac017-1) × treatment (Control / Antimycin A),
3 replicates each = 12 DRS libraries. ANAC017 is the retrograde-signalling transcription
factor under study; Antimycin A triggers the mitochondrial retrograde response.

**Two analysis axes.**

- **Transcriptome axis** (primary, epitranscriptomic): reads mapped to FLAIR/StringTie
  isoforms; m6A detection and transcript-level differential expression.
- **Genome axis** (alternative QC/DE): reads mapped to TAIR10 with spliced alignment.

**Servers and locations.**

| Item | Location |
|------|----------|
| Compute server | `ssh <server>` (UniOvi VPN required) |
| Pipeline repo | `<conda-prefix>/kchopore-anac017-drs` |
| Working directory | `<data-dir>/run_transcriptome` |
| Conda env (pipeline) | `kchopore` (via `<conda-prefix>/etc/profile.d/conda.sh`) |
| Conda env (figures) | `viz` |
| System R | `Rscript` (has DESeq2 / ggplot2 / patchwork / ragg) |
| NAS data (read-only) | `<data-dir>/nas/Comun`, `<data-dir>/nas/HTData` (rclone mounts) |
| Final figures | `docs/figures/` |
| Result tables | `docs/tables/` |

**The single editable file is `config/config_transcriptome.yml`.** The `Snakefile` rules
are never touched for a new run.

---

## 1. Install / Docker

Two equivalent ways to run the QC + alignment pipeline: native conda, or the reproducible
Docker image (FAIR path).

### Conda (native)

```bash
ssh <server>
source <conda-prefix>/etc/profile.d/conda.sh
conda activate kchopore
```

### Docker (reproducible, FAIR)

```bash
cd <conda-prefix>/kchopore-anac017-drs
docker build -t kchopore-anac017-drs:latest .          # image already built
docker run --rm \
  -v $PWD:/workspace -v <data-dir>/nas:/nas -w /workspace \
  kchopore-anac017-drs:latest \
  snakemake --configfile config/config_transcriptome.yml \
            --cores 12 --rerun-triggers mtime --keep-going
```

The container mounts the repo at `/workspace` and the NAS at `/nas`, so the same Snakemake
invocation runs inside or outside Docker.

---

## 2. Configuration

### Mount the NAS (read-only, no sudo)

Data live on the NAS and are read **in place** — no copying.

```bash
rclone mount nas:Comun          <data-dir>/nas/Comun  --read-only --vfs-cache-mode off --dir-cache-time 720h --daemon
rclone mount nas:HTData_and_DBs <data-dir>/nas/HTData --read-only --vfs-cache-mode off --dir-cache-time 720h --daemon
ls <data-dir>/nas/Comun/Chus   # sanity check
```

### Prepare the workdir

```bash
conda activate kchopore
python3 scripts/setup_transcriptome.py   # builds run_transcriptome with fastq/summaries/BAMs/ref symlinked
```

### The config file

`config/config_transcriptome.yml` is the only file you edit for a new run. It controls:

| Key | Meaning |
|-----|---------|
| `samples:` | list of sample IDs |
| `conditions:` | control / treatment label per sample |
| `input_files:` | paths to fastq, summaries, reference |
| `params: run_*:` | which modules to run (`true` / `false`) |

The `Snakefile` rules are not modified.

---

## 3. Run the pipeline

### Golden rule

```
-n  (see the plan)   →   if it looks right, drop -n and launch
```

Always keep `--rerun-triggers mtime` so already-computed heavy steps are reused, not recomputed.

### 3.1 QC + alignment (Snakemake)

Run everything from the workdir:

```bash
cd <data-dir>/run_transcriptome
```

Dry-run first (always):

```bash
snakemake -n \
  --snakefile Snakefile \
  --configfile config/config_transcriptome.yml \
  --rerun-triggers mtime
```

`Nothing to be done` means everything is up to date (reuse, no recompute).

Execute for real:

```bash
snakemake \
  --snakefile Snakefile \
  --configfile config/config_transcriptome.yml \
  --cores 12 \
  --rerun-triggers mtime \
  --keep-going
# -> results/{nanoplot,nanocomp,quality_analysis,samtools_stats,multiqc}
```

Force one rule (to watch it run / demo it):

```bash
snakemake --snakefile Snakefile --configfile config/config_transcriptome.yml \
  --cores 4 --rerun-triggers mtime \
  --forcerun multiqc results/multiqc/multiqc_report.html
```

Request one specific target (Snakemake does only what's needed to produce it):

```bash
snakemake --snakefile Snakefile --configfile config/config_transcriptome.yml \
  --cores 4 --rerun-triggers mtime \
  results/nanoplot/WT_C_R1/NanoStats.txt
```

#### Snakemake flag reference

| Flag | What it does |
|------|--------------|
| `-n` | dry-run: shows what it would do, changes nothing |
| `--configfile` | the file that controls everything (samples, modules, paths) |
| `--rerun-triggers mtime` | decide what to redo by file date only → no recompute of finished heavy steps |
| `--cores N` / `-j N` | cores in parallel (40 available; 12–16 is plenty) |
| `--keep-going` | if a rule fails, continue with the rest |
| `--forcerun <rule>` | re-run a rule even if already done |
| `-p` / `--printshellcmds` | print the real shell command of each step |
| `-r` / `--reason` | explain why each rule runs |
| `--unlock` | unlock the directory after an interruption (then relaunch) |
| `--dag \| dot -Tpng > dag.png` | draw the dependency graph |

### 3.2 Differential expression (DESeq2, 2×2)

```bash
Rscript scripts/run_deseq2.R
# idxstats of the 12 BAMs -> counts -> DESeq2 -> out (PCA, volcanoes, heatmap)
```

### 3.3 GO / KEGG enrichment (g:Profiler)

```bash
conda activate kchopore
Rscript scripts/go_enrich.R      # -> go
```

### 3.4 m6A figures (ELIGOS2 already run upstream)

```bash
conda run -n viz python scripts/eligos_m6a_figures_v2.py
# -> eligos_figs_v2 (motif, normalised rate, volcano, heatmap)
# baseExt0 files already extracted in eligos_results
```

### 3.5 m6anet (orthogonal validation, subset by disk space)

```bash
# copy one eventalign locally and process with a local FIFO (faster than SMB):
cp "<data-dir>/nas/Comun/Chus/Chus_DRS_Nanopolish_eventalign_m6Anet_NO BORRAR/anac017-1_C_R3_eventalign.txt.gz" /tmp/
bash scripts/m6anet_stream.sh anac017-1_C_R3   # (adjust to read from /tmp)
```

---

## 4. Per-tool settings

All settings are verbatim from `METHODS_TOOLS.md`. The transcriptome axis is the primary
(epitranscriptomic) track; the genome axis is the alternative QC/DE track.

| # | Tool | Step | Settings / filters | Output | Plots |
|---|------|------|--------------------|--------|-------|
| 1 | **Guppy 6.2.1** (Chus) | Basecalling FAST5→FASTQ | `rna_r9.4.1_70bps_hac.cfg`, `--reverse_sequence 1 --u_substitution 1`, `--min_qscore 5` | fastq | — |
| 2 | **NanoFilt** | Read filtering | config `q≥7`, `length≥200` (transcriptome run used Chus fastq **without** re-filtering) | filtered fastq | — |
| 3 | **NanoPlot** | Per-sample QC | defaults | `NanoStats.txt` | ~15 plots/sample (L vs Q, yield, histograms) |
| 4 | **NanoComp** | Comparative QC | defaults, 12 samples | HTML report | violins length/quality/N50/reads, overlay histograms |
| 5 | **minimap2** | Mapping | **transcriptome:** `-ax map-ont -L`; **genome:** `-ax splice -k14 --secondary=no --MD` | SAM | — |
| 6 | **samtools** | Filter / sort | `view -F 2324 -q 10` (primary, MAPQ≥10), `sort`, `index` | BAM + bai | — |
| 7 | **samtools stats** | Alignment stats | `flagstat` + `stats` | txt | → MultiQC |
| 8 | **pycoQC** | Signal + align QC | `--min_pass_qual 7` (summary + BAM) | interactive HTML | reads/channels/time QC |
| 9 | **FLAIR / StringTie** | Isoforms | FLAIR collapse (support≥3, stringent); StringTie `min_cov 2.5` | isoforms.bed/gtf | (reused from Chus) |
| 10 | **ELIGOS2** | m6A by base-error | `rna_mod`/`pair_diff_mod` `--pval 0.05 --oddR 5 --esb 0.2`; on FLAIR transcriptome; **vs 5 writer-mutants** (mta/mtb/fip37/vir/hakai) | baseExt0/1/2, combine, `.A.filtered`, BedGraph (IGV) | native: BedGraph + QC; **ours**: RRACH motif, volcano, UpSet, heatmap |
| 11 | **m6anet** | m6A by signal | `dataprep --n_processes 1` (FIFO), `inference --num_iterations 1000` | `data.site_proba.csv` | (own downstream) |
| 12 | **xpore** | Differential m6A | **pending** (needs eventalign for 2 conditions) | `diffmod.table` | — |
| 13 | **DESeq2** | Differential expression | counts = `samtools idxstats` (transcript-level); design `~genotype+treatment+genotype:treatment`; filter `rowSums≥10`; significant `padj<0.05 & |LFC|>1` | CSV per contrast | PCA, 3× volcano, heatmap, dispersion |
| 14 | **g:Profiler** (gprofiler2) | GO / KEGG | `organism=athaliana`, `significant=TRUE`, GO:BP/MF/CC + KEGG; on DE genes (AT IDs) | `GO.csv` | barplot of top terms |
| 15 | **MultiQC** | QC aggregation | defaults | HTML report | aggregate |

### m6A calling criteria (ELIGOS2, our downstream)

- **Candidate m6A site:** `ref==A` & `pval<0.05` & `oddR<1` & `ESB_ctrl>ESB_test` (the sample
  has more base-error than the writer-mutant → loses m6A when the writer is removed).
- **Consensus rule (≥3/5 writers):** a site is m6A in a condition if it holds in **≥3 of 5** writer mutants.
- **Cross-condition normalisation:** **rate = sites / 1,000 A tested** (corrects for depth;
  important because WT_C had lower coverage).

### Open review points (from `METHODS_TOOLS.md`)

1. **NanoFilt not actually applied** in the transcriptome axis (Chus fastq used without
   re-filtering) → decide whether to re-filter `q≥7 l≥200`.
2. **oddR:** the ELIGOS config asks `oddR 5`, but Chus's `.A.filtered` appears pre-filtered to
   `pval<0.05` without a strict oddR threshold; our consensus uses `oddR<1` (m6A direction).
   Consider tightening (minimum `|log2 oddR|`).
3. **DESeq2 at transcript level** (not gene): de-novo FLAIR transcripts lack an AT gene → those
   are lost for GO. Alternative: aggregate by gene via overlap with AtRTD.
4. **m6anet** only feasible for a subset (disk space). Define scope.
5. **xpore** not yet run.

---

## 5. Figures and tables

Final figures are in `docs/figures/` (PNG, most with a matching PDF). The result tables that
back them are in `docs/tables/`. The composite paper figures (`fig1`–`fig5`) assemble the
single-panel figures listed below.

### Figure provenance

Each figure maps to its file in `docs/figures/` and the table in `docs/tables/` it is built
from.

| Figure | File | Table | Content |
|--------|------|-------|---------|
| PCA | `docs/figures/pca.png` | `docs/tables/tidy_de_genotype.csv` | transcript-level PCA (PC1 / PC2) of the 12 libraries |
| Volcanoes | `docs/figures/volcanoes.png` | `docs/tables/tidy_de_AA_in_WT.csv`, `docs/tables/tidy_de_genotype.csv`, `docs/tables/tidy_de_interaction.csv` | AA-in-WT / genotype / interaction contrasts |
| Heatmap top 40 | `docs/figures/heatmap-top40.png` | `docs/tables/tidy_de_interaction.csv` | z-score vst heatmap of the top 40 DE transcripts |
| Differential expression (composite) | `docs/figures/fig2-differential-expression.png` | `docs/tables/tidy_de_AA_in_WT.csv`, `docs/tables/tidy_de_genotype.csv`, `docs/tables/tidy_de_interaction.csv` | PCA + volcanoes + heatmap assembled |
| m6A rate (normalised) | `docs/figures/m6a-rate-normalised.png` | `docs/tables/m6a_summary.csv` | m6A sites per 1,000 A tested, per condition |
| 86% ANAC017-dependent m6A | `docs/figures/m6a-86pct-anac017-dependent.png` | `docs/tables/m6a_summary.csv` | stacked bars: AA-gained m6A that needs ANAC017 |
| m6A and ANAC017 (composite) | `docs/figures/fig3-m6a-anac017.png` | `docs/tables/m6a_summary.csv`, `docs/tables/m6a_motif_kmers.csv` | rate + motif + ANAC017 dependence assembled |
| m6A identity / metagene | `docs/figures/fig4-m6a-identity-metagene.png` | `docs/tables/taskA_mta_identity.csv`, `docs/tables/taskB_mapping_coverage.csv` | m6A identity (mta) and mapping coverage |
| AOX1a track | `docs/figures/aox1a-track.png` | — | ESB track over AOX1a (Chr3:7,906,521-7,908,746), 4 conditions |
| GO compareCluster | `docs/figures/go-compareCluster.png` | `docs/tables/compareCluster_GO_BP_interaction.csv` | GO:BP terms across contrasts |
| Function / GO (composite) | `docs/figures/fig5-function-go.png` | `docs/tables/ora_GO_BP_interaction.csv`, `docs/tables/kegg_interaction.csv` | GO / KEGG enrichment of the interaction set |
| QC + design (composite) | `docs/figures/fig1-qc-design.png` | — | NanoComp QC summary and the 2x2 design |

The RRACH motif comes from `docs/tables/m6a_motif_kmers.csv`. The DE contrasts hold 151
(AA-in-WT) / 93 (genotype) / 19 (interaction) significant transcripts at `padj<0.05 &
|LFC|>1`. The GO / KEGG enrichment tables (`ora_*`, `gseGO_*`, `kegg_*`, `compareCluster_*`)
are split by contrast (`AA_in_WT`, `genotype`, `interaction`).

---

## 6. Reproducibility

### NAS access

- Read-only rclone mounts, no sudo (see §2): `nas:Comun` → `<data-dir>/nas/Comun`,
  `nas:HTData_and_DBs` → `<data-dir>/nas/HTData`. Data are read in place; nothing is copied.
- Source DRS data and the precomputed ELIGOS2 / eventalign products live under
  `<data-dir>/nas/Comun`.

### Provenance

- Every figure is traced to its file in `docs/figures/` and the table in `docs/tables/` it is
  built from (see §5).
- The pipeline is config-driven: a run is fully described by
  `config/config_transcriptome.yml` plus the immutable `Snakefile`.

### Versions and determinism

- Pinned tool versions: **Guppy 6.2.1** (basecalling), ELIGOS2, m6anet, DESeq2,
  gprofiler2; the full toolchain is pinned in the Docker image
  `kchopore-anac017-drs:latest`.
- **m6anet inference:** `--num_iterations 1000` (fixed iteration count).
- **m6anet dataprep:** `--n_processes 1` over a FIFO (deterministic, low-memory).
- **ELIGOS2:** fixed thresholds `--pval 0.05 --oddR 5 --esb 0.2`; consensus ≥3/5 writer mutants.
- **DESeq2:** fixed design `~genotype+treatment+genotype:treatment`, filter `rowSums≥10`,
  significance `padj<0.05 & |LFC|>1`.

> No explicit random seed is recorded in the source documents. The deterministic knobs
> above (fixed thresholds, fixed iteration count, single-process dataprep) are what the
> sources specify; set/record an RNG seed if any stochastic step is added.

---

## Appendix A — Output glossary

| Output | Produced by | Meaning |
|--------|-------------|---------|
| `NanoStats.txt` | NanoPlot | per-sample read summary (no per-read dump) |
| NanoComp HTML | NanoComp | comparative QC across 12 libraries |
| `*.bam` / `*.bai` | samtools | sorted/indexed primary alignments, MAPQ≥10 (`-F 2324`) |
| samtools stats txt | samtools stats | flagstat + stats → MultiQC |
| pycoQC HTML | pycoQC | interactive signal + alignment QC |
| isoforms.bed/gtf | FLAIR / StringTie | collapsed isoform models |
| `*_baseExt0/1/2.txt` | ELIGOS2 | per-base error tables (extension windows) |
| `.A.filtered` | ELIGOS2 | A-site error calls filtered to `pval<0.05` |
| `*.A.ESB_ctrl.bdg` | ELIGOS2 | BedGraph of ESB at A sites (IGV / tracks) |
| `m6A_master_table.csv` | `scripts/eligos_m6a_figures_v2.py` | consolidated m6A calls (figure source) |
| `m6a_volcano_percond.csv` | `scripts/eligos_m6a_figures_v2.py` | per-condition A-site volcano (log2 oddR, −log10 p, hit) |
| `data.site_proba.csv` | m6anet | per-site m6A probability (signal-based) |
| `diffmod.table` | xpore | differential m6A (pending) |
| DESeq2 CSV per contrast | `scripts/run_deseq2.R` | transcript-level DE (AA, genotype, interaction) |
| `*_GO.csv` | g:Profiler | GO:BP/MF/CC + KEGG enrichment of DE genes |
| MultiQC HTML | MultiQC | aggregate QC report |
| result CSVs in `docs/tables/` | DESeq2 / GO / m6A steps | per-figure plotting and result tables |

---

## Appendix B — Pending cleanup

- **xpore** not yet run (needs eventalign for both conditions).
- **m6anet** scope limited by disk space — only a subset of libraries processed.
- **NanoFilt** not actually applied on the transcriptome axis; decide whether to re-filter
  `q≥7 l≥200`.
- **ELIGOS2 oddR threshold** inconsistency to resolve (config `oddR 5` vs Chus `.A.filtered`
  pre-filtered on `pval` only vs our `oddR<1` consensus).
- **DESeq2 transcript-level** loses de-novo FLAIR transcripts for GO; consider gene-level
  aggregation via AtRTD overlap.
- Figures to **improve**: `m6a/fig2_rate_normalised`, `m6a/fig6_aa_anac017` (barplots),
  `m6a/fig3_upset_overlap`, `deseq2/F2_volcano_panel`, `deseq2/F3_MA_panel`.
- Figures to **discard**: `go/interaction_ANAC017dep_barplot` (use the dotplot instead).
- `integration/fig_integration_expr_m6a` is descriptive only — do **not** present it as
  expression–m6A convergence.
- No RNG seed recorded; add one if any stochastic step is introduced.
