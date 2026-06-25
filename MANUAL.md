# K-CHOPORE · anac017 DRS — Consolidated Code Manual

Direct RNA sequencing (DRS, Oxford Nanopore) analysis of the *Arabidopsis thaliana*
**WT / anac017-1 × Control / Antimycin A** 2×2 design (3 replicates per group, 12 libraries).
The workflow covers two parallel axes: the **transcriptome** axis (FLAIR isoforms — the
epitranscriptomic, primary track) and the **genome** axis (TAIR10 — alternative QC/DE track).

> This manual consolidates `METHODS_TOOLS.md`, `RUN_CLI.md`, `COMANDOS_SNAKEMAKE.md`,
> `FIGS_PAPER/INVENTARIO_PROVENANCE.md` and `POSTER/FIGURAS_4_BLUEPRINT.md`. All tool
> settings, filters and file paths are taken verbatim from those sources.

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
| Compute server (EpiPower) | `ssh usuario2@156.35.42.17` (UniOvi VPN required) |
| Pipeline repo | `~/kchopore-anac017-drs` |
| Working directory | `/media/usuario2/ssd4TB1/kchopore_arabidopsis/run_transcriptome` |
| Conda env (pipeline) | `kchopore` (via `~/miniconda/etc/profile.d/conda.sh`) |
| Conda env (figures) | `viz` |
| System R | `/usr/bin/Rscript` (has DESeq2 / ggplot2 / patchwork / ragg) |
| NAS data (read-only) | `~/nas/Comun`, `~/nas/HTData` (rclone mounts) |
| Figure code | `FIGS_PAPER/code/` |
| Figure tidy data | `FIGS_PAPER/data/` |
| Poster assets | `POSTER/` |

**The single editable file is `config/config_transcriptome.yml`.** The `Snakefile` rules
are never touched for a new run.

---

## 1. Install / Docker

Two equivalent ways to run the QC + alignment pipeline: native conda, or the reproducible
Docker image (FAIR path).

### Conda (native)

```bash
ssh usuario2@156.35.42.17
source ~/miniconda/etc/profile.d/conda.sh
conda activate kchopore
```

### Docker (reproducible, FAIR)

```bash
cd ~/kchopore-anac017-drs
docker build -t kchopore-anac017-drs:latest .          # image already built
docker run --rm \
  -v $PWD:/workspace -v ~/nas:/nas -w /workspace \
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
~/bin/rclone mount nas:Comun          ~/nas/Comun  --read-only --vfs-cache-mode off --dir-cache-time 720h --daemon
~/bin/rclone mount nas:HTData_and_DBs ~/nas/HTData --read-only --vfs-cache-mode off --dir-cache-time 720h --daemon
ls ~/nas/Comun/Chus   # sanity check
```

### Prepare the workdir

```bash
conda activate kchopore
python3 ~/setup_transcriptome.py   # builds run_transcriptome with fastq/summaries/BAMs/ref symlinked
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
cd /media/usuario2/ssd4TB1/kchopore_arabidopsis/run_transcriptome
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
bash ~/run_deseq2.sh
# idxstats of the 12 BAMs -> counts -> DESeq2 -> ~/deseq2/out (PCA, volcanos, heatmap)
```

### 3.3 GO / KEGG enrichment (g:Profiler)

```bash
conda activate kchopore
Rscript ~/go_enrich.R      # -> ~/deseq2/go
```

### 3.4 m6A figures (ELIGOS2 already run by Chus)

```bash
conda run -n viz python ~/eligos_m6a_figures_v2.py
# -> ~/eligos_figs_v2 (motif, normalised rate, volcano, heatmap)
# baseExt0 files already extracted in ~/eligos_results
```

### 3.5 m6anet (orthogonal validation, subset by disk space)

```bash
# copy one eventalign locally and process with a local FIFO (faster than SMB):
cp "~/nas/Comun/Chus/Chus_DRS_Nanopolish_eventalign_m6Anet_NO BORRAR/anac017-1_C_R3_eventalign.txt.gz" /tmp/
bash ~/m6anet_stream.sh anac017-1_C_R3   # (adjust to read from /tmp)
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
- **Robust consensus:** a site is m6A in a condition if it holds in **≥3 of 5** writer mutants.
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

## 5. Figures and tidy data

Figure code lives in `FIGS_PAPER/code/`; tidy data in `FIGS_PAPER/data/` unless noted.
Final poster figures: rate = **1a lollipop**, AA/ANAC017 = **2a stacked**. The whole
biology scenario (Scenario B) is in R; the RRACH motif is in R too.

### How to re-run the figure scripts

| Script | Reads | Produces |
|--------|-------|----------|
| `code/deseq2_paper.R` | `data/transcript_counts_matrix.csv` | re-runs DESeq2 (vst/PCA), PNG+PDF; PCA legend outside |
| `code/m6a_figures.R` (run from `FIGS_PAPER/m6a/`) | `m6a_summary.csv`, `m6A_master_table.csv` | 1a lollipop (FINAL), 1c dotplot, 2a stacked (FINAL), 2b lollipop, 2c waffle, fig5 heatmap. Env vars `M6A_SUMMARY/M6A_MASTER/M6A_OUT` |
| `code/motif_figures.R` (run from `FIGS_PAPER/m6a/`) | `m6a_motif_kmers.csv` | ggseqlogo 2×2 |
| `code/go_figures.R` | `~/deseq2/go/*_GO.csv` (server) | 3 dotplots (interaction = dotplot, not barplot) |
| `code/track_figures.R` | `eligos_bedgraph/*.A.ESB_ctrl.bdg` | QC violin + AOX1a (Chr3:7906521-7908746) + CYP81D8 (Chr4:17569740-17571743) |
| `code/deseq2_paper.R` → `F0_dispersion` | `data/transcript_counts_matrix.csv` | dispersion (model health) |
| `code/panel_nanoplot.py` | 12-library NanoPlot images | PIL 4×3 montage; var `PLOT` for another view |
| `eligos_m6a_paper.py` | ELIGOS2 `*_baseExt0.txt` (server) | m6A source + dumps CSVs |
| `bdg_trackplots.py` | `eligos_bedgraph/*.bdg` | Python tracks |
| `go_paper.py`, `fig_integration.py` | see below | GO / integration panels |

### Figure → script → data (provenance)

#### Scenario B — anac017 biology

| Figure | Script | Data | Lang | Type / status |
|--------|--------|------|------|---------------|
| `deseq2/F1_PCA` | `deseq2_paper.R` (L147-164) | `data/transcript_counts_matrix.csv` | R | scatter · legend outside ✅ |
| `deseq2/F2_volcano_interaction` | `deseq2_paper.R` (L169-208) | idem | R | volcano (19 ANAC017-dependent tx) ✅ |
| `deseq2/F2_volcano_panel` | `deseq2_paper.R` (L206-208) | idem | R | panel (A) AA (B) genotype (C) interaction · polish |
| `deseq2/F3_MA_panel` | `deseq2_paper.R` (L237-238) | idem | R | panel · polish |
| `deseq2/F4_heatmap_top40` | `deseq2_paper.R` (L243-281) | idem | R | heatmap (z-score vst) ✅ |
| `m6a/fig1_motif_rrach` | `eligos_m6a_paper.py` | ELIGOS2 `*_baseExt0.txt` → `m6a/m6A_master_table.csv` | Py | logo ✅ |
| `m6a/fig2_rate_normalised` | `eligos_m6a_paper.py` | idem | Py | barplot (sites per 1,000 A) · improve |
| `m6a/fig3_upset_overlap` | `eligos_m6a_paper.py` | idem | Py | upset · polish |
| `m6a/fig6_aa_anac017` | `eligos_m6a_paper.py` | idem | Py | barplot (~86% AA-gained need ANAC017) · improve |
| `bdg_tracks/track_AT3G22370` | `bdg_trackplots.py` | `eligos_bedgraph/*.bdg` | Py | track AOX1a, 4 conditions ✅ |
| `bdg_tracks/track_AT4G37370` | `bdg_trackplots.py` | idem | Py | track CYP81D8, 4 conditions ✅ |
| `go/interaction_ANAC017dep_dotplot` | `go_paper.py` | `RESULTADOS/03_GO/interaction_ANAC017dep_GO.csv` | Py | dotplot (stress & hypoxia) ✅ — prefer over barplot |
| `go/interaction_ANAC017dep_barplot` | `go_paper.py` | idem | Py | barplot · discard, use dotplot |
| `integration/fig_integration_expr_m6a` | `fig_integration.py` | `RESULTADOS/02_DESeq2/D5_volcano_interaction.csv` + `m6a/m6A_master_table.csv` | Py | panel · ⚠️ descriptive only, not convergence |

#### Scenario A — K-CHOPORE works

| Figure | Script | Data | Lang | Action |
|--------|--------|------|------|--------|
| Hero architecture | (prompt v3, sci-diagram) | — | — | generate |
| QC `RESULTADOS/01_QC/nanocomp` | NanoComp (pipeline) | BAMs / fastq | — | pick 1 ✅ |
| QR repo | `POSTER/make_qr.py` | URL `github.com/biopelayo/...` | Py | ✅ |

### Poster figure blueprint (tidyplots, 4 figures)

**Global style.** Okabe-Ito condition palette: `WT_C=#0072B2 · WT_AA=#56B4E9 ·
anac017-1_C=#D55E00 · anac017-1_AA=#E69F00`. Significant `#BB5566` / non-sig `#BBBBBB`.
GO sources `GO:BP=#3B6CB7 · GO:MF=#2E8B6B · GO:CC=#C9892B · KEGG=#8E5AA6`. White background,
Arial sans-serif, top/right spines off, faint grid `#E8E8E8`. Export PNG 300 dpi + PDF with
`save_plot(..., bg="white")`. Panel labels A/B/C bold 11 pt, top-left; axis captions ≤6 words.

| Fig | Panels (tidy data) | Status |
|-----|--------------------|--------|
| **F1 · DRS quality + design** | A QC `qc_nanostats_summary.csv` (scatter or NanoComp image — NanoPlot dumps no per-read data) · B PCA `tidy_pca_coords.csv` (PC1 41%, PC2 30%) · C dispersion `tidy_dispersion.csv` | ✅ |
| **F2 · Differential expression** | A/B/C volcanos `tidy_de_{AA_in_WT,genotype,interaction}.csv` (151 / 93 / 19) · MA from same CSVs · D heatmap `tidy_heatmap_top40_long.csv` | ✅ |
| **F3 · m6A detection** | A ESB violin `eligos_bedgraph/*.A.ESB_ctrl.bdg` (pivot long) · B motif `m6a_motif_kmers.csv` (ggseqlogo) · C per-condition m6A volcano `m6a_volcano_percond.csv` | ✅ |
| **F4 · ANAC017-dependent m6A + function** | A 86% stacked bars `data/m6a_summary.csv` (aa_gain_wt=7202, anac017_dep=6225) · B track AOX1a (`*.bdg` subset Chr3:7,906,221-7,909,046) · C GO dotplot `RESULTADOS/03_GO/interaction_ANAC017dep_GO.csv` | ✅ |

All 12 panels have turnkey tidy data in `FIGS_PAPER/data/`.

---

## 6. Reproducibility

### NAS access

- Read-only rclone mounts, no sudo (see §2): `nas:Comun` → `~/nas/Comun`,
  `nas:HTData_and_DBs` → `~/nas/HTData`. Data are read in place; nothing is copied.
- Source DRS data and Chus's precomputed ELIGOS2 / eventalign products live under
  `~/nas/Comun/Chus`.

### Provenance

- Every poster figure is traced to its **script** and its **data file** in
  `FIGS_PAPER/INVENTARIO_PROVENANCE.md` (consolidated in §5). Figure code in
  `FIGS_PAPER/code/`, data in `FIGS_PAPER/data/`.
- The pipeline is config-driven: a run is fully described by
  `config/config_transcriptome.yml` plus the immutable `Snakefile`.

### GitHub + QR

- Repo mirrored at `github.com/biopelayo/...`; the poster QR is generated by
  `POSTER/make_qr.py` and links to the repo (Scenario A "Reproducible" panel).

### Versions and determinism

- Pinned tool versions: **Guppy 6.2.1** (basecalling), ELIGOS2, m6anet, DESeq2,
  gprofiler2; the full toolchain is captured in the Docker image
  `kchopore-anac017-drs:latest` for FAIR re-execution.
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
| `m6A_master_table.csv` | `eligos_m6a_paper.py` | consolidated m6A calls (figure source) |
| `m6a_volcano_percond.csv` | `eligos_m6a_paper.py` | per-condition A-site volcano (log2 oddR, −log10 p, hit) |
| `data.site_proba.csv` | m6anet | per-site m6A probability (signal-based) |
| `diffmod.table` | xpore | differential m6A (pending) |
| DESeq2 CSV per contrast | `run_deseq2.sh` | transcript-level DE (AA, genotype, interaction) |
| `*_GO.csv` | g:Profiler | GO:BP/MF/CC + KEGG enrichment of DE genes |
| MultiQC HTML | MultiQC | aggregate QC report |
| tidy CSVs in `FIGS_PAPER/data/` | figure scripts | per-panel plotting data (turnkey) |

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
