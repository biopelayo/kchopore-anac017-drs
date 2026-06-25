# K-CHOPORE · anac017 DRS

Reproducible **Direct RNA nanopore (DRS)** pipeline for the **m6A epitranscriptome and transcriptome of the mitochondrial retrograde response** in *Arabidopsis thaliana*.

**Design (2×2):** genotype (WT vs *anac017-1*) × treatment (Control vs Antimycin A), 3 biological replicates. Antimycin A inhibits mitochondrial complex III and triggers retrograde signalling, of which **ANAC017** is the master regulator.

## Pipeline (Snakemake)
Basecalling (Guppy) → splice/transcriptome mapping (minimap2) → QC (NanoPlot, NanoComp, pycoQC, MultiQC) → isoforms (FLAIR, StringTie) → **m6A detection** by ELIGOS2 (error-based, validated against m6A-writer mutants *mta/mtb/fip37/vir/hakai*) and m6anet (signal-level) → differential expression (DESeq2) → functional enrichment (g:Profiler). Containerised (Docker) for FAIR reproducibility.

## Repo layout
```
Snakefile                 workflow rules
config/                   config_*.yml (the single control point)
scripts/                  analysis + figure scripts (DESeq2, ELIGOS2/m6A, GO, setup)
Dockerfile                full software stack
requirements.txt          python deps (installed in the image)
METHODS_TOOLS.md          per-tool settings, filters and outputs (review sheet)
RUN_CLI.md                step-by-step CLI runbook
```

## Quick start
```bash
docker build -t kchopore-anac017-drs:latest .
docker run --rm -v $PWD:/workspace -v /path/to/data:/data -w /workspace kchopore-anac017-drs:latest \
  snakemake --configfile config/config_transcriptome.yml --cores 12 --rerun-triggers mtime --keep-going
```
Edit only `config/config_*.yml` (samples, conditions, data paths, `run_*` module flags) for a new run. See `RUN_CLI.md`.

## Data
Raw DRS data and intermediate signal files (FAST5/eventalign) are large and are **not** stored here; see `RUN_CLI.md` for mounting the lab NAS or fetching from a public repository (ENA/GEO).

Raw DRS data: ENA PRJEBXXXXXX / GEO GSEXXXXXX (to be assigned).

## Author
Pelayo González de Lena Rodríguez · Cancer Epigenetics & Nanomedicine Lab (FINBA) / Systems Biology Lab (Univ. Oviedo).
Follow-up study; methods build on the K-CHOPORE pipeline.

## License
MIT — see [`LICENSE`](LICENSE). If you use this pipeline, please cite it via [`CITATION.cff`](CITATION.cff).

## Environment / reproducibility
Exact software versions from the working analysis environment are frozen in
[`envs/frozen/`](envs/frozen/) (one file per conda env: `kchopore`, `m6anet`, `xpore`,
`viz`, `pycoqc_env`). The primary pipeline env is specified in
[`envs/environment.yml`](envs/environment.yml); the full software stack is built by the
`Dockerfile`.
