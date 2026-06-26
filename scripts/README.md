# Scripts

Helper scripts called by the Snakemake pipeline or run by hand for downstream analysis.
The pipeline itself is driven by the top-level `Snakefile`; the only one of these scripts it
invokes directly is `run_deseq2.R`.

## Top level

| Script | Purpose |
|--------|---------|
| `run_deseq2.R` | DESeq2 2x2 differential expression (called by the `run_deseq2` rule). Counts from `samtools idxstats`, design `~genotype+treatment+genotype:treatment`. Produces PCA, volcanoes, heatmap. |
| `go_enrich.R` | GO / KEGG enrichment of DE genes with gprofiler2 (organism `athaliana`). |
| `eligos_m6a_figures_v2.py` | m6A figures from ELIGOS2 output: RRACH motif, normalised rate, volcano, heatmap. |
| `m6anet_stream.sh` | Run m6anet on one eventalign file through a local FIFO (low memory). |
| `setup_transcriptome.py` | Build the working directory: symlink fastq, summaries, BAMs and reference for a transcriptome run. |
| `setup_directories.sh` | Create the expected output directory tree. |
| `generate_samples.sh` | Generate the sample list for the config from the input files. |
| `map_with_minimap2.py` | minimap2 mapping wrapper (transcriptome `map-ont`, genome `splice`). |
| `make_bams.sh` | Filter, sort and index alignments into BAM + bai. |
| `manual_kchopo.sh` | End-to-end manual run, mirroring the Snakemake steps. |
| `probar_pipe.sh` | Smoke-test driver for the pipeline. |

## Subfolders

- `alignment/` — `minimap2_mapping.py`, the mapping step used by the alignment rules.
- `pipelines/` — standalone shell drivers (`00_pipeline.sh`, `01_pipeline.sh`) and
  `sort_and_index_bam.py`.
