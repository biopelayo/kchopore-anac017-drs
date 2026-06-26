# Configuration

YAML config files for the pipeline. A run is fully described by one config file plus the
immutable `Snakefile`; each file lists the samples, their control/treatment labels, the input
file paths, and the `run_*` flags that switch modules on or off.

| File | Role |
|------|------|
| `config_transcriptome.yml` | Canonical config. The transcriptome (FLAIR isoform) axis used for the anac017 analysis. This is the file to edit for a new run. |
| `config_anac017.yml` | Earlier genome-axis config for the same anac017 samples. Kept for reference. |
| `config.yml` | Annotated template with every key documented and the full WT / anac017 sample set. Copy from here when starting a new design. |

Edit only the canonical config; the `Snakefile` rules are not touched.
