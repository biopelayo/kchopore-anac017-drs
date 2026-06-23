# K-CHOPORE · Reporte de herramientas, settings y outputs (revisión)

Diseño 2×2: **WT / anac017-1** × **Control / Antimycin A**, 3 réplicas. DRS nanopore.
Dos vías: **transcriptoma** (FLAIR isoforms — eje epitranscriptómico, oficial) y **genoma** TAIR10 (QC/DE alternativo).
> Objetivo de este documento: que Pelayo revise settings/filtros de cada paso (bien / regular / mal).

| # | Tool | Paso | Settings / filtros usados | Output | Plots |
|---|------|------|---------------------------|--------|-------|
| 1 | **Guppy 6.2.1** (Chus) | Basecalling FAST5→FASTQ | `rna_r9.4.1_70bps_hac.cfg`, `--reverse_sequence 1 --u_substitution 1`, `--min_qscore 5` | fastq | — |
| 2 | **NanoFilt** | Filtrado lectura | config: `q≥7`, `length≥200` (en la corrida transcriptoma se usó fastq sin re-filtrar) | fastq filtrado | — |
| 3 | **NanoPlot** | QC por muestra | defaults | NanoStats.txt | ~15 plots/muestra (L vs Q, yield, histogramas) |
| 4 | **NanoComp** | QC comparativo | defaults, 12 muestras | report HTML | violins length/quality/N50/reads, histogramas overlay |
| 5 | **minimap2** | Mapeo | **transcriptoma**: `-ax map-ont -L`; **genoma**: `-ax splice -k14 --secondary=no --MD` | SAM | — |
| 6 | **samtools** | Filtrar/sort | `view -F 2324 -q 10` (primarios, MAPQ≥10), `sort`, `index` | BAM+bai | — |
| 7 | **samtools stats** | Stats align | flagstat + stats | txt | → MultiQC |
| 8 | **pycoQC** | QC señal+align | `--min_pass_qual 7` (summary + BAM) | HTML interactivo | QC reads/canales/tiempo |
| 9 | **FLAIR / StringTie** | Isoformas | FLAIR collapse (support≥3, stringent); StringTie min_cov 2.5 | isoforms.bed/gtf | (reutilizado de Chus) |
| 10 | **ELIGOS2** | m6A por error | `rna_mod`/`pair_diff_mod` `--pval 0.05 --oddR 5 --esb 0.2`; sobre transcriptoma flair; **vs 5 writer-mutants** (mta/mtb/fip37/vir/hakai) | baseExt0/1/2, combine, .A.filtered, BedGraph(IGV) | nativo: BedGraph + QC; **nuestros**: motivo RRACH, volcano, UpSet, heatmap |
| 11 | **m6anet** | m6A por señal | `dataprep --n_processes 1` (FIFO), `inference --num_iterations 1000` | data.site_proba.csv | (downstream propio) |
| 12 | **xpore** | m6A diferencial | pendiente (necesita eventalign 2 condiciones) | diffmod.table | — |
| 13 | **DESeq2** | Expresión diferencial | counts = `samtools idxstats` (transcript-level); diseño `~genotype+treatment+genotype:treatment`; filtro `rowSums≥10`; sig `padj<0.05 & |LFC|>1` | CSV por contraste | PCA, volcano×3, heatmap, dispersión |
| 14 | **g:Profiler** (gprofiler2) | GO/KEGG | `organism=athaliana`, `significant=TRUE`, GO:BP/MF/CC+KEGG; sobre genes DE (AT IDs) | GO.csv | barplot top términos |
| 15 | **MultiQC** | Agregado QC | defaults | report HTML | agregado |

## Criterios de llamada m6A (ELIGOS2, nuestro downstream)
- Sitio candidato m6A: `ref==A` & `pval<0.05` & `oddR<1` & `ESB_ctrl>ESB_test` (la muestra tiene más error que el writer-mutante → pierde m6A al quitar el writer).
- **Consenso robusto**: sitio m6A en una condición si cumple en **≥3 de 5** writers.
- Normalización entre condiciones: **tasa = sitios / 1.000 A testeadas** (corrige profundidad; clave porque WT_C tenía menos cobertura).

## Puntos a revisar (posibles mejoras / dudas)
1. **NanoFilt no aplicado realmente** en la vía transcriptoma (usamos fastq de Chus sin re-filtrar) → decidir si re-filtrar `q≥7 l≥200`.
2. **oddR**: el config ELIGOS pide `oddR 5` pero el `.A.filtered` de Chus parece pre-filtrado a `pval<0.05` sin umbral oddR estricto; nuestro consenso usa `oddR<1` (dirección m6A). Revisar si queremos endurecer (`|log2 oddR|` mínimo).
3. **DESeq2 a nivel transcrito** (no gen): los transcritos FLAIR de novo no tienen AT gene → para GO se pierden esos. Alternativa: agregar por gen vía solapamiento con AtRTD.
4. **m6anet** solo factible para subconjunto (espacio). Definir alcance.
5. **xpore** aún no corrido.
