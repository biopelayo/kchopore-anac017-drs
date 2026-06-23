# 🧬 K-CHOPORE · Comandos Snakemake

Guía rápida para correr el pipeline **a mano** (lo que `scripts/probar_pipe.sh` hace por dentro).
Workdir de referencia: `…/kchopore_arabidopsis/run_transcriptome`.

---

## 0 · Preparar la sesión (una vez)

```bash
ssh usuario2@156.35.42.17
source ~/miniconda/etc/profile.d/conda.sh
conda activate kchopore
cd /media/usuario2/ssd4TB1/kchopore_arabidopsis/run_transcriptome
```

> Todo lo demás se lanza **desde este directorio**.

---

## 1 · Dry-run (plan, NO ejecuta) — SIEMPRE primero

```bash
snakemake -n \
  --snakefile Snakefile \
  --configfile config/config_transcriptome.yml \
  --rerun-triggers mtime
```

| Flag | Qué hace |
|------|----------|
| `-n` | *dry-run*: muestra qué haría, sin tocar nada |
| `--configfile` | el fichero que controla todo (muestras, módulos, rutas) |
| `--rerun-triggers mtime` | decide qué rehacer **solo por fecha** → no recomputa lo pesado ya hecho |

✅ Si dice `Nothing to be done` = todo está al día (reutiliza, no recomputa).

---

## 2 · Ejecutar de verdad

```bash
snakemake \
  --snakefile Snakefile \
  --configfile config/config_transcriptome.yml \
  --cores 12 \
  --rerun-triggers mtime \
  --keep-going
```

| Flag | Qué hace |
|------|----------|
| *(sin `-n`)* | ejecuta de verdad |
| `--cores 12` | núcleos en paralelo (hay 40; 12–16 va sobrado). Alias: `-j 12` |
| `--keep-going` | si una regla falla, sigue con las demás |

---

## 3 · Forzar una regla para VERLA correr

```bash
snakemake \
  --snakefile Snakefile \
  --configfile config/config_transcriptome.yml \
  --cores 4 --rerun-triggers mtime \
  --forcerun multiqc results/multiqc/multiqc_report.html
```

`--forcerun <regla>` re-ejecuta esa regla aunque ya esté hecha (útil para demostrar el pipeline).

---

## 4 · Pedir un resultado concreto (un *target*)

```bash
snakemake … --cores 4 --rerun-triggers mtime \
  results/nanoplot/WT_C_R1/NanoStats.txt
```

Pones como argumento el **fichero de salida** que quieres → Snakemake hace solo lo necesario para producirlo.

---

## 🔧 Flags útiles

| Flag | Para qué |
|------|----------|
| `-p` / `--printshellcmds` | imprime el comando real de cada paso |
| `-r` / `--reason` | explica *por qué* ejecuta cada regla |
| `--unlock` | desbloquea el directorio tras un corte (luego relanzas) |
| `--dag \| dot -Tpng > dag.png` | dibuja el grafo de dependencias |
| `-j N` | alias de `--cores N` |

---

## ✏️ Lo único que editas para una corrida nueva

El archivo **`config/config_transcriptome.yml`**:

- `samples:` — lista de muestras
- `conditions:` — control / treatment por muestra
- `input_files:` — rutas de fastq, summaries, referencia
- `params: run_*:` — qué módulos correr (`true`/`false`)

Las reglas del `Snakefile` **no se tocan**.

---

## 🥇 Regla de oro

```
-n  (ver el plan)   →   si te convence, quitas el -n y lanzas
```

Deja siempre `--rerun-triggers mtime` para reutilizar lo ya calculado y no rehacer los pasos pesados.
