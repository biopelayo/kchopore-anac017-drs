# FAIR_TODO — K-CHOPORE anac017 DRS

Plan de remediación FAIR. No ejecuta builds: lista cada cambio con la línea exacta a tocar y su prioridad.

> **Actualización 2026-06-25 — versiones reales capturadas del servidor.**
> Se capturaron las versiones exactas de los entornos conda del servidor de análisis
> (`usuario2@156.35.42.17`) y se guardaron en `envs/frozen/` (kchopore, m6anet, xpore,
> viz, pycoqc_env). El pipeline usa un entorno **primario `kchopore`** (Snakemake, mapeo,
> QC, isoformas, DESeq2/GO) más entornos aislados para `m6anet`, `xpore`, `pycoqc_env` y
> figuras (`viz`).
>
> **Resolución plotly / NanoPlot (corrige el supuesto de C5.1/Dockerfile):** NO hay
> conflicto real en el entorno que funciona. En `kchopore` coexisten **NanoPlot 1.46.1 +
> plotly 6.3.1 + pycoQC 2.5.2** sin romper (pycoQC se instala `--no-deps`, tolera plotly 6
> en ejecución). El pin antiguo `plotly==4.1.0` solo vive en el entorno aislado `pycoqc_env`
> como respaldo, sin compartir proceso con NanoPlot. Por eso `requirements.txt` y el
> Dockerfile fijan **plotly 6.3.1**, no 4.1.0.
>
> Estado por ítem marcado abajo con `[HECHO]`, `[PARCIAL]` o `[PENDIENTE-PELAYO]`.

Leyenda de prioridad:

- **P0** bloqueante de reproducibilidad (sin esto el build no es determinista).
- **P1** importante para FAIR (licencia, accesiones, semillas).
- **P2** recomendable, no bloqueante.

Ficheros afectados:

- `requirements.txt`
- `Dockerfile`
- `README` (placeholder de accesión)
- `LICENSE` (a crear)
- Scripts de R y de m6anet (semillas)

---

## C5 — Pinear versiones (reproducibilidad del build)

### C5.1 — `requirements.txt`: pasar todo de `>=` a `==` (P0) — `[HECHO]`

Hecho con las versiones REALES del entorno `kchopore` (no las mínimas de abajo):
numpy 2.3.4, pandas 2.3.3, scipy 1.16.3, mappy 2.24, h5py 3.15.1, pybedtools 0.12.0,
ncls 0.0.70, **plotly 6.3.1** (ver nota de resolución arriba), Jinja2 3.1.6,
MarkupSafe 3.0.3, tqdm 4.67.1, python-dateutil 2.9.0.post0, pytz 2025.2, tzdata 2025.2,
retrying 1.4.2, tabulate 0.9.0, setuptools 80.9.0, wheel 0.45.1.
`Cython` y `rpy2` quedan SIN pin: no aparecían instalados en ningún entorno capturado
(Cython es solo build; rpy2 vive en el stack pip del Docker). `[PENDIENTE-PELAYO]`: fijarlos
desde un build de imagen validado.

Tabla original (versiones mínimas, ya superada por las reales de arriba):

Cada línea con `>=` deja que el resolutor coja la última versión publicada en el momento del build, así que dos builds en fechas distintas dan entornos distintos. Hay que congelar a la versión concreta que se valide. Las versiones objetivo de abajo son las mínimas declaradas hoy; **antes de fijar, capturar las reales con `pip freeze` dentro de un build que funcione** y usar esas. Líneas exactas a cambiar:

```
setuptools>=69.0.0      ->  setuptools==69.0.0
numpy>=1.23.5           ->  numpy==1.23.5
pandas>=2.0.0           ->  pandas==2.0.0
scipy>=1.10.0           ->  scipy==1.10.0
mappy>=2.24             ->  mappy==2.24
h5py>=3.11.0            ->  h5py==3.11.0
pybedtools>=0.10.0      ->  pybedtools==0.10.0
ncls>=0.0.68            ->  ncls==0.0.68
plotly>=4.1.0           ->  plotly==4.1.0
Jinja2>=3.1.4           ->  Jinja2==3.1.4
MarkupSafe>=2.1.5       ->  MarkupSafe==2.1.5
tqdm>=4.66.5            ->  tqdm==4.66.5
rpy2>=3.5.0             ->  rpy2==3.5.0
python-dateutil>=2.9.0  ->  python-dateutil==2.9.0
pytz>=2024.1            ->  pytz==2024.1
tzdata>=2024.1          ->  tzdata==2024.1
retrying>=1.3.4         ->  retrying==1.3.4
tabulate>=0.9.0         ->  tabulate==0.9.0
```

Notas:

- `wheel` y `Cython` (líneas 11-12) no llevan versión: fijarlos también, p. ej. `wheel==0.43.0` y `Cython==3.0.10`.
- `plotly==4.1.0` choca con NanoPlot (que pide `plotly>=6.1.1`) y por eso pycoQC se instala con `--no-deps` en el Dockerfile (líneas 146-151). Mantener la separación: el `plotly` de `requirements.txt` queda para pycoQC; verificar que no rompe el entorno final.

### C5.2 — `Dockerfile`: pinear los pip install sin `==` (P0) — `[HECHO / PARCIAL]`

Fijado con las versiones reales observadas en los entornos del servidor (`envs/frozen/`):
m6anet **2.1.0** (no 2.0.1), xpore 2.1, pycoQC 2.5.2, NanoPlot 1.46.1, NanoComp 1.25.6,
NanoFilt 2.8.0, nanoget 1.19.3, nanomath 1.4.0, multiqc 1.32, pod5 0.3.39, pysam 0.23.3,
flair-brookslab 2.0.0, ont-bonito 1.1.0, snakemake **9.13.4** (con pulp 2.8.0),
logomaker 0.8.7, upsetplot 0.9.0, seaborn 0.13.2.
`[PENDIENTE-PELAYO]`: `torch`, `scikit-learn`, `rpy2`, `ont-fast5-api` quedan SIN pin porque
no se observó una versión instalada compatible con el Python 3.10 de la imagen (el torch del
entorno m6anet es 1.6.0, atado a Python 3.8). Fijarlos desde un build de imagen validado.
La línea del `sed list_solvers→listSolvers` queda como no-op (snakemake 9.x ya usa
`listSolvers`); marcada para eliminar tras validar el build.

Mismas versiones objetivo (notas originales, ya aplicadas con versiones reales):

- **Línea 117** `pip install --upgrade snakemake` (instala HEAD del índice). Cambiar a versión fija, p. ej. `pip install snakemake==7.32.4`. Ojo: el `sed` de la línea 118 parchea `list_solvers`→`listSolvers`, dependiente de versión; al fijar snakemake, revalidar que ese parche sigue aplicando o eliminarlo si la versión ya usa `listSolvers`.
- **Línea 123** `pip install ont-bonito || true` -> `pip install ont-bonito==0.8.1 || true` (ajustar a la validada).
- **Línea 159** `pip install --no-cache-dir torch scikit-learn` -> `torch==2.2.2 scikit-learn==1.4.2` (ajustar). `torch` sin pin es especialmente grave: cambia entre builds y arrastra CUDA.
- **Línea 158** `pip install --no-cache-dir --no-deps m6anet` -> `m6anet==2.0.1` (el comentario ya cita 2.0.1).
- **Líneas 164-165** `xpore` (dos intentos) -> `xpore==2.1` en ambos.
- **Línea 178** `pip install --no-cache-dir multiqc` -> `multiqc==1.21` (ajustar).
- **Líneas 139-144** NanoPlot/NanoComp/NanoFilt/nanoget/nanomath sin pin -> fijar cada uno (`NanoPlot==1.42.0`, etc.).
- **Línea 183** `pip install --no-cache-dir pod5` -> `pod5==0.3.10` (ajustar).
- **Líneas 188-192** `pysam`, `rpy2`, `ont-fast5-api` sin pin -> `pysam==0.22.0`, `rpy2==3.5.16`, `ont-fast5-api==4.1.3` (ajustar). `tabulate==0.9.0` ya está fijado, dejar.
- **Línea 204** `pip install --no-cache-dir flair-brookslab` -> `flair-brookslab==2.0.0` (ajustar).
- **Línea 151** pycoQC `--no-deps` sin pin -> `pycoQC==2.5.2`.
- **Línea 254** `logomaker upsetplot seaborn "kaleido==0.2.1"` -> fijar `logomaker==0.8`, `upsetplot==0.9.0`, `seaborn==0.13.2` (kaleido ya pineado).

### C5.3 — Pin de commit/tag para los clones que cogen HEAD (P0) — `[HECHO / PARCIAL]`

- **eligos2: `[HECHO]`** fijado al SHA real del servidor
  `b205a5cec34cd2974e5ddae78f7fcf7beb49d9c8` (`git describe` → `v2.1.0-8-gb205a5c`,
  commit 2024-01-29), que es el que usa el análisis (`/home/usuario2/eligos2`).
- **nanopolish: `[PARCIAL / PENDIENTE-PELAYO]`** se añadió `git checkout v0.14.0` +
  `submodule update`, pero NO se observó binario de nanopolish en ningún entorno conda del
  servidor del que leer la versión real (corre dentro del Docker). v0.14.0 es el tag estable
  documentado, no un install observado. Confirmar con `nanopolish --version` del binario real.

Detalle original:

Ahora mismo `git clone` sin `--branch` ni `checkout` trae la rama por defecto en su estado actual: el build deja de ser reproducible en cuanto el upstream haga un commit.

- **Línea 101** nanopolish:
  ```
  RUN git clone --recursive https://github.com/jts/nanopolish.git /opt/nanopolish && \
      cd /opt/nanopolish && \
  ```
  añadir checkout a un tag fijo tras el clone, p. ej.:
  ```
  RUN git clone --recursive https://github.com/jts/nanopolish.git /opt/nanopolish && \
      cd /opt/nanopolish && \
      git checkout v0.14.0 && git submodule update --init --recursive && \
      make -j1 && \
  ```
  (usar el tag/commit que se valide; v0.14.0 es la última release estable conocida).

- **Línea 170** eligos2:
  ```
  RUN git clone https://gitlab.com/piroonj/eligos2.git /home/eligos2 && \
  ```
  eligos2 no publica tags fiables: fijar por commit SHA. Cambiar a:
  ```
  RUN git clone https://gitlab.com/piroonj/eligos2.git /home/eligos2 && \
      cd /home/eligos2 && git checkout <SHA_VALIDADO> && \
  ```
  (anotar el SHA del `HEAD` que usa el análisis actual: `git -C /home/eligos2 rev-parse HEAD`).

### C5.4 — Binarios descargados por URL (P2, ya están pineados, solo verificar) — `[HECHO / nota]`

minimap2 vía apt: añadida nota en el Dockerfile (línea ~59). El servidor corre minimap2 2.30;
apt no permite fijar versión. `[PENDIENTE-PELAYO, P2, opcional]`: sustituir el `apt install`
por una descarga de release de minimap2 2.30 si se quiere determinismo total. El resto
(samtools 1.19, picard 2.25.7, dorado 0.8.0, guppy 6.1.5, stringtie 2.2.1, yq 4.9.8) ya van
pineados por URL.

Detalle original:

samtools 1.19, minimap2 (apt, sin pin de versión — línea 59), picard 2.25.7, dorado 0.8.0, guppy 6.1.5, stringtie 2.2.1, yq 4.9.8. Todos llevan versión salvo `minimap2` vía `apt` (línea 59): si se quiere determinismo total, instalar minimap2 por release con versión fija en vez de `apt-get install minimap2`.

---

## M13 — Metadatos FAIR

### M13.1 — LICENSE MIT (P1) — `[HECHO]`

`LICENSE` (MIT, 2026, titulares Pelayo G. de Lena, Jesús Pascual, Mario F. Fraga, Luis
Valledor / Univ. de Oviedo) ya existe en el repo, y el README ya lleva su sección de licencia.
Añadido además `CITATION.cff` con `license: MIT` y los mismos autores (versión, fecha de
release y DOI marcados como TODO hasta cortar release).

Detalle original:

### M13.2 — Semillas deterministas (P0 reproducibilidad de resultados, P1 FAIR) — `[PENDIENTE-PELAYO]`

`[PENDIENTE-PELAYO]`: NO tocado en este pase. Fijar semillas implica editar código de análisis
(scripts de R, llamada a `m6anet inference`, scripts numpy/torch propios), que queda fuera de
alcance de esta tarea de pinning (no se toca código de análisis ni figuras). Aplicar
`set.seed(42)` / `--seed 42` / `np.random.seed(42)` / `torch.manual_seed(42)` como detalla abajo.

Detalle original:

- **R**: en cada script de R que muestree, modele o use cualquier RNG (DESeq2/apeglm usan optimización con componente aleatorio; `EnhancedVolcano`/heatmaps no, pero `samplesizeCMH` y cualquier bootstrap sí), añadir al inicio:
  ```r
  set.seed(42)
  ```
  Donde haya funciones con argumento de aleatoriedad, pasar `seed=TRUE` / el valor de semilla explícito. Localizar los `.R` en `scripts/` y aplicarlo en la cabecera de cada uno.

- **m6anet**: en la llamada a `m6anet inference` (o `m6anet-run_inference`) del Snakefile/scripts, añadir el flag de semilla:
  ```
  m6anet inference ... --seed 42
  ```
  (m6anet expone `--seed`; fijarlo para que el muestreo del modelo sea reproducible).

- Revisar también xpore y cualquier `numpy`/`torch` en scripts propios: fijar `np.random.seed(42)` y `torch.manual_seed(42)` donde haya estocasticidad.

### M13.3 — Placeholder de accesión ENA/GEO en README (P1) — `[HECHO el placeholder / PENDIENTE-PELAYO el número]`

El README ya incluye la línea de disponibilidad de datos con placeholder visible:
`ENA PRJEBXXXXXX / GEO GSEXXXXXX (to be assigned)`. `[PENDIENTE-PELAYO]`: sustituir por las
accesiones reales antes de publicar (requiere número real, no se inventa).

Detalle original:

Añadir en el README un bloque de disponibilidad de datos con placeholder visible hasta tener el número real, p. ej.:

```
## Data availability

Raw Direct RNA-seq data: ENA accession `PRJEBXXXXXX` (TBD).
Processed counts / modification tables: GEO accession `GSEXXXXXX` (TBD).
```

Dejar el marcador `XXXXXX`/`TBD` explícito para que no se olvide rellenarlo antes de publicar. Recordatorio de tono: el proyecto es follow-up, no usar "first/primer" en la redacción del README.

---

## Resumen de prioridades

| Item | Prioridad | Fichero | Estado |
|------|-----------|---------|--------|
| C5.1 requirements `>=`→`==` | P0 | requirements.txt | HECHO (versiones reales; Cython/rpy2 pendientes) |
| C5.2 pip sin pin en Dockerfile | P0 | Dockerfile | HECHO/PARCIAL (torch/scikit-learn/rpy2/ont-fast5-api sin pin) |
| C5.3 pin commit/tag nanopolish + eligos2 | P0 | Dockerfile | eligos2 HECHO (SHA real); nanopolish PARCIAL (tag v0.14.0, sin binario que verificar) |
| C5.4 minimap2 apt sin versión | P2 | Dockerfile | NOTA añadida; opcional |
| M13.1 LICENSE MIT | P1 | LICENSE / README / CITATION.cff | HECHO |
| M13.2 set.seed / --seed | P0 | scripts R, Snakefile/m6anet | PENDIENTE-PELAYO (toca código de análisis) |
| M13.3 placeholder ENA/GEO | P1 | README | placeholder HECHO; número PENDIENTE-PELAYO |
| Captura de versiones reales | P0 | envs/frozen/*.txt | HECHO (kchopore, m6anet, xpore, viz, pycoqc_env) |
| environment.yml (name+channels+pins) | P1 | envs/environment.yml | HECHO |

### Pendiente real de Pelayo (decisiones / datos que no se pueden inventar)
1. **Accesiones ENA/GEO reales** (M13.3) — sustituir `PRJEBXXXXXX` / `GSEXXXXXX`.
2. **Tag de release + DOI** — fijar `version` y `date-released` en `CITATION.cff`, cortar release y archivar (Zenodo) para el DOI; confirmar la URL canónica del repo en `CITATION.cff`.
3. **Semillas deterministas** (M13.2) — editar scripts de R / m6anet (fuera del alcance de pinning).
4. **Pins residuales desde un build de imagen validado**: `torch`, `scikit-learn`, `rpy2`, `ont-fast5-api`, `Cython`, y confirmar `nanopolish` (binario) y minimap2 2.30 por release si se quiere determinismo total.
