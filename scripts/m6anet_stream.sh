#!/usr/bin/env bash
# m6anet dataprep+inference leyendo el eventalign COMPLETO del NAS por streaming (FIFO),
# sin descargar el .gz gigante ni re-correr nanopolish.
# Uso: m6anet_stream.sh <sample> [gz_basename]
set -uo pipefail
S="$1"
GZBASE="${2:-${S}_eventalign.txt.gz}"
EVDIR="$HOME/nas/Comun/Chus/Chus_DRS_Nanopolish_eventalign_m6Anet_NO BORRAR"
GZ="$EVDIR/$GZBASE"
OUT="$HOME/m6anet_anac017/$S"
FIFO="$HOME/m6anet_anac017/${S}.fifo"
mkdir -p "$OUT"
rm -f "$FIFO"; mkfifo "$FIFO"
source /home/usuario2/miniconda/etc/profile.d/conda.sh

echo "[m6anet] $(date) inicio $S"
echo "[m6anet] eventalign: $GZ"
if [ ! -f "$GZ" ]; then echo "ERROR: no existe $GZ"; exit 1; fi

# productor: descomprime el eventalign del NAS al FIFO
zcat "$GZ" > "$FIFO" &
ZPID=$!

# consumidor: m6anet dataprep (n_processes 1 -> lectura secuencial, compatible con FIFO)
conda run -n m6anet m6anet dataprep --eventalign "$FIFO" --out_dir "$OUT" --n_processes 1
RC=$?
wait "$ZPID" 2>/dev/null
rm -f "$FIFO"
echo "M6ANET_DATAPREP_EXIT=$RC $(date)"

if [ "$RC" -eq 0 ]; then
  conda run -n m6anet m6anet inference --input_dir "$OUT" --out_dir "$OUT" --n_processes 4 --num_iterations 1000
  echo "M6ANET_INFERENCE_EXIT=$? $(date)"
  echo "Resultado: $OUT/data.site_proba.csv"
fi
