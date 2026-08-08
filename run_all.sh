#!/usr/bin/env bash
# Pipeline completo. Ejecutar desde la raiz del repositorio.
# Requiere que data/ contenga los ficheros de ADNI (ver data/README.md).
set -e
mkdir -p output figures
echo "== 00 extraer ADNIMERGE2 ==" ; Rscript src/00_extract_adnimerge2.R
echo "== 01 fusion ==="            ; python3 src/01_merge.py
echo "== 02 bootstrap B=1000 (~75 min) ==="
for i in 1 2 3 4 5; do python3 src/02_bootstrap.py 200; done
echo "== 03 agregacion ==="        ; python3 src/03_aggregate.py
echo "== 04 evidencia de cola ===" ; python3 src/04_tail_evidence.py
echo "== 05 estratificacion ==="   ; python3 src/05_stratified.py
echo "== 06 evidencia estratificada ===" ; python3 src/06_stratified_evidence.py
echo "== 07 solo-jitter ==="       ; python3 src/07_jitter_only.py 110
echo "== 08 figuras ==="           ; python3 src/08_figures.py
echo "== GCGM benchmark ==="       ; Rscript gcgm_benchmark.R
echo; echo "HECHO. Resultados en output/, figuras en figures/"
