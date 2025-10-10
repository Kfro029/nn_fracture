#!/bin/bash
# generature_random.sh - updated
# Генерация набора vtk + json метаданных для каждого созданного файла

set -euo pipefail

BASE_DIR="/home/faki/one_fracture_chimera"
FRACTURES_DIR="$BASE_DIR/fractures"

# Количество файлов для генерации
NUM_FILES=2

# Общие параметры генерации
MIN_N=1
MAX_N=10
VTK_TYPE="structured"   # или "polydata"
SPACING="1.0"
SAMPLES=1               # если >1 fracture_generator_vtk создаст суффиксы _0000 и т.д.

echo "=== ГЕНЕРАЦИЯ $NUM_FILES ФАЙЛОВ ТРЕЩИН (VTK + JSON META) ==="

mkdir -p "$FRACTURES_DIR"

for ((i=1; i<=$NUM_FILES; i++)); do
    RANDOM_N=$(( MIN_N + RANDOM % (MAX_N - MIN_N + 1) ))
    RANDOM_ID=$((1000 + RANDOM % 9000))
    BASE_NAME="fractures_${RANDOM_ID}_${RANDOM_N}tr"

    OUT_FILE="$FRACTURES_DIR/${BASE_NAME}.vtk"
    META_FILE="$FRACTURES_DIR/${BASE_NAME}.json"

    echo "[$i/$NUM_FILES] Генерация: $OUT_FILE (N=$RANDOM_N)"
    
    # передаём --out и --meta. при необходимости --seed и др.
    python3 fracture_generator_vtk.py \
        --xmin -200 --xmax 200 \
        --ymin -200 --ymax 200 \
        --N "$RANDOM_N" \
        --out "$OUT_FILE" \
        --meta "$META_FILE" \
        --vtk-type "$VTK_TYPE" \
        --spacing "$SPACING" \
        --samples "$SAMPLES"

    echo "Сохранено: $OUT_FILE"
    echo "Мета: $META_FILE"
    echo "----------------------------------------"
done

echo "Генерация завершена! Создано $NUM_FILES файлов в $FRACTURES_DIR/"

