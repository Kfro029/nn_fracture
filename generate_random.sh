set -euo pipefail

BASE_DIR="/home/faki/one_fracture_chimera"
FRACTURES_DIR="$BASE_DIR/fractures_1500"

# Количество файлов для генерации
NUM_FILES=1500

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
    # RANDOM_ID=$((1000 + RANDOM % 9000))
    
    # Генерация случайных параметров для названия файла
    MEAN_ANGLE=$(( RANDOM % 180 ))           # средний угол наклона трещин (0-179)
    ANGLE_SPREAD=$(( 10 + RANDOM % 80 ))     # угол разброса трещин (10-89)
    MEAN_LENGTH=$(( 10 + RANDOM % 100 ))     # средняя длина трещин (10-109)
    LENGTH_VARIANCE=$(( 1 + RANDOM % 20 ))   # дисперсия средней длины (1-20)

    BASE_NAME="fractures_${i}_${RANDOM_N}_${MEAN_ANGLE}_${ANGLE_SPREAD}_${MEAN_LENGTH}_${LENGTH_VARIANCE}"

    OUT_FILE="$FRACTURES_DIR/${BASE_NAME}.txt"
    META_FILE="$FRACTURES_DIR/${BASE_NAME}.json"

    # echo "[$i/$NUM_FILES] Генерация: $OUT_FILE (N=$RANDOM_N)"
    
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

    # echo "Сохранено: $OUT_FILE"
    # echo "Мета: $META_FILE"
    # echo "----------------------------------------"
done

#echo "Генерация завершена! Создано $NUM_FILES файлов в $FRACTURES_DIR/"
